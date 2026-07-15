/**
 * @file   bin_chunk_loader.cpp
 * @brief  File selection, bin loading, and chunk-slicing for the
 *         noise-marking GUI.
 */

#include "gui_handler.h"
#include "chart_utils.hpp"
#include "gui_handler.h"

#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QProgressDialog>
#include <QApplication>
#include <QEventLoop>
#include <algorithm>
#include <cstring>
#include <limits>
#include <vector>

void noise_marking_gui::setFileSource(const QString& filePath) {
    loadSelectedFile(filePath);
}

void noise_marking_gui::loadSelectedFile(const QString& filePath) {
    if (!m_binFilePath.isEmpty()) {
        m_genExc.filePath = m_binFilePath;
        m_fileMarkings[m_binFilePath] = m_genExc;
    }

    m_binFilePath = filePath;
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly)) return;

    // Header layout: 1 scalar (sleep epoch length) + 36 upsampled sizes + 36 raw sizes + 36 native rates + 36 upsampled rates + 1 sleep count.
    constexpr int n_fields = 1 + 4 * NUM_CHANNELS + 1;
    static_assert(FILE_HEADER_SIZE == n_fields * 4, "FILE_HEADER_SIZE must match 1 + 4*N + 1 slot count");

    uint32_t raw32[n_fields] = {};
    file.read(reinterpret_cast<char*>(raw32), sizeof(raw32));
    file.close();

    double sleepEpoch = static_cast<double>(raw32[0]);
    m_sleepSR = (sleepEpoch > 0) ? (1.0 / sleepEpoch) : 0;

    constexpr int kSizesRawBase = 1 + NUM_CHANNELS;
    constexpr int kNativeRatesBase = kSizesRawBase + NUM_CHANNELS;
    constexpr int kUpRatesBase = kNativeRatesBase + NUM_CHANNELS;
    constexpr int kSleepCountIdx = kUpRatesBase + NUM_CHANNELS;

    for (int i = 0; i < NUM_CHANNELS; ++i) {
        upsampled_channel_sizes[i] = raw32[1 + i];
        raw_channel_sizes[i] = raw32[kSizesRawBase + i];
        std::memcpy(&channel_native_rates[i], &raw32[kNativeRatesBase + i], sizeof(float));
        std::memcpy(&channel_upsampled_rates[i], &raw32[kUpRatesBase + i], sizeof(float));
    }
    total_sleep_samples = raw32[kSleepCountIdx];

    if (m_fileMarkings.contains(filePath)) {
        m_genExc = m_fileMarkings[filePath];
        m_noiseManager = std::make_unique<annotation_handler>();
        for (int i = 0; i < m_genExc.noiseExc.size(); ++i) {
            double sr = sampleRateForSignal(m_genExc.data_type[i]);
            m_noiseManager->addSegment(
                static_cast<int>(m_genExc.noiseExc[i].first * sr),
                static_cast<int>(m_genExc.noiseExc[i].second * sr),
                m_genExc.data_type[i].toStdString(),
                m_genExc.marking_type[i].toStdString(),
                sr);
        }
    }
    else {
        m_genExc = GenExcStruct();
        m_genExc.filePath = filePath;
        m_noiseManager = std::make_unique<annotation_handler>();
    }

    current_start_time = 0.0;
    m_markArmed = false;
    for (const QString& lbl : markableChannelLabels()) cancelMarking(lbl);

    setWindowTitle("Marking: " + QFileInfo(filePath).fileName());
    loadChunkFromFile(0);
}

void noise_marking_gui::handleBrowseFile() {
    QString startDir;
    if (!m_binFilePath.isEmpty())
        startDir = QFileInfo(m_binFilePath).absolutePath();
    else if (!m_cfg.bin_file_path.empty())
        startDir = QString::fromStdString(m_cfg.bin_file_path);

    QString binPath = QFileDialog::getOpenFileName(
        this, "Select Bin File", startDir,
        "Converted bin files (*.bin);;All files (*)");

    if (binPath.isEmpty() || binPath == m_binFilePath) return;

    if (!QFileInfo(binPath).isReadable()) {
        QMessageBox::warning(this, "Cannot open file",
            QString("File is not readable:\n%1").arg(binPath));
        return;
    }

    QWidget* prevFocus = QApplication::focusWidget();
    QProgressDialog progress(this);
    progress.setWindowTitle("Loading");
    progress.setLabelText(QString("Loading %1\u2026").arg(QFileInfo(binPath).fileName()));
    progress.setRange(0, 0); progress.setCancelButton(nullptr);
    progress.setMinimumDuration(0); progress.setWindowModality(Qt::WindowModal);
    progress.setAutoClose(false); progress.setAutoReset(false);
    progress.setFocusPolicy(Qt::NoFocus); progress.show();
    QApplication::processEvents(QEventLoop::ExcludeUserInputEvents);

    loadSelectedFile(binPath);

    progress.close();
    if (prevFocus && !qobject_cast<QLineEdit*>(prevFocus))
        prevFocus->setFocus(Qt::OtherFocusReason);
    else
        this->setFocus(Qt::OtherFocusReason);
    QApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

bool noise_marking_gui::loadChunkFromFile(uint64_t chunkIndex) {
    QFile file(m_binFilePath);
    if (!file.open(QIODevice::ReadOnly)) return false;
    current_chunk_index = chunkIndex;

    uint64_t chanUpOffset[NUM_CHANNELS], chanRawOffset[NUM_CHANNELS];
    uint64_t running = 0;
    for (int i = 0; i < NUM_CHANNELS; ++i) {
        chanUpOffset[i] = running; running += upsampled_channel_sizes[i];
        chanRawOffset[i] = running; running += raw_channel_sizes[i] * 2;
    }
    const uint64_t sleepByteOffset = running;


    auto loadSignal = [&](QVector<double>& dest, int chIdx) {
        const double rate = channel_upsampled_rates[chIdx];
        uint64_t totalSamples = upsampled_channel_sizes[chIdx];
        uint64_t perChunk = static_cast<uint64_t>(seconds_in_memory_at_once * rate);
        uint64_t start = chunkIndex * perChunk;
        uint64_t count = (totalSamples > start)
            ? std::min(perChunk, totalSamples - start) : 0;
        dest.resize(static_cast<int>(count));
        file.seek(FILE_HEADER_SIZE + (chanUpOffset[chIdx] + start) * sizeof(double));
        file.read(reinterpret_cast<char*>(dest.data()), count * sizeof(double));
        };

    auto loadRaw = [&](QVector<QPointF>& dest, int chIdx) {
        dest.clear();
        const uint64_t totalPairs = raw_channel_sizes[chIdx];
        if (totalPairs <= 1) {
            if (totalPairs == 1) {
                double pair[2] = { -1.0, -1.0 };
                file.seek(FILE_HEADER_SIZE + chanRawOffset[chIdx] * sizeof(double));
                file.read(reinterpret_cast<char*>(pair), 2 * sizeof(double));
                dest.append(QPointF(pair[0], pair[1]));
            }
            return;
        }
        const float nativeHz = channel_native_rates[chIdx];
        if (nativeHz <= 0.0f) return;

        // Slice the raw block by TIME, not by a pair count. The raw pairs are
        // stored with x = row-index time (row index / grid rate) and are
        // monotonic in x, but they are NOT uniformly dense at nativeHz -- gaps
        // and 2x-packing corrections mean a chunk holds fewer pairs than
        // 8h * nativeHz. The old "firstPair = chunkIndex * 8h * nativeHz"
        // slicing assumed density and so overshot on every chunk after the
        // first (landing past the real chunk data -> no scatter). Instead,
        // load exactly the pairs whose stored time falls in this chunk's
        // [chunkStart, chunkStart + 8h) window. This matches loadSignal's
        // time-based chunking so raw and upsampled always cover the same span.
        const double chunkStartT =
            static_cast<double>(chunkIndex) * seconds_in_memory_at_once;
        const double chunkEndT = chunkStartT + seconds_in_memory_at_once;

        const qint64 baseBytes = FILE_HEADER_SIZE
            + static_cast<qint64>(chanRawOffset[chIdx]) * sizeof(double);

        // Binary-search the first pair with t >= chunkStartT. Pairs are (t, v)
        // doubles, 16 bytes each, monotonic in t.
        auto pairTimeAt = [&](uint64_t idx) -> double {
            double t = 0.0;
            file.seek(baseBytes + static_cast<qint64>(idx) * 16);
            file.read(reinterpret_cast<char*>(&t), sizeof(double));
            return t;
            };
        uint64_t lo = 0, hi = totalPairs;
        while (lo < hi) {
            const uint64_t mid = (lo + hi) / 2;
            if (pairTimeAt(mid) < chunkStartT) lo = mid + 1; else hi = mid;
        }
        const uint64_t firstPair = lo;
        if (firstPair >= totalPairs) return;

        // Read forward in a bounded block, keeping pairs until t >= chunkEndT.
        // Cap the read so a huge chunk doesn't allocate unboundedly; loop if
        // the window spans more than one block.
        if (!file.seek(baseBytes + static_cast<qint64>(firstPair) * 16)) return;
        constexpr uint64_t kBlockPairs = 1u << 20;   // 1M pairs = 16 MB per read
        uint64_t remaining = totalPairs - firstPair;
        bool done = false;
        dest.reserve(static_cast<int>(std::min<uint64_t>(remaining,
            static_cast<uint64_t>(std::numeric_limits<int>::max()))));
        std::vector<double> buf;
        while (remaining > 0 && !done) {
            const uint64_t thisBlock = std::min(kBlockPairs, remaining);
            try { buf.resize(thisBlock * 2); }
            catch (const std::bad_alloc&) { return; }
            const qint64 got = file.read(reinterpret_cast<char*>(buf.data()),
                static_cast<qint64>(thisBlock) * 16);
            if (got <= 0) break;
            const uint64_t gotPairs = static_cast<uint64_t>(got) / 16;
            for (uint64_t k = 0; k < gotPairs; ++k) {
                const double t = buf[k * 2];
                if (t >= chunkEndT) { done = true; break; }
                // Store chunk-local x (subtract chunk start), matching the
                // upsampled block, which the renderer plots at index/rate =
                // chunk-local seconds. current_start_time resets to 0 per chunk.
                dest.append(QPointF(t - chunkStartT, buf[k * 2 + 1]));
            }
            remaining -= gotPairs;
            if (gotPairs < thisBlock) break;   // short read = EOF
        }
        };

    loadSignal(m_ecg1, CH_ECG1);   loadSignal(m_ecg2, CH_ECG2);
    loadSignal(m_ecg3, CH_ECG3);   loadSignal(m_ppg, CH_PPG);
    loadSignal(m_accelX, CH_ACCEL_X); loadSignal(m_accelY, CH_ACCEL_Y);
    loadSignal(m_accelZ, CH_ACCEL_Z); loadSignal(m_cvp, CH_CVP);
    loadSignal(m_resp, CH_RESP);    loadSignal(m_abp, CH_ABP);
    loadSignal(m_art, CH_ART); loadSignal(m_artPulm, CH_ART_PULM);
    loadSignal(m_temp, CH_TEMP);    loadSignal(m_marker, CH_MARKER);
    loadSignal(m_pacemaker, CH_PACEMAKER_EVENT);

    loadRaw(m_ecg1Raw, CH_ECG1);   loadRaw(m_ecg2Raw, CH_ECG2);
    loadRaw(m_ecg3Raw, CH_ECG3);   loadRaw(m_ppgRaw, CH_PPG);
    loadRaw(m_abpRaw, CH_ABP);    loadRaw(m_accelXRaw, CH_ACCEL_X);
    loadRaw(m_accelYRaw, CH_ACCEL_Y); loadRaw(m_accelZRaw, CH_ACCEL_Z);
    loadRaw(m_respRaw, CH_RESP);   loadRaw(m_cvpRaw, CH_CVP);
    loadRaw(m_artRaw, CH_ART); loadRaw(m_artPulmRaw, CH_ART_PULM);
    loadRaw(m_tempRaw, CH_TEMP);   loadRaw(m_markerRaw, CH_MARKER);
    loadRaw(m_pacemakerRaw, CH_PACEMAKER_EVENT);

    // The raw (t, v) block is stored with x = row-index time (row index /
    // grid rate), the same clock the upsampled block uses, so the raw scatter
    // already lands exactly on the upsampled trace -- no re-timing needed. (The
    // old code rewrote x to per-channel ordinal time here, which desynced the
    // two blocks at every gap; that rewrite, and the gap_indicator it fed, have
    // been removed.)

    {
        uint64_t perChunk = static_cast<uint64_t>(seconds_in_memory_at_once * m_sleepSR);
        uint64_t start = chunkIndex * perChunk;
        uint64_t count = (total_sleep_samples > start)
            ? std::min(perChunk, total_sleep_samples - start) : 0;
        m_sleepStages.resize(static_cast<int>(count));
        file.seek(FILE_HEADER_SIZE + (sleepByteOffset + start) * sizeof(double));
        file.read(reinterpret_cast<char*>(m_sleepStages.data()), count * sizeof(double));
    }
    file.close();
    current_start_time = 0;

    m_activeChannels.clear();
    auto markActive = [this](const QString& label, const QVector<double>& data) {
        bool missing = is_missing_signal(data);
        if (auto* cv = chartViewForSignalLabel(label)) cv->setVisible(!missing);
        if (!missing) m_activeChannels.insert(label);
        };
    const bool bittium = (m_cfg.dataset_type == "BITTIUM");

    markActive("ECG1", m_ecg1);
    markActive("ECG2", m_ecg2);
    markActive("ECG3", m_ecg3);
    markActive("PPG", m_ppg);
    markActive("ACCEL", m_accelX);
    markActive("ART", m_art);
    markActive("ART_PULM", m_artPulm);


    bool anyAccel = !is_missing_signal(m_accelX) || !is_missing_signal(m_accelY) || !is_missing_signal(m_accelZ);

    if (ui->abp_axis)
        ui->abp_axis->setVisible(!bittium && !is_missing_signal(m_abp));
    if (ui->ppg_ampogram_axis)
        ui->ppg_ampogram_axis->setVisible(!is_missing_signal(m_ppg));

    if (bittium) {
        if (ui->cvp_eeg_axis)
            ui->cvp_eeg_axis->setVisible(!is_missing_signal(m_temp));
        if (ui->hyp_resp_axis)
            ui->hyp_resp_axis->setVisible(!is_missing_signal(m_marker));
        if (ui->pacemaker_axis)
            ui->pacemaker_axis->setVisible(!is_missing_signal(m_pacemaker));
    }
    else {
        if (!anyAccel) markActive("ABP", m_abp);
        if (ui->pacemaker_axis)
            ui->pacemaker_axis->setVisible(false);   // BITTIUM-only chart
        if (ui->cvp_eeg_axis)
            ui->cvp_eeg_axis->setVisible(m_cvpRaw.size() >= 2);
        if (ui->hyp_resp_axis) {
            ui->hyp_resp_axis->setVisible(
                sleep_data_present(m_sleepStages)
                || !is_missing_signal(m_resp) || anyAccel);
        }
    }

    updateMarkingButtons();
    ampogram();
    handle_data_plot();
    setupHypnogram();
    updateAmpogramCursor();

    uint64_t ecgPerChunk = static_cast<uint64_t>(seconds_in_memory_at_once * channel_upsampled_rates[CH_ECG1]);
    ui->prev8hours->setEnabled(chunkIndex > 0);
    ui->next8hours->setEnabled(
        (chunkIndex * ecgPerChunk + m_ecg1.size()) < upsampled_channel_sizes[CH_ECG1]);
    return true;
}

void noise_marking_gui::on_next8hours_clicked() {
    /*
        If the 'w' key or the next 8 hours button is clicked, see if there is another 8 hour chunk of data to load, if so, load it.
        If there is additional data less than 8 hours, that data is loaded.
    */
    const uint64_t ecgPerChunk = static_cast<uint64_t>(seconds_in_memory_at_once * channel_upsampled_rates[CH_ECG1]);
    const uint64_t nextStart = (current_chunk_index + 1) * ecgPerChunk;
    if (nextStart >= upsampled_channel_sizes[CH_ECG1]) return;   // no data ahead
    resetUnpinnedGains();
    loadChunkFromFile(current_chunk_index + 1);
}
void noise_marking_gui::on_prev8hours_clicked() {
    if (current_chunk_index > 0) {
        resetUnpinnedGains(); loadChunkFromFile(current_chunk_index - 1);
    }
}