/**
 * @file   bin_chunk_loader.cpp
 * @brief  File selection, bin loading, and chunk-slicing for the
 *         noise-marking GUI.
 */

#include "gui_handler.h"
#include "chart_utils.hpp"

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
        const uint64_t perChunk = static_cast<uint64_t>(seconds_in_memory_at_once * (double)nativeHz);
        const uint64_t firstPair = chunkIndex * perChunk;
        if (firstPair >= totalPairs) return;
        const uint64_t count = std::min(perChunk, totalPairs - firstPair);
        if (count == 0) return;

        const qint64 baseBytes = FILE_HEADER_SIZE
            + static_cast<qint64>(chanRawOffset[chIdx]) * sizeof(double);
        if (!file.seek(baseBytes + static_cast<qint64>(firstPair) * 16)) return;

        std::vector<double> buf;
        try { buf.resize(count * 2); }
        catch (const std::bad_alloc&) { return; }
        const qint64 got = file.read(reinterpret_cast<char*>(buf.data()),
            static_cast<qint64>(count) * 16);
        if (got <= 0) return;
        const uint64_t gotPairs = static_cast<uint64_t>(got) / 16;
        dest.reserve(static_cast<int>(
            std::min<uint64_t>(gotPairs,
                static_cast<uint64_t>(std::numeric_limits<int>::max()))));
        for (uint64_t k = 0; k < gotPairs; ++k)
            dest.append(QPointF(buf[k * 2], buf[k * 2 + 1]));
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

    if (m_gapIndicator) m_gapIndicator->rescan();

    auto rewriteRawToIndexTime = [&](QVector<QPointF>& raw, int chIdx) {
        if (raw.size() < 2) return;
        const float nativeHz = channel_native_rates[chIdx];
        if (nativeHz <= 0.0f) return;
        const double dt = 1.0 / nativeHz;
        for (int i = 0; i < raw.size(); ++i) raw[i].setX(i * dt);
        };
    rewriteRawToIndexTime(m_ecg1Raw, CH_ECG1);
    rewriteRawToIndexTime(m_ecg2Raw, CH_ECG2);
    rewriteRawToIndexTime(m_ecg3Raw, CH_ECG3);
    rewriteRawToIndexTime(m_ppgRaw, CH_PPG);
    rewriteRawToIndexTime(m_abpRaw, CH_ABP);
    rewriteRawToIndexTime(m_tempRaw, CH_TEMP);
    rewriteRawToIndexTime(m_markerRaw, CH_MARKER);
    rewriteRawToIndexTime(m_accelXRaw, CH_ACCEL_X);
    rewriteRawToIndexTime(m_accelYRaw, CH_ACCEL_Y);
    rewriteRawToIndexTime(m_accelZRaw, CH_ACCEL_Z);
    rewriteRawToIndexTime(m_respRaw, CH_RESP);
    rewriteRawToIndexTime(m_cvpRaw, CH_CVP);
    rewriteRawToIndexTime(m_pacemakerRaw, CH_PACEMAKER_EVENT);
    rewriteRawToIndexTime(m_artRaw, CH_ART);
    rewriteRawToIndexTime(m_artPulmRaw, CH_ART_PULM);

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
    restoreMarkingMarkers();

    uint64_t ecgPerChunk = static_cast<uint64_t>(seconds_in_memory_at_once * channel_upsampled_rates[CH_ECG1]);
    ui->prev8hours->setEnabled(chunkIndex > 0);
    ui->next8hours->setEnabled(
        (chunkIndex * ecgPerChunk + m_ecg1.size()) < upsampled_channel_sizes[CH_ECG1]);
    return true;
}

void noise_marking_gui::on_next8hours_clicked() {
    resetUnpinnedGains(); loadChunkFromFile(current_chunk_index + 1);
}
void noise_marking_gui::on_prev8hours_clicked() {
    if (current_chunk_index > 0) { resetUnpinnedGains(); loadChunkFromFile(current_chunk_index - 1); }
}