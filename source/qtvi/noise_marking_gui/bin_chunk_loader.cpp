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

    constexpr int kNumChannels = NUM_CHANNELS;
    constexpr int kNumFields = 4 + 3 * kNumChannels + 1;
    static_assert(FILE_HEADER_SIZE == kNumFields * 4,
        "FILE_HEADER_SIZE must match 4 + 3*N + 1 slot count");

    uint32_t raw32[kNumFields] = {};
    file.read(reinterpret_cast<char*>(raw32), sizeof(raw32));
    file.close();

    m_ecgSR = static_cast<double>(raw32[0]);
    m_boolSR = static_cast<double>(raw32[1]);
    m_ppgSR = m_ecgSR;
    double sleepEpoch = static_cast<double>(raw32[3]);
    m_sleepSR = (sleepEpoch > 0) ? (1.0 / sleepEpoch) : 0;

    constexpr int kSizesUpBase = 4;
    constexpr int kSizesRawBase = kSizesUpBase + kNumChannels;
    constexpr int kNativeRatesBase = kSizesRawBase + kNumChannels;
    constexpr int kSleepCountIdx = kNativeRatesBase + kNumChannels;

    for (int i = 0; i < kNumChannels; ++i) {
        upsampled_channel_sizes[i] = raw32[kSizesUpBase + i];
        raw_channel_sizes[i] = raw32[kSizesRawBase + i];
        std::memcpy(&channel_native_rates[i], &raw32[kNativeRatesBase + i], sizeof(float));
    }
    total_sleep_samples = raw32[kSleepCountIdx];

    if (m_fileMarkings.contains(filePath)) {
        m_genExc = m_fileMarkings[filePath];
        m_noiseManager = std::make_unique<annotation_handler>(m_ecgSR);
        for (int i = 0; i < m_genExc.noiseExc.size(); ++i) {
            double sr = sampleRateForSignal(m_genExc.data_type[i]);
            m_noiseManager->addSegment(
                static_cast<int>(m_genExc.noiseExc[i].first * sr),
                static_cast<int>(m_genExc.noiseExc[i].second * sr),
                m_genExc.data_type[i].toStdString(),
                m_genExc.marking_type[i].toStdString());
        }
    }
    else {
        m_genExc = GenExcStruct();
        m_genExc.filePath = filePath;
        m_noiseManager = std::make_unique<annotation_handler>(m_ecgSR);
    }

    current_start_time = 0.0;
    m_markAllMode = MarkAllMode::None;
    single_ecg_marker_clicked = false;
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

    auto rateForChannel = [this](int chIdx) -> double {
        if (chIdx == CH_MARKER || chIdx == CH_TEMP || chIdx == CH_PACEMAKER)
            return m_boolSR;
        if (chIdx >= CH_EKG_OFF && chIdx <= CH_EEG3_OFF) return m_boolSR;
        if (chIdx == CH_OXSTATUS || chIdx == CH_SPO2 || chIdx == CH_HR)
            return m_boolSR;
        return m_ecgSR;
        };

    auto loadSignal = [&](QVector<double>& dest, int chIdx) {
        double   sr = rateForChannel(chIdx);
        uint64_t totalSamples = upsampled_channel_sizes[chIdx];
        uint64_t perChunk = static_cast<uint64_t>(seconds_in_memory_at_once * sr);
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
    loadSignal(m_accelZ, CH_ACCEL_Z); loadSignal(m_cvp, CH_PRES);
    loadSignal(m_resp, CH_RESP);    loadSignal(m_abp, CH_ABP);

    loadRaw(m_ecg1Raw, CH_ECG1);   loadRaw(m_ecg2Raw, CH_ECG2);
    loadRaw(m_ecg3Raw, CH_ECG3);   loadRaw(m_ppgRaw, CH_PPG);
    loadRaw(m_abpRaw, CH_ABP);    loadRaw(m_accelXRaw, CH_ACCEL_X);
    loadRaw(m_accelYRaw, CH_ACCEL_Y); loadRaw(m_accelZRaw, CH_ACCEL_Z);
    loadRaw(m_respRaw, CH_RESP);   loadRaw(m_cvpRaw, CH_PRES);

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
    rewriteRawToIndexTime(m_accelXRaw, CH_ACCEL_X);
    rewriteRawToIndexTime(m_accelYRaw, CH_ACCEL_Y);
    rewriteRawToIndexTime(m_accelZRaw, CH_ACCEL_Z);
    rewriteRawToIndexTime(m_respRaw, CH_RESP);
    rewriteRawToIndexTime(m_cvpRaw, CH_PRES);

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
        bool missing = isMissingSignal(data);
        if (auto* cv = chartViewForSignalLabel(label)) cv->setVisible(!missing);
        if (!missing) m_activeChannels.insert(label);
        };
    markActive("ECG1", m_ecg1); markActive("ECG2", m_ecg2);
    markActive("ECG3", m_ecg3); markActive("PPG", m_ppg);

    bool anyAccel = !isMissingSignal(m_accelX)
        || !isMissingSignal(m_accelY) || !isMissingSignal(m_accelZ);
    if (!anyAccel) markActive("ABP", m_abp);

    if (ui->accel_or_abp_axis)
        ui->accel_or_abp_axis->setVisible(!isMissingSignal(m_abp));
    if (ui->ppg_ampogram_axis)
        ui->ppg_ampogram_axis->setVisible(!isMissingSignal(m_ppg));
    if (ui->cvp_axis)
        ui->cvp_axis->setVisible(!isMissingSignal(m_cvp));
    if (ui->hyp_accel_resp_axis) {
        ui->hyp_accel_resp_axis->setVisible(
            sleepDataPresent(m_sleepStages)
            || !isMissingSignal(m_resp) || anyAccel);
    }

    updateAllChannelButtonStates();
    ampogram();
    handle_data_plot();
    setupHypnogram();
    updateAmpogramCursor();
    restoreMarkingMarkers();

    uint64_t ecgPerChunk = static_cast<uint64_t>(seconds_in_memory_at_once * m_ecgSR);
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