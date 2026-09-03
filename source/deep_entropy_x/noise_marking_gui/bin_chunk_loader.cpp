/**
 * @file   bin_chunk_loader.cpp
 * @brief  File selection, bin loading, and chunk-slicing for the
 *         noise-marking GUI.
 */

#include "gui_handler.h"
#include "chart_utils.hpp"
#include "annotation_types.hpp"
#include "notch_filter.hpp"
#include "user_annotation_handler.h"
#include "vcg_lead.hpp"

#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QProgressDialog>
#include <QApplication>
#include <QEventLoop>
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <limits>
#include <vector>

void noise_marking_gui::setFileSource(const QString& filePath) {
    loadSelectedFile(filePath);
}

namespace {
    // Read a noise-markings .bin. ONE format: magic, version, count, then rows
    // of noise_markings::kColumns doubles indexed by the Column enum. See
    // noise_markings in user_annotation_handler.h, which the writer shares.
    //
    // A file without the magic is REFUSED, not parsed. Those are pre-versioned
    // files whose rows are six bare doubles and which carry no threshold or
    // blanking values, so reading one would restore every parameter-edit span at
    // the config defaults and move the R peaks inside it -- looking entirely
    // correct while being wrong about the one thing the span existed to record.
    GenExcStruct readNoiseMarkingsBin(const std::filesystem::path& path,
        const QString& filePath) {
        namespace nm = noise_markings;
        GenExcStruct g;
        g.filePath = filePath;
        std::ifstream f(path, std::ios::binary);
        if (!f.is_open()) return g;

        char magic[sizeof(nm::kMagic)] = {};
        f.read(magic, sizeof(magic));
        if (!f || std::memcmp(magic, nm::kMagic, sizeof(magic)) != 0) {
            std::fprintf(stderr,
                "[noise-markings] %s has no magic header, so it predates "
                "parameter storage; refusing to read it. Those files carry no "
                "threshold or blanking values -- re-save from the CSV or "
                "re-mark.\n", path.string().c_str());
            return g;
        }

        uint32_t version = 0;
        f.read(reinterpret_cast<char*>(&version), sizeof(version));
        if (!f) return g;
        // A version mismatch is a hard stop, not a best-effort read: guessing
        // the row stride is how spans land at the wrong times while looking
        // plausible.
        if (version != nm::kVersion) {
            std::fprintf(stderr,
                "[noise-markings] %s is version %u; this build reads version %u "
                "only.\n", path.string().c_str(), version, nm::kVersion);
            return g;
        }

        uint64_t count = 0;
        f.read(reinterpret_cast<char*>(&count), sizeof(count));
        if (!f) return g;

        for (uint64_t i = 0; i < count; ++i) {
            double row[nm::kColumns];
            f.read(reinterpret_cast<char*>(row), sizeof(row));
            if (!f) {
                std::fprintf(stderr,
                    "[noise-markings] %s ended after %llu of %llu rows\n",
                    path.string().c_str(),
                    static_cast<unsigned long long>(i),
                    static_cast<unsigned long long>(count));
                break;
            }

            // Both lookups read the tables the writer wrote from, so a miss here
            // means the file names something this build does not know. Worth a
            // line each rather than a silent skip -- silent skipping is what hid
            // the writer/reader channel-map divergence before the tables were
            // unified.
            const char* chan = nm::channel_for_code(
                static_cast<uint8_t>(row[nm::kChannelCode]));
            if (!chan) {
                std::fprintf(stderr, "[noise-markings] row %llu: unknown channel "
                    "code %d, skipped\n", static_cast<unsigned long long>(i),
                    static_cast<int>(row[nm::kChannelCode]));
                continue;
            }
            const char* type = nullptr;
            for (const auto& t : annotation_types::noise_types)
                if (t.code == static_cast<int>(row[nm::kAnnotationCode])) {
                    type = t.label; break;
                }
            if (!type) {
                std::fprintf(stderr, "[noise-markings] row %llu: unknown marking "
                    "code %d, skipped\n", static_cast<unsigned long long>(i),
                    static_cast<int>(row[nm::kAnnotationCode]));
                continue;
            }

            g.appendMarking(row[nm::kStartSec], row[nm::kEndSec],
                QString::fromLatin1(chan), QString::fromLatin1(type),
                row[nm::kThreshold], row[nm::kBlankingMs]);
        }
        return g;
    }
}   // namespace

void noise_marking_gui::loadSelectedFile(const QString& filePath) {
    if (!m_binFilePath.isEmpty()) {
        m_genExc.filePath = m_binFilePath;
        m_fileMarkings[m_binFilePath] = m_genExc;
    }

    m_binFilePath = filePath;
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly)) return;

    // Header layout (v1, 592 bytes): version + n_channels + sleep-epoch length,
    // then 36 upsampled sizes, 36 raw sizes, 36 native rates, 36 upsampled rates,
    // and 1 sleep count. 'version' and 'n_channels' were added at the front, so
    // every per-channel array sits 2 slots (8 bytes) past the old 584-byte layout.
    constexpr int n_fields = 2 + 1 + 4 * NUM_CHANNELS + 1;   // version,n_channels + old layout = 148
    static_assert(FILE_HEADER_SIZE == n_fields * 4,
        "FILE_HEADER_SIZE must match version + n_channels + sleep + 4*N + sleep_count");

    uint32_t raw32[n_fields] = {};
    file.read(reinterpret_cast<char*>(raw32), sizeof(raw32));
    file.close();

    const uint32_t header_version = raw32[0];               // offset 0
    const uint32_t n_channels = raw32[1];               // offset 4
    if (header_version != BIN_HEADER_VERSION || static_cast<int>(n_channels) != NUM_CHANNELS) {
        // Legacy (pre-version) bins began directly with the sleep-epoch length
        // and have no version/n_channels, so their fields would be misread.
        // Refuse loudly rather than plot garbage; regenerate with file_to_bin.
        QMessageBox::warning(this, "Unsupported .bin format",
            QString("Expected header version %1 with %2 channels, but got "
                "version=%3, n_channels=%4.\nThis .bin predates the header "
                "change - re-run file_to_bin to regenerate it.")
            .arg(BIN_HEADER_VERSION).arg(NUM_CHANNELS)
            .arg(header_version).arg(n_channels));
        return;
    }

    double sleepEpoch = static_cast<double>(raw32[2]);      // offset 8
    m_sleepSR = (sleepEpoch > 0) ? (1.0 / sleepEpoch) : 0;

    constexpr int kSizesUpBase = 3;                     // after version, n_channels, sleep_epoch
    constexpr int kSizesRawBase = kSizesUpBase + NUM_CHANNELS;
    constexpr int kNativeRatesBase = kSizesRawBase + NUM_CHANNELS;
    constexpr int kUpRatesBase = kNativeRatesBase + NUM_CHANNELS;
    constexpr int kSleepCountIdx = kUpRatesBase + NUM_CHANNELS;

    for (int i = 0; i < NUM_CHANNELS; ++i) {
        upsampled_channel_sizes[i] = raw32[kSizesUpBase + i];
        raw_channel_sizes[i] = raw32[kSizesRawBase + i];
        std::memcpy(&channel_native_rates[i], &raw32[kNativeRatesBase + i], sizeof(float));
        std::memcpy(&channel_upsampled_rates[i], &raw32[kUpRatesBase + i], sizeof(float));
    }
    total_sleep_samples = raw32[kSleepCountIdx];

    // First time we see this file this session: if a saved noise-markings
    // .bin exists on disk, load it so prior markings are restored. In-session
    // edits (already in m_fileMarkings) take precedence and are not clobbered.
    if (!m_fileMarkings.contains(filePath) && !m_cfg.noise_data_path.empty()) {
        const std::filesystem::path nb =
            std::filesystem::path(m_cfg.noise_data_path)
            / (QFileInfo(filePath).completeBaseName().toStdString() + "_noise_markings.bin");
        if (std::filesystem::exists(nb)) {
            GenExcStruct g = readNoiseMarkingsBin(nb, filePath);
            if (!g.noiseExc.isEmpty()) m_fileMarkings[filePath] = g;
        }
    }

    // Per-file VCG state: clear the cached basis so the next chunk rebuilds
    // it for THIS recording, and tell vcg_lead where to file it. Both belong
    // here and not in loadChunkFromFile -- a basis rebuilt per chunk gives a
    // markable channel whose meaning drifts as the operator scrolls.
    m_vcgCfg.ortho = vcg::OrthoBasis{};
    m_vcgCfg.orthoAcc.reset();
    // Same reasoning for the cached lead-polarity correction: a sign flip
    // decided for the PREVIOUS file must not leak into this one.
    m_vcgCfg.leadSign[0] = m_vcgCfg.leadSign[1] = m_vcgCfg.leadSign[2] = 1;
    {
        const QString stem = QFileInfo(filePath).completeBaseName();
        const QString dir = QString::fromStdString(m_cfg.vcg_output);
        QDir().mkpath(dir);   // nothing upstream creates vcg_output/
        m_vcgCfg.basisCsvPath = (dir + stem + "_vcg_basis.csv").toStdString();
        m_vcgCfg.basisCsvSubject = stem.toStdString();
    }

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

    // Measured, not asked: now that m_ecg1/m_ecg2/m_ecg3 hold this file's
    // first chunk, check whether any of the three is polarity-inverted
    // relative to the other two and pre-set the corresponding checkbox.
    // Once per FILE (here), not once per chunk (loadChunkFromFile) --
    // lead polarity is a property of the recording, not of a time window.
    autoDetectLeadPolarity();
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

bool noise_marking_gui::loadChunkFromFile(uint64_t chunkIndex, bool resetScroll) {
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

    // Raw (t,v) x-values are absolute Unix-epoch milliseconds (channel 0,
    // sample 0 is the recording start). Convert to seconds-from-start so the
    // raw scatter shares the upsampled block's index/rate time axis. Legacy
    // bins (x already in seconds, no epoch anchor) pass through unchanged.
    double recStartEpochMs = 0.0;
    if (upsampled_channel_sizes[CH_TIMESTAMP] > 0) {
        file.seek(FILE_HEADER_SIZE + chanUpOffset[CH_TIMESTAMP] * sizeof(double));
        file.read(reinterpret_cast<char*>(&recStartEpochMs), sizeof(double));
    }
    const bool rawIsEpochMs = (recStartEpochMs > 1.0e9);   // plausible epoch(ms) => new format
    auto rawToLocalSec = [recStartEpochMs, rawIsEpochMs](double t) -> double {
        return rawIsEpochMs ? (t - recStartEpochMs) / 1000.0 : t;
        };


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
            if (rawToLocalSec(pairTimeAt(mid)) < chunkStartT) lo = mid + 1; else hi = mid;
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
                const double ts = rawToLocalSec(buf[k * 2]);   // epoch-ms -> s from start
                if (ts >= chunkEndT) { done = true; break; }
                // Store chunk-local x (subtract chunk start), matching the
                // upsampled block, which the renderer plots at index/rate =
                // chunk-local seconds. current_start_time resets to 0 per chunk.
                dest.append(QPointF(ts - chunkStartT, buf[k * 2 + 1]));
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
    loadSignal(m_flow, CH_FLOW);    loadSignal(m_thor, CH_THOR);
    loadSignal(m_abdo, CH_ABDO);    loadSignal(m_spo2, CH_SPO2);

    // Powerline notch (config-driven Hz, toggled by the repurposed "Notch
    // Filter" checkbox). Applied once per chunk load, zero-phase, so both the
    // display and downstream peak detection see the same cleaned signal.
    if (m_notchFilterEnabled && m_cfg.notch_filter_hz != 0) {
        auto notch = [&](QVector<double>& sig, double rateHz) {
            if (is_missing_signal(sig)) return;
            std::vector<double> tmp(sig.begin(), sig.end());
            tmp = notch_filter::apply(tmp, rateHz, m_cfg.notch_filter_hz);
            sig = QVector<double>(tmp.begin(), tmp.end());
            };
        notch(m_ecg1, channel_upsampled_rates[CH_ECG1]);
        notch(m_ecg2, channel_upsampled_rates[CH_ECG2]);
        notch(m_ecg3, channel_upsampled_rates[CH_ECG3]);
        notch(m_ppg, channel_upsampled_rates[CH_PPG]);
        notch(m_abp, channel_upsampled_rates[CH_ABP]);
        notch(m_art, channel_upsampled_rates[CH_ART]);
        notch(m_artPulm, channel_upsampled_rates[CH_ART_PULM]);
    }

    loadRaw(m_ecg1Raw, CH_ECG1);   loadRaw(m_ecg2Raw, CH_ECG2);
    loadRaw(m_ecg3Raw, CH_ECG3);   loadRaw(m_ppgRaw, CH_PPG);
    loadRaw(m_abpRaw, CH_ABP);    loadRaw(m_accelXRaw, CH_ACCEL_X);
    loadRaw(m_accelYRaw, CH_ACCEL_Y); loadRaw(m_accelZRaw, CH_ACCEL_Z);
    loadRaw(m_respRaw, CH_RESP);   loadRaw(m_cvpRaw, CH_CVP);
    loadRaw(m_artRaw, CH_ART); loadRaw(m_artPulmRaw, CH_ART_PULM);
    loadRaw(m_tempRaw, CH_TEMP);   loadRaw(m_markerRaw, CH_MARKER);
    loadRaw(m_pacemakerRaw, CH_PACEMAKER_EVENT);
    loadRaw(m_flowRaw, CH_FLOW);   loadRaw(m_thorRaw, CH_THOR);
    loadRaw(m_abdoRaw, CH_ABDO);   loadRaw(m_spo2Raw, CH_SPO2);

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
    if (resetScroll) current_start_time = 0;
    m_activeChannels.clear();
    auto markActive = [this](const QString& label, const QVector<double>& data) {
        bool missing = is_missing_signal(data);
        if (auto* cv = chartViewForSignalLabel(label)) cv->setVisible(!missing);
        if (!missing) m_activeChannels.insert(label);
        };
    const bool bittium = (m_cfg.dataset_type == "BITTIUM");
    const bool shhs = (m_cfg.dataset_type == "SHHS");

    markActive("ECG1", m_ecg1);
    markActive("ECG2", m_ecg2);
    markActive("ECG3", m_ecg3);
    refreshVcgFromLeadFlags();   // shared with the ecg_N_reverse toggled handler; sets VCG's chart visibility itself
    markActive("PPG", m_ppg);
    markActive("ACCEL", m_accelX);
    markActive("ART", m_art);
    markActive("ART_PULM", m_artPulm);


    bool anyAccel = !is_missing_signal(m_accelX) || !is_missing_signal(m_accelY) || !is_missing_signal(m_accelZ);

    if (ui->abp_axis)
        ui->abp_axis->setVisible(!bittium && !is_missing_signal(m_abp));
    // SaO2 has its own chart at the bottom of the main plot column. It is set
    // here rather than inside the per-dataset branches below because it is the
    // one non-markable chart that is NOT a shared slot: no other dataset has
    // anything to put in it, so leaving it visible would cost a stretch slot in
    // main_plots and show an empty chart.
    if (ui->spo2_shhs_plot)
        ui->spo2_shhs_plot->setVisible(shhs && !is_missing_signal(m_spo2));
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
    else if (shhs) {
        // SHHS has no CVP, RESP, marker or ABP channel, and its sleep staging
        // owns hyp_resp_axis, so the two remaining non-markable slots carry the
        // PSG context a reviewer wants while marking ECG noise: cvp_eeg_axis
        // takes AIRFLOW/THOR/ABDO together, pacemaker_axis takes SaO2. The chart
        // titles come from determine_which_nonmarkable_charts_to_plot, so the
        // slot NAMES are the only thing still saying CVP/EEG and PACEMAKER.
        const bool anyResp = !is_missing_signal(m_flow)
            || !is_missing_signal(m_thor) || !is_missing_signal(m_abdo);
        if (ui->cvp_eeg_axis)   ui->cvp_eeg_axis->setVisible(anyResp);
        if (ui->pacemaker_axis) ui->pacemaker_axis->setVisible(false);   // BITTIUM-only chart
        if (ui->hyp_resp_axis)  ui->hyp_resp_axis->setVisible(sleep_data_present(m_sleepStages));
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
    if (resetScroll) handle_data_plot();
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