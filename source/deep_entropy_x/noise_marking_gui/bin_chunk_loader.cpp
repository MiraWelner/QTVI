/**
 * @file   bin_chunk_loader.cpp
 * @brief  File selection, bin loading, and chunk-slicing for the
 *         noise-marking GUI.
 */

#include "gui_handler.h"
#include "chart_utils.hpp"
#include "annotation_types.hpp"
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
#include <cmath>
#include <fstream>
#include <filesystem>
#include <limits>
#include <vector>

void noise_marking_gui::setFileSource(const QString& filePath) {
    loadSelectedFile(filePath);
}

namespace {
    GenExcStruct read_noise_markings_bin(const std::filesystem::path& path, const QString& filePath) {
        namespace nm = noise_markings;
        GenExcStruct g;
        g.filePath = filePath;

        const nm::RowsResult rr = nm::loadRows(path.string());
        if (!rr.read) {
            if (!rr.error.empty())
                std::fprintf(stderr, "[noise-markings] %s: %s\n",
                    path.string().c_str(), rr.error.c_str());
            return g;
        }
        if (rr.legacy)
            std::fprintf(stderr,  "[noise-markings] %s is old lecacy file\n",  path.string().c_str());

        for (std::size_t i = 0; i < rr.rows.size(); ++i) {
            const nm::Row& row = rr.rows[i];
            const char* chan = nm::channel_for_code(row.channel_code);
            if (!chan) {
                std::fprintf(stderr, "[noise-markings] %s row %zu: unknown "
                    "channel code %d, skipped\n", path.string().c_str(), i,
                    static_cast<int>(row.channel_code));
                continue;
            }
            const char* type = nullptr;
            for (const auto& t : annotation_types::noise_types)
                if (t.code == static_cast<int>(row.annotation_code)) {
                    type = t.label; break;
                }
            if (!type) {
                std::fprintf(stderr, "[noise-markings] %s row %zu: unknown "
                    "marking code %d, skipped\n", path.string().c_str(), i,
                    static_cast<int>(row.annotation_code));
                continue;
            }

            g.appendMarking(row.start_sec, row.end_sec, chan, type,
                row.threshold, row.blanking_ms);
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
            GenExcStruct g = read_noise_markings_bin(nb, filePath);
            if (!g.marks.isEmpty()) m_fileMarkings[filePath] = g;
        }
    }

    // markable channel whose meaning drifts as the operator scrolls.
    m_vcgCfg.ortho = vcg::OrthoBasis{};
    m_vcgCfg.orthoAcc.reset();
    // Same reasoning for the cached lead-polarity correction: a sign flip
    // decided for the PREVIOUS file must not leak into this one.
    m_vcgCfg.leadSign[0] = m_vcgCfg.leadSign[1] = m_vcgCfg.leadSign[2] = 1;
    {
        const QString stem = QFileInfo(filePath).completeBaseName();
        const QString dir = QString::fromStdString(m_cfg.vcg_output);
        m_vcgCfg.basisCsvPath = (dir + stem + "_vcg_basis.csv").toStdString();
        m_vcgCfg.basisCsvSubject = stem.toStdString();
    }

    if (m_fileMarkings.contains(filePath)) {
        // No replay into a second store any more: m_genExc IS the store.
        m_genExc = m_fileMarkings[filePath];
        rehydrateParamOverrides();
    }
    else {
        m_genExc = GenExcStruct();
        m_genExc.filePath = filePath;
        m_thresholdOverrides.clear();
        m_blankingOverrides.clear();
        m_invertOverrides.clear();
    }

    current_start_time = 0.0;
    m_markArmed = false;
    for (const QString& lbl : markableChannelLabels()) cancelMarking(lbl);

    setWindowTitle("Marking: " + QFileInfo(filePath).fileName());
    loadChunkFromFile(0);
    autoDetectLeadPolarity(); //technically you should be able to figure out lead polarity via vcg but in practice this has never actually worked
}

void noise_marking_gui::rehydrateParamOverrides() {
    m_thresholdOverrides.clear();
    m_blankingOverrides.clear();
    m_invertOverrides.clear();

    // ParamOverride still keys on QString because the renderer compares it
    // against chart labels; that conversion happens once per override here,
    // not once per marking per frame.
    for (const Marking& m : m_genExc.marks) {
        if (m.type == annotation_types::kParamEditLabel) {
            const QString ch = QString::fromStdString(m.channel);
            if (!std::isnan(m.threshold))
                m_thresholdOverrides.append(ParamOverride{ ch, m.start, m.end, m.threshold });
            if (!std::isnan(m.blanking))
                m_blankingOverrides.append(ParamOverride{ ch, m.start, m.end, m.blanking });
        }
        else if (m.type == annotation_types::kInvertEditLabel) {
            m_invertOverrides.append(ParamOverride{
                QString::fromStdString(m.channel), m.start, m.end, 1.0 });
        }
    }
}

void noise_marking_gui::handleBrowseFile() {
    QString startDir;
    if (!m_binFilePath.isEmpty())
        startDir = QFileInfo(m_binFilePath).absolutePath();
    else if (!m_cfg.input_path.empty())
        startDir = QString::fromStdString(m_cfg.input_path);

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
    progress.setLabelText(QString("Loading %1…").arg(QFileInfo(binPath).fileName()));
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

        const double chunkStartT = static_cast<double>(chunkIndex) * seconds_in_memory_at_once;
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
        // Reserve for THIS chunk, not for the rest of the channel: remaining
        // counts every pair from here to EOF, so on chunk 0 of a long
        // recording reserve(remaining) asked for the whole raw block up front.
        // nativeHz * chunk seconds bounds what the window can hold.
        const uint64_t expectedPairs = static_cast<uint64_t>(
            static_cast<double>(nativeHz) * seconds_in_memory_at_once) + 1;
        dest.reserve(static_cast<int>(std::min<uint64_t>(
            std::min(remaining, expectedPairs),
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

    // One table, so the upsampled and raw loads cannot drift apart. These were
    // two hand-maintained lists of the same 19 channels in different orders;
    // adding a channel to one and forgetting the other gives a trace with no
    // scatter, or a scatter with no trace.
    struct ChannelLoad { QVector<double>* up; QVector<QPointF>* raw; int ch; };
    const ChannelLoad kChannels[] = {
        { &m_ecg1,      &m_ecg1Raw,      CH_ECG1            },
        { &m_ecg2,      &m_ecg2Raw,      CH_ECG2            },
        { &m_ecg3,      &m_ecg3Raw,      CH_ECG3            },
        { &m_ppg,       &m_ppgRaw,       CH_PPG             },
        { &m_accelX,    &m_accelXRaw,    CH_ACCEL_X         },
        { &m_accelY,    &m_accelYRaw,    CH_ACCEL_Y         },
        { &m_accelZ,    &m_accelZRaw,    CH_ACCEL_Z         },
        { &m_cvp,       &m_cvpRaw,       CH_CVP             },
        { &m_resp,      &m_respRaw,      CH_RESP            },
        { &m_abp,       &m_abpRaw,       CH_ABP             },
        { &m_art,       &m_artRaw,       CH_ART             },
        { &m_artPulm,   &m_artPulmRaw,   CH_ART_PULM        },
        { &m_temp,      &m_tempRaw,      CH_TEMP            },
        { &m_marker,    &m_markerRaw,    CH_MARKER          },
        { &m_pacemaker, &m_pacemakerRaw, CH_PACEMAKER_EVENT },
        { &m_flow,      &m_flowRaw,      CH_FLOW            },
        { &m_thor,      &m_thorRaw,      CH_THOR            },
        { &m_abdo,      &m_abdoRaw,      CH_ABDO            },
        { &m_spo2,      &m_spo2Raw,      CH_SPO2            },
    };
    for (const ChannelLoad& c : kChannels) loadSignal(*c.up, c.ch);
    for (const ChannelLoad& c : kChannels) loadRaw(*c.raw, c.ch);

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
    //load subsequent chunk of data
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