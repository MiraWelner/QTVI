/**
 * @file   file_to_bin.cpp
 * @brief  Implementation of the source-to-bin conversion declared in
 *         file_to_bin.hpp. EDF and .dat paths share the same on-disk
 *         channel layout but pull samples from very different shapes of
 *         input: EDF channels are uniform native-rate streams; .dat
 *         columns are sparse populated cells inside a uniform monitor
 *         row grid. Each channel is resampled to its OWN target rate
 *         (cfg.<chan>_upsample_rate) for the upsampled block, with a
 *         (t, v) raw block preserving the channel's native sampling.
 *
 *         All dataset-specific channel names live in config_entry and
 *         are populated by apply_dataset_specific_channel_labels(). This
 *         file knows nothing about MESA / Bittium / CHAOS conventions --
 *         if a label field is empty (or a channel has no configured
 *         rate) that slot writes a missing-channel placeholder.
 */

#include "file_to_bin.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

extern "C" {
#include "file_format_parsing/edflib.h"
}
#include "file_format_parsing/pugixml.hpp"

// Polyphase resampler + zero-phase filtering helpers (single source of truth;
// this file previously carried its own copy of the resampler).
#include "filter_utils.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// File-private helpers
// ============================================================================

namespace {

    // ---------- string utilities ----------

    std::vector<std::string> parse_csv_row(const std::string& line) {
        std::vector<std::string> fields;
        std::string cur;
        for (char c : line) {
            if (c == ',') { fields.push_back(cur); cur.clear(); }
            else cur += c;
        }
        fields.push_back(cur);
        for (auto& f : fields) {
            size_t first = f.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) f = "";
            else {
                size_t last = f.find_last_not_of(" \t\r\n");
                f = f.substr(first, last - first + 1);
            }
        }
        return fields;
    }

    bool contains(std::string haystack, std::string needle) {
        if (needle.empty()) return false;
        std::transform(haystack.begin(), haystack.end(), haystack.begin(), ::toupper);
        std::transform(needle.begin(), needle.end(), needle.begin(), ::toupper);
        return haystack.find(needle) != std::string::npos;
    }

    // Case-insensitive equality on trimmed strings. Use this when channel names
    // share a common prefix (NLS_NOM_PRESS_BLD_ART vs ..._ART_PULM vs ..._ART_ABP).
    bool equals_ci(std::string a, std::string b) {
        std::transform(a.begin(), a.end(), a.begin(), ::toupper);
        std::transform(b.begin(), b.end(), b.begin(), ::toupper);
        return a == b;
    }

    // ---------- .dat column reading ----------

    std::vector<std::string> find_real_header(std::istream& in) {
        std::string line;
        while (std::getline(in, line)) {
            if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
            bool looksLikeHeader =
                contains(line, "NLS_NOM_") ||
                contains(line, "NLS_EEG_") ||
                (contains(line, "Index") && contains(line, "TimeStamp"));
            if (looksLikeHeader) return parse_csv_row(line);
        }
        return {};
    }

    // infer_row_rate was removed; the .dat row rate now comes from
    // cfg.ecg_raw_rate directly. One source of truth.


    // One-pass column store for a CHAOS .dat file. Reading the same file 80+
    // times to extract 40 channels is the main bottleneck in make_binfile_dat.
    // This struct caches every column's populated values + their row indices,
    // keyed by header-column index. Downstream readers do O(1) lookups.
    struct PrescannedDat {
        std::vector<std::string> headers;
        size_t totalRows = 0;
        std::vector<std::vector<size_t>> rowIdxPerCol;
        std::vector<std::vector<double>> valuePerCol;
    };

    PrescannedDat prescan_dat_columns(const std::filesystem::path& path) {
        PrescannedDat out;
        std::ifstream in(path);
        if (!in) return out;

        out.headers = find_real_header(in);
        if (out.headers.empty()) return out;

        const size_t nCols = out.headers.size();
        out.rowIdxPerCol.assign(nCols, {});
        out.valuePerCol.assign(nCols, {});
        for (size_t c = 0; c < nCols; ++c) {
            out.rowIdxPerCol[c].reserve(1 << 16);
            out.valuePerCol[c].reserve(1 << 16);
        }

        std::string line;
        size_t rowIdx = 0;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::vector<std::string> row = parse_csv_row(line);
            const size_t lim = std::min(row.size(), nCols);
            for (size_t c = 0; c < lim; ++c) {
                if (row[c].empty()) continue;
                try {
                    double v = std::stod(row[c]);
                    out.rowIdxPerCol[c].push_back(rowIdx);
                    out.valuePerCol[c].push_back(v);
                }
                catch (...) { /* ignore malformed cell */ }
            }
            ++rowIdx;
        }
        out.totalRows = rowIdx;
        return out;
    }

    int find_column_index(const PrescannedDat& dat,
        const std::string& label, bool exactMatch)
    {
        if (label.empty()) return -1;
        for (int i = 0; i < (int)dat.headers.size(); ++i) {
            bool match = exactMatch ? equals_ci(dat.headers[i], label)
                : contains(dat.headers[i], label);
            if (match) return i;
        }
        return -1;
    }

    void column_raw_with_indices(const PrescannedDat& dat,
        const std::string& label, bool exactMatch,
        std::vector<double>& outValues,
        std::vector<size_t>& outRowIndices)
    {
        outValues.clear();
        outRowIndices.clear();
        int colIdx = find_column_index(dat, label, exactMatch);
        if (colIdx < 0) return;
        outValues = dat.valuePerCol[colIdx];
        outRowIndices = dat.rowIdxPerCol[colIdx];
    }

    // ---------- output writers ----------

    // Write a placeholder for a missing channel: a single -1.0 for the upsampled
    // slot and a single (-1.0, -1.0) sentinel pair for the raw slot. Native rate
    // and upsample rate = 0.
    void write_missing(std::ofstream& out,
        uint32_t& sizeUpOut, uint32_t& sizeRawOut,
        float& nativeRateOut, float& upRateOut)
    {
        double v = -1.0;
        out.write(reinterpret_cast<const char*>(&v), 8);
        sizeUpOut = 1;

        double sentinel[2] = { -1.0, -1.0 };
        out.write(reinterpret_cast<const char*>(sentinel), 16);
        sizeRawOut = 1;

        nativeRateOut = 0.0f;
        upRateOut = 0.0f;
    }

    // Write an EDF channel as a (upsampled, raw-pairs) pair. `finalSamplingRate`
    // is this channel's own target rate from config; it is stamped into the
    // per-channel upsample-rate header block via upRateOut.
    // Read a channel over the sample window [winStartSamp, winEndSamp) on this
    // channel's OWN native grid (or the whole channel when winEndSamp < 0), then
    // upsample and write the (upsampled, raw-pairs) blocks. The raw block's
    // epoch axis is anchored at winStartEpochMs so a segment's absolute time is
    // correct regardless of which 8-hour window it is.
    void edf_to_bin(int handle, int idx, long long n,
        double old_rate, double finalSamplingRate, double winStartEpochMs,
        std::ofstream& out,
        uint32_t& sizeUpOut, uint32_t& sizeRawOut,
        float& nativeRateOut, float& upRateOut,
        long long winStartSamp = 0, long long winEndSamp = -1)
    {
        if (idx < 0 || n <= 0 || finalSamplingRate <= 0.0) {
            write_missing(out, sizeUpOut, sizeRawOut, nativeRateOut, upRateOut);
            return;
        }

        // Clamp the window to the channel's available samples.
        long long lo = (winStartSamp < 0) ? 0 : winStartSamp;
        long long hi = (winEndSamp < 0) ? n : std::min(winEndSamp, n);
        if (lo >= hi) {                       // nothing in this window
            write_missing(out, sizeUpOut, sizeRawOut, nativeRateOut, upRateOut);
            return;
        }
        long long count = hi - lo;

        // Position this signal's read cursor at the window start, then read only
        // the window. edfseek is per-signal; EDFSEEK_SET == 0 in edflib.
        edfseek(handle, idx, lo, EDFSEEK_SET);
        std::vector<double> raw(count);
        edfread_physical_samples(handle, idx, (int)count, raw.data());

        std::vector<double> up = filterutils::upsample(raw, old_rate, finalSamplingRate);

        out.write(reinterpret_cast<const char*>(up.data()), up.size() * 8);
        sizeUpOut = (uint32_t)up.size();

        const double dtMs = (old_rate > 0.0) ? (1000.0 / old_rate) : 0.0;
        for (size_t k = 0; k < raw.size(); ++k) {
            double pair[2] = { winStartEpochMs + static_cast<double>(k) * dtMs, raw[k] };
            out.write(reinterpret_cast<const char*>(pair), 16);
        }
        sizeRawOut = (uint32_t)raw.size();
        nativeRateOut = static_cast<float>(old_rate);
        upRateOut = static_cast<float>(finalSamplingRate);
    }

    void write_synthetic_timestamp(std::ofstream& out,
        double durationSec, double nativeRate,
        double finalSamplingRate, double startEpochMs,
        uint32_t& sizeUpOut, uint32_t& sizeRawOut,
        float& nativeRateOut, float& upRateOut)
    {
        if (durationSec <= 0.0 || nativeRate <= 0.0 || finalSamplingRate <= 0.0) {
            write_missing(out, sizeUpOut, sizeRawOut, nativeRateOut, upRateOut);
            return;
        }

        // Upsampled block: absolute Unix-epoch milliseconds on the ECG target grid.
        const size_t upLen = (size_t)std::ceil(durationSec * finalSamplingRate);
        std::vector<double> up(upLen);
        const double dtUpMs = 1000.0 / finalSamplingRate;
        for (size_t k = 0; k < upLen; ++k) up[k] = startEpochMs + (double)k * dtUpMs;
        out.write(reinterpret_cast<const char*>(up.data()), up.size() * 8);
        sizeUpOut = (uint32_t)up.size();

        // Raw block: same epoch-ms axis at the native grid.
        const size_t rawLen = (size_t)std::floor(durationSec * nativeRate) + 1;
        const double dtRawMs = 1000.0 / nativeRate;
        for (size_t k = 0; k < rawLen; ++k) {
            const double t = startEpochMs + (double)k * dtRawMs;
            double pair[2] = { t, t };
            out.write(reinterpret_cast<const char*>(pair), 16);
        }
        sizeRawOut = (uint32_t)rawLen;
        nativeRateOut = static_cast<float>(nativeRate);
        upRateOut = static_cast<float>(finalSamplingRate);
    }

    void write_header_and_close(std::ofstream& out,
        const config_entry& cfg,
        const uint32_t sizes_up[NUM_CHANNELS],
        const uint32_t sizes_raw[NUM_CHANNELS],
        const float    native_rates[NUM_CHANNELS],
        const float    up_rates[NUM_CHANNELS],
        uint32_t sleep_size)
    {
        out.seekp(0);
        uint32_t header_version = BIN_HEADER_VERSION;   // offset 0
        uint32_t n_channels = (uint32_t)NUM_CHANNELS; // offset 4
        out.write(reinterpret_cast<const char*>(&header_version), 4);
        out.write(reinterpret_cast<const char*>(&n_channels), 4);
        uint32_t sleep_state_len = (uint32_t)cfg.sleepstate_length;
        out.write(reinterpret_cast<const char*>(&sleep_state_len), 4);
        out.write(reinterpret_cast<const char*>(sizes_up), NUM_CHANNELS * 4);
        out.write(reinterpret_cast<const char*>(sizes_raw), NUM_CHANNELS * 4);
        out.write(reinterpret_cast<const char*>(native_rates), NUM_CHANNELS * 4);
        out.write(reinterpret_cast<const char*>(up_rates), NUM_CHANNELS * 4);
        out.write(reinterpret_cast<const char*>(&sleep_size), 4);
        out.close();
    }

    // ---------- EDF helpers ----------

    // EDF channel map: edf_signal_idx[ChannelIdx] = signal index in the EDF file
    // (or -1 if absent). All slot labels come from config_entry -- if a field
    // is empty for this dataset, that slot stays -1 and write_missing is used.
    using EdfSignalMap = std::array<int, NUM_CHANNELS>;

    EdfSignalMap build_edf_channel_map(const edf_hdr_struct* hdr,
        const config_entry& cfg)
    {
        EdfSignalMap m;
        m.fill(-1);
        std::set<int> used;

        auto find = [&](const std::string& label) -> int {
            if (label.empty()) return -1;
            // Exact (case-insensitive) match first: channel names in this family
            // share prefixes (NLS_NOM_PRESS_BLD_ART vs ..._ART_PULM vs ..._ART_ABP),
            // so a substring test would let the shorter name steal the longer
            // column. equals_ci pins it to the exact column.
            for (int i = 0; i < hdr->edfsignals; ++i) {
                if (!used.count(i) && equals_ci(hdr->signalparam[i].label, label)) {
                    used.insert(i);
                    return i;
                }
            }
            // Fall back to substring for datasets whose labels aren't exact
            // (e.g. trailing units/whitespace in EDF labels).
            for (int i = 0; i < hdr->edfsignals; ++i) {
                if (!used.count(i) && contains(hdr->signalparam[i].label, label)) {
                    used.insert(i);
                    return i;
                }
            }
            return -1;
            };

        m[CH_ECG1] = find(cfg.ecg_1_label);
        m[CH_ECG2] = find(cfg.ecg_2_label);
        m[CH_ECG3] = find(cfg.ecg_3_label);
        m[CH_PPG] = find(cfg.ppg_label);
        m[CH_ACCEL_X] = find(cfg.accel_x_label);
        m[CH_ACCEL_Y] = find(cfg.accel_y_label);
        m[CH_ACCEL_Z] = find(cfg.accel_z_label);
        m[CH_EEG1] = find(cfg.eeg_1_label);
        m[CH_EEG2] = find(cfg.eeg_2_label);
        m[CH_EEG3] = find(cfg.eeg_3_label);
        m[CH_EEG4] = find(cfg.eeg_4_label);
        m[CH_CVP] = find(cfg.cvp_label);
        m[CH_RESP] = find(cfg.resp_label);
        m[CH_MARKER] = find(cfg.marker_label);
        m[CH_PACEMAKER_EVENT] = find(cfg.pacemaker_label);
        m[CH_TEMP] = find(cfg.temp_label);
        m[CH_ABP] = find(cfg.abp_label);
        m[CH_ART] = find(cfg.art_label);
        m[CH_ART_PULM] = find(cfg.art_pulm_label);
        m[CH_PRES] = find(cfg.pres_label);
        m[CH_FLOW] = find(cfg.flow_label);
        m[CH_SNORE] = find(cfg.snore_label);
        m[CH_THOR] = find(cfg.thor_label);
        m[CH_ABDO] = find(cfg.abdo_label);
        m[CH_LEG] = find(cfg.leg_label);
        m[CH_OXSTATUS] = find(cfg.oxstatus_label);
        m[CH_SPO2] = find(cfg.spo2_label);
        m[CH_HR] = find(cfg.hr_label);
        m[CH_DHR] = find(cfg.dhr_label);
        m[CH_EOG_L] = find(cfg.eog_l_label);
        m[CH_EOG_R] = find(cfg.eog_r_label);
        m[CH_EMG] = find(cfg.emg_label);
        m[CH_AUXAC] = find(cfg.auxac_label);
        m[CH_THERM] = find(cfg.therm_label);
        m[CH_POS] = find(cfg.pos_label);
        return m;
    }

    void repair_edf_startdate(const std::filesystem::path& path) {
        /* SHHS EDFs exported by Compumedics wrote the start date as mm.dd.yy, so
        "05.29.01" (29 May 2001) fails edflib's dd.mm.yy validation with
        "startdate of recording is invalid: expected dd.mm.yy where mm should be
        more than 00 and less than 13". edfopen_file_readonly refuses the file, so
        nothing downstream ever runs.*/
        std::fstream f(path, std::ios::in | std::ios::out | std::ios::binary);
        if (!f) return;   // unreadable or read-only; edfopen will say so

        char raw[8] = {};
        f.seekg(168);
        f.read(raw, 8);
        if (!f) return;

        int a = 0, b = 0, c = 0;
        if (std::sscanf(std::string(raw, 8).c_str(), "%d.%d.%d", &a, &b, &c) != 3)
            return;

        const bool ddmm = (a >= 1 && a <= 31) && (b >= 1 && b <= 12);
        const bool mmdd = (b >= 1 && b <= 31) && (a >= 1 && a <= 12);

        if (ddmm) return;    // already conforming (covers the ambiguous case too)
        if (!mmdd) return;   // not a date either way; let edflib report it

        char fixed[16];
        std::snprintf(fixed, sizeof(fixed), "%02d.%02d.%02d", b, a, c);
        f.seekp(168);
        f.write(fixed, 8);
        f.flush();

        std::cerr << "  [edf] repaired mm.dd.yy start date in "
            << path.filename().string() << ": "
            << std::string(raw, 8) << " -> " << std::string(fixed, 8) << "\n";
    }


    double edf_channel_rate(const edf_hdr_struct* hdr, int idx) {
        if (idx < 0) return 0.0;
        return (double)hdr->signalparam[idx].smp_in_datarecord /
            ((double)hdr->datarecord_duration / 10000000.0);
    }

    long long edf_samples(const edf_hdr_struct* hdr, int idx) {
        return (idx < 0) ? 0 : hdr->signalparam[idx].smp_in_file;
    }

    std::filesystem::path make_out_path(const std::filesystem::path& src,
        const config_entry& cfg)
    {
        // Channels no longer share one upsample rate, so the file name no
        // longer carries a rate suffix -- just <stem>.bin.
        return std::filesystem::path(cfg.output_path) /
            (src.stem().string() + ".bin");
    }


    // Find a sleep-stage XML next to the source file when sleepExt is set.
    std::filesystem::path findSleepXml(const std::filesystem::path& src,
        const std::string& sleepExt) {
        if (sleepExt.empty()) return {};
        std::string want = sleepExt;
        std::transform(want.begin(), want.end(), want.begin(), ::toupper);
        std::string stem = src.stem().string();

        std::filesystem::path profusion, other;
        for (const auto& f :
            std::filesystem::directory_iterator(src.parent_path())) {
            std::string e = f.path().extension().string();
            std::transform(e.begin(), e.end(), e.begin(), ::toupper);
            if (e != want) continue;
            if (f.path().stem().string().find(stem) == std::string::npos) continue;

            if (contains(f.path().stem().string(), "profusion")) {
                profusion = f.path();           // the one the parser can read
            }
            else if (!contains(f.path().stem().string(), "nsrr")) {
                // Neither Profusion nor NSRR: a dataset with one unqualified XML
                // (MESA), which is Profusion-schema anyway. Kept as the fallback
                // so this stays a no-op for every existing dataset.
                if (other.empty()) other = f.path();
            }
        }

        if (!profusion.empty()) return profusion;
        if (!other.empty())     return other;

        // Only an -nsrr.xml exists. Say so rather than hand back a file whose
        // schema the caller cannot parse and let it look like "no stages".
        std::cerr << "WARNING: only an NSRR-schema XML found for "
            << src.filename().string()
            << "; its <ScoredEvent> stages are not readable by the "
            "//SleepStage parser, so no hypnogram will be written\n";
        return {};
    }

    // Parse the wall-clock timestamp column of a CHAOS .dat into a per-row
    // vector of "seconds since the first row's timestamp." Returns a vector
    // of size totalRows; rows that fail to parse get a sentinel of -1.0
    // (caller can then fall back to the synthetic row_idx/row_rate timeline
    // for those rows).
    //
    // Prefers "System TimeStamp UTC" (wall-clock, has the gaps) over
    // "Monitor TimeStamp" (device-relative). std::stod can't parse
    // "YYYYMMDD HH:MM:SS.mmm" so the prescan's per-column double vectors
    // are useless for this; we do a dedicated text pass here.
    std::vector<double> parse_dat_timestamps(const std::filesystem::path& path,
        size_t totalRows)
    {
        std::vector<double> times(totalRows, -1.0);
        if (totalRows == 0) return times;

        std::ifstream in(path);
        if (!in) return times;

        std::vector<std::string> hdrs = find_real_header(in);
        if (hdrs.empty()) return times;

        int tsCol = -1;
        for (int i = 0; i < (int)hdrs.size(); ++i) {
            if (contains(hdrs[i], "System TimeStamp UTC")) { tsCol = i; break; }
        }
        if (tsCol < 0) {
            for (int i = 0; i < (int)hdrs.size(); ++i) {
                if (contains(hdrs[i], "Monitor TimeStamp")) { tsCol = i; break; }
            }
        }
        if (tsCol < 0) return times;

        auto daysFromCivil = [](int y, int m, int d) -> long long {
            y -= (m <= 2);
            const int era = (y >= 0 ? y : y - 399) / 400;
            const unsigned yoe = static_cast<unsigned>(y - era * 400);
            const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
            const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
            return (long long)era * 146097 + (long long)doe - 719468;
            };

        auto parseAbsSec = [&](const std::string& s, double& outSec) -> bool {
            int y = 0, mo = 0, d = 0, hh = 0, mm = 0, ss = 0, ms = 0;
            if (std::sscanf(s.c_str(), "%4d%2d%2d %d:%d:%d.%d",
                &y, &mo, &d, &hh, &mm, &ss, &ms) >= 6 ||
                std::sscanf(s.c_str(), "%4d%2d%2d %d:%d:%d",
                    &y, &mo, &d, &hh, &mm, &ss) == 6)
            {
                const long long days = daysFromCivil(y, mo, d);
                outSec = (double)days * 86400.0
                    + hh * 3600.0 + mm * 60.0 + ss + ms / 1000.0;
                return true;
            }
            return false;
            };

        double firstAbs = 0.0;
        bool   firstSet = false;
        std::string line;
        size_t rowIdx = 0;
        while (std::getline(in, line) && rowIdx < totalRows) {
            if (line.empty()) { ++rowIdx; continue; }
            std::vector<std::string> cells = parse_csv_row(line);
            if (tsCol < (int)cells.size() && !cells[tsCol].empty()) {
                double abs = 0.0;
                if (parseAbsSec(cells[tsCol], abs)) {
                    if (!firstSet) { firstAbs = abs; firstSet = true; }
                    times[rowIdx] = abs - firstAbs;
                }
            }
            ++rowIdx;
        }

        // Backfill any unparsed rows by linear interpolation between the
        // nearest valid neighbors. Edge unparsed runs (start / end) carry
        // the nearest valid value.
        if (firstSet) {
            int n = (int)times.size();
            int firstValid = -1;
            for (int i = 0; i < n; ++i) if (times[i] >= 0.0) { firstValid = i; break; }
            if (firstValid < 0) return times;
            for (int i = 0; i < firstValid; ++i) times[i] = times[firstValid];

            int i = firstValid;
            while (i < n) {
                int j = i + 1;
                while (j < n && times[j] < 0.0) ++j;
                if (j >= n) {
                    for (int k = i + 1; k < n; ++k) times[k] = times[i];
                    break;
                }
                if (j > i + 1) {
                    const double t0 = times[i], t1 = times[j];
                    const double step = (t1 - t0) / (j - i);
                    for (int k = i + 1; k < j; ++k) times[k] = t0 + step * (k - i);
                }
                i = j;
            }
        }

        return times;
    }

    // Correct "double-rate then blank" anomalies in a sparse .dat channel.
    //
    // Some recordings contain runs where a channel was sampled at ~2x its
    // configured rate (e.g. PPG populated every 2 rows instead of every 4 on a
    // 500 Hz grid), each run immediately followed by a blank of equal length
    // where that channel has no samples while ECG / the row grid continue. The
    // per-sample timestamps in these files are unreliable, so detection is
    // purely on row-index stride: a run whose spacing is below the channel's
    // normal stride is the 2x region. Re-spacing those samples back onto the
    // normal stride spreads them across the trailing blank and closes it, so
    // both the upsampled and raw blocks (which share the row-index clock) see
    // one continuous stream.
    //
    // Operates on row indices only (values untouched). Handles many runs per
    // channel. No-op for grid-rate channels (ECG, normalStride <= 1).
    //
    // Two guards prevent a tail-shift: a run is only re-spaced if it has a
    // trailing anchor (a normal sample after it -- otherwise there is no blank
    // to fill and nothing to bound the expansion, so it would run past the real
    // end of data), and if it is at least kMinRun samples long (a one- or
    // two-sample short delta is sensor jitter, not a 2x region, and re-spacing
    // it would permanently offset every later sample). A genuine gap -- normal
    // stride before it, no 2x compression -- never triggers and stays open.
    void fix_double_rate_regions(std::vector<size_t>& rawRowIdx,
        double rowRate, double nativeHz)
    {
        if (rawRowIdx.size() < 2 || rowRate <= 0.0 || nativeHz <= 0.0) return;

        const double normalStride = rowRate / nativeHz;        // PPG: 500/125 = 4
        if (normalStride <= 1.0) return;                       // ECG: grid rate
        const double compressedThresh = normalStride * 0.75;   // delta below => 2x
        constexpr size_t kMinRun = 4;                          // ignore short blips

        const size_t n = rawRowIdx.size();
        size_t i = 0;
        while (i + 1 < n) {
            size_t j = i;
            while (j + 1 < n &&
                (double)(rawRowIdx[j + 1] - rawRowIdx[j]) < compressedThresh) {
                ++j;
            }
            const size_t runLen = j - i + 1;
            if (j > i && runLen >= kMinRun && j + 1 < n) {
                const size_t anchor = rawRowIdx[i];
                const size_t nextRow = rawRowIdx[j + 1];
                for (size_t m = 1; m <= j - i; ++m) {
                    size_t newIdx = anchor +
                        (size_t)std::llround((double)m * normalStride);
                    if (newIdx <= rawRowIdx[i + m - 1]) newIdx = rawRowIdx[i + m - 1] + 1;
                    if (newIdx >= nextRow)              newIdx = nextRow - 1;
                    rawRowIdx[i + m] = newIdx;
                }
                i = j;
            }
            else {
                ++i;
            }
        }
    }
    // Build a dense NATIVE-rate signal from a blank-padded .dat column. Each real
    // sample is placed on the native grid at slot = round(rowIdx / stride). Over-
    // dense (2x-packed) samples that would collide are pushed to the next slot
    // (spread); native slots with no real sample stay NaN (absent -> a gap, never
    // interpolated). The caller fills these gaps before polyphase upsampling.
    std::vector<double> depad_column_to_native(
        const std::vector<double>& vals, const std::vector<size_t>& rowIdx,
        size_t totalRows, double rowRate, double nativeHz)
    {
        if (vals.empty() || rowRate <= 0.0 || nativeHz <= 0.0) return {};
        const double stride = rowRate / nativeHz;                 // rows per native sample
        const size_t nOut = (size_t)std::ceil((double)totalRows / stride);
        if (nOut == 0) return {};
        std::vector<double> out(nOut, std::numeric_limits<double>::quiet_NaN());
        size_t lastSlot = (size_t)-1;
        for (size_t k = 0; k < vals.size(); ++k) {
            size_t slot = (size_t)std::llround((double)rowIdx[k] / stride);
            if (lastSlot != (size_t)-1 && slot <= lastSlot) slot = lastSlot + 1;  // spread 2x
            if (slot >= nOut) break;
            out[slot] = vals[k];
            lastSlot = slot;
        }
        return out;
    }

    // ---------- absolute time anchoring (Unix epoch milliseconds) ----------

    // Days since 1970-01-01 for a civil date (proleptic Gregorian).
    long long days_from_civil(int y, int m, int d) {
        y -= (m <= 2);
        const int era = (y >= 0 ? y : y - 399) / 400;
        const unsigned yoe = (unsigned)(y - era * 400);
        const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
        const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
        return (long long)era * 146097 + (long long)doe - 719468;
    }

    // Recording start as Unix-epoch MILLISECONDS from the EDF header.
    // NOTE: EDF stores wall-clock with no timezone; we treat it as UTC. If the
    // EDFs are in local time the epoch is offset by the TZ. subsecond is
    // edflib's units of 100 ns (0..9,999,999).
    double edf_start_epoch_ms(const edf_hdr_struct* hdr) {
        if (hdr->startdate_year < 1970) {
            std::cerr << "WARNING: EDF start year " << hdr->startdate_year
                << " looks invalid; epoch anchor set to 0\n";
            return 0.0;
        }
        const long long days = days_from_civil(hdr->startdate_year,
            hdr->startdate_month,
            hdr->startdate_day);
        const double sec = (double)days * 86400.0
            + hdr->starttime_hour * 3600.0
            + hdr->starttime_minute * 60.0
            + hdr->starttime_second
            + (double)hdr->starttime_subsecond / 1e7;
        return sec * 1000.0;
    }

    // First parseable "System TimeStamp UTC" (fallback "Monitor TimeStamp") of a
    // .dat as Unix-epoch ms. Returns 0.0 (and warns) if none is found -- the file
    // then carries an elapsed-from-zero axis rather than an absolute one.
    double dat_start_epoch_ms(const std::filesystem::path& path) {
        std::ifstream in(path);
        if (!in) return 0.0;
        std::vector<std::string> hdrs = find_real_header(in);
        if (hdrs.empty()) return 0.0;

        int tsCol = -1;
        for (int i = 0; i < (int)hdrs.size(); ++i)
            if (contains(hdrs[i], "System TimeStamp UTC")) { tsCol = i; break; }
        if (tsCol < 0)
            for (int i = 0; i < (int)hdrs.size(); ++i)
                if (contains(hdrs[i], "Monitor TimeStamp")) { tsCol = i; break; }
        if (tsCol < 0) {
            std::cerr << "WARNING: no UTC timestamp column in "
                << path.filename().string() << "; epoch anchor set to 0\n";
            return 0.0;
        }

        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::vector<std::string> cells = parse_csv_row(line);
            if (tsCol >= (int)cells.size() || cells[tsCol].empty()) continue;
            int d = 0, mo = 0, y = 0, hh = 0, mm = 0, ss = 0, ms = 0;
            // CHAOS "System TimeStamp UTC" comes in two forms across files:
            //   dashed  "DD-MM-YYYY HH:MM:SS.mmm"  (day first)
            //   packed  "YYYYMMDD HH:MM:SS.mmm"    (year first, no separators)
            // Detect by the '-': dashed => d,mo,y; packed => %4d%2d%2d => y,mo,d.
            const std::string& ts = cells[tsCol];
            bool ok;
            if (ts.find('-') != std::string::npos)
                ok = std::sscanf(ts.c_str(), "%d-%d-%d %d:%d:%d.%d",
                    &d, &mo, &y, &hh, &mm, &ss, &ms) >= 6;
            else
                ok = std::sscanf(ts.c_str(), "%4d%2d%2d %d:%d:%d.%d",
                    &y, &mo, &d, &hh, &mm, &ss, &ms) >= 6;
            if (ok) {
                const long long days = days_from_civil(y, mo, d);
                const double sec = (double)days * 86400.0 + hh * 3600.0 + mm * 60.0 + ss + ms / 1000.0;
                return sec * 1000.0;
            }
        }
        std::cerr << "WARNING: could not parse any UTC timestamp in "
            << path.filename().string() << "; epoch anchor set to 0\n";
        return 0.0;
    }

}   // anonymous namespace

// ============================================================================
// Public entry points
// ============================================================================

// Full-recording conversion: window covers the entire file.
void make_binfile_edf(const std::filesystem::path& path, const config_entry& cfg)
{
    make_binfile_edf_window(path, cfg, make_out_path(path, cfg), 0.0, -1.0);
}

// Windowed conversion. [winStartSec, winEndSec) bounds the output on the
// recording's own seconds-from-start clock; winEndSec < 0 means "to the end".
// outPath lets a caller name the segment (e.g. <stem>_0.bin). Every channel is
// sliced to the same wall-clock window on its own native grid, so the two time
// axes (native raw / upsampled) stay identical across the segment set.
void make_binfile_edf_window(const std::filesystem::path& path, const config_entry& cfg, const std::filesystem::path& outPath,  double winStartSec, double winEndSec)
{
    char filebuf[1 << 16];
    std::ofstream out;
    out.rdbuf()->pubsetbuf(filebuf, sizeof(filebuf));
    out.open(outPath, std::ios::binary | std::ios::trunc);
    if (!out.is_open()) {
        std::cerr << "ERROR: cannot create " << outPath << "\n";
        return;
    }

    std::vector<char> zeroes(HEADER_SIZE, 0);
    out.write(zeroes.data(), HEADER_SIZE);

    auto hdr = std::make_unique<edf_hdr_struct>();
    repair_edf_startdate(path);
    if (edfopen_file_readonly(path.string().c_str(), hdr.get(),
        EDFLIB_READ_ALL_ANNOTATIONS)) {
        std::cerr << "ERROR: cannot open EDF " << path
            << " (edflib code " << hdr->filetype << ")\n";
        out.close();
        std::filesystem::remove(outPath);
        return;
    }

    EdfSignalMap sigmap = build_edf_channel_map(hdr.get(), cfg);

    // The recording's absolute start; a windowed segment's own start is offset
    // by winStartSec so its timestamps remain correct wall-clock epoch ms.
    const double recStartEpochMs = edf_start_epoch_ms(hdr.get());
    const double segStartEpochMs = recStartEpochMs
        + ((winStartSec > 0.0) ? winStartSec * 1000.0 : 0.0);

    uint32_t sizes_up[NUM_CHANNELS] = {};
    uint32_t sizes_raw[NUM_CHANNELS] = {};
    float    native_rates[NUM_CHANNELS] = {};
    float    up_rates[NUM_CHANNELS] = {};

    // Convert a [winStartSec, winEndSec) span into a [lo, hi) sample window on a
    // channel of the given native rate. winEndSec < 0 => to end (hi = -1).
    auto sampleWindow = [&](double rate, long long& lo, long long& hi) {
        lo = (winStartSec > 0.0) ? (long long)std::floor(winStartSec * rate) : 0;
        hi = (winEndSec >= 0.0) ? (long long)std::floor(winEndSec * rate) : -1;
        };

    // Channel 0: synthetic timestamp (seconds from start). EDF has no monitor
    // clock, so anchor the raw block to the primary ECG channel's rate and the
    // upsampled block to the ECG target rate. For a window, the duration is the
    // window length (clamped to what the ECG channel actually has).
    {
        double tsRate = (cfg.ecg_raw_rate > 0.0)
            ? cfg.ecg_raw_rate
            : edf_channel_rate(hdr.get(), sigmap[CH_ECG1]);
        double tsDur = 0.0;
        if (sigmap[CH_ECG1] >= 0 && tsRate > 0.0) {
            long long total = edf_samples(hdr.get(), sigmap[CH_ECG1]);
            long long lo, hi; sampleWindow(tsRate, lo, hi);
            if (lo < 0) lo = 0;
            if (hi < 0 || hi > total) hi = total;
            long long cnt = (hi > lo) ? (hi - lo) : 0;
            tsDur = (double)cnt / tsRate;
        }
        write_synthetic_timestamp(out, tsDur, tsRate, cfg.ecg_upsample_rate, segStartEpochMs,
            sizes_up[CH_TIMESTAMP], sizes_raw[CH_TIMESTAMP],
            native_rates[CH_TIMESTAMP], up_rates[CH_TIMESTAMP]);
    }

    // Each channel is upsampled to its own configured target rate. Unmapped
    // slots (-1 in sigmap) or channels with no configured rate become
    // missing-channel placeholders via edf_to_bin -> write_missing. The window
    // is recomputed per channel because each has its own native rate.
    auto write_signal_to_bin = [&](ChannelIdx ch, double rawRate = 0.0,
        double upRate = 0.0) {
            int chIdx = sigmap[ch];
            long long n = (chIdx < 0) ? 0 : edf_samples(hdr.get(), chIdx);
            long long lo = 0, hi = -1;
            if (chIdx >= 0 && rawRate > 0.0) sampleWindow(rawRate, lo, hi);
            edf_to_bin(hdr->handle, chIdx, n, rawRate, upRate, segStartEpochMs, out,
                sizes_up[ch], sizes_raw[ch], native_rates[ch], up_rates[ch],
                lo, hi);
        };

    write_signal_to_bin(CH_ECG1, cfg.ecg_raw_rate, cfg.ecg_upsample_rate);
    write_signal_to_bin(CH_ECG2, cfg.ecg_raw_rate, cfg.ecg_upsample_rate);
    write_signal_to_bin(CH_ECG3, cfg.ecg_raw_rate, cfg.ecg_upsample_rate);
    write_signal_to_bin(CH_PPG, cfg.ppg_raw_rate, cfg.ppg_upsample_rate);

    write_signal_to_bin(CH_ACCEL_X, cfg.accel_raw_rate, cfg.accel_upsample_rate);
    write_signal_to_bin(CH_ACCEL_Y, cfg.accel_raw_rate, cfg.accel_upsample_rate);
    write_signal_to_bin(CH_ACCEL_Z, cfg.accel_raw_rate, cfg.accel_upsample_rate);
    write_signal_to_bin(CH_MARKER, cfg.marker_raw_rate, cfg.marker_upsample_rate);
    write_signal_to_bin(CH_TEMP, cfg.temp_raw_rate, cfg.temp_upsample_rate);
    write_signal_to_bin(CH_PACEMAKER_EVENT, cfg.pacemaker_raw_rate, cfg.pacemaker_upsample_rate);

    write_signal_to_bin(CH_EOG_L, cfg.eog_l_raw_rate, cfg.eog_l_upsample_rate);
    write_signal_to_bin(CH_EOG_R, cfg.eog_r_raw_rate, cfg.eog_r_upsample_rate);
    write_signal_to_bin(CH_EMG, cfg.emg_raw_rate, cfg.emg_upsample_rate);
    write_signal_to_bin(CH_EEG1, cfg.eeg_raw_rate, cfg.eeg_upsample_rate);
    write_signal_to_bin(CH_EEG2, cfg.eeg_raw_rate, cfg.eeg_upsample_rate);
    write_signal_to_bin(CH_EEG3, cfg.eeg_raw_rate, cfg.eeg_upsample_rate);
    write_signal_to_bin(CH_EEG4, cfg.eeg_raw_rate, cfg.eeg_upsample_rate);

    write_signal_to_bin(CH_CVP, cfg.cvp_raw_rate, cfg.cvp_upsample_rate);
    write_signal_to_bin(CH_PRES, cfg.pres_raw_rate, cfg.pres_upsample_rate);
    write_signal_to_bin(CH_FLOW, cfg.flow_raw_rate, cfg.flow_upsample_rate);
    write_signal_to_bin(CH_SNORE, cfg.snore_raw_rate, cfg.snore_upsample_rate);
    write_signal_to_bin(CH_THOR, cfg.thor_raw_rate, cfg.thor_upsample_rate);
    write_signal_to_bin(CH_ABDO, cfg.abdo_raw_rate, cfg.abdo_upsample_rate);
    write_signal_to_bin(CH_LEG, cfg.leg_raw_rate, cfg.leg_upsample_rate);
    write_signal_to_bin(CH_AUXAC, cfg.auxac_raw_rate, cfg.auxac_upsample_rate);
    write_signal_to_bin(CH_THERM, cfg.therm_raw_rate, cfg.therm_upsample_rate);
    write_signal_to_bin(CH_POS, cfg.pos_raw_rate, cfg.pos_upsample_rate);

    write_signal_to_bin(CH_OXSTATUS, cfg.oxstatus_raw_rate, cfg.oxstatus_upsample_rate);
    write_signal_to_bin(CH_SPO2, cfg.spo2_raw_rate, cfg.spo2_upsample_rate);
    write_signal_to_bin(CH_HR, cfg.hr_raw_rate, cfg.hr_upsample_rate);
    write_signal_to_bin(CH_DHR, cfg.dhr_raw_rate, cfg.dhr_upsample_rate);

    write_signal_to_bin(CH_RESP, cfg.resp_raw_rate, cfg.resp_upsample_rate);
    write_signal_to_bin(CH_ABP, cfg.abp_raw_rate, cfg.abp_upsample_rate);
    write_signal_to_bin(CH_ART, cfg.art_raw_rate, cfg.art_upsample_rate);
    write_signal_to_bin(CH_ART_PULM, cfg.art_pulm_raw_rate, cfg.art_pulm_upsample_rate);

    edfclose_file(hdr->handle);

    auto sleep_path = findSleepXml(path, cfg.sleep_file_extention);

    std::vector<double> stages;
    if (!cfg.sleep_file_extention.empty() && !sleep_path.empty()
        && std::filesystem::exists(sleep_path)) {
        pugi::xml_document doc;
        if (doc.load_file(sleep_path.string().c_str())) {
            for (auto node : doc.select_nodes("//SleepStage")) {
                double v = node.node().text().as_double();
                stages.push_back(v == 5.0 ? 4.0 : v);
            }
        }
    }
    if (stages.empty()) stages.push_back(-1.0);
    uint32_t sleep_size = (uint32_t)stages.size();
    out.write(reinterpret_cast<const char*>(stages.data()),
        sleep_size * sizeof(double));

    write_header_and_close(out, cfg, sizes_up, sizes_raw, native_rates, up_rates, sleep_size);
}

void make_binfile_dat(const std::filesystem::path& path,
    const config_entry& cfg)
{
    std::filesystem::path outPath = make_out_path(path, cfg);

    char filebuf[1 << 16];
    std::ofstream out;
    out.rdbuf()->pubsetbuf(filebuf, sizeof(filebuf));
    out.open(outPath, std::ios::binary | std::ios::trunc);
    if (!out.is_open()) {
        std::cerr << "ERROR: cannot create " << outPath << "\n";
        return;
    }

    std::vector<char> zeroes(HEADER_SIZE, 0);
    out.write(zeroes.data(), HEADER_SIZE);

    uint32_t sizes_up[NUM_CHANNELS] = {};
    uint32_t sizes_raw[NUM_CHANNELS] = {};
    float    native_rates[NUM_CHANNELS] = {};
    float    up_rates[NUM_CHANNELS] = {};

    // The .dat has one row per ECG sample, so the row rate IS the ECG
    // native rate from the config. No inference -- config is the single
    // source of truth.
    const double row_rate = cfg.ecg_raw_rate;
    if (row_rate <= 0.0) {
        std::cerr << "ERROR: cfg.ecg_raw_rate is 0; aborting conversion of "
            << path.filename().string() << "\n";
        return;
    }
    PrescannedDat prescan = prescan_dat_columns(path);

    const double startEpochMs = dat_start_epoch_ms(path);



    // Write one column from the .dat. Empty label, column-not-found,
    // nativeHz<=0, or upHz<=0 all produce a missing-channel placeholder.
    //
    // `exact` controls header matching: substring (false) for most channels,
    // exact (true) when names share a common prefix (e.g. ART vs ART_PULM
    // vs ART_ABP) and substring would alias.
    // `nativeHz` is stamped into the .bin's native-rates header; `upHz` is
    // this channel's target rate, stamped into the upsample-rates header.
    auto writeCol = [&](const std::string& label, bool exact, ChannelIdx ch,
        double nativeHz, double upHz) {
            if (label.empty() || nativeHz <= 0.0 || upHz <= 0.0) {
                write_missing(out, sizes_up[ch], sizes_raw[ch], native_rates[ch], up_rates[ch]);
                return;
            }

            std::vector<double> rawValues;
            std::vector<size_t> rawRowIdx;
            column_raw_with_indices(prescan, label, exact, rawValues, rawRowIdx);

            if (rawValues.empty()) {
                write_missing(out, sizes_up[ch], sizes_raw[ch], native_rates[ch], up_rates[ch]);
                return;
            }

            // De-pad the blank-padded column onto its native grid (config rate;
            // 2x-packed runs spread; absent stretches left as NaN gaps), then
            // polyphase native -> target -- the same upsample() the EDF path uses.
            fix_double_rate_regions(rawRowIdx, row_rate, nativeHz);
            std::vector<double> dense =
                depad_column_to_native(rawValues, rawRowIdx, prescan.totalRows, row_rate, nativeHz);
            // Fill absent (NaN) slots by linear hold so the polyphase input is
            // continuous: the upsampled LINE runs straight through gaps (OpenGL
            // can't break on NaN, which is fine -- the raw scatter, empty over
            // gaps, is the real "data present" cue). No NaN reaches the block.
            for (size_t nn = 0; nn < dense.size(); ++nn) {
                if (!std::isnan(dense[nn])) continue;
                size_t p = nn; while (p > 0 && std::isnan(dense[p - 1])) --p;
                size_t q = nn; while (q + 1 < dense.size() && std::isnan(dense[q + 1])) ++q;
                const double a = (p > 0) ? dense[p - 1] : 0.0;
                const double b = (q + 1 < dense.size()) ? dense[q + 1] : a;
                for (size_t m = p; m <= q; ++m) {
                    const double f = (q + 1 > p) ? (double)(m - p + 1) / (double)(q - p + 2) : 0.0;
                    dense[m] = a + (b - a) * f;
                }
                nn = q;
            }
            std::vector<double> up = filterutils::upsample(dense, nativeHz, upHz);

            if (up.empty()) {
                double v = -1.0;
                out.write(reinterpret_cast<const char*>(&v), 8);
                sizes_up[ch] = 1;
            }
            else {
                out.write(reinterpret_cast<const char*>(up.data()), up.size() * 8);
                sizes_up[ch] = (uint32_t)up.size();
            }

            // Raw block x = corrected row-index time (row index / grid rate) --
            // the SAME clock resample_from_sparse uses for the upsampled block.
            // The .dat's own timestamps are unreliable, and fix_double_rate_regions
            // has already corrected rawRowIdx for 2x packing while leaving genuine
            // gaps open, so row index is the single source of truth: the raw
            // scatter now lands exactly on the upsampled trace everywhere.
            const double dtMs = (row_rate > 0.0) ? (1000.0 / row_rate) : 0.0;
            for (size_t k = 0; k < rawValues.size(); ++k) {
                const double t = startEpochMs + (double)rawRowIdx[k] * dtMs;
                double pair[2] = { t, rawValues[k] };
                out.write(reinterpret_cast<const char*>(pair), 16);
            }
            sizes_raw[ch] = (uint32_t)rawValues.size();
            native_rates[ch] = (float)nativeHz;
            up_rates[ch] = (float)upHz;
        };

    auto writeMissing = [&](ChannelIdx ch) {
        write_missing(out, sizes_up[ch], sizes_raw[ch], native_rates[ch], up_rates[ch]);
        };

    // Channel 0: synthetic timestamp. Upsampled at the ECG target rate.
    {
        const double durationSec = (row_rate > 0.0)
            ? (double)prescan.totalRows / row_rate : 0.0;
        write_synthetic_timestamp(out, durationSec, row_rate, cfg.ecg_upsample_rate, startEpochMs,
            sizes_up[CH_TIMESTAMP], sizes_raw[CH_TIMESTAMP],
            native_rates[CH_TIMESTAMP], up_rates[CH_TIMESTAMP]);
    }

    // Algorithm-facing channels.
    writeCol(cfg.ecg_1_label, false, CH_ECG1, cfg.ecg_raw_rate, cfg.ecg_upsample_rate);
    writeCol(cfg.ecg_2_label, false, CH_ECG2, cfg.ecg_raw_rate, cfg.ecg_upsample_rate);
    writeCol(cfg.ecg_3_label, false, CH_ECG3, cfg.ecg_raw_rate, cfg.ecg_upsample_rate);
    writeCol(cfg.ppg_label, false, CH_PPG, cfg.ppg_raw_rate, cfg.ppg_upsample_rate);

    // Accel / marker / temp / pacemaker / EOG / EMG aren't in CHAOS .dat
    // files. Emit placeholders so the 35-slot layout stays consistent.
    for (int ch = CH_ACCEL_X; ch <= CH_EMG; ++ch)
        writeMissing((ChannelIdx)ch);

    // EEG columns exist in config.csv now (populated for MESA), but CHAOS
    // .dat files don't carry EEG, so these stay placeholders here.
    writeMissing(CH_EEG1);
    writeMissing(CH_EEG2);
    writeMissing(CH_EEG3);
    writeMissing(CH_EEG4);

    // Central venous pressure. `pres` is a separate slot filled only for
    // MESA (an EDF dataset), so it's a placeholder on the CHAOS .dat path.
    writeCol(cfg.cvp_label, false, CH_CVP, cfg.cvp_raw_rate, cfg.cvp_upsample_rate);
    writeMissing(CH_PRES);

    // Sleep apnea slots (flow, snore, thor, abdo, leg, therm, pos) and the
    // SpO2 family (oxstatus, spo2, HR, DHR) aren't in CHAOS .dat.
    for (int ch = CH_FLOW; ch <= CH_DHR; ++ch)
        writeMissing((ChannelIdx)ch);

    // Respiration + arterial pressures.
    //
    // ABP / ART / ART_PULM use exact-match because their labels share the
    // prefix NLS_NOM_PRESS_BLD_ART and substring matching would alias them.
    // ART / ART_PULM are still emitted as placeholders here -- see the note
    // in the summary if you want them populated from the new art_* config
    // columns.
    writeCol(cfg.resp_label, false, CH_RESP, cfg.resp_raw_rate, cfg.resp_upsample_rate);
    writeCol(cfg.abp_label, true, CH_ABP, cfg.abp_raw_rate, cfg.abp_upsample_rate);
    writeCol(cfg.art_label, true, CH_ART, cfg.art_raw_rate, cfg.art_upsample_rate);
    writeCol(cfg.art_pulm_label, true, CH_ART_PULM, cfg.art_pulm_raw_rate, cfg.art_pulm_upsample_rate);

    double placeholder = -1.0;
    out.write(reinterpret_cast<const char*>(&placeholder), sizeof(double));
    uint32_t sleep_size = 1;

    write_header_and_close(out, cfg, sizes_up, sizes_raw, native_rates, up_rates, sleep_size);
}

namespace {
    // Recording duration in seconds, taken from the primary ECG channel's
    // sample count and configured native rate. Returns 0 on any failure.
    double edf_duration_seconds(const std::filesystem::path& path,
        const config_entry& cfg)
    {
        auto hdr = std::make_unique<edf_hdr_struct>();
        if (edfopen_file_readonly(path.string().c_str(), hdr.get(),
            EDFLIB_READ_ALL_ANNOTATIONS)) {
            return 0.0;
        }
        EdfSignalMap sigmap = build_edf_channel_map(hdr.get(), cfg);
        double rate = (cfg.ecg_raw_rate > 0.0)
            ? cfg.ecg_raw_rate
            : edf_channel_rate(hdr.get(), sigmap[CH_ECG1]);
        long long n = (sigmap[CH_ECG1] >= 0)
            ? edf_samples(hdr.get(), sigmap[CH_ECG1]) : 0;
        edfclose_file(hdr->handle);
        return (rate > 0.0 && n > 0) ? (double)n / rate : 0.0;
    }

    // Seconds in one output segment. 8 hours.
    constexpr double SEGMENT_SECONDS = 8.0 * 3600.0;
}

std::filesystem::path make_binfile(const std::filesystem::path& path, const config_entry& cfg)
{
    /*
        creates output path, and calls the edf or dat specific function to write the bin file.
    */
    std::filesystem::path out = make_out_path(path, cfg);

    std::string ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);

    if (ext == ".EDF") {
        // BITTIUM: split the recording into fixed 8-hour segments, each written
        // as its own valid .bin named <stem>_0.bin, <stem>_1.bin, ... The final
        // segment shorter than 8 h is kept as-is. All other EDF datasets write
        // a single whole-recording .bin.
        if (cfg.dataset_type == "BITTIUM") {
            double dur = edf_duration_seconds(path, cfg);
            if (dur <= 0.0) {
                std::cerr << "ERROR: could not determine duration of "
                    << path.filename().string() << "; writing single .bin\n";
                make_binfile_edf(path, cfg);
                return out;
            }
            int nSeg = (int)std::ceil(dur / SEGMENT_SECONDS);
            if (nSeg < 1) nSeg = 1;
            const std::string stem = path.stem().string();
            std::filesystem::path firstOut;
            for (int s = 0; s < nSeg; ++s) {
                double t0 = s * SEGMENT_SECONDS;
                double t1 = (s == nSeg - 1) ? -1.0 : (s + 1) * SEGMENT_SECONDS;
                std::filesystem::path segOut =
                    std::filesystem::path(cfg.output_path) /
                    (stem + "_" + std::to_string(s) + ".bin");
                if (s == 0) firstOut = segOut;
                make_binfile_edf_window(path, cfg, segOut, t0, t1);
            }
            return firstOut;
        }
        make_binfile_edf(path, cfg);
    }
    else if (ext == ".DAT") {
        make_binfile_dat(path, cfg);
    }
    else {
        std::cerr << "ERROR: unsupported file type " << ext
            << " for " << path << "\n";
        return {};
    }
    return out;
}