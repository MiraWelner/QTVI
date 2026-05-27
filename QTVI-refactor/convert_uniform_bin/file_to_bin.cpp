/**
 * @file   file_to_bin.cpp
 * @brief  Implementation of the source-to-bin conversion declared in
 *         file_to_bin.hpp. EDF and .dat paths share the same on-disk
 *         channel layout but pull samples from very different shapes of
 *         input: EDF channels are uniform native-rate streams; .dat
 *         columns are sparse populated cells inside a uniform monitor
 *         row grid. Both end up resampled to cfg.finalSamplingRate for
 *         the upsampled block, with a (t, v) raw block preserving the
 *         channel's native sampling.
 *
 *         All dataset-specific channel names live in config_entry and
 *         are populated by applyDefaultChannelLabels(). This file knows
 *         nothing about MESA / Bittium / CHAOS conventions -- if a
 *         label field is empty, that slot writes a missing-channel
 *         placeholder.
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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// File-private helpers
// ============================================================================

namespace {

    // ---------- polyphase resampler ----------
    //
    // The resampler builds a filter bank of small filters (one per
    // fractional phase between input samples), then for each output sample
    // dot-products the right sub-filter against a window of inputs. Output is
    // split into a leading boundary (filter hangs off the left edge), a multi-
    // threaded interior (no bounds checks), and a trailing boundary.

    int greatest_common_divisor(int a, int b) {
        a = std::abs(a); b = std::abs(b);
        while (b) { int t = b; b = a % b; a = t; }
        return a;
    }

    std::vector<std::vector<double>> buildPolyphaseBank(int P, int Q, int halfLobes) {
        int maxPQ = std::max(P, Q);
        int numTaps = 2 * halfLobes * maxPQ + 1;
        double fc = 1.0 / static_cast<double>(maxPQ);
        int M = numTaps - 1;
        double halfM = M / 2.0;

        std::vector<double> h(numTaps);
        for (int n = 0; n < numTaps; ++n) {
            double x = n - halfM;
            double sinc = (std::abs(x) < 1e-12) ? 1.0
                : std::sin(M_PI * fc * x) / (M_PI * x);
            double w = 0.42 - 0.5 * std::cos(2.0 * M_PI * n / M)
                + 0.08 * std::cos(4.0 * M_PI * n / M);
            h[n] = sinc * w;
        }

        int subLen = (numTaps + P - 1) / P;
        std::vector<std::vector<double>> bank(P, std::vector<double>(subLen, 0.0));
        for (int i = 0; i < numTaps; ++i) {
            bank[i % P][i / P] = h[i];
        }

        for (int p = 0; p < P; ++p) {
            double s = 0.0;
            for (int k = 0; k < subLen; ++k) s += bank[p][k];
            if (std::abs(s) > 1e-15) {
                for (int k = 0; k < subLen; ++k) bank[p][k] /= s;
            }
        }
        return bank;
    }

    void processRangeBoundary(
        const double* inPtr, long long inLen,
        double* outPtr, long long mStart, long long mEnd,
        const std::vector<const double*>& bankPtrs,
        int subLen, int filterCenter, int P, int Q)
    {
        for (long long m = mStart; m < mEnd; ++m) {
            long long upsampledIdx = m * Q;
            int phase = static_cast<int>(upsampledIdx % P);
            long long baseInput = upsampledIdx / P;
            const double* sub = bankPtrs[phase];

            double sum = 0.0;
            for (int k = 0; k < subLen; ++k) {
                long long inIdx = baseInput - k + filterCenter;
                if (inIdx >= 0 && inIdx < inLen) {
                    sum += sub[k] * inPtr[inIdx];
                }
            }
            outPtr[m] = sum;
        }
    }

    void processRangeInterior(
        const double* inPtr,
        double* outPtr, long long mStart, long long mEnd,
        const std::vector<const double*>& bankPtrs,
        int subLen, int filterCenter, int P, int Q)
    {
        for (long long m = mStart; m < mEnd; ++m) {
            long long upsampledIdx = m * Q;
            int phase = static_cast<int>(upsampledIdx % P);
            long long baseInput = upsampledIdx / P;
            const double* sub = bankPtrs[phase];
            const double* inBase = inPtr + baseInput + filterCenter;

            double sum = 0.0;
            for (int k = 0; k < subLen; ++k) {
                sum += sub[k] * inBase[-k];
            }
            outPtr[m] = sum;
        }
    }

    std::vector<double> polyphase_resample(const std::vector<double>& input, int P, int Q) {
        if (input.empty()) return {};
        if (P == 1 && Q == 1) return input;

        int halfLobes = std::max(16, std::max(P, Q) / 2);
        auto bank = buildPolyphaseBank(P, Q, halfLobes);
        int subLen = static_cast<int>(bank[0].size());

        long long inLen = static_cast<long long>(input.size());
        long long outLen = static_cast<long long>(
            std::ceil(static_cast<double>(inLen) * P / Q));
        std::vector<double> output(outLen);

        int maxPQ = std::max(P, Q);
        int filterCenter = halfLobes * maxPQ / P;

        const double* inPtr = input.data();
        double* outPtr = output.data();

        std::vector<const double*> bankPtrs(P);
        for (int p = 0; p < P; ++p) bankPtrs[p] = bank[p].data();

        long long safeStartM = 0;
        long long safeEndM = 0;
        for (long long m = 0; m < outLen; ++m) {
            long long baseInput = (m * Q) / P;
            if (baseInput - subLen + 1 + filterCenter >= 0 &&
                baseInput + filterCenter < inLen) {
                safeStartM = m;
                break;
            }
        }
        for (long long m = outLen - 1; m >= safeStartM; --m) {
            long long baseInput = (m * Q) / P;
            if (baseInput - subLen + 1 + filterCenter >= 0 &&
                baseInput + filterCenter < inLen) {
                safeEndM = m + 1;
                break;
            }
        }

        processRangeBoundary(inPtr, inLen, outPtr, 0, safeStartM,
            bankPtrs, subLen, filterCenter, P, Q);
        processRangeBoundary(inPtr, inLen, outPtr, safeEndM, outLen,
            bankPtrs, subLen, filterCenter, P, Q);

        long long interiorLen = safeEndM - safeStartM;
        if (interiorLen <= 0) return output;

        unsigned int numThreads = std::thread::hardware_concurrency();
        if (numThreads == 0) numThreads = 4;
        if (interiorLen < 100000) numThreads = 1;
        numThreads = std::min(numThreads, static_cast<unsigned int>(interiorLen));

        if (numThreads == 1) {
            processRangeInterior(inPtr, outPtr, safeStartM, safeEndM,
                bankPtrs, subLen, filterCenter, P, Q);
        }
        else {
            std::vector<std::thread> threads;
            threads.reserve(numThreads);
            long long chunkSize = interiorLen / numThreads;

            for (unsigned int t = 0; t < numThreads; ++t) {
                long long start = safeStartM + t * chunkSize;
                long long end = (t == numThreads - 1) ? safeEndM : start + chunkSize;

                threads.emplace_back(processRangeInterior,
                    inPtr, outPtr, start, end,
                    std::cref(bankPtrs), subLen, filterCenter, P, Q);
            }
            for (auto& th : threads) th.join();
        }
        return output;
    }

    std::vector<double> upsample(const std::vector<double>& input,
        double sourceRate, double targetRate) {
        if (input.empty()) return {};
        if (sourceRate == targetRate) return input;

        int gcd = greatest_common_divisor((int)targetRate, (int)sourceRate);
        int P = (int)targetRate / gcd;
        int Q = (int)sourceRate / gcd;

        if (P > 1000 || Q > 1000) {
            throw std::runtime_error(
                "Resampling ratio " + std::to_string(P) + "/" +
                std::to_string(Q) + " too large - source rate: " +
                std::to_string(sourceRate));
        }
        return polyphase_resample(input, P, Q);
    }

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
    // cfg.ecgRate directly. One source of truth.


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
    // slot and a single (-1.0, -1.0) sentinel pair for the raw slot. Native rate = 0.
    void write_missing(std::ofstream& out,
        uint32_t& sizeUpOut, uint32_t& sizeRawOut,
        float& nativeRateOut)
    {
        double v = -1.0;
        out.write(reinterpret_cast<const char*>(&v), 8);
        sizeUpOut = 1;

        double sentinel[2] = { -1.0, -1.0 };
        out.write(reinterpret_cast<const char*>(sentinel), 16);
        sizeRawOut = 1;

        nativeRateOut = 0.0f;
    }

    // Write an EDF channel as a (upsampled, raw-pairs) pair.
    void edf_to_bin(int handle, int idx, long long n,
        double old_rate, double finalSamplingRate,
        std::ofstream& out,
        uint32_t& sizeUpOut, uint32_t& sizeRawOut,
        float& nativeRateOut)
    {
        if (idx < 0 || n <= 0) {
            write_missing(out, sizeUpOut, sizeRawOut, nativeRateOut);
            return;
        }

        std::vector<double> raw(n);
        edfread_physical_samples(handle, idx, (int)n, raw.data());

        const bool skip_resample = (old_rate <= BOOLEAN_RATE);
        if (skip_resample) {
            out.write(reinterpret_cast<const char*>(raw.data()), raw.size() * 8);
            sizeUpOut = (uint32_t)raw.size();
        }
        else {
            std::vector<double> up = upsample(raw, old_rate, finalSamplingRate);
            out.write(reinterpret_cast<const char*>(up.data()), up.size() * 8);
            sizeUpOut = (uint32_t)up.size();
        }

        const double dt = (old_rate > 0.0) ? (1.0 / old_rate) : 0.0;
        for (size_t k = 0; k < raw.size(); ++k) {
            double pair[2] = { static_cast<double>(k) * dt, raw[k] };
            out.write(reinterpret_cast<const char*>(pair), 16);
        }
        sizeRawOut = (uint32_t)raw.size();
        nativeRateOut = static_cast<float>(old_rate);
    }

    void write_synthetic_timestamp(std::ofstream& out,
        double durationSec, double nativeRate,
        double finalSamplingRate,
        uint32_t& sizeUpOut, uint32_t& sizeRawOut,
        float& nativeRateOut)
    {
        if (durationSec <= 0.0 || nativeRate <= 0.0) {
            write_missing(out, sizeUpOut, sizeRawOut, nativeRateOut);
            return;
        }

        const size_t upLen = (size_t)std::ceil(durationSec * finalSamplingRate);
        std::vector<double> up(upLen);
        const double dtUp = 1.0 / finalSamplingRate;
        for (size_t k = 0; k < upLen; ++k) up[k] = (double)k * dtUp;
        out.write(reinterpret_cast<const char*>(up.data()), up.size() * 8);
        sizeUpOut = (uint32_t)up.size();

        const size_t rawLen = (size_t)std::floor(durationSec * nativeRate) + 1;
        const double dtRaw = 1.0 / nativeRate;
        for (size_t k = 0; k < rawLen; ++k) {
            const double t = (double)k * dtRaw;
            double pair[2] = { t, t };
            out.write(reinterpret_cast<const char*>(pair), 16);
        }
        sizeRawOut = (uint32_t)rawLen;
        nativeRateOut = static_cast<float>(nativeRate);
    }

    void write_header_and_close(std::ofstream& out,
        double finalSamplingRate,
        const uint32_t sizes_up[NUM_CHANNELS],
        const uint32_t sizes_raw[NUM_CHANNELS],
        const float    native_rates[NUM_CHANNELS],
        uint32_t sleep_size)
    {
        out.seekp(0);
        uint32_t scalars[4] = {
            (uint32_t)finalSamplingRate,
            (uint32_t)BOOLEAN_RATE,
            PACEMAKER_RATE,
            (uint32_t)SLEEP_STATE_LENGTH
        };
        out.write(reinterpret_cast<const char*>(scalars), sizeof(scalars));
        out.write(reinterpret_cast<const char*>(sizes_up), NUM_CHANNELS * 4);
        out.write(reinterpret_cast<const char*>(sizes_raw), NUM_CHANNELS * 4);
        out.write(reinterpret_cast<const char*>(native_rates), NUM_CHANNELS * 4);
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
            for (int i = 0; i < hdr->edfsignals; ++i) {
                if (!used.count(i) && contains(hdr->signalparam[i].label, label)) {
                    used.insert(i);
                    return i;
                }
            }
            return -1;
            };

        // Algorithm-facing channels.
        m[CH_ECG1] = find(cfg.ecg1Label);
        m[CH_ECG2] = find(cfg.ecg2Label);
        m[CH_ECG3] = find(cfg.ecg3Label);
        m[CH_PPG] = find(cfg.ppgLabel);

        // Accelerometer.
        m[CH_ACCEL_X] = find(cfg.accelXLabel);
        m[CH_ACCEL_Y] = find(cfg.accelYLabel);
        m[CH_ACCEL_Z] = find(cfg.accelZLabel);

        // EEG.
        m[CH_EEG1] = find(cfg.eeg1Label);
        m[CH_EEG2] = find(cfg.eeg2Label);
        m[CH_EEG3] = find(cfg.eeg3Label);
        m[CH_EEG4] = find(cfg.eeg4Label);

        // Pressures / respiratory.
        m[CH_PRES] = find(cfg.cvpLabel);
        m[CH_RESP] = find(cfg.respLabel);
        m[CH_ABP] = find(cfg.abpLabel);
        m[CH_ART] = find(cfg.artLabel);
        m[CH_ART_PULM] = find(cfg.artPulmLabel);

        // Slots not yet promoted to config fields stay -1 and get
        // write_missing automatically:
        //   CH_MARKER, CH_TEMP, CH_PACEMAKER,
        //   CH_EOG_L, CH_EOG_R, CH_EMG,
        //   CH_FLOW, CH_THOR, CH_ABDO, CH_LEG, CH_THERM, CH_POS,
        //   CH_EKG_OFF, CH_EOG_L_OFF, CH_EOG_R_OFF, CH_EMG_OFF,
        //   CH_EEG1_OFF, CH_EEG2_OFF, CH_EEG3_OFF,
        //   CH_OXSTATUS, CH_SPO2, CH_HR, CH_DHR

        return m;
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
        return std::filesystem::path(cfg.output_path) /
            (src.stem().string() + "_" +
                std::to_string((int)cfg.finalSamplingRate) + ".bin");
    }


    // Find a sleep-stage XML next to the source file when sleepExt is set.
    std::filesystem::path findSleepXml(const std::filesystem::path& src,
        const std::string& sleepExt) {
        if (sleepExt.empty()) return {};
        std::string want = sleepExt;
        std::transform(want.begin(), want.end(), want.begin(), ::toupper);
        std::string stem = src.stem().string();
        for (const auto& f :
            std::filesystem::directory_iterator(src.parent_path())) {
            std::string e = f.path().extension().string();
            std::transform(e.begin(), e.end(), e.begin(), ::toupper);
            if (e == want && f.path().stem().string().find(stem) != std::string::npos)
                return f.path();
        }
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

}   // anonymous namespace

// ============================================================================
// Public entry points
// ============================================================================

void make_binfile_edf(const std::filesystem::path& path, const config_entry& cfg)
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

    auto hdr = std::make_unique<edf_hdr_struct>();
    if (edfopen_file_readonly(path.string().c_str(), hdr.get(),
        EDFLIB_READ_ALL_ANNOTATIONS)) {
        std::cerr << "ERROR: cannot open EDF " << path << "\n";
        out.close();
        std::filesystem::remove(outPath);
        return;
    }

    EdfSignalMap sigmap = build_edf_channel_map(hdr.get(), cfg);

    uint32_t sizes_up[NUM_CHANNELS] = {};
    uint32_t sizes_raw[NUM_CHANNELS] = {};
    float    native_rates[NUM_CHANNELS] = {};

    // Channel 0: synthetic timestamp (seconds from start). EDF has no monitor
    // clock, so anchor the raw block to the primary ECG channel's rate.
    {
        double tsRate = (cfg.ecgRate > 0.0)
            ? cfg.ecgRate
            : edf_channel_rate(hdr.get(), sigmap[CH_ECG1]);
        double tsDur = 0.0;
        if (sigmap[CH_ECG1] >= 0 && tsRate > 0.0) {
            tsDur = (double)edf_samples(hdr.get(), sigmap[CH_ECG1]) / tsRate;
        }
        write_synthetic_timestamp(out, tsDur, tsRate, cfg.finalSamplingRate,
            sizes_up[CH_TIMESTAMP], sizes_raw[CH_TIMESTAMP],
            native_rates[CH_TIMESTAMP]);
    }

    auto writeChannel = [&](ChannelIdx ch, double rateOverride = 0.0) {
        int chIdx = sigmap[ch];
        double rate;
        if (chIdx < 0) {
            rate = 0.0;   // doesn't matter; edf_to_bin will write_missing
        }
        else if (rateOverride > 0.0) {
            rate = rateOverride;
        }
        else {
            rate = edf_channel_rate(hdr.get(), chIdx);
        }
        long long n = (chIdx < 0) ? 0 : edf_samples(hdr.get(), chIdx);
        edf_to_bin(hdr->handle, chIdx, n, rate, cfg.finalSamplingRate, out,
            sizes_up[ch], sizes_raw[ch], native_rates[ch]);
        };

    // Channels in fixed slot order. Unmapped slots (-1 in sigmap) become
    // missing-channel placeholders via edf_to_bin -> write_missing.
    writeChannel(CH_ECG1, cfg.ecgRate);
    writeChannel(CH_ECG2, cfg.ecgRate);
    writeChannel(CH_ECG3, cfg.ecgRate);
    writeChannel(CH_PPG, cfg.ppgRate);

    writeChannel(CH_ACCEL_X);
    writeChannel(CH_ACCEL_Y);
    writeChannel(CH_ACCEL_Z);

    writeChannel(CH_MARKER);
    writeChannel(CH_TEMP);
    writeChannel(CH_PACEMAKER);

    writeChannel(CH_EOG_L);
    writeChannel(CH_EOG_R);
    writeChannel(CH_EMG);
    writeChannel(CH_EEG1);
    writeChannel(CH_EEG2);
    writeChannel(CH_EEG3);
    writeChannel(CH_EEG4);

    writeChannel(CH_PRES);
    writeChannel(CH_FLOW);
    writeChannel(CH_THOR);
    writeChannel(CH_ABDO);
    writeChannel(CH_LEG);
    writeChannel(CH_THERM);
    writeChannel(CH_POS);

    writeChannel(CH_EKG_OFF);
    writeChannel(CH_EOG_L_OFF);
    writeChannel(CH_EOG_R_OFF);
    writeChannel(CH_EMG_OFF);
    writeChannel(CH_EEG1_OFF);
    writeChannel(CH_EEG2_OFF);
    writeChannel(CH_EEG3_OFF);

    writeChannel(CH_OXSTATUS);
    writeChannel(CH_SPO2);
    writeChannel(CH_HR);
    writeChannel(CH_DHR);

    writeChannel(CH_RESP);
    writeChannel(CH_ABP);
    writeChannel(CH_ART);
    writeChannel(CH_ART_PULM);

    edfclose_file(hdr->handle);

    auto sleep_path = findSleepXml(path, cfg.sleepExt);

    std::vector<double> stages;
    if (!cfg.sleepExt.empty() && !sleep_path.empty()
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

    write_header_and_close(out, cfg.finalSamplingRate,
        sizes_up, sizes_raw, native_rates, sleep_size);
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

    // The .dat has one row per ECG sample, so the row rate IS the ECG
    // native rate from the config. No inference -- config is the single
    // source of truth.
    const double row_rate = cfg.ecgRate;
    if (row_rate <= 0.0) {
        std::cerr << "ERROR: cfg.ecgRate is 0; aborting conversion of "
            << path.filename().string() << "\n";
        return;
    }
    PrescannedDat prescan = prescan_dat_columns(path);

    // Real per-row wall-clock timestamps (seconds since first row's
    // System TimeStamp UTC). Used by writeCol() to put each raw (t, v)
    // pair at its true time -- otherwise gaps in the source recording
    // get silently filled by the row_idx/row_rate timeline below.
    std::vector<double> rowTimestamps = parse_dat_timestamps(path, prescan.totalRows);
    const bool haveRealTimestamps = !rowTimestamps.empty()
        && std::any_of(rowTimestamps.begin(), rowTimestamps.end(),
            [](double t) { return t > 0.0; });

    // Resample a channel to finalSamplingRate using ONLY its real populated
    // samples, placed at their true times. Linear interpolation between
    // consecutive real samples. Avoids the "fake plateau" bug from treating
    // a sparse column (populated every Nth row) as a dense signal.
    const double finalSR = cfg.finalSamplingRate;
    auto resample_from_sparse = [&](const std::vector<double>& rawValues,
        const std::vector<size_t>& rawRowIdx,
        size_t totalRows) -> std::vector<double>
        {
            if (rawValues.size() < 2 || row_rate <= 0.0) return {};

            const double totalDur = (double)totalRows / row_rate;
            const size_t outLen = (size_t)std::ceil(totalDur * finalSR);
            if (outLen == 0) return {};

            std::vector<double> result(outLen);
            const double inv_row_rate = 1.0 / row_rate;
            const size_t N = rawValues.size();

            size_t k = 0;
            for (size_t m = 0; m < outLen; ++m) {
                const double t = (double)m / finalSR;
                while (k + 1 < N &&
                    (double)rawRowIdx[k + 1] * inv_row_rate <= t) {
                    ++k;
                }
                if (k + 1 >= N) {
                    result[m] = rawValues[N - 1];
                    continue;
                }
                const double t0 = (double)rawRowIdx[k] * inv_row_rate;
                const double t1 = (double)rawRowIdx[k + 1] * inv_row_rate;
                if (t <= t0) { result[m] = rawValues[k]; continue; }
                const double span = t1 - t0;
                const double f = (span > 0.0) ? (t - t0) / span : 0.0;
                result[m] = rawValues[k] * (1.0 - f) + rawValues[k + 1] * f;
            }
            return result;
        };

    // Write one column from the .dat. Empty label, column-not-found, or
    // nativeHz<=0 all produce a missing-channel placeholder.
    //
    // `exact` controls header matching: substring (false) for most channels,
    // exact (true) when names share a common prefix (e.g. ART vs ART_PULM
    // vs ART_ABP) and substring would alias.
    // `nativeHz` is stamped directly into the .bin's native-rates header.
    auto writeCol = [&](const std::string& label, bool exact, ChannelIdx ch,
        double nativeHz) {
            if (label.empty() || nativeHz <= 0.0) {
                write_missing(out, sizes_up[ch], sizes_raw[ch], native_rates[ch]);
                return;
            }

            std::vector<double> rawValues;
            std::vector<size_t> rawRowIdx;
            column_raw_with_indices(prescan, label, exact, rawValues, rawRowIdx);

            if (rawValues.empty()) {
                write_missing(out, sizes_up[ch], sizes_raw[ch], native_rates[ch]);
                return;
            }

            std::vector<double> up =
                resample_from_sparse(rawValues, rawRowIdx, prescan.totalRows);

            if (up.empty()) {
                double v = -1.0;
                out.write(reinterpret_cast<const char*>(&v), 8);
                sizes_up[ch] = 1;
            }
            else {
                out.write(reinterpret_cast<const char*>(up.data()), up.size() * 8);
                sizes_up[ch] = (uint32_t)up.size();
            }

            // Raw block: real wall-clock timestamps when we have them, so
            // gaps in the source recording survive into the .bin and the
            // gap_indicator can find them. Falls back to row_idx*dt only for
            // rows whose timestamp didn't parse.
            const double dt = (row_rate > 0.0) ? (1.0 / row_rate) : 0.0;
            for (size_t k = 0; k < rawValues.size(); ++k) {
                const size_t rIdx = rawRowIdx[k];
                double t;
                if (haveRealTimestamps && rIdx < rowTimestamps.size()
                    && rowTimestamps[rIdx] >= 0.0) {
                    t = rowTimestamps[rIdx];
                }
                else {
                    t = (double)rIdx * dt;
                }
                double pair[2] = { t, rawValues[k] };
                out.write(reinterpret_cast<const char*>(pair), 16);
            }
            sizes_raw[ch] = (uint32_t)rawValues.size();
            native_rates[ch] = (float)nativeHz;
        };

    auto writeMissing = [&](ChannelIdx ch) {
        write_missing(out, sizes_up[ch], sizes_raw[ch], native_rates[ch]);
        };

    // Channel 0: synthetic timestamp.
    {
        const double durationSec = (row_rate > 0.0)
            ? (double)prescan.totalRows / row_rate : 0.0;
        write_synthetic_timestamp(out, durationSec, row_rate, finalSR,
            sizes_up[CH_TIMESTAMP], sizes_raw[CH_TIMESTAMP],
            native_rates[CH_TIMESTAMP]);
    }

    // Algorithm-facing channels.
    writeCol(cfg.ecg1Label, false, CH_ECG1, cfg.ecgRate);
    writeCol(cfg.ecg2Label, false, CH_ECG2, cfg.ecgRate);
    writeCol(cfg.ecg3Label, false, CH_ECG3, cfg.ecgRate);
    writeCol(cfg.ppgLabel, false, CH_PPG, cfg.ppgRate);

    // Accel / marker / temp / pacemaker / EOG / EMG aren't in CHAOS .dat
    // files. Emit placeholders so the 40-slot layout stays consistent.
    for (int ch = CH_ACCEL_X; ch <= CH_EMG; ++ch)
        writeMissing((ChannelIdx)ch);

    // EEG: no rate in config.csv, so passes 0 -> writeCol emits a missing
    // placeholder. Add an eeg_rate column to config.csv + cfg.eegRate to
    // config_entry if you ever need to populate these.
    writeCol(cfg.eeg1Label, false, CH_EEG1, 0.0);
    writeCol(cfg.eeg2Label, false, CH_EEG2, 0.0);
    writeCol(cfg.eeg3Label, false, CH_EEG3, 0.0);
    writeCol(cfg.eeg4Label, false, CH_EEG4, 0.0);

    // Central venous pressure.
    writeCol(cfg.cvpLabel, false, CH_PRES, cfg.cvpRate);

    // Sleep apnea slots and SpO2 family aren't in CHAOS .dat.
    for (int ch = CH_FLOW; ch <= CH_DHR; ++ch)
        writeMissing((ChannelIdx)ch);

    // Respiration + arterial pressures.
    //
    // ABP / ART / ART_PULM use exact-match because their labels share the
    // prefix NLS_NOM_PRESS_BLD_ART and substring matching would alias them.
    // ART and ART_PULM use abpRate -- they're variants on the same
    // arterial-pressure sensor and share its native rate.
    writeCol(cfg.respLabel, false, CH_RESP, cfg.respRate);
    writeCol(cfg.abpLabel, true, CH_ABP, cfg.abpRate);
    writeCol(cfg.artLabel, true, CH_ART, cfg.abpRate);
    writeCol(cfg.artPulmLabel, true, CH_ART_PULM, cfg.abpRate);

    double placeholder = -1.0;
    out.write(reinterpret_cast<const char*>(&placeholder), sizeof(double));
    uint32_t sleep_size = 1;

    write_header_and_close(out, finalSR,
        sizes_up, sizes_raw, native_rates, sleep_size);
}

std::filesystem::path make_binfile(const std::filesystem::path& path, const config_entry& cfg)
{
    /*
        creates output path, and calls the edf or dat specific function to write the bin file.
    */
    std::filesystem::path out = std::filesystem::path(cfg.output_path) /
        (path.stem().string() + "_" + std::to_string((int)cfg.finalSamplingRate) + ".bin");

    std::string ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);

    if (ext == ".EDF") {
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