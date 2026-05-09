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
 *         Anything in here that another .cpp would call lives in the
 *         header. Everything else is in an anonymous namespace.
 */

#include "file_to_bin.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

extern "C" {
#include "edflib.h"
}
#include "pugixml.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// File-private helpers
// ============================================================================

namespace {

    // ---------- polyphase resampler ----------
    //
    // Used to be in resample.hpp; folded in here because file_to_bin is the only
    // consumer. The resampler builds a filter bank of small filters (one per
    // fractional phase between input samples), then for each output sample
    // dot-products the right sub-filter against a window of inputs. Output is
    // split into a leading boundary (filter hangs off the left edge), a multi-
    // threaded interior (no bounds checks), and a trailing boundary.

    int greatest_common_divisor(int a, int b) {
        a = std::abs(a); b = std::abs(b);
        while (b) { int t = b; b = a % b; a = t; }
        return a;
    }

    // Build the filter bank. Designs a windowed-sinc lowpass (Blackman window),
    // then deals its taps round-robin into P sub-filters (the polyphase
    // decomposition). Each sub-filter handles one fractional interpolation
    // phase. Sub-filters are normalized so their taps sum to 1.0, ensuring
    // every output sample has correct amplitude.
    //
    // Returns bank[phase][tap]. During resampling, each output sample picks
    // its phase, grabs that row, and dot-products with the surrounding inputs.
    //
    // halfLobes: number of sinc lobes on each side of the filter center.
    //   More lobes = longer filter = better quality but slower. Caller picks
    //   max(16, max(P, Q) / 2).
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

    // Process output samples [mStart, mEnd) with per-tap bounds checks.
    // Used near the input boundaries where some filter taps would read out
    // of range.
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

    // Interior version: no bounds checks, raw pointer arithmetic. ~99% of the
    // output lands here, and it's the part that gets multithreaded.
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

        // Find safe interior range (no bounds checks needed).
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

        // Leading + trailing boundaries (small, single-threaded).
        processRangeBoundary(inPtr, inLen, outPtr, 0, safeStartM,
            bankPtrs, subLen, filterCenter, P, Q);
        processRangeBoundary(inPtr, inLen, outPtr, safeEndM, outLen,
            bankPtrs, subLen, filterCenter, P, Q);

        // Interior: split across threads.
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

    // Resample a signal from sourceRate to targetRate using polyphase filtering.
    // Picks P (upsample factor) and Q (downsample factor) from the GCD of the
    // two rates, then dispatches to polyphase_resample.
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

    double infer_row_rate(const std::filesystem::path& path,
        const std::string& tsColumnName)
    {
        std::ifstream in(path);
        if (!in) return 0.0;

        std::vector<std::string> hdrs = find_real_header(in);
        if (hdrs.empty()) return 0.0;

        int tsCol = -1;
        for (int i = 0; i < (int)hdrs.size(); ++i) {
            if (contains(hdrs[i], tsColumnName)) { tsCol = i; break; }
        }
        if (tsCol < 0) return 0.0;

        auto parseMs = [](const std::string& s) -> long long {
            size_t sp = s.find_last_of(' ');
            std::string t = (sp == std::string::npos) ? s : s.substr(sp + 1);
            int hh = 0, mm = 0, ss = 0, ms = 0;
            if (std::sscanf(t.c_str(), "%d:%d:%d.%d", &hh, &mm, &ss, &ms) < 3)
                return -1;
            return ((long long)hh * 3600 + mm * 60 + ss) * 1000LL + ms;
            };

        long long firstMs = -1, secondMs = -1;
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::vector<std::string> cells = parse_csv_row(line);
            if (tsCol >= (int)cells.size() || cells[tsCol].empty()) continue;
            long long ms = parseMs(cells[tsCol]);
            if (ms < 0) continue;

            if (firstMs < 0) { firstMs = ms; continue; }
            if (secondMs < 0) {
                if (ms == firstMs) continue;
                secondMs = ms;
                break;
            }
        }
        if (firstMs < 0 || secondMs < 0 || secondMs <= firstMs) return 0.0;

        double dtSec = (secondMs - firstMs) / 1000.0;
        return (dtSec > 0.0) ? (1.0 / dtSec) : 0.0;
    }

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
        out.write((char*)&v, 8);
        sizeUpOut = 1;

        double sentinel[2] = { -1.0, -1.0 };
        out.write((char*)sentinel, 16);
        sizeRawOut = 1;

        nativeRateOut = 0.0f;
    }

    // Write an EDF channel as a (upsampled, raw-pairs) pair.
    //   - Upsampled block: native samples resampled to finalSamplingRate, OR
    //     a copy of raw if the channel is at/below BOOLEAN_RATE.
    //   - Raw block: (timestamp, value) pairs interleaved. EDF samples are
    //     uniformly spaced at `old_rate`, so t_k = k / old_rate.
    // sizeRawOut counts PAIRS (not doubles); byte length = 16 * sizeRawOut.
    // If the channel index is invalid, writes a missing-channel placeholder.
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
            out.write((char*)raw.data(), raw.size() * 8);
            sizeUpOut = (uint32_t)raw.size();
        }
        else {
            std::vector<double> up = upsample(raw, old_rate, finalSamplingRate);
            out.write((char*)up.data(), up.size() * 8);
            sizeUpOut = (uint32_t)up.size();
        }

        const double dt = (old_rate > 0.0) ? (1.0 / old_rate) : 0.0;
        for (size_t k = 0; k < raw.size(); ++k) {
            double pair[2] = { static_cast<double>(k) * dt, raw[k] };
            out.write((char*)pair, 16);
        }
        sizeRawOut = (uint32_t)raw.size();
        nativeRateOut = static_cast<float>(old_rate);
    }

    // Synthesize and write a timestamp channel. Upsampled block is t_k = k /
    // finalSamplingRate; raw block is (k / nativeRate, k / nativeRate) pairs.
    // If durationSec or nativeRate <= 0, writes a missing-channel placeholder.
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
        out.write((char*)up.data(), up.size() * 8);
        sizeUpOut = (uint32_t)up.size();

        const size_t rawLen = (size_t)std::floor(durationSec * nativeRate) + 1;
        const double dtRaw = 1.0 / nativeRate;
        for (size_t k = 0; k < rawLen; ++k) {
            const double t = (double)k * dtRaw;
            double pair[2] = { t, t };
            out.write((char*)pair, 16);
        }
        sizeRawOut = (uint32_t)rawLen;
        nativeRateOut = static_cast<float>(nativeRate);
    }

    // Final step shared by both .edf and .dat paths: seek to start, write the
    // 500-byte header, append the sleep-stage block size, close the file.
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
        out.write((char*)scalars, sizeof(scalars));
        out.write((char*)sizes_up, NUM_CHANNELS * 4);
        out.write((char*)sizes_raw, NUM_CHANNELS * 4);
        out.write((char*)native_rates, NUM_CHANNELS * 4);
        out.write((char*)&sleep_size, 4);
        out.close();
    }

    // ---------- EDF helpers ----------

    // EDF channel map: edf_signal_idx[ChannelIdx] = signal index in the EDF file
    // (or -1 if absent). One label per slot is attempted; a few slots try
    // multiple fallback labels (pacemaker, spo2, resp).
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
        auto findAny = [&](std::initializer_list<const char*> labels) -> int {
            for (const char* l : labels) {
                int idx = find(l);
                if (idx >= 0) return idx;
            }
            return -1;
            };

        m[CH_ECG1] = find(cfg.ecg1Label);
        m[CH_ECG2] = find(cfg.ecg2Label);
        m[CH_ECG3] = find(cfg.ecg3Label);
        m[CH_PPG] = find(cfg.ppgLabel);

        m[CH_ACCEL_X] = find(cfg.accelXLabel);
        m[CH_ACCEL_Y] = find(cfg.accelYLabel);
        m[CH_ACCEL_Z] = find(cfg.accelZLabel);

        m[CH_MARKER] = find("Marker");
        m[CH_TEMP] = find("DEV_Temperature");
        m[CH_PACEMAKER] = findAny({ "Pacemaker", "Pace_Event", "Pace" });

        m[CH_EOG_L] = find("EOG-L");
        m[CH_EOG_R] = find("EOG-R");
        m[CH_EMG] = find("EMG");
        m[CH_EEG1] = find("EEG1");
        m[CH_EEG2] = find("EEG2");
        m[CH_EEG3] = find("EEG3");
        m[CH_EEG4] = find("EEG4");

        m[CH_PRES] = find("Pres");
        m[CH_FLOW] = find("Flow");
        m[CH_THOR] = find("Thor");
        m[CH_ABDO] = find("Abdo");
        m[CH_LEG] = find("Leg");
        m[CH_THERM] = find("Therm");
        m[CH_POS] = find("Pos");

        m[CH_EKG_OFF] = find("EKG_Off");
        m[CH_EOG_L_OFF] = find("EOG-L_Off");
        m[CH_EOG_R_OFF] = find("EOG-R_Off");
        m[CH_EMG_OFF] = find("EMG_Off");
        m[CH_EEG1_OFF] = find("EEG1_Off");
        m[CH_EEG2_OFF] = find("EEG2_Off");
        m[CH_EEG3_OFF] = find("EEG3_Off");

        m[CH_OXSTATUS] = find("OxStatus");
        m[CH_SPO2] = findAny({ "SpO2", "Sp02" });
        m[CH_HR] = find("HR");
        m[CH_DHR] = find("DHR");

        m[CH_RESP] = findAny({ "NLS_NOM_RESP", "Resp" });

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

}   // anonymous namespace

// ============================================================================
// Public entry points
// ============================================================================

void make_binfile_edf(const std::filesystem::path& path,
    const std::filesystem::path& xmlPath,
    const config_entry& cfg)
{
    std::filesystem::path outPath =
        std::filesystem::path(cfg.bin_file_path) /
        (path.stem().string() + "_" +
            std::to_string((int)cfg.finalSamplingRate) + "_" + std::format("{:03d}", static_cast<int>(cfg.bin_length_minutes)) + ".bin");

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
        double rate = (rateOverride > 0.0)
            ? rateOverride
            : edf_channel_rate(hdr.get(), chIdx);
        edf_to_bin(hdr->handle, chIdx, edf_samples(hdr.get(), chIdx),
            rate, cfg.finalSamplingRate, out,
            sizes_up[ch], sizes_raw[ch], native_rates[ch]);
        };

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

    std::vector<double> stages;
    if (!cfg.sleepExt.empty() && !xmlPath.empty()
        && std::filesystem::exists(xmlPath)) {
        pugi::xml_document doc;
        if (doc.load_file(xmlPath.string().c_str())) {
            for (auto node : doc.select_nodes("//SleepStage")) {
                double v = node.node().text().as_double();
                stages.push_back(v == 5.0 ? 4.0 : v);
            }
        }
    }
    if (stages.empty()) stages.push_back(-1.0);
    uint32_t sleep_size = (uint32_t)stages.size();
    out.write((char*)stages.data(), sleep_size * sizeof(double));

    write_header_and_close(out, cfg.finalSamplingRate,
        sizes_up, sizes_raw, native_rates, sleep_size);
}

void make_binfile_dat(const std::filesystem::path& path,
    const config_entry& cfg)
{
    std::filesystem::path outPath =
        std::filesystem::path(cfg.bin_file_path) /
        (path.stem().string() + "_" +
            std::to_string((int)cfg.finalSamplingRate) + "_" + std::format("{:03d}", static_cast<int>(cfg.bin_length_minutes)) + ".bin");

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

    double row_rate = infer_row_rate(path, "Monitor TimeStamp");
    if (row_rate <= 0.0) row_rate = infer_row_rate(path, "System TimeStamp UTC");
    if (row_rate <= 0.0) {
        row_rate = (cfg.ecgRate > 0.0) ? cfg.ecgRate : 500.0;
        std::cout << "  [warn] couldn't infer row rate; using fallback "
            << row_rate << " Hz\n";
    }
    PrescannedDat prescan = prescan_dat_columns(path);

    // Resample a channel to finalSamplingRate using ONLY its real populated
    // samples, placed at their true times. This avoids the "fake plateau" bug
    // where a sparse channel (e.g. PPG populated every 4th row) is treated
    // as a dense signal. Linear interpolation between consecutive real samples.
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

    auto infer_native_rate = [&](const std::vector<size_t>& rawRowIdx) -> double {
        if (rawRowIdx.size() < 2 || row_rate <= 0.0) return 0.0;
        const size_t span = rawRowIdx.back() - rawRowIdx.front();
        if (span == 0) return 0.0;
        const double rowsPerSample = (double)span / (rawRowIdx.size() - 1);
        return (rowsPerSample > 0.0) ? (row_rate / rowsPerSample) : 0.0;
        };

    auto writeCol = [&](const std::string& label, bool exact, ChannelIdx ch) {
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
            out.write((char*)&v, 8);
            sizes_up[ch] = 1;
        }
        else {
            out.write((char*)up.data(), up.size() * 8);
            sizes_up[ch] = (uint32_t)up.size();
        }

        const double dt = (row_rate > 0.0) ? (1.0 / row_rate) : 0.0;
        for (size_t k = 0; k < rawValues.size(); ++k) {
            double pair[2] = { (double)rawRowIdx[k] * dt, rawValues[k] };
            out.write((char*)pair, 16);
        }
        sizes_raw[ch] = (uint32_t)rawValues.size();
        native_rates[ch] = (float)infer_native_rate(rawRowIdx);
        };

    auto writeMissing = [&](ChannelIdx ch) {
        write_missing(out, sizes_up[ch], sizes_raw[ch], native_rates[ch]);
        };

    {
        const double durationSec = (row_rate > 0.0)
            ? (double)prescan.totalRows / row_rate : 0.0;
        write_synthetic_timestamp(out, durationSec, row_rate, finalSR,
            sizes_up[CH_TIMESTAMP], sizes_raw[CH_TIMESTAMP],
            native_rates[CH_TIMESTAMP]);
    }

    writeCol(cfg.ecg1Label, false, CH_ECG1);
    writeCol(cfg.ecg2Label, false, CH_ECG2);
    writeCol(cfg.ecg3Label, false, CH_ECG3);
    writeCol(cfg.ppgLabel, false, CH_PPG);

    for (int ch = CH_ACCEL_X; ch <= CH_EMG; ++ch)
        writeMissing((ChannelIdx)ch);

    writeCol("NLS_EEG_NAMES_EEG_CHAN1", false, CH_EEG1);
    writeCol("NLS_EEG_NAMES_EEG_CHAN2", false, CH_EEG2);
    writeCol("NLS_EEG_NAMES_EEG_CHAN3", false, CH_EEG3);
    writeCol("NLS_EEG_NAMES_EEG_CHAN4", false, CH_EEG4);

    writeCol("NLS_NOM_PRESS_BLD_VEN_CENT", false, CH_PRES);

    for (int ch = CH_FLOW; ch <= CH_DHR; ++ch)
        writeMissing((ChannelIdx)ch);

    writeCol("NLS_NOM_RESP", false, CH_RESP);
    writeCol("NLS_NOM_PRESS_BLD_ART_ABP", true, CH_ABP);
    writeCol("NLS_NOM_PRESS_BLD_ART", true, CH_ART);
    writeCol("NLS_NOM_PRESS_BLD_ART_PULM", true, CH_ART_PULM);

    double placeholder = -1.0;
    out.write((char*)&placeholder, sizeof(double));
    uint32_t sleep_size = 1;

    write_header_and_close(out, finalSR,
        sizes_up, sizes_raw, native_rates, sleep_size);
}

void make_binfile(const std::filesystem::path& path,
    const std::filesystem::path& xmlPath,
    const config_entry& cfg)
{
    std::string ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);

    if (ext == ".EDF") {
        make_binfile_edf(path, xmlPath, cfg);
    }
    else if (ext == ".DAT" || ext == ".CSV") {
        make_binfile_dat(path, cfg);
    }
    else {
        std::cerr << "ERROR: unsupported file type " << ext
            << " for " << path << "\n";
    }
}