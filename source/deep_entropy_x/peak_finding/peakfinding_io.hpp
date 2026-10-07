/**
 * @file   peakfinding_io.hpp
 * @brief  the io for the peak locations binfile and csv
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-20
 */
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <filesystem>
#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include <cmath>
#include <limits>
#include <array>
#include <cstdio>

#include "peakfinding_structs.hpp"

inline constexpr char     annealed_data_magic[8] = { 'A','N','N','L','S','E','G','S' };
inline constexpr uint32_t annealed_data_version = 0;


inline constexpr char     peak_data_magic[8] = { 'R','P','E','A','K','L','O','C' };
inline constexpr uint32_t peak_data_version = 0;

// Shared by the three readers below: consume an 8-byte magic and a uint32
// version, or throw naming the file. ONE PLACE, so a reader cannot check the
// magic and forget the version.
inline void read_and_check_header(std::ifstream& f, const std::string& path, const char(&magic)[8], uint32_t expectVersion, const char* what)
{
    char header[8] = {};
    uint32_t version = 0;
    if (!f.read(header, sizeof(header))
        || !f.read(reinterpret_cast<char*>(&version), 4))
        throw std::runtime_error(std::string("truncated ") + what
            + " header: " + path);
    for (std::size_t i = 0; i < sizeof(header); ++i)
        if (header[i] != magic[i])
            throw std::runtime_error(std::string("not a ") + what
                + " file (bad magic), or written before the format had one: "
                + path);
    if (version != expectVersion)
        throw std::runtime_error(std::string(what) + " version "
            + std::to_string(version) + " != "
            + std::to_string(expectVersion) + " -- regenerate: " + path);
}
/**
 * @brief  Read an annealed-segments .bin file produced by the upstream
 *         preparation step.
 *
 * On-disk layout (all values little-endian, native widths):
 *
 *   Header:
 *     char[8]  magic = kAnnealedMagic ("ANNLSEGS")
 *     uint32   version = kAnnealedVersion (0, the only accepted value)
 *     uint64   numBins                        // segment count
 *     double   filePpgSR                      // PPG sample rate (consumed, not retained)
 *     double   fileEcgSR                      // ECG sample rate (consumed, not retained)
 *     double   scoringEpoch                   // (consumed, not retained)
 *     uint32   nChannels                      // pass-through channel count
 *     uint32   nativeSR[nChannels]            // per-channel native rate (skipped)
 *     uint8    ecg1_inverted                  // "Inverted Lead?" checkbox, CH1 (0/1)
 *     uint8    ecg2_inverted                  // same, CH2
 *     uint8    ecg3_inverted                  // same, CH3
 *
 *   Per bin (numBins of these):
 *     uint64   nPpgPairs
 *     (uint64,uint64) ppg_bin_indexs[nPpgPairs]
 *     uint64   nEcgPairs
 *     (uint64,uint64) ecg_bin_indexs[nEcgPairs]
 *     For each of {ppg_signal, ecg_signal_1, ecg_signal_2, ecg_signal_3, sleep_state_signal}:
 *       uint64   N;  double samples[N]
 *     For each of nChannels pass-through slots:
 *       uint64   nUp;     double upsampled[nUp]
 *       uint64   nPairs;  double raw_tv_interleaved[2 * nPairs]   // (t,v,t,v,...)
 *
 *  The pass-through channels are not consumed by the peak-detection
 *  pipeline. They were previously routed through to the wave_markings
 *  file unchanged; that copy has been dropped, but the read path here
 *  still loads them so the re-hydrating reader in peakfinding_io.hpp
 *  can populate any field a downstream consumer needs.
 */
inline AnnealedData read_input_binfile(const std::string& path) {
    AnnealedData data;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Could not open: " + path);

    // Buffered I/O for throughput.
    char read_buf[1 << 16];
    file.rdbuf()->pubsetbuf(read_buf, sizeof(read_buf));

    read_and_check_header(file, path, annealed_data_magic, annealed_data_version,
        "annealed-segments");

    uint64_t numBins = 0;
    file.read(reinterpret_cast<char*>(&numBins), 8);

    double filePpgSR, fileEcgSR, scoringEpoch;
    file.read(reinterpret_cast<char*>(&filePpgSR), 8);
    file.read(reinterpret_cast<char*>(&fileEcgSR), 8);
    file.read(reinterpret_cast<char*>(&scoringEpoch), 8);

    // nChannels is the per-bin pass-through slot count. It's followed by
    // nChannels uint32 native sample rates which we skip -- the writer
    // recovers the per-segment slot count from bin.all_upsampled.size().
    uint32_t nChannels = 0;
    file.read(reinterpret_cast<char*>(&nChannels), 4);
    if (nChannels > 0)
        file.seekg(static_cast<std::streamoff>(nChannels) * 4, std::ios::cur);

    // Per-channel "Inverted Lead?" checkbox state, one value for the whole
    // file (not per-bin, not auto-detected). Written by the annealed .bin
    // writer immediately after the native-rate array above; read here in
    // the same order.
    uint8_t inv1 = 0, inv2 = 0, inv3 = 0;
    file.read(reinterpret_cast<char*>(&inv1), 1);
    file.read(reinterpret_cast<char*>(&inv2), 1);
    file.read(reinterpret_cast<char*>(&inv3), 1);
    data.ecg1_inverted = inv1 != 0;
    data.ecg2_inverted = inv2 != 0;
    data.ecg3_inverted = inv3 != 0;

    data.bins.resize(numBins);

    for (uint64_t i = 0; i < numBins; ++i) {
        auto& bin = data.bins[i];

        // PPG bin-index pairs: count + (uint64, uint64) records, 16 bytes each.
        uint64_t nPpgPairs;
        if (!file.read(reinterpret_cast<char*>(&nPpgPairs), 8)) break;
        bin.ppg_bin_indexs.resize(nPpgPairs);
        if (nPpgPairs > 0)
            file.read(reinterpret_cast<char*>(bin.ppg_bin_indexs.data()), nPpgPairs * 16);

        // ECG bin-index pairs: same layout as PPG.
        uint64_t nEcgPairs;
        if (!file.read(reinterpret_cast<char*>(&nEcgPairs), 8)) break;
        bin.ecg_bin_indexs.resize(nEcgPairs);
        if (nEcgPairs > 0)
            file.read(reinterpret_cast<char*>(bin.ecg_bin_indexs.data()), nEcgPairs * 16);

        // Helper: read a (uint64 count, double samples[count]) block.
        auto readDoubleArray = [&](std::vector<double>& vec) -> bool {
            uint64_t sz;
            if (!file.read(reinterpret_cast<char*>(&sz), 8)) return false;
            vec.resize(sz);
            if (sz > 0) file.read(reinterpret_cast<char*>(vec.data()), sz * 8);
            return true;
            };

        if (!readDoubleArray(bin.ppg_signal)) break;
        if (!readDoubleArray(bin.ecg_signal_1)) break;
        if (!readDoubleArray(bin.ecg_signal_2)) break;
        if (!readDoubleArray(bin.ecg_signal_3)) break;
        if (!readDoubleArray(bin.sleep_state_signal)) break;

        // Per-segment pass-through: nChannels slots, each as
        //   uint64 nUp;     double upsampled[nUp]
        //   uint64 nPairs;  double raw_tv[2 * nPairs]   // interleaved (t,v)
        // The (t,v) doubles are kept interleaved here so the writer can
        // re-emit them with a single bulk write per slot.
        bin.all_upsampled.resize(nChannels);
        bin.all_raw_pairs_flat.resize(nChannels);
        bool ok = true;
        for (uint32_t ch = 0; ch < nChannels; ++ch) {
            if (!readDoubleArray(bin.all_upsampled[ch])) { ok = false; break; }
            uint64_t nPairs;
            if (!file.read(reinterpret_cast<char*>(&nPairs), 8)) { ok = false; break; }
            bin.all_raw_pairs_flat[ch].resize(nPairs * 2);
            if (nPairs > 0)
                file.read(reinterpret_cast<char*>(bin.all_raw_pairs_flat[ch].data()),
                    nPairs * 2 * sizeof(double));
        }
        if (!ok) break;
    }
    return data;
}

/**
 * @brief  Write per-segment R-peak / PPG-pairing results to a .bin file.
 *
 * On-disk layout (little-endian, native widths). Only fields that are
 * actually computed by the peak-finding step are written; signals and
 * pass-through channels live in the annealed .bin and are reconstructed
 * from there by peakfinding_io.hpp's re-hydrating reader.
 *
 *   uint64   numBins
 *
 *   For each bin (in order):
 *
 *     // R-peak indices: 3 channels x 3 preprocessing methods = 9 arrays.
 *     // Each: uint64 count + count uint64 indices (1-based on disk).
 *     ch1.raw, ch1.squared, ch1.absval
 *     ch2.raw, ch2.squared, ch2.absval
 *     ch3.raw, ch3.squared, ch3.absval
 *
 *     // PPG event indices (same uint64 + 1-based layout as above):
 *     ppgMaxAmps
 *     ppgMinAmps
 *
 *     // Preprocessed signals (uint64 count + count doubles): squared
 *     // then absval, per channel. These are easily recomputed pointwise
 *     // from the raw ECG (x*x and |x|), but we still write them so the
 *     // rebuild path doesn't have to redo that work on every load.
 *     ch1.squared_signal, ch1.absval_signal
 *     ch2.squared_signal, ch2.absval_signal
 *     ch3.squared_signal, ch3.absval_signal
 *
 *     // Pairs (PPG-valley index, ECG-R-peak index):
 *     uint64   numPairs
 *     int64    pairBuf[2 * numPairs]    // interleaved (ppg, ecg); 1-based with -1 sentinel
 *
 * Index 1-basing summary:
 *   1-based on disk (writer adds 1):  R-peak arrays, ppgMaxAmps, ppgMinAmps, pairBuf entries
 *
 * The -1 sentinel in pairBuf marks an unpaired side; NaN or negative
 * doubles in bin.pairs are mapped to -1 on write.
 */
inline void write_output_binfile(const std::string& path, const std::vector<output_binfile_data>& results) {
    std::ofstream file(path, std::ios::binary);
    if (!file.is_open()) return;

    char write_buf[1 << 16];
    file.rdbuf()->pubsetbuf(write_buf, sizeof(write_buf));

    file.write(peak_data_magic, sizeof(peak_data_magic));
    file.write(reinterpret_cast<const char*>(&peak_data_version), 4);

    uint64_t numBins = results.size();
    file.write(reinterpret_cast<const char*>(&numBins), 8);

    for (const auto& bin : results) {

        // Write a size_t index array as uint64, adding 1 for MATLAB 1-based indexing.
        auto writeIdx = [&](const std::vector<std::size_t>& v) {
            uint64_t sz = v.size();
            file.write(reinterpret_cast<const char*>(&sz), 8);
            if (sz > 0) {
                std::vector<uint64_t> tmp(sz);
                for (uint64_t i = 0; i < sz; ++i)
                    tmp[i] = static_cast<uint64_t>(v[i]) + 1;
                file.write(reinterpret_cast<const char*>(tmp.data()), sz * 8);
            }
            };

        // Write a (uint64 count, count doubles) block.
        auto writeSignal = [&](const std::vector<double>& sig) {
            uint64_t sz = sig.size();
            file.write(reinterpret_cast<const char*>(&sz), 8);
            if (sz > 0) file.write(reinterpret_cast<const char*>(sig.data()), sz * 8);
            };

        /* R-peak indices: 3 channels x 3 methods = 9 arrays (all 1-based). */
        writeIdx(bin.ch1.raw);
        writeIdx(bin.ch1.squared);
        writeIdx(bin.ch1.absval);

        writeIdx(bin.ch2.raw);
        writeIdx(bin.ch2.squared);
        writeIdx(bin.ch2.absval);

        writeIdx(bin.ch3.raw);
        writeIdx(bin.ch3.squared);
        writeIdx(bin.ch3.absval);

        /* PPG event indices (1-based). Order: maxAmps then minAmps. */
        writeIdx(bin.ppgMaxAmps);
        writeIdx(bin.ppgMinAmps);

        /* Preprocessed signals: squared then absval, per channel. Kept on
           disk (rather than recomputed on read) to avoid redoing pointwise
           work on every template rebuild. */
        writeSignal(bin.ch1.squared_signal);
        writeSignal(bin.ch1.absval_signal);
        writeSignal(bin.ch2.squared_signal);
        writeSignal(bin.ch2.absval_signal);
        writeSignal(bin.ch3.squared_signal);
        writeSignal(bin.ch3.absval_signal);

        /* Pairs: uint64 count, then count interleaved (int64 ppg, int64 ecg).
           1-based with -1 sentinel for the unpaired side. */
        uint64_t numPairs = bin.pairs.size();
        file.write(reinterpret_cast<const char*>(&numPairs), 8);
        if (numPairs > 0) {
            std::vector<int64_t> pairBuf(numPairs * 2);
            for (uint64_t i = 0; i < numPairs; ++i) {
                const auto& p = bin.pairs[i];
                pairBuf[i * 2] = (std::isnan(p[0]) || p[0] < -0.1) ? -1 : static_cast<int64_t>(std::round(p[0])) + 1;
                pairBuf[i * 2 + 1] = (std::isnan(p[1]) || p[1] < -0.1) ? -1 : static_cast<int64_t>(std::round(p[1])) + 1;
            }
            file.write(reinterpret_cast<const char*>(pairBuf.data()), numPairs * 16);
        }
    }
}

// Per-beat columns, placed right after the CH1 R group. values[bin][slice] is
// that beat's row of cells, parallel to names; NaN or a missing row is blank.
// slice is the R-pair ordinal: the beat cut at ch1.raw[slice], which is row
// `slice` of the ch1_r column. Built by the caller
// (analysis_job::buildBeatMoveColumns), so this writer needs to know nothing
// about what the columns mean.
//
// TEXT COLUMNS. textSlot[col] >= 0 marks a column as text: its cell is
// text[bin][slice][textSlot[col]], blank when empty. Numeric columns read
// values; textSlot shorter than names means the rest are numeric.
struct BeatColumns {
    std::vector<std::string> names;
    std::vector<int> textSlot;
    std::vector<std::vector<std::vector<double>>> values;
    std::vector<std::vector<std::vector<std::string>>> text;
};

inline void write_output_csvfile(const std::string& path, const std::vector<output_binfile_data>& bins, const std::string& fileID, double sampleRateHz = 0.0, double ppgRateHz = 0.0,
    const BeatColumns* beatCols = nullptr)
{
    std::ofstream f(path);
    if (!f.is_open())
        throw std::runtime_error("cannot open for write: " + path);

    constexpr std::size_t COLS_PER_ROW = 11;
    const bool haveRate = sampleRateHz > 0.0;
    const double secPerSample = haveRate ? (1.0 / sampleRateHz) : 0.0;
    // PPG indices are on the PPG clock, which need not be the ECG's. 0 = same as ECG.
    const double ppgSecPerSample = haveRate
        ? ((ppgRateHz > 0.0) ? 1.0 / ppgRateHz : secPerSample) : 0.0;

    // Local PR-segment baseline: median of samples [peak-100, peak-40)
    // (~40-100 ms before the peak at 1000 Hz -- the flat PR interval).
    constexpr std::size_t kBaselineLo = 100;
    constexpr std::size_t kBaselineHi = 40;
    auto localBaseline = [](const std::vector<double>& sig, std::size_t peak) -> double {
        if (sig.empty()) return 0.0;
        const std::size_t lo = (peak > kBaselineLo) ? (peak - kBaselineLo) : 0;
        const std::size_t hi = (peak > kBaselineHi) ? (peak - kBaselineHi) : 0;
        if (hi <= lo) return sig[peak];
        std::vector<double> w(sig.begin() + lo, sig.begin() + hi);
        std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end());
        return w[w.size() / 2];
        };

    // Most recent PPG foot before a given peak sample.
    auto precedingFoot = [](const std::vector<std::size_t>& mins,
        std::size_t peakSample) -> long long {
            long long best = -1;
            for (std::size_t m : mins) {
                if (m < peakSample) best = static_cast<long long>(m);
                else break;
            }
            return best;
        };

    // First PPG onset after a given peak sample: the end of that pulse.
    auto followingOnset = [](const std::vector<std::size_t>& mins,
        std::size_t peakSample) -> long long {
            for (std::size_t m : mins)
                if (m > peakSample) return static_cast<long long>(m);
            return -1;
        };

    // The 50% point: where the upstroke first reaches halfway from onset to
    // peak, linearly interpolated between samples. -1 when it can't be found.
    auto t50Between = [](const std::vector<double>& sig, std::size_t on,
        std::size_t pk) -> double {
            if (pk <= on || pk >= sig.size()) return -1.0;
            const double lo = sig[on], hi = sig[pk];
            if (!std::isfinite(lo) || !std::isfinite(hi) || !(hi > lo)) return -1.0;
            const double level = lo + 0.5 * (hi - lo);
            for (std::size_t k = on + 1; k <= pk; ++k) {
                const double a = sig[k - 1], b = sig[k];
                if (!std::isfinite(a) || !std::isfinite(b) || b < level) continue;
                return (b == a) ? static_cast<double>(k)
                    : static_cast<double>(k - 1) + (level - a) / (b - a);
            }
            return -1.0;
        };

    // t80 and pw80, defined exactly as FeatureMarks defines them for the
    // templates (amplitude_crossing / crossing_at_level in feature_marks.cpp),
    // so the per-beat and per-template numbers mean the same thing.
    //
    // amplitudeCrossing: where the trace first reaches `frac` of the way from
    // v[a] to v[b], interpolated between the straddling samples; with no
    // straddle, the closest sample, as FeatureMarks falls back. -1 on bad input.
    auto amplitudeCrossing = [](const std::vector<double>& v, long long a,
        long long b, double frac) -> double {
            const long long N = static_cast<long long>(v.size());
            if (a < 0 || b <= a || b >= N) return -1.0;
            const double va = v[a], vb = v[b];
            if (std::isnan(va) || std::isnan(vb)) return -1.0;
            const double target = va + frac * (vb - va);
            const bool rising = (vb >= va);
            for (long long i = a + 1; i <= b; ++i) {
                if (std::isnan(v[i]) || std::isnan(v[i - 1])) continue;
                const bool crossed = rising ? (v[i] >= target && v[i - 1] < target)
                    : (v[i] <= target && v[i - 1] > target);
                if (!crossed) continue;
                const double den = v[i] - v[i - 1];
                const double f = (den != 0.0) ? (target - v[i - 1]) / den : 0.0;
                return static_cast<double>(i - 1) + std::clamp(f, 0.0, 1.0);
            }
            long long best = a; double bestDiff = std::numeric_limits<double>::infinity();
            for (long long i = a; i <= b; ++i) {
                if (std::isnan(v[i])) continue;
                const double d = std::abs(v[i] - target);
                if (d < bestDiff) { bestDiff = d; best = i; }
            }
            return static_cast<double>(best);
        };
    // crossingAtLevel: first crossing of an absolute level between a and b,
    // interpolated; -1 when there is none.
    auto crossingAtLevel = [](const std::vector<double>& v, long long a,
        long long b, double target) -> double {
            const long long N = static_cast<long long>(v.size());
            if (a < 0 || b <= a || b >= N) return -1.0;
            const double va = v[a], vb = v[b];
            if (std::isnan(va) || std::isnan(vb) || std::isnan(target)) return -1.0;
            const bool rising = (vb >= va);
            for (long long i = a + 1; i <= b; ++i) {
                if (std::isnan(v[i]) || std::isnan(v[i - 1])) continue;
                const bool crossed = rising ? (v[i] >= target && v[i - 1] < target)
                    : (v[i] <= target && v[i - 1] > target);
                if (!crossed) continue;
                const double den = v[i] - v[i - 1];
                const double f = (den != 0.0) ? (target - v[i - 1]) / den : 0.0;
                return static_cast<double>(i - 1) + std::clamp(f, 0.0, 1.0);
            }
            return -1.0;
        };

    static const char* const column_suffix[COLS_PER_ROW] = {
        "ch1_r", "ch1_squared_r", "ch1_absval_r",
        "ch2_r", "ch2_squared_r", "ch2_absval_r",
        "ch3_r", "ch3_squared_r", "ch3_absval_r",
        "ppg_onset", "ppg_peak"
    };
    static const char* const interval_column_name[COLS_PER_ROW] = {
        "ch1_rr_interval_ms", "ch1_squared_rr_interval_ms", "ch1_absval_rr_interval_ms",
        "ch2_rr_interval_ms", "ch2_squared_rr_interval_ms", "ch2_absval_rr_interval_ms",
        "ch3_rr_interval_ms", "ch3_squared_rr_interval_ms", "ch3_absval_rr_interval_ms",
        "ppg_onset_interval_ms", "ppg_peak_interval_ms"
    };
    // The two ways every time is written, header and cells alike.
    static const char* const kTimeSuffix[2] = {
        "_ms_from_bin_start", "_ms_from_file_start" };
    auto timeHeader = [&](const char* prefix) {
        for (const char* s : kTimeSuffix) f << ',' << prefix << s;
        };

    // Header: file_id, bin, then 4 cells per column, then the pulse end and
    // onset -> 50% point next to the last group, ppg_peak. The CH1 R group
    // is followed by the per-beat block (see BeatColumns above).
    f << "file_id,bin";
    for (std::size_t c = 0; c < COLS_PER_ROW; ++c) {
        timeHeader(column_suffix[c]);
        // A height per ECG R (mV) and per PPG peak (sensor units, so no unit
        // in the name). None for the onset: it is the baseline end of the
        // peak's amplitude, not a height of its own.
        if (c < 9)        f << ',' << column_suffix[c] << "_peak_height_mv";
        else if (c == 10) f << ",ppg_peak_amplitude";
        f << ',' << interval_column_name[c];
        if (c == 0 && beatCols)
            for (const std::string& nm : beatCols->names) f << ',' << nm;
    }
    timeHeader("ppg_end");
    f << ",ppg_onset_to_t50_ms,ppg_onset_to_t80_ms,ppg_pw80_ms\n";

    struct ColBundle {
        const std::vector<std::size_t>* peaks;
        const std::vector<double>* signal;
    };

    // Running "seconds elapsed in prior bins" -- the start time of the
    // current bin relative to the whole recording. Advances at the end of
    // each bin by that bin's ECG length in seconds.
    double binStartSec = 0.0;

    // The two time cells, in kTimeSuffix's order; blank when !has.
    auto timeCells = [&](bool has, double binSec) {
        const double v[2] = { binSec * 1000.0, (binStartSec + binSec) * 1000.0 };
        for (const double x : v) { f << ','; if (has) f << x; }
        };

    for (std::size_t b = 0; b < bins.size(); ++b) {
        const auto& bin = bins[b];

        std::array<ColBundle, COLS_PER_ROW> cols = { {
            {&bin.ch1.raw,     &bin.ecgSignal},
            {&bin.ch1.squared, &bin.ecgSignal},
            {&bin.ch1.absval,  &bin.ecgSignal},
            {&bin.ch2.raw,     &bin.ecgSignal2},
            {&bin.ch2.squared, &bin.ecgSignal2},
            {&bin.ch2.absval,  &bin.ecgSignal2},
            {&bin.ch3.raw,     &bin.ecgSignal3},
            {&bin.ch3.squared, &bin.ecgSignal3},
            {&bin.ch3.absval,  &bin.ecgSignal3},
            {&bin.ppgMinAmps,  &bin.ppgSignal},
            {&bin.ppgMaxAmps,  &bin.ppgSignal}
        } };

        std::size_t maxLen = 0;
        for (const auto& c : cols)
            if (c.peaks->size() > maxLen) maxLen = c.peaks->size();

        for (std::size_t r = 0; r < maxLen; ++r) {
            f << fileID << ',' << b;

            for (std::size_t c = 0; c < COLS_PER_ROW; ++c) {
                const auto& bundle = cols[c];
                const auto& peaks = *bundle.peaks;
                const bool has = r < peaks.size();
                const std::size_t idx0 = has ? peaks[r] : 0;
                const double sps = (c >= 9) ? ppgSecPerSample : secPerSample;

                const double binSec = has ? (idx0 * sps) : 0.0;
                timeCells(has && haveRate, binSec);

                // peak_height -- no cell at all for the onset (see header)
                //   ECG columns: signal[peak] minus local PR-segment baseline.
                //   ppg_peak:    ppgSignal[peak] - ppgSignal[preceding onset].
                if (c != 9) f << ',';
                if (c != 9 && has && bundle.signal && idx0 < bundle.signal->size()) {
                    if (c == 10) {
                        const long long foot = precedingFoot(bin.ppgMinAmps, idx0);
                        if (foot >= 0 &&
                            static_cast<std::size_t>(foot) < bundle.signal->size()) {
                            f << ((*bundle.signal)[idx0]
                                - (*bundle.signal)[static_cast<std::size_t>(foot)]);
                        }
                    }
                    else {
                        const double baseline = localBaseline(*bundle.signal, idx0);
                        f << ((*bundle.signal)[idx0] - baseline);
                    }
                }

                // rinterval_ms: this peak minus previous in same column, in ms.
                f << ',';
                if (has && haveRate && r > 0) {
                    const std::size_t prev0 = peaks[r - 1];
                    f << ((static_cast<double>(idx0) - static_cast<double>(prev0))
                        * sps * 1000.0);
                }

                // CH1 R's beat: the caller's per-beat block, same key.
                if (c == 0) {
                    if (beatCols) {
                        const std::vector<double>* cells = nullptr;
                        if (has && b < beatCols->values.size() && r < beatCols->values[b].size()
                            && !beatCols->values[b][r].empty())
                            cells = &beatCols->values[b][r];
                        const std::vector<std::string>* texts = nullptr;
                        if (has && b < beatCols->text.size() && r < beatCols->text[b].size()
                            && !beatCols->text[b][r].empty())
                            texts = &beatCols->text[b][r];
                        for (std::size_t k = 0; k < beatCols->names.size(); ++k) {
                            f << ',';
                            const int ts = (k < beatCols->textSlot.size()) ? beatCols->textSlot[k] : -1;
                            if (ts >= 0) {
                                if (texts && static_cast<std::size_t>(ts) < texts->size()) f << (*texts)[ts];
                            }
                            else if (cells && k < cells->size() && std::isfinite((*cells)[k])) f << (*cells)[k];
                        }
                    }
                }
            }

            // This row's pulse: its end (the next onset after the peak) and
            // onset -> t50 (onset = the onset preceding the peak).
            {
                const auto& pk = bin.ppgMaxAmps;
                const bool hasPk = haveRate && r < pk.size();
                const long long end = hasPk ? followingOnset(bin.ppgMinAmps, pk[r]) : -1;
                const double endBin = (end >= 0) ? static_cast<double>(end) * ppgSecPerSample : 0.0;
                timeCells(end >= 0, endBin);
                f << ',';
                const long long on = hasPk ? precedingFoot(bin.ppgMinAmps, pk[r]) : -1;
                if (hasPk) {
                    const double t = (on >= 0)
                        ? t50Between(bin.ppgSignal, static_cast<std::size_t>(on), pk[r]) : -1.0;
                    if (t >= 0.0)
                        f << (t - static_cast<double>(on)) * ppgSecPerSample * 1000.0;
                }

                // t80: where the downstroke has fallen 80% of the way from the
                // peak to the end. Written as its distance from the onset.
                // pw80: the pulse's width at t80's level -- t80 minus where the
                // upstroke first crosses that same level.
                const std::vector<double>& sig = bin.ppgSignal;
                const long long ipk = hasPk ? static_cast<long long>(pk[r]) : -1;
                const double t80 = (hasPk && end > ipk)
                    ? amplitudeCrossing(sig, ipk, end, 0.80) : -1.0;
                f << ',';
                if (t80 >= 0.0 && on >= 0)
                    f << (t80 - static_cast<double>(on)) * ppgSecPerSample * 1000.0;
                f << ',';
                if (t80 >= 0.0 && on >= 0 && ipk > on
                    && static_cast<std::size_t>(end) < sig.size()) {
                    const double level = sig[ipk] + 0.80 * (sig[end] - sig[ipk]);
                    const double rise = crossingAtLevel(sig, on, ipk, level);
                    if (rise >= 0.0 && t80 > rise)
                        f << (t80 - rise) * ppgSecPerSample * 1000.0;
                }
            }
            f << '\n';
        }

        // Advance the running total for the next bin. Use the deepest of
        // the three ECG channels so a missing ch1 doesn't zero the step.
        if (haveRate) {
            const std::size_t nSamples = std::max({
                bin.ecgSignal.size(),
                bin.ecgSignal2.size(),
                bin.ecgSignal3.size() });
            binStartSec += static_cast<double>(nSamples) * secPerSample;
        }
    }
}

// ===========================================================================
// MOVED IN FROM THE OLD peakfinding_io.hpp.
//
// read_input_binfile above is what the two-arg overload needs, which is why
// these belong here: in the old layout it had to be forward-declared across
// the header boundary. The one-arg overload is defined before the two-arg one
// because the latter calls it.
// ===========================================================================

/**
 * @brief  Read a wave_markings .bin (R-peaks-only layout).
 *
 *         On-disk layout (see write_output_binfile above in this file
 *         for the authoritative spec):
 *
 *           char[8] magic = kPeaksMagic ("RPEAKLOC")
 *           uint32  version = kPeaksVersion (0)
 *           uint64 numBins
 *           For each bin:
 *             9 index arrays (ch1/2/3 x raw/squared/absval), each
 *               uint64 count + count * uint64 indices (1-based on disk)
 *             ppgMaxAmps, ppgMinAmps   (same uint64-count + 1-based layout)
 *             6 preprocessed signals (ch1/2/3 x squared/absval), each
 *               uint64 count + count * double samples
 *             uint64 numPairs
 *             int64 pairBuf[2 * numPairs]    (interleaved ppg, ecg; -1 sentinel)
 *
 *         The raw signals, bin-index ranges and pass-through channels are
 *         NOT in this file. Use the overload below that takes the annealed
 *         path if the caller needs them.
 */
inline std::vector<output_binfile_data> read_output_binfile(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("cannot open: " + path);

    char rbuf[1 << 16];
    f.rdbuf()->pubsetbuf(rbuf, sizeof(rbuf));

    read_and_check_header(f, path, peak_data_magic, peak_data_version, "R-peak");

    uint64_t numBins = 0;
    f.read(reinterpret_cast<char*>(&numBins), 8);
    std::vector<output_binfile_data> bins(numBins);

    auto readIdx = [&](std::vector<std::size_t>& v) {
        uint64_t sz;
        f.read(reinterpret_cast<char*>(&sz), 8);
        v.resize(sz);
        if (sz > 0) {
            std::vector<uint64_t> tmp(sz);
            f.read(reinterpret_cast<char*>(tmp.data()), sz * 8);
            // Writer added 1 for MATLAB compatibility; subtract back to 0-based.
            for (uint64_t i = 0; i < sz; ++i)
                v[i] = static_cast<std::size_t>(tmp[i]) - 1;
        }
        };
    auto readSig = [&](std::vector<double>& v) {
        uint64_t sz;
        f.read(reinterpret_cast<char*>(&sz), 8);
        v.resize(sz);
        if (sz > 0) f.read(reinterpret_cast<char*>(v.data()), sz * 8);
        };

    for (uint64_t i = 0; i < numBins; ++i) {
        auto& b = bins[i];
        b.index = i;

        // 9 R-peak index arrays (3 channels x 3 methods)
        readIdx(b.ch1.raw);     readIdx(b.ch1.squared);   readIdx(b.ch1.absval);
        readIdx(b.ch2.raw);     readIdx(b.ch2.squared);   readIdx(b.ch2.absval);
        readIdx(b.ch3.raw);     readIdx(b.ch3.squared);   readIdx(b.ch3.absval);

        // PPG indices
        readIdx(b.ppgMaxAmps);  readIdx(b.ppgMinAmps);

        // 6 preprocessed signals (squared/absval per channel)
        readSig(b.ch1.squared_signal); readSig(b.ch1.absval_signal);
        readSig(b.ch2.squared_signal); readSig(b.ch2.absval_signal);
        readSig(b.ch3.squared_signal); readSig(b.ch3.absval_signal);

        // pairs (PPG - ECG matching), written as int64 with -1 sentinel
        uint64_t numPairs;
        f.read(reinterpret_cast<char*>(&numPairs), 8);
        b.pairs.resize(numPairs, std::vector<double>(2));
        if (numPairs > 0) {
            std::vector<int64_t> tmp(numPairs * 2);
            f.read(reinterpret_cast<char*>(tmp.data()), numPairs * 16);
            for (uint64_t k = 0; k < numPairs; ++k) {
                b.pairs[k][0] = (tmp[k * 2] == -1) ? -1.0
                    : static_cast<double>(tmp[k * 2] - 1);
                b.pairs[k][1] = (tmp[k * 2 + 1] == -1) ? -1.0
                    : static_cast<double>(tmp[k * 2 + 1] - 1);
            }
        }
    }
    return bins;
}


/**
 * @brief  Read a wave_markings .bin AND re-hydrate the raw signals from
 *         the matching annealed .bin.
 *
 *         The wave_markings file holds R-peak indices, PPG event indices,
 *         the preprocessed (squared/absval) ECG channels and pairs.
 *         Template generation also needs the raw ECG/PPG signals
 *         (to extract beats around each R-peak) and the bin-index ranges;
 *         those live in the annealed .bin. This overload reads both
 *         files and stitches the fields together so the returned vector
 *         matches what make_beats.hpp's create_ecg_ppg_pairs_raw()
 *         produces in memory.
 *
 *         Per-bin re-hydration:
 *           ecgSignal     <- annealed.ecg_signal_1
 *           ecgSignal2    <- annealed.ecg_signal_2
 *           ecgSignal3    <- annealed.ecg_signal_3
 *           ppgSignal     <- annealed.ppg_signal
 *           ppg_bin_indexs / ecg_bin_indexs <- annealed
 *
 *         The preprocessed squared_signal / absval_signal are NOT
 *         recomputed here: they come straight off the wave_markings .bin
 *         (already loaded by the single-arg overload above).
 *
 *         The bin count in the two files must match; otherwise this
 *         throws.
 *
 * @param  wavePath       Path to the wave_markings .bin.
 * @param  annealedPath   Path to the matching annealed .bin.
 */
inline std::vector<output_binfile_data> read_output_binfile(const std::string& wavePath, const std::string& annealedPath)
{
    std::vector<output_binfile_data> bins = read_output_binfile(wavePath);
    AnnealedData ann = read_input_binfile(annealedPath);

    if (ann.bins.size() != bins.size()) {
        throw std::runtime_error(
            "bin-count mismatch between wave_markings and annealed .bin: "
            + std::to_string(bins.size()) + " vs "
            + std::to_string(ann.bins.size()));
    }

    for (std::size_t i = 0; i < bins.size(); ++i) {
        auto& b = bins[i];
        auto& a = ann.bins[i];

        b.ppgSignal = std::move(a.ppg_signal);
        b.ecgSignal = std::move(a.ecg_signal_1);
        b.ecgSignal2 = std::move(a.ecg_signal_2);
        b.ecgSignal3 = std::move(a.ecg_signal_3);

        b.ppg_bin_indexs = std::move(a.ppg_bin_indexs);
        b.ecg_bin_indexs = std::move(a.ecg_bin_indexs);

        // Arterial pass-through channels, pulled from the annealed slot set
        // for use as background-context traces in the template viewer. Slot
        // indices match the file_to_bin / gui_handler channel layout:
        //   CH_ABP = 33, CH_ART = 34, CH_ART_PULM = 35.
        // A missing/short slot yields an empty vector (=> not drawn).
        auto slot = [&](std::size_t ch) -> std::vector<double> {
            return (ch < a.all_upsampled.size())
                ? std::move(a.all_upsampled[ch]) : std::vector<double>{};
            };
        b.abpSignal = slot(33);
        b.artSignal = slot(34);
        b.artPulmSignal = slot(35);

        // WAS qDebug(). Converted to stderr so this header needs no Qt:
        // <QDebug> was the only Qt dependency in the peak-finding subsystem,
        // and it reached build_bins.hpp and the whole template pipeline
        // through the struct definitions. Same tag, same fields.
        if (i < 3)
            std::fprintf(stderr,
                "  [REHYDRATE] bin %zu all_upsampled.size %zu abp(33) %zu "
                "art(34) %zu artPulm(35) %zu\n",
                static_cast<std::size_t>(i), a.all_upsampled.size(),
                b.abpSignal.size(), b.artSignal.size(), b.artPulmSignal.size());

        // Pass-through channels stay empty: template generation does not
        // read them (the previous reader was already seeking past them).
    }
    return bins;
}