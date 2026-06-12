/**
 * @file   peakfinding_io.hpp
 * @brief  Structs for input/output .bin files and their I/O functions.
 *         Each ECG channel stores R-peak indices from three detection
 *         methods: raw, squared, and absolute-value preprocessing,
 *         along with the preprocessed signals themselves.
 *
 *         The wave_markings .bin (written by write_output_binfile in
 *         run_find_r_peaks.hpp) was historically large because it
 *         duplicated the raw signals, the PPG/ECG bin-index ranges, and
 *         the full 40-slot pass-through channel set that the annealed
 *         .bin already contains.
 *
 *         The current on-disk format keeps everything the peakfinding
 *         step actually computes -- R-peak indices, PPG min/max-amp
 *         indices, the squared/absval preprocessed ECG channels, noise
 *         flags, and PPG-to-ECG pairs -- but drops the raw signals and
 *         pass-through channels. When a downstream consumer needs the
 *         raw signals or bin-index ranges (e.g. template generation
 *         rebuilding from disk), call read_output_binfile() with the
 *         annealed .bin path to re-hydrate them from there.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-30
 */
#pragma once

#include <vector>
#include <string>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <cmath>


struct AnnealedSegment {
    std::vector<double> ppg_signal;
    std::vector<double> ecg_signal_1;
    std::vector<double> ecg_signal_2;
    std::vector<double> ecg_signal_3;
    std::vector<double> sleep_state_signal;
    std::vector<std::pair<uint64_t, uint64_t>> ppg_bin_indexs;
    std::vector<std::pair<uint64_t, uint64_t>> ecg_bin_indexs;

    // Pass-through: full set of input channels carried alongside the
    // algorithm-facing signals above. The peakfinding algorithm doesn't
    // touch these -- they're just routed from the annealed input through
    // (formerly into the wave_markings output, now dropped from there
    // since they're already on disk in the annealed .bin). Indexed by the
    // same 41-slot layout as step 3's annealed bin (slot 0 = timestamp,
    // 1 = ECG1, ..., 4 = PPG, 5..40 = other channels).
    std::vector<std::vector<double>> all_upsampled;          // per-slot upsampled samples
    std::vector<std::vector<double>> all_raw_pairs_flat;     // per-slot interleaved (t, v, t, v, ...)
};


struct AnnealedData {
    std::vector<AnnealedSegment> bins;
};

/**
 * @brief  R-peak results and preprocessed signals for a single channel
 *         across three preprocessing methods.
 */
struct ChannelRPeaks {
    std::vector<std::size_t> raw;       ///< R-peaks from unmodified signal
    std::vector<std::size_t> squared;   ///< R-peaks from squared signal
    std::vector<std::size_t> absval;    ///< R-peaks from abs-value signal

    std::vector<double> squared_signal; ///< The squared version of the ECG channel
    std::vector<double> absval_signal;  ///< The abs-value version of the ECG channel

    bool raw_noisy = false;
    bool squared_noisy = false;
    bool absval_noisy = false;
};

struct output_binfile_data {
    std::vector<std::vector<double>> pairs;
    bool bad_segment = false;

    // Raw signals (always stored unmodified in-memory; not serialized to
    // wave_markings .bin -- the annealed .bin is the source of truth for
    // these, and read_output_binfile() can re-hydrate them from there).
    std::vector<double> ecgSignal;
    std::vector<double> ecgSignal2;
    std::vector<double> ecgSignal3;
    std::vector<double> ppgSignal;

    // Per-channel R-peaks + preprocessed signals from 3 methods. The
    // R-peak indices, noise flags, and preprocessed (squared/absval)
    // signals are all serialized to the wave_markings .bin. The raw ECG
    // they're derived from lives in the annealed .bin.
    ChannelRPeaks ch1;
    ChannelRPeaks ch2;
    ChannelRPeaks ch3;

    std::vector<std::size_t> ppgMinAmps;
    std::vector<std::size_t> ppgMaxAmps;
    std::size_t index = 0;

    // PPG/ECG bin-index ranges live in the annealed .bin; populated here
    // only when read_output_binfile() is given the annealed path.
    std::vector<std::pair<uint64_t, uint64_t>> ppg_bin_indexs;
    std::vector<std::pair<uint64_t, uint64_t>> ecg_bin_indexs;

    // Pass-through channels live exclusively in the annealed .bin. These
    // fields stay in the struct for the in-memory (peakResultsInMemory)
    // path that hands peakResults straight to template generation
    // without a disk round-trip; they are never populated by the
    // wave_markings reader.
    std::vector<std::vector<double>> all_upsampled;
    std::vector<std::vector<double>> all_raw_pairs_flat;
};


/**
 * @brief  Read a wave_markings .bin (R-peaks-only layout).
 *
 *         On-disk layout (see write_output_binfile in run_find_r_peaks.hpp
 *         for the authoritative spec):
 *
 *           uint64 numBins
 *           For each bin:
 *             9 index arrays (ch1/2/3 x raw/squared/absval), each
 *               uint64 count + count * uint64 indices (1-based on disk)
 *             ppgMaxAmps, ppgMinAmps   (same uint64-count + 1-based layout)
 *             6 preprocessed signals (ch1/2/3 x squared/absval), each
 *               uint64 count + count * double samples
 *             uint8 flags[9]           (noise flags, channel-major then method)
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

        // 9 noise-flag bytes (3 channels x 3 methods)
        uint8_t flags[9];
        f.read(reinterpret_cast<char*>(flags), 9);
        b.ch1.raw_noisy = flags[0]; b.ch1.squared_noisy = flags[1]; b.ch1.absval_noisy = flags[2];
        b.ch2.raw_noisy = flags[3]; b.ch2.squared_noisy = flags[4]; b.ch2.absval_noisy = flags[5];
        b.ch3.raw_noisy = flags[6]; b.ch3.squared_noisy = flags[7]; b.ch3.absval_noisy = flags[8];

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


// Forward-declared; defined in run_find_r_peaks.hpp. Reads an annealed
// segments .bin into AnnealedData. Declared here so the re-hydrating
// reader below can use it without inverting the header dependency.
inline AnnealedData read_input_binfile(const std::string& path);


/**
 * @brief  Read a wave_markings .bin AND re-hydrate the raw signals from
 *         the matching annealed .bin.
 *
 *         The wave_markings file holds R-peak indices, PPG event indices,
 *         the preprocessed (squared/absval) ECG channels, noise flags and
 *         pairs. Template generation also needs the raw ECG/PPG signals
 *         (to extract beats around each R-peak) and the bin-index ranges;
 *         those live in the annealed .bin. This overload reads both
 *         files and stitches the fields together so the returned vector
 *         matches what create_ecg_ppg_pairs() produces in memory.
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
 *         throws. bad_segment is not on disk in either file -- it stays
 *         false on the rebuild path, matching the previous behaviour.
 *
 * @param  wavePath       Path to the wave_markings .bin.
 * @param  annealedPath   Path to the matching annealed .bin.
 */
inline std::vector<output_binfile_data> read_output_binfile(
    const std::string& wavePath,
    const std::string& annealedPath)
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

        // Pass-through channels stay empty: template generation does not
        // read them (the previous reader was already seeking past them).
    }
    return bins;
}