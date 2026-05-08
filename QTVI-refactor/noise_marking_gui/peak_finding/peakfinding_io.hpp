/**
 * @file   peakfinding_io.hpp
 * @brief  Structs for input/output .bin files and their I/O functions.
 *         Each ECG channel stores R-peak indices from three detection
 *         methods: raw, squared, and absolute-value preprocessing,
 *         along with the preprocessed signals themselves.
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
    // to the wave_markings output. Indexed by the same 41-slot layout as
    // step 3's annealed bin (slot 0 = timestamp, 1 = ECG1, ..., 4 = PPG,
    // 5..40 = other channels).
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

    // Raw signals (always stored unmodified)
    std::vector<double> ecgSignal;
    std::vector<double> ecgSignal2;
    std::vector<double> ecgSignal3;
    std::vector<double> ppgSignal;

    // Per-channel R-peaks + preprocessed signals from 3 methods
    ChannelRPeaks ch1;
    ChannelRPeaks ch2;
    ChannelRPeaks ch3;

    std::vector<std::size_t> ppgMinAmps;
    std::vector<std::size_t> ppgMaxAmps;
    std::size_t index = 0;
    std::vector<std::pair<uint64_t, uint64_t>> ppg_bin_indexs;
    std::vector<std::pair<uint64_t, uint64_t>> ecg_bin_indexs;

    // Pass-through channels routed from the annealed input through to the
    // wave_markings output without modification. Same 41-slot layout as
    // AnnealedSegment.all_upsampled / all_raw_pairs_flat.
    std::vector<std::vector<double>> all_upsampled;
    std::vector<std::vector<double>> all_raw_pairs_flat;
};


inline std::vector<output_binfile_data> read_output_binfile(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("cannot open: " + path);

    char rbuf[1 << 16];
    f.rdbuf()->pubsetbuf(rbuf, sizeof(rbuf));

    uint64_t numBins = 0;
    f.read(reinterpret_cast<char*>(&numBins), 8);
    std::vector<output_binfile_data> bins(numBins);

    auto readIdx = [&](std::vector<std::size_t>& v) {
        uint64_t sz; f.read(reinterpret_cast<char*>(&sz), 8);
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
        uint64_t sz; f.read(reinterpret_cast<char*>(&sz), 8);
        v.resize(sz);
        if (sz > 0) f.read(reinterpret_cast<char*>(v.data()), sz * 8);
        };
    auto readPairVec = [&](std::vector<std::pair<uint64_t, uint64_t>>& v) {
        uint64_t sz; f.read(reinterpret_cast<char*>(&sz), 8);
        v.resize(sz);
        if (sz > 0) f.read(reinterpret_cast<char*>(v.data()), sz * 16);
        };

    for (uint64_t i = 0; i < numBins; ++i) {
        auto& b = bins[i];
        b.index = i;

        // 9 R-peak index arrays (3 channels × 3 methods)
        readIdx(b.ch1.raw);     readIdx(b.ch1.squared);   readIdx(b.ch1.absval);
        readIdx(b.ch2.raw);     readIdx(b.ch2.squared);   readIdx(b.ch2.absval);
        readIdx(b.ch3.raw);     readIdx(b.ch3.squared);   readIdx(b.ch3.absval);

        // PPG indices
        readIdx(b.ppgMaxAmps);  readIdx(b.ppgMinAmps);

        // Raw signals
        readSig(b.ppgSignal); readSig(b.ecgSignal);
        readSig(b.ecgSignal2); readSig(b.ecgSignal3);

        // 6 preprocessed signals (squared/absval per channel)
        readSig(b.ch1.squared_signal); readSig(b.ch1.absval_signal);
        readSig(b.ch2.squared_signal); readSig(b.ch2.absval_signal);
        readSig(b.ch3.squared_signal); readSig(b.ch3.absval_signal);

        // 9 noise-flag bytes
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

        readPairVec(b.ppg_bin_indexs);
        readPairVec(b.ecg_bin_indexs);

        // 40 pass-through channel slices: each is (upsampled doubles, then
        // raw t,v pair doubles). buildTemplatesFromPeakResults doesn't read
        // these, so we just seek past them.
        constexpr int PASS_NCH = 40;
        for (int ch = 0; ch < PASS_NCH; ++ch) {
            uint64_t nUp;
            f.read(reinterpret_cast<char*>(&nUp), 8);
            f.seekg(static_cast<std::streamoff>(nUp * 8), std::ios::cur);
            uint64_t nP;
            f.read(reinterpret_cast<char*>(&nP), 8);
            f.seekg(static_cast<std::streamoff>(nP * 16), std::ios::cur);
        }
    }
    return bins;
}