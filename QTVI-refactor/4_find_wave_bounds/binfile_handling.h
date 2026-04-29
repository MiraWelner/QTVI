/**
 * @file   binfile_handling.h
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

AnnealedData read_input_binfile(const std::string& path);
void write_output_binfile(const std::string& path, const std::vector<output_binfile_data>& results);