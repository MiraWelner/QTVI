/**
 * @file   peakfinding_structs.hpp
 * @brief  The data types the R-peak / PPG-pairing stage passes around.
 *         TYPES ONLY -- no reading, no writing, no Qt.
 *
 *         Was peakfinding_io.hpp, which held these structs AND two
 *         read_output_binfile overloads. Those have moved to
 *         peakfinding_io.hpp (formerly run_find_r_peaks.hpp), which is
 *         where every other reader and writer already lived. Three
 *         consequences, and they are the reason for the split:
 *
 *           - The forward declaration of read_input_binfile is gone. It
 *             existed because the two-arg read_output_binfile needed a
 *             function defined in the other file, while that file included
 *             this one -- a cycle broken with a forward declaration of an
 *             inline function. Both readers now sit beside it instead.
 *
 *           - <QDebug> is gone. qDebug() was used in exactly one place, the
 *             re-hydration trace inside the two-arg reader. Because
 *             build_bins.hpp includes this header for output_binfile_data,
 *             that one line pulled Qt into the entire template pipeline.
 *
 *           - <fstream> and <stdexcept> are gone with the readers.
 *
 *         So make_beats.hpp and build_bins.hpp now depend on the types and
 *         nothing else; only analysis_job.hpp needs the I/O header.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-30
 */
#pragma once

#include <vector>
#include <string>
#include <utility>
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
    // since they're already on disk in the annealed .bin). 
    std::vector<std::vector<double>> all_upsampled;          // per-slot upsampled samples
    std::vector<std::vector<double>> all_raw_pairs_flat;     // per-slot interleaved (t, v, t, v, ...)
};


struct AnnealedData {
    std::vector<AnnealedSegment> bins;

    // Per-channel "Inverted Lead?" checkbox state (noise_marking_gui.ui:
    // ecg_1_reverse / ecg_2_reverse / ecg_3_reverse), one value for the
    // whole file -- not auto-detected, not per-bin. Read from the annealed
    // .bin header by read_input_binfile() in peakfinding_io.hpp; must be
    // written at the matching header position by the annealed .bin writer.
    bool ecg1_inverted = false;
    bool ecg2_inverted = false;
    bool ecg3_inverted = false;
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

    // Raw arterial signals, for display as faint background-context traces
    // in the template viewer. Present only when the dataset carried them;
    // re-hydrated from the annealed .bin's pass-through slots by the
    // two-arg read_output_binfile below. Empty otherwise.
    std::vector<double> abpSignal;
    std::vector<double> artSignal;
    std::vector<double> artPulmSignal;

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