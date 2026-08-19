#pragma once
/**
 * @file   file_to_bin.hpp
 * @brief  Uses polyphase upsampling with a filter bank for speed and efficiency to convert an .edf or .dat file to a
 *         .bin file. Every channel is upsampled to its OWN target rate,
 *         taken from the per-channel *_upsampled_rate columns of
 *         config.csv (there is no longer a single global rate or a <=1 Hz
 *         boolean tier). The per-channel target rates are recorded in the
 *         header so a consumer can reconstruct each channel's upsampled
 *         time axis.
 *
 *   592-byte header (148 x 32-bit fields):
 *
 *     Offset   0: header_version  (uint32)  = BIN_HEADER_VERSION (currently 1)
 *     Offset   4: n_channels      (uint32)  = NUM_CHANNELS (currently 36)
 *     Offset   8: sleep_state_len (uint32)  sleep-stage epoch length in seconds (e.g. 30 s)
 *
 *     Offset  12: upsampled sizes       36 x uint32  (size_<chan>)
 *     Offset 156: raw-pair sizes        36 x uint32  (size_<chan>_raw, counts PAIRS)
 *     Offset 300: native sampling rates 36 x float32 (Hz)
 *                 0.0 = channel absent; negative values never used.
 *     Offset 444: upsampled rates       36 x float32 (Hz)
 *                 0.0 = channel absent (matches a missing-channel placeholder).
 *     Offset 588: size_sleep            (uint32)
 *
 *   Header size check: 2 + 1 + 36 + 36 + 36 + 36 + 1 = 148 fields x 4 bytes = 592 bytes
 *   (2 = version + n_channels)
 *
 *   Channel index order (35 slots, identical across upsampled/raw/native-rate/upsample-rate blocks):
 *      0: Unix epoch milliseconds (absolute)           1: ecg_1
 *      2: ecg_2                                     3: ecg_3
 *      4: ppg                                       5: accel_x
 *      6: accel_y                                   7: accel_z
 *      8: marker                                    9: temp
 *     10: pacemaker_event                          11: eog_l
 *     12: eog_r                                    13: emg
 *     14: eeg_1                                    15: eeg_2
 *     16: eeg_3                                    17: eeg_4
 *     18: cvp                                      19: pres
 *     20: flow                                     21: snore
 *     22: thor                                     23: abdo
 *     24: leg                                      25: auxac
 *     26: therm                                    27: pos
 *     28: oxstatus                                 29: spo2
 *     30: HR                                       31: DHR
 *     32: resp                                     33: aBP
 *     34: art                                      35: art_pulm
 *
 * @author Mira Welner
 * @email MEW386@pitt.edu
 * @date   2026-07-02
 */

#include <cstdint>
#include <filesystem>
#include "config_entry.hpp"   // config_entry

 // ============================================================================
 // Public constants
 // ============================================================================
inline constexpr int          NUM_CHANNELS = 36;
// On-disk header format version, written as a uint32 at offset 0 so future
// readers can detect layout changes (e.g. a larger slot count). Bump when
// the header or channel layout changes incompatibly.
inline constexpr uint32_t     BIN_HEADER_VERSION = 1;
// Scalar fields: version, n_channels, sleep_state_len, size_sleep (4),
// plus 4 per-channel arrays (upsampled sizes, raw sizes, native rates,
// upsampled rates). 4 + 4*36 = 148 fields x 4 bytes = 592 bytes.
inline constexpr int          NUM_HEADER_FIELDS = 4 + 4 * NUM_CHANNELS; // = 148
inline constexpr std::streamoff HEADER_SIZE = NUM_HEADER_FIELDS * 4;    // = 592

// Channel indices (one source of truth, used throughout file_to_bin and
// any consumer that needs to address channels by name).
enum ChannelIdx {
    CH_TIMESTAMP = 0,
    CH_ECG1, CH_ECG2, CH_ECG3, CH_PPG,
    CH_ACCEL_X, CH_ACCEL_Y, CH_ACCEL_Z,
    CH_MARKER, CH_TEMP, CH_PACEMAKER_EVENT,
    CH_EOG_L, CH_EOG_R, CH_EMG,
    CH_EEG1, CH_EEG2, CH_EEG3, CH_EEG4,
    CH_CVP, CH_PRES, CH_FLOW, CH_SNORE, CH_THOR, CH_ABDO,
    CH_LEG, CH_AUXAC, CH_THERM, CH_POS,
    CH_OXSTATUS, CH_SPO2, CH_HR, CH_DHR,
    CH_RESP, CH_ABP,
    CH_ART, CH_ART_PULM
};

// ============================================================================
// Public conversion entry points
// ============================================================================

/**
 * @brief Convert one EDF recording (with optional sleep-stage XML sidecar)
 *        into a 35-channel .bin file in cfg.output_path.
 *
 * @param path     Path to the source .edf
 * @param cfg      Dataset config (paths, rates, channel labels)
 */
void make_binfile_edf(const std::filesystem::path& path, const config_entry& cfg);

/**
 * @brief Convert a time window [winStartSec, winEndSec) of one EDF recording
 *        into a .bin at outPath. winEndSec < 0 means "to the end of the file".
 *        Used to split a recording into fixed-length segments (e.g. the
 *        8-hour Bittium segments written as <stem>_0.bin, <stem>_1.bin, ...).
 *        Each channel is sliced to the same wall-clock window on its own
 *        native grid, so the native and upsampled time axes stay identical.
 */
void make_binfile_edf_window(const std::filesystem::path& path,
    const config_entry& cfg,
    const std::filesystem::path& outPath,
    double winStartSec, double winEndSec);

/**
 * @brief Convert one CHAOS / Bittium .dat recording into a 35-channel
 *        .bin file in cfg.output_path. See make_binfile_edf for details.
 */
void make_binfile_dat(const std::filesystem::path& path, const config_entry& cfg);

/**
 * @brief Dispatch to make_binfile_edf or make_binfile_dat based on `path`'s
 *        extension. Logs an error and does nothing for unsupported types.
 */
std::filesystem::path make_binfile(const std::filesystem::path& path, const config_entry& cfg);