#pragma once
/**
 * @file   file_to_bin.hpp
 * @brief  The 'step zero' of the QTVi processing
           Take in a MESA, Bittium, or CHAOS file and convert it to a .bin file with 40 channels.
 *         Some channels will be empty depending on the file being read - for example, when loading
 *         CHAOS, sleep state will be empty.
 *         Uses polyphase upsampling with a filter bank for speed and efficiency
 *
 *   500-byte header (125 x 32-bit fields):
 *
 *     Offset   0: signal_rate           (uint32)  common upsampled rate (1000 Hz)
 *     Offset   4: boolean_rate          (uint32)  shared rate for <=1 Hz channels (1 Hz)
 *     Offset   8: pacemaker_event_rate  (uint32)  pacemaker epoch rate (8 Hz)
 *     Offset  12: sleep_state_rate      (uint32)  sleep-stage epoch length in seconds (30 s)
 *
 *     Offset  16: upsampled sizes       40 x uint32  (size_<chan>)
 *     Offset 176: raw-pair sizes        40 x uint32  (size_<chan>_raw, counts PAIRS)
 *     Offset 336: native sampling rates 40 x float32 (Hz)
 *                 0.0 = channel absent; negative values never used.
 *     Offset 496: size_sleep            (uint32)
 *
 *   Channel index order (40 slots, identical across upsampled/raw/native-rate blocks):
 *      0: seconds from start of recording           1: ecg_1
 *      2: ecg_2                                     3: ecg_3
 *      4: ppg                                       5: accel_x
 *      6: accel_y                                   7: accel_z
 *      8: marker                                    9: temp
 *     10: pacemaker                                11: eog_l
 *     12: eog_r                                    13: emg
 *     14: eeg_1                                    15: eeg_2
 *     16: eeg_3                                    17: eeg_4
 *     18: pres                                     19: flow
 *     20: thor                                     21: abdo
 *     22: leg                                      23: therm
 *     24: pos                                      25: ekg_off
 *     26: eog_l_off                                27: eog_r_off
 *     28: emg_off                                  29: eeg1_off
 *     30: eeg2_off                                 31: eeg3_off
 *     32: oxstatus                                 33: spo2
 *     34: HR                                       35: DHR
 *     36: resp                                     37: aBP
 *     38: art                                      39: art_pulm
 *
 *   Header size check: 4 + 40 + 40 + 40 + 1 = 125 fields x 4 bytes = 500 bytes
 *
 * @author Mira Welner
 * @email MEW386@pitt.edu
 * @date   2026-05-04
 */

#include <filesystem>
#include "config_entry.hpp"   // ConfigEntry

 // ============================================================================
 // Public constants
 // ============================================================================
inline constexpr int          NUM_CHANNELS = 40;
inline constexpr int          NUM_HEADER_FIELDS = 4 + 3 * NUM_CHANNELS + 1;  // = 125
inline constexpr std::streamoff HEADER_SIZE = NUM_HEADER_FIELDS * 4;     // = 500
inline constexpr double       SLEEP_STATE_LENGTH = 30.0;
inline constexpr double       BOOLEAN_RATE = 1.0;   // <=1 Hz channels are not upsampled
inline constexpr uint32_t     PACEMAKER_RATE = 8;

// Channel indices (one source of truth, used throughout file_to_bin and
// any consumer that needs to address channels by name).
enum ChannelIdx {
    CH_TIMESTAMP = 0,
    CH_ECG1, CH_ECG2, CH_ECG3, CH_PPG,
    CH_ACCEL_X, CH_ACCEL_Y, CH_ACCEL_Z,
    CH_MARKER, CH_TEMP, CH_PACEMAKER,
    CH_EOG_L, CH_EOG_R, CH_EMG,
    CH_EEG1, CH_EEG2, CH_EEG3, CH_EEG4,
    CH_PRES, CH_FLOW, CH_THOR, CH_ABDO,
    CH_LEG, CH_THERM, CH_POS,
    CH_EKG_OFF, CH_EOG_L_OFF, CH_EOG_R_OFF, CH_EMG_OFF,
    CH_EEG1_OFF, CH_EEG2_OFF, CH_EEG3_OFF,
    CH_OXSTATUS, CH_SPO2, CH_HR, CH_DHR,
    CH_RESP, CH_ABP,
    CH_ART, CH_ART_PULM
};

// ============================================================================
// Public conversion entry points
// ============================================================================

/**
 * @brief Convert one EDF recording (with optional sleep-stage XML sidecar)
 *        into a 40-channel .bin file at cfg.binFilePath.
 *
 *        Output filename is `<source-stem>_<rate>.bin` where rate is
 *        cfg.finalSamplingRate. Every channel is resampled to that rate
 *        for the upsampled block; the raw block keeps the channel's
 *        native (t, v) samples unchanged.
 *
 * @param path     Path to the source .edf
 * @param xmlPath  Path to a matching sleep-stage XML, or empty
 * @param cfg      Dataset config (paths, rates, channel labels)
 */
void make_binfile_edf(const std::filesystem::path& path,
    const std::filesystem::path& xmlPath,
    const config_entry& cfg);

/**
 * @brief Convert one CHAOS / Bittium .dat (CSV) recording into a 40-channel
 *        .bin file at cfg.binFilePath. See make_binfile_edf for details.
 */
void make_binfile_dat(const std::filesystem::path& path,
    const config_entry& cfg);

/**
 * @brief Dispatch to make_binfile_edf or make_binfile_dat based on `path`'s
 *        extension. Logs an error and does nothing for unsupported types.
 */
void make_binfile(const std::filesystem::path& path,
    const std::filesystem::path& xmlPath,
    const config_entry& cfg);

