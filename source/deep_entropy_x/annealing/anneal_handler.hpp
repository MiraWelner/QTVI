#pragma once
/**
 * @file   anneal_handler.hpp
 * @brief  Public interface to step 3 of the QTVi pipeline. Takes a v2
 *         data .bin (from file_to_bin) plus a noise-markings .bin (from
 *         the marking GUI), runs AnnealSegments, writes the output .bin
 *         with all 40 channels preserved.
 *
 *         The marking GUI's "Process Output" button calls
 *         annealAllForDataset(cfg), which walks every .bin in
 *         cfg.binFilePath, finds the matching noise file in
 *         cfg.noiseDataPath, and annealOneFile()s each pair.
 */

#include <filesystem>
#include "config_file_handling/config_entry.hpp"

 /**
  * @brief Anneal one recording. Reads `binPath` and (optionally) `noisePath`,
  *        runs AnnealSegments, writes to `outPath`. If `noisePath` doesn't
  *        exist, runs with empty noise markings (no exclusions). Returns
  *        true on success.
  *
  * @param ecg1_inverted  "Inverted Lead?" checkbox state for CH1
  *                       (noise_marking_gui.ui: ecg_1_reverse), written into
  *                       the output .bin header so downstream peak detection
  *                       (create_ecg_ppg_pairs.hpp) knows the true R-peak is
  *                       a local min rather than a local max for this
  *                       channel. Not auto-detected -- caller must supply
  *                       the actual checkbox value.
  * @param ecg2_inverted  Same, for CH2.
  * @param ecg3_inverted  Same, for CH3.
  */
bool annealOneFile(const std::filesystem::path& binPath,
    const std::filesystem::path& noisePath,
    const std::filesystem::path& outPath,
    double binLengthMin,
    bool ecg1_inverted, bool ecg2_inverted, bool ecg3_inverted);