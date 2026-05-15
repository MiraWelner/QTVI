/**
 * @file   run_find_r_peaks.hpp
 * @brief  Binary I/O for the R-peak / PPG-pairing pipeline.
 *         Provides read_input_binfile() to load an annealed-segments .bin
 *         produced upstream, and write_output_binfile() to serialize the
 *         per-segment results (R-peaks, preprocessed signals, PPG pairing)
 *         to disk. The peak detection itself lives in create_ecg_ppg_pairs().
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

#include "peakfinding_io.hpp"
#include "create_ecg_ppg_pairs.hpp"

 /**
  * @brief  Read an annealed-segments .bin file produced by the upstream
  *         preparation step.
  *
  * On-disk layout (all values little-endian, native widths):
  *
  *   Header:
  *     uint64   numBins                        // segment count
  *     double   filePpgSR                      // PPG sample rate (consumed, not retained)
  *     double   fileEcgSR                      // ECG sample rate (consumed, not retained)
  *     double   scoringEpoch                   // (consumed, not retained)
  *     uint32   nChannels                      // pass-through channel count
  *     uint32   nativeSR[nChannels]            // per-channel native rate (skipped)
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
  * The pass-through channels are not consumed by the peak-detection pipeline;
  * they are routed through to write_output_binfile() unchanged. Keeping (t,v)
  * interleaved on disk lets the writer pass them through with a single bulk
  * write per slot.
  */
inline AnnealedData read_input_binfile(const std::string& path) {
    AnnealedData data;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Could not open: " + path);

    // Buffered I/O for throughput.
    char read_buf[1 << 16];
    file.rdbuf()->pubsetbuf(read_buf, sizeof(read_buf));

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
 * On-disk layout (little-endian, native widths):
 *
 *   uint64   numBins
 *
 *   For each bin (in order):
 *
 *     // R-peak indices: 3 channels x 3 preprocessing methods = 9 arrays.
 *     // Each: uint64 count + count uint64 indices (1-based, see note below).
 *     ch1.raw, ch1.squared, ch1.absval
 *     ch2.raw, ch2.squared, ch2.absval
 *     ch3.raw, ch3.squared, ch3.absval
 *
 *     // PPG event indices (same uint64 + 1-based layout as above):
 *     ppgMaxAmps
 *     ppgMinAmps
 *
 *     // Raw signals (uint64 count + count doubles, unmodified):
 *     ppgSignal, ecgSignal, ecgSignal2, ecgSignal3
 *
 *     // Preprocessed signals (uint64 count + count doubles):
 *     ch1.squared_signal, ch1.absval_signal
 *     ch2.squared_signal, ch2.absval_signal
 *     ch3.squared_signal, ch3.absval_signal
 *
 *     // 9 noise flags, one byte each, in the same channel x method order
 *     // as the index arrays above:
 *     uint8    flags[9]
 *
 *     // Pairs (PPG-valley index, ECG-R-peak index):
 *     uint64   numPairs
 *     int64    pairBuf[2 * numPairs]    // interleaved (ppg, ecg); 1-based with -1 sentinel
 *
 *     // PPG / ECG bin-index ranges (uint64 pairs, written raw / 0-based):
 *     uint64   nPpgBinIdx;  (uint64,uint64) ppg_bin_indexs[nPpgBinIdx]
 *     uint64   nEcgBinIdx;  (uint64,uint64) ecg_bin_indexs[nEcgBinIdx]
 *
 *     // Pass-through channels: bin.all_upsampled.size() slots, same order
 *     // as the input file. Each slot:
 *     //   uint64 nUp;     double upsampled[nUp]
 *     //   uint64 nPairs;  double raw_tv[2 * nPairs]    // interleaved (t,v)
 *     // The slot count is recovered from bin.all_upsampled.size() at read
 *     // time -- it is NOT written into the file as a per-segment header.
 *
 * Index 1-basing summary:
 *   1-based (writer adds 1):  R-peak arrays, ppgMaxAmps, ppgMinAmps, pairBuf entries
 *   0-based (raw passthrough): ppg_bin_indexs, ecg_bin_indexs, nUp/nPairs counts,
 *                               and the raw (t, v) doubles in pass-through slots
 *
 * The -1 sentinel in pairBuf marks an unpaired side; NaN or negative
 * doubles in bin.pairs are mapped to -1 on write.
 */
inline void write_output_binfile(const std::string& path, const std::vector<output_binfile_data>& results) {
    std::ofstream file(path, std::ios::binary);
    if (!file.is_open()) return;

    char write_buf[1 << 16];
    file.rdbuf()->pubsetbuf(write_buf, sizeof(write_buf));

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

        // Write a (uint64 count, count (uint64,uint64) pairs) block, raw (no offset).
        auto writePairVec = [&](const std::vector<std::pair<uint64_t, uint64_t>>& v) {
            uint64_t sz = v.size();
            file.write(reinterpret_cast<const char*>(&sz), 8);
            if (sz > 0) file.write(reinterpret_cast<const char*>(v.data()), sz * 16);
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

        /* Raw signals (unmodified). */
        writeSignal(bin.ppgSignal);
        writeSignal(bin.ecgSignal);
        writeSignal(bin.ecgSignal2);
        writeSignal(bin.ecgSignal3);

        /* Preprocessed signals: squared then absval, per channel. */
        writeSignal(bin.ch1.squared_signal);
        writeSignal(bin.ch1.absval_signal);
        writeSignal(bin.ch2.squared_signal);
        writeSignal(bin.ch2.absval_signal);
        writeSignal(bin.ch3.squared_signal);
        writeSignal(bin.ch3.absval_signal);

        /* 9 noise flags: 3 channels x 3 methods (raw/squared/absval). */
        uint8_t flags[9] = {
            static_cast<uint8_t>(bin.ch1.raw_noisy),
            static_cast<uint8_t>(bin.ch1.squared_noisy),
            static_cast<uint8_t>(bin.ch1.absval_noisy),
            static_cast<uint8_t>(bin.ch2.raw_noisy),
            static_cast<uint8_t>(bin.ch2.squared_noisy),
            static_cast<uint8_t>(bin.ch2.absval_noisy),
            static_cast<uint8_t>(bin.ch3.raw_noisy),
            static_cast<uint8_t>(bin.ch3.squared_noisy),
            static_cast<uint8_t>(bin.ch3.absval_noisy),
        };
        file.write(reinterpret_cast<const char*>(flags), 9);

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

        /* Bin-index range arrays: written raw as uint64 pairs (NOT 1-based). */
        writePairVec(bin.ppg_bin_indexs);
        writePairVec(bin.ecg_bin_indexs);

        /* Pass-through: per-slot { uint64 nUp + nUp doubles, uint64 nPairs +
           2*nPairs interleaved (t,v) doubles }. Slot order matches the
           input file. No per-segment slot-count header is written -- the
           reader sizes both vectors identically and recovers the count
           from bin.all_upsampled.size(). */
        for (size_t ch = 0; ch < bin.all_upsampled.size(); ++ch) {
            writeSignal(bin.all_upsampled[ch]);

            // Raw (t,v) block: uint64 pair-count, then 2*pair-count doubles.
            // all_raw_pairs_flat[ch] is already interleaved on disk-layout,
            // so we emit count + a single bulk-write.
            const auto& rwSrc = bin.all_raw_pairs_flat[ch];
            uint64_t nPairs = rwSrc.size() / 2;
            file.write(reinterpret_cast<const char*>(&nPairs), 8);
            if (nPairs > 0)
                file.write(reinterpret_cast<const char*>(rwSrc.data()),
                    nPairs * 2 * sizeof(double));
        }
    }
}