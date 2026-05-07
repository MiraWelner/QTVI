/**
 * @file   main.cpp
 * @brief  Identify the R peaks on a ECG/PPG signal and pair them with PPG diastolic dips
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-20
 */

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

#include "binfile_handling.hpp"
#include "create_ecg_ppg_pairs.hpp"

struct ConfigSettings {
    std::string dataType;
    std::string annealedPath;
    std::string wavePath;
};

AnnealedData read_input_binfile(const std::string& path) {
    AnnealedData data;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Could not open: " + path);

    // Read buffer for better I/O performance
    char read_buf[1 << 16];
    file.rdbuf()->pubsetbuf(read_buf, sizeof(read_buf));

    uint64_t numBins = 0;
    file.read(reinterpret_cast<char*>(&numBins), 8);

    double filePpgSR, fileEcgSR, scoringEpoch;
    file.read(reinterpret_cast<char*>(&filePpgSR), 8);
    file.read(reinterpret_cast<char*>(&fileEcgSR), 8);
    file.read(reinterpret_cast<char*>(&scoringEpoch), 8);

    // Channel count + native rates: consumed but not retained. The writer
    // recovers per-segment channel count from bin.all_upsampled.size().
    uint32_t nChannels = 0;
    file.read(reinterpret_cast<char*>(&nChannels), 4);
    if (nChannels > 0)
        file.seekg(static_cast<std::streamoff>(nChannels) * 4, std::ios::cur);

    data.bins.resize(numBins);

    for (uint64_t i = 0; i < numBins; ++i) {
        auto& bin = data.bins[i];

        // PPG index pairs - bulk read
        uint64_t nPpgPairs;
        if (!file.read(reinterpret_cast<char*>(&nPpgPairs), 8)) break;
        bin.ppg_bin_indexs.resize(nPpgPairs);
        if (nPpgPairs > 0)
            file.read(reinterpret_cast<char*>(bin.ppg_bin_indexs.data()), nPpgPairs * 16);

        // ECG index pairs - bulk read
        uint64_t nEcgPairs;
        if (!file.read(reinterpret_cast<char*>(&nEcgPairs), 8)) break;
        bin.ecg_bin_indexs.resize(nEcgPairs);
        if (nEcgPairs > 0)
            file.read(reinterpret_cast<char*>(bin.ecg_bin_indexs.data()), nEcgPairs * 16);

        // Helper: read a double array preceded by its uint64 count
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

        // Per-segment 41 x {upsampled samples, raw (t, v) pairs}. Each
        // upsampled block is preceded by its sample count; each raw block
        // by its pair count (block byte length = 2 * pair_count * 8). We
        // read both into flat vectors -- keeping (t, v) interleaved
        // matches the on-disk layout exactly, so the writer can pass it
        // through with one bulk write.
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
 * @brief  Write output bin file with R-peaks, preprocessed signals,
 *         pre-bandpass signals, and raw signals from 3 methods x 3 channels,
 *         plus pass-through copies of every input channel.
 *
 * Layout per bin:
 *   - 9 index arrays (ch1 raw/sq/abs, ch2 raw/sq/abs, ch3 raw/sq/abs)
 *   - 2 PPG index arrays (maxAmps, minAmps)
 *   - 4 raw signal arrays (ppg, ecg1, ecg2, ecg3)
 *   - 6 preprocessed signal arrays (ch1 sq/abs, ch2 sq/abs, ch3 sq/abs)
 *   - 9 noise flags (1 byte each)
 *   - pairs array
 *   - ppg_bin_indexs
 *   - ecg_bin_indexs
 *   - bin.all_upsampled.size() x { upsampled doubles + raw (t, v) pair doubles }
 *     (pass-through from the annealed input; same slot order as step 3)
 *
 * All index arrays are written 1-based for MATLAB compatibility.
 */
void write_output_binfile(const std::string& path, const std::vector<output_binfile_data>& results) {
    std::ofstream file(path, std::ios::binary);
    if (!file.is_open()) return;

    char write_buf[1 << 16];
    file.rdbuf()->pubsetbuf(write_buf, sizeof(write_buf));

    uint64_t numBins = results.size();
    file.write(reinterpret_cast<const char*>(&numBins), 8);

    for (const auto& bin : results) {

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

        auto writeSignal = [&](const std::vector<double>& sig) {
            uint64_t sz = sig.size();
            file.write(reinterpret_cast<const char*>(&sz), 8);
            if (sz > 0) file.write(reinterpret_cast<const char*>(sig.data()), sz * 8);
            };

        auto writePairVec = [&](const std::vector<std::pair<uint64_t, uint64_t>>& v) {
            uint64_t sz = v.size();
            file.write(reinterpret_cast<const char*>(&sz), 8);
            if (sz > 0) file.write(reinterpret_cast<const char*>(v.data()), sz * 16);
            };

        /* R-peak indices: 3 methods x 3 channels = 9 arrays */
        writeIdx(bin.ch1.raw);
        writeIdx(bin.ch1.squared);
        writeIdx(bin.ch1.absval);

        writeIdx(bin.ch2.raw);
        writeIdx(bin.ch2.squared);
        writeIdx(bin.ch2.absval);

        writeIdx(bin.ch3.raw);
        writeIdx(bin.ch3.squared);
        writeIdx(bin.ch3.absval);

        /* PPG indices */
        writeIdx(bin.ppgMaxAmps);
        writeIdx(bin.ppgMinAmps);

        /* Raw signals (unmodified) */
        writeSignal(bin.ppgSignal);
        writeSignal(bin.ecgSignal);
        writeSignal(bin.ecgSignal2);
        writeSignal(bin.ecgSignal3);

        /* Preprocessed signals: squared and absval per channel */
        writeSignal(bin.ch1.squared_signal);
        writeSignal(bin.ch1.absval_signal);
        writeSignal(bin.ch2.squared_signal);
        writeSignal(bin.ch2.absval_signal);
        writeSignal(bin.ch3.squared_signal);
        writeSignal(bin.ch3.absval_signal);

        /* 9 noise flags: 3 methods x 3 channels */
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

        /* Pairs */
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

        writePairVec(bin.ppg_bin_indexs);
        writePairVec(bin.ecg_bin_indexs);

        /* Pass-through: per-channel { upsampled samples, raw (t, v) pairs },
           in the same slot order as step 3. The reader sized both vectors
           identically per segment, so they have the same length and that
           length serves as the segment's channel count. */
        for (size_t ch = 0; ch < bin.all_upsampled.size(); ++ch) {
            writeSignal(bin.all_upsampled[ch]);

            // Raw block: pair-count followed by 2*pair-count doubles
            // (interleaved t, v). bin.all_raw_pairs_flat[ch] is already in
            // interleaved form so just emit count + bulk-write.
            const auto& rwSrc = bin.all_raw_pairs_flat[ch];
            uint64_t nPairs = rwSrc.size() / 2;
            file.write(reinterpret_cast<const char*>(&nPairs), 8);
            if (nPairs > 0)
                file.write(reinterpret_cast<const char*>(rwSrc.data()),
                    nPairs * 2 * sizeof(double));
        }
    }
}
// ============================================================================
// Config
// ============================================================================

ConfigSettings parseConfig(const std::string& configPath, int type_row) {
    std::ifstream file(configPath);
    if (!file.is_open()) throw std::runtime_error("Could not open config.csv");

    std::string line;
    std::getline(file, line); // Skip header

    int currentIdx = 0;
    while (std::getline(file, line)) {
        if (currentIdx++ == type_row - 1) {
            std::stringstream ss(line);
            std::vector<std::string> row;
            std::string val;
            while (std::getline(ss, val, ',')) row.push_back(val);

            return { row[0], row[6], row[7] };
        }
    }
    throw std::runtime_error("Config type row not found");
}

// Processes one already-annealed .bin file. Returns true on success.
bool processWaveMarkingsForFile(const std::string& annealedBinPath,
    const std::string& wavePath,
    const std::string& currentID)
{
    try {
        AnnealedData annealedData = read_input_binfile(annealedBinPath);

        std::vector<output_binfile_data> results = create_ecg_ppg_pairs(
            std::move(annealedData.bins), 0, true, currentID);

        std::string outputPath = wavePath + "/" + currentID + "_wave_markings.bin";
        write_output_binfile(outputPath, results);
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "Wave-marking failed for " << currentID << ": " << e.what() << std::endl;
        return false;
    }
}