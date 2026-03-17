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

#include "FindWaveBounds.h"
#include "FindWaveBounds_EKGandPPG.h"

namespace fs = std::filesystem;

struct ConfigSettings {
    std::string dataType;
    std::string annealedPath;
    std::string wavePath;
    double ecgFs;
    double ppgFs;
};

AnnealedData readCppBin(const std::string& path) {
    AnnealedData data;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Could not open: " + path);

    // --- Header ---
    uint64_t numBins = 0;
    file.read(reinterpret_cast<char*>(&numBins), 8);

    double filePpgSR, fileEcgSR, scoringEpoch;
    file.read(reinterpret_cast<char*>(&filePpgSR), 8);
    file.read(reinterpret_cast<char*>(&fileEcgSR), 8);
    file.read(reinterpret_cast<char*>(&scoringEpoch), 8);

    data.bins.resize(numBins);

    for (uint64_t i = 0; i < numBins; ++i) {
        // --- PPG index pairs ---
        uint64_t nPpgPairs;
        if (!file.read(reinterpret_cast<char*>(&nPpgPairs), 8)) break;
        data.bins[i].ppg_bin_indexs.resize(nPpgPairs);
        for (uint64_t j = 0; j < nPpgPairs; ++j) {
            file.read(reinterpret_cast<char*>(&data.bins[i].ppg_bin_indexs[j].first), 8);
            file.read(reinterpret_cast<char*>(&data.bins[i].ppg_bin_indexs[j].second), 8);
        }

        // --- ECG index pairs ---
        uint64_t nEcgPairs;
        if (!file.read(reinterpret_cast<char*>(&nEcgPairs), 8)) break;
        data.bins[i].ecg_bin_indexs.resize(nEcgPairs);
        for (uint64_t j = 0; j < nEcgPairs; ++j) {
            file.read(reinterpret_cast<char*>(&data.bins[i].ecg_bin_indexs[j].first), 8);
            file.read(reinterpret_cast<char*>(&data.bins[i].ecg_bin_indexs[j].second), 8);
        }

        // --- PPG samples ---
        uint64_t pS;
        if (!file.read(reinterpret_cast<char*>(&pS), 8)) break;
        data.bins[i].po.resize(pS);
        if (pS > 0) file.read(reinterpret_cast<char*>(data.bins[i].po.data()), pS * 8);

        // --- ECG channel 1 ---
        uint64_t e1S;
        if (!file.read(reinterpret_cast<char*>(&e1S), 8)) break;
        data.bins[i].ecg.resize(e1S);
        if (e1S > 0) file.read(reinterpret_cast<char*>(data.bins[i].ecg.data()), e1S * 8);

        // --- ECG channel 2 ---
        uint64_t e2S;
        if (!file.read(reinterpret_cast<char*>(&e2S), 8)) break;
        data.bins[i].ecg2.resize(e2S);
        if (e2S > 0) file.read(reinterpret_cast<char*>(data.bins[i].ecg2.data()), e2S * 8);

        // --- ECG channel 3 ---
        uint64_t e3S;
        if (!file.read(reinterpret_cast<char*>(&e3S), 8)) break;
        data.bins[i].ecg3.resize(e3S);
        if (e3S > 0) file.read(reinterpret_cast<char*>(data.bins[i].ecg3.data()), e3S * 8);

        // --- Sleep stages ---
        uint64_t sS;
        if (!file.read(reinterpret_cast<char*>(&sS), 8)) break;
        data.bins[i].sleepStages.resize(sS);
        if (sS > 0) file.read(reinterpret_cast<char*>(data.bins[i].sleepStages.data()), sS * 8);

        data.bins[i].ecgSampleRate = fileEcgSR;
        data.bins[i].ppgSampleRate = filePpgSR;
    }
    return data;
}

// --- Binary Saving Logic (updated for 3 ECG channels) ---
void saveWaveData(const std::string& path, const std::vector<WaveData>& results) {
    std::ofstream file(path, std::ios::binary);
    if (!file.is_open()) return;

    uint64_t numBins = results.size();
    file.write(reinterpret_cast<const char*>(&numBins), 8);

    for (const auto& bin : results) {
        file.write(reinterpret_cast<const char*>(&bin.ecgSamplingRate), 8);
        file.write(reinterpret_cast<const char*>(&bin.ppgSamplingRate), 8);

        // Lambda: write a size_t index array (converted to 1-based for MATLAB)
        auto writeIdx = [&](const std::vector<std::size_t>& v) {
            uint64_t sz = v.size();
            file.write(reinterpret_cast<const char*>(&sz), 8);
            for (std::size_t val : v) {
                uint64_t out = static_cast<uint64_t>(val) + 1; // MATLAB 1-based
                file.write(reinterpret_cast<const char*>(&out), 8);
            }
            };

        // Lambda: write a double signal array
        auto writeSignal = [&](const std::vector<double>& sig) {
            uint64_t sz = sig.size();
            file.write(reinterpret_cast<const char*>(&sz), 8);
            if (sz > 0) file.write(reinterpret_cast<const char*>(sig.data()), sz * 8);
            };

        // --- Channel 1 R-indices ---
        writeIdx(bin.ecgRIndex);
        // --- Channel 2 R-indices ---
        writeIdx(bin.ecgRIndex2);
        // --- Channel 3 R-indices ---
        writeIdx(bin.ecgRIndex3);

        // --- PPG indices (same as before) ---
        writeIdx(bin.ppgMaxAmps);
        writeIdx(bin.ppgMinAmps);

        // --- PPG signal ---
        writeSignal(bin.ppgSignal);

        // --- ECG channel 1 signal ---
        writeSignal(bin.ecgSignal);
        // --- ECG channel 2 signal ---
        writeSignal(bin.ecgSignal2);
        // --- ECG channel 3 signal ---
        writeSignal(bin.ecgSignal3);

        // --- Pairs (channel 1 pairing, same as before) ---
        uint64_t numPairs = bin.pairs.size();
        file.write(reinterpret_cast<const char*>(&numPairs), 8);
        for (const auto& p : bin.pairs) {
            int64_t ppg_out = (std::isnan(p[0]) || p[0] < -0.1) ? -1 : static_cast<int64_t>(std::round(p[0])) + 1;
            int64_t ecg_out = (std::isnan(p[1]) || p[1] < -0.1) ? -1 : static_cast<int64_t>(std::round(p[1])) + 1;
            file.write(reinterpret_cast<char*>(&ppg_out), 8);
            file.write(reinterpret_cast<char*>(&ecg_out), 8);
        }
    }
}

// --- Config Parser ---
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

            return { row[0], row[6], row[7], std::stod(row[12]),
                   (row.size() > 13 && !row[13].empty()) ? std::stod(row[13]) : std::stod(row[12]) };
        }
    }
    throw std::runtime_error("Config type row not found");
}

int main() {
    try {
        int choice;
        std::cout << "Enter Data Type:\n1) MESA\n2) Bittium\n3) CHAOS" << std::endl;
        if (!(std::cin >> choice)) return 1;

        ConfigSettings cfg = parseConfig("config.csv", choice);

        // Standard C++17 filesystem call
        if (!fs::exists(cfg.wavePath)) {
            fs::create_directories(cfg.wavePath);
        }

        // --- SINGLE LOOP: Process each file exactly once ---
        for (const auto& entry : fs::directory_iterator(cfg.annealedPath)) {
            if (entry.path().extension() == ".bin") {
                std::string currentID = entry.path().stem().string();
                std::cout << "Processing: " << currentID << "..." << std::endl;

                // 1. Read the binary file into separate 30s bins
                AnnealedData annealedData = readCppBin(entry.path().string());

                // 2. Process the bins
                std::string detrendedPath = R"(D:\USERS\MiraWelner\QTVI\QTVI-refactor\peakfind_output\)" + currentID + "_detrended.csv";
                if (std::filesystem::exists(detrendedPath)) {
                    std::filesystem::remove(detrendedPath);
                }
                std::vector<WaveData> results = FindWaveBounds_EKGandPPG(annealedData.bins, 0, true, currentID);

                // 3. Save the results to the output path
                std::string outputPath = cfg.wavePath + "/" + currentID + "_wave_data.bin";
                saveWaveData(outputPath, results);
            }
        }
        std::cout << "Processing Complete." << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "CRITICAL ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}