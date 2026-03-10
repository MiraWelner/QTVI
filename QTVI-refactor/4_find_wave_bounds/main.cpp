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

// --- Binary Loading Logic ---
AnnealedData readCppBin(const std::string& path, double ecgFs, double ppgFs) {
    AnnealedData data;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Could not open: " + path);

    uint64_t numBins = 0;
    if (!file.read(reinterpret_cast<char*>(&numBins), 8)) return data;
    data.bins.resize(numBins);

    for (uint64_t i = 0; i < numBins; ++i) {
        uint64_t pS, eS, sS;

        if (!file.read(reinterpret_cast<char*>(&pS), 8)) break;
        data.bins[i].po.resize(pS);
        file.read(reinterpret_cast<char*>(data.bins[i].po.data()), pS * 8);

        if (!file.read(reinterpret_cast<char*>(&eS), 8)) break;
        data.bins[i].ecg.resize(eS);
        file.read(reinterpret_cast<char*>(data.bins[i].ecg.data()), eS * 8);

        data.bins[i].ecgSampleRate = ecgFs;
        data.bins[i].ppgSampleRate = ppgFs;

        // Skip sleep data as it's not used in this step
        if (!file.read(reinterpret_cast<char*>(&sS), 8)) break;
        if (sS > 0) file.seekg(sS * 8, std::ios::cur);
    }
    return data;
}

// --- Binary Saving Logic ---
void saveWaveData(const std::string& path, const std::vector<WaveData>& results) {
    std::ofstream file(path, std::ios::binary);
    if (!file.is_open()) return;

    uint64_t numBins = results.size();
    file.write(reinterpret_cast<const char*>(&numBins), 8);

    for (const auto& bin : results) {
        file.write(reinterpret_cast<const char*>(&bin.ecgSamplingRate), 8);
        file.write(reinterpret_cast<const char*>(&bin.ppgSamplingRate), 8);

        auto writeIdx = [&](const std::vector<std::size_t>& v) {
            uint64_t sz = v.size();
            file.write(reinterpret_cast<const char*>(&sz), 8);
            for (std::size_t val : v) {
                uint64_t out = static_cast<uint64_t>(val) + 1; // MATLAB 1-based
                file.write(reinterpret_cast<const char*>(&out), 8);
            }
            };

        writeIdx(bin.ecgRIndex);
        writeIdx(bin.ppgMaxAmps);
        writeIdx(bin.ppgMinAmps);

        auto writeSignal = [&](const std::vector<double>& sig) {
            uint64_t sz = sig.size();
            file.write(reinterpret_cast<const char*>(&sz), 8);
            if (sz > 0) file.write(reinterpret_cast<const char*>(sig.data()), sz * 8);
            };

        writeSignal(bin.ppgSignal);
        writeSignal(bin.ecgSignal);

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
                AnnealedData annealedData = readCppBin(entry.path().string(), cfg.ecgFs, cfg.ppgFs);

                if (annealedData.bins.size() > 501) {
                    auto& b = annealedData.bins[501];
                    std::cout << "Bin 501 ECG size: " << b.ecg.size()
                        << " min: " << *std::min_element(b.ecg.begin(), b.ecg.end())
                        << " max: " << *std::max_element(b.ecg.begin(), b.ecg.end())
                        << std::endl;
                }

                // 2. Process the bins directly (DO NOT consolidate into one giant signal)
                // This matches the MATLAB 'for i = 1:length(annealedSegments)' logic
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
