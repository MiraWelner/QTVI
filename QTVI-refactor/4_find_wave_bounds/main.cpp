#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <filesystem>
#include <algorithm>
#include <stdexcept>

// This include brings in the structures automatically
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

        if (!file.read(reinterpret_cast<char*>(&sS), 8)) break;
        if (sS > 0) file.seekg(sS * 8, std::ios_base::cur); // Skip broken sleep data
    }
    return data;
}

// --- Binary Saving Logic ---
void saveWaveData(const std::string& path, const std::vector<WaveData>& results) {
    std::ofstream file(path, std::ios::binary);
    if (!file.is_open()) return;

    uint64_t numBins = results.size();
    file.write(reinterpret_cast<char*>(&numBins), 8);

    for (const auto& bin : results) {
        // Save the raw R-peaks count
        uint64_t numRPeaks = bin.ecgRIndex.size();
        file.write(reinterpret_cast<char*>(&numRPeaks), 8);

        // Save each R-peak index as a double (to match MATLAB's precision)
        for (size_t rIdx : bin.ecgRIndex) {
            double val = static_cast<double>(rIdx);
            file.write(reinterpret_cast<char*>(&val), 8);
        }
    }
}

// --- Config Parser ---
ConfigSettings parseConfig(const std::string& configPath, const int type_row) {
    std::ifstream file(configPath);
    std::string line; std::getline(file, line);
    int currentIdx = 0; // Track the current data row index
    while (std::getline(file, line)) {
        // If the current index matches the requested type_row (0-indexed)
        if (currentIdx++ == type_row-1) {
            std::stringstream ss(line);
            std::vector<std::string> row;
            std::string val;
            while (std::getline(ss, val, ','))
                row.push_back(val);

            return { row[0], row[6], row[7], std::stod(row[12]),
                   (row.size() > 13 && !row[13].empty()) ? std::stod(row[13]) : std::stod(row[12]) };
        }
    }

    throw std::runtime_error("Type not found");
}

int main() {
    try {
        int choice;
        std::cout << "Enter Data Type:\nMESA: 1\nBittium: 2\nCHAOS: 3\n"; std::cin >> choice;
        ConfigSettings cfg = parseConfig("config.csv", choice);

        if (!fs::exists(cfg.wavePath)) fs::create_directories(cfg.wavePath);

        for (const auto& entry : fs::directory_iterator(cfg.annealedPath)) {
            if (entry.path().extension() == ".bin") {

                AnnealedData annealedData = readCppBin(entry.path().string(), cfg.ecgFs, cfg.ppgFs);

                // Process logic
                auto results = FindWaveBounds_EKGandPPG(annealedData.bins, cfg.ecgFs, true);

                saveWaveData(cfg.wavePath + "/" + entry.path().stem().string() + "_wave_data.bin", results);
                std::cout << "Processed: " << entry.path().filename() << std::endl;
            }
        }
    }
    catch (const std::exception& e) { std::cerr << e.what() << std::endl; }
    return 0;
}
