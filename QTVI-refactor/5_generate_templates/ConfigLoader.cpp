#include "ConfigLoader.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

ConfigLoader::ConfigLoader() {}

bool ConfigLoader::loadConfig(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string line;
    if (!std::getline(file, line)) return false; // Skip header

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> columns;
        while (std::getline(ss, cell, ',')) columns.push_back(trim(cell));

        if (columns.size() >= 9) {
            ConfigRow row;
            row.dataType = columns[0];
            row.waveBinPath = columns[7];
            row.generateTemplatePath = columns[8];
            rows.push_back(row);
        }
    }
    return !rows.empty();
}

std::vector<BinFile> ConfigLoader::getBinFiles(const std::string& directoryPath) const {
    std::vector<BinFile> fileList;
    try {
        if (fs::exists(directoryPath) && fs::is_directory(directoryPath)) {
            // Recursive search (like Matlab's **/*.bin)
            for (const auto& entry : fs::recursive_directory_iterator(directoryPath)) {
                if (entry.path().extension() == ".bin") {
                    fileList.push_back({ entry.path().filename().string(), entry.path().string() });
                }
            }
        }
    }
    catch (...) {}
    return fileList;
}

ConfigPaths ConfigLoader::getPathsFromUserSelection() const {
    ConfigPaths result;
    if (rows.empty()) return result;

    int selection = 0;
    std::cout << "Select Data Type (1-MESA, 2-Bittium, 3-CHAOS): ";
    if (std::cin >> selection && selection >= 1 && selection <= (int)rows.size()) {
        const ConfigRow& r = rows[selection - 1];
        result.sourcePath = r.waveBinPath;
        result.destinationPath = r.generateTemplatePath;
        result.dataType = r.dataType;
    }
    return result;
}

std::string ConfigLoader::trim(const std::string& s) const {
    size_t first = s.find_first_not_of(" \t\r\n");
    if (std::string::npos == first) return "";
    size_t last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, (last - first + 1));
}

// Static method to load the binary wave data
std::vector<WaveSegment> ConfigLoader::loadWaveData(const std::string& filePath) {
    std::ifstream file(filePath, std::ios::binary);
    std::vector<WaveSegment> segments;
    if (!file.is_open()) return segments;

    int32_t numSegs = 0;
    file.read(reinterpret_cast<char*>(&numSegs), sizeof(int32_t));

    for (int i = 0; i < numSegs; ++i) {
        WaveSegment seg;
        seg.segmentIndex = i;
        file.read(reinterpret_cast<char*>(&seg.ppgSamplingRate), sizeof(double));
        file.read(reinterpret_cast<char*>(&seg.ecgSamplingRate), sizeof(double));
        file.read(reinterpret_cast<char*>(&seg.isBadSegment), sizeof(bool));

        int32_t ppgN, ecgN;
        file.read(reinterpret_cast<char*>(&ppgN), sizeof(int32_t));
        seg.ppgData.resize(ppgN);
        file.read(reinterpret_cast<char*>(seg.ppgData.data()), ppgN * sizeof(double));

        file.read(reinterpret_cast<char*>(&ecgN), sizeof(int32_t));
        seg.ecgData.resize(ecgN);
        file.read(reinterpret_cast<char*>(seg.ecgData.data()), ecgN * sizeof(double));

        segments.push_back(seg);
    }
    return segments;
}
