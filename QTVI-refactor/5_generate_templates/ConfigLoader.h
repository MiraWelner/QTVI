#ifndef CONFIG_LOADER_H
#define CONFIG_LOADER_H

#include <string>
#include <vector>
#include <filesystem>
#include <Eigen/Dense>

namespace fs = std::filesystem;

// Metadata for a file found on disk
struct BinFile {
    std::string fileName;
    std::string fullPath;
};

// Result of user selection
struct ConfigPaths {
    std::string sourcePath;
    std::string destinationPath;
    std::string dataType;
};

// Structure for a single segment (Matlab wave_data{i})
struct WaveSegment {
    int segmentIndex;
    double ppgSamplingRate;
    double ecgSamplingRate;
    bool isBadSegment;
    Eigen::VectorXd ppgData;
    Eigen::VectorXd ecgData;
};

// Internal row from config.csv
struct ConfigRow {
    std::string dataType;
    std::string mainFileExt;
    std::string sleepFileExt;
    std::string originalFilePath;
    std::string binFilePath;
    std::string noiseDataPath;
    std::string annealedBinPath;
    std::string waveBinPath;
    std::string generateTemplatePath;
};

class ConfigLoader {
public:
    ConfigLoader();

    bool loadConfig(const std::string& filename = "config.csv");
    ConfigPaths getPathsFromUserSelection() const;

    // NEW: Finds .bin files in a directory (like Matlab's dir command)
    std::vector<BinFile> getBinFiles(const std::string& directoryPath) const;

    // NEW: Loads structured binary data (like Matlab's load command)
    static std::vector<WaveSegment> loadWaveData(const std::string& filePath);

private:
    std::vector<ConfigRow> rows;
    std::string trim(const std::string& s) const;
};

#endif
