
#include <Eigen/Dense>
#include <iostream>
#include "ConfigLoader.h"

int main() {
    ConfigLoader loader;
    if (loader.loadConfig("config.csv")) {
        ConfigPaths paths = loader.getPathsFromUserSelection();

        // 1. Get the list of files (the C++ version of Matlab's 'dir')
        std::vector<BinFile> files = loader.getBinFiles(paths.sourcePath);

        if (!files.empty()) {
            std::cout << "Found " << files.size() << " binary files." << std::endl;

            // 2. Load data from the first file (the C++ version of Matlab's 'load')
            std::vector<WaveSegment> wave_data = ConfigLoader::loadWaveData(files[0].fullPath);

            if (!wave_data.empty()) {
                std::cout << "Loaded file: " << files[0].fileName << std::endl;
                std::cout << "Total segments: " << wave_data.size() << std::endl;
            }
        }
    }
    return 0;
}
