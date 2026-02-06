#include <iostream>
#include <string>
#include <vector>
#include <filesystem> // Requires C++17
#include <map>

// Include your previously converted headers
#include "readProps.h"
#include "FindWaveBounds.h"

namespace fs = std::filesystem;

/**
 * Main function to orchestrate the processing of ECG/PPG data files.
 * Replicates the behavior of FindWaveBoundRunner.m
 */
int main(int argc, char* argv[]) {
    std::string configFilePath = "config.txt";
    
    // Allow overriding config file via command line
    if (argc > 1) {
        configFilePath = argv[1];
    }

    try {
        std::cout << "--- Starting Wave Boundary Detection Runner ---" << std::endl;

        // 1. Read configuration properties
        std::map<std::string, std::string> props = readProps(configFilePath);
        
        std::string inputPathStr = props["FWB_input_path"];
        std::string outputPathStr = props["FWB_output_path"];
        bool skipExisting = (props["Skip_Existing"] == "1");

        fs::path inputPath(inputPathStr);
        fs::path outputPath(outputPathStr);

        // 2. Validate input directory
        if (!fs::exists(inputPath) || !fs::is_directory(inputPath)) {
            std::cerr << "[Error] Input directory does not exist: " << inputPathStr << std::endl;
            return 1;
        }

        // 3. Create output directory if it doesn't exist (Matlab: mkdir)
        if (!fs::exists(outputPath)) {
            std::cout << "[Info] Creating output directory: " << outputPathStr << std::endl;
            fs::create_directories(outputPath);
        }

        // 4. Iterate through directory (Matlab: files = dir(fullfile(inputPath, '*_annealedSegments.mat')))
        int processedCount = 0;
        int skippedCount = 0;
        int errorCount = 0;

        std::cout << "[Info] Scanning directory: " << inputPathStr << std::endl;

        for (const auto& entry : fs::directory_iterator(inputPath)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().filename().string();

                // Filter for files containing "_annealedSegments" (equivalent to Matlab pattern matching)
                if (filename.find("_annealedSegments") != std::string::npos) {
                    
                    // Construct output filename
                    // Example: Subject01_annealedSegments.mat -> Subject01_wave_data
                    std::string baseName = filename.substr(0, filename.find("_annealedSegments"));
                    std::string outFileName = baseName + "_wave_data"; 
                    fs::path fullOutPath = outputPath / outFileName;

                    // 5. Check if we should skip existing files
                    if (skipExisting && fs::exists(fullOutPath)) {
                        std::cout << "[Skip] Already exists: " << filename << std::endl;
                        skippedCount++;
                        continue;
                    }

                    std::cout << "[Process] Processing file: " << filename << std::endl;

                    // 6. Call the processing function
                    int result = FindWaveBounds(entry.path().string(), outputPath.string());

                    if (result == 1) {
                        processedCount++;
                    } else {
                        std::cerr << "[Fail] Processing failed for: " << filename << std::endl;
                        errorCount++;
                    }
                }
            }
        }

        // Summary report
        std::cout << "\n--- Processing Summary ---" << std::endl;
        std::cout << "Successfully Processed: " << processedCount << std::endl;
        std::cout << "Skipped (Existing):     " << skippedCount << std::endl;
        std::cout << "Errors Encountered:     " << errorCount << std::endl;
        std::cout << "--------------------------" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "[Critical Error] " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
