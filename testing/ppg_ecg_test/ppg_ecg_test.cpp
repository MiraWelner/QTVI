#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdint>

/**
 * This script reads the wave_data.bin file and exports peak locations
 * to a CSV file for comparison.
 *
 * CSV Format: SegmentIndex, PeakType, PeakLocation
 */

int main() {
    std::string binPath = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\cpp\mesa_files\3010155_20110511_annealed_wave_data.bin)"; // Update this to your filename
    std::string csvPath = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\cpp\mesa_files\3010155_20110511_peaks_output.csv)";


    std::ifstream binFile(binPath, std::ios::binary);
    if (!binFile.is_open()) {
        std::cerr << "Error: Could not open binary file: " << binPath << std::endl;
        return 1;
    }

    std::ofstream csvFile(csvPath);
    if (!csvFile.is_open()) {
        std::cerr << "Error: Could not create CSV file: " << csvPath << std::endl;
        return 1;
    }

    // Write CSV Header
    csvFile << "SegmentIndex,ECG_R_Peak_1Based,PPG_Systolic_Peak_1Based\n";

    uint64_t numBins = 0;
    if (!binFile.read(reinterpret_cast<char*>(&numBins), 8)) {
        std::cerr << "Error: Failed to read file header." << std::endl;
        return 1;
    }

    std::cout << "Processing " << numBins << " segments..." << std::endl;

    for (uint64_t i = 0; i < numBins; ++i) {

        // 1. Load ECG R-peaks for this segment into a vector
        uint64_t numRPeaks = 0;
        binFile.read(reinterpret_cast<char*>(&numRPeaks), 8);
        std::vector<double> rPeaks(numRPeaks);
        for (uint64_t j = 0; j < numRPeaks; ++j) {
            binFile.read(reinterpret_cast<char*>(&rPeaks[j]), 8);
        }

        // 2. Load PPG Systolic Peaks for this segment into a vector
        uint64_t numPPGPeaks = 0;
        binFile.read(reinterpret_cast<char*>(&numPPGPeaks), 8);
        std::vector<double> ppgPeaks(numPPGPeaks);
        for (uint64_t j = 0; j < numPPGPeaks; ++j) {
            binFile.read(reinterpret_cast<char*>(&ppgPeaks[j]), 8);
        }

        // 3. Determine how many rows we need for this segment
        size_t maxRows = std::max(rPeaks.size(), ppgPeaks.size());

        // 4. Write to CSV
        for (size_t row = 0; row < maxRows; ++row) {
            csvFile << i << ","; // Column 1: Segment Index

            // Column 2: ECG R Peak (+1 for MATLAB compatibility)
            if (row < rPeaks.size()) {
                csvFile << (rPeaks[row] + 1);
            }
            csvFile << ",";

            // Column 3: PPG Systolic Peak (+1 for MATLAB compatibility)
            if (row < ppgPeaks.size()) {
                csvFile << (ppgPeaks[row] + 1);
            }
            csvFile << "\n";
        }
    }

    std::cout << "Success! Data saved to: " << csvPath << std::endl;

    binFile.close();
    csvFile.close();

    return 0;
}