#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <matio.h> // Make sure matio library is linked
#include <cmath>
#include <algorithm> // For std::min

// Helper function to print bytes (for debugging)
void printBytes(const char* data, size_t size, const std::string& label) {
    std::cerr << "DEBUG: " << label << " raw bytes (hex): ";
    for (size_t i = 0; i < size; ++i) {
        // Cast to unsigned char before int to ensure proper hex printing of byte values 0-255
        std::cerr << std::hex << (static_cast<unsigned int>(static_cast<unsigned char>(data[i]))) << " ";
    }
    std::cerr << std::dec << std::endl; // Switch back to decimal for subsequent prints
}

// 1. Reads the eIdx (R-peaks) from the existing 'pairs' in the BIN file
std::vector<std::vector<double>> readRFromBinPairs(const std::string& path) {
    std::vector<std::vector<double>> binR;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return binR;

    uint64_t numBins = 0;
    file.read(reinterpret_cast<char*>(&numBins), 8);
    binR.resize(numBins);

    for (uint64_t i = 0; i < numBins; ++i) {
        uint64_t numRPeaks = 0;
        file.read(reinterpret_cast<char*>(&numRPeaks), 8);
        for (uint64_t j = 0; j < numRPeaks; ++j) {
            double val;
            file.read(reinterpret_cast<char*>(&val), 8); // Read only one double
            binR[i].push_back(val);
        }
    }
    return binR;
}

// 2. Reads the ecgRIndex from the MAT file
std::vector<std::vector<double>> readRFromMat(const std::string& path) {
    std::vector<std::vector<double>> matR;
    mat_t* matfp = Mat_Open(path.c_str(), MAT_ACC_RDONLY);
    if (!matfp) {
        std::cerr << "Error: Could not open MAT file: " << path << std::endl;
        return matR;
    }
    matvar_t* wave_data = Mat_VarRead(matfp, "wave_data");
    if (!wave_data) {
        std::cerr << "Error: Could not read 'wave_data' from MAT file: " << path << std::endl;
        Mat_Close(matfp);
        return matR;
    }

    size_t numBins = wave_data->dims[0] * wave_data->dims[1];
    matR.resize(numBins);

    for (size_t i = 0; i < numBins; ++i) {
        matvar_t* cell = Mat_VarGetCell(wave_data, i);
        if (!cell) {
            std::cerr << "Warning: Could not get cell " << i << " from 'wave_data' in MAT file." << std::endl;
            continue;
        }
        matvar_t* r_var = Mat_VarGetStructFieldByName(cell, "ecgRIndex", 0);
        if (r_var && r_var->data) {
            size_t n = 1;
            if (r_var->rank >= 1) n = r_var->dims[0];
            if (r_var->rank >= 2) n *= r_var->dims[1];

            double* data = (double*)r_var->data;
            for (size_t j = 0; j < n; ++j) matR[i].push_back(data[j]);
        }
    }
    Mat_VarFree(wave_data);
    Mat_Close(matfp);
    return matR;
}

int main(int argc, char** argv) {

    // Use raw string literals for paths for better readability and to avoid backslash escaping issues
    std::string binPath = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\cpp\mesa_files\3010155_20110511_annealed_wave_data.bin)";
    std::string matPath = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\matlab\3010155_20110511_wave_data.mat)";

    auto binR = readRFromBinPairs(binPath);
    auto matR = readRFromMat(matPath);

    if (binR.empty() && matR.empty()) {
        std::cerr << "Warning: No data loaded from either BIN or MAT files." << std::endl;
        return 1;
    }

    size_t compare_limit = 0;
    if (!binR.empty() || !matR.empty()) {
        compare_limit = std::min(binR.size(), matR.size());
        if (compare_limit == 0 && (!binR.empty() || !matR.empty())) {
            // If one is empty and the other isn't, compare up to 5 bins of the non-empty one
            compare_limit = (binR.empty() ? matR.size() : binR.size());
            compare_limit = std::min(compare_limit, (size_t)5);
        }
        else {
            compare_limit = std::min(compare_limit, (size_t)5); // Compare first 5, or fewer if less available
        }
    }


    for (size_t i = 0; i < compare_limit; ++i) {
        std::cout << "--- Bin " << i << " ---" << std::endl;

        size_t bin_count = (i < binR.size() ? binR[i].size() : 0);
        size_t mat_count = (i < matR.size() ? matR[i].size() : 0);

        std::cout << "  BIN count: " << bin_count << " | MAT count: " << mat_count << std::endl;

        if (bin_count > 0 && binR[i][0] == -1) {
            std::cout << "  [STATUS] C++ found peaks but failed to pair them (values are -1)." << std::endl;
        }
        else if (bin_count > 0 && mat_count > 0) {
            std::cout << "  First Peak: BIN=" << binR[i][0] << ", MAT=" << matR[i][0] << std::endl;
        }
        else if (bin_count == 0 && mat_count > 0) {
            std::cout << "  [STATUS] BIN has no peaks for this bin, MAT has " << mat_count << " peaks." << std::endl;
        }
        else if (bin_count > 0 && mat_count == 0) {
            std::cout << "  [STATUS] BIN has " << bin_count << " peaks, MAT has no peaks for this bin." << std::endl;
        }
        else {
            std::cout << "  [STATUS] Both BIN and MAT have no peaks for this bin." << std::endl;
        }
    }
    return 0;
}
