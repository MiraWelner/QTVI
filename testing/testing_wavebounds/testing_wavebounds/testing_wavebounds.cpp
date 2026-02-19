#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <matio.h> 
#include <cmath>
#include <algorithm>
#include <cstdint> // Added for int64_t

struct ComparisonData {
    std::vector<double> ecgRIndex;
    std::vector<double> ppgMaxAmps;
    std::vector<double> ppgMinAmps;
    std::vector<std::pair<double, double>> pairs;
};

// Helper to read a vector of doubles from binary
void readVecFromBin(std::ifstream& file, std::vector<double>& vec) {
    uint64_t count = 0;
    if (!file.read(reinterpret_cast<char*>(&count), 8)) return;
    vec.resize(count);
    for (uint64_t i = 0; i < count; ++i) {
        int64_t val; // CHANGED: Read as signed int64_t to correctly handle -1
        file.read(reinterpret_cast<char*>(&val), 8);
        vec[i] = static_cast<double>(val);
    }
}

// 1. Reads ALL sections from the BIN file
std::vector<ComparisonData> readAllFromBin(const std::string& path) {
    std::vector<ComparisonData> results;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return results;

    uint64_t numBins = 0;
    file.read(reinterpret_cast<char*>(&numBins), 8);
    results.resize(numBins);

    for (uint64_t i = 0; i < numBins; ++i) {
        // Section 1: ECG R-Peaks
        readVecFromBin(file, results[i].ecgRIndex);

        // Section 2: PPG Max Amps
        readVecFromBin(file, results[i].ppgMaxAmps);

        // Section 3: PPG Min Amps
        readVecFromBin(file, results[i].ppgMinAmps);

        // Section 4: Pairs
        uint64_t numPairs = 0;
        file.read(reinterpret_cast<char*>(&numPairs), 8);
        for (uint64_t j = 0; j < numPairs; ++j) {
            int64_t ppgIdx, ecgIdx; // CHANGED: Read as signed int64_t
            file.read(reinterpret_cast<char*>(&ppgIdx), 8);
            file.read(reinterpret_cast<char*>(&ecgIdx), 8);
            results[i].pairs.push_back({ static_cast<double>(ppgIdx), static_cast<double>(ecgIdx) });
        }
    }
    return results;
}

// Helper to read a field from a MAT struct cell
std::vector<double> getMatField(matvar_t* cell, const std::string& fieldName) {
    std::vector<double> vec;
    matvar_t* field = Mat_VarGetStructFieldByName(cell, fieldName.c_str(), 0);
    if (field && field->data) {
        size_t n = 1;
        for (int i = 0; i < field->rank; ++i) n *= field->dims[i];
        double* data = (double*)field->data;
        for (size_t j = 0; j < n; ++j) vec.push_back(data[j]);
    }
    return vec;
}

// 2. Reads ALL sections from the MAT file
std::vector<ComparisonData> readAllFromMat(const std::string& path) {
    std::vector<ComparisonData> results;
    mat_t* matfp = Mat_Open(path.c_str(), MAT_ACC_RDONLY);
    if (!matfp) return results;

    matvar_t* wave_data = Mat_VarRead(matfp, "wave_data");
    if (!wave_data) {
        Mat_Close(matfp);
        return results;
    }

    size_t numBins = wave_data->dims[0] * wave_data->dims[1];
    results.resize(numBins);

    for (size_t i = 0; i < numBins; ++i) {
        matvar_t* cell = Mat_VarGetCell(wave_data, i);
        if (!cell) continue;

        results[i].ecgRIndex = getMatField(cell, "ecgRIndex");
        results[i].ppgMaxAmps = getMatField(cell, "ppgMaxAmps");
        results[i].ppgMinAmps = getMatField(cell, "ppgMinAmps");

        matvar_t* pair_var = Mat_VarGetStructFieldByName(cell, "pairs", 0);
        if (pair_var && pair_var->data && pair_var->rank == 2) {
            size_t rows = pair_var->dims[0];
            double* data = (double*)pair_var->data;
            for (size_t r = 0; r < rows; ++r) {
                results[i].pairs.push_back({ data[r], data[r + rows] });
            }
        }
    }
    Mat_VarFree(wave_data);
    Mat_Close(matfp);
    return results;
}

void compareVectors(const std::string& label, const std::vector<double>& bin, const std::vector<double>& mat) {
    std::cout << "  " << label << " -> C++: " << bin.size() << " | MAT: " << mat.size();
    if (bin.size() != mat.size()) {
        std::cout << " [MISMATCH]";
    }
    std::cout << std::endl;

    size_t limit = std::min(bin.size(), mat.size());
    for (size_t k = 0; k < limit; ++k) {
        // Checking against (bin[k] + 1) because MATLAB is 1-indexed and C++ is 0-indexed
        if ((bin[k] + 1 != mat[k]) && !(bin[k] == -1 && mat[k] == -1)) {
            // Note: If bin[k] is -1, bin[k]+1 is 0. If MATLAB shows -1, this will correctly flag it.
            std::cout << "    ! Offset at index " << k << ": BIN=" << bin[k] << ", MAT=" << mat[k] << std::endl;
        }
    }
}

int main() {
    std::string binPath = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\cpp\mesa_files\3010023_20110817_annealed_wave_data.bin)";
    std::string matPath = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\matlab\3010023_20110817_wave_data.mat)";

    auto binData = readAllFromBin(binPath);
    auto matData = readAllFromMat(matPath);

    size_t compare_limit = std::min({ (size_t)5, binData.size(), matData.size() });

    for (size_t i = 0; i < compare_limit; ++i) {
        std::cout << "\n--- Bin " << i << " ---" << std::endl;

        compareVectors("ECG R-Peaks", binData[i].ecgRIndex, matData[i].ecgRIndex);
        compareVectors("PPG Max (Systolic)", binData[i].ppgMaxAmps, matData[i].ppgMaxAmps);
        compareVectors("PPG Min (Valleys) ", binData[i].ppgMinAmps, matData[i].ppgMinAmps);

        std::cout << "  Pairs -> C++: " << binData[i].pairs.size() << " | MAT: " << matData[i].pairs.size();
        if (binData[i].pairs.size() != matData[i].pairs.size()) std::cout << " [MISMATCH]";
        std::cout << std::endl;

        size_t pLimit = std::min(binData[i].pairs.size(), matData[i].pairs.size());
        for (size_t k = 0; k < pLimit; ++k) {
            double binP = binData[i].pairs[k].first;
            double binE = binData[i].pairs[k].second;
            double matP = matData[i].pairs[k].first;
            double matE = matData[i].pairs[k].second;

            // Handle the +1 offset for comparison, but if value is -1, it stays -1 in comparison logic 
            if (binP + 1 != matP || (binE + 1 !=  matE && !(binE == -1 && matE == -1))) {
                std::cout << "    ! Pair " << k << " mismatch: BIN=[" << binP << "," << binE
                    << "], MAT=[" << matP << "," << matE << "]" << std::endl;
            }
        }
    }
    return 0;
}
