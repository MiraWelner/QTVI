#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <matio.h> 
#include <cmath>
#include <algorithm>
#include <cstdint>

struct ComparisonData {
    std::vector<double> ecgRIndex;
    std::vector<double> ppgMaxAmps;
    std::vector<double> ppgMinAmps;
    std::vector<double> ppgSignal;
    std::vector<std::pair<double, double>> pairs;
};

void readVecFromBin(std::ifstream& file, std::vector<double>& vec) {
    uint64_t count = 0;
    if (!file.read(reinterpret_cast<char*>(&count), 8)) return;
    vec.resize(count);
    for (uint64_t i = 0; i < count; ++i) {
        int64_t val;
        file.read(reinterpret_cast<char*>(&val), 8);
        vec[i] = static_cast<double>(val);
    }
}

void readDoubleVecFromBin(std::ifstream& file, std::vector<double>& vec) {
    uint64_t count = 0;
    if (!file.read(reinterpret_cast<char*>(&count), 8)) return;
    vec.resize(count);
    file.read(reinterpret_cast<char*>(vec.data()), count * 8);
}

std::vector<ComparisonData> readAllFromBin(const std::string& path) {
    std::vector<ComparisonData> results;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return results;

    uint64_t numBins = 0;
    file.read(reinterpret_cast<char*>(&numBins), 8);
    results.resize(numBins);

    for (uint64_t i = 0; i < numBins; ++i) {
        readVecFromBin(file, results[i].ecgRIndex);
        readVecFromBin(file, results[i].ppgMaxAmps);
        readVecFromBin(file, results[i].ppgMinAmps);
        readDoubleVecFromBin(file, results[i].ppgSignal);

        uint64_t numPairs = 0;
        file.read(reinterpret_cast<char*>(&numPairs), 8);
        for (uint64_t j = 0; j < numPairs; ++j) {
            int64_t ppgIdx, ecgIdx;
            file.read(reinterpret_cast<char*>(&ppgIdx), 8);
            file.read(reinterpret_cast<char*>(&ecgIdx), 8);
            results[i].pairs.push_back({ static_cast<double>(ppgIdx), static_cast<double>(ecgIdx) });
        }
    }
    return results;
}

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
        results[i].ppgSignal = getMatField(cell, "ppgSeg");

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

void compareVectors(std::ostream& out, const std::string& label, const std::vector<double>& bin, const std::vector<double>& mat) {
    out << "  " << label << " -> C++: " << bin.size() << " | MAT: " << mat.size();
    if (bin.size() != mat.size()) {
        out << " [MISMATCH]";
    }
    out << std::endl;

    size_t binIdx = 0, matIdx = 0;
    while (binIdx < bin.size() && matIdx < mat.size()) {
        if (bin[binIdx] == mat[matIdx]) {
            binIdx++;
            matIdx++;
        }
        else {
            if (binIdx + 1 < bin.size() && bin[binIdx + 1] == mat[matIdx]) {
                out << "    ! Extra element in C++ at index " << binIdx << ": " << bin[binIdx] << std::endl;
                binIdx++;
            }
            else if (matIdx + 1 < mat.size() && bin[binIdx] == mat[matIdx + 1]) {
                out << "    ! Extra element in MAT at index " << matIdx << ": " << mat[matIdx] << std::endl;
                matIdx++;
            }
            else {
                out << "    ! Value mismatch at C++ idx " << binIdx << " / MAT idx " << matIdx
                    << ": C++=" << bin[binIdx] << ", MAT=" << mat[matIdx] << std::endl;
                binIdx++;
                matIdx++;
            }
        }
    }

    while (binIdx < bin.size()) {
        out << "    ! Extra element in C++ at index " << binIdx << ": " << bin[binIdx] << std::endl;
        binIdx++;
    }
    while (matIdx < mat.size()) {
        out << "    ! Extra element in MAT at index " << matIdx << ": " << mat[matIdx] << std::endl;
        matIdx++;
    }
}

void comparePairs(std::ostream& out, const std::vector<std::pair<double, double>>& bin, const std::vector<std::pair<double, double>>& mat) {
    out << "  Pairs -> C++: " << bin.size() << " | MAT: " << mat.size();
    if (bin.size() != mat.size()) out << " [MISMATCH]";
    out << std::endl;

    size_t binIdx = 0, matIdx = 0;
    while (binIdx < bin.size() && matIdx < mat.size()) {
        if (bin[binIdx] == mat[matIdx]) {
            binIdx++;
            matIdx++;
        }
        else {
            if (binIdx + 1 < bin.size() && bin[binIdx + 1] == mat[matIdx]) {
                out << "    ! Extra pair in C++ at index " << binIdx << ": ["
                    << bin[binIdx].first << "," << bin[binIdx].second << "]" << std::endl;
                binIdx++;
            }
            else if (matIdx + 1 < mat.size() && bin[binIdx] == mat[matIdx + 1]) {
                out << "    ! Extra pair in MAT at index " << matIdx << ": ["
                    << mat[matIdx].first << "," << mat[matIdx].second << "]" << std::endl;
                matIdx++;
            }
            else {
                out << "    ! Pair mismatch at C++ idx " << binIdx << " / MAT idx " << matIdx
                    << ": C++=[" << bin[binIdx].first << "," << bin[binIdx].second
                    << "], MAT=[" << mat[matIdx].first << "," << mat[matIdx].second << "]" << std::endl;
                binIdx++;
                matIdx++;
            }
        }
    }

    while (binIdx < bin.size()) {
        out << "    ! Extra pair in C++ at index " << binIdx << ": ["
            << bin[binIdx].first << "," << bin[binIdx].second << "]" << std::endl;
        binIdx++;
    }
    while (matIdx < mat.size()) {
        out << "    ! Extra pair in MAT at index " << matIdx << ": ["
            << mat[matIdx].first << "," << mat[matIdx].second << "]" << std::endl;
        matIdx++;
    }
}

int main() {
    std::string ID = "3010023_20110817";
    std::string binPath = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\cpp\mesa_files\)" + ID + "_annealed_wave_data.bin";
    std::string matPath = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\matlab\)" + ID + "_wave_data.mat";
    std::string outPath = ID + ".txt";

    std::ofstream outFile(outPath);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open output file: " << outPath << std::endl;
        return 1;
    }

    auto binData = readAllFromBin(binPath);
    auto matData = readAllFromMat(matPath);

    size_t max_bins = std::min(binData.size(), matData.size());

    for (size_t i = 0; i < max_bins; ++i) {
        outFile << "\n--- Bin " << i << " ---" << std::endl;

        compareVectors(outFile, "ECG R-Peaks", binData[i].ecgRIndex, matData[i].ecgRIndex);
        compareVectors(outFile, "PPG Max (Systolic)", binData[i].ppgMaxAmps, matData[i].ppgMaxAmps);
        compareVectors(outFile, "PPG Min (Valleys) ", binData[i].ppgMinAmps, matData[i].ppgMinAmps);

        comparePairs(outFile, binData[i].pairs, matData[i].pairs);
    }

    outFile.close();
    std::cout << "Comparison complete. Results written to: " << outPath << std::endl;

    return 0;
}
