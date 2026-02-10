#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <matio.h>

struct SegmentData {
    std::vector<double> ppg;
    std::vector<double> ecg;
    std::vector<double> sleep;
};

// ============================================================================
// LOAD FROM .MAT FILE (Directly from MATLAB Output)
// ============================================================================
std::vector<SegmentData> LoadFromMat(const std::string& path) {
    std::vector<SegmentData> segments;
    mat_t* matfp = Mat_Open(path.c_str(), MAT_ACC_RDONLY);
    if (!matfp) throw std::runtime_error("Could not open MAT file: " + path);

    matvar_t* annealedVar = Mat_VarRead(matfp, "annealedSegments");
    if (!annealedVar) {
        Mat_Close(matfp);
        throw std::runtime_error("Variable 'annealedSegments' not found in MAT file.");
    }

    size_t nBins = 1;
    for (int i = 0; i < annealedVar->rank; i++) nBins *= annealedVar->dims[i];
    segments.resize(nBins);

    for (size_t i = 0; i < nBins; ++i) {
        matvar_t* cell = Mat_VarGetCell(annealedVar, (int)i);
        if (cell && cell->class_type == MAT_C_STRUCT) {
            auto loadField = [&](const char* name, std::vector<double>& dest) {
                matvar_t* f = Mat_VarGetStructFieldByName(cell, name, 0);
                if (f && f->data) {
                    double* d = (double*)f->data;
                    size_t n = 1;
                    for (int k = 0; k < f->rank; k++) n *= f->dims[k];
                    dest.assign(d, d + n);
                }
                };
            loadField("po", segments[i].ppg);
            loadField("ecg", segments[i].ecg);
            loadField("sleep_stages", segments[i].sleep);
        }
    }
    Mat_VarFree(annealedVar);
    Mat_Close(matfp);
    return segments;
}

// ============================================================================
// LOAD FROM C++ .BIN (Updated for 8-byte Headers and 8-byte Sleep)
// ============================================================================
std::vector<SegmentData> LoadFromCppBin(const std::string& path) {
    std::vector<SegmentData> segments;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return segments;

    uint64_t nBins;
    file.read((char*)&nBins, 8);
    segments.resize(nBins);

    for (uint64_t i = 0; i < nBins; ++i) {
        uint64_t pS, eS, sS; // Updated to 8-byte headers

        // Read PPG
        file.read((char*)&pS, 8);
        segments[i].ppg.resize(pS);
        file.read((char*)segments[i].ppg.data(), pS * 8);

        // Read ECG
        file.read((char*)&eS, 8);
        segments[i].ecg.resize(eS);
        file.read((char*)segments[i].ecg.data(), eS * 8);

        // Read Sleep (Now stored as 8-byte double)
        file.read((char*)&sS, 8);
        segments[i].sleep.resize(sS);
        file.read((char*)segments[i].sleep.data(), sS * 8);
    }
    return segments;
}

// ============================================================================
// COMPARISON LOGIC (STRICT)
// ============================================================================
void Verify(int bin_num, const std::string& label, const std::vector<double>& ref, const std::vector<double>& test) {
    if (ref.size() != test.size()) {
        std::cout << "\nBin #" << bin_num << ": [FAIL] " << label << " Size Mismatch! Mat: " << ref.size() << " C++: " << test.size() << std::endl;
        return;
    }
    double max_err = 0;
    for (size_t i = 0; i < ref.size(); ++i) {
        double err = std::abs(ref[i] - test[i]);
        if (err > max_err) max_err = err;
    }
    if (max_err >= 1e-9) {
        std::cout << "\nBin #" << bin_num << " [FAIL] " << label << " (Max Diff: " << max_err << ")" << std::endl;
    }
}

int main() {
    std::string mat_file = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\3_annealed_files\matlab\3010155_20110511_annealedSegments.mat)";
    std::string bin_file = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\3_annealed_files\cpp\mesa_bins\3010155_20110511_annealed.bin)";

    try {
        auto matData = LoadFromMat(mat_file);
        auto cppData = LoadFromCppBin(bin_file);

        size_t commonBins = std::min(matData.size(), cppData.size());
        for (size_t b = 0; b < commonBins; ++b) {
            Verify(b + 1, "PPG", matData[b].ppg, cppData[b].ppg);
            Verify(b + 1, "ECG", matData[b].ecg, cppData[b].ecg);
            Verify(b + 1, "Sleep", matData[b].sleep, cppData[b].sleep);
        }
    }
    catch (const std::exception& e) { std::cerr << "Error: " << e.what() << std::endl; }
    return 0;
}
