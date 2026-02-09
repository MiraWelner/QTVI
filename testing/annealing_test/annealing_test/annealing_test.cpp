#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <matio.h> // Ensure you have matio installed and linked

struct SegmentData {
    std::vector<double> ppg;
    std::vector<double> ecg;
    std::vector<double> sleep; // MATLAB stores as double
};

// ============================================================================
// LOAD FROM .MAT FILE (Directly from MATLAB Output)
// ============================================================================
std::vector<SegmentData> LoadFromMat(const std::string& path) {
    std::vector<SegmentData> segments;
    mat_t* matfp = Mat_Open(path.c_str(), MAT_ACC_RDONLY);
    if (!matfp) throw std::runtime_error("Could not open MAT file: " + path);

    // Read the 'annealedSegments' cell array
    matvar_t* annealedVar = Mat_VarRead(matfp, "annealedSegments");
    if (!annealedVar) {
        Mat_Close(matfp);
        throw std::runtime_error("Variable 'annealedSegments' not found in MAT file.");
    }

    // FIX: Calculate total number of elements from dimensions
    size_t nBins = 1;
    for (int i = 0; i < annealedVar->rank; i++) {
        nBins *= annealedVar->dims[i];
    }

    segments.resize(nBins);

    for (size_t i = 0; i < nBins; ++i) {
        // Mat_VarGetCell retrieves the struct inside each cell
        matvar_t* cell = Mat_VarGetCell(annealedVar, (int)i);
        if (cell && cell->class_type == MAT_C_STRUCT) {

            // 1. Get PPG ('po' field)
            matvar_t* poVar = Mat_VarGetStructFieldByName(cell, "po", 0);
            if (poVar && poVar->data) {
                double* data = (double*)poVar->data;
                size_t n = 1;
                for (int d = 0; d < poVar->rank; d++) n *= poVar->dims[d];
                segments[i].ppg.assign(data, data + n);
            }

            // 2. Get ECG ('ecg' field)
            matvar_t* ecgVar = Mat_VarGetStructFieldByName(cell, "ecg", 0);
            if (ecgVar && ecgVar->data) {
                double* data = (double*)ecgVar->data;
                size_t n = 1;
                for (int d = 0; d < ecgVar->rank; d++) n *= ecgVar->dims[d];
                segments[i].ecg.assign(data, data + n);
            }

            // 3. Get Sleep ('sleep_stages' field)
            matvar_t* slpVar = Mat_VarGetStructFieldByName(cell, "sleep_stages", 0);
            if (slpVar && slpVar->data) {
                double* data = (double*)slpVar->data;
                size_t n = 1;
                for (int d = 0; d < slpVar->rank; d++) n *= slpVar->dims[d];
                segments[i].sleep.assign(data, data + n);
            }
        }
    }

    Mat_VarFree(annealedVar);
    Mat_Close(matfp);
    return segments;
}

// ============================================================================
// LOAD FROM C++ .BIN (Our Application Output)
// ============================================================================
std::vector<SegmentData> LoadFromCppBin(const std::string& path) {
    std::vector<SegmentData> segments;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return segments;

    uint64_t nBins;
    file.read((char*)&nBins, 8);
    segments.resize(nBins);

    for (uint64_t i = 0; i < nBins; ++i) {
        uint32_t pS, eS, sS;
        // Read PPG
        file.read((char*)&pS, 4);
        segments[i].ppg.resize(pS);
        file.read((char*)segments[i].ppg.data(), (uint64_t)pS * 8);
        // Read ECG
        file.read((char*)&eS, 4);
        segments[i].ecg.resize(eS);
        file.read((char*)segments[i].ecg.data(), (uint64_t)eS * 8);
        // Read Sleep (Stored as int in C++, convert to double for comparison)
        file.read((char*)&sS, 4);
        std::vector<int> tempSleep(sS);
        file.read((char*)tempSleep.data(), (uint64_t)sS * 4);
        for (int val : tempSleep) segments[i].sleep.push_back((double)val);
    }
    return segments;
}

// ============================================================================
// COMPARISON LOGIC
// ============================================================================
void Verify(const std::string& label, const std::vector<double>& ref, const std::vector<double>& test) {
    if (ref.size() != test.size()) {
        std::cout << "  [FAIL] " << label << " Size Mismatch! Mat: " << ref.size() << " C++: " << test.size() << std::endl;
        return;
    }
    double max_err = 0;
    for (size_t i = 0; i < ref.size(); ++i) {
        double err = std::abs(ref[i] - test[i]);
        if (err > max_err) max_err = err;
    }
    if (max_err >= 1e-9) {
        std::cout << "[FAIL]" <<  label << " (Max Diff: " << max_err << ")" << std::endl;
    }
}

int main() {
    std::string mat_file = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\3_annealed_files\matlab\3010155_20110511_annealedSegments.mat)";
    std::string bin_file = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\3_annealed_files\cpp\mesa_bins\3010155_20110511_annealed.bin)";

    try {
        std::cout << "Loading MATLAB .mat file..." << std::endl;
        auto matData = LoadFromMat(mat_file);

        std::cout << "Loading C++ .bin file..." << std::endl;
        auto cppData = LoadFromCppBin(bin_file);

        size_t commonBins = std::min(matData.size(), cppData.size());
        for (size_t b = 0; b < commonBins; ++b) {
            std::cout << "\nBin #" << b + 1 << ":" << std::endl;
            Verify("PPG", matData[b].ppg, cppData[b].ppg);
            Verify("ECG", matData[b].ecg, cppData[b].ecg);
            Verify("Sleep", matData[b].sleep, cppData[b].sleep);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}
