#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <matio.h>

// ============================================================================
// DATA STRUCTURES
// ============================================================================

struct AnnealedData {
    struct Bin {
        std::vector<double> ppg;
        std::vector<double> ecg;
        std::vector<double> sleep;
    };
    std::vector<Bin> bins;
};

// ============================================================================
// HELPERS
// ============================================================================

// Helper to manually calculate element count from dims (fixes C2039 error)
size_t getMatElementCount(const matvar_t* var) {
    if (!var || var->rank == 0) return 0;
    size_t count = 1;
    for (int i = 0; i < var->rank; i++) {
        count *= var->dims[i];
    }
    return count;
}

// ============================================================================
// C++ BINARY READER
// ============================================================================

AnnealedData readCppBin(const std::string& path) {
    AnnealedData data;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Could not open C++ bin: " + path);

    uint64_t numBins;
    if (!file.read((char*)&numBins, 8)) return data;
    data.bins.resize(numBins);

    for (uint64_t i = 0; i < numBins; ++i) {
        uint64_t pS, eS, sS;

        // Read PPG
        file.read((char*)&pS, 8);
        data.bins[i].ppg.resize(pS);
        file.read((char*)data.bins[i].ppg.data(), pS * 8);

        // Read ECG
        file.read((char*)&eS, 8);
        data.bins[i].ecg.resize(eS);
        file.read((char*)data.bins[i].ecg.data(), eS * 8);

        // Read Sleep
        file.read((char*)&sS, 8);
        data.bins[i].sleep.resize(sS);
        file.read((char*)data.bins[i].sleep.data(), sS * 8);
    }
    return data;
}

// ============================================================================
// MATLAB .MAT READER (USING MATIO)
// ============================================================================

AnnealedData readMatlabMat(const std::string& path) {
    AnnealedData data;
    mat_t* matfp = Mat_Open(path.c_str(), MAT_ACC_RDONLY);
    if (!matfp) throw std::runtime_error("Could not open MATLAB mat: " + path);

    // Variable name in .mat is 'annealedSegments'
    matvar_t* annSegs = Mat_VarRead(matfp, "annealedSegments");
    if (!annSegs) {
        Mat_Close(matfp);
        throw std::runtime_error("Variable 'annealedSegments' not found in " + path);
    }

    size_t numBins = getMatElementCount(annSegs);
    data.bins.resize(numBins);

    for (size_t i = 0; i < numBins; ++i) {
        matvar_t* cell = Mat_VarGetCell(annSegs, (int)i);
        if (!cell) continue;

        // 1. Get PPG (Field: 'po')
        matvar_t* ppgVar = Mat_VarGetStructFieldByName(cell, "po", 0);
        if (ppgVar && ppgVar->data) {
            size_t n = getMatElementCount(ppgVar);
            data.bins[i].ppg.assign((double*)ppgVar->data, (double*)ppgVar->data + n);
        }

        // 2. Get ECG (Field: 'ecg')
        matvar_t* ecgVar = Mat_VarGetStructFieldByName(cell, "ecg", 0);
        if (ecgVar && ecgVar->data) {
            size_t n = getMatElementCount(ecgVar);
            data.bins[i].ecg.assign((double*)ecgVar->data, (double*)ecgVar->data + n);
        }

        // 3. Get Sleep Stages (Field: 'sleep_stages')
        matvar_t* sleepVar = Mat_VarGetStructFieldByName(cell, "sleep_stages", 0);
        if (sleepVar && sleepVar->data) {
            size_t n = getMatElementCount(sleepVar);
            data.bins[i].sleep.assign((double*)sleepVar->data, (double*)sleepVar->data + n);
        }
    }

    Mat_VarFree(annSegs);
    Mat_Close(matfp);
    return data;
}

// ============================================================================
// COMPARISON ENGINE
// ============================================================================

void runComparison(const std::string& cppPath, const std::string& matPath) {
    std::cout << "====================================================" << std::endl;
    std::cout << "        ANNEALING COMPARISON (MAT vs BIN)" << std::endl;
    std::cout << "====================================================" << std::endl;

    try {
        std::cout << "Loading C++ File: " << cppPath << "..." << std::endl;
        AnnealedData cpp = readCppBin(cppPath);

        std::cout << "Loading MATLAB File: " << matPath << "..." << std::endl;
        AnnealedData mat = readMatlabMat(matPath);

        std::cout << "\n--- Header Data ---" << std::endl;
        std::cout << "C++ Bin Count:    " << cpp.bins.size() << std::endl;
        std::cout << "MAT Bin Count:    " << mat.bins.size() << std::endl;

        if (cpp.bins.size() != mat.bins.size()) {
            std::cout << "[WARN] Bin count mismatch!" << std::endl;
        }

        size_t totalBins = std::min(cpp.bins.size(), mat.bins.size());
        double maxPPGDiff = 0, maxECGDiff = 0;
        int sleepErrorCount = 0;

        for (size_t i = 0; i < totalBins; ++i) {
            // Signal Check
            size_t pLen = std::min(cpp.bins[i].ppg.size(), mat.bins[i].ppg.size());
            for (size_t p = 0; p < pLen; ++p)
                maxPPGDiff = std::max(maxPPGDiff, std::abs(cpp.bins[i].ppg[p] - mat.bins[i].ppg[p]));

            size_t eLen = std::min(cpp.bins[i].ecg.size(), mat.bins[i].ecg.size());
            for (size_t e = 0; e < eLen; ++e)
                maxECGDiff = std::max(maxECGDiff, std::abs(cpp.bins[i].ecg[e] - mat.bins[i].ecg[e]));

            // Sleep Check
            bool sleepMismatch = (cpp.bins[i].sleep.size() != mat.bins[i].sleep.size());
            if (!sleepMismatch) {
                for (size_t s = 0; s < cpp.bins[i].sleep.size(); ++s) {
                    if (std::abs(cpp.bins[i].sleep[s] - mat.bins[i].sleep[s]) > 0.001) {
                        sleepMismatch = true; break;
                    }
                }
            }

            if (sleepMismatch) {
                sleepErrorCount++;
                std::cout << "\n[FAIL] Sleep Mismatch in Bin #" << i << std::endl;
                std::cout << "C++ (" << cpp.bins[i].sleep.size() << "): [ ";
                for (double v : cpp.bins[i].sleep) std::cout << (int)v << " ";
                std::cout << "]" << std::endl;

                std::cout << "MAT (" << mat.bins[i].sleep.size() << "): [ ";
                for (double v : mat.bins[i].sleep) std::cout << (int)v << " ";
                std::cout << "]" << std::endl;
            }
        }

        std::cout << "\n--- Final Comparison Stats ---" << std::endl;
        std::cout << "Max PPG Abs Diff:  " << maxPPGDiff << std::endl;
        std::cout << "Max ECG Abs Diff:  " << maxECGDiff << std::endl;
        std::cout << "Total Sleep Error Bins: " << sleepErrorCount << " / " << totalBins << std::endl;

        if (sleepErrorCount == 0 && maxPPGDiff < 1e-7 && cpp.bins.size() == mat.bins.size()) {
            std::cout << "\nRESULT: PASS - Files are logically identical." << std::endl;
        }
        else {
            std::cout << "\nRESULT: FAIL - Discrepancies found." << std::endl;
        }

    }
    catch (const std::exception& e) {
        std::cerr << "\nFATAL TEST ERROR: " << e.what() << std::endl;
    }
}

// ============================================================================
// MAIN
// ============================================================================

int main() {
    // Paths provided by user
    std::string cppAnnealedFile = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\3_annealed_files\cpp\mesa_bins\3010155_20110511_annealed.bin)";
    std::string matlabAnnealedFile = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\3_annealed_files\matlab\3010155_20110511_annealedSegments.mat)";

    runComparison(cppAnnealedFile, matlabAnnealedFile);

    std::cout << "\nPress Enter to exit..." << std::endl;
    std::cin.get();
    return 0;
}
