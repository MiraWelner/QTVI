#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <matio.h>

#include <numeric>
#include <algorithm>


// Simple Progress Bar Helper
void printProgress(double progress) {
    int barWidth = 50;
    std::cout << "[";
    int pos = (int)(barWidth * progress);
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

bool compareDoubles(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) return true;
    return std::abs(a - b) < 1e-9;
}

int main() {
    // HARDCODED PATHS
    std::string matPath = "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\3_annealed_files\\matlab\\3010155_20110511_annealedSegments.mat";
    std::string binPath = "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\3_annealed_files\\cpp\\mesa_bins\\3010155_20110511_annealed.bin";

    // 1. LOAD C++ BIN DATA
    std::cout << "Loading C++ Bin Data..." << std::endl;
    std::ifstream binFile(binPath, std::ios::binary);
    if (!binFile.is_open()) { std::cerr << "Error opening bin file: " << binPath << std::endl; return 1; }

    uint64_t binNumSegments;
    binFile.read(reinterpret_cast<char*>(&binNumSegments), sizeof(uint64_t));

    std::vector<std::vector<double>> binResults(binNumSegments);
    for (uint64_t i = 0; i < binNumSegments; ++i) {
        uint32_t pSize;
        binFile.read(reinterpret_cast<char*>(&pSize), sizeof(uint32_t));
        binResults[i].resize(pSize);
        binFile.read(reinterpret_cast<char*>(binResults[i].data()), pSize * sizeof(double));
    }
    std::cout << "Loaded " << binNumSegments << " segments from C++ Bin." << std::endl;

    // 2. LOAD MATLAB MAT DATA
    std::cout << "Loading MATLAB Mat Data..." << std::endl;
    mat_t* matfp = Mat_Open(matPath.c_str(), MAT_ACC_RDONLY);
    if (!matfp) { std::cerr << "Error opening mat file: " << matPath << std::endl; return 1; }

    matvar_t* annSegs = Mat_VarRead(matfp, "annealedSegments");
    if (!annSegs) {
        std::cerr << "Variable 'annealedSegments' not found in MAT file." << std::endl;
        Mat_Close(matfp); return 1;
    }

    // MANUAL ELEMENT COUNT FOR CELL ARRAY
    size_t matNumSegments = 1;
    for (int i = 0; i < annSegs->rank; i++) {
        matNumSegments *= annSegs->dims[i];
    }

    if (matNumSegments != binNumSegments) {
        std::cout << "\nCRITICAL FAILURE: Segment count mismatch!" << std::endl;
        std::cout << "MATLAB: " << matNumSegments << " vs C++: " << binNumSegments << std::endl;
        Mat_VarFree(annSegs); Mat_Close(matfp); return 1;
    }

    std::cout << "Starting statistical comparison (Mean, Median, Range)..." << std::endl;
    bool allMatch = true;

    for (size_t i = 0; i < matNumSegments; ++i) {
        matvar_t* cell = Mat_VarGetCell(annSegs, (int)i);
        if (!cell) { allMatch = false; break; }

        matvar_t* poField = Mat_VarGetStructFieldByName(cell, "po", 0);
        if (!poField) { allMatch = false; break; }

        // 1. Get MATLAB data and size
        double* matData = static_cast<double*>(poField->data);
        size_t matDataSize = 1;
        for (int k = 0; k < poField->rank; k++) matDataSize *= poField->dims[k];

        // 2. Get C++ data and size
        const std::vector<double>& cppData = binResults[i]; // Assuming binResults[i] is std::vector<double>
        size_t cppDataSize = cppData.size();

        // 3. Helper to calculate stats (lambda function)
        auto getStats = [](const double* data, size_t size) {
            if (size == 0) return std::make_tuple(0.0, 0.0, 0.0, 0.0); // Mean, Median, Min, Max

            // Mean and Range
            double sum = 0;
            double minVal = data[0];
            double maxVal = data[0];
            for (size_t i = 0; i < size; ++i) {
                sum += data[i];
                if (data[i] < minVal) minVal = data[i];
                if (data[i] > maxVal) maxVal = data[i];
            }
            double mean = sum / size;

            // Median (requires sorting a copy)
            std::vector<double> sortedCopy(data, data + size);
            std::sort(sortedCopy.begin(), sortedCopy.end());
            double median = (size % 2 == 0)
                ? (sortedCopy[size / 2 - 1] + sortedCopy[size / 2]) / 2.0
                : sortedCopy[size / 2];

            return std::make_tuple(mean, median, minVal, maxVal);
            };

        // 4. Calculate Stats
        auto [mMean, mMedian, mMin, mMax] = getStats(matData, matDataSize);
        auto [cMean, cMedian, cMin, cMax] = getStats(cppData.data(), cppDataSize);

        // 5. Compare (using a very small epsilon for floating point mean/median)
        auto isNear = [](double a, double b) { return std::abs(a - b) < 1e-5; };

        bool match = isNear(mMean, cMean) &&
            isNear(mMedian, cMedian) &&
            isNear(mMin, cMin) &&
            isNear(mMax, cMax);

        if (!match) {
            std::cout << "\nFAILED: Statistical mismatch at segment " << i << std::endl;
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "        |   MEAN   |  MEDIAN  |   MIN    |   MAX    |" << std::endl;
            std::cout << "MATLAB: | " << mMean << " | " << mMedian << " | " << mMin << " | " << mMax << " |" << std::endl;
            std::cout << "C++:    | " << cMean << " | " << cMedian << " | " << cMin << " | " << cMax << " |" << std::endl;
            std::cout << "Size:   MAT(" << matDataSize << ") vs C++(" << cppDataSize << ")" << std::endl;
            allMatch = false;
            break;
        }

        printProgress(static_cast<double>(i + 1) / matNumSegments);
    }

    if (allMatch) {
        std::cout << "\nVERIFICATION SUCCESSFUL: All segments match statistically." << std::endl;
    }

    Mat_VarFree(annSegs);
    Mat_Close(matfp);
    return allMatch ? 0 : 1;
}
