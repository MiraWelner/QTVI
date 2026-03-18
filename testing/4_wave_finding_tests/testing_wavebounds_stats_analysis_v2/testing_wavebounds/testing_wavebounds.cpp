#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <matio.h>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <cstdint>
#include <functional>

// ============================================================================
// Data structures
// ============================================================================

struct ComparisonData {
    double ecgSamplingRate = 0.0;
    double ppgSamplingRate = 0.0;
    std::vector<double> ecgRIndex;
    std::vector<double> ecgRIndex2;
    std::vector<double> ecgRIndex3;
    std::vector<double> ppgMaxAmps;
    std::vector<double> ppgMinAmps;
    std::vector<double> ppgSignal;
    std::vector<double> ecgSignal;
    std::vector<double> ecgSignal2;
    std::vector<double> ecgSignal3;
    std::vector<std::pair<double, double>> pairs;
};

struct DatasetStats {
    double mean = 0.0;
    double stdDev = 0.0;
    long long total = 0;
};

struct FileResult {
    std::string id;
    long long rPeaksBin = 0;
    long long rPeaksMat = 0;
    int binCountBin = 0;
    int binCountMat = 0;
    size_t binsWithZeroSSD = 0;
    size_t totalBeatsBin = 0;
    size_t totalBeatsMat = 0;
    double ssdMean = 0.0;
    double ssdStd = 0.0;
};

// ============================================================================
// Statistics helpers
// ============================================================================

template <typename Projection>
DatasetStats calculateStats(const std::vector<ComparisonData>& data, Projection proj) {
    if (data.empty()) return {};

    double sum = 0.0, sumSq = 0.0;
    for (const auto& item : data) {
        double val = static_cast<double>(proj(item));
        sum += val;
        sumSq += val * val;
    }

    double n = static_cast<double>(data.size());
    double mean = sum / n;
    double variance = std::max(0.0, (sumSq / n) - (mean * mean));
    return { mean, std::sqrt(variance), static_cast<long long>(sum) };
}

double computeRPeakLocationSSD(const ComparisonData& bin, const ComparisonData& mat) {
    size_t nBin = bin.ecgRIndex.size();
    size_t nMat = mat.ecgRIndex.size();
    size_t nMax = std::max(nBin, nMat);
    size_t nMin = std::min(nBin, nMat);

    double binLength = static_cast<double>(bin.ecgSignal.size());

    double ssd = 0.0;
    for (size_t i = 0; i < nMin; ++i) {
        double diff = bin.ecgRIndex[i] - mat.ecgRIndex[i];
        ssd += diff * diff;
    }

    size_t nExtra = nMax - nMin;
    ssd += nExtra * binLength * binLength;

    return ssd;
}

// ============================================================================
// Binary file I/O
// ============================================================================

void readInt64VecAsDouble(std::ifstream& file, std::vector<double>& vec) {
    uint64_t count = 0;
    if (!file.read(reinterpret_cast<char*>(&count), 8)) return;
    vec.resize(count);
    for (uint64_t i = 0; i < count; ++i) {
        int64_t val;
        file.read(reinterpret_cast<char*>(&val), 8);
        vec[i] = static_cast<double>(val);
    }
}

void readDoubleVec(std::ifstream& file, std::vector<double>& vec) {
    uint64_t count = 0;
    if (!file.read(reinterpret_cast<char*>(&count), 8)) return;
    vec.resize(count);
    if (count > 0)
        file.read(reinterpret_cast<char*>(vec.data()), count * 8);
}

std::vector<ComparisonData> readAllFromBin(const std::string& path) {
    std::vector<ComparisonData> results;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return results;

    uint64_t numBins = 0;
    if (!file.read(reinterpret_cast<char*>(&numBins), 8)) return results;
    results.resize(numBins);

    for (uint64_t i = 0; i < numBins; ++i) {
        auto& r = results[i];
        file.read(reinterpret_cast<char*>(&r.ecgSamplingRate), 8);
        file.read(reinterpret_cast<char*>(&r.ppgSamplingRate), 8);

        readInt64VecAsDouble(file, r.ecgRIndex);
        readInt64VecAsDouble(file, r.ecgRIndex2);
        readInt64VecAsDouble(file, r.ecgRIndex3);
        readInt64VecAsDouble(file, r.ppgMaxAmps);
        readInt64VecAsDouble(file, r.ppgMinAmps);

        readDoubleVec(file, r.ppgSignal);
        readDoubleVec(file, r.ecgSignal);
        readDoubleVec(file, r.ecgSignal2);
        readDoubleVec(file, r.ecgSignal3);

        uint64_t numPairs = 0;
        if (file.read(reinterpret_cast<char*>(&numPairs), 8) && numPairs < 100000) {
            r.pairs.reserve(numPairs);
            for (uint64_t j = 0; j < numPairs; ++j) {
                int64_t ppgIdx, ecgIdx;
                if (!file.read(reinterpret_cast<char*>(&ppgIdx), 8)) break;
                if (!file.read(reinterpret_cast<char*>(&ecgIdx), 8)) break;
                r.pairs.push_back({ static_cast<double>(ppgIdx), static_cast<double>(ecgIdx) });
            }
        }
    }
    return results;
}

// ============================================================================
// MAT file I/O
// ============================================================================

std::vector<double> getMatField(matvar_t* cell, const std::string& fieldName) {
    std::vector<double> vec;
    matvar_t* field = Mat_VarGetStructFieldByName(cell, fieldName.c_str(), 0);

    if (!field && fieldName == "ppgMaxAmps") {
        field = Mat_VarGetStructFieldByName(cell, "ppgMaxIdx", 0);
        if (!field) field = Mat_VarGetStructFieldByName(cell, "maxAmps", 0);
    }
    if (!field && fieldName == "ppgMinAmps")
        field = Mat_VarGetStructFieldByName(cell, "ppgMinIdx", 0);

    if (!field || !field->data) return vec;

    size_t n = 1;
    for (int i = 0; i < field->rank; ++i) n *= field->dims[i];

    vec.reserve(n);
    for (size_t j = 0; j < n; ++j) {
        if (field->class_type == MAT_C_DOUBLE)
            vec.push_back(static_cast<double*>(field->data)[j]);
        else if (field->class_type == MAT_C_UINT64)
            vec.push_back(static_cast<double>(static_cast<uint64_t*>(field->data)[j]));
        else if (field->class_type == MAT_C_INT64)
            vec.push_back(static_cast<double>(static_cast<int64_t*>(field->data)[j]));
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

        auto& r = results[i];
        auto eRate = getMatField(cell, "ecgSamplingRate");
        r.ecgSamplingRate = eRate.empty() ? 0.0 : eRate[0];

        auto pRate = getMatField(cell, "ppgSamplingRate");
        r.ppgSamplingRate = pRate.empty() ? 0.0 : pRate[0];

        r.ecgRIndex = getMatField(cell, "ecgRIndex");
        r.ppgMaxAmps = getMatField(cell, "ppgMaxAmps");
        r.ppgMinAmps = getMatField(cell, "ppgMinAmps");
        r.ppgSignal = getMatField(cell, "ppgSignal");
        r.ecgSignal = getMatField(cell, "ecgSignal");

        matvar_t* pair_var = Mat_VarGetStructFieldByName(cell, "pairs", 0);
        if (pair_var && pair_var->data && pair_var->rank == 2) {
            size_t rows = pair_var->dims[0];
            double* data = static_cast<double*>(pair_var->data);
            for (size_t row = 0; row < rows; ++row)
                r.pairs.push_back({ data[row], data[row + rows] });
        }
    }

    Mat_VarFree(wave_data);
    Mat_Close(matfp);
    return results;
}

// ============================================================================
// SVG generation
// ============================================================================

void generateSVGPlot(const std::string& filename,
    const std::vector<ComparisonData>& binData,
    const std::vector<ComparisonData>& matData,
    const std::string& label,
    std::function<double(const ComparisonData&)> selector) {

    size_t n = std::min(binData.size(), matData.size());
    if (n == 0) return;

    double maxVal = 0.0;
    for (size_t i = 0; i < n; ++i)
        maxVal = std::max({ maxVal, selector(binData[i]), selector(matData[i]) });
    maxVal = (maxVal == 0.0) ? 10.0 : maxVal * 1.1;

    constexpr int size = 600, padding = 80;
    constexpr int chartSize = size - padding * 2;

    std::ofstream svg(filename);
    if (!svg.is_open()) return;

    auto mapX = [&](double v) { return padding + v / maxVal * chartSize; };
    auto mapY = [&](double v) { return (size - padding) - v / maxVal * chartSize; };

    svg << "<svg width=\"" << size << "\" height=\"" << size
        << "\" xmlns=\"http://www.w3.org/2000/svg\">\n"
        << "<rect width=\"100%\" height=\"100%\" fill=\"#fcfcfc\"/>\n";

    svg << "<line x1=\"" << mapX(0) << "\" y1=\"" << mapY(0)
        << "\" x2=\"" << mapX(maxVal) << "\" y2=\"" << mapY(maxVal)
        << "\" stroke=\"#999\" stroke-width=\"1\" stroke-dasharray=\"5,5\"/>\n";

    svg << "<line x1=\"" << padding << "\" y1=\"" << (size - padding)
        << "\" x2=\"" << (size - padding) << "\" y2=\"" << (size - padding)
        << "\" stroke=\"black\" stroke-width=\"2\"/>\n";
    svg << "<line x1=\"" << padding << "\" y1=\"" << padding
        << "\" x2=\"" << padding << "\" y2=\"" << (size - padding)
        << "\" stroke=\"black\" stroke-width=\"2\"/>\n";

    svg << "<text x=\"" << size / 2 << "\" y=\"" << size - 20
        << "\" text-anchor=\"middle\" font-family=\"sans-serif\">Bin " << label << "</text>\n";
    svg << "<text x=\"25\" y=\"" << size / 2
        << "\" text-anchor=\"middle\" transform=\"rotate(-90 25," << size / 2
        << ")\" font-family=\"sans-serif\">Mat " << label << "</text>\n";
    svg << "<text x=\"" << size / 2 << "\" y=\"40\" text-anchor=\"middle\" "
        << "font-weight=\"bold\" font-family=\"sans-serif\">" << label << " Comparison</text>\n";

    for (size_t i = 0; i < n; ++i) {
        double x = selector(binData[i]);
        double y = selector(matData[i]);
        const char* color = (std::abs(x - y) < 0.0001) ? "#3498db" : "#e74c3c";

        svg << "<circle cx=\"" << mapX(x) << "\" cy=\"" << mapY(y)
            << "\" r=\"3.5\" fill=\"" << color << "\" fill-opacity=\"0.7\">\n"
            << "  <title>Bin: " << x << ", Mat: " << y << "</title>\n"
            << "</circle>\n";
    }

    svg << "</svg>";
}

// ============================================================================
// Per-file comparison
// ============================================================================

FileResult compareFile(const std::string& id,
    const std::vector<ComparisonData>& binData,
    const std::vector<ComparisonData>& matData) {

    size_t maxBins = std::min(binData.size(), matData.size());
    std::string outPath = id + ".txt";
    std::ofstream out(outPath);
    if (!out.is_open()) {
        std::cerr << "Failed to open output file: " << outPath << std::endl;
        return {};
    }

    auto rPeakProj = [](const ComparisonData& c) { return (double)c.ecgRIndex.size(); };
    auto ppgPeakProj = [](const ComparisonData& c) { return (double)c.ppgMaxAmps.size(); };
    auto ppgDipProj = [](const ComparisonData& c) { return (double)c.ppgMinAmps.size(); };
    auto pairProj = [](const ComparisonData& c) { return (double)c.pairs.size(); };

    DatasetStats binR = calculateStats(binData, rPeakProj);
    DatasetStats matR = calculateStats(matData, rPeakProj);
    DatasetStats binPPG = calculateStats(binData, ppgPeakProj);
    DatasetStats matPPG = calculateStats(matData, ppgPeakProj);
    DatasetStats binDip = calculateStats(binData, ppgDipProj);
    DatasetStats matDip = calculateStats(matData, ppgDipProj);
    DatasetStats binPrs = calculateStats(binData, pairProj);
    DatasetStats matPrs = calculateStats(matData, pairProj);

    size_t binsWithZeroSSD = 0;
    size_t totalBeatsBin = 0, totalBeatsMat = 0;
    std::vector<double> ssdValues;
    ssdValues.reserve(maxBins);

    for (size_t i = 0; i < maxBins; ++i) {
        double ssd = computeRPeakLocationSSD(binData[i], matData[i]);
        ssdValues.push_back(ssd);
        if (ssd == 0.0) ++binsWithZeroSSD;
        totalBeatsBin += binData[i].ecgRIndex.size();
        totalBeatsMat += matData[i].ecgRIndex.size();
    }

    // Compute SSD mean and std
    double ssdMean = 0.0, ssdStd = 0.0;
    if (!ssdValues.empty()) {
        double sum = 0.0, sumSq = 0.0;
        for (double v : ssdValues) {
            sum += v;
            sumSq += v * v;
        }
        double n = static_cast<double>(ssdValues.size());
        ssdMean = sum / n;
        ssdStd = (n > 1) ? std::sqrt((sumSq - sum * sum / n) / (n - 1)) : 0.0;
    }

    // --- Write summary ---
    out << "Data taken from MESA file " << id << "\n\n";
    out << "=== Summary ===\n";
    out << "Total bins: " << binData.size() << "\n";
    out << "Bins with R peaks not in same place: " << (binData.size() - binsWithZeroSSD) << "\n";
    out << "R-peak index SSD — Mean: " << std::fixed << std::setprecision(4) << ssdMean
        << "  Std: " << ssdStd << "\n\n";

    out << "\t\tTotal in .bin\tTotal in .mat\tDiff\n";
    auto writeRow = [&](const char* label, DatasetStats b, DatasetStats m) {
        out << label << "\t" << b.total << "\t\t" << m.total
            << "\t\t" << (m.total - b.total) << "\n";
        };
    writeRow("R Peaks\t", binR, matR);
    writeRow("PPG Peaks", binPPG, matPPG);
    writeRow("PPG Dips", binDip, matDip);
    writeRow("Pairs\t", binPrs, matPrs);
    out << "\n";

    // --- Per-bin detail ---
    out << "bin #\tBeats(.bin)\tBeats(.mat)\tMean R idx(.bin)\tMean R idx(.mat)\tR-peak location SSD\n";

    for (size_t i = 0; i < maxBins; ++i) {
        size_t nBin = binData[i].ecgRIndex.size();
        size_t nMat = matData[i].ecgRIndex.size();

        auto meanIdx = [](const std::vector<double>& v) {
            if (v.empty()) return 0.0;
            return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
            };

        double meanBin = meanIdx(binData[i].ecgRIndex);
        double meanMat = meanIdx(matData[i].ecgRIndex);

        out << std::setw(3) << std::setfill('0') << i << std::setfill(' ')
            << "\t" << nBin << "\t\t" << nMat << "\t\t"
            << std::fixed << std::setprecision(1) << meanBin << "\t\t\t"
            << meanMat << "\t\t\t"
            << std::setprecision(4) << ssdValues[i] << "\n";
    }

    return { id, binR.total, matR.total,
             static_cast<int>(binData.size()), static_cast<int>(matData.size()),
             binsWithZeroSSD, totalBeatsBin, totalBeatsMat, ssdMean, ssdStd };
}

// ============================================================================
// All-files summary
// ============================================================================

void writeAllFilesSummary(const std::vector<FileResult>& results) {
    std::ofstream out("allfiles.txt");
    if (!out) return;

    long long sumBin = 0, sumMat = 0;
    int sumBinsBin = 0, sumBinsMat = 0;
    size_t sumZeroSSD = 0;

    for (const auto& r : results) {
        sumBin += r.rPeaksBin;
        sumMat += r.rPeaksMat;
        sumBinsBin += r.binCountBin;
        sumBinsMat += r.binCountMat;
        sumZeroSSD += r.binsWithZeroSSD;
    }

    out << "=== R-Peak Count Summary (All Files) ===\n";
    out << "\t\t\t\tR Peaks (.bin)\tR Peaks (.mat)\tDiff\tBins w/ same R index\t"
        << "Total Bins\t% Bins w/ same R index\tSSD Mean\tSSD Std\n";

    out << "All Files Sum\t\t\t"
        << sumBin << "\t\t" << sumMat << "\t\t" << (sumMat - sumBin) << "\t\t"
        << sumZeroSSD << "\t\t"
        << sumBinsBin << "\t\t"
        << std::fixed << std::setprecision(4)
        << static_cast<double>(sumZeroSSD) / sumBinsBin << "\n\n";

    for (const auto& r : results) {
        out << "Mesa ID=" << r.id << "\t"
            << r.rPeaksBin << "\t\t" << r.rPeaksMat << "\t\t"
            << (r.rPeaksMat - r.rPeaksBin) << "\t\t"
            << r.binsWithZeroSSD << "\t\t"
            << r.binCountBin << "\t\t"
            << std::fixed << std::setprecision(4)
            << static_cast<double>(r.binsWithZeroSSD) / r.binCountBin << "\t\t"
            << r.ssdMean << "\t\t" << r.ssdStd << "\n";
    }
}

// ============================================================================
// SVG helper lambdas
// ============================================================================

static auto rPeakCount = [](const ComparisonData& d) { return (double)d.ecgRIndex.size(); };
static auto ppgMaxCount = [](const ComparisonData& d) { return (double)d.ppgMaxAmps.size(); };
static auto ppgMinCount = [](const ComparisonData& d) { return (double)d.ppgMinAmps.size(); };
static auto pairCount = [](const ComparisonData& d) { return (double)d.pairs.size(); };

// ============================================================================
// Main
// ============================================================================

int main() {
    const std::string binDir = "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\4_wave_bound_files\\mesa_files\\";
    const std::string matDir = "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\4_wave_bound_files\\matlab\\";

    const std::vector<std::string> IDs = {
        "3010023_20110817", "3010104_20111027", "3010112_20110725",
        "3010139_20110210", "3010201_20120320", "3010228_20110426",
        "3010317_20110413", "3010457_20111109", "3010660_20120322",
        "3010724_20110811", "3010740_20110303"
    };

    std::vector<FileResult> results;
    std::vector<ComparisonData> allBinData, allMatData;

    for (const auto& id : IDs) {
        std::string binPath = binDir + id + "_annealed_wave_data.bin";
        std::string matPath = matDir + id + "_wave_data.mat";

        auto binData = readAllFromBin(binPath);
        auto matData = readAllFromMat(matPath);

        if (binData.empty()) std::cerr << "WARNING: No data from " << binPath << "\n";
        if (matData.empty()) std::cerr << "WARNING: No data from " << matPath << "\n";
        if (binData.empty() || matData.empty()) {
            std::cerr << "Skipping " << id << " due to missing data.\n";
            continue;
        }

        results.push_back(compareFile(id, binData, matData));

        generateSVGPlot(id + "_ecg_comparison.svg", binData, matData, "R-Peak Count", rPeakCount);
        generateSVGPlot(id + "_ppg_max_comparison.svg", binData, matData, "PPG Max Amps Count", ppgMaxCount);
        generateSVGPlot(id + "_ppg_min_comparison.svg", binData, matData, "PPG Min Amps Count", ppgMinCount);
        generateSVGPlot(id + "_pairs_comparison.svg", binData, matData, "Detected Pairs", pairCount);

        size_t maxBins = std::min(binData.size(), matData.size());
        for (size_t i = 0; i < maxBins; ++i) {
            allBinData.push_back(binData[i]);
            allMatData.push_back(matData[i]);
        }

        std::cout << "Comparison complete for: " << id << "\n";
    }

    generateSVGPlot("ALL_ecg_comparison.svg", allBinData, allMatData, "R-Peak Count (All Files)", rPeakCount);
    generateSVGPlot("ALL_ppg_max_comparison.svg", allBinData, allMatData, "PPG Max Amps Count (All Files)", ppgMaxCount);
    generateSVGPlot("ALL_ppg_min_comparison.svg", allBinData, allMatData, "PPG Min Amps Count (All Files)", ppgMinCount);
    generateSVGPlot("ALL_pairs_comparison.svg", allBinData, allMatData, "Detected Pairs (All Files)", pairCount);

    writeAllFilesSummary(results);

    return 0;
}