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
#include <sstream>

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
    bool ecg1_noisy = false;
    bool ecg2_noisy = false;
    bool ecg3_noisy = false;
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
    int noisyBins = 0;
    int cleanBins = 0;
    size_t totalBeatsBin = 0;
    size_t totalBeatsMat = 0;
    double timeSsdMean = 0.0;
    std::vector<double> perBinSsd;       // clean bins
    std::vector<double> perBinSsdNoisy;  // noisy bins
};

// ============================================================================
// Index-to-time conversion
// ============================================================================

std::vector<double> indicesToSeconds(const std::vector<double>& indices, double sampleRate) {
    std::vector<double> times(indices.size());
    for (size_t i = 0; i < indices.size(); i++) {
        times[i] = (indices[i] - 1.0) / sampleRate;
    }
    return times;
}

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

// ============================================================================
// Time-based R-peak comparison
// ============================================================================

struct TimedSSDResult {
    double ssd_seconds;
    size_t matchedCount;
    size_t unmatchedBin;
    size_t unmatchedMat;
    double binRMeanSec;
    double binRStdSec;
    double matRMeanSec;
    double matRStdSec;
};

TimedSSDResult computeTimedRPeakSSD(
    const ComparisonData& bin, double binSR,
    const ComparisonData& mat, double matSR,
    double toleranceSec = 0.150)
{
    auto binTimes = indicesToSeconds(bin.ecgRIndex, binSR);
    auto matTimes = indicesToSeconds(mat.ecgRIndex, matSR);

    TimedSSDResult result = { 0.0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0 };

    auto meanStd = [](const std::vector<double>& times) -> std::pair<double, double> {
        if (times.empty()) return { 0.0, 0.0 };
        double sum = 0.0, sumSq = 0.0;
        for (double t : times) { sum += t; sumSq += t * t; }
        double n = static_cast<double>(times.size());
        double mean = sum / n;
        double var = (n > 1) ? (sumSq - sum * sum / n) / (n - 1) : 0.0;
        return { mean, std::sqrt(std::max(0.0, var)) };
        };

    auto [binMean, binStd] = meanStd(binTimes);
    auto [matMean, matStd] = meanStd(matTimes);
    result.binRMeanSec = binMean;
    result.binRStdSec = binStd;
    result.matRMeanSec = matMean;
    result.matRStdSec = matStd;

    if (binTimes.empty() && matTimes.empty()) {
        return result;
    }
    if (binTimes.empty() || matTimes.empty()) {
        result.unmatchedBin = binTimes.size();
        result.unmatchedMat = matTimes.size();
        double duration = std::max(
            bin.ecgSignal.empty() ? 0.0 : (double)bin.ecgSignal.size() / binSR,
            mat.ecgSignal.empty() ? 0.0 : (double)mat.ecgSignal.size() / matSR);
        result.ssd_seconds = (result.unmatchedBin + result.unmatchedMat) * duration * duration;
        return result;
    }

    std::vector<bool> matUsed(matTimes.size(), false);
    double sumSqDiff = 0.0;
    size_t matched = 0;

    for (size_t i = 0; i < binTimes.size(); i++) {
        double bestDiff = toleranceSec + 1.0;
        int bestJ = -1;

        for (size_t j = 0; j < matTimes.size(); j++) {
            if (matUsed[j]) continue;
            double diff = std::abs(binTimes[i] - matTimes[j]);
            if (diff < bestDiff) {
                bestDiff = diff;
                bestJ = static_cast<int>(j);
            }
            if (matTimes[j] > binTimes[i] + toleranceSec) break;
        }

        if (bestJ >= 0 && bestDiff <= toleranceSec) {
            matUsed[bestJ] = true;
            double diff = binTimes[i] - matTimes[bestJ];
            sumSqDiff += diff * diff;
            matched++;
        }
    }

    size_t unmatchedMat = 0;
    for (bool u : matUsed) if (!u) unmatchedMat++;

    double duration = std::max(
        bin.ecgSignal.empty() ? 0.0 : (double)bin.ecgSignal.size() / binSR,
        mat.ecgSignal.empty() ? 0.0 : (double)mat.ecgSignal.size() / matSR);
    size_t unmatchedBin = binTimes.size() - matched;
    double penalty = (unmatchedBin + unmatchedMat) * duration * duration;

    result.ssd_seconds = sumSqDiff + penalty;
    result.matchedCount = matched;
    result.unmatchedBin = unmatchedBin;
    result.unmatchedMat = unmatchedMat;

    return result;
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

        readInt64VecAsDouble(file, r.ecgRIndex);
        readInt64VecAsDouble(file, r.ecgRIndex2);
        readInt64VecAsDouble(file, r.ecgRIndex3);
        readInt64VecAsDouble(file, r.ppgMaxAmps);
        readInt64VecAsDouble(file, r.ppgMinAmps);

        readDoubleVec(file, r.ppgSignal);
        readDoubleVec(file, r.ecgSignal);
        readDoubleVec(file, r.ecgSignal2);
        readDoubleVec(file, r.ecgSignal3);

        // Noise flags
        uint8_t n1 = 0, n2 = 0, n3 = 0;
        file.read(reinterpret_cast<char*>(&n1), 1);
        file.read(reinterpret_cast<char*>(&n2), 1);
        file.read(reinterpret_cast<char*>(&n3), 1);
        r.ecg1_noisy = (n1 != 0);
        r.ecg2_noisy = (n2 != 0);
        r.ecg3_noisy = (n3 != 0);

        // Pairs
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

        r.ecgSamplingRate = 2000.0;
        r.ppgSamplingRate = 2000.0;
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
        << "\" text-anchor=\"middle\" font-family=\"sans-serif\">C++ (2000 Hz) " << label << "</text>\n";
    svg << "<text x=\"25\" y=\"" << size / 2
        << "\" text-anchor=\"middle\" transform=\"rotate(-90 25," << size / 2
        << ")\" font-family=\"sans-serif\">MATLAB (256 Hz) " << label << "</text>\n";
    svg << "<text x=\"" << size / 2 << "\" y=\"40\" text-anchor=\"middle\" "
        << "font-weight=\"bold\" font-family=\"sans-serif\">" << label << " Comparison</text>\n";

    for (size_t i = 0; i < n; ++i) {
        double x = selector(binData[i]);
        double y = selector(matData[i]);
        bool noisy = binData[i].ecg1_noisy;
        const char* color = noisy ? "#cccccc" : (std::abs(x - y) < 0.0001) ? "#3498db" : "#e74c3c";
        double opacity = noisy ? 0.4 : 0.7;

        svg << "<circle cx=\"" << mapX(x) << "\" cy=\"" << mapY(y)
            << "\" r=\"3.5\" fill=\"" << color << "\" fill-opacity=\"" << opacity << "\">\n"
            << "  <title>Bin " << i << (noisy ? " (NOISY)" : "") << " C++: " << x << ", MATLAB: " << y << "</title>\n"
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

    // Separate clean vs noisy bins
    int noisyCount = 0;
    for (size_t i = 0; i < maxBins; ++i) {
        if (binData[i].ecg1_noisy) noisyCount++;
    }
    int cleanCount = static_cast<int>(maxBins) - noisyCount;

    // Stats on clean bins only
    std::vector<ComparisonData> cleanBin, cleanMat;
    for (size_t i = 0; i < maxBins; ++i) {
        if (!binData[i].ecg1_noisy) {
            cleanBin.push_back(binData[i]);
            cleanMat.push_back(matData[i]);
        }
    }

    auto rPeakProj = [](const ComparisonData& c) { return (double)c.ecgRIndex.size(); };
    auto ppgPeakProj = [](const ComparisonData& c) { return (double)c.ppgMaxAmps.size(); };
    auto ppgDipProj = [](const ComparisonData& c) { return (double)c.ppgMinAmps.size(); };
    auto pairProj = [](const ComparisonData& c) { return (double)c.pairs.size(); };

    DatasetStats binR = calculateStats(cleanBin, rPeakProj);
    DatasetStats matR = calculateStats(cleanMat, rPeakProj);
    DatasetStats binPPG = calculateStats(cleanBin, ppgPeakProj);
    DatasetStats matPPG = calculateStats(cleanMat, ppgPeakProj);
    DatasetStats binDip = calculateStats(cleanBin, ppgDipProj);
    DatasetStats matDip = calculateStats(cleanMat, ppgDipProj);
    DatasetStats binPrs = calculateStats(cleanBin, pairProj);
    DatasetStats matPrs = calculateStats(cleanMat, pairProj);

    size_t totalBeatsBin = 0, totalBeatsMat = 0;
    std::vector<TimedSSDResult> timedResults(maxBins);
    std::vector<double> perBinSsd;
    std::vector<double> perBinSsdNoisy;

    double sumTimeSsd = 0;

    for (size_t i = 0; i < maxBins; ++i) {
        double binSR = binData[i].ecgSamplingRate;
        double matSR = matData[i].ecgSamplingRate;
        if (binSR <= 0) binSR = 2000.0;
        if (matSR <= 0) matSR = 256.0;

        timedResults[i] = computeTimedRPeakSSD(binData[i], binSR, matData[i], matSR);

        if (!binData[i].ecg1_noisy) {
            totalBeatsBin += binData[i].ecgRIndex.size();
            totalBeatsMat += matData[i].ecgRIndex.size();
            sumTimeSsd += timedResults[i].ssd_seconds;
            perBinSsd.push_back(timedResults[i].ssd_seconds);
        }
        else {
            perBinSsdNoisy.push_back(timedResults[i].ssd_seconds);
        }
    }

    double nClean = static_cast<double>(cleanCount);
    double timeSsdMean = nClean > 0 ? sumTimeSsd / nClean : 0.0;

    out << "Data taken from MESA file " << id << "\n";

    out << "=== Summary ===\n";
    out << "Total bins: " << maxBins << "\n";
    out << "Noisy bins (excluded): " << noisyCount << "\n";
    out << "Clean bins (used for stats): " << cleanCount << "\n\n";

    out << "--- Time-based R-peak comparison ---\n";

    out << "\t\tTotal in C++\tTotal in MATLAB\tDiff\n";
    auto writeRow = [&](const char* label, DatasetStats b, DatasetStats m) {
        out << label << "\t" << b.total << "\t\t" << m.total
            << "\t\t" << (m.total - b.total) << "\n";
        };
    writeRow("R Peaks\t", binR, matR);
    writeRow("PPG Peaks", binPPG, matPPG);
    writeRow("PPG Dips", binDip, matDip);
    writeRow("Pairs\t", binPrs, matPrs);
    out << "\n";

    out << "bin#\tNoisy\tBeats(C++)\tBeats(MAT)\tMean Beat Time (C++)\tMean Beat Time (MAT)\tSTD Beat Time (C++)\tSTD Beat Time (MAT)\tSum of Squared Differences\n";

    for (size_t i = 0; i < maxBins; ++i) {
        const auto& tr = timedResults[i];
        bool noisy = binData[i].ecg1_noisy;
        out << std::setw(3) << std::setfill('0') << i << std::setfill(' ')
            << "\t" << (noisy ? "Y" : "N")
            << "\t" << binData[i].ecgRIndex.size()
            << "\t\t" << matData[i].ecgRIndex.size()
            << "\t\t" << std::fixed << std::setprecision(4) << tr.binRMeanSec
            << "\t\t\t" << tr.matRMeanSec
            << "\t\t\t" << tr.binRStdSec
            << "\t\t\t" << tr.matRStdSec
            << "\t\t\t" << std::setprecision(6) << tr.ssd_seconds
            << "\n";
    }

    return { id, binR.total, matR.total,
             static_cast<int>(binData.size()), static_cast<int>(matData.size()),
             noisyCount, cleanCount,
             totalBeatsBin, totalBeatsMat, timeSsdMean,
             std::move(perBinSsd), std::move(perBinSsdNoisy) };
}

// ============================================================================
// All-files summary: single SVG with one SSD histogram per file
// ============================================================================

void writeAllFilesSummary(const std::vector<FileResult>& results,
    int numBuckets = 50)
{
    if (results.empty()) return;

    int numFiles = static_cast<int>(results.size());
    int cols = 3;
    int rows = (numFiles + cols - 1) / cols;

    constexpr int cellW = 360, cellH = 260;
    constexpr int padL = 55, padR = 45, padT = 45, padB = 45;
    int chartW = cellW - padL - padR;
    int chartH = cellH - padT - padB;

    int svgW = cols * cellW;
    int svgH = rows * cellH + 50;

    std::ofstream svg("allfiles_ssd_histograms.svg");
    if (!svg.is_open()) return;

    svg << std::fixed;
    svg << "<svg width=\"" << svgW << "\" height=\"" << svgH
        << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
    svg << "<style>\n"
        << "  text { font-family: Consolas, 'Courier New', monospace; }\n"
        << "  .title { font-size: 16px; font-weight: bold; }\n"
        << "  .subtitle { font-size: 10px; font-weight: bold; }\n"
        << "  .stats { font-size: 12px; fill: #555; }\n"
        << "  .tick { font-size: 10px; }\n"
        << "  .axis-label { font-size: 12px; }\n"
        << "  .y-axis-label { font-size: 10px; }\n"
        << "  .x-tick { font-size: 10px; }\n"
        << "  .clipped { font-size: 10px; fill: #e74c3c; font-weight: bold; }\n"
        << "  .empty { font-size: 10px; fill: #999; }\n"
        << "</style>\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"#fcfcfc\"/>\n";

    svg << "<text x=\"" << svgW / 2 << "\" y=\"30\" text-anchor=\"middle\" "
        << "class=\"title\">"
        << "Per-Bin R-Peak SSD Histograms</text>\n";

    for (int f = 0; f < numFiles; ++f) {
        const auto& fr = results[f];
        int col = f % cols;
        int row = f / cols;
        int ox = col * cellW;
        int oy = row * cellH + 50;

        // Subtitle
        svg << "<text x=\"" << ox + cellW / 2 << "\" y=\"" << oy + 16
            << "\" text-anchor=\"middle\" class=\"subtitle\">"
            << fr.id << "</text>\n";

        // Merge clean + noisy for range calculation
        std::vector<double> allSsd;
        allSsd.insert(allSsd.end(), fr.perBinSsd.begin(), fr.perBinSsd.end());
        allSsd.insert(allSsd.end(), fr.perBinSsdNoisy.begin(), fr.perBinSsdNoisy.end());

        if (allSsd.empty()) {
            svg << "<text x=\"" << ox + cellW / 2 << "\" y=\"" << oy + cellH / 2
                << "\" text-anchor=\"middle\" class=\"empty\">"
                << "(no bins)</text>\n";
            continue;
        }

        // Stats (on clean bins only, as before)
        double sum = 0.0, sumSq = 0.0;
        double nClean = static_cast<double>(fr.perBinSsd.size());
        for (double v : fr.perBinSsd) { sum += v; sumSq += v * v; }
        double mean = nClean > 0 ? sum / nClean : 0.0;

        std::vector<double> sorted = fr.perBinSsd;
        std::sort(sorted.begin(), sorted.end());
        double median = 0.0;
        if (!sorted.empty()) {
            median = (sorted.size() % 2 == 0)
                ? (sorted[sorted.size() / 2 - 1] + sorted[sorted.size() / 2]) / 2.0
                : sorted[sorted.size() / 2];
        }

        svg << "<text x=\"" << ox + cellW / 2 << "\" y=\"" << oy + 28
            << "\" text-anchor=\"middle\" class=\"stats\">"
            << "clean=" << std::setprecision(0) << nClean
            << " noisy=" << fr.perBinSsdNoisy.size()
            << " mean=" << std::setprecision(4) << mean
            << " med=" << median << "</text>\n";

        // Histogram range from ALL bins (clean + noisy)
        std::sort(allSsd.begin(), allSsd.end());
        double minVal = allSsd.front();
        double maxVal = allSsd.back();
        if (maxVal - minVal < 1e-15) maxVal = minVal + 1.0;
        double bucketWidth = (maxVal - minVal) / numBuckets;

        // Separate counts for clean and noisy
        std::vector<int> cleanCounts(numBuckets, 0);
        std::vector<int> noisyCounts(numBuckets, 0);

        for (double v : fr.perBinSsd) {
            int idx = static_cast<int>((v - minVal) / bucketWidth);
            if (idx >= numBuckets) idx = numBuckets - 1;
            cleanCounts[idx]++;
        }
        for (double v : fr.perBinSsdNoisy) {
            int idx = static_cast<int>((v - minVal) / bucketWidth);
            if (idx >= numBuckets) idx = numBuckets - 1;
            noisyCounts[idx]++;
        }

        // Stacked total per bucket
        std::vector<int> totalCounts(numBuckets);
        for (int b = 0; b < numBuckets; ++b)
            totalCounts[b] = cleanCounts[b] + noisyCounts[b];

        int maxCount = *std::max_element(totalCounts.begin(), totalCounts.end());
        if (maxCount == 0) maxCount = 1;

        // Y-axis cutoff (same logic as before, but on total counts)
        std::vector<int> sortedCounts = totalCounts;
        std::sort(sortedCounts.begin(), sortedCounts.end());
        int yAxisMax = maxCount;
        for (int k = (int)sortedCounts.size() - 2; k >= 0; --k) {
            if (sortedCounts[k] < maxCount && sortedCounts[k] > 0) {
                yAxisMax = static_cast<int>(std::ceil(sortedCounts[k] * 1.2));
                break;
            }
        }
        if (yAxisMax <= 0) yAxisMax = maxCount;

        int cx = ox + padL;
        int cy = oy + padT;
        double barW = static_cast<double>(chartW) / numBuckets;

        // Axes
        svg << "<line x1=\"" << cx << "\" y1=\"" << cy
            << "\" x2=\"" << cx << "\" y2=\"" << cy + chartH
            << "\" stroke=\"#333\" stroke-width=\"1\"/>\n";
        svg << "<line x1=\"" << cx << "\" y1=\"" << cy + chartH
            << "\" x2=\"" << cx + chartW << "\" y2=\"" << cy + chartH
            << "\" stroke=\"#333\" stroke-width=\"1\"/>\n";

        // Y-axis label (rotated)
        svg << "<text x=\"" << ox + 22 << "\" y=\"" << cy + chartH / 2
            << "\" text-anchor=\"middle\" class=\"y-axis-label\""
            << " transform=\"rotate(-90 " << ox + 12 << "," << cy + chartH / 2 << ")\">"
            << "# bins with this SSD</text>\n";

        // Y-axis ticks
        for (int t = 0; t <= 3; ++t) {
            int yVal = static_cast<int>(std::round(yAxisMax * t / 3.0));
            double yPos = cy + chartH - (static_cast<double>(yVal) / yAxisMax) * chartH;
            svg << "<text x=\"" << cx - 4 << "\" y=\"" << yPos + 3
                << "\" text-anchor=\"end\" class=\"tick\">"
                << yVal << "</text>\n";
            if (t > 0) {
                svg << "<line x1=\"" << cx + 1 << "\" y1=\"" << yPos
                    << "\" x2=\"" << cx + chartW << "\" y2=\"" << yPos
                    << "\" stroke=\"#eee\" stroke-width=\"0.5\"/>\n";
            }
        }

        // Stacked bars: clean on bottom (blue), noisy on top (orange)
        for (int b = 0; b < numBuckets; ++b) {
            int total = totalCounts[b];
            if (total == 0) continue;

            bool clipped = total > yAxisMax;
            double bx = cx + b * barW;
            double rangeStart = minVal + b * bucketWidth;
            double rangeEnd = minVal + (b + 1) * bucketWidth;

            // Clean portion (bottom)
            double cleanDisplay = std::min(static_cast<double>(cleanCounts[b]),
                static_cast<double>(yAxisMax));
            double cleanBarH = (cleanDisplay / yAxisMax) * chartH;
            if (cleanCounts[b] > 0) {
                double by = cy + chartH - cleanBarH;
                svg << "<rect x=\"" << bx + 0.5 << "\" y=\"" << by
                    << "\" width=\"" << barW - 1 << "\" height=\"" << cleanBarH
                    << "\" fill=\"#3498db\" fill-opacity=\"0.8\">"
                    << "<title>[" << std::setprecision(4) << rangeStart
                    << ", " << rangeEnd << "): clean=" << cleanCounts[b]
                    << " noisy=" << noisyCounts[b] << "</title></rect>\n";
            }

            // Noisy portion (stacked on top of clean)
            if (noisyCounts[b] > 0) {
                double stackBase = cleanDisplay;
                double noisyDisplay = std::min(
                    static_cast<double>(noisyCounts[b]),
                    static_cast<double>(yAxisMax) - stackBase);
                if (noisyDisplay > 0) {
                    double noisyBarH = (noisyDisplay / yAxisMax) * chartH;
                    double by = cy + chartH - (stackBase / yAxisMax) * chartH - noisyBarH;
                    svg << "<rect x=\"" << bx + 0.5 << "\" y=\"" << by
                        << "\" width=\"" << barW - 1 << "\" height=\"" << noisyBarH
                        << "\" fill=\"#e67e22\" fill-opacity=\"0.7\">"
                        << "<title>[" << std::setprecision(4) << rangeStart
                        << ", " << rangeEnd << "): clean=" << cleanCounts[b]
                        << " noisy=" << noisyCounts[b] << "</title></rect>\n";
                }
            }

            // Annotate clipped bars
            if (clipped) {
                double topY = cy; // top of chart
                svg << "<text x=\"" << bx + barW / 2.0 << "\" y=\"" << topY - 2
                    << "\" text-anchor=\"middle\" class=\"clipped\">"
                    << total << "</text>\n";
            }
        }

        // X-axis labels
        for (int t = 0; t <= 2; ++t) {
            double val = minVal + (maxVal - minVal) * t / 2.0;
            double xPos = cx + static_cast<double>(chartW) * t / 2.0;
            svg << "<text x=\"" << xPos << "\" y=\"" << cy + chartH + 12
                << "\" text-anchor=\"middle\" class=\"x-tick\">"
                << std::setprecision(4) << val << "</text>\n";
        }

        // X label
        svg << "<text x=\"" << cx + chartW / 2 << "\" y=\"" << cy + chartH + 25
            << "\" text-anchor=\"middle\" class=\"axis-label\">"
            << "SSD Between my and Daniels R Peaks (s)</text>\n";
    }

    svg << "</svg>\n";
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

        std::cout << "  Done: " << id << "\n";
    }

    generateSVGPlot("ALL_ecg_comparison.svg", allBinData, allMatData, "R-Peak Count (All Files)", rPeakCount);
    generateSVGPlot("ALL_ppg_max_comparison.svg", allBinData, allMatData, "PPG Max Amps Count (All Files)", ppgMaxCount);
    generateSVGPlot("ALL_ppg_min_comparison.svg", allBinData, allMatData, "PPG Min Amps Count (All Files)", ppgMinCount);
    generateSVGPlot("ALL_pairs_comparison.svg", allBinData, allMatData, "Detected Pairs (All Files)", pairCount);

    writeAllFilesSummary(results);

    std::cout << "\nAll comparisons complete. See allfiles_ssd_histograms.svg for summary.\n";
    return 0;
}