/**
 * @file   compare_wave_bounds.cpp
 * @brief  Compare wave bound outputs between C++ (_wave_markings.bin)
 *         and MATLAB (_wave_data.mat). Only compares raw (ch1) ECG.
 *         Outputs per-file CSVs and a summary CSV.
 *
 * @author Mira Welner
 * @date   2026-03-31
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <functional>

#include <matio.h>

 // ============================================================================
 // Config
 // ============================================================================

static const std::string BIN_DIR =
"D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\4_wave_bound_files\\mesa_files\\";
static const std::string MAT_DIR =
"D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\4_wave_bound_files\\matlab\\";
static const std::string OUTPUT_DIR =
"D:\\USERS\\MiraWelner\\QTVI\\testing\\4_wave_bounds\\results\\";

static constexpr double BIN_SR = 2000.0;
static constexpr double MAT_SR = 256.0;

// ============================================================================
// Data structures
// ============================================================================

struct BinData {
    std::vector<size_t> ecgRIndex;      // ch1 raw only
    std::vector<size_t> ppgMaxAmps;
    std::vector<size_t> ppgMinAmps;
    std::vector<double> ppgSignal;
    std::vector<double> ecgSignal;
    std::vector<std::vector<double>> pairs;
    std::vector<std::pair<uint64_t, uint64_t>> ppg_bin_indexs;
    std::vector<std::pair<uint64_t, uint64_t>> ecg_bin_indexs;
    bool noisy = false;
    bool bad_segment = false;
};

struct PerBinResult {
    int bin;
    bool noisy;
    bool bad_segment;
    size_t rPeaksBin;
    size_t rPeaksMat;
    size_t ppgMaxBin;
    size_t ppgMaxMat;
    size_t ppgMinBin;
    size_t ppgMinMat;
    size_t pairsBin;
    size_t pairsMat;
    double ssd_seconds;
    size_t matched;
    size_t unmatchedBin;
    size_t unmatchedMat;
};

struct FileResult {
    std::string id;
    int totalBins;
    int noisyBins;
    int cleanBins;
    std::vector<PerBinResult> bins;
};

// ============================================================================
// Statistics
// ============================================================================

struct Stats {
    double mean = 0.0;
    double median = 0.0;
    double q1 = 0.0;
    double q3 = 0.0;
    size_t count = 0;
};

Stats computeStats(const std::vector<double>& vals) {
    Stats s;
    if (vals.empty()) return s;
    std::vector<double> sorted = vals;
    std::sort(sorted.begin(), sorted.end());
    s.count = sorted.size();
    double sum = 0.0;
    for (double v : sorted) sum += v;
    s.mean = sum / s.count;
    s.median = (s.count % 2 == 0)
        ? (sorted[s.count / 2 - 1] + sorted[s.count / 2]) / 2.0
        : sorted[s.count / 2];
    s.q1 = sorted[s.count / 4];
    s.q3 = sorted[s.count * 3 / 4];
    return s;
}

// ============================================================================
// Time-based R-peak SSD
// ============================================================================

struct SSDResult {
    double ssd_seconds = 0.0;
    size_t matched = 0;
    size_t unmatchedBin = 0;
    size_t unmatchedMat = 0;
};

std::vector<double> indicesToSeconds(const std::vector<size_t>& indices, double sr) {
    std::vector<double> times(indices.size());
    for (size_t i = 0; i < indices.size(); ++i)
        times[i] = static_cast<double>(indices[i]) / sr;
    return times;
}

SSDResult computeSSD(const std::vector<size_t>& binIdx, double binSR,
    const std::vector<size_t>& matIdx, double matSR,
    double toleranceSec = 0.150)
{
    auto binTimes = indicesToSeconds(binIdx, binSR);
    auto matTimes = indicesToSeconds(matIdx, matSR);
    SSDResult res;

    if (binTimes.empty() && matTimes.empty()) return res;
    if (binTimes.empty() || matTimes.empty()) {
        res.unmatchedBin = binTimes.size();
        res.unmatchedMat = matTimes.size();
        return res;
    }

    std::vector<bool> matUsed(matTimes.size(), false);
    double sumSqDiff = 0.0;

    for (size_t i = 0; i < binTimes.size(); ++i) {
        double bestDiff = toleranceSec + 1.0;
        int bestJ = -1;
        for (size_t j = 0; j < matTimes.size(); ++j) {
            if (matUsed[j]) continue;
            double diff = std::abs(binTimes[i] - matTimes[j]);
            if (diff < bestDiff) { bestDiff = diff; bestJ = static_cast<int>(j); }
            if (matTimes[j] > binTimes[i] + toleranceSec) break;
        }
        if (bestJ >= 0 && bestDiff <= toleranceSec) {
            matUsed[bestJ] = true;
            double diff = binTimes[i] - matTimes[bestJ];
            sumSqDiff += diff * diff;
            res.matched++;
        }
    }

    for (bool u : matUsed) if (!u) res.unmatchedMat++;
    res.unmatchedBin = binTimes.size() - res.matched;
    res.ssd_seconds = sumSqDiff;
    return res;
}

// ============================================================================
// Read C++ _wave_markings.bin
// ============================================================================

std::vector<BinData> readBinFile(const std::string& path) {
    std::vector<BinData> results;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return results;

    constexpr uint64_t MAX_SANE = 50000000;

    uint64_t numBins = 0;
    file.read(reinterpret_cast<char*>(&numBins), 8);
    if (numBins > 100000) return results;
    results.resize(numBins);

    auto readIdx = [&](std::vector<size_t>& v) {
        uint64_t sz = 0;
        if (!file.read(reinterpret_cast<char*>(&sz), 8)) { v.clear(); return; }
        if (sz > MAX_SANE) { v.clear(); return; }
        v.resize(sz);
        if (sz > 0) {
            std::vector<uint64_t> tmp(sz);
            file.read(reinterpret_cast<char*>(tmp.data()), sz * 8);
            for (uint64_t j = 0; j < sz; ++j)
                v[j] = static_cast<size_t>(tmp[j] > 0 ? tmp[j] - 1 : 0);
        }
        };

    auto skipIdx = [&]() {
        uint64_t sz = 0;
        if (!file.read(reinterpret_cast<char*>(&sz), 8)) return;
        if (sz > MAX_SANE) return;
        file.seekg(sz * 8, std::ios::cur);
        };

    auto readSignal = [&](std::vector<double>& sig) {
        uint64_t sz = 0;
        if (!file.read(reinterpret_cast<char*>(&sz), 8)) { sig.clear(); return; }
        if (sz > MAX_SANE) { sig.clear(); return; }
        sig.resize(sz);
        if (sz > 0) file.read(reinterpret_cast<char*>(sig.data()), sz * 8);
        };

    auto skipSignal = [&]() {
        uint64_t sz = 0;
        if (!file.read(reinterpret_cast<char*>(&sz), 8)) return;
        if (sz > MAX_SANE) return;
        file.seekg(sz * 8, std::ios::cur);
        };

    auto readPairVec = [&](std::vector<std::pair<uint64_t, uint64_t>>& v) {
        uint64_t sz = 0;
        if (!file.read(reinterpret_cast<char*>(&sz), 8)) { v.clear(); return; }
        if (sz > MAX_SANE) { v.clear(); return; }
        v.resize(sz);
        if (sz > 0) file.read(reinterpret_cast<char*>(v.data()), sz * 16);
        };

    for (uint64_t i = 0; i < numBins; ++i) {
        auto& bin = results[i];

        // 9 R-peak index arrays: ch1(raw,sq,abs), ch2(raw,sq,abs), ch3(raw,sq,abs)
        readIdx(bin.ecgRIndex);  // ch1 raw
        skipIdx();               // ch1 squared
        skipIdx();               // ch1 absval
        skipIdx();               // ch2 raw
        skipIdx();               // ch2 squared
        skipIdx();               // ch2 absval
        skipIdx();               // ch3 raw
        skipIdx();               // ch3 squared
        skipIdx();               // ch3 absval

        // 2 PPG index arrays
        readIdx(bin.ppgMaxAmps);
        readIdx(bin.ppgMinAmps);

        // 4 raw signals
        readSignal(bin.ppgSignal);
        readSignal(bin.ecgSignal);
        skipSignal();  // ecgSignal2
        skipSignal();  // ecgSignal3

        // 6 preprocessed signals
        skipSignal();  // ch1 squared
        skipSignal();  // ch1 absval
        skipSignal();  // ch2 squared
        skipSignal();  // ch2 absval
        skipSignal();  // ch3 squared
        skipSignal();  // ch3 absval

        // 9 noise flags
        uint8_t flags[9] = {};
        file.read(reinterpret_cast<char*>(flags), 9);
        bin.noisy = (flags[0] != 0);  // ch1 raw noisy

        // Pairs
        uint64_t numPairs = 0;
        file.read(reinterpret_cast<char*>(&numPairs), 8);
        if (numPairs < MAX_SANE && numPairs > 0) {
            bin.pairs.resize(numPairs);
            std::vector<int64_t> pairBuf(numPairs * 2);
            file.read(reinterpret_cast<char*>(pairBuf.data()), numPairs * 16);
            for (uint64_t j = 0; j < numPairs; ++j) {
                bin.pairs[j].resize(2);
                bin.pairs[j][0] = (pairBuf[j * 2] == -1)
                    ? -1.0 : static_cast<double>(pairBuf[j * 2] - 1);
                bin.pairs[j][1] = (pairBuf[j * 2 + 1] == -1)
                    ? -1.0 : static_cast<double>(pairBuf[j * 2 + 1] - 1);
            }
        }

        // ppg/ecg bin index pairs
        readPairVec(bin.ppg_bin_indexs);
        readPairVec(bin.ecg_bin_indexs);

        // bad_segment
        bool has_any_ecg = !bin.ecgSignal.empty();
        bool has_any_peaks = !bin.ecgRIndex.empty();
        bin.bad_segment = !has_any_ecg || (!has_any_peaks && bin.ppgMinAmps.empty());
    }
    return results;
}

// ============================================================================
// Read MATLAB _wave_data.mat
// ============================================================================

std::vector<double> getMatField(matvar_t* cell, const char* name) {
    std::vector<double> vec;
    matvar_t* field = Mat_VarGetStructFieldByName(cell, name, 0);
    if (!field || !field->data) return vec;
    size_t n = 1;
    for (int d = 0; d < field->rank; ++d) n *= field->dims[d];
    vec.reserve(n);
    for (size_t j = 0; j < n; ++j) {
        if (field->class_type == MAT_C_DOUBLE)
            vec.push_back(static_cast<double*>(field->data)[j]);
        else if (field->class_type == MAT_C_UINT64)
            vec.push_back(static_cast<double>(static_cast<uint64_t*>(field->data)[j]));
        else if (field->class_type == MAT_C_INT64)
            vec.push_back(static_cast<double>(static_cast<int64_t*>(field->data)[j]));
        else if (field->class_type == MAT_C_SINGLE)
            vec.push_back(static_cast<double>(static_cast<float*>(field->data)[j]));
        else if (field->class_type == MAT_C_INT32)
            vec.push_back(static_cast<double>(static_cast<int32_t*>(field->data)[j]));
        else if (field->class_type == MAT_C_UINT32)
            vec.push_back(static_cast<double>(static_cast<uint32_t*>(field->data)[j]));
    }
    return vec;
}

struct MatBinData {
    std::vector<size_t> ecgRIndex;
    std::vector<size_t> ppgMaxAmps;
    std::vector<size_t> ppgMinAmps;
    size_t pairCount = 0;
    double ecgSamplingRate = 0.0;
    double ppgSamplingRate = 0.0;
};

std::vector<MatBinData> readMatFile(const std::string& path) {
    std::vector<MatBinData> results;
    mat_t* matfp = Mat_Open(path.c_str(), MAT_ACC_RDONLY);
    if (!matfp) return results;

    matvar_t* wave_data = Mat_VarRead(matfp, "wave_data");
    if (!wave_data) { Mat_Close(matfp); return results; }

    size_t numBins = wave_data->dims[0] * wave_data->dims[1];
    results.resize(numBins);

    for (size_t i = 0; i < numBins; ++i) {
        matvar_t* cell = Mat_VarGetCell(wave_data, i);
        if (!cell) continue;
        auto& r = results[i];

        auto eRate = getMatField(cell, "ecgSamplingRate");
        r.ecgSamplingRate = eRate.empty() ? MAT_SR : eRate[0];
        auto pRate = getMatField(cell, "ppgSamplingRate");
        r.ppgSamplingRate = pRate.empty() ? MAT_SR : pRate[0];

        auto rIdx = getMatField(cell, "ecgRIndex");
        r.ecgRIndex.reserve(rIdx.size());
        for (double v : rIdx) r.ecgRIndex.push_back(static_cast<size_t>(v));

        auto maxA = getMatField(cell, "ppgMaxAmps");
        r.ppgMaxAmps.reserve(maxA.size());
        for (double v : maxA) r.ppgMaxAmps.push_back(static_cast<size_t>(v));

        auto minA = getMatField(cell, "ppgMinAmps");
        r.ppgMinAmps.reserve(minA.size());
        for (double v : minA) r.ppgMinAmps.push_back(static_cast<size_t>(v));

        matvar_t* pair_var = Mat_VarGetStructFieldByName(cell, "pairs", 0);
        if (pair_var && pair_var->data && pair_var->rank == 2)
            r.pairCount = pair_var->dims[0];
    }

    Mat_VarFree(wave_data);
    Mat_Close(matfp);
    return results;
}

// ============================================================================
// Compare one file
// ============================================================================

FileResult compareFile(const std::string& id,
    const std::vector<BinData>& binData,
    const std::vector<MatBinData>& matData) {
    FileResult fr;
    fr.id = id;
    fr.totalBins = static_cast<int>(std::min(binData.size(), matData.size()));
    fr.noisyBins = 0;
    fr.cleanBins = 0;

    for (int i = 0; i < fr.totalBins; ++i) {
        const auto& b = binData[i];
        const auto& m = matData[i];

        PerBinResult pbr;
        pbr.bin = i;
        pbr.noisy = b.noisy;
        pbr.bad_segment = b.bad_segment;
        pbr.rPeaksBin = b.ecgRIndex.size();
        pbr.rPeaksMat = m.ecgRIndex.size();
        pbr.ppgMaxBin = b.ppgMaxAmps.size();
        pbr.ppgMaxMat = m.ppgMaxAmps.size();
        pbr.ppgMinBin = b.ppgMinAmps.size();
        pbr.ppgMinMat = m.ppgMinAmps.size();
        pbr.pairsBin = b.pairs.size();
        pbr.pairsMat = m.pairCount;

        double binSR = BIN_SR;
        double matSR = m.ecgSamplingRate > 0 ? m.ecgSamplingRate : MAT_SR;
        auto ssd = computeSSD(b.ecgRIndex, binSR, m.ecgRIndex, matSR);
        pbr.ssd_seconds = ssd.ssd_seconds;
        pbr.matched = ssd.matched;
        pbr.unmatchedBin = ssd.unmatchedBin;
        pbr.unmatchedMat = ssd.unmatchedMat;

        if (b.noisy) fr.noisyBins++;
        else fr.cleanBins++;

        fr.bins.push_back(pbr);
    }
    return fr;
}

// ============================================================================
// CSV output
// ============================================================================

void writePerFileCSV(const std::string& path, const FileResult& fr) {
    std::ofstream out(path);
    if (!out.is_open()) return;

    // Header
    out << "bin,noisy,bad_segment,"
        << "rpeaks_cpp,rpeaks_mat,rpeaks_diff,"
        << "ppg_max_cpp,ppg_max_mat,ppg_max_diff,"
        << "ppg_min_cpp,ppg_min_mat,ppg_min_diff,"
        << "pairs_cpp,pairs_mat,pairs_diff,"
        << "ssd_seconds,matched,unmatched_cpp,unmatched_mat\n";

    // Stats rows (clean bins only)
    std::vector<double> ssdVals, rDiffVals, ppgMaxDiffVals, ppgMinDiffVals, pairDiffVals;
    for (const auto& b : fr.bins) {
        if (b.noisy) continue;
        ssdVals.push_back(b.ssd_seconds);
        rDiffVals.push_back(static_cast<double>(b.rPeaksBin) - b.rPeaksMat);
        ppgMaxDiffVals.push_back(static_cast<double>(b.ppgMaxBin) - b.ppgMaxMat);
        ppgMinDiffVals.push_back(static_cast<double>(b.ppgMinBin) - b.ppgMinMat);
        pairDiffVals.push_back(static_cast<double>(b.pairsBin) - b.pairsMat);
    }

    auto ssdS = computeStats(ssdVals);
    auto rS = computeStats(rDiffVals);
    auto maxS = computeStats(ppgMaxDiffVals);
    auto minS = computeStats(ppgMinDiffVals);
    auto pairS = computeStats(pairDiffVals);

    out << std::fixed << std::setprecision(6);

    auto writeStatRow = [&](const char* label) {
        std::function<double(const Stats&)> getter;
        if (std::string(label) == "MEAN") getter = [](const Stats& s) { return s.mean; };
        else if (std::string(label) == "MEDIAN") getter = [](const Stats& s) { return s.median; };
        else if (std::string(label) == "Q1") getter = [](const Stats& s) { return s.q1; };
        else getter = [](const Stats& s) { return s.q3; };

        out << label << ",,,"
            << ",," << getter(rS) << ","
            << ",," << getter(maxS) << ","
            << ",," << getter(minS) << ","
            << ",," << getter(pairS) << ","
            << getter(ssdS) << ",,,\n";
        };

    writeStatRow("MEAN");
    writeStatRow("MEDIAN");
    writeStatRow("Q1");
    writeStatRow("Q3");

    // Per-bin rows
    for (const auto& b : fr.bins) {
        out << b.bin << ","
            << (b.noisy ? "Y" : "N") << ","
            << (b.bad_segment ? "Y" : "N") << ","
            << b.rPeaksBin << "," << b.rPeaksMat << ","
            << (static_cast<int>(b.rPeaksBin) - static_cast<int>(b.rPeaksMat)) << ","
            << b.ppgMaxBin << "," << b.ppgMaxMat << ","
            << (static_cast<int>(b.ppgMaxBin) - static_cast<int>(b.ppgMaxMat)) << ","
            << b.ppgMinBin << "," << b.ppgMinMat << ","
            << (static_cast<int>(b.ppgMinBin) - static_cast<int>(b.ppgMinMat)) << ","
            << b.pairsBin << "," << b.pairsMat << ","
            << (static_cast<int>(b.pairsBin) - static_cast<int>(b.pairsMat)) << ","
            << b.ssd_seconds << ","
            << b.matched << "," << b.unmatchedBin << "," << b.unmatchedMat << "\n";
    }
}

void writeSummaryCSV(const std::string& path, const std::vector<FileResult>& results) {
    std::ofstream out(path);
    if (!out.is_open()) return;

    out << "subject,total_bins,clean_bins,noisy_bins,"
        << "ssd_mean,ssd_median,ssd_q1,ssd_q3,"
        << "rpeak_diff_mean,rpeak_diff_median,rpeak_diff_q1,rpeak_diff_q3,"
        << "ppg_min_diff_mean,ppg_min_diff_median,ppg_min_diff_q1,ppg_min_diff_q3,"
        << "pairs_diff_mean,pairs_diff_median,pairs_diff_q1,pairs_diff_q3\n";

    out << std::fixed << std::setprecision(6);

    // Collect all clean-bin stats across all files for grand summary
    std::vector<double> allSsd, allRDiff, allMinDiff, allPairDiff;

    for (const auto& fr : results) {
        std::vector<double> ssdVals, rDiffVals, minDiffVals, pairDiffVals;
        for (const auto& b : fr.bins) {
            if (b.noisy) continue;
            ssdVals.push_back(b.ssd_seconds);
            rDiffVals.push_back(static_cast<double>(b.rPeaksBin) - b.rPeaksMat);
            minDiffVals.push_back(static_cast<double>(b.ppgMinBin) - b.ppgMinMat);
            pairDiffVals.push_back(static_cast<double>(b.pairsBin) - b.pairsMat);
        }

        allSsd.insert(allSsd.end(), ssdVals.begin(), ssdVals.end());
        allRDiff.insert(allRDiff.end(), rDiffVals.begin(), rDiffVals.end());
        allMinDiff.insert(allMinDiff.end(), minDiffVals.begin(), minDiffVals.end());
        allPairDiff.insert(allPairDiff.end(), pairDiffVals.begin(), pairDiffVals.end());

        auto ssdS = computeStats(ssdVals);
        auto rS = computeStats(rDiffVals);
        auto minS = computeStats(minDiffVals);
        auto pairS = computeStats(pairDiffVals);

        out << fr.id << "," << fr.totalBins << "," << fr.cleanBins << "," << fr.noisyBins << ","
            << ssdS.mean << "," << ssdS.median << "," << ssdS.q1 << "," << ssdS.q3 << ","
            << rS.mean << "," << rS.median << "," << rS.q1 << "," << rS.q3 << ","
            << minS.mean << "," << minS.median << "," << minS.q1 << "," << minS.q3 << ","
            << pairS.mean << "," << pairS.median << "," << pairS.q1 << "," << pairS.q3 << "\n";
    }

    // Grand average row
    auto gSsd = computeStats(allSsd);
    auto gR = computeStats(allRDiff);
    auto gMin = computeStats(allMinDiff);
    auto gPair = computeStats(allPairDiff);

    out << "OVERALL,,,,"
        << gSsd.mean << "," << gSsd.median << "," << gSsd.q1 << "," << gSsd.q3 << ","
        << gR.mean << "," << gR.median << "," << gR.q1 << "," << gR.q3 << ","
        << gMin.mean << "," << gMin.median << "," << gMin.q1 << "," << gMin.q3 << ","
        << gPair.mean << "," << gPair.median << "," << gPair.q1 << "," << gPair.q3 << "\n";
}

// ============================================================================
// SVG generation
// ============================================================================

void generateSVGPlot(const std::string& filename,
    const std::vector<BinData>& binData,
    const std::vector<MatBinData>& matData,
    const std::string& label,
    std::function<double(const BinData&)> binSelector,
    std::function<double(const MatBinData&)> matSelector) {

    size_t n = std::min(binData.size(), matData.size());
    if (n == 0) return;

    double maxVal = 0.0;
    for (size_t i = 0; i < n; ++i)
        maxVal = std::max({ maxVal, binSelector(binData[i]), matSelector(matData[i]) });
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
        double x = binSelector(binData[i]);
        double y = matSelector(matData[i]);
        bool noisy = binData[i].noisy;
        const char* color = noisy ? "#cccccc" : (std::abs(x - y) < 0.0001) ? "#3498db" : "#e74c3c";
        double opacity = noisy ? 0.4 : 0.7;

        svg << "<circle cx=\"" << mapX(x) << "\" cy=\"" << mapY(y)
            << "\" r=\"3.5\" fill=\"" << color << "\" fill-opacity=\"" << opacity << "\">\n"
            << "  <title>Bin " << i << (noisy ? " (NOISY)" : "") << " C++: " << x << ", MATLAB: " << y << "</title>\n"
            << "</circle>\n";
    }

    svg << "</svg>";
}

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

    std::string svgPath = OUTPUT_DIR + "allfiles_ssd_histograms.svg";
    std::ofstream svg(svgPath);
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

        svg << "<text x=\"" << ox + cellW / 2 << "\" y=\"" << oy + 16
            << "\" text-anchor=\"middle\" class=\"subtitle\">"
            << fr.id << "</text>\n";

        // Collect clean SSD values
        std::vector<double> cleanSsd, noisySsd;
        for (const auto& b : fr.bins) {
            if (b.noisy) noisySsd.push_back(b.ssd_seconds);
            else cleanSsd.push_back(b.ssd_seconds);
        }

        std::vector<double> allSsd;
        allSsd.insert(allSsd.end(), cleanSsd.begin(), cleanSsd.end());
        allSsd.insert(allSsd.end(), noisySsd.begin(), noisySsd.end());

        if (allSsd.empty()) {
            svg << "<text x=\"" << ox + cellW / 2 << "\" y=\"" << oy + cellH / 2
                << "\" text-anchor=\"middle\" class=\"empty\">"
                << "(no bins)</text>\n";
            continue;
        }

        auto s = computeStats(cleanSsd);
        svg << "<text x=\"" << ox + cellW / 2 << "\" y=\"" << oy + 28
            << "\" text-anchor=\"middle\" class=\"stats\">"
            << "clean=" << std::setprecision(0) << cleanSsd.size()
            << " noisy=" << noisySsd.size()
            << " mean=" << std::setprecision(4) << s.mean
            << " med=" << s.median << "</text>\n";

        std::sort(allSsd.begin(), allSsd.end());
        double minVal = allSsd.front();
        double maxVal = allSsd.back();
        if (maxVal - minVal < 1e-15) maxVal = minVal + 1.0;
        double bucketWidth = (maxVal - minVal) / numBuckets;

        std::vector<int> cleanCounts(numBuckets, 0);
        std::vector<int> noisyCounts(numBuckets, 0);

        for (double v : cleanSsd) {
            int idx = static_cast<int>((v - minVal) / bucketWidth);
            if (idx >= numBuckets) idx = numBuckets - 1;
            cleanCounts[idx]++;
        }
        for (double v : noisySsd) {
            int idx = static_cast<int>((v - minVal) / bucketWidth);
            if (idx >= numBuckets) idx = numBuckets - 1;
            noisyCounts[idx]++;
        }

        std::vector<int> totalCounts(numBuckets);
        for (int b = 0; b < numBuckets; ++b)
            totalCounts[b] = cleanCounts[b] + noisyCounts[b];

        int maxCount = *std::max_element(totalCounts.begin(), totalCounts.end());
        if (maxCount == 0) maxCount = 1;

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

        svg << "<line x1=\"" << cx << "\" y1=\"" << cy
            << "\" x2=\"" << cx << "\" y2=\"" << cy + chartH
            << "\" stroke=\"#333\" stroke-width=\"1\"/>\n";
        svg << "<line x1=\"" << cx << "\" y1=\"" << cy + chartH
            << "\" x2=\"" << cx + chartW << "\" y2=\"" << cy + chartH
            << "\" stroke=\"#333\" stroke-width=\"1\"/>\n";

        svg << "<text x=\"" << ox + 22 << "\" y=\"" << cy + chartH / 2
            << "\" text-anchor=\"middle\" class=\"y-axis-label\""
            << " transform=\"rotate(-90 " << ox + 12 << "," << cy + chartH / 2 << ")\">"
            << "# bins</text>\n";

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

        for (int b = 0; b < numBuckets; ++b) {
            int total = totalCounts[b];
            if (total == 0) continue;

            bool clipped = total > yAxisMax;
            double bx = cx + b * barW;
            double rangeStart = minVal + b * bucketWidth;
            double rangeEnd = minVal + (b + 1) * bucketWidth;

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

            if (clipped) {
                double topY = cy;
                svg << "<text x=\"" << bx + barW / 2.0 << "\" y=\"" << topY - 2
                    << "\" text-anchor=\"middle\" class=\"clipped\">"
                    << total << "</text>\n";
            }
        }

        for (int t = 0; t <= 2; ++t) {
            double val = minVal + (maxVal - minVal) * t / 2.0;
            double xPos = cx + static_cast<double>(chartW) * t / 2.0;
            svg << "<text x=\"" << xPos << "\" y=\"" << cy + chartH + 12
                << "\" text-anchor=\"middle\" class=\"x-tick\">"
                << std::setprecision(4) << val << "</text>\n";
        }

        svg << "<text x=\"" << cx + chartW / 2 << "\" y=\"" << cy + chartH + 25
            << "\" text-anchor=\"middle\" class=\"axis-label\">"
            << "SSD (s)</text>\n";
    }

    svg << "</svg>\n";
}

// ============================================================================
// Main
// ============================================================================

int main() {
    std::filesystem::create_directories(OUTPUT_DIR);

    // Find matching files
    std::vector<std::string> ids;
    for (const auto& entry : std::filesystem::directory_iterator(BIN_DIR)) {
        if (entry.path().extension() != ".bin") continue;
        std::string stem = entry.path().stem().string();
        auto pos = stem.find("_wave_markings");
        if (pos == std::string::npos) continue;
        std::string id = stem.substr(0, pos);

        std::string matPath = MAT_DIR + id + "_wave_data.mat";
        if (std::filesystem::exists(matPath))
            ids.push_back(id);
    }
    std::sort(ids.begin(), ids.end());

    std::cout << "======================================================================\n";
    std::cout << "WAVE BOUNDS COMPARISON: C++ (2000 Hz) vs MATLAB (256 Hz)\n";
    std::cout << "======================================================================\n";
    std::cout << "C++ dir:    " << BIN_DIR << "\n";
    std::cout << "MATLAB dir: " << MAT_DIR << "\n";
    std::cout << "Output dir: " << OUTPUT_DIR << "\n";
    std::cout << "Matched:    " << ids.size() << " files\n\n";

    std::vector<FileResult> results;
    std::vector<BinData> allBinData;
    std::vector<MatBinData> allMatData;

    for (const auto& id : ids) {
        std::cout << "  " << id << "..." << std::flush;

        std::string binPath = BIN_DIR + id + "_wave_markings.bin";
        std::string matPath = MAT_DIR + id + "_wave_data.mat";

        auto binData = readBinFile(binPath);
        auto matData = readMatFile(matPath);

        if (binData.empty() || matData.empty()) {
            std::cout << " SKIPPED (empty data)\n";
            continue;
        }

        auto fr = compareFile(id, binData, matData);
        results.push_back(fr);

        // Write per-file CSV
        std::string perFilePath = OUTPUT_DIR + id + "_wave_comparison.csv";
        writePerFileCSV(perFilePath, fr);

        // Per-file SVGs
        auto rPeakBin = [](const BinData& d) { return (double)d.ecgRIndex.size(); };
        auto rPeakMat = [](const MatBinData& d) { return (double)d.ecgRIndex.size(); };
        auto ppgMaxBin = [](const BinData& d) { return (double)d.ppgMaxAmps.size(); };
        auto ppgMaxMat = [](const MatBinData& d) { return (double)d.ppgMaxAmps.size(); };
        auto ppgMinBin = [](const BinData& d) { return (double)d.ppgMinAmps.size(); };
        auto ppgMinMat = [](const MatBinData& d) { return (double)d.ppgMinAmps.size(); };
        auto pairBin = [](const BinData& d) { return (double)d.pairs.size(); };
        auto pairMat = [](const MatBinData& d) { return (double)d.pairCount; };

        generateSVGPlot(OUTPUT_DIR + id + "_ecg_comparison.svg", binData, matData, "R-Peak Count", rPeakBin, rPeakMat);
        generateSVGPlot(OUTPUT_DIR + id + "_ppg_max_comparison.svg", binData, matData, "PPG Max Amps Count", ppgMaxBin, ppgMaxMat);
        generateSVGPlot(OUTPUT_DIR + id + "_ppg_min_comparison.svg", binData, matData, "PPG Min Amps Count", ppgMinBin, ppgMinMat);
        generateSVGPlot(OUTPUT_DIR + id + "_pairs_comparison.svg", binData, matData, "Detected Pairs", pairBin, pairMat);

        // Accumulate for all-files SVGs
        size_t maxBins = std::min(binData.size(), matData.size());
        for (size_t i = 0; i < maxBins; ++i) {
            allBinData.push_back(binData[i]);
            allMatData.push_back(matData[i]);
        }

        // Print summary line
        std::vector<double> cleanSsd;
        for (const auto& b : fr.bins)
            if (!b.noisy) cleanSsd.push_back(b.ssd_seconds);
        auto s = computeStats(cleanSsd);

        std::cout << " " << fr.totalBins << " bins"
            << " (clean=" << fr.cleanBins << " noisy=" << fr.noisyBins << ")"
            << " SSD[mean=" << std::fixed << std::setprecision(4) << s.mean
            << " med=" << s.median
            << " IQR=(" << s.q1 << "," << s.q3 << ")]\n";
    }

    // All-files SVGs
    {
        auto rPeakBin = [](const BinData& d) { return (double)d.ecgRIndex.size(); };
        auto rPeakMat = [](const MatBinData& d) { return (double)d.ecgRIndex.size(); };
        auto ppgMaxBin = [](const BinData& d) { return (double)d.ppgMaxAmps.size(); };
        auto ppgMaxMat = [](const MatBinData& d) { return (double)d.ppgMaxAmps.size(); };
        auto ppgMinBin = [](const BinData& d) { return (double)d.ppgMinAmps.size(); };
        auto ppgMinMat = [](const MatBinData& d) { return (double)d.ppgMinAmps.size(); };
        auto pairBin = [](const BinData& d) { return (double)d.pairs.size(); };
        auto pairMat = [](const MatBinData& d) { return (double)d.pairCount; };

        generateSVGPlot(OUTPUT_DIR + "ALL_ecg_comparison.svg", allBinData, allMatData, "R-Peak Count (All Files)", rPeakBin, rPeakMat);
        generateSVGPlot(OUTPUT_DIR + "ALL_ppg_max_comparison.svg", allBinData, allMatData, "PPG Max Amps Count (All Files)", ppgMaxBin, ppgMaxMat);
        generateSVGPlot(OUTPUT_DIR + "ALL_ppg_min_comparison.svg", allBinData, allMatData, "PPG Min Amps Count (All Files)", ppgMinBin, ppgMinMat);
        generateSVGPlot(OUTPUT_DIR + "ALL_pairs_comparison.svg", allBinData, allMatData, "Detected Pairs (All Files)", pairBin, pairMat);
    }

    // SSD histogram summary SVG
    writeAllFilesSummary(results);

    // Write summary CSV
    std::string summaryPath = OUTPUT_DIR + "summary.csv";
    writeSummaryCSV(summaryPath, results);

    // Print summary table
    std::cout << "\n======================================================================\n";
    std::cout << "SUMMARY\n";
    std::cout << "======================================================================\n";
    std::cout << std::left << std::setw(25) << "Subject"
        << std::right << std::setw(6) << "Bins"
        << std::setw(7) << "Clean"
        << std::setw(10) << "SSD mean"
        << std::setw(10) << "SSD med"
        << std::setw(10) << "SSD Q1"
        << std::setw(10) << "SSD Q3" << "\n";
    std::cout << std::string(78, '-') << "\n";

    for (const auto& fr : results) {
        std::vector<double> cleanSsd;
        for (const auto& b : fr.bins)
            if (!b.noisy) cleanSsd.push_back(b.ssd_seconds);
        auto s = computeStats(cleanSsd);

        std::cout << std::left << std::setw(25) << fr.id
            << std::right << std::setw(6) << fr.totalBins
            << std::setw(7) << fr.cleanBins
            << std::fixed << std::setprecision(4)
            << std::setw(10) << s.mean
            << std::setw(10) << s.median
            << std::setw(10) << s.q1
            << std::setw(10) << s.q3 << "\n";
    }

    std::cout << std::string(78, '-') << "\n";
    std::cout << "\nPer-file CSVs:  " << OUTPUT_DIR << "<subject>_wave_comparison.csv\n";
    std::cout << "Summary CSV:    " << summaryPath << "\n";
    std::cout << "SSD Histograms: " << OUTPUT_DIR << "allfiles_ssd_histograms.svg\n";
    std::cout << "Scatter plots:  " << OUTPUT_DIR << "<subject>_*_comparison.svg\n";

    return 0;
}