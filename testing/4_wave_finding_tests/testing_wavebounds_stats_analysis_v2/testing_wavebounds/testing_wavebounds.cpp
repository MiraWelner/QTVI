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
    std::vector<size_t> ecgRIndex;
    std::vector<size_t> ppgMaxAmps;
    std::vector<size_t> ppgMinAmps;
    std::vector<double> ppgSignal;
    std::vector<double> ecgSignal;
    std::vector<std::vector<double>> pairs;
    std::vector<std::pair<uint64_t, uint64_t>> ppg_bin_indexs;
    std::vector<std::pair<uint64_t, uint64_t>> ecg_bin_indexs;
    bool bad_segment = false;
};

struct MatBinData {
    std::vector<size_t> ecgRIndex;
    std::vector<size_t> ppgMaxAmps;
    std::vector<size_t> ppgMinAmps;
    size_t pairCount = 0;
    double ecgSamplingRate = 0.0;
    double ppgSamplingRate = 0.0;
};

struct PerBinResult {
    int bin;
    size_t rPeaksBin;
    size_t rPeaksMat;
    size_t ppgMinBin;
    size_t ppgMinMat;
    double ecg_ssd_seconds;
    double ppg_ssd_seconds;
};

struct FileResult {
    std::string id;
    int totalBins;
    std::vector<PerBinResult> bins;
};

// ============================================================================
// Statistics
// ============================================================================

struct Stats {
    double mean = 0.0;
    double median = 0.0;
    double std_dev = 0.0;
    double sum = 0.0;
    size_t count = 0;
};

Stats computeStats(const std::vector<double>& vals) {
    Stats s;
    if (vals.empty()) return s;
    std::vector<double> sorted = vals;
    std::sort(sorted.begin(), sorted.end());
    s.count = sorted.size();
    for (double v : sorted) s.sum += v;
    s.mean = s.sum / s.count;
    double sumSq = 0.0;
    for (double v : sorted) sumSq += (v - s.mean) * (v - s.mean);
    s.std_dev = s.count > 1 ? std::sqrt(sumSq / (s.count - 1)) : 0.0;
    s.median = (s.count % 2 == 0)
        ? (sorted[s.count / 2 - 1] + sorted[s.count / 2]) / 2.0
        : sorted[s.count / 2];
    return s;
}

// ============================================================================
// Time-based SSD
// ============================================================================

std::vector<double> indicesToSeconds(const std::vector<size_t>& indices, double sr) {
    std::vector<double> times(indices.size());
    for (size_t i = 0; i < indices.size(); ++i)
        times[i] = static_cast<double>(indices[i]) / sr;
    return times;
}

double computeSSD(const std::vector<size_t>& binIdx, double binSR,
    const std::vector<size_t>& matIdx, double matSR,
    double toleranceSec = 0.150)
{
    auto binTimes = indicesToSeconds(binIdx, binSR);
    auto matTimes = indicesToSeconds(matIdx, matSR);

    if (binTimes.empty() || matTimes.empty()) return 0.0;

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
        }
    }

    return sumSqDiff;
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

        readIdx(bin.ecgRIndex);
        skipIdx(); skipIdx();
        skipIdx(); skipIdx(); skipIdx();
        skipIdx(); skipIdx(); skipIdx();

        readIdx(bin.ppgMaxAmps);
        readIdx(bin.ppgMinAmps);

        readSignal(bin.ppgSignal);
        readSignal(bin.ecgSignal);
        skipSignal(); skipSignal();

        skipSignal(); skipSignal();
        skipSignal(); skipSignal();
        skipSignal(); skipSignal();

        uint8_t flags[9] = {};
        file.read(reinterpret_cast<char*>(flags), 9);

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

        readPairVec(bin.ppg_bin_indexs);
        readPairVec(bin.ecg_bin_indexs);

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

    for (int i = 0; i < fr.totalBins; ++i) {
        const auto& b = binData[i];
        const auto& m = matData[i];

        PerBinResult pbr;
        pbr.bin = i;
        pbr.rPeaksBin = b.ecgRIndex.size();
        pbr.rPeaksMat = m.ecgRIndex.size();
        pbr.ppgMinBin = b.ppgMinAmps.size();
        pbr.ppgMinMat = m.ppgMinAmps.size();

        double matEcgSR = m.ecgSamplingRate > 0 ? m.ecgSamplingRate : MAT_SR;
        double matPpgSR = m.ppgSamplingRate > 0 ? m.ppgSamplingRate : MAT_SR;

        pbr.ecg_ssd_seconds = computeSSD(b.ecgRIndex, BIN_SR, m.ecgRIndex, matEcgSR);
        pbr.ppg_ssd_seconds = computeSSD(b.ppgMinAmps, BIN_SR, m.ppgMinAmps, matPpgSR);

        fr.bins.push_back(pbr);
    }
    return fr;
}

// ============================================================================
// Location stats for per-bin CSV
// ============================================================================

struct LocStats {
    size_t count = 0;
    double rate = 0.0;
    double mean = 0.0;
    double median = 0.0;
    double std_dev = 0.0;
};

LocStats computeLocStats(const std::vector<size_t>& indices, double sampleRate) {
    LocStats ls;
    if (indices.empty()) return ls;
    ls.count = indices.size();
    std::vector<double> secs(indices.size());
    for (size_t i = 0; i < indices.size(); ++i)
        secs[i] = static_cast<double>(indices[i]) / sampleRate;
    double duration = secs.back();
    ls.rate = duration > 0 ? (ls.count / duration) * 60.0 : 0.0;
    double sum = 0.0;
    for (double v : secs) sum += v;
    ls.mean = sum / ls.count;
    double sumSq = 0.0;
    for (double v : secs) sumSq += (v - ls.mean) * (v - ls.mean);
    ls.std_dev = ls.count > 1 ? std::sqrt(sumSq / (ls.count - 1)) : 0.0;
    std::sort(secs.begin(), secs.end());
    ls.median = (ls.count % 2 == 0)
        ? (secs[ls.count / 2 - 1] + secs[ls.count / 2]) / 2.0
        : secs[ls.count / 2];
    return ls;
}

// ============================================================================
// Per-file CSV
// ============================================================================

void writePerFileCSV(const std::string& path, const FileResult& fr,
    const std::vector<BinData>& binData,
    const std::vector<MatBinData>& matData) {
    std::ofstream out(path);
    if (!out.is_open()) return;

    size_t n = static_cast<size_t>(fr.totalBins);

    // Header
    out << "bin,"
        << "mat_ecg_count,cpp_ecg_count,mat_ecg_rate,cpp_ecg_rate,mat_ecg_mean,cpp_ecg_mean,mat_ecg_median,cpp_ecg_median,mat_ecg_std,cpp_ecg_std,"
        << "mat_ppg_count,cpp_ppg_count,mat_ppg_rate,cpp_ppg_rate,mat_ppg_mean,cpp_ppg_mean,mat_ppg_median,cpp_ppg_median,mat_ppg_std,cpp_ppg_std,"
        << "ecg_ssd_seconds,ppg_ssd_seconds\n";

    out << std::fixed << std::setprecision(6);

    // Collect all values for summary rows
    std::vector<double> matEcgCounts, cppEcgCounts, matEcgRates, cppEcgRates;
    std::vector<double> matEcgMeans, cppEcgMeans, matEcgMedians, cppEcgMedians;
    std::vector<double> matEcgStds, cppEcgStds;
    std::vector<double> matPpgCounts, cppPpgCounts, matPpgRates, cppPpgRates;
    std::vector<double> matPpgMeans, cppPpgMeans, matPpgMedians, cppPpgMedians;
    std::vector<double> matPpgStds, cppPpgStds;
    std::vector<double> ecgSsds, ppgSsds;

    struct BinRow {
        LocStats matEcg, cppEcg, matPpg, cppPpg;
        double ecgSsd, ppgSsd;
    };
    std::vector<BinRow> rows(n);

    for (size_t i = 0; i < n; ++i) {
        const auto& b = binData[i];
        const auto& m = matData[i];
        const auto& pb = fr.bins[i];

        rows[i].matEcg = computeLocStats(m.ecgRIndex, m.ecgSamplingRate > 0 ? m.ecgSamplingRate : MAT_SR);
        rows[i].cppEcg = computeLocStats(b.ecgRIndex, BIN_SR);
        rows[i].matPpg = computeLocStats(m.ppgMinAmps, m.ppgSamplingRate > 0 ? m.ppgSamplingRate : MAT_SR);
        rows[i].cppPpg = computeLocStats(b.ppgMinAmps, BIN_SR);
        rows[i].ecgSsd = pb.ecg_ssd_seconds;
        rows[i].ppgSsd = pb.ppg_ssd_seconds;

        matEcgCounts.push_back(rows[i].matEcg.count);
        cppEcgCounts.push_back(rows[i].cppEcg.count);
        matEcgRates.push_back(rows[i].matEcg.rate);
        cppEcgRates.push_back(rows[i].cppEcg.rate);
        matEcgMeans.push_back(rows[i].matEcg.mean);
        cppEcgMeans.push_back(rows[i].cppEcg.mean);
        matEcgMedians.push_back(rows[i].matEcg.median);
        cppEcgMedians.push_back(rows[i].cppEcg.median);
        matEcgStds.push_back(rows[i].matEcg.std_dev);
        cppEcgStds.push_back(rows[i].cppEcg.std_dev);
        matPpgCounts.push_back(rows[i].matPpg.count);
        cppPpgCounts.push_back(rows[i].cppPpg.count);
        matPpgRates.push_back(rows[i].matPpg.rate);
        cppPpgRates.push_back(rows[i].cppPpg.rate);
        matPpgMeans.push_back(rows[i].matPpg.mean);
        cppPpgMeans.push_back(rows[i].cppPpg.mean);
        matPpgMedians.push_back(rows[i].matPpg.median);
        cppPpgMedians.push_back(rows[i].cppPpg.median);
        matPpgStds.push_back(rows[i].matPpg.std_dev);
        cppPpgStds.push_back(rows[i].cppPpg.std_dev);
        ecgSsds.push_back(rows[i].ecgSsd);
        ppgSsds.push_back(rows[i].ppgSsd);
    }

    // MEAN row
    auto ms = [](const std::vector<double>& v) { return computeStats(v).mean; };
    out << "MEAN,"
        << ms(matEcgCounts) << "," << ms(cppEcgCounts) << ","
        << ms(matEcgRates) << "," << ms(cppEcgRates) << ","
        << ms(matEcgMeans) << "," << ms(cppEcgMeans) << ","
        << ms(matEcgMedians) << "," << ms(cppEcgMedians) << ","
        << ms(matEcgStds) << "," << ms(cppEcgStds) << ","
        << ms(matPpgCounts) << "," << ms(cppPpgCounts) << ","
        << ms(matPpgRates) << "," << ms(cppPpgRates) << ","
        << ms(matPpgMeans) << "," << ms(cppPpgMeans) << ","
        << ms(matPpgMedians) << "," << ms(cppPpgMedians) << ","
        << ms(matPpgStds) << "," << ms(cppPpgStds) << ","
        << ms(ecgSsds) << "," << ms(ppgSsds) << "\n";

    // SUM row
    auto ss = [](const std::vector<double>& v) { return computeStats(v).sum; };
    out << "SUM,"
        << ss(matEcgCounts) << "," << ss(cppEcgCounts) << ","
        << ss(matEcgRates) << "," << ss(cppEcgRates) << ","
        << ss(matEcgMeans) << "," << ss(cppEcgMeans) << ","
        << ss(matEcgMedians) << "," << ss(cppEcgMedians) << ","
        << ss(matEcgStds) << "," << ss(cppEcgStds) << ","
        << ss(matPpgCounts) << "," << ss(cppPpgCounts) << ","
        << ss(matPpgRates) << "," << ss(cppPpgRates) << ","
        << ss(matPpgMeans) << "," << ss(cppPpgMeans) << ","
        << ss(matPpgMedians) << "," << ss(cppPpgMedians) << ","
        << ss(matPpgStds) << "," << ss(cppPpgStds) << ","
        << ss(ecgSsds) << "," << ss(ppgSsds) << "\n";

    // Per-bin rows
    for (size_t i = 0; i < n; ++i) {
        const auto& r = rows[i];
        out << i << ","
            << r.matEcg.count << "," << r.cppEcg.count << ","
            << r.matEcg.rate << "," << r.cppEcg.rate << ","
            << r.matEcg.mean << "," << r.cppEcg.mean << ","
            << r.matEcg.median << "," << r.cppEcg.median << ","
            << r.matEcg.std_dev << "," << r.cppEcg.std_dev << ","
            << r.matPpg.count << "," << r.cppPpg.count << ","
            << r.matPpg.rate << "," << r.cppPpg.rate << ","
            << r.matPpg.mean << "," << r.cppPpg.mean << ","
            << r.matPpg.median << "," << r.cppPpg.median << ","
            << r.matPpg.std_dev << "," << r.cppPpg.std_dev << ","
            << r.ecgSsd << "," << r.ppgSsd << "\n";
    }
}

// ============================================================================
// Summary CSV
// Rows: Total Subjects, Total Beats (N), Mean RR Interval, Mean Beats/Person, Median Beats/Person
// Columns: mat_ecg, cpp_ecg, mat_ppg, cpp_ppg
// ============================================================================

void writeSummaryCSV(const std::string& path, const std::vector<FileResult>& results,
    const std::vector<std::vector<BinData>>& allBinDataPerFile,
    const std::vector<std::vector<MatBinData>>& allMatDataPerFile) {
    std::ofstream out(path);
    if (!out.is_open()) return;

    out << std::fixed << std::setprecision(6);

    // Collect per-subject totals
    size_t nSubjects = results.size();

    size_t totalMatEcgBeats = 0, totalCppEcgBeats = 0;
    size_t totalMatPpgBeats = 0, totalCppPpgBeats = 0;

    std::vector<double> perSubjMatEcgBeats, perSubjCppEcgBeats;
    std::vector<double> perSubjMatPpgBeats, perSubjCppPpgBeats;

    // For mean RR interval: collect all RR intervals across all subjects
    std::vector<double> allMatEcgRR, allCppEcgRR;
    std::vector<double> allMatPpgRR, allCppPpgRR;

    for (size_t f = 0; f < results.size(); ++f) {
        const auto& binData = allBinDataPerFile[f];
        const auto& matData = allMatDataPerFile[f];
        size_t n = static_cast<size_t>(results[f].totalBins);

        size_t subjMatEcg = 0, subjCppEcg = 0;
        size_t subjMatPpg = 0, subjCppPpg = 0;

        for (size_t i = 0; i < n; ++i) {
            const auto& b = binData[i];
            const auto& m = matData[i];

            subjMatEcg += m.ecgRIndex.size();
            subjCppEcg += b.ecgRIndex.size();
            subjMatPpg += m.ppgMinAmps.size();
            subjCppPpg += b.ppgMinAmps.size();

            double matEcgSR = m.ecgSamplingRate > 0 ? m.ecgSamplingRate : MAT_SR;
            double matPpgSR = m.ppgSamplingRate > 0 ? m.ppgSamplingRate : MAT_SR;

            // RR intervals for mat ECG
            for (size_t k = 1; k < m.ecgRIndex.size(); ++k)
                allMatEcgRR.push_back(static_cast<double>(m.ecgRIndex[k] - m.ecgRIndex[k - 1]) / matEcgSR);

            // RR intervals for cpp ECG
            for (size_t k = 1; k < b.ecgRIndex.size(); ++k)
                allCppEcgRR.push_back(static_cast<double>(b.ecgRIndex[k] - b.ecgRIndex[k - 1]) / BIN_SR);

            // PP intervals for mat PPG
            for (size_t k = 1; k < m.ppgMinAmps.size(); ++k)
                allMatPpgRR.push_back(static_cast<double>(m.ppgMinAmps[k] - m.ppgMinAmps[k - 1]) / matPpgSR);

            // PP intervals for cpp PPG
            for (size_t k = 1; k < b.ppgMinAmps.size(); ++k)
                allCppPpgRR.push_back(static_cast<double>(b.ppgMinAmps[k] - b.ppgMinAmps[k - 1]) / BIN_SR);
        }

        totalMatEcgBeats += subjMatEcg;
        totalCppEcgBeats += subjCppEcg;
        totalMatPpgBeats += subjMatPpg;
        totalCppPpgBeats += subjCppPpg;

        perSubjMatEcgBeats.push_back(static_cast<double>(subjMatEcg));
        perSubjCppEcgBeats.push_back(static_cast<double>(subjCppEcg));
        perSubjMatPpgBeats.push_back(static_cast<double>(subjMatPpg));
        perSubjCppPpgBeats.push_back(static_cast<double>(subjCppPpg));
    }

    auto matEcgBeatStats = computeStats(perSubjMatEcgBeats);
    auto cppEcgBeatStats = computeStats(perSubjCppEcgBeats);
    auto matPpgBeatStats = computeStats(perSubjMatPpgBeats);
    auto cppPpgBeatStats = computeStats(perSubjCppPpgBeats);

    auto matEcgRRStats = computeStats(allMatEcgRR);
    auto cppEcgRRStats = computeStats(allCppEcgRR);
    auto matPpgRRStats = computeStats(allMatPpgRR);
    auto cppPpgRRStats = computeStats(allCppPpgRR);

    // Mean rate (beats per minute) = 60 / mean RR interval
    double matEcgBPM = matEcgRRStats.mean > 0 ? 60.0 / matEcgRRStats.mean : 0.0;
    double cppEcgBPM = cppEcgRRStats.mean > 0 ? 60.0 / cppEcgRRStats.mean : 0.0;
    double matPpgBPM = matPpgRRStats.mean > 0 ? 60.0 / matPpgRRStats.mean : 0.0;
    double cppPpgBPM = cppPpgRRStats.mean > 0 ? 60.0 / cppPpgRRStats.mean : 0.0;

    // Write
    out << ",mat_ecg,cpp_ecg,mat_ppg,cpp_ppg\n";

    out << "Total Subjects,"
        << nSubjects << "," << nSubjects << ","
        << nSubjects << "," << nSubjects << "\n";

    out << "Total Beats (N),"
        << totalMatEcgBeats << "," << totalCppEcgBeats << ","
        << totalMatPpgBeats << "," << totalCppPpgBeats << "\n";

    out << "Mean Rate (beats/min),"
        << matEcgBPM << "," << cppEcgBPM << ","
        << matPpgBPM << "," << cppPpgBPM << "\n";

    out << "Mean Beats Per Person,"
        << matEcgBeatStats.mean << "," << cppEcgBeatStats.mean << ","
        << matPpgBeatStats.mean << "," << cppPpgBeatStats.mean << "\n";
}

// ============================================================================
// Main
// ============================================================================

int main() {
    std::filesystem::create_directories(OUTPUT_DIR);

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
    std::vector<std::vector<BinData>> allBinDataPerFile;
    std::vector<std::vector<MatBinData>> allMatDataPerFile;

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
        allBinDataPerFile.push_back(binData);
        allMatDataPerFile.push_back(matData);

        std::string perFilePath = OUTPUT_DIR + id + "_wave_comparison.csv";
        writePerFileCSV(perFilePath, fr, binData, matData);

        std::vector<double> ecgSsdVec;
        for (const auto& b : fr.bins)
            ecgSsdVec.push_back(b.ecg_ssd_seconds);
        auto s = computeStats(ecgSsdVec);

        std::cout << " " << fr.totalBins << " bins"
            << " ECG SSD[mean=" << std::fixed << std::setprecision(4) << s.mean
            << " med=" << s.median
            << " std=" << s.std_dev << "]\n";
    }

    std::string summaryPath = OUTPUT_DIR + "summary.csv";
    writeSummaryCSV(summaryPath, results, allBinDataPerFile, allMatDataPerFile);

    std::cout << "\n======================================================================\n";
    std::cout << "SUMMARY\n";
    std::cout << "======================================================================\n";
    std::cout << std::left << std::setw(25) << "Subject"
        << std::right << std::setw(6) << "Bins"
        << std::setw(12) << "ECG SSD mn"
        << std::setw(12) << "ECG SSD md"
        << std::setw(12) << "ECG SSD sd" << "\n";
    std::cout << std::string(67, '-') << "\n";

    for (const auto& fr : results) {
        std::vector<double> ecgSsdVec;
        for (const auto& b : fr.bins)
            ecgSsdVec.push_back(b.ecg_ssd_seconds);
        auto s = computeStats(ecgSsdVec);

        std::cout << std::left << std::setw(25) << fr.id
            << std::right << std::setw(6) << fr.totalBins
            << std::fixed << std::setprecision(4)
            << std::setw(12) << s.mean
            << std::setw(12) << s.median
            << std::setw(12) << s.std_dev << "\n";
    }

    std::cout << std::string(67, '-') << "\n";
    std::cout << "\nPer-file CSVs: " << OUTPUT_DIR << "<subject>_wave_comparison.csv\n";
    std::cout << "Summary CSV:   " << summaryPath << "\n";

    return 0;
}