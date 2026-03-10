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

struct ComparisonData {
    double ecgSamplingRate;
    double ppgSamplingRate;
    std::vector<double> ecgRIndex;
    std::vector<double> ppgMaxAmps;
    std::vector<double> ppgMinAmps;
    std::vector<double> ppgSignal;
    std::vector<double> ecgSignal;
    std::vector<std::pair<double, double>> pairs;
};

struct DatasetStats {
    double mean = 0.0;
    double stdDev = 0.0;
    long long total = 0;
};

/**
 * @brief Generates an SVG 45-degree identity plot comparing specific attributes
 */
void generateSVGPlot(const std::string& filename,
    const std::vector<ComparisonData>& binData,
    const std::vector<ComparisonData>& matData,
    const std::string& label,
    std::function<double(const ComparisonData&)> selector) {

    size_t n = std::min(binData.size(), matData.size());
    if (n == 0) return;

    double maxVal = 0;
    for (size_t i = 0; i < n; ++i) {
        maxVal = std::max({ maxVal, selector(binData[i]), selector(matData[i]) });
    }
    maxVal = (maxVal == 0) ? 10 : maxVal * 1.1;

    int size = 600;
    int padding = 80;
    int chartSize = size - (padding * 2);

    std::ofstream svg(filename);
    if (!svg.is_open()) return;

    svg << "<svg width=\"" << size << "\" height=\"" << size << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"#fcfcfc\"/>\n";

    auto mapX = [&](double val) { return padding + (val / maxVal * chartSize); };
    auto mapY = [&](double val) { return (size - padding) - (val / maxVal * chartSize); };

    svg << "<line x1=\"" << mapX(0) << "\" y1=\"" << mapY(0)
        << "\" x2=\"" << mapX(maxVal) << "\" y2=\"" << mapY(maxVal)
        << "\" stroke=\"#999\" stroke-width=\"1\" stroke-dasharray=\"5,5\" />\n";

    svg << "<line x1=\"" << padding << "\" y1=\"" << (size - padding) << "\" x2=\"" << (size - padding) << "\" y2=\"" << (size - padding) << "\" stroke=\"black\" stroke-width=\"2\"/>\n";
    svg << "<line x1=\"" << padding << "\" y1=\"" << padding << "\" x2=\"" << padding << "\" y2=\"" << (size - padding) << "\" stroke=\"black\" stroke-width=\"2\"/>\n";

    svg << "<text x=\"" << size / 2 << "\" y=\"" << size - 20 << "\" text-anchor=\"middle\" font-family=\"sans-serif\">Bin " << label << " (X)</text>\n";
    svg << "<text x=\"" << 25 << "\" y=\"" << size / 2 << "\" text-anchor=\"middle\" transform=\"rotate(-90 25," << size / 2 << ")\" font-family=\"sans-serif\">Mat " << label << " (Y)</text>\n";
    svg << "<text x=\"" << size / 2 << "\" y=\"" << 40 << "\" text-anchor=\"middle\" font-weight=\"bold\" font-family=\"sans-serif\">" << label << " Comparison</text>\n";

    for (size_t i = 0; i < n; ++i) {
        double x = selector(binData[i]);
        double y = selector(matData[i]);
        std::string color = (std::abs(x - y) < 0.0001) ? "#3498db" : "#e74c3c";

        svg << "<circle cx=\"" << mapX(x) << "\" cy=\"" << mapY(y)
            << "\" r=\"3.5\" fill=\"" << color << "\" fill-opacity=\"0.7\">\n"
            << "  <title>Bin: " << x << ", Mat: " << y << "</title>\n"
            << "</circle>\n";
    }

    svg << "</svg>";
    svg.close();
}

template <typename Projection>
DatasetStats calculateStats(const std::vector<ComparisonData>& data, Projection projector) {
    if (data.empty()) return {};

    struct Accumulator {
        double sum = 0.0;
        double sumSq = 0.0;
    };

    Accumulator acc = std::accumulate(data.begin(), data.end(), Accumulator{},
        [&](Accumulator a, const ComparisonData& item) {
            double val = static_cast<double>(projector(item));
            return Accumulator{ a.sum + val, a.sumSq + (val * val) };
        });

    double n = static_cast<double>(data.size());
    double mean = acc.sum / n;
    double variance = std::max(0.0, (acc.sumSq / n) - (mean * mean));
    double stdDev = std::sqrt(variance);

    return { mean, stdDev, static_cast<long long>(acc.sum) };
}

/**
 * @brief Reads a vector of doubles that were stored as int64_t in the binary file.
 * Typically used for indices.
 */
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

/**
 * @brief Reads a vector of doubles directly from binary.
 * Typically used for amplitudes and signals.
 */
void readDoubleVecFromBin(std::ifstream& file, std::vector<double>& vec) {
    uint64_t count = 0;
    if (!file.read(reinterpret_cast<char*>(&count), 8)) return;
    vec.resize(count);
    if (count > 0) {
        file.read(reinterpret_cast<char*>(vec.data()), count * 8);
    }
}

std::vector<ComparisonData> readAllFromBin(const std::string& path) {
    std::vector<ComparisonData> results;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return results;

    uint64_t numBins = 0;
    if (!file.read(reinterpret_cast<char*>(&numBins), 8)) return results;
    results.resize(numBins);

    for (uint64_t i = 0; i < numBins; ++i) {
        // 1. Read Sampling Rates
        file.read(reinterpret_cast<char*>(&results[i].ecgSamplingRate), 8);
        file.read(reinterpret_cast<char*>(&results[i].ppgSamplingRate), 8);

        // 2. Read Indices and Amplitudes
        // Check if your binary writer saves PPG Peaks/Dips as doubles (amplitudes) or int64 (indices)
        // If they are amplitudes, use readDoubleVecFromBin. 
        // If they are indices, use readInt64VecAsDouble.
        readInt64VecAsDouble(file, results[i].ecgRIndex);
        readDoubleVecFromBin(file, results[i].ppgMaxAmps); // Fixed: Use double reader for Amps
        readDoubleVecFromBin(file, results[i].ppgMinAmps); // Fixed: Use double reader for Amps

        // 3. Read Raw Signals
        readDoubleVecFromBin(file, results[i].ppgSignal);
        readDoubleVecFromBin(file, results[i].ecgSignal);

        // 4. Read Pairs
        uint64_t numPairs = 0;
        if (file.read(reinterpret_cast<char*>(&numPairs), 8)) {
            for (uint64_t j = 0; j < numPairs; ++j) {
                int64_t ppgIdx, ecgIdx;
                file.read(reinterpret_cast<char*>(&ppgIdx), 8);
                file.read(reinterpret_cast<char*>(&ecgIdx), 8);
                results[i].pairs.push_back({ static_cast<double>(ppgIdx), static_cast<double>(ecgIdx) });
            }
        }
    }
    return results;
}

std::vector<double> getMatField(matvar_t* cell, const std::string& fieldName) {
    std::vector<double> vec;
    // Try primary name, then common alternatives if NULL
    matvar_t* field = Mat_VarGetStructFieldByName(cell, fieldName.c_str(), 0);
    if (!field) {
        if (fieldName == "ppgMaxAmps") field = Mat_VarGetStructFieldByName(cell, "ppgMaxIdx", 0);
        if (fieldName == "ppgMaxAmps" && !field) field = Mat_VarGetStructFieldByName(cell, "maxAmps", 0);
        if (fieldName == "ppgMinAmps") field = Mat_VarGetStructFieldByName(cell, "ppgMinIdx", 0);
    }

    if (field && field->data) {
        size_t n = 1;
        for (int i = 0; i < field->rank; ++i) n *= field->dims[i];

        // Handle different data types safely
        if (field->class_type == MAT_C_DOUBLE) {
            double* data = (double*)field->data;
            for (size_t j = 0; j < n; ++j) vec.push_back(data[j]);
        }
        else if (field->class_type == MAT_C_UINT64) {
            uint64_t* data = (uint64_t*)field->data;
            for (size_t j = 0; j < n; ++j) vec.push_back(static_cast<double>(data[j]));
        }
        else if (field->class_type == MAT_C_INT64) {
            int64_t* data = (int64_t*)field->data;
            for (size_t j = 0; j < n; ++j) vec.push_back(static_cast<double>(data[j]));
        }
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

        std::vector<double> eRate = getMatField(cell, "ecgSamplingRate");
        results[i].ecgSamplingRate = eRate.empty() ? 0 : eRate[0];

        std::vector<double> pRate = getMatField(cell, "ppgSamplingRate");
        results[i].ppgSamplingRate = pRate.empty() ? 0 : pRate[0];

        results[i].ecgRIndex = getMatField(cell, "ecgRIndex");
        results[i].ppgMaxAmps = getMatField(cell, "ppgMaxAmps");
        results[i].ppgMinAmps = getMatField(cell, "ppgMinAmps");
        results[i].ppgSignal = getMatField(cell, "ppgSignal");
        results[i].ecgSignal = getMatField(cell, "ecgSignal");

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

void compareIndividual(std::ostream& out, const std::vector<ComparisonData>& bin_files, const std::vector<ComparisonData>& mat_files, const size_t max_bins) {
    out << "The n for each 60s bin" << std::endl;
    out << "bin #\tR Peak(.bin)\tR Peak(.mat)\tR Peak diff\tPPG Peak(.bin)\tPPG Peak(.mat)\tPPG Peak diff\tPPG Dip(.bin)\tPPG Dip(.mat)\tPPG Dip diff\tPair(.bin)\tPair(.mat)\tPair diff" << std::endl;
    for (size_t i = 0; i < max_bins; ++i) {
        out << std::setw(3) << std::setfill('0') << i << std::setfill(' ') << "\t"
            << bin_files[i].ecgRIndex.size() << "\t\t" << mat_files[i].ecgRIndex.size() << "\t\t"
            << static_cast<int>(mat_files[i].ecgRIndex.size()) - static_cast<int>(bin_files[i].ecgRIndex.size()) << "\t\t"
            << bin_files[i].ppgMaxAmps.size() << "\t\t" << mat_files[i].ppgMaxAmps.size() << "\t\t"
            << static_cast<int>(mat_files[i].ppgMaxAmps.size()) - static_cast<int>(bin_files[i].ppgMaxAmps.size()) << "\t\t"
            << bin_files[i].ppgMinAmps.size() << "\t\t" << mat_files[i].ppgMinAmps.size() << "\t\t"
            << static_cast<int>(mat_files[i].ppgMinAmps.size()) - static_cast<int>(bin_files[i].ppgMinAmps.size()) << "\t\t"
            << bin_files[i].pairs.size() << "\t\t" << mat_files[i].pairs.size() << "\t\t"
            << static_cast<int>(mat_files[i].pairs.size()) - static_cast<int>(bin_files[i].pairs.size()) << "\t\t"
            << std::endl;
    }
}

int main() {
    std::vector<std::string> IDs = {"3010023_20110817", "3010104_20111027", "3010112_20110725", "3010139_20110210", "3010201_20120320","3010228_20110426", "3010317_20110413", "3010457_20111109", "3010660_20120322", "3010724_20110811", "3010740_20110303" };
    std::vector<long long> R_peaks_all_files_bins;
    std::vector<long long> R_peaks_all_files_mat;
    std::vector<int> bins_binfile;
    std::vector<int> bins_matfile;

    for (const std::string& ID : IDs) {
        std::string binPath = "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\4_wave_bound_files\\cpp\\mesa_files\\" + ID + "_annealed_wave_data.bin";
        std::string matPath = "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\4_wave_bound_files\\matlab\\" + ID + "_wave_data.mat";
        std::string outPath = ID + ".txt";

        std::ofstream outFile(outPath);
        if (!outFile.is_open()) {
            std::cerr << "Failed to open output file: " << outPath << std::endl;
            continue;
        }

        auto binData = readAllFromBin(binPath);
        auto matData = readAllFromMat(matPath);

        size_t max_bins = std::min(binData.size(), matData.size());
        outFile << "Data taken from MESA file " << ID << std::endl;
        outFile << "Bins in .bin file: " << binData.size() << std::endl << "Bins in .mat file: " << matData.size() << std::endl << std::endl;

        auto rPeakProjector = [](const ComparisonData& c) { return (double)c.ecgRIndex.size(); };
        auto ppgPeakProjector = [](const ComparisonData& c) { return (double)c.ppgMaxAmps.size(); };
        auto ppgDipProjector = [](const ComparisonData& c) { return (double)c.ppgMinAmps.size(); };
        auto pairProjector = [](const ComparisonData& c) { return (double)c.pairs.size(); };

        DatasetStats binRStats = calculateStats(binData, rPeakProjector);
        DatasetStats matRStats = calculateStats(matData, rPeakProjector);
        DatasetStats binPPGPeakStats = calculateStats(binData, ppgPeakProjector);
        DatasetStats matPPGPeakStats = calculateStats(matData, ppgPeakProjector);
        DatasetStats binPPGDipStats = calculateStats(binData, ppgDipProjector);
        DatasetStats matPPGDipStats = calculateStats(matData, ppgDipProjector);
        DatasetStats binPairStats = calculateStats(binData, pairProjector);
        DatasetStats matPairStats = calculateStats(matData, pairProjector);

        int bins_r_error_5 = 0, bins_r_error_10 = 0;
        int bins_ppg_peak_error_5 = 0, bins_ppg_peak_error_10 = 0;
        int bins_ppg_dip_error_5 = 0, bins_ppg_dip_error_10 = 0;
        int bins_pair_error_5 = 0, bins_pair_error_10 = 0;

        for (size_t i = 0; i < max_bins; ++i) {
            auto checkErr = [&](size_t b, size_t m, int& e5, int& e10) {
                int diff = std::abs(static_cast<int>(m) - static_cast<int>(b));
                if (diff > 5) e5++;
                if (diff > 50) e10++;
                };
            checkErr(binData[i].ecgRIndex.size(), matData[i].ecgRIndex.size(), bins_r_error_5, bins_r_error_10);
            checkErr(binData[i].ppgMaxAmps.size(), matData[i].ppgMaxAmps.size(), bins_ppg_peak_error_5, bins_ppg_peak_error_10);
            checkErr(binData[i].ppgMinAmps.size(), matData[i].ppgMinAmps.size(), bins_ppg_dip_error_5, bins_ppg_dip_error_10);
            checkErr(binData[i].pairs.size(), matData[i].pairs.size(), bins_pair_error_5, bins_pair_error_10);
        }

        outFile << "\t\tTotal in .bin\tTotal in .mat\tTotal diff\tSTD in .bin\tSTD in .mat\tSTD diff\t% error\t\t>5 errors\t>10 error" << std::endl;

        auto writeRow = [&](std::string label, DatasetStats b, DatasetStats m, int e25, int e50) {
            double errPct = (m.total == 0) ? 0 : (m.total - b.total) / static_cast<double>(m.total);
            outFile << label << std::fixed << std::setprecision(3) << "\t" << b.total << "\t\t" << m.total << "\t\t" << m.total - b.total
                << "\t\t" << b.stdDev << "\t\t" << m.stdDev << "\t\t" << m.stdDev - b.stdDev << "\t\t"
                << errPct << "\t\t" << e25 << "\t\t" << e50 << std::endl << std::endl;
            };

        writeRow("R Peaks\t", binRStats, matRStats, bins_r_error_5, bins_r_error_10);
        writeRow("PPG Peaks", binPPGPeakStats, matPPGPeakStats, bins_ppg_peak_error_5, bins_ppg_peak_error_10);
        writeRow("PPG Dips", binPPGDipStats, matPPGDipStats, bins_ppg_dip_error_5, bins_ppg_dip_error_10);
        writeRow("Pairs\t", binPairStats, matPairStats, bins_pair_error_5, bins_pair_error_10);

        compareIndividual(outFile, binData, matData, max_bins);
        outFile.close();

        R_peaks_all_files_bins.push_back(binRStats.total);
        R_peaks_all_files_mat.push_back(matRStats.total);
        bins_binfile.push_back((int)binData.size());
        bins_matfile.push_back((int)matData.size());

        generateSVGPlot(ID + "_ecg_comparison.svg", binData, matData, "R-Peak Count", [](const ComparisonData& d) { return (double)d.ecgRIndex.size(); });
        generateSVGPlot(ID + "_ppg_max_comparison.svg", binData, matData, "PPG Max Amps Count", [](const ComparisonData& d) { return (double)d.ppgMaxAmps.size(); });
        generateSVGPlot(ID + "_ppg_min_comparison.svg", binData, matData, "PPG Min Amps Count", [](const ComparisonData& d) { return (double)d.ppgMinAmps.size(); });
        generateSVGPlot(ID + "_pairs_comparison.svg", binData, matData, "Detected Pairs", [](const ComparisonData& d) { return (double)d.pairs.size(); });

        std::cout << "Comparison complete for: " << ID << std::endl;
    }

    std::ofstream total_outfile("allfiles.txt");
    if (total_outfile) {
        long long sum_R_bins = std::accumulate(R_peaks_all_files_bins.begin(), R_peaks_all_files_bins.end(), 0LL);
        long long sum_R_mat = std::accumulate(R_peaks_all_files_mat.begin(), R_peaks_all_files_mat.end(), 0LL);
        int sum_bins_bin = std::accumulate(bins_binfile.begin(), bins_binfile.end(), 0);
        int sum_bins_mat = std::accumulate(bins_matfile.begin(), bins_matfile.end(), 0);

        total_outfile << std::fixed << std::setprecision(3) << "\t\t\t\t" << "Total R Peaks (.bin)" << "\t" << "Total R Peaks (.mat)" << "\t" << "R Peaks Diff" << "\t\t" << "% error" << "\t\t\t" << "Total Bins (.bin)" << "\t" << "Total Bins (.mat)" << "\n";
        total_outfile << "All Tested Files Sum" << "\t\t" << sum_R_bins << "\t\t\t" << sum_R_mat << "\t\t\t" << sum_R_mat - sum_R_bins << "\t\t\t" << (sum_R_mat == 0 ? 0 : (sum_R_mat - sum_R_bins) / static_cast<double>(sum_R_mat)) << "\t\t\t" << sum_bins_bin << "\t\t\t" << sum_bins_mat << "\n\n";

        for (size_t i = 0; i < IDs.size(); ++i) {
            double err = (R_peaks_all_files_mat[i] == 0) ? 0 : (R_peaks_all_files_mat[i] - R_peaks_all_files_bins[i]) / static_cast<double>(R_peaks_all_files_mat[i]);
            total_outfile << "Mesa ID=" << IDs[i] << "\t" << R_peaks_all_files_bins[i] << "\t\t\t" << R_peaks_all_files_mat[i] << "\t\t\t" << R_peaks_all_files_mat[i] - R_peaks_all_files_bins[i] << "\t\t\t" << err << "\t\t\t" << bins_binfile[i] << "\t\t\t" << bins_matfile[i] << "\n";
        }
        total_outfile.close();
    }

    return 0;
}

