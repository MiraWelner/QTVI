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
 * between two datasets.
 *
 * @param filename Output SVG file path.
 * @param binData The "Reference" or first dataset.
 * @param matData The "Comparison" or second dataset.
 * @param label A descriptive string for the attribute being compared.
 * @param selector A lambda function that extracts the numeric value to plot from a ComparisonData object.
 */
void generateSVGPlot(const std::string& filename,
    const std::vector<ComparisonData>& binData,
    const std::vector<ComparisonData>& matData,
    const std::string& label,
    std::function<double(const ComparisonData&)> selector) {

    size_t n = std::min(binData.size(), matData.size());
    if (n == 0) return;

    // 1. Find the maximum value to set the scale dynamically
    double maxVal = 0;
    for (size_t i = 0; i < n; ++i) {
        maxVal = std::max({ maxVal, selector(binData[i]), selector(matData[i]) });
    }
    // Add 10% padding for visual clarity
    maxVal = (maxVal == 0) ? 10 : maxVal * 1.1;

    // 2. SVG Configuration
    int size = 600;
    int padding = 80; // Increased padding for longer axis labels
    int chartSize = size - (padding * 2);

    std::ofstream svg(filename);
    if (!svg.is_open()) return;

    svg << "<svg width=\"" << size << "\" height=\"" << size << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"#fcfcfc\"/>\n";

    auto mapX = [&](double val) { return padding + (val / maxVal * chartSize); };
    auto mapY = [&](double val) { return (size - padding) - (val / maxVal * chartSize); };

    // 3. Identity Line (Ideal match)
    svg << "<line x1=\"" << mapX(0) << "\" y1=\"" << mapY(0)
        << "\" x2=\"" << mapX(maxVal) << "\" y2=\"" << mapY(maxVal)
        << "\" stroke=\"#999\" stroke-width=\"1\" stroke-dasharray=\"5,5\" />\n";

    // 4. Draw Axes
    svg << "<line x1=\"" << padding << "\" y1=\"" << (size - padding) << "\" x2=\"" << (size - padding) << "\" y2=\"" << (size - padding) << "\" stroke=\"black\" stroke-width=\"2\"/>\n";
    svg << "<line x1=\"" << padding << "\" y1=\"" << padding << "\" x2=\"" << padding << "\" y2=\"" << (size - padding) << "\" stroke=\"black\" stroke-width=\"2\"/>\n";

    // Labels
    svg << "<text x=\"" << size / 2 << "\" y=\"" << size - 20 << "\" text-anchor=\"middle\" font-family=\"sans-serif\">Bin " << label << " (X)</text>\n";
    svg << "<text x=\"" << 25 << "\" y=\"" << size / 2 << "\" text-anchor=\"middle\" transform=\"rotate(-90 25," << size / 2 << ")\" font-family=\"sans-serif\">Mat " << label << " (Y)</text>\n";
    svg << "<text x=\"" << size / 2 << "\" y=\"" << 40 << "\" text-anchor=\"middle\" font-weight=\"bold\" font-family=\"sans-serif\">" << label << " Comparison</text>\n";

    // 5. Plot Data Points
    for (size_t i = 0; i < n; ++i) {
        double x = selector(binData[i]);
        double y = selector(matData[i]);

        // Points are blue if they land exactly on the line (within epsilon), red if they mismatch
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

    // Cast the sum to long long here
    return { mean, stdDev, static_cast<long long>(acc.sum) };
}

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
        // 1. Read Sampling Rates (First)
        file.read(reinterpret_cast<char*>(&results[i].ecgSamplingRate), 8);
        file.read(reinterpret_cast<char*>(&results[i].ppgSamplingRate), 8);

        // 2. Read Indices
        readVecFromBin(file, results[i].ecgRIndex);
        readVecFromBin(file, results[i].ppgMaxAmps);
        readVecFromBin(file, results[i].ppgMinAmps);

        // 3. Read Signals
        readDoubleVecFromBin(file, results[i].ppgSignal);
        readDoubleVecFromBin(file, results[i].ecgSignal);

        // 4. Read Pairs
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

        // Fetch Metadata
        std::vector<double> eRate = getMatField(cell, "ecgSamplingRate");
        results[i].ecgSamplingRate = eRate.empty() ? 0 : eRate[0];

        std::vector<double> pRate = getMatField(cell, "ppgSamplingRate");
        results[i].ppgSamplingRate = pRate.empty() ? 0 : pRate[0];

        // Fetch Indices
        results[i].ecgRIndex = getMatField(cell, "ecgRIndex");
        results[i].ppgMaxAmps = getMatField(cell, "ppgMaxAmps");
        results[i].ppgMinAmps = getMatField(cell, "ppgMinAmps");

        // Fetch Raw Signals
        results[i].ppgSignal = getMatField(cell, "ppgSeg");
        results[i].ecgSignal = getMatField(cell, "ecgSeg"); // Added

        // Fetch Pairs (Existing Logic)
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

void compareIndividual(std::ostream& out, auto bin_files, auto mat_files, const int max_bins) {
    out << "The n for each 60s bin" << std::endl;
    out << "bin #\tR Peak(.bin)\tR Peak(.mat)\tR Peak diff\tPPG Peak(.bin)\tPPG Peak(.mat)\tPPG Peak diff\tPPG Dip(.bin)\tPPG Dip(.mat)\tPPG Dip diff\tPair(.bin)\tPair(.mat)\tPair diff"<< std::endl;
    for (size_t i = 0; i < max_bins; ++i) {

        out << std::setw(3) << std::setfill('0') << i  << std::setfill(' ') << "\t" 

            << bin_files[i].ecgRIndex.size() << "\t\t" << mat_files[i].ecgRIndex.size() << "\t\t" 
            << static_cast<int>(mat_files[i].ecgRIndex.size()) - static_cast<int>(bin_files[i].ecgRIndex.size()) << "\t\t"

            << bin_files[i].ppgMaxAmps.size() << "\t\t" << mat_files[i].ppgMaxAmps.size() << "\t\t"
            << static_cast<int>(mat_files[i].ppgMaxAmps.size()) - static_cast<int>(bin_files[i].ppgMaxAmps.size()) << "\t\t"

            << bin_files[i].ppgMinAmps.size() << "\t\t" << mat_files[i].ppgMinAmps.size() << "\t\t"
            << static_cast<int>(mat_files[i].ppgMinAmps.size()) - static_cast<int>(bin_files[i].ppgMinAmps.size()) << "\t\t"

            << bin_files[i].pairs.size() << "\t\t" << mat_files[i].pairs.size() << "\t\t"
            << static_cast<int>(mat_files[i].pairs.size()) - static_cast<int>(bin_files[i].pairs.size()) << "\t\t"

            << std::endl <<std::endl;
    }
}


int main() {


    std::vector<std::string> IDs = { "3010023_20110817","3010104_20111027","3010112_20110725", "3010139_20110210","3010155_20110511", "3010201_20120320","3010228_20110426", "3010317_20110413", "3010457_20111109", "3010660_20120322", "3010724_20110811", "3010740_20110303", "3010767_20120423"};
    std::vector<int> R_peaks_all_files_bins = {};
    std::vector<int> R_peaks_all_files_mat = {};
    std::vector<int> bins_binfile = {};
    std::vector<int> bins_matfile = {};

    for (int z = 0; z < IDs.size(); z++) {
        std::string ID = IDs[z];
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
        outFile << "Data taken from MESA file " << ID << std::endl;
        outFile << "Bins in .bin file: " << binData.size() << std::endl << "Bins in .mat file: " << matData.size() << std::endl << std::endl;


        //stats among all bins
        auto rPeakProjector = [](const ComparisonData& c) { return c.ecgRIndex.size(); };
        auto ppgPeakProjector = [](const ComparisonData& c) { return c.ppgMaxAmps.size(); };
        auto ppgDipProjector = [](const ComparisonData& c) { return c.ppgMinAmps.size(); };
        auto pairProjector = [](const ComparisonData& c) { return c.pairs.size(); };
        
        // Calculate metrics
        DatasetStats binRStats = calculateStats(binData, rPeakProjector);
        DatasetStats matRStats = calculateStats(matData, rPeakProjector);
        DatasetStats binPPGPeakStats = calculateStats(binData, ppgPeakProjector);
        DatasetStats matPPGPeakStats = calculateStats(matData, ppgPeakProjector);
        DatasetStats binPPGDipStats = calculateStats(binData, ppgDipProjector);
        DatasetStats matPPGDipStats = calculateStats(matData, ppgDipProjector);
        DatasetStats binPairStats = calculateStats(binData, pairProjector);
        DatasetStats matPairStats = calculateStats(matData, pairProjector);


        //get number of bins with certain error
        int bins_r_error_25 = 0;
        int bins_r_error_50 = 0;
        int bins_ppg_peak_error_25 = 0;
        int bins_ppg_peak_error_50 = 0;
        int bins_ppg_dip_error_25 = 0;
        int bins_ppg_dip_error_50 = 0;
        int bins_pair_error_25 = 0;
        int bins_pair_error_50 = 0;

        for (size_t i = 0; i < max_bins; ++i) {
                
            if (std::abs(static_cast<int>(matData[i].ecgRIndex.size()) - static_cast<int>(binData[i].ecgRIndex.size()) > 25)) {
                bins_r_error_25 += 1;
            }
            if (std::abs(static_cast<int>(matData[i].ecgRIndex.size()) - static_cast<int>(binData[i].ecgRIndex.size()) > 50)) {
                bins_r_error_50 += 1;
            }

            if (std::abs(static_cast<int>(matData[i].ppgMaxAmps.size()) - static_cast<int>(binData[i].ppgMaxAmps.size()) > 25)) {
                bins_ppg_peak_error_25 += 1;
            }
            if (std::abs(static_cast<int>(matData[i].ppgMaxAmps.size()) - static_cast<int>(binData[i].ppgMaxAmps.size()) > 50)) {
                bins_ppg_peak_error_50 += 1;
            }
            if (std::abs(static_cast<int>(matData[i].ppgMinAmps.size()) - static_cast<int>(binData[i].ppgMinAmps.size()) > 25)) {
                bins_ppg_dip_error_25 += 1;
            }
            if (std::abs(static_cast<int>(matData[i].ppgMinAmps.size()) - static_cast<int>(binData[i].ppgMinAmps.size()) > 50)) {
                bins_ppg_dip_error_50 += 1;
            }
            if (std::abs(static_cast<int>(matData[i].pairs.size()) - static_cast<int>(binData[i].pairs.size()) > 25)) {
                bins_pair_error_25 += 1;
            }
            if (std::abs(static_cast<int>(matData[i].pairs.size()) - static_cast<int>(binData[i].pairs.size()) > 50)) {
                bins_pair_error_50 += 1;
            }
        }



        outFile << "\t\tTotal in .bin\tTotal in .mat\tTotal diff\tSTD in .bin\tSTD in .mat\tSTD diff\t%  error\tbins >25 errors\tbins >50 error" << std::endl;

        outFile << "R Peaks" << std::fixed << std::setprecision(3) << "\t\t" << binRStats.total << "\t\t" << matRStats.total << "\t\t" << matRStats.total - binRStats.total
            << "\t\t" << binRStats.stdDev << "\t\t" << matRStats.stdDev << "\t\t" << matRStats.stdDev - binRStats.stdDev << "\t\t" 
           << (matRStats.total-binRStats.total)/ static_cast<double>(matRStats.total) << "\t\t" << bins_r_error_25 << "\t\t" << bins_r_error_50 << std::endl << std::endl;

        outFile << "PPG Peaks" << std::fixed << std::setprecision(3) << "\t" << binPPGPeakStats.total << "\t\t" << matPPGPeakStats.total << "\t\t" << matPPGPeakStats.total - binPPGPeakStats.total
            << "\t\t" << binPPGPeakStats.stdDev << "\t\t" << matPPGPeakStats.stdDev << "\t\t" << matPPGPeakStats.stdDev - binPPGPeakStats.stdDev << "\t\t"
            << (matPPGPeakStats.total - binPPGPeakStats.total) / static_cast<double>(matPPGPeakStats.total) << "\t\t" << bins_ppg_peak_error_25 << "\t\t" << bins_ppg_peak_error_50 << std::endl << std::endl;

        outFile << "PPG Dips" << std::fixed << std::setprecision(3) << "\t" << binPPGDipStats.total << "\t\t" << matPPGDipStats.total << "\t\t" << matPPGDipStats.total - binPPGDipStats.total
            << "\t\t" << binPPGDipStats.stdDev << "\t\t" << matPPGDipStats.stdDev << "\t\t" << matPPGDipStats.stdDev - binPPGDipStats.stdDev << "\t\t"
            << (matPPGDipStats.total - binPPGDipStats.total) / static_cast<double>(matPPGDipStats.total) <<
            "\t\t" << bins_ppg_dip_error_25 << "\t\t" << bins_ppg_dip_error_50 << std::endl << std::endl;

        outFile << "Pairs" << std::fixed << std::setprecision(3) << "\t\t" << binPairStats.total << "\t\t" << matPairStats.total << "\t\t" << matPairStats.total - binPairStats.total
            << "\t\t" << binPairStats.stdDev << "\t\t" << matPairStats.stdDev << "\t\t" << matPairStats.stdDev - binPairStats.stdDev << "\t\t"
             << (matPairStats.total - binPairStats.total) / static_cast<double>(matPairStats.total) << 
            "\t\t" << bins_pair_error_25 << "\t\t" << bins_pair_error_50 << std::endl << std::endl;


        compareIndividual(outFile, binData, matData, max_bins);


        outFile.close();

        R_peaks_all_files_bins.push_back(binRStats.total);
        R_peaks_all_files_mat.push_back(matRStats.total);
        bins_binfile.push_back(binData.size());
        bins_matfile.push_back(matData.size());

        // Generate the plot
        generateSVGPlot(ID+"_ecg_comparison.svg", binData, matData, "R-Peak Count",
            [](const ComparisonData& d) { return (double)d.ecgRIndex.size(); });

        generateSVGPlot(ID+"_ppg_max_comparison.svg", binData, matData, "PPG Max Amps Count",
            [](const ComparisonData& d) { return (double)d.ppgMaxAmps.size(); });

        generateSVGPlot(ID+"_ppg_min_comparison.svg", binData, matData, "PPG Min Amps Count",
            [](const ComparisonData& d) { return (double)d.ppgMinAmps.size(); });

        generateSVGPlot(ID+"_pairs_comparison.svg", binData, matData, "Detected Pairs",
            [](const ComparisonData& d) { return (double)d.pairs.size(); });

        std::cout << "Comparison complete. Results written to: " << outPath << std::endl;
    }

    std::ofstream total_outfile("allfiles.txt");

    if (!total_outfile) {
        std::cerr << "Error opening file.\n";
        return 1;
    }

    // Compute sums
    int sum_R_bins = std::accumulate(R_peaks_all_files_bins.begin(),
        R_peaks_all_files_bins.end(), 0);

    int sum_R_mat = std::accumulate(R_peaks_all_files_mat.begin(),
        R_peaks_all_files_mat.end(), 0);

    int sum_bins_bin = std::accumulate(bins_binfile.begin(),
        bins_binfile.end(), 0);

    int sum_bins_mat = std::accumulate(bins_matfile.begin(),
        bins_matfile.end(), 0);

    // Write header
    total_outfile << std::fixed << std::setprecision(3) <<  "\t\t\t\t" <<  "Total R Peaks (.bin)" << "\t" << "Total R Peaks (.mat)" << "\t" << "R Peaks Diff" << "\t\t" << "% error"  << "\t\t\t" << "Total Bins (.bin)" << "\t" << "Total Bins (.mat)" << "\n";

    // Write sums row
    total_outfile << "All Tested Files Sum" << "\t\t" << sum_R_bins << "\t\t\t" << sum_R_mat << "\t\t\t" << sum_R_mat - sum_R_bins << "\t\t\t" << (sum_R_mat - sum_R_bins)/static_cast<double>(sum_R_mat) <<  "\t\t\t" << sum_bins_bin << "\t\t\t" << sum_bins_mat << "\n\n";

    // Write per-file values
    size_t n = IDs.size();
    for (size_t i = 0; i < n; ++i) {
        total_outfile << "Mesa ID=" << IDs[i] << "\t" << R_peaks_all_files_bins[i] << "\t\t\t" << R_peaks_all_files_mat[i] << "\t\t\t" 
            << R_peaks_all_files_mat[i] - R_peaks_all_files_bins[i]<< "\t\t\t"<< (R_peaks_all_files_mat[i] - R_peaks_all_files_bins[i]) / static_cast<double>(R_peaks_all_files_mat[i]) << "\t\t\t" << bins_binfile[i] << "\t\t\t" << bins_matfile[i] << "\n";
    }

    total_outfile.close();
    return 0;
}
