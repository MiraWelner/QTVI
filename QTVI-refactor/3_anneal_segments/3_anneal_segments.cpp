#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <cstdint>

namespace fs = std::filesystem;

// ============================================================================
// DATA STRUCTURES
// ============================================================================

struct ProjectConfig {
    std::string dataType, binPath, noisePath, annealedPath;
    double ecgSR, ppgSR;
};

struct RawData {
    std::vector<double> ppg, ecg, sleepStages;
    double ppgSR, ecgSR, scoringEpochSec;
};

struct FinalSegment {
    std::vector<double> ppg;
    std::vector<double> ecg;
    std::vector<int> sleep_stages;
};

// ============================================================================
// HELPERS
// ============================================================================

std::string cleanField(std::string s) {
    s.erase(0, s.find_first_not_of(" \t\r\n\""));
    s.erase(s.find_last_not_of(" \t\r\n\"") + 1);
    return s;
}

int timeToSleepIdx(uint64_t ppgIdx, double ppgSR, double epochSec) {
    double timeSec = (double)(ppgIdx - 1) / ppgSR;
    return static_cast<int>(std::floor(timeSec / epochSec));
}

// ============================================================================
// FILE I/O (BINARY & CSV)
// ============================================================================

RawData readStructuredBin(const std::string& path) {
    RawData data;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Cannot open bin: " + path);

    uint64_t nEcg1, nEcg2, nEcg3, nPpg, nSleep;
    file.read((char*)&data.ecgSR, 8);
    file.read((char*)&data.ppgSR, 8);
    file.read((char*)&data.scoringEpochSec, 8);
    file.read((char*)&nEcg1, 8); file.read((char*)&nEcg2, 8); file.read((char*)&nEcg3, 8);
    file.read((char*)&nPpg, 8);
    file.read((char*)&nSleep, 8);

    data.ecg.resize(nEcg1 + nEcg2 + nEcg3);
    file.read((char*)data.ecg.data(), data.ecg.size() * 8);
    data.ppg.resize(nPpg);
    file.read((char*)data.ppg.data(), nPpg * 8);
    data.sleepStages.resize(nSleep);
    file.read((char*)data.sleepStages.data(), nSleep * 8);
    return data;
}

std::vector<std::pair<double, double>> readNoiseCSV(const std::string& path) {
    std::vector<std::pair<double, double>> noise;
    std::ifstream file(path);
    if (!file.is_open()) return noise; // Return empty if no noise file exists
    std::string line; std::getline(file, line); // skip header
    while (std::getline(file, line)) {
        std::stringstream ss(line); std::string col; std::vector<std::string> cols;
        while (std::getline(ss, col, ',')) cols.push_back(cleanField(col));
        if (cols.size() >= 4) noise.push_back({ std::stod(cols[2]), std::stod(cols[3]) });
    }
    return noise;
}

// ============================================================================
// CORE ANNEALING ALGORITHM
// ============================================================================

#include <algorithm>
#include <cmath>
#include <vector>

// Helper to merge adjacent segments (MATLAB's MergeSegments)
std::vector<std::pair<uint64_t, uint64_t>> MergeSegments(std::vector<std::pair<uint64_t, uint64_t>>& segs) {
    if (segs.empty()) return {};
    std::sort(segs.begin(), segs.end());
    std::vector<std::pair<uint64_t, uint64_t>> merged;
    merged.push_back(segs[0]);
    for (size_t i = 1; i < segs.size(); ++i) {
        if (segs[i].first <= merged.back().second) {
            merged.back().second = std::max(merged.back().second, segs[i].second);
        }
        else {
            merged.push_back(segs[i]);
        }
    }
    return merged;
}

std::vector<FinalSegment> AnnealSegments(const RawData& data,
    const std::vector<std::pair<double, double>>& noise,
    double targetLenMins)
{
    const double ppgSR = data.ppgSR;
    const uint64_t bin_samples = static_cast<uint64_t>(ppgSR * 60.0 * targetLenMins);
    const double min_bin_size_samples = (targetLenMins / 2.0) * 60.0 * ppgSR;
    const uint64_t total_samples = data.ppg.size();

    // 1. Ideal Bin Breaks (MATLAB: bin_breaks)
    size_t bin_count = total_samples / bin_samples;
    if (total_samples % bin_samples != 0) bin_count++;

    std::vector<uint64_t> bin_breaks;
    for (size_t i = 1; i <= bin_count; ++i) {
        uint64_t brk = i * bin_samples;
        bin_breaks.push_back(std::min(brk, total_samples));
    }

    // 2. Extract "Good" (Non-Noise) Fragments
    std::vector<std::pair<uint64_t, uint64_t>> good_frags;
    uint64_t cur = 0;
    auto sorted_noise = noise;
    std::sort(sorted_noise.begin(), sorted_noise.end());

    for (const auto& n : sorted_noise) {
        uint64_t n_start = static_cast<uint64_t>(n.first * ppgSR);
        uint64_t n_end = static_cast<uint64_t>(n.second * ppgSR);
        if (n_start > cur) good_frags.push_back({ cur, n_start });
        cur = std::max(cur, n_end);
    }
    if (cur < total_samples) good_frags.push_back({ cur, total_samples });

    // 3. Assign Fragments to Bins based on Midpoint
    // This structure holds the sample ranges assigned to each bin
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> bin_assignments(bin_count);

    for (auto& frag : good_frags) {
        double midpoint = (frag.first + frag.second) / 2.0;
        int target_bin = static_cast<int>(std::floor(midpoint / bin_samples));

        if (target_bin >= 0 && target_bin < (int)bin_count) {
            bin_assignments[target_bin].push_back(frag);
        }
    }

    // 4. Construct Final Segments
    std::vector<FinalSegment> results(bin_count);
    double ratio = data.ecgSR / ppgSR;

    for (size_t i = 0; i < bin_count; ++i) {
        auto merged = MergeSegments(bin_assignments[i]);

        for (auto& seg : merged) {
            // PPG
            for (uint64_t k = seg.first; k < seg.second; ++k) {
                results[i].ppg.push_back(data.ppg[k]);
            }
            // ECG
            uint64_t ecg_s = static_cast<uint64_t>(seg.first * ratio);
            uint64_t ecg_e = static_cast<uint64_t>(seg.second * ratio);
            for (uint64_t k = ecg_s; k < ecg_e; ++k) {
                if (k < data.ecg.size()) results[i].ecg.push_back(data.ecg[k]);
            }
            // Sleep (Expanded)
            for (uint64_t k = seg.first; k < seg.second; ++k) {
                int epoch = static_cast<int>(std::floor((double)k / ppgSR / data.scoringEpochSec));
                if (epoch >= 0 && epoch < (int)data.sleepStages.size())
                    results[i].sleep_stages.push_back((int)data.sleepStages[epoch]);
                else
                    results[i].sleep_stages.push_back(-1);
            }
        }

        // --- THE SHAVE LOGIC ---
        // If this bin ends exactly where the next one starts, subtract 1 
        // to match MATLAB's inclusive boundary adjustment
        if (i < bin_count - 1 && !results[i].ppg.empty()) {
            // Only shave if the next bin actually has data and follows immediately
            if (!bin_assignments[i + 1].empty() && merged.back().second == bin_assignments[i + 1][0].first) {
                if (results[i].ppg.size() > 0) results[i].ppg.pop_back();
                if (results[i].ecg.size() > 0) results[i].ecg.pop_back();
                if (results[i].sleep_stages.size() > 0) results[i].sleep_stages.pop_back();
            }
        }
    }

    return results;
}

// ============================================================================
// MAIN APPLICATION
// ============================================================================

int main() {
    std::vector<ProjectConfig> projects;
    std::ifstream cfgFile("config.csv");
    if (!cfgFile.is_open()) { std::cerr << "Config.csv missing!" << std::endl; return 1; }

    std::string line; std::getline(cfgFile, line); // skip header
    while (std::getline(cfgFile, line)) {
        std::stringstream ss(line); std::string col; std::vector<std::string> cols;
        while (std::getline(ss, col, ',')) cols.push_back(cleanField(col));
        if (cols.size() >= 13) projects.push_back({ cols[0], cols[4], cols[5], cols[6],
                                                   std::stod(cols[11]), cols[12].empty() ? 0 : std::stod(cols[12]) });
    }

    std::cout << "Select Dataset: \n";
    for (size_t i = 0; i < projects.size(); ++i) std::cout << i + 1 << ". " << projects[i].dataType << "\n";
    int choice; std::cin >> choice;
    if (choice < 1 || choice >(int)projects.size()) return 1;

    ProjectConfig sel = projects[choice - 1];
    if (!fs::exists(sel.annealedPath)) fs::create_directories(sel.annealedPath);

    for (const auto& entry : fs::directory_iterator(sel.binPath)) {
        if (entry.path().extension() == ".bin") {
            std::string id = entry.path().stem().string();
            std::cout << "Processing: " << id << std::endl;

            try {
                RawData raw = readStructuredBin(entry.path().string());
                auto noise = readNoiseCSV(sel.noisePath + "/" + id + "_noise_markings.csv");
                auto results = AnnealSegments(raw, noise, 1.0);



                std::ofstream out(sel.annealedPath + "/" + id + "_annealed.bin", std::ios::binary);
                uint64_t nBins = results.size();
                out.write((char*)&nBins, 8);
                for (auto& s : results) {
                    uint32_t pS = s.ppg.size(), eS = s.ecg.size(), sS = s.sleep_stages.size();
                    out.write((char*)&pS, 4); out.write((char*)s.ppg.data(), (uint64_t)pS * 8);
                    out.write((char*)&eS, 4); out.write((char*)s.ecg.data(), (uint64_t)eS * 8);
                    out.write((char*)&sS, 4); out.write((char*)s.sleep_stages.data(), (uint64_t)sS * 4);
                }
                std::cout << "  -> Saved " << nBins << " bins." << std::endl;
            }
            catch (const std::exception& e) { std::cerr << "  -> Error: " << e.what() << std::endl; }
        }
    }
    return 0;
}
