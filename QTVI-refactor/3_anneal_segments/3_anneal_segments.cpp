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
    std::vector<double> sleep_stages;
};

// ============================================================================
// HELPERS
// ============================================================================

std::string cleanField(std::string s) {
    s.erase(0, s.find_first_not_of(" \t\r\n\""));
    s.erase(s.find_last_not_of(" \t\r\n\"") + 1);
    return s;
}

uint64_t closest_idx(double target_time, double sr) {
    return static_cast<uint64_t>(std::round(target_time * sr)) + 1;
}

std::vector<std::pair<uint64_t, uint64_t>> MergeSegments(std::vector<std::pair<uint64_t, uint64_t>> segs) {
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

// ============================================================================
// CORE ANNEALING (REPLICA OF AnnealSegments.m)
// ============================================================================

std::vector<FinalSegment> AnnealSegments(const RawData& data, const std::vector<std::pair<uint64_t, uint64_t>>& noise, double targetLenMins) {
    const double ppgSR = data.ppgSR;
    const uint64_t bin_size_samples = static_cast<uint64_t>(ppgSR * 60.0 * targetLenMins);
    const double min_bin_size_mins = targetLenMins / 2.0;
    const uint64_t total_len = data.ppg.size();

    // 1. Bin Count Logic (Exact MATLAB Replica)
    int bin_count;
    double remainder_mins = (double)(total_len % bin_size_samples) / ppgSR / 60.0;
    if (remainder_mins < min_bin_size_mins) {
        bin_count = (int)std::floor((double)total_len / bin_size_samples);
    }
    else {
        bin_count = (int)std::ceil((double)total_len / bin_size_samples);
    }

    // 2. Bin Breaks (Exact MATLAB Replica)
    std::vector<uint64_t> bin_breaks;
    for (uint64_t b = bin_size_samples + 1; b <= total_len; b += bin_size_samples) {
        bin_breaks.push_back(b);
    }
    if (bin_breaks.size() < (size_t)bin_count) {
        bin_breaks.push_back(total_len);
    }
    else if (!bin_breaks.empty()) {
        bin_breaks.back() = total_len;
    }

    // 3. Final Indices (Exact MATLAB Replica)
    struct BinIdx { std::vector<std::pair<uint64_t, uint64_t>> po; };
    std::vector<BinIdx> final_bin_idx(bin_count);

    for (int b = 1; b <= bin_count; ++b) {
        uint64_t b_end = bin_breaks[b - 1];
        uint64_t b_start = b_end - bin_size_samples;
        if (b == 1) b_start = 1;
        final_bin_idx[b - 1].po.push_back({ b_start, b_end });
    }

    // 4. Correct Overlaps (MATLAB: final_bin_idx{i}.po(end,2) = ... - 1)
    for (int i = 0; i < bin_count - 1; ++i) {
        if (!final_bin_idx[i].po.empty() && !final_bin_idx[i + 1].po.empty()) {
            if (final_bin_idx[i].po.back().second == final_bin_idx[i + 1].po.front().first) {
                final_bin_idx[i].po.back().second--;
            }
        }
    }

    // 5. Data Collection
    std::vector<FinalSegment> results(bin_count);
    for (int i = 0; i < bin_count; ++i) {
        for (auto& seg : final_bin_idx[i].po) {
            // PPG
            for (uint64_t k = seg.first; k <= seg.second; ++k) {
                results[i].ppg.push_back(data.ppg[k - 1]);
            }

            // ECG (using closest_idx like MATLAB)
            uint64_t e_s = closest_idx((double)(seg.first - 1) / ppgSR, data.ecgSR);
            uint64_t e_e = closest_idx((double)(seg.second - 1) / ppgSR, data.ecgSR);
            for (uint64_t k = e_s; k <= e_e && k <= data.ecg.size(); ++k) {
                results[i].ecg.push_back(data.ecg[k - 1]);
            }

            // Sleep Stages (Exact MATLAB boundary filter)
            double time_start = (double)(seg.first - 1) / ppgSR;
            double time_end = (double)(seg.second - 1) / ppgSR;
            for (size_t s = 0; s < data.sleepStages.size(); ++s) {
                double epoch_end_time = (s + 1) * data.scoringEpochSec;
                if (epoch_end_time >= time_start && epoch_end_time <= time_end) {
                    results[i].sleep_stages.push_back(data.sleepStages[s]);
                }
            }
        }
    }

    return results;
}

// ============================================================================
// FILE I/O
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

int main() {
    std::vector<ProjectConfig> projects;
    std::ifstream cfg("config.csv");
    if (!cfg.is_open()) {
        std::cerr << "Could not open config.csv" << std::endl;
        return 1;
    }

    std::string line;
    std::getline(cfg, line); // Header
    while (std::getline(cfg, line)) {
        std::stringstream ss(line);
        std::string col;
        std::vector<std::string> cols;
        while (std::getline(ss, col, ',')) cols.push_back(cleanField(col));
        if (cols.size() >= 13) {
            projects.push_back({ cols[0], cols[4], cols[5], cols[6], std::stod(cols[11]), std::stod(cols[12]) });
        }
    }

    std::cout << "Select Dataset:\n";
    for (size_t i = 0; i < projects.size(); ++i) std::cout << i + 1 << ". " << projects[i].dataType << "\n";

    int choice;
    std::cin >> choice;
    if (choice < 1 || choice >(int)projects.size()) return 1;
    ProjectConfig sel = projects[choice - 1];

    if (!fs::exists(sel.annealedPath)) fs::create_directories(sel.annealedPath);

    for (const auto& entry : fs::directory_iterator(sel.binPath)) {
        if (entry.path().extension() == ".bin") {
            std::string id = entry.path().stem().string();
            try {
                RawData raw = readStructuredBin(entry.path().string());
                auto results = AnnealSegments(raw, {}, 1.0); // No noise expansion needed if noise empty

                std::ofstream out(sel.annealedPath + "/" + id + "_annealed.bin", std::ios::binary);
                uint64_t nB = results.size();
                out.write((char*)&nB, 8);

                for (auto& s : results) {
                    uint64_t pS = s.ppg.size();
                    uint64_t eS = s.ecg.size();
                    uint64_t sS = s.sleep_stages.size();

                    out.write((char*)&pS, 8);
                    out.write((char*)s.ppg.data(), pS * 8);

                    out.write((char*)&eS, 8);
                    out.write((char*)s.ecg.data(), eS * 8);

                    out.write((char*)&sS, 8);
                    out.write((char*)s.sleep_stages.data(), sS * 8);
                }
                std::cout << "Processed: " << id << " (" << nB << " bins)" << std::endl;
            }
            catch (const std::exception& e) {
                std::cerr << "Error processing " << id << ": " << e.what() << std::endl;
            }
        }
    }

    return 0;
}
