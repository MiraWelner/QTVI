#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <filesystem>
#include <cstdint>

namespace fs = std::filesystem;

// ============================================================================
// HARDCODED CONFIGURATION VALUES
// ============================================================================
const double CONF_TARGET_LENGTH_MINS = 1.0;
const double CONF_NOISE_EXPANSION_SEC = 15.0;
const double CONF_MIN_EXCLUSION_SEC = 5.0;

// ============================================================================
// STRUCTURES
// ============================================================================

struct ConfigData {
    std::string data_name;
    std::string bin_path;
    std::string noise_path;
    std::string out_path;
};

struct SegmentIndices {
    uint64_t start_idx;
    uint64_t end_idx;
};

struct Exclusion {
    size_t begin_idx;
    size_t end_idx;
    size_t bin_num;
};

struct GoodSection {
    size_t beg;
    size_t end;
    int move_dir;
    int move_flag;
    size_t original_bin_idx; 
};

struct FinalBin {
    std::vector<SegmentIndices> ppg_indices;
    std::vector<SegmentIndices> ecg_indices;
    std::vector<double> ppg_data;
    std::vector<double> ecg_data;
    std::vector<int> sleep_stages;
    double ppgSR;
    double ecgSR;
};

struct BinData {
    double ecgSR;
    double ppgSR;
    double scoring_epoch;
    std::vector<double> ecg;
    std::vector<double> ppg;
    std::vector<double> sleepStages;
};


// ============================================================================
// HELPERS
// ============================================================================

std::string trim(const std::string& s) {
    auto first = s.find_first_not_of(" \t\r\n\"");
    if (first == std::string::npos) return "";
    auto last = s.find_last_not_of(" \t\r\n\"");
    return s.substr(first, (last - first + 1));
}

std::vector<ConfigData> readConfigCSV(const std::string& filename) {
    std::vector<ConfigData> configs;
    std::ifstream file(filename);
    if (!file.is_open()) return configs;

    std::string line;
    std::getline(file, line); // Skip Header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string col;
        std::vector<std::string> row;
        while (std::getline(ss, col, ',')) row.push_back(trim(col));

        if (row.size() >= 7) {
            ConfigData cd;
            cd.data_name = row[0];
            cd.bin_path = row[4];
            cd.noise_path = row[5];
            cd.out_path = row[6];
            configs.push_back(cd);
        }
    }
    return configs;
}

size_t closest_idx(const std::vector<double>& time_vec, double target) {
    auto it = std::lower_bound(time_vec.begin(), time_vec.end(), target);
    if (it == time_vec.end()) return time_vec.size() - 1;
    if (it == time_vec.begin()) return 0;
    if (std::abs(*it - target) < std::abs(*(it - 1) - target))
        return std::distance(time_vec.begin(), it);
    return std::distance(time_vec.begin(), it - 1);
}

std::vector<SegmentIndices> MergeSegments(std::vector<SegmentIndices> segs) {
    if (segs.empty()) return {};
    std::sort(segs.begin(), segs.end(), [](const SegmentIndices& a, const SegmentIndices& b) {
        return a.start_idx < b.start_idx;
        });
    std::vector<SegmentIndices> merged;
    SegmentIndices curr = segs[0];
    for (size_t i = 1; i < segs.size(); ++i) {
        if (segs[i].start_idx <= curr.end_idx + 1) {
            curr.end_idx = std::max(curr.end_idx, segs[i].end_idx);
        }
        else {
            merged.push_back(curr);
            curr = segs[i];
        }
    }
    merged.push_back(curr);
    return merged;
}

// ============================================================================
// FILE I/O
// ============================================================================

BinData readStructuredBin(const std::string& path) {
    BinData data;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("File error");

    double dummy;
    uint64_t nEcg, nEcg2, nEcg3, nPpg, nSleep;

    file.read((char*)&dummy, 8);         // ecgSR (ignore)
    file.read((char*)&data.ppgSR, 8);    // ppgSR (save)
    file.read((char*)&dummy, 8);         // scoring_epoch (ignore)

    file.read((char*)&nEcg, 8);
    file.read((char*)&nEcg2, 8);
    file.read((char*)&nEcg3, 8);
    file.read((char*)&nPpg, 8);
    file.read((char*)&nSleep, 8);

    // CRITICAL: Skip exactly (nEcg + nEcg2 + nEcg3) doubles
    file.seekg((nEcg + nEcg2 + nEcg3) * 8, std::ios::cur);

    data.ppg.resize(nPpg);
    file.read((char*)data.ppg.data(), nPpg * 8);
    std::cout << data.ppg.size() << std::endl;
    return data;
}

std::vector<std::pair<double, double>> readNoiseCSV(const std::string& path) {
    std::vector<std::pair<double, double>> noise;
    std::ifstream file(path);
    if (!file.is_open()) return noise;
    std::string line;
    std::getline(file, line);
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string col;
        std::vector<std::string> row;
        while (std::getline(ss, col, ',')) row.push_back(trim(col));
        if (row.size() >= 4) {
            try { noise.push_back({ std::stod(row[2]), std::stod(row[3]) }); }
            catch (...) {}
        }
    }
    return noise;
}

// ============================================================================
// ANNEALING
// ============================================================================

std::vector<FinalBin> AnnealSegments(const BinData& data, const std::vector<std::pair<double, double>>& rawNoise) {
    const double targetLength = 1.0;
    const double min_bin_size_mins = 0.5;
    const double min_exclusion_sec = 5.0;

    double ppgSR = data.ppgSR;
    // MATLAB: bin_size_samples = ppgSampleRate * 60 * targetLength;
    size_t bin_size_samples = (size_t)(ppgSR * 60 * targetLength);

    // 1. Bin Count
    size_t bin_count;
    double remainder_mins = (double)(data.ppg.size() % bin_size_samples) / ppgSR / 60.0;
    if (remainder_mins < min_bin_size_mins) {
        bin_count = (size_t)std::floor((double)data.ppg.size() / (double)bin_size_samples);
    }
    else {
        bin_count = (size_t)std::ceil((double)data.ppg.size() / (double)bin_size_samples);
    }

    // 2. Exact MATLAB bin_breaks translation
    // MATLAB: bin_breaks = (bin_size_samples + 1:bin_size_samples:length(ppg));
    std::vector<uint64_t> bin_breaks;
    for (uint64_t v = bin_size_samples + 1; v <= data.ppg.size(); v += bin_size_samples) {
        bin_breaks.push_back(v);
    }

    // MATLAB: if length(bin_breaks) < bin_count ... [bin_breaks length(ppg)]
    if (bin_breaks.size() < bin_count) {
        bin_breaks.push_back(data.ppg.size());
    }
    else {
        bin_breaks.back() = data.ppg.size();
    }

    // 3. Exclusions (Keep only if duration >= 5s)
    std::vector<std::pair<double, double>> exclusions_sec;
    for (auto& n : rawNoise) {
        if ((n.second - n.first) >= min_exclusion_sec) {
            exclusions_sec.push_back(n);
        }
    }

    // 4. Annealing Loop
    std::vector<std::vector<SegmentIndices>> final_bin_idx(bin_count);
    std::vector<SegmentIndices> temp_bin;
    size_t current_bin = 0;

    for (size_t b = 0; b < bin_count; ++b) {
        // MATLAB: bin_end = bin_breaks(cur_bin); bin_begin = bin_end - bin_size_samples;
        uint64_t m_end = bin_breaks[b];               // 1-based end (e.g., 15361)
        uint64_t m_begin = m_end - bin_size_samples;  // 1-based start (e.g., 1)

        // Convert to 0-based for C++ vectors
        uint64_t b_start = m_begin - 1;
        uint64_t b_end = m_end - 1;

        // MATLAB: bin_half = bin_end - (bin_size_samples / 2);
        double bin_half = (double)m_end - ((double)bin_size_samples / 2.0);

        // Find noise indices for this bin
        std::vector<std::pair<uint64_t, uint64_t>> bin_noise;
        for (auto& ex : exclusions_sec) {
            // MATLAB uses closest_idx which is similar to round(t * SR)
            uint64_t ns = (uint64_t)std::round(ex.first * ppgSR);
            uint64_t ne = (uint64_t)std::round(ex.second * ppgSR);

            uint64_t is = std::max(ns, b_start);
            uint64_t ie = std::min(ne, b_end);
            if (is <= ie) bin_noise.push_back({ is, ie });
        }
        std::sort(bin_noise.begin(), bin_noise.end());

        // Good sections in bin
        std::vector<SegmentIndices> good_frags;
        uint64_t curr = b_start;
        for (auto& n : bin_noise) {
            if (n.first > curr) good_frags.push_back({ curr, n.first - 1 });
            curr = n.second + 1;
        }
        if (curr <= b_end) good_frags.push_back({ curr, b_end });

        for (auto& seg : good_frags) {
            double len_m = (double)(seg.end_idx - seg.start_idx + 1) / ppgSR / 60.0;
            if (len_m < min_bin_size_mins) {
                // MATLAB uses 1-based index for midpoint comparison: seg.end_idx + 1
                int move_dir = ((double)seg.end_idx + 1 < bin_half) ? 1 : 2;
                if (b == 0) move_dir = 2;
                if (b == bin_count - 1) move_dir = 1;

                if (move_dir == 1) {
                    if (current_bin > 0) final_bin_idx[current_bin - 1].push_back(seg);
                }
                else {
                    temp_bin.push_back(seg);
                }
            }
            else {
                for (auto& t : temp_bin) final_bin_idx[current_bin].push_back(t);
                temp_bin.clear();
                final_bin_idx[current_bin].push_back(seg);
                current_bin++;
            }
        }
    }

    // 5. Final Extraction
    std::vector<FinalBin> results(bin_count);
    for (size_t i = 0; i < bin_count; ++i) {
        results[i].ppgSR = ppgSR;
        for (auto& idx : final_bin_idx[i]) {
            // Correctly inserts [start, end] inclusive
            results[i].ppg_data.insert(results[i].ppg_data.end(),
                data.ppg.begin() + idx.start_idx,
                data.ppg.begin() + idx.end_idx + 1);
        }
    }
    return results;
}


// ============================================================================
// MAIN EXECUTION
// ============================================================================

int main() {
    auto configs = readConfigCSV("config.csv");
    if (configs.empty()) {
        std::cerr << "Error: config.csv not found or empty." << std::endl;
        return 1;
    }

    std::cout << "Available Datasets:\n";
    for (size_t i = 0; i < configs.size(); ++i) std::cout << i + 1 << ". " << configs[i].data_name << "\n";

    int choice;
    std::cout << "Select Dataset: "; std::cin >> choice;
    if (choice < 1 || choice >(int)configs.size()) return 1;
    ConfigData sel = configs[choice - 1];

    if (!fs::exists(sel.out_path)) fs::create_directories(sel.out_path);

    std::cout << "\nScanning for .bin files in: " << sel.bin_path << "...\n";

    for (const auto& entry : fs::directory_iterator(sel.bin_path)) {
        if (entry.path().extension() == ".bin") {
            std::string id = entry.path().stem().string();
            std::string nFile = sel.noise_path + "/" + id + "_noise_markings.csv";
            std::cout << "Processing ID: " << id << std::endl;
            try {
                BinData data = readStructuredBin(entry.path().string());
                auto results = AnnealSegments(data, readNoiseCSV(nFile));
                std::ofstream out(sel.out_path + "/" + id + "_annealed.bin", std::ios::binary);
                uint64_t nBins = results.size();
                out.write(reinterpret_cast<char*>(&nBins), sizeof(uint64_t));
                for (auto& b : results) {
                    uint32_t pS = (uint32_t)b.ppg_data.size();
                    out.write(reinterpret_cast<char*>(&pS), sizeof(uint32_t));
                    out.write(reinterpret_cast<char*>(b.ppg_data.data()), pS * sizeof(double));
                }
                std::cout << "  -> Saved " << nBins << " annealed segments." << std::endl;
            }
            catch (const std::exception& e) { std::cerr << "  -> Failed: " << e.what() << std::endl; }
        }
    }
    std::cout << "\nDone." << std::endl;
    return 0;
}
