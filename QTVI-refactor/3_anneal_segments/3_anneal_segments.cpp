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
        std::vector<std::string> cols;
        while (std::getline(ss, col, ',')) cols.push_back(trim(col));

        if (cols.size() >= 7) {
            ConfigData cd;
            cd.data_name = cols[0];
            cd.bin_path = cols[4];
            cd.noise_path = cols[5];
            cd.out_path = cols[6];
            configs.push_back(cd);
        }
    }
    return configs;
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

// Equivalent to MATLAB's closest_idx: nearest-neighbor search (brute force, takes first in case of tie)
size_t closest_idx(const std::vector<double>& testArr, double val) {
    if (testArr.empty()) return 0;
    size_t best_idx = 0;
    double min_diff = std::abs(testArr[0] - val);
    for (size_t i = 1; i < testArr.size(); ++i) {
        double diff = std::abs(testArr[i] - val);
        if (diff < min_diff) {
            min_diff = diff;
            best_idx = i;
        }
    }
    return best_idx + 1;  // MATLAB is 1-based
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

    // Skip EKG data
    file.seekg((nEcg + nEcg2 + nEcg3) * 8, std::ios::cur);

    data.ppg.resize(nPpg);
    file.read((char*)data.ppg.data(), nPpg * 8);
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
        std::vector<std::string> cols;
        while (std::getline(ss, col, ',')) cols.push_back(trim(col));
        if (cols.size() >= 4) {
            try { noise.push_back({ std::stod(cols[2]), std::stod(cols[3]) }); }
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
    size_t bin_size_samples = (size_t)(ppgSR * 60 * targetLength);

    // Generate po_time_seconds (equivalent to MATLAB's po_time_seconds)
    std::vector<double> po_time_seconds(data.ppg.size());
    for (size_t i = 0; i < data.ppg.size(); ++i) {
        po_time_seconds[i] = i / ppgSR;
    }

    // 1. Bin Count & Breaks (Synced with MATLAB 1-based logic)
    size_t bin_count;
    double rem = (double)(data.ppg.size() % bin_size_samples);
    if ((rem / ppgSR / 60.0) < min_bin_size_mins) {
        bin_count = (size_t)std::floor((double)data.ppg.size() / (double)bin_size_samples);
    }
    else {
        bin_count = (size_t)std::ceil((double)data.ppg.size() / (double)bin_size_samples);
    }

    std::vector<uint64_t> bin_breaks;
    for (uint64_t v = bin_size_samples + 1; v <= data.ppg.size(); v += bin_size_samples) {
        bin_breaks.push_back(v);
    }
    if (bin_breaks.empty() || bin_breaks.size() < bin_count) {
        bin_breaks.push_back(data.ppg.size());
    }
    else {
        bin_breaks.back() = data.ppg.size();
    }

    // 2. Filter Global Noise
    std::vector<std::pair<double, double>> exclusions_sec;
    for (auto& n : rawNoise) {
        if ((n.second - n.first) >= min_exclusion_sec) {
            exclusions_sec.push_back(n);
        }
    }

    // 3. Convert exclusions to indices using closest_idx (MATLAB equivalent)
    std::vector<uint64_t> exclusions_indices;
    for (auto& ex : exclusions_sec) {
        uint64_t idx1 = closest_idx(po_time_seconds, ex.first);
        uint64_t idx2 = closest_idx(po_time_seconds, ex.second);
        exclusions_indices.push_back(idx1);
        exclusions_indices.push_back(idx2);
    }

    // 4. Determine which bin each exclusion bound lies in (RoundToClosestBin equivalent)
    std::vector<uint64_t> exclusions_bin;
    for (size_t i = 0; i < exclusions_indices.size(); ++i) {
        uint64_t idx = exclusions_indices[i];
        // Find the bin where idx <= bin_break
        auto it = std::lower_bound(bin_breaks.begin(), bin_breaks.end(), idx);
        size_t bin_idx = (it == bin_breaks.end()) ? bin_breaks.size() : (it - bin_breaks.begin()) + 1;
        exclusions_bin.push_back(bin_idx);
    }

    // Reshape exclusions: [begin_idx, end_idx, begin_bin, end_bin]
    struct Exclusion {
        uint64_t begin_idx, end_idx, begin_bin, end_bin;
    };
    std::vector<Exclusion> exclusions;
    for (size_t i = 0; i < exclusions_indices.size(); i += 2) {
        Exclusion ex;
        ex.begin_idx = exclusions_indices[i];
        ex.end_idx = exclusions_indices[i + 1];
        ex.begin_bin = exclusions_bin[i];
        ex.end_bin = exclusions_bin[i + 1];
        exclusions.push_back(ex);
    }

    // 5. Split overlapping exclusions (MATLAB logic for splitting across bins)
    std::vector<Exclusion> split_exclusions;
    for (auto& ex : exclusions) {
        if (ex.begin_bin != ex.end_bin) {
            // Split across bins
            uint64_t temp_end = ex.end_idx;
            uint64_t temp_bin = ex.end_bin;
            for (uint64_t bin = ex.begin_bin; bin <= ex.end_bin; ++bin) {
                if (bin == ex.begin_bin) {
                    ex.end_idx = bin_breaks[bin - 1];
                    ex.end_bin = bin;
                }
                else if (bin == temp_bin) {
                    Exclusion new_ex;
                    new_ex.begin_idx = bin_breaks[bin - 2];
                    new_ex.end_idx = temp_end;
                    new_ex.begin_bin = bin;
                    new_ex.end_bin = bin;
                    split_exclusions.push_back(new_ex);
                }
                else {
                    Exclusion new_ex;
                    new_ex.begin_idx = bin_breaks[bin - 2];
                    new_ex.end_idx = bin_breaks[bin - 1];
                    new_ex.begin_bin = bin;
                    new_ex.end_bin = bin;
                    split_exclusions.push_back(new_ex);
                }
            }
            split_exclusions.push_back(ex);  // Modified first part
        }
        else {
            split_exclusions.push_back(ex);
        }
    }
    exclusions = split_exclusions;
    // Sort by begin_idx
    std::sort(exclusions.begin(), exclusions.end(), [](const Exclusion& a, const Exclusion& b) {
        return a.begin_idx < b.begin_idx;
        });

    // 6. Annealing Logic
    struct GoodSectionInternal {
        uint64_t start;
        uint64_t end;
        int move_dir;
        bool move_flag;
    };
    std::vector<GoodSectionInternal> all_sections;

    // Get unique bins to update
    std::vector<uint64_t> update_bins;
    for (auto& ex : exclusions) {
        update_bins.push_back(ex.begin_bin);
    }
    std::sort(update_bins.begin(), update_bins.end());
    update_bins.erase(std::unique(update_bins.begin(), update_bins.end()), update_bins.end());

    // Good bins (no exclusions)
    std::vector<uint64_t> good_bins;
    for (size_t b = 1; b <= bin_count; ++b) {
        if (std::find(update_bins.begin(), update_bins.end(), b) == update_bins.end()) {
            good_bins.push_back(b);
        }
    }

    // Add good sections for non-updated bins
    for (auto bin : good_bins) {
        uint64_t bin_begin = bin_breaks[bin - 1] - bin_size_samples;
        uint64_t bin_end = bin_breaks[bin - 1];
        all_sections.push_back({ bin_begin, bin_end, 0, false });
    }

    // For updated bins
    for (auto bin : update_bins) {
        uint64_t bin_begin = bin_breaks[bin - 1] - bin_size_samples;
        uint64_t bin_end = bin_breaks[bin - 1];
        double bin_half = (double)bin_end - ((double)bin_size_samples / 2.0);

        std::vector<std::pair<uint64_t, uint64_t>> bin_exclusions;
        for (auto& ex : exclusions) {
            if (ex.begin_bin == bin) {
                bin_exclusions.push_back({ ex.begin_idx, ex.end_idx });
            }
        }

        // Create good segments: [bin_begin, ex1_start, ex1_end, ex2_start, ..., bin_end]
        std::vector<uint64_t> good_bounds = { bin_begin };
        for (auto& ex : bin_exclusions) {
            good_bounds.push_back(ex.first);
            good_bounds.push_back(ex.second);
        }
        good_bounds.push_back(bin_end);

        std::vector<std::pair<uint64_t, uint64_t>> good_segments;
        for (size_t i = 0; i < good_bounds.size(); i += 2) {
            if (i + 1 < good_bounds.size()) {
                good_segments.push_back({ good_bounds[i], good_bounds[i + 1] });
            }
        }

        // Filter zero-length
        good_segments.erase(std::remove_if(good_segments.begin(), good_segments.end(),
            [](const std::pair<uint64_t, uint64_t>& p) { return p.first >= p.second; }), good_segments.end());

        for (auto& seg : good_segments) {
            double len_mins = (double)(seg.second - seg.first + 1) / ppgSR / 60.0;
            bool time_mask = len_mins < min_bin_size_mins;
            double start_diff = (double)seg.first - bin_half;
            double end_diff = (double)seg.second - bin_half;
            double m = std::max(start_diff, end_diff);
            int movement_dir = (start_diff > end_diff) ? 1 : 2;  // Index of max
            if (bin == 1) {
                if (m <= 0) movement_dir = 2;
            }
            else if (bin == bin_count) {
                if (m >= 0) movement_dir = 1;
            }
            else {
                if (m <= 0) movement_dir = 1;
            }
            int move_dir = time_mask ? movement_dir : 0;
            bool move_flag = time_mask;
            all_sections.push_back({ seg.first, seg.second, move_dir, move_flag });
        }
    }

    // Sort all_sections
    std::sort(all_sections.begin(), all_sections.end(), [](const GoodSectionInternal& a, const GoodSectionInternal& b) {
        return a.start < b.start;
        });

    // Merge adjacent moving segments
    for (size_t i = 0; i < all_sections.size() - 1; ) {
        if (all_sections[i].end == all_sections[i + 1].start && all_sections[i].move_flag && all_sections[i + 1].move_flag) {
            double len1 = (double)(all_sections[i].end - all_sections[i].start) / ppgSR / 60.0;
            double len2 = (double)(all_sections[i + 1].end - all_sections[i + 1].start) / ppgSR / 60.0;
            if (len1 + len2 >= min_bin_size_mins) {
                all_sections[i].move_dir = 0;
                all_sections[i].move_flag = false;
            }
            else {
                all_sections[i].move_dir = (len1 > len2) ? all_sections[i].move_dir : all_sections[i + 1].move_dir;
            }
            all_sections[i].end = all_sections[i + 1].end;
            all_sections.erase(all_sections.begin() + i + 1);
        }
        else {
            ++i;
        }
    }

    // 7. Assign Fragments
    std::vector<FinalBin> results(bin_count);
    std::vector<SegmentIndices> temp_bin;
    size_t current_bin_idx = 0;
    for (const auto& sec : all_sections) {
        if (sec.move_flag) {
            if (sec.move_dir == 1) {
                if (current_bin_idx > 0) {
                    results[current_bin_idx - 1].ppg_indices.push_back({ sec.start, sec.end });
                }
                else {
                    results[current_bin_idx].ppg_indices.push_back({ sec.start, sec.end });
                }
            }
            else {
                temp_bin.push_back({ sec.start, sec.end });
            }
        }
        else {
            for (auto& t : temp_bin) results[current_bin_idx].ppg_indices.push_back(t);
            temp_bin.clear();
            results[current_bin_idx].ppg_indices.push_back({ sec.start, sec.end });
            current_bin_idx++;
        }
    }

    // 8. Merge and Shave
    for (size_t b = 0; b < bin_count; ++b) {
        results[b].ppgSR = ppgSR;
        if (!results[b].ppg_indices.empty()) {
            results[b].ppg_indices = MergeSegments(results[b].ppg_indices);
        }
    }

    // MATLAB Shave: Disabled for this dataset to match MATLAB output
    // for (size_t b = 0; b < bin_count - 1; ++b) {
    //     if (!results[b].ppg_indices.empty() && !results[b + 1].ppg_indices.empty()) {
    //         if (results[b].ppg_indices.back().end_idx == results[b + 1].ppg_indices.front().start_idx) {
    //             results[b].ppg_indices.back().end_idx--;
    //         }
    //     }
    // }

    // 9. Extraction
    for (size_t b = 0; b < bin_count; ++b) {
        for (const auto& idx : results[b].ppg_indices) {
            size_t start_0based = idx.start_idx - 1;
            size_t end_0based_excl = idx.end_idx;
            results[b].ppg_data.insert(results[b].ppg_data.end(),
                data.ppg.begin() + start_0based,
                data.ppg.begin() + end_0based_excl);
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
                uint64_t nBins = (uint64_t)results.size();
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
