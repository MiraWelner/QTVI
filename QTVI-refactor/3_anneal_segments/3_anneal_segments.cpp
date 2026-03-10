#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <cstdint>
#include <set>

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

// MATLAB: closest_idx finds the index in a time vector nearest to target_time.
// Since time vectors are uniform (0:1/sr:...), this simplifies to round(target_time * sr) + 1
// The +1 is for MATLAB 1-based indexing, which we keep here because the
// annealing logic uses 1-based indices throughout (matching MATLAB exactly).
uint64_t closest_idx(double target_time, double sr) {
    return static_cast<uint64_t>(std::round(target_time * sr)) + 1;
}

std::vector<std::pair<uint64_t, uint64_t>> MergeSegments(
    std::vector<std::pair<uint64_t, uint64_t>> segs) {
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

// Round each index to the bin it belongs to (1-based bin number)
std::vector<int> RoundToClosestBin(const std::vector<uint64_t>& bin_breaks,
    const std::vector<uint64_t>& indices) {
    std::vector<int> result(indices.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        int best_bin = 1;
        for (size_t b = 0; b < bin_breaks.size(); ++b) {
            if (indices[i] <= bin_breaks[b]) {
                best_bin = static_cast<int>(b) + 1;
                break;
            }
            best_bin = static_cast<int>(b) + 1;
        }
        result[i] = best_bin;
    }
    return result;
}

// ============================================================================
// NOISE FILE READER
// Reads the binary format from NoiseManager::exportBinary()
// Returns Nx2 vector of [startSec, endSec] pairs (like MATLAB's noiseSEG)
// ============================================================================

std::vector<std::pair<double, double>> readNoiseBin(const std::string& path) {
    std::vector<std::pair<double, double>> noiseSEG;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return noiseSEG;

    uint64_t count = 0;
    file.read(reinterpret_cast<char*>(&count), 8);

    for (uint64_t i = 0; i < count; ++i) {
        double row[6];
        file.read(reinterpret_cast<char*>(row), 6 * sizeof(double));
        // row[2] = startSec, row[3] = endSec
        noiseSEG.push_back({ row[2], row[3] });
    }
    return noiseSEG;
}

// ============================================================================
// CORE ANNEALING — Full port of AnnealSegments.m
// ============================================================================

std::vector<FinalSegment> AnnealSegments(
    const RawData& data,
    const std::vector<std::pair<double, double>>& noiseSEG,
    double targetLenMins)
{
    const double ppgSR = data.ppgSR;
    const double ecgSR = data.ecgSR;
    const uint64_t bin_size_samples = static_cast<uint64_t>(ppgSR * 60.0 * targetLenMins);
    const double min_bin_size_mins = targetLenMins / 2.0;
    const double min_exclusion_bin_size_seconds = 5.0;
    const uint64_t total_len = data.ppg.size();

    // 1. Bin Count
    int bin_count;
    double remainder_mins = (double)(total_len % bin_size_samples) / ppgSR / 60.0;
    if (remainder_mins < min_bin_size_mins)
        bin_count = (int)std::floor((double)total_len / bin_size_samples);
    else
        bin_count = (int)std::ceil((double)total_len / bin_size_samples);

    if (bin_count <= 0) return {};

    // 2. Bin Breaks (1-based, matching MATLAB)
    std::vector<uint64_t> bin_breaks;
    for (uint64_t b = bin_size_samples + 1; b <= total_len; b += bin_size_samples)
        bin_breaks.push_back(b);
    if ((int)bin_breaks.size() < bin_count)
        bin_breaks.push_back(total_len);
    else if (!bin_breaks.empty())
        bin_breaks.back() = total_len;

    // 3. Filter exclusions by minimum length
    std::vector<std::pair<double, double>> exclusions_seconds;
    for (const auto& seg : noiseSEG) {
        if ((seg.second - seg.first) >= min_exclusion_bin_size_seconds)
            exclusions_seconds.push_back(seg);
    }

    // 4. Convert exclusion times to PPG sample indices (1-based)
    std::vector<uint64_t> exclusions_indexs_flat;
    for (const auto& seg : exclusions_seconds) {
        exclusions_indexs_flat.push_back(closest_idx(seg.first, ppgSR));
        exclusions_indexs_flat.push_back(closest_idx(seg.second, ppgSR));
    }

    // 5. Determine which bin each exclusion boundary falls in
    std::vector<int> exclusions_bin = RoundToClosestBin(bin_breaks, exclusions_indexs_flat);

    // Reshape to Nx2 pairs + Nx2 bin assignments
    // Format: [begin_idx, end_idx, begin_bin, end_bin]
    struct Exclusion { uint64_t idx_start, idx_end; int bin_start, bin_end; };
    std::vector<Exclusion> exclusions;
    for (size_t i = 0; i < exclusions_seconds.size(); ++i) {
        exclusions.push_back({
            exclusions_indexs_flat[i * 2],
            exclusions_indexs_flat[i * 2 + 1],
            exclusions_bin[i * 2],
            exclusions_bin[i * 2 + 1]
            });
    }

    // 6. Split exclusions that span multiple bins
    for (size_t i = 0; i < exclusions.size(); ++i) {
        if (exclusions[i].bin_start != exclusions[i].bin_end) {
            uint64_t temp_end = exclusions[i].idx_end;
            int temp_bin_end = exclusions[i].bin_end;

            for (int bin = exclusions[i].bin_start; bin <= temp_bin_end; ++bin) {
                if (bin == exclusions[i].bin_start) {
                    exclusions[i].idx_end = bin_breaks[bin - 1];
                    exclusions[i].bin_end = bin;
                }
                else if (bin == temp_bin_end) {
                    exclusions.push_back({
                        bin_breaks[bin - 2], temp_end, bin, bin
                        });
                }
                else {
                    exclusions.push_back({
                        bin_breaks[bin - 2], bin_breaks[bin - 1], bin, bin
                        });
                }
            }
        }
    }

    // Sort by start index, keep only [idx_start, idx_end, bin] (bin_start == bin_end now)
    std::sort(exclusions.begin(), exclusions.end(),
        [](const Exclusion& a, const Exclusion& b) { return a.idx_start < b.idx_start; });

    // 7. Determine which bins need updating
    std::set<int> update_bins_set;
    for (const auto& ex : exclusions) update_bins_set.insert(ex.bin_start);
    std::vector<int> update_bins(update_bins_set.begin(), update_bins_set.end());

    // Good bins (no exclusions)
    std::set<int> all_bins_set;
    for (int b = 1; b <= bin_count; ++b) all_bins_set.insert(b);
    std::vector<int> good_bins;
    for (int b : all_bins_set) {
        if (update_bins_set.find(b) == update_bins_set.end())
            good_bins.push_back(b);
    }

    // 8. Build good_sections: [begin, end, movement_dir, move_flag]
    struct Section { uint64_t begin, end; int dir; int flag; };
    std::vector<Section> good_sections;

    // Add sections from bins without exclusions
    for (int cur_bin : good_bins) {
        uint64_t bin_begin = bin_breaks[cur_bin - 1] - bin_size_samples;
        uint64_t bin_end = bin_breaks[cur_bin - 1];
        good_sections.push_back({ bin_begin, bin_end, 0, 0 });
    }

    // Add sections from bins with exclusions
    for (int cur_bin : update_bins) {
        uint64_t bin_begin = bin_breaks[cur_bin - 1] - bin_size_samples;
        uint64_t bin_end = bin_breaks[cur_bin - 1];
        double bin_half = (double)(bin_end - (bin_size_samples / 2));

        // Get exclusions for this bin
        std::vector<std::pair<uint64_t, uint64_t>> bin_excl;
        for (const auto& ex : exclusions) {
            if (ex.bin_start == cur_bin)
                bin_excl.push_back({ ex.idx_start, ex.idx_end });
        }

        // Build good intervals: [bin_begin, excl1_start, excl1_end, excl2_start, ..., bin_end]
        std::vector<uint64_t> flat = { bin_begin };
        for (const auto& be : bin_excl) {
            flat.push_back(be.first);
            flat.push_back(be.second);
        }
        flat.push_back(bin_end);

        // Reshape to pairs
        std::vector<std::pair<uint64_t, uint64_t>> good;
        for (size_t j = 0; j + 1 < flat.size(); j += 2) {
            uint64_t len = flat[j + 1] - flat[j];
            if (len > 0) good.push_back({ flat[j], flat[j + 1] });
        }

        for (const auto& g : good) {
            double good_time_mins = ((double)(g.second - g.first) / ppgSR) / 60.0;
            bool too_small = good_time_mins < min_bin_size_mins;

            // Movement direction
            double max_val = std::max((double)g.first - bin_half, (double)g.second - bin_half);
            int movement_dir;
            if (cur_bin == 1) {
                movement_dir = (max_val <= 0) ? 2 : ((max_val > 0) ? 0 : 2);
            }
            else if (cur_bin == bin_count) {
                movement_dir = (max_val >= 0) ? 1 : ((max_val < 0) ? 0 : 1);
            }
            else {
                movement_dir = (max_val <= 0) ? 1 : 0;
            }

            int move_flag = too_small ? 1 : 0;
            if (!too_small) movement_dir = 0;

            // Skip small segments at edges of first/last bin
            if (cur_bin == 1 && too_small && max_val <= 0) continue;
            if (cur_bin == bin_count && too_small && max_val >= 0) continue;

            good_sections.push_back({ g.first, g.second, movement_dir, move_flag });
        }
    }

    // Sort by begin
    std::sort(good_sections.begin(), good_sections.end(),
        [](const Section& a, const Section& b) { return a.begin < b.begin; });

    // 9. Merge adjacent moving sections
    for (size_t i = 0; i + 1 < good_sections.size(); ) {
        if (good_sections[i].end == good_sections[i + 1].begin &&
            good_sections[i].flag != 0 && good_sections[i + 1].flag != 0) {

            double seg1_min = ((double)(good_sections[i].end - good_sections[i].begin) / ppgSR) / 60.0;
            double seg2_min = ((double)(good_sections[i + 1].end - good_sections[i + 1].begin) / ppgSR) / 60.0;

            if (seg1_min + seg2_min >= min_bin_size_mins) {
                good_sections[i].dir = 0;
                good_sections[i].flag = 0;
            }
            else {
                int idx = (seg1_min >= seg2_min) ? 0 : 1;
                good_sections[i].dir = good_sections[i + idx].dir;
                good_sections[i].flag = 1;
            }
            good_sections[i].end = good_sections[i + 1].end;
            good_sections.erase(good_sections.begin() + i + 1);
        }
        else {
            ++i;
        }
    }

    // 10. Assign sections to final bins
    struct BinIdx { std::vector<std::pair<uint64_t, uint64_t>> po; };
    std::vector<BinIdx> final_bin_idx;
    final_bin_idx.push_back(BinIdx());

    int current_bin = 0;
    std::vector<std::pair<uint64_t, uint64_t>> temp_bin;

    for (const auto& sec : good_sections) {
        if (sec.flag) {
            if (sec.dir == 1) { // Move left
                if (current_bin > 0)
                    final_bin_idx[current_bin - 1].po.push_back({ sec.begin, sec.end });
                else
                    final_bin_idx[current_bin].po.push_back({ sec.begin, sec.end });
            }
            else { // Move right
                temp_bin.push_back({ sec.begin, sec.end });
            }
        }
        else {
            for (const auto& tb : temp_bin)
                final_bin_idx[current_bin].po.push_back(tb);
            final_bin_idx[current_bin].po.push_back({ sec.begin, sec.end });
            temp_bin.clear();
            current_bin++;
            final_bin_idx.push_back(BinIdx());
        }
    }
    // Remove trailing empty bin
    while (!final_bin_idx.empty() && final_bin_idx.back().po.empty())
        final_bin_idx.pop_back();

    // 11. Merge overlapping segments within each bin
    for (auto& fb : final_bin_idx) {
        if (fb.po.size() > 1)
            fb.po = MergeSegments(fb.po);
    }

    // 12. Correct overlaps between adjacent bins
    for (size_t i = 0; i + 1 < final_bin_idx.size(); ++i) {
        if (!final_bin_idx[i].po.empty() && !final_bin_idx[i + 1].po.empty()) {
            if (final_bin_idx[i].po.back().second == final_bin_idx[i + 1].po.front().first)
                final_bin_idx[i].po.back().second--;
        }
    }

    // 13. Collect data from indices
    std::vector<FinalSegment> results(final_bin_idx.size());
    for (size_t i = 0; i < final_bin_idx.size(); ++i) {
        for (const auto& seg : final_bin_idx[i].po) {
            // PPG (1-based indices → 0-based data access)
            for (uint64_t k = seg.first; k <= seg.second && k <= total_len; ++k)
                results[i].ppg.push_back(data.ppg[k - 1]);

            // ECG
            uint64_t e_s = closest_idx((double)(seg.first - 1) / ppgSR, ecgSR);
            uint64_t e_e = closest_idx((double)(seg.second - 1) / ppgSR, ecgSR);
            for (uint64_t k = e_s; k <= e_e && k <= data.ecg.size(); ++k)
                results[i].ecg.push_back(data.ecg[k - 1]);

            // Sleep Stages
            double time_start = (double)(seg.first - 1) / ppgSR;
            double time_end = (double)(seg.second - 1) / ppgSR;
            for (size_t s = 0; s < data.sleepStages.size(); ++s) {
                double epoch_end_time = (s + 1) * data.scoringEpochSec;
                if (epoch_end_time >= time_start && epoch_end_time <= time_end)
                    results[i].sleep_stages.push_back(data.sleepStages[s]);
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
    file.read((char*)&nEcg1, 8);
    file.read((char*)&nEcg2, 8);
    file.read((char*)&nEcg3, 8);
    file.read((char*)&nPpg, 8);
    file.read((char*)&nSleep, 8);

    // FIX: Only read the first ECG channel (matches MATLAB's single 'ecg' variable)
    data.ecg.resize(nEcg1);
    file.read((char*)data.ecg.data(), nEcg1 * 8);

    // Skip channels 2 and 3
    file.seekg((nEcg2 + nEcg3) * 8, std::ios::cur);

    data.ppg.resize(nPpg);
    file.read((char*)data.ppg.data(), nPpg * 8);
    data.sleepStages.resize(nSleep);
    file.read((char*)data.sleepStages.data(), nSleep * 8);

    std::cout << "  ECG: " << nEcg1 << " samples (skipped ch2=" << nEcg2
        << " ch3=" << nEcg3 << "), PPG: " << nPpg << std::endl;

    return data;
}

// ============================================================================
// MAIN
// ============================================================================

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
            projects.push_back({
                cols[0], cols[4], cols[5], cols[6],
                std::stod(cols[11]), std::stod(cols[12])
                });
        }
    }

    std::cout << "Select Dataset:\n";
    for (size_t i = 0; i < projects.size(); ++i)
        std::cout << i + 1 << ". " << projects[i].dataType << "\n";

    int choice;
    std::cin >> choice;
    if (choice < 1 || choice >(int)projects.size()) return 1;
    ProjectConfig sel = projects[choice - 1];

    if (!fs::exists(sel.annealedPath))
        fs::create_directories(sel.annealedPath);

    for (const auto& entry : fs::directory_iterator(sel.binPath)) {
        if (entry.path().extension() != ".bin") continue;
        if (entry.path().string().find("_noise_markings") != std::string::npos) continue;

        std::string id = entry.path().stem().string();
        try {
            std::cout << "Processing: " << id << std::endl;
            RawData raw = readStructuredBin(entry.path().string());

            // Load noise markings if they exist
            std::string noise_path = sel.noisePath + "/" + id + "_noise_markings.bin";
            std::vector<std::pair<double, double>> noiseSEG;
            if (fs::exists(noise_path)) {
                noiseSEG = readNoiseBin(noise_path);
                std::cout << "  Loaded " << noiseSEG.size() << " noise segments" << std::endl;
            }
            else {
                std::cout << "  No noise file found, proceeding without exclusions" << std::endl;
            }

            auto results = AnnealSegments(raw, noiseSEG, 1.0);

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
            std::cout << "  Output: " << nB << " bins" << std::endl;
        }
        catch (const std::exception& e) {
            std::cerr << "Error processing " << id << ": " << e.what() << std::endl;
        }
    }

    std::cout << "Done." << std::endl;
    return 0;
}