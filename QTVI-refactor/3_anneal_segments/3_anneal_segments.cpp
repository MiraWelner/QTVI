#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <filesystem>
#include <map>
#include <cstdint>  // For uint64_t
#include <tuple>    // For std::tuple

namespace fs = std::filesystem;

// ============================================================================
// STRUCTURES (Matching Original Config)
// ============================================================================

struct ConfigData {
    std::string data_type;             // Col 0
    std::string bin_file_path;         // Col 4 (Where the .bin files are)
    std::string noise_data_path;       // Col 5 (Where the _noise_markings.csv are)
    std::string annealed_bin_path;     // Col 6 (Where the output goes)
    double ecg_sampling_rate;          // Col 11 (Used as default if not in .bin)
};

// ============================================================================
// STRUCTURES FOR LOADING FROM .BIN (Matching GUI's setFileSource Format)
// ============================================================================

// .bin structure from GUI code:
// - double ecgSR
// - double ppgSR
// - double scoring_epoch
// - uint64_t totalEcgSamples (Signal 1)
// - uint64_t totalSignal2Samples (Signal 2, ECG2)
// - uint64_t totalSignal3Samples (Signal 3, ECG3)
// - uint64_t totalPpgSamples
// - uint64_t totalSleepSamples
// Then data: ECG1, ECG2, ECG3, PPG, Sleep stages (all as doubles)

struct BinData {
    double ecgSR;
    double ppgSR;
    double scoring_epoch_sec;
    std::vector<double> ecg;       // Signal 1
    std::vector<double> ecg2;      // Signal 2
    std::vector<double> ecg3;      // Signal 3
    std::vector<double> ppg;
    std::vector<double> sleepStages;  // As doubles, cast to int if needed
};

// ============================================================================
// STRUCTURES FOR ANNEALING OUTPUT (Matching MATLAB annealedSegments)
// ============================================================================

struct SegmentIndices {
    size_t start_idx;
    size_t end_idx;
};

struct Exclusion {
    size_t begin_idx;
    size_t end_idx;
    size_t begin_bin;
    size_t end_bin;
};

struct GoodSection {
    size_t beg;
    size_t end;
    size_t move_dir;
    size_t move_flag;
};

struct FinalBin {
    std::vector<SegmentIndices> po_indices;
    std::vector<SegmentIndices> ecg_indices;
    std::vector<double> po_data;
    std::vector<double> ecg_data;
    std::vector<int> sleep_stages;
    std::vector<double> r_peaks;
    double ppgSampleRate;
    double ecgSampleRate;
    std::vector<SegmentIndices> ppg_bin_indexs;
    std::vector<SegmentIndices> ecg_bin_indexs;
};

// ============================================================================
// HELPERS (Matching Original C++)
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
    if (!file.is_open()) throw std::runtime_error("Could not open config: " + filename);

    std::string line;
    std::getline(file, line); // Skip Header

    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::string col;
        std::vector<std::string> row;
        while (std::getline(ss, col, ',')) row.push_back(trim(col));

        if (row.size() >= 12) {
            ConfigData cd;
            cd.data_type = row[0];
            cd.bin_file_path = row[4];
            cd.noise_data_path = row[5];
            cd.annealed_bin_path = row[6];
            try {
                cd.ecg_sampling_rate = std::stod(row[11]);
            }
            catch (...) {
                cd.ecg_sampling_rate = 256.0;  // Default
            }
            configs.push_back(cd);
        }
    }
    return configs;
}

std::vector<std::pair<double, double>> readNoiseCSV(const std::string& path) {
    std::vector<std::pair<double, double>> noiseSEG;
    std::ifstream file(path);
    if (!file.is_open()) return noiseSEG;

    std::string line;
    std::getline(file, line); // Skip Header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string col;
        std::vector<std::string> row;
        while (std::getline(ss, col, ',')) row.push_back(trim(col));

        if (row.size() >= 4) {
            try {
                double start = std::stod(row[2]);
                double end = std::stod(row[3]);
                noiseSEG.push_back({ start, end });
            }
            catch (...) {}
        }
    }
    return noiseSEG;
}

// Load structured .bin file matching GUI's setFileSource
BinData readStructuredBin(const std::string& path) {
    BinData data;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open .bin file: " + path);
    }

    const size_t MAX_SAMPLES = 100000000;  // 100 million samples max (adjust for your data)
    const size_t MAX_SLEEP_STAGES = 100000;  // Reasonable limit for sleep stages

    // Read header
    file.read(reinterpret_cast<char*>(&data.ecgSR), sizeof(double));
    file.read(reinterpret_cast<char*>(&data.ppgSR), sizeof(double));
    file.read(reinterpret_cast<char*>(&data.scoring_epoch_sec), sizeof(double));

    uint64_t totalEcgSamples, totalSignal2Samples, totalSignal3Samples, totalPpgSamples, totalSleepSamples;
    file.read(reinterpret_cast<char*>(&totalEcgSamples), sizeof(uint64_t));
    file.read(reinterpret_cast<char*>(&totalSignal2Samples), sizeof(uint64_t));
    file.read(reinterpret_cast<char*>(&totalSignal3Samples), sizeof(uint64_t));
    file.read(reinterpret_cast<char*>(&totalPpgSamples), sizeof(uint64_t));
    file.read(reinterpret_cast<char*>(&totalSleepSamples), sizeof(uint64_t));

    // Validate sizes
    if (totalEcgSamples > MAX_SAMPLES || totalSignal2Samples > MAX_SAMPLES || totalSignal3Samples > MAX_SAMPLES ||
        totalPpgSamples > MAX_SAMPLES || totalSleepSamples > MAX_SAMPLES) {
        throw std::runtime_error("Invalid sample counts in .bin file (too large)");
    }

    // Read data
    data.ecg.resize(totalEcgSamples);
    file.read(reinterpret_cast<char*>(data.ecg.data()), totalEcgSamples * sizeof(double));

    data.ecg2.resize(totalSignal2Samples);
    file.read(reinterpret_cast<char*>(data.ecg2.data()), totalSignal2Samples * sizeof(double));

    data.ecg3.resize(totalSignal3Samples);
    file.read(reinterpret_cast<char*>(data.ecg3.data()), totalSignal3Samples * sizeof(double));

    data.ppg.resize(totalPpgSamples);
    file.read(reinterpret_cast<char*>(data.ppg.data()), totalPpgSamples * sizeof(double));

    data.sleepStages.resize(totalSleepSamples);
    file.read(reinterpret_cast<char*>(data.sleepStages.data()), totalSleepSamples * sizeof(double));

    return data;
}

// ============================================================================
// MATLAB-LIKE ANNEALING LOGIC (Exact Translation)
// ============================================================================

size_t closest_idx(const std::vector<double>& time_seconds, double target) {
    auto it = std::lower_bound(time_seconds.begin(), time_seconds.end(), target);
    if (it == time_seconds.end()) return time_seconds.size() - 1;
    if (it == time_seconds.begin()) return 0;
    if (std::abs(*it - target) < std::abs(*(it - 1) - target)) return it - time_seconds.begin();
    return (it - 1) - time_seconds.begin();
}

std::vector<SegmentIndices> MergeSegments(const std::vector<SegmentIndices>& segments) {
    if (segments.empty()) return {};
    std::vector<SegmentIndices> merged = segments;
    std::sort(merged.begin(), merged.end(), [](const SegmentIndices& a, const SegmentIndices& b) {
        return a.start_idx < b.start_idx;
        });
    std::vector<SegmentIndices> result;
    SegmentIndices current = merged[0];
    for (size_t i = 1; i < merged.size(); ++i) {
        if (merged[i].start_idx <= current.end_idx + 1) {
            current.end_idx = std::max(current.end_idx, merged[i].end_idx);
        }
        else {
            result.push_back(current);
            current = merged[i];
        }
    }
    result.push_back(current);
    return result;
}

size_t RoundToClosestBin(const std::vector<size_t>& bin_breaks, size_t idx) {
    auto it = std::lower_bound(bin_breaks.begin(), bin_breaks.end(), idx);
    if (it == bin_breaks.end()) return bin_breaks.size();
    if (it == bin_breaks.begin()) return 1;
    if (std::abs(static_cast<int>(*it) - static_cast<int>(idx)) < std::abs(static_cast<int>(*(it - 1)) - static_cast<int>(idx))) {
        return it - bin_breaks.begin() + 1;
    }
    return (it - bin_breaks.begin());
}

std::vector<FinalBin> AnnealSegments(
    const std::vector<double>& ppg, double ppgSampleRate,
    const std::vector<double>& ecg, double ecgSampleRate,
    const std::vector<std::pair<double, double>>& noiseSEG,
    double scoring_epoch_size_sec,
    const std::vector<double>& sleepStages_double,
    double targetLength
) {
    // Convert sleepStages_double to int
    std::vector<int> sleepStages(sleepStages_double.begin(), sleepStages_double.end());

    // Constants from MATLAB
    double min_exclusion_bin_size_seconds = 5.0;
    double min_bin_size_mins = targetLength / 2.0;

    // Time vectors
    std::vector<double> ecg_time_seconds(ecg.size());
    for (size_t i = 0; i < ecg.size(); ++i) {
        ecg_time_seconds[i] = static_cast<double>(i) / ecgSampleRate;
    }
    std::vector<double> po_time_seconds(ppg.size());
    for (size_t i = 0; i < ppg.size(); ++i) {
        po_time_seconds[i] = static_cast<double>(i) / ppgSampleRate;
    }
    size_t bin_size_samples = static_cast<size_t>(ppgSampleRate * 60.0 * targetLength);

    // Bin count
    size_t remainder = ppg.size() % bin_size_samples;
    double remainder_mins = static_cast<double>(remainder) / ppgSampleRate / 60.0;
    size_t bin_count;
    if (remainder_mins < min_bin_size_mins) {
        bin_count = ppg.size() / bin_size_samples;
    }
    else {
        bin_count = (ppg.size() + bin_size_samples - 1) / bin_size_samples;
    }

    // Bin breaks
    std::vector<size_t> bin_breaks;
    for (size_t i = bin_size_samples; i < ppg.size(); i += bin_size_samples) {
        bin_breaks.push_back(i);
    }
    if (bin_breaks.size() < bin_count) {
        bin_breaks.push_back(ppg.size() - 1);
    }
    else {
        bin_breaks.back() = ppg.size() - 1;
    }

    // Filter exclusions
    std::vector<std::pair<double, double>> exclusions_seconds;
    for (const auto& n : noiseSEG) {
        if (n.second - n.first >= min_exclusion_bin_size_seconds) {
            exclusions_seconds.push_back(n);
        }
    }

    // Map exclusions to indices
    std::vector<size_t> exclusions_indices;
    for (const auto& e : exclusions_seconds) {
        exclusions_indices.push_back(closest_idx(po_time_seconds, e.first));
        exclusions_indices.push_back(closest_idx(po_time_seconds, e.second));
    }

    // Determine bins for exclusions
    std::vector<size_t> exclusions_bin;
    for (size_t idx : exclusions_indices) {
        exclusions_bin.push_back(RoundToClosestBin(bin_breaks, idx));
    }

    // Reshape exclusions: [begin_idx, end_idx, begin_bin, end_bin]
    std::vector<Exclusion> exclusions;
    for (size_t i = 0; i < exclusions_indices.size(); i += 2) {
        size_t begin_idx = exclusions_indices[i];
        size_t end_idx = exclusions_indices[i + 1];
        size_t begin_bin = exclusions_bin[i];
        size_t end_bin = exclusions_bin[i + 1];
        exclusions.push_back({ begin_idx, end_idx, begin_bin, end_bin });
    }

    // Split exclusions across bins
    for (size_t i = 0; i < exclusions.size(); ++i) {
        if (exclusions[i].begin_bin != exclusions[i].end_bin) {
            size_t temp_bin_end = exclusions[i].end_idx;
            size_t temp_bin_num = exclusions[i].end_bin;
            for (size_t bin = exclusions[i].begin_bin; bin <= exclusions[i].end_bin; ++bin) {
                if (bin == exclusions[i].begin_bin) {
                    exclusions[i].end_idx = bin_breaks[bin - 1];
                    exclusions[i].end_bin = bin;
                }
                else if (bin == temp_bin_num) {
                    exclusions.push_back({ bin_breaks[bin - 2] + 1, temp_bin_end, bin, bin });
                }
                else {
                    exclusions.push_back({ bin_breaks[bin - 2] + 1, bin_breaks[bin - 1], bin, bin });
                }
            }
        }
    }

    // Sort exclusions and drop last column
    std::sort(exclusions.begin(), exclusions.end(), [](const Exclusion& a, const Exclusion& b) {
        return a.begin_idx < b.begin_idx;
        });
    std::vector<std::tuple<size_t, size_t, size_t>> exclusions_trimmed;
    for (const auto& e : exclusions) {
        exclusions_trimmed.push_back(std::make_tuple(e.begin_idx, e.end_idx, e.begin_bin));
    }

    // Unique update bins
    std::vector<size_t> update_bins;
    for (const auto& e : exclusions_trimmed) {
        if (std::find(update_bins.begin(), update_bins.end(), std::get<2>(e)) == update_bins.end()) {
            update_bins.push_back(std::get<2>(e));
        }
    }

    // Good bins
    std::vector<size_t> good_bins;
    for (size_t i = 1; i <= bin_count; ++i) {
        if (std::find(update_bins.begin(), update_bins.end(), i) == update_bins.end()) {
            good_bins.push_back(i);
        }
    }

    // Good sections: [beg, end, move_dir, move_flag]
    std::vector<GoodSection> good_sections;
    for (size_t bin : good_bins) {
        size_t bin_begin = (bin == 1) ? 0 : bin_breaks[bin - 2] + 1;
        size_t bin_end = bin_breaks[bin - 1];
        good_sections.push_back({ bin_begin, bin_end, 0, 0 });
    }

    // For update bins
    for (size_t bin : update_bins) {
        size_t bin_begin = (bin == 1) ? 0 : bin_breaks[bin - 2] + 1;
        size_t bin_end = bin_breaks[bin - 1];
        size_t bin_half = bin_end - bin_size_samples / 2;

        std::vector<std::pair<size_t, size_t>> bin_exclusions;
        for (const auto& e : exclusions_trimmed) {
            if (std::get<2>(e) == bin) {
                bin_exclusions.push_back({ std::get<0>(e), std::get<1>(e) });
            }
        }

        std::vector<size_t> good_flat = { bin_begin };
        for (const auto& e : bin_exclusions) {
            good_flat.push_back(e.first);
            good_flat.push_back(e.second);
        }
        good_flat.push_back(bin_end);

        std::vector<std::pair<size_t, size_t>> good;
        for (size_t i = 0; i < good_flat.size(); i += 2) {
            good.push_back({ good_flat[i], good_flat[i + 1] });
        }

        // Remove zero-length
        good.erase(std::remove_if(good.begin(), good.end(), [](const std::pair<size_t, size_t>& p) {
            return p.second - p.first == 0;
            }), good.end());

        std::vector<double> good_lengths;
        for (const auto& g : good) {
            good_lengths.push_back(static_cast<double>(g.second - g.first) / ppgSampleRate / 60.0);
        }

        std::vector<size_t> time_mask;
        for (double len : good_lengths) {
            time_mask.push_back(len < min_bin_size_mins ? 1 : 0);
        }

        // Movement direction
        std::vector<size_t> movement_dir(good.size());
        for (size_t i = 0; i < good.size(); ++i) {
            double center = static_cast<double>(good[i].first + good[i].second) / 2.0;
            if (center <= bin_half) {
                movement_dir[i] = (bin == 1) ? 2 : 1;
            }
            else {
                movement_dir[i] = (bin == bin_count) ? 1 : 2;
            }
        }
        for (size_t i = 0; i < movement_dir.size(); ++i) {
            if (time_mask[i]) movement_dir[i] = 0;
        }

        std::vector<size_t> move_flag(movement_dir.size(), 1);
        for (size_t i = 0; i < move_flag.size(); ++i) {
            if (time_mask[i]) move_flag[i] = 0;
        }

        // Potential additions
        std::vector<GoodSection> potential_additions;
        for (size_t i = 0; i < good.size(); ++i) {
            potential_additions.push_back({ good[i].first, good[i].second, movement_dir[i], move_flag[i] });
        }

        // Special for first/last bins
        if (bin == 1) {
            potential_additions.erase(std::remove_if(potential_additions.begin(), potential_additions.end(),
                [bin_half](const GoodSection& p) {
                    double center = static_cast<double>(p.beg + p.end) / 2.0;
                    return center <= 0;
                }), potential_additions.end());
        }
        else if (bin == bin_count) {
            potential_additions.erase(std::remove_if(potential_additions.begin(), potential_additions.end(),
                [bin_half](const GoodSection& p) {
                    double center = static_cast<double>(p.beg + p.end) / 2.0;
                    return center >= bin_half * 2;
                }), potential_additions.end());
        }

        good_sections.insert(good_sections.end(), potential_additions.begin(), potential_additions.end());
    }

    // Sort good_sections
    std::sort(good_sections.begin(), good_sections.end(), [](const GoodSection& a, const GoodSection& b) {
        return a.beg < b.beg;
        });

    // Merge adjacent moving segments
    size_t i = 0;
    while (i < good_sections.size() - 1) {
        if (good_sections[i].end == good_sections[i + 1].beg && good_sections[i].move_flag != 0 && good_sections[i + 1].move_flag != 0) {
            double seg1_size = static_cast<double>(good_sections[i].end - good_sections[i].beg) / ppgSampleRate / 60.0;
            double seg2_size = static_cast<double>(good_sections[i + 1].end - good_sections[i + 1].beg) / ppgSampleRate / 60.0;
            if (seg1_size + seg2_size >= min_bin_size_mins) {
                good_sections[i].move_dir = 0;
                good_sections[i].move_flag = 0;
            }
            else {
                size_t idx = (seg1_size > seg2_size) ? 0 : 1;
                good_sections[i].move_dir = good_sections[i + idx].move_dir;
                good_sections[i].move_flag = 1;
            }
            good_sections[i].end = good_sections[i + 1].end;
            good_sections.erase(good_sections.begin() + i + 1);
        }
        else {
            ++i;
        }
    }

    // Build final_bin_idx
    std::vector<std::vector<SegmentIndices>> final_bin_idx;
    size_t current_bin = 0;
    std::vector<SegmentIndices> temp_bin;

    for (const auto& gs : good_sections) {
        if (gs.move_flag) {
            if (gs.move_dir == 1) {  // left
                if (current_bin > 0) {
                    final_bin_idx[current_bin - 1].push_back({ gs.beg, gs.end });
                }
                else {
                    final_bin_idx.push_back({ {gs.beg, gs.end} });
                    ++current_bin;
                }
            }
            else {  // right
                temp_bin.push_back({ gs.beg, gs.end });
            }
        }
        else {
            std::vector<SegmentIndices> new_bin = temp_bin;
            new_bin.push_back({ gs.beg, gs.end });
            final_bin_idx.push_back(new_bin);
            ++current_bin;
            temp_bin.clear();
        }
    }

    // Merge segments in final_bin_idx
    for (auto& bin : final_bin_idx) {
        if (bin.size() > 1) {
            bin = MergeSegments(bin);
        }
    }

    // Adjust overlapping ends
    for (size_t j = 0; j < final_bin_idx.size() - 1; ++j) {
        if (!final_bin_idx[j].empty() && !final_bin_idx[j + 1].empty()) {
            if (final_bin_idx[j].back().end_idx == final_bin_idx[j + 1][0].start_idx) {
                final_bin_idx[j].back().end_idx -= 1;
            }
        }
    }

    // Sleep stage times
    std::vector<double> sleep_stage_times(sleepStages.size());
    for (size_t j = 0; j < sleepStages.size(); ++j) {
        sleep_stage_times[j] = (j + 1) * scoring_epoch_size_sec;
    }

    // Fill data
    std::vector<FinalBin> annealedSegments;
    for (size_t k = 0; k < final_bin_idx.size(); ++k) {
        const auto& bin = final_bin_idx[k];
        FinalBin fb;
        fb.ppgSampleRate = ppgSampleRate;
        fb.ecgSampleRate = ecgSampleRate;
        fb.ppg_bin_indexs = bin;  // PO indices

        // Compute ecg_indices for this bin
        std::vector<SegmentIndices> ecg_indices;
        for (const auto& po_seg : bin) {
            double start_time = static_cast<double>(po_seg.start_idx) / ppgSampleRate;
            double end_time = static_cast<double>(po_seg.end_idx) / ppgSampleRate;
            size_t ecg_start = closest_idx(ecg_time_seconds, start_time);
            size_t ecg_end = closest_idx(ecg_time_seconds, end_time);
            ecg_indices.push_back({ ecg_start, ecg_end });
        }
        fb.ecg_bin_indexs = ecg_indices;

        for (size_t w = 0; w < bin.size(); ++w) {
            const auto& seg = bin[w];
            const auto& ecg_seg = ecg_indices[w];
            double start_time = static_cast<double>(seg.start_idx) / ppgSampleRate;
            double end_time = static_cast<double>(seg.end_idx) / ppgSampleRate;

            // Sleep stages
            for (size_t t = 0; t < sleep_stage_times.size(); ++t) {
                if (sleep_stage_times[t] >= start_time && sleep_stage_times[t] <= end_time) {
                    fb.sleep_stages.push_back(sleepStages[t]);
                }
            }

            // R-peaks (empty since rs removed)
            // fb.r_peaks remains empty

            // PPG data
            if (seg.start_idx < ppg.size() && seg.end_idx < ppg.size()) {
                fb.po_data.insert(fb.po_data.end(), ppg.begin() + seg.start_idx, ppg.begin() + seg.end_idx + 1);
            }

            // ECG data
            if (ecg_seg.start_idx < ecg.size() && ecg_seg.end_idx < ecg.size()) {
                fb.ecg_data.insert(fb.ecg_data.end(), ecg.begin() + ecg_seg.start_idx, ecg.begin() + ecg_seg.end_idx + 1);
            }
        }
        annealedSegments.push_back(fb);
    }

    return annealedSegments;
}

// ============================================================================
// PROCESSING (With Metadata in Output .bin)
// ============================================================================

void processAnnealing(const ConfigData& sel, const std::string& id, const BinData& binData,
    const std::vector<std::pair<double, double>>& noise) {

    // Use ECG signal 1 as the main ECG (adjust if needed)
    auto annealedSegments = AnnealSegments(binData.ppg, binData.ppgSR, binData.ecg, binData.ecgSR,
        noise, binData.scoring_epoch_sec, binData.sleepStages, 15.0);

    // Output to structured .bin (equivalent to MATLAB .mat metadata)
    // Format:
    // uint64_t num_segments
    // For each segment:
    //   uint32_t ppg_bin_indexs_size
    //   For each: uint64_t start_idx, uint64_t end_idx
    //   uint32_t ecg_bin_indexs_size
    //   For each: uint64_t start_idx, uint64_t end_idx
    //   uint32_t po_data_size
    //   double po_data[po_data_size]
    //   uint32_t ecg_data_size
    //   double ecg_data[ecg_data_size]
    //   uint32_t sleep_stages_size
    //   int32_t sleep_stages[sleep_stages_size]
    //   uint32_t r_peaks_size
    //   double r_peaks[r_peaks_size]
    //   double ppgSampleRate
    //   double ecgSampleRate

    std::string outputFile = sel.annealed_bin_path + "/" + id + "_annealed.bin";
    std::ofstream out(outputFile, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Error: Could not create output file: " << outputFile << std::endl;
        return;
    }

    uint64_t num_segments = annealedSegments.size();
    out.write(reinterpret_cast<const char*>(&num_segments), sizeof(uint64_t));

    for (const auto& seg : annealedSegments) {
        // ppg_bin_indexs
        uint32_t ppg_idx_size = seg.ppg_bin_indexs.size();
        out.write(reinterpret_cast<const char*>(&ppg_idx_size), sizeof(uint32_t));
        for (const auto& idx : seg.ppg_bin_indexs) {
            uint64_t start = idx.start_idx;
            uint64_t end = idx.end_idx;
            out.write(reinterpret_cast<const char*>(&start), sizeof(uint64_t));
            out.write(reinterpret_cast<const char*>(&end), sizeof(uint64_t));
        }

        // ecg_bin_indexs
        uint32_t ecg_idx_size = seg.ecg_bin_indexs.size();
        out.write(reinterpret_cast<const char*>(&ecg_idx_size), sizeof(uint32_t));
        for (const auto& idx : seg.ecg_bin_indexs) {
            uint64_t start = idx.start_idx;
            uint64_t end = idx.end_idx;
            out.write(reinterpret_cast<const char*>(&start), sizeof(uint64_t));
            out.write(reinterpret_cast<const char*>(&end), sizeof(uint64_t));
        }

        // po_data
        uint32_t po_size = seg.po_data.size();
        out.write(reinterpret_cast<const char*>(&po_size), sizeof(uint32_t));
        out.write(reinterpret_cast<const char*>(seg.po_data.data()), po_size * sizeof(double));

        // ecg_data
        uint32_t ecg_size = seg.ecg_data.size();
        out.write(reinterpret_cast<const char*>(&ecg_size), sizeof(uint32_t));
        out.write(reinterpret_cast<const char*>(seg.ecg_data.data()), ecg_size * sizeof(double));

        // sleep_stages
        uint32_t sleep_size = seg.sleep_stages.size();
        out.write(reinterpret_cast<const char*>(&sleep_size), sizeof(uint32_t));
        for (int stage : seg.sleep_stages) {
            int32_t s = stage;
            out.write(reinterpret_cast<const char*>(&s), sizeof(int32_t));
        }

        // r_peaks
        uint32_t r_size = seg.r_peaks.size();
        out.write(reinterpret_cast<const char*>(&r_size), sizeof(uint32_t));
        out.write(reinterpret_cast<const char*>(seg.r_peaks.data()), r_size * sizeof(double));

        // Rates
        out.write(reinterpret_cast<const char*>(&seg.ppgSampleRate), sizeof(double));
        out.write(reinterpret_cast<const char*>(&seg.ecgSampleRate), sizeof(double));
    }
    out.close();

    std::cout << "Created " << annealedSegments.size() << " annealed segments with full metadata." << std::endl;
}

// ============================================================================
// MAIN (Loads .bin Files, Loads Noise .csv, Runs Annealing, Outputs Annealed .bin)
// ============================================================================

int main() {
    try {
        auto configs = readConfigCSV("config.csv");

        std::cout << "Select Dataset:\n";
        for (size_t i = 0; i < configs.size(); ++i) {
            std::cout << i + 1 << ". " << configs[i].data_type << "\n";
        }

        int choice;
        std::cin >> choice;
        if (choice < 1 || choice >(int)configs.size()) return 0;
        ConfigData sel = configs[choice - 1];

        // Ensure output directory exists
        if (!fs::exists(sel.annealed_bin_path)) {
            fs::create_directories(sel.annealed_bin_path);
        }

        std::cout << "\nScanning for .bin files in: " << sel.bin_file_path << "...\n";

        bool foundFiles = false;
        for (const auto& entry : fs::directory_iterator(sel.bin_file_path)) {
            if (entry.path().extension() == ".bin") {
                foundFiles = true;
                std::string id = entry.path().stem().string();

                // Load all data from .bin
                BinData binData = readStructuredBin(entry.path().string());

                // Load noise from CSV
                std::string noiseFile = sel.noise_data_path + "/" + id + "_noise_markings.csv";
                auto noiseSEG = readNoiseCSV(noiseFile);

                std::cout << "Processing ID: " << id << std::endl;
                processAnnealing(sel, id, binData, noiseSEG);
            }
        }

        if (!foundFiles) {
            std::cout << "No .bin files found in " << sel.bin_file_path << "\n";
        }

    }
    catch (const std::exception& e) {
        std::cerr << "Fatal Error: " << e.what() << std::endl;
    }

    std::cout << "\nDone. Press Enter to exit...";
    std::cin.ignore(); std::cin.get();
    return 0;
}
