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

// MATLAB's closest_idx logic
uint64_t closest_idx(double target_time, double sr) {
    // MATLAB: closest_idx(time_vector, target) -> round(target * SR) + 1
    return static_cast<uint64_t>(std::round(target_time * sr)) + 1;
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

std::vector<std::pair<uint64_t, uint64_t>> readNoiseCSV(const std::string& path) {
    std::vector<std::pair<uint64_t, uint64_t>> noise;
    std::ifstream file(path);
    if (!file.is_open()) return noise;
    std::string line; std::getline(file, line);
    while (std::getline(file, line)) {
        std::stringstream ss(line); std::string col; std::vector<std::string> cols;
        while (std::getline(ss, col, ',')) cols.push_back(cleanField(col));
        if (cols.size() >= 2) {
            noise.push_back({ std::stoull(cols[0]), std::stoull(cols[1]) });
        }
    }
    return noise;
}

// ============================================================================
// CORE ANNEALING (EXACT MATLAB REPLICA)
// ============================================================================

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

struct Exclusion { uint64_t start, end; int bin; };
struct GoodSection { uint64_t start, end; int dir, flag; };

std::vector<FinalSegment> AnnealSegments(const RawData& data, const std::vector<std::pair<uint64_t, uint64_t>>& noise, double targetLenMins) {
    const double ppgSR = data.ppgSR;
    const uint64_t bin_size = static_cast<uint64_t>(ppgSR * 60.0 * targetLenMins);
    const double min_bin_size_mins = targetLenMins / 2.0;
    const uint64_t total_len = data.ppg.size();

    // 1. Ideal Bin Breaks (MATLAB logic)
    int bin_count;
    double rem_mins = (double)(total_len % bin_size) / ppgSR / 60.0;
    if (rem_mins < min_bin_size_mins) bin_count = (int)std::floor(total_len / (double)bin_size);
    else bin_count = (int)std::ceil(total_len / (double)bin_size);

    std::vector<uint64_t> bin_breaks;
    for (uint64_t b = bin_size + 1; b <= total_len; b += bin_size) bin_breaks.push_back(b);
    if (bin_breaks.size() < (size_t)bin_count) bin_breaks.push_back(total_len);
    else bin_breaks.back() = total_len;

    // 2. Map Exclusions (1-based)
    std::vector<Exclusion> exclusions;
    for (auto n : noise) {
        if ((double)(n.second - n.first) / ppgSR < 5.0) continue;
        int b1 = -1, b2 = -1;
        for (int i = 0; i < (int)bin_breaks.size(); ++i) {
            if (b1 == -1 && n.first <= bin_breaks[i]) b1 = i + 1;
            if (b2 == -1 && n.second <= bin_breaks[i]) b2 = i + 1;
        }
        if (b1 == -1) b1 = bin_count; if (b2 == -1) b2 = bin_count;

        for (int b = b1; b <= b2; ++b) {
            uint64_t b_end = bin_breaks[b - 1];
            uint64_t b_start = (b == 1) ? 1 : bin_breaks[b - 2];
            exclusions.push_back({ std::max(n.first, b_start), std::min(n.second, b_end), b });
        }
    }

    std::vector<int> update_bins;
    for (auto& ex : exclusions) if (std::find(update_bins.begin(), update_bins.end(), ex.bin) == update_bins.end()) update_bins.push_back(ex.bin);
    // After collecting update_bins from exclusions
    std::sort(update_bins.begin(), update_bins.end());


    // 3. Extract Good Sections
    std::vector<GoodSection> good_sections;
    for (int b = 1; b <= bin_count; ++b) {
        uint64_t b_end = bin_breaks[b - 1];
        uint64_t b_start = b_end - bin_size; if (b == 1) b_start = 1;

        if (std::find(update_bins.begin(), update_bins.end(), b) == update_bins.end()) {
            good_sections.push_back({ b_start, b_end, 0, 0 });
        }
        else {
            uint64_t b_half = b_end - (bin_size / 2);
            std::vector<std::pair<uint64_t, uint64_t>> bin_ex;
            for (auto& ex : exclusions) if (ex.bin == b) bin_ex.push_back({ ex.start, ex.end });
            bin_ex = MergeSegments(bin_ex);

            uint64_t cur = b_start;
            auto add_seg = [&](uint64_t s, uint64_t e) {
                if (e <= s) return;
                int dir = 0, flag = 0;
                if ((double)(e - s) / ppgSR / 60.0 < min_bin_size_mins) {
                    flag = 1;
                    if (b == 1) { if (e > b_half) dir = 2; }
                    else if (b == bin_count) { if (s < b_half) dir = 1; }
                    else { dir = ((s + e) / 2 < b_half) ? 1 : 2; }
                }
                if (flag == 0 || dir != 0) good_sections.push_back({ s, e, dir, flag });
                };
            for (auto& ex : bin_ex) { add_seg(cur, ex.first); cur = ex.second; }
            add_seg(cur, b_end);
        }
    }
    std::sort(good_sections.begin(), good_sections.end(), [](const GoodSection& a, const GoodSection& b) { return a.start < b.start; });

    // 4. Merge Moving
    for (size_t i = 0; i < good_sections.size() && good_sections.size() > 1 && i < good_sections.size() - 1; ) {
        if (good_sections[i].end == good_sections[i + 1].start && good_sections[i].flag && good_sections[i + 1].flag) {
            double s1 = (double)(good_sections[i].end - good_sections[i].start) / ppgSR / 60.0;
            double s2 = (double)(good_sections[i + 1].end - good_sections[i + 1].start) / ppgSR / 60.0;
            if (s1 + s2 >= min_bin_size_mins) { good_sections[i].dir = 0; good_sections[i].flag = 0; }
            else { good_sections[i].dir = (s1 > s2) ? good_sections[i].dir : good_sections[i + 1].dir; }
            good_sections[i].end = good_sections[i + 1].end;
            good_sections.erase(good_sections.begin() + i + 1);
        }
        else i++;
    }

    // 5. Assignment
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> final_bins_idx(bin_count);
    std::vector<std::pair<uint64_t, uint64_t>> temp_moving;
    int cur_bin = 0;
    int cur_b_idx = 0;
    for (auto& gs : good_sections) {
        if (gs.flag) {
            if (gs.dir == 1) {
                // Move Left
                int target = (cur_b_idx > 0) ? cur_b_idx - 1 : cur_b_idx;
                final_bins_idx[target].push_back({ gs.start, gs.end });
            }
            else {
                // Move Right
                temp_moving.push_back({ gs.start, gs.end });
            }
        }
        else {
            // Stay: Append any segments waiting to move right, then add current
            for (auto& t : temp_moving) final_bins_idx[cur_b_idx].push_back(t);
            temp_moving.clear();
            final_bins_idx[cur_b_idx].push_back({ gs.start, gs.end });
            cur_b_idx++;
        }
    }

    for (auto& b : final_bins_idx) b = MergeSegments(b);
    for (int i = 0; i < bin_count - 1; ++i) {
        if (!final_bins_idx[i].empty() && !final_bins_idx[i + 1].empty())
            if (final_bins_idx[i].back().second == final_bins_idx[i + 1].front().first) final_bins_idx[i].back().second--;
    }

    // 6. Data Collection (Exact MATLAB condition)
    std::vector<FinalSegment> results(bin_count);
    for (int i = 0; i < bin_count; ++i) {
        for (auto& seg : final_bins_idx[i]) {
            // PPG
            for (uint64_t k = seg.first; k <= seg.second && k <= total_len; ++k)
                results[i].ppg.push_back(data.ppg[k - 1]);

            // ECG
            uint64_t e_s = closest_idx((double)(seg.first - 1) / ppgSR, data.ecgSR);
            uint64_t e_e = closest_idx((double)(seg.second - 1) / ppgSR, data.ecgSR);
            for (uint64_t k = e_s; k <= e_e && k <= data.ecg.size(); ++k)
                results[i].ecg.push_back(data.ecg[k - 1]);

            // Sleep Extraction
            double t1 = (double)(seg.first - 1) / ppgSR;
            double tend = (double)(seg.second - 1) / ppgSR;

            for (size_t j = 0; j < data.sleepStages.size(); ++j) {
                double stime = (double)(j + 1) * data.scoringEpochSec;
                if (stime >= t1 && stime <= tend) {
                    results[i].sleep_stages.push_back(data.sleepStages[j]);
                }
            }


        }
    }
    return results;
}

int main() {
    std::vector<ProjectConfig> projects;
    std::ifstream cfg("config.csv");
    if (!cfg.is_open()) return 1;
    std::string line; std::getline(cfg, line);
    while (std::getline(cfg, line)) {
        std::stringstream ss(line); std::string col; std::vector<std::string> cols;
        while (std::getline(ss, col, ',')) cols.push_back(cleanField(col));
        if (cols.size() >= 13) projects.push_back({ cols[0], cols[4], cols[5], cols[6], std::stod(cols[11]), std::stod(cols[12]) });
    }
    std::cout << "Select Dataset:\n";
    for (size_t i = 0; i < projects.size(); ++i) std::cout << i + 1 << ". " << projects[i].dataType << "\n";
    int choice; std::cin >> choice;
    ProjectConfig sel = projects[choice - 1];
    if (!fs::exists(sel.annealedPath)) fs::create_directories(sel.annealedPath);
    for (const auto& entry : fs::directory_iterator(sel.binPath)) {
        if (entry.path().extension() == ".bin") {
            std::string id = entry.path().stem().string();
            try {
                RawData raw = readStructuredBin(entry.path().string());
                auto results = AnnealSegments(raw, readNoiseCSV(sel.noisePath + "/" + id + "_noise_markings.csv"), 1.0);
                std::ofstream out(sel.annealedPath + "/" + id + "_annealed.bin", std::ios::binary);
                uint64_t nB = results.size(); out.write((char*)&nB, 8);
                for (auto& s : results) {
                    // Headers should be 8 bytes (uint64_t) to match the reader
                    uint64_t pS = s.ppg.size();
                    uint64_t eS = s.ecg.size();
                    uint64_t sS = s.sleep_stages.size();

                    // Write PPG (8-byte size header + 8-byte doubles)
                    out.write((char*)&pS, 8);
                    out.write((char*)s.ppg.data(), pS * 8);

                    // Write ECG (8-byte size header + 8-byte doubles)
                    out.write((char*)&eS, 8);
                    out.write((char*)s.ecg.data(), eS * 8);

                    // Write Sleep (8-byte size header + 8-byte doubles)
                    out.write((char*)&sS, 8);
                    out.write((char*)s.sleep_stages.data(), sS * 8);
                }

                std::cout << "Processed: " << id << std::endl;
            }
            catch (...) {}
        }
    }
    return 0;
}

