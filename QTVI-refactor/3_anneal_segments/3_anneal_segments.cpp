/**
 * @file   3_anneal_segments.cpp
 * @brief  Takes a raw data .bin (88-byte header) and a noise markings .bin,
 *         outputs 1-minute segments with marked noise removed.
 *
 *         Noise exclusion logic:
 *           - ECG noise is only excluded if ALL 3 ECG channels have
 *             overlapping noise marks for a given time region.
 *           - PPG noise is excluded independently.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-22
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include "AnnealSegments.hpp"

namespace fs = std::filesystem;
double bin_length = 1.0;

// ============================================================================
// Binary I/O
// ============================================================================

/**
 * @brief Read a noise markings .bin produced by step 2 (noise_marking_gui).
 *
 * Format: [uint64 count] then count × [6 doubles per row]
 *   Row: startSample, endSample, startSec, endSec, labelId, typeId
 *   labelId: 0=unknown, 1=PPG, 2=ECG1, 3=ECG2, 4=ECG3
 */
static NoiseMarkings read_noise_bin(const std::string& path) {
    NoiseMarkings m;
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) return m;

    uint64_t count = 0;
    f.read(reinterpret_cast<char*>(&count), 8);

    for (uint64_t i = 0; i < count; ++i) {
        double row[6];
        f.read(reinterpret_cast<char*>(row), 48);

        std::pair<double, double> iv = { row[2], row[3] };
        switch (static_cast<int>(row[4])) {
        case 1: m.ppg.push_back(iv);  break;
        case 2: m.ecg1.push_back(iv); break;
        case 3: m.ecg2.push_back(iv); break;
        case 4: m.ecg3.push_back(iv); break;
        default: break;
        }
    }
    return m;
}

/**
 * @brief Read the 88-byte header data .bin from file_to_bin (step 1).
 *
 * 88-byte header:
 *   [double] ecgSR, ppgSR, scoringEpochSec
 *   [uint64] nEcg1, nEcg2, nEcg3, nPpg, nSleep, nAbs1, nAbs2, nAbs3
 *
 * Signal order: ECG1, ECG2, ECG3, PPG, Sleep, |ECG1|, |ECG2|, |ECG3|
 * Absolute-value channels are skipped.
 */
static RawData read_data_bin(const std::string& path) {
    RawData d;
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("Cannot open bin: " + path);

    f.seekg(0, std::ios::end);
    const uint64_t fileSize = static_cast<uint64_t>(f.tellg());
    f.seekg(0, std::ios::beg);

    double ecgSR, ppgSR, epochSec;
    uint64_t nEcg1, nEcg2, nEcg3, nPpg, nSleep;

    f.read(reinterpret_cast<char*>(&ecgSR), 8);
    f.read(reinterpret_cast<char*>(&ppgSR), 8);
    f.read(reinterpret_cast<char*>(&epochSec), 8);
    f.read(reinterpret_cast<char*>(&nEcg1), 8);
    f.read(reinterpret_cast<char*>(&nEcg2), 8);
    f.read(reinterpret_cast<char*>(&nEcg3), 8);
    f.read(reinterpret_cast<char*>(&nPpg), 8);
    f.read(reinterpret_cast<char*>(&nSleep), 8);

    // Skip nAbs1, nAbs2, nAbs3 (3 × uint64)
    f.seekg(24, std::ios::cur);

    d.ecgSR = ecgSR;
    d.ppgSR = ppgSR;
    d.scoringEpochSec = epochSec;

    auto safeRead = [&](std::vector<double>& dest, uint64_t count) {
        uint64_t pos = static_cast<uint64_t>(f.tellg());
        uint64_t remaining = (fileSize > pos) ? (fileSize - pos) / 8 : 0;
        uint64_t actual = std::min(count, remaining);
        dest.resize(actual);
        if (actual > 0)
            f.read(reinterpret_cast<char*>(dest.data()), actual * 8);
        };

    safeRead(d.ecg1, nEcg1);
    safeRead(d.ecg2, nEcg2);
    safeRead(d.ecg3, nEcg3);
    safeRead(d.ppg, nPpg);
    safeRead(d.sleepStages, nSleep);

    return d;
}

/**
 * @brief Write the annealed output .bin.
 *
 * Header: [uint64 nSegments] [double ppgSR] [double ecgSR] [double epochSec]
 * Per segment: ppg_bin_indexs, ecg_bin_indexs, ppg, ecg1, ecg2, ecg3, sleep
 */
static void write_output_bin(const std::string& path,
    const std::vector<FinalSegment>& segs)
{
    std::ofstream out(path, std::ios::binary);
    uint64_t n = segs.size();
    out.write(reinterpret_cast<char*>(&n), 8);

    if (n > 0) {
        out.write(reinterpret_cast<const char*>(&segs[0].ppgSampleRate), 8);
        out.write(reinterpret_cast<const char*>(&segs[0].ecgSampleRate), 8);
        out.write(reinterpret_cast<const char*>(&segs[0].scoring_epoch_size_sec), 8);
    }

    auto writePairs = [&](const std::vector<std::pair<uint64_t, uint64_t>>& v) {
        uint64_t sz = v.size();
        out.write(reinterpret_cast<char*>(&sz), 8);
        for (const auto& p : v) {
            out.write(reinterpret_cast<const char*>(&p.first), 8);
            out.write(reinterpret_cast<const char*>(&p.second), 8);
        }
        };

    auto writeVec = [&](const std::vector<double>& v) {
        uint64_t sz = v.size();
        out.write(reinterpret_cast<char*>(&sz), 8);
        if (!v.empty())
            out.write(reinterpret_cast<const char*>(v.data()), sz * 8);
        };

    for (const auto& s : segs) {
        writePairs(s.ppg_bin_indexs);
        writePairs(s.ecg_bin_indexs);
        writeVec(s.ppg);
        writeVec(s.ecg1);
        writeVec(s.ecg2);
        writeVec(s.ecg3);
        writeVec(s.sleep_stages);
    }
}

// ============================================================================
// Config parsing
// ============================================================================

static std::vector<ProjectConfig> readConfig(const std::string& path) {
    std::vector<ProjectConfig> projects;
    std::ifstream cfg(path);
    if (!cfg.is_open()) {
        std::cerr << "Could not open " << path << std::endl;
        return projects;
    }

    std::string line;
    std::getline(cfg, line); // skip header
    while (std::getline(cfg, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string col;
        std::vector<std::string> cols;
        while (std::getline(ss, col, ',')) cols.push_back(col);

        ProjectConfig pc;
        pc.dataType = cols.size() > 0 ? cols[0] : "unknown";
        pc.binPath = cols.size() > 4 ? cols[4] : "";
        pc.noisePath = cols.size() > 5 ? cols[5] : "";
        pc.annealedPath = cols.size() > 6 ? cols[6] : "";
        projects.push_back(pc);
    }
    return projects;
}

// ============================================================================
// Main
// ============================================================================

int main() {
    auto projects = readConfig("config.csv");
    if (projects.empty()) return 1;

    std::cout << "Select Dataset:\n";
    for (size_t i = 0; i < projects.size(); ++i)
        std::cout << i + 1 << ". " << projects[i].dataType << "\n";

    int choice;
    std::cin >> choice;
    if (choice < 1 || choice >(int)projects.size()) return 1;
    const auto& sel = projects[choice - 1];

    if (!fs::exists(sel.annealedPath))
        fs::create_directories(sel.annealedPath);

    int processed = 0;
    for (const auto& entry : fs::recursive_directory_iterator(sel.binPath)) {
        if (!entry.is_regular_file()) continue;

        std::string ext = entry.path().extension().string();
        std::string name = entry.path().filename().string();
        if (name.find("_noise_markings") != std::string::npos) continue;
        if (ext != ".bin" && ext != ".BIN") continue;

        std::string id = entry.path().stem().string();
        try {
            RawData raw = read_data_bin(entry.path().string());

            NoiseMarkings noise;
            std::string npath = sel.noisePath + "/" + id + "_noise_markings.bin";
            if (fs::exists(npath))
                noise = read_noise_bin(npath);

            auto results = AnnealSegments(raw, noise, bin_length);
            write_output_bin(sel.annealedPath + "/" + id + ".bin", results);

            ++processed;
            std::cerr << "  [" << processed << "] " << id
                << " -> " << results.size() << " bins\n";
        }
        catch (const std::exception& e) {
            std::cerr << "Error processing " << id << ": " << e.what() << "\n";
        }
    }

    std::cerr << "Done. Processed " << processed << " files.\n";
    return 0;
}