#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>

#include "AnnealSegments.hpp"

namespace fs = std::filesystem;

// ============================================================================
// HELPERS
// ============================================================================

std::string cleanField(std::string s) {
    s.erase(0, s.find_first_not_of(" \t\r\n\""));
    s.erase(s.find_last_not_of(" \t\r\n\"") + 1);
    return s;
}

// ============================================================================
// NOISE FILE READER
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
        noiseSEG.push_back({ row[2], row[3] });
    }
    return noiseSEG;
}

// ============================================================================
// DATA FILE READER
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

    data.ecg1.resize(nEcg1);
    file.read((char*)data.ecg1.data(), nEcg1 * 8);

    data.ecg2.resize(nEcg2);
    file.read((char*)data.ecg2.data(), nEcg2 * 8);

    data.ecg3.resize(nEcg3);
    file.read((char*)data.ecg3.data(), nEcg3 * 8);

    data.ppg.resize(nPpg);
    file.read((char*)data.ppg.data(), nPpg * 8);
    data.sleepStages.resize(nSleep);
    file.read((char*)data.sleepStages.data(), nSleep * 8);

    std::cout << "  ECG ch1: " << nEcg1 << " ch2: " << nEcg2
        << " ch3: " << nEcg3 << ", PPG: " << nPpg
        << ", Sleep: " << nSleep << std::endl;

    return data;
}


void writeAnnealedBin(const std::string& path, const std::vector<FinalSegment>& segments) {
    /*
        This function outputs a .bin file containing the annealed data. The format is as such:

        Header (written once):
                uint64 — number of segments
                double — PPG sample rate
                double — ECG sample rate
                double — scoring epoch size (seconds)

        Per segment (repeated for each segment):
            uint64 — number of PPG index pairs, followed by that many (uint64 start, uint64 end) pairs
            uint64 — number of ECG index pairs, followed by that many (uint64 start, uint64 end) pairs
            uint64 nPpg, then double[nPpg] — PPG samples
            uint64 nEcg1, then double[nEcg1] — ECG channel 1
            uint64 nEcg2, then double[nEcg2] — ECG channel 2
            uint64 nEcg3, then double[nEcg3] — ECG channel 3
            uint64 nSleep, then double[nSleep] — sleep stage values
    
    */
    std::ofstream out(path, std::ios::binary);
    uint64_t nB = segments.size();
    out.write((char*)&nB, 8);

    if (nB > 0) {
        out.write((char*)&segments[0].ppgSampleRate, 8);
        out.write((char*)&segments[0].ecgSampleRate, 8);
        out.write((char*)&segments[0].scoring_epoch_size_sec, 8);
    }

    for (const auto& s : segments) {
        uint64_t nPairs = s.ppg_bin_indexs.size();
        out.write((char*)&nPairs, 8);
        for (const auto& p : s.ppg_bin_indexs) {
            out.write((char*)&p.first, 8);
            out.write((char*)&p.second, 8);
        }

        uint64_t ePairs = s.ecg_bin_indexs.size();
        out.write((char*)&ePairs, 8);
        for (const auto& p : s.ecg_bin_indexs) {
            out.write((char*)&p.first, 8);
            out.write((char*)&p.second, 8);
        }

        uint64_t pS = s.ppg.size();
        out.write((char*)&pS, 8);
        out.write((char*)s.ppg.data(), pS * 8);

        uint64_t e1S = s.ecg1.size();
        out.write((char*)&e1S, 8);
        out.write((char*)s.ecg1.data(), e1S * 8);

        uint64_t e2S = s.ecg2.size();
        out.write((char*)&e2S, 8);
        out.write((char*)s.ecg2.data(), e2S * 8);

        uint64_t e3S = s.ecg3.size();
        out.write((char*)&e3S, 8);
        out.write((char*)s.ecg3.data(), e3S * 8);

        uint64_t sS = s.sleep_stages.size();
        out.write((char*)&sS, 8);
        out.write((char*)s.sleep_stages.data(), sS * 8);
    }
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
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string col;
        std::vector<std::string> cols;
        while (std::getline(ss, col, ',')) cols.push_back(cleanField(col));

        std::cout << "Config row: " << (cols.empty() ? "(empty)" : cols[0])
            << " | columns: " << cols.size() << std::endl;

        // Attempt to read whatever we can from the row
        ProjectConfig pc;
        pc.dataType = cols.size() > 0 ? cols[0] : "unknown";
        pc.binPath = cols.size() > 4 ? cols[4] : "";
        pc.noisePath = cols.size() > 5 ? cols[5] : "";
        pc.annealedPath = cols.size() > 6 ? cols[6] : "";

        if (pc.binPath.empty()) {
            std::cout << "  Skipping: no binPath" << std::endl;
            continue;
        }
        if (!fs::exists(pc.binPath)) {
            std::cout << "  Warning: binPath does not exist: [" << pc.binPath << "]" << std::endl;
        }

        projects.push_back(pc);
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

    std::cout << "Scanning: [" << sel.binPath << "]" << std::endl;
    int fileCount = 0;
    for (const auto& entry : fs::recursive_directory_iterator(sel.binPath)) {
        if (!entry.is_regular_file()) continue;

        std::string ext = entry.path().extension().string();
        std::string name = entry.path().filename().string();

        // Skip noise marking files
        if (name.find("_noise_markings") != std::string::npos) continue;
        // Skip non-bin files
        if (ext != ".bin" && ext != ".BIN") continue;

        fileCount++;
        std::string id = entry.path().stem().string();
        std::cout << "  Found: " << entry.path() << std::endl;
        try {
            std::cout << "Processing: " << id << std::endl;
            RawData raw = readStructuredBin(entry.path().string());

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

            writeAnnealedBin(sel.annealedPath + "/" + id + "_annealed.bin", results);
            std::cout << "  Output: " << results.size() << " bins" << std::endl;
        }
        catch (const std::exception& e) {
            std::cerr << "Error processing " << id << ": " << e.what() << std::endl;
        }
    }

    if (fileCount == 0) {
        std::cout << "No .bin files found in [" << sel.binPath << "]" << std::endl;
    }
    std::cout << "Done. Processed " << fileCount << " files." << std::endl;
    return 0;
}