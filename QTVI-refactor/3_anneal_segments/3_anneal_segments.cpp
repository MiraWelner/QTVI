/**
 * @file   anneal_segments.cpp
 * @brief  Takes in a .bin file representing the data (88-byte header format from file_to_bin),
 *         and a noise markings .bin, and outputs 1-minute segments with marked noise removed.
 *
 *         Noise exclusion logic:
 *           - ECG noise is only excluded if ALL 3 ECG channels have overlapping noise marks
 *             for a given time region. If only 1 or 2 channels are marked, the data is kept.
 *           - PPG noise is excluded independently.
 *
 *         The output .bin file is a vector of the following format:
 *         Header:
 *                The number of segments (uint64)
 *                PPG sample rate (float64)
 *                ECG sample rate (float64)
 *                Scoring epoch size in seconds
 *         Then for each segment in the file you get:
 *                Number of PPG bins (uint64)
 *                Start and end indices (uint64) of each PPG bin
 *                Number of ECG bins (uint64)
 *                Start and end indices (uint64) of each ECG bin
 *                The number of samples in the PPG signal (uint64)
 *                The PPG values (float64)
 *                The number of samples in the ECG1 signal (uint64)
 *                The ECG1 values (float64)
 *                The number of samples in the ECG2 signal (uint64)
 *                The ECG2 values (float64)
 *                The number of samples in the ECG3 signal (uint64)
 *                The ECG3 values (float64)
 *
 * @author Mira Welner
 * @email MEW386@pitt.edu
 * @date   2026-03-22
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include "AnnealSegments.hpp"


 /**
  * @brief Reads a noise markings .bin produced by step 2 (noise_marking_gui).
  *
  * Format: [uint64 count] then count x [6 doubles per row]
  *   Row: startSample, endSample, startSec, endSec, labelId, typeId
  *   labelId: 0=unknown, 1=PPG, 2=ECG1, 3=ECG2, 4=ECG3
  *
  * @return  A NoiseMarkings struct with separate interval lists per channel
  */
NoiseMarkings read_noise_bin(const std::string& path) {
    NoiseMarkings markings;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return markings;

    uint64_t count = 0;
    file.read(reinterpret_cast<char*>(&count), 8);

    for (uint64_t i = 0; i < count; ++i) {
        double row[6];
        file.read(reinterpret_cast<char*>(row), 6 * sizeof(double));

        double startSec = row[2];
        double endSec = row[3];
        int labelId = static_cast<int>(row[4]);

        std::pair<double, double> interval = { startSec, endSec };

        switch (labelId) {
        case 1: markings.ppg.push_back(interval);  break;
        case 2: markings.ecg1.push_back(interval); break;
        case 3: markings.ecg2.push_back(interval); break;
        case 4: markings.ecg3.push_back(interval); break;
        default: break;
        }
    }
    return markings;
}


/**
 * @brief Reads the 88-byte header data .bin from file_to_bin (step 1).
 *
 * 88-byte header:
 *   [double] ecgSR, ppgSR, scoringEpochSec
 *   [uint64] nEcg1, nEcg2, nEcg3, nPpg, nSleep, nAbs1, nAbs2, nAbs3
 *
 * Signal data order: ECG1, ECG2, ECG3, PPG, Sleep, |ECG1|, |ECG2|, |ECG3|
 * Absolute value channels are skipped.
 *
 * @return A RawData object containing the signals
 */
RawData read_data_bin(const std::string& path) {
    RawData data;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Cannot open bin: " + path);

    file.seekg(0, std::ios::end);
    const uint64_t fileSize = static_cast<uint64_t>(file.tellg());
    file.seekg(0, std::ios::beg);

    double ecgSR, ppgSR, epochSec;
    uint64_t nEcg1, nEcg2, nEcg3, nPpg, nSleep;

    file.read((char*)&ecgSR, 8);
    file.read((char*)&ppgSR, 8);
    file.read((char*)&epochSec, 8);
    file.read((char*)&nEcg1, 8);
    file.read((char*)&nEcg2, 8);
    file.read((char*)&nEcg3, 8);
    file.read((char*)&nPpg, 8);
    file.read((char*)&nSleep, 8);

    data.ecgSR = ecgSR;
    data.ppgSR = ppgSR;
    data.scoringEpochSec = epochSec;

    // Safe read: clamp to remaining file bytes
    auto safeRead = [&](std::vector<double>& dest, uint64_t count, const char* name) {
        uint64_t pos = static_cast<uint64_t>(file.tellg());
        uint64_t remaining = (fileSize > pos) ? (fileSize - pos) / 8 : 0;
        uint64_t actual = std::min(count, remaining);
        dest.resize(actual);
        if (actual > 0)
            file.read((char*)dest.data(), actual * 8);
        };

    safeRead(data.ecg1, nEcg1, "ecg1");
    safeRead(data.ecg2, nEcg2, "ecg2");
    safeRead(data.ecg3, nEcg3, "ecg3");
    safeRead(data.ppg, nPpg, "ppg");
    safeRead(data.sleepStages, nSleep, "sleep");

    return data;
}


void write_output_bin(const std::string& path, const std::vector<FinalSegment>& segments) {
    /**
    * @brief  Writes the bin of the format described in the header string
    * @param path       The place to write the output
    * @param segments   A vector representing the segments being written
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
        while (std::getline(ss, col, ',')) cols.push_back(col);

        ProjectConfig pc;
        pc.dataType = cols.size() > 0 ? cols[0] : "unknown";
        pc.binPath = cols.size() > 4 ? cols[4] : "";
        pc.noisePath = cols.size() > 5 ? cols[5] : "";
        pc.annealedPath = cols.size() > 6 ? cols[6] : "";
        projects.push_back(pc);
    }

    std::cout << "Select Dataset:\n";
    for (size_t i = 0; i < projects.size(); ++i)
        std::cout << i + 1 << ". " << projects[i].dataType << "\n";

    int choice;
    std::cin >> choice;
    if (choice < 1 || choice >(int)projects.size()) return 1;
    ProjectConfig sel = projects[choice - 1];

    if (!std::filesystem::exists(sel.annealedPath))
        std::filesystem::create_directories(sel.annealedPath);

    int fileCount = 0;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(sel.binPath)) {
        if (!entry.is_regular_file()) continue;

        std::string ext = entry.path().extension().string();
        std::string name = entry.path().filename().string();

        if (name.find("_noise_markings") != std::string::npos) continue;
        if (ext != ".bin" && ext != ".BIN") continue;

        fileCount++;
        std::string id = entry.path().stem().string();
        try {
            RawData raw = read_data_bin(entry.path().string());

            std::string noise_path = sel.noisePath + "/" + id + "_noise_markings.bin";
            NoiseMarkings noiseMarkings;
            if (std::filesystem::exists(noise_path)) {
                noiseMarkings = read_noise_bin(noise_path);
            }

            auto results = AnnealSegments(raw, noiseMarkings, 1.0);

            write_output_bin(sel.annealedPath + "/" + id + ".bin", results);
        }
        catch (const std::exception& e) {
            std::cerr << "Error processing " << id << ": " << e.what() << std::endl;
        }
    }
    return 0;
}