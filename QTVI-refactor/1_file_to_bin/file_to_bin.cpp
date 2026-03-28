/**
 * @file   file_to_bin.cpp
 * @brief  Take in a MESA, Bittium, or CHAOS file and converts it to a .bin file of uniform format
 *
 * The format is the following:
 *  * 64-byte header:
 *   Offset  0: ecgRate       (double)  — ECG sampling rate in Hz - this will ALWAYS be 2000.0 because it is resampled
 *   Offset  8: ppgRate       (double)  — PPG sampling rate in Hz this will ALWAYS be 2000.0 because it is resampled
 *   Offset 16: epochSize     (double)  — sleep stage epoch duration in seconds - this is 30 seconds in MESA which is the only format (so far) with sleep
 *   Offset 24: size1         (uint64)  — number of samples in ECG channel 1
 *   Offset 32: size2         (uint64)  — number of samples in ECG channel 2
 *   Offset 40: size3         (uint64)  — number of samples in ECG channel 3
 *   Offset 48: sizeP         (uint64)  — number of samples in PPG channel
 *   Offset 56: sizeS         (uint64)  — number of sleep stage values
 *
 * Signal data (contiguous doubles, immediately after header):
 *   [size1 doubles]  ECG channel 1
 *   [size2 doubles]  ECG channel 2
 *   [size3 doubles]  ECG channel 3
 *   [sizeP doubles]  PPG
 *   [sizeS doubles]  Sleep stages (one value per epoch: 0=Wake, 1=N1, 2=N2, 3=N3, 4=REM, -1=unknown)
 *
 * Missing channels are stored as a single -1.0 with their size field set to 1.
 *
 * @author Mira Welner
 * @email MEW386@pitt.edu
 * @date   2026-03-18
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <set>
#include <cmath>

extern "C" {
#include "edflib.h"
}
#include "pugixml.hpp"
#include "resample.hpp"

static const std::streamoff HEADER_SIZE = 64;
static const double SLEEP_STATE_LENGTH = 30.0;
static const std::string CONFIG_PATH = "config.csv";
static const double final_sampling_rate = 2000.0;

struct config_csv_data {
    std::string dataType, mainExt, sleepExt, inputPath, outputPath;
    std::string ecg1Label, ecg2Label, ecg3Label, ppgLabel;
    double ecgRate, ppgRate;
};

std::vector<std::string> parse_csv_row(const std::string& line) {
    std::vector<std::string> fields;
    std::string cur;
    for (size_t i = 0; i < line.length(); ++i) {
        if (line[i] == ',') { fields.push_back(cur); cur = ""; }
        else cur += line[i];
    }
    fields.push_back(cur);
    for (auto& f : fields) {
        size_t first = f.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) f = "";
        else {
            size_t last = f.find_last_not_of(" \t\r\n");
            f = f.substr(first, last - first + 1);
        }
    }
    return fields;
}

bool contains(std::string search_string, std::string substring) {
    if (substring.empty()) return false;
    std::transform(search_string.begin(), search_string.end(), search_string.begin(), ::toupper);
    std::transform(substring.begin(), substring.end(), substring.begin(), ::toupper);
    return search_string.find(substring) != std::string::npos;
}

void edf_to_bin(int handle, int idx, long long n, double old_rate, std::ofstream& out, uint64_t& sizeOut) {
    if (idx < 0 || n <= 0) {
        double dummy = -1.0;
        out.write((char*)&dummy, 8);
        sizeOut = 1;
        return;
    }

    std::vector<double> buf(n);
    edfread_physical_samples(handle, idx, (int)n, buf.data());

    buf = upsample(buf, old_rate);
    out.write((char*)buf.data(), buf.size() * 8);
    sizeOut = (uint64_t)buf.size();
}

static void dat_to_bin(const std::filesystem::path& path, std::string label, double old_rate, std::ofstream& out, uint64_t& sizeOut) {
    std::ifstream in(path);
    if (!in || label.empty()) {
        double v = -1.0; out.write((char*)&v, 8); sizeOut = 1; return;
    }

    std::string line; int colIdx = -1; bool headerFound = false;
    while (std::getline(in, line)) {
        if (contains(line, "Index") || contains(line, label)) {
            std::vector<std::string> hdrs = parse_csv_row(line);
            for (int i = 0; i < (int)hdrs.size(); ++i) {
                if (contains(hdrs[i], label))
                {
                    colIdx = i;
                    headerFound = true;
                    break;
                }
            }
            if (headerFound) break;
        }
    }

    if (!headerFound || colIdx == -1) {
        double v = -1.0; out.write((char*)&v, 8); sizeOut = 1; return;
    }

    uint64_t count = 0;
    std::vector<double> samples;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = parse_csv_row(line);
        if (colIdx < (int)row.size() && !row[colIdx].empty()) {
            try { samples.push_back(std::stod(row[colIdx])); }
            catch (...) {}
        }
    }
    if (samples.empty()) { double v = -1.0; out.write((char*)&v, 8); sizeOut = 1; }
    else {
        samples = upsample(samples, old_rate);
        out.write((char*)samples.data(), samples.size() * 8);
        sizeOut = (uint64_t)samples.size();
    }
}

static bool load_config(int data_type, config_csv_data& out) {
    std::ifstream file(CONFIG_PATH);
    if (!file.is_open()) return false;

    std::string target = (data_type == 1) ? "MESA" : (data_type == 2) ? "BITTIUM" : (data_type == 3) ? "CHAOS" : "";
    std::string line;
    std::getline(file, line);

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = parse_csv_row(line);

        if (row.size() < 12) continue;

        std::string rType = row[0];
        std::transform(rType.begin(), rType.end(), rType.begin(), ::toupper);

        if (rType == target) {
            out.dataType = row[0];
            out.mainExt = row[1];
            out.sleepExt = row[2];
            out.inputPath = row[3];
            out.outputPath = row[4];

            out.ecg1Label = row[6];
            out.ecg2Label = row[7];
            out.ecg3Label = row[8];
            out.ppgLabel = row[9];

            try {
                out.ecgRate = (!row[10].empty()) ? std::stod(row[10]) : 0.0;
                out.ppgRate = (!row[11].empty()) ? std::stod(row[11]) : 0.0;
            }
            catch (...) {
                out.ecgRate = 256.0;
                out.ppgRate = 256.0;
            }
            return true;
        }
    }
    return false;
}

static void make_binfile(const std::filesystem::path& path, const std::filesystem::path& xmlPath, const config_csv_data& cfg) {
    std::filesystem::path outPath = std::filesystem::path(cfg.outputPath) / (path.stem().string() + ".bin");
    std::ofstream out(outPath, std::ios::binary | std::ios::trunc);
    if (!out.is_open()) {
        std::cerr << "ERROR: cannot create " << outPath << "\n";
        return;
    }

    std::cout << "Processing: " << path.filename().string() << std::endl;

    char zeroes[HEADER_SIZE] = {};
    out.write(zeroes, HEADER_SIZE);

    uint64_t s1 = 0, s2 = 0, s3 = 0, sp = 0;

    std::string ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);

    if (ext == ".EDF") {
        auto hdr = std::make_unique<edf_hdr_struct>();
        if (edfopen_file_readonly(path.string().c_str(), hdr.get(), EDFLIB_READ_ALL_ANNOTATIONS)) {
            std::cerr << "ERROR: cannot open EDF " << path << "\n";
            out.close();
            std::filesystem::remove(outPath);
            return;
        }

        std::set<int> used;
        auto findChannel = [&](const std::string& label) -> int {
            if (label.empty()) return -1;
            for (int i = 0; i < hdr->edfsignals; ++i) {
                if (!used.count(i) && contains(hdr->signalparam[i].label, label)) {
                    used.insert(i);
                    return i;
                }
            }
            return -1;
            };

        int i1 = findChannel(cfg.ecg1Label);
        int i2 = findChannel(cfg.ecg2Label);
        int i3 = findChannel(cfg.ecg3Label);
        int ip = findChannel(cfg.ppgLabel);

        edf_to_bin(hdr->handle, i1, (i1 >= 0) ? hdr->signalparam[i1].smp_in_file : 0, cfg.ecgRate, out, s1);
        edf_to_bin(hdr->handle, i2, (i2 >= 0) ? hdr->signalparam[i2].smp_in_file : 0, cfg.ecgRate, out, s2);
        edf_to_bin(hdr->handle, i3, (i3 >= 0) ? hdr->signalparam[i3].smp_in_file : 0, cfg.ecgRate, out, s3);
        edf_to_bin(hdr->handle, ip, (ip >= 0) ? hdr->signalparam[ip].smp_in_file : 0, cfg.ppgRate, out, sp);

        edfclose_file(hdr->handle);
    }
    else {
        dat_to_bin(path, cfg.ecg1Label, cfg.ecgRate, out, s1);
        dat_to_bin(path, cfg.ecg2Label, cfg.ecgRate, out, s2);
        dat_to_bin(path, cfg.ecg3Label, cfg.ecgRate, out, s3);
        dat_to_bin(path, cfg.ppgLabel, cfg.ppgRate, out, sp);
    }

    std::vector<double> stages;

    if (!cfg.sleepExt.empty() && !xmlPath.empty() && std::filesystem::exists(xmlPath)) {
        pugi::xml_document doc;
        if (doc.load_file(xmlPath.string().c_str())) {
            for (auto node : doc.select_nodes("//SleepStage")) {
                double v = node.node().text().as_double();
                stages.push_back(v == 5.0 ? 4.0 : v);
            }
        }
    }
    if (stages.empty()) {
        stages.push_back(-1.0);
    }

    uint64_t ss = static_cast<uint64_t>(stages.size());
    out.write(reinterpret_cast<const char*>(stages.data()), ss * sizeof(double));

    // Write header
    out.seekp(0);

    auto writeField = [&](const auto& val) {
        out.write(reinterpret_cast<const char*>(&val), sizeof(val));
        };

    writeField(final_sampling_rate);
    writeField(final_sampling_rate);
    writeField(SLEEP_STATE_LENGTH);
    writeField(s1);
    writeField(s2);
    writeField(s3);
    writeField(sp);
    writeField(ss);

    out.close();
}


int main(int argc, char* argv[]) {
    std::cout << "FILE TO BIN\n";
    std::cout << "Select Dataset:\n1: MESA\n2: Bittium\n3: CHAOS\nChoice: ";
    int choice;
    std::cin >> choice;

    config_csv_data cfg;
    if (!load_config(choice, cfg)) {
        std::cerr << "Error: Could not find configuration for selection " << choice << " in config.csv" << std::endl;
        return 1;
    }

    std::filesystem::create_directories(cfg.outputPath);
    std::string tExt = cfg.mainExt;
    std::transform(tExt.begin(), tExt.end(), tExt.begin(), ::toupper);

    for (const auto& entry : std::filesystem::recursive_directory_iterator(cfg.inputPath)) {
        if (!entry.is_regular_file()) continue;

        std::string fExt = entry.path().extension().string();
        std::transform(fExt.begin(), fExt.end(), fExt.begin(), ::toupper);

        if (fExt == tExt) {
            std::filesystem::path xml;
            if (!cfg.sleepExt.empty()) {
                std::string stem = entry.path().stem().string();
                for (const auto& f : std::filesystem::directory_iterator(entry.path().parent_path())) {
                    std::string cExt = f.path().extension().string();
                    std::transform(cExt.begin(), cExt.end(), cExt.begin(), ::toupper);
                    if (cExt == cfg.sleepExt && f.path().stem().string().find(stem) != std::string::npos) {
                        xml = f.path();
                        break;
                    }
                }
            }
            make_binfile(entry.path(), xml, cfg);
        }
    }
    std::cout << "Processing Complete." << std::endl;
    return 0;
}