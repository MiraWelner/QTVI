#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <set>

extern "C" {
#include "edflib.h"
}
#include "pugixml.hpp"

namespace fs = std::filesystem;

struct Config {
    std::string dataType, mainExt, sleepExt, inputPath, outputPath;
    std::string ecg1Label, ecg2Label, ecg3Label, ppgLabel;
    double ecgRate, ppgRate;
};

// Helper to parse CSV lines correctly, handling quotes and commas
std::vector<std::string> parseCsvLine(const std::string& line) {
    std::vector<std::string> fields;
    std::string cur;
    bool inQuotes = false;
    for (size_t i = 0; i < line.length(); ++i) {
        char c = line[i];
        if (c == '\"') inQuotes = !inQuotes;
        else if (c == ',' && !inQuotes) { fields.push_back(cur); cur = ""; }
        else cur += c;
    }
    fields.push_back(cur);
    for (auto& f : fields) {
        size_t first = f.find_first_not_of(" \t\r\n\"");
        if (first == std::string::npos) f = "";
        else {
            size_t last = f.find_last_not_of(" \t\r\n\"");
            f = f.substr(first, last - first + 1);
        }
    }
    return fields;
}

// Case-insensitive substring check
bool contains(std::string haystack, std::string needle) {
    if (needle.empty()) return false;
    std::transform(haystack.begin(), haystack.end(), haystack.begin(), ::toupper);
    std::transform(needle.begin(), needle.end(), needle.begin(), ::toupper);
    return haystack.find(needle) != std::string::npos;
}

// Function to write EDF signal data to the binary file
// Added "const edf_param_struct& param" as the 4th argument
void writeEdfSignal(int handle, int idx, long long n, const edf_param_struct& param, std::ofstream& out, uint64_t& sizeOut) {
    if (idx < 0 || n <= 0) {
        double dummy = -1.0;
        out.write((char*)&dummy, 8);
        sizeOut = 1;
        return;
    }

    std::vector<int> digitalBuf(n);
    std::vector<double> physicalBuf(n);

    // 1. Read the raw digital bits. 
    // If edflib still clamps here, we will see it in the next step.
    edfread_digital_samples(handle, idx, (int)n, digitalBuf.data());

    // 2. Extract scaling factors
    double pMax = param.phys_max;
    double pMin = param.phys_min;
    double dMax = (double)param.dig_max;
    double dMin = (double)param.dig_min;

    // This is the "True" scaling factor used by MATLAB
    double scale = (pMax - pMin) / (dMax - dMin);

    for (long long i = 0; i < n; ++i) {
        // 3. THE FIX: Ignore the 'param' limits. 
        // We treat the digitalBuf[i] as a raw 16-bit value.
        // This allows the result to exceed 312.5.
        physicalBuf[i] = ((double)digitalBuf[i] - dMin) * scale + pMin;
    }

    out.write((char*)physicalBuf.data(), n * 8);
    sizeOut = (uint64_t)n;
}

// Function to write Text/DAT signal data to the binary file
void writeTextSignalByLabel(const fs::path& path, std::string label, std::ofstream& out, uint64_t& sizeOut) {
    std::ifstream in(path);
    if (!in || label.empty()) {
        double v = -1.0; out.write((char*)&v, 8); sizeOut = 1; return;
    }

    std::string line; int colIdx = -1; bool headerFound = false;
    while (std::getline(in, line)) {
        if (contains(line, "Index") || contains(line, label)) {
            std::vector<std::string> hdrs = parseCsvLine(line);
            for (int i = 0; i < (int)hdrs.size(); ++i) {
                if (contains(hdrs[i], label)) { colIdx = i; headerFound = true; break; }
            }
            if (headerFound) break;
        }
    }

    if (!headerFound || colIdx == -1) {
        double v = -1.0; out.write((char*)&v, 8); sizeOut = 1; return;
    }

    uint64_t count = 0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = parseCsvLine(line);
        if (colIdx < (int)row.size() && !row[colIdx].empty()) {
            try {
                double val = std::stod(row[colIdx]);
                out.write((char*)&val, 8);
                count++;
            }
            catch (...) {}
        }
    }

    if (count == 0) { double v = -1.0; out.write((char*)&v, 8); sizeOut = 1; }
    else { sizeOut = count; }
}

// LOADS DYNAMIC CONFIG FROM YOUR 12-COLUMN CSV
bool loadConfig(const std::string& filename, int choice, Config& out) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string target = (choice == 1) ? "MESA" : (choice == 2) ? "BITTIUM" : (choice == 3) ? "CHAOS" : "";
    std::string line;
    std::getline(file, line); // Skip the header row

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = parseCsvLine(line);

        // Expecting 12 columns based on your new structure
        if (row.size() < 12) continue;

        std::string rType = row[0];
        std::transform(rType.begin(), rType.end(), rType.begin(), ::toupper);

        if (rType == target) {
            out.dataType = row[0];
            out.mainExt = row[1];
            out.sleepExt = row[2];
            out.inputPath = row[3];
            out.outputPath = row[4];

            // Map labels from columns 6, 7, 8, 9
            out.ecg1Label = row[6]; // EKG (MESA) or ECG_1 (Bittium)
            out.ecg2Label = row[7];
            out.ecg3Label = row[8];
            out.ppgLabel = row[9]; // Pleth (MESA)

            // Map sampling rates from columns 10, 11
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

void processFile(const fs::path& path, const fs::path& xmlPath, const Config& cfg) {
    fs::path outPath = fs::path(cfg.outputPath) / (path.stem().string() + ".bin");
    std::fstream out(outPath, std::ios::binary | std::ios::out | std::ios::in | std::ios::trunc);
    if (!out.is_open()) {
        std::ofstream create(outPath, std::ios::binary);
        create.close();
        out.open(outPath, std::ios::binary | std::ios::out | std::ios::in);
    }

    std::cout << "Processing: " << path.filename().string() << std::endl;

    uint64_t s1 = 0, s2 = 0, s3 = 0, sp = 0, ss = 0;
    out.seekp(64); // Skip space for the 64-byte binary header

    std::string ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);

    if (ext == ".EDF" || ext == ".BDF") {
        edf_hdr_struct hdr;
        if (edfopen_file_readonly(path.string().c_str(), &hdr, EDFLIB_READ_ALL_ANNOTATIONS)) { out.close(); return; }

        int i1 = -1, i2 = -1, i3 = -1, ip = -1;
        std::set<int> used;

        auto find = [&](std::string l) {
            if (l.empty()) return -1;
            for (int i = 0; i < hdr.edfsignals; ++i) {
                if (!used.count(i) && contains(hdr.signalparam[i].label, l)) {
                    used.insert(i);
                    return i;
                }
            }
            return -1;
            };

        i1 = find(cfg.ecg1Label);
        i2 = find(cfg.ecg2Label);
        i3 = find(cfg.ecg3Label);
        ip = find(cfg.ppgLabel);

        std::ofstream wrapper;
        wrapper.basic_ios<char>::rdbuf(out.rdbuf());

        // NOTE: Passing hdr.signalparam[i1] as the NEW 4th argument
        writeEdfSignal(hdr.handle, i1, (i1 >= 0) ? hdr.signalparam[i1].smp_in_file : 0, hdr.signalparam[i1], wrapper, s1);
        writeEdfSignal(hdr.handle, i2, (i2 >= 0) ? hdr.signalparam[i2].smp_in_file : 0, hdr.signalparam[i2], wrapper, s2);
        writeEdfSignal(hdr.handle, i3, (i3 >= 0) ? hdr.signalparam[i3].smp_in_file : 0, hdr.signalparam[i3], wrapper, s3);
        writeEdfSignal(hdr.handle, ip, (ip >= 0) ? hdr.signalparam[ip].smp_in_file : 0, hdr.signalparam[ip], wrapper, sp);

        edfclose_file(hdr.handle);
    }
    else {
        // Handle Text/DAT files (CHAOS)
        std::ofstream wrapper;
        wrapper.basic_ios<char>::rdbuf(out.rdbuf());
        writeTextSignalByLabel(path, cfg.ecg1Label, wrapper, s1);
        writeTextSignalByLabel(path, cfg.ecg2Label, wrapper, s2);
        writeTextSignalByLabel(path, cfg.ecg3Label, wrapper, s3);
        writeTextSignalByLabel(path, cfg.ppgLabel, wrapper, sp);
    }

    // Process Sleep Stages from XML
    std::vector<double> ds;
    if (!cfg.sleepExt.empty() && !xmlPath.empty() && fs::exists(xmlPath)) {
        pugi::xml_document doc;
        if (doc.load_file(xmlPath.string().c_str())) {
            for (auto n : doc.select_nodes("//SleepStage")) {
                double v = n.node().text().as_double();
                ds.push_back(v == 5.0 ? 4.0 : v); // Map stage 5 to 4 if needed
            }
        }
    }
    if (ds.empty()) ds.push_back(-1.0);
    ss = (uint64_t)ds.size();
    out.write((char*)ds.data(), ss * 8);

    // Write the 64-byte Header
    out.seekp(0);
    double er = cfg.ecgRate, pr = cfg.ppgRate, ep = 30.0;
    out.write((char*)&er, 8); // ECG Sample Rate
    out.write((char*)&pr, 8); // PPG Sample Rate
    out.write((char*)&ep, 8); // Epoch Length
    out.write((char*)&s1, 8); // Count ECG1
    out.write((char*)&s2, 8); // Count ECG2
    out.write((char*)&s3, 8); // Count ECG3
    out.write((char*)&sp, 8); // Count PPG
    out.write((char*)&ss, 8); // Count Sleep Stages
    out.close();
}

int main(int argc, char* argv[]) {
    std::cout << "Select Dataset:\n1: MESA\n2: Bittium\n3: CHAOS\nChoice: ";
    int choice;
    if (!(std::cin >> choice)) return 1;

    Config cfg;
    if (!loadConfig("config.csv", choice, cfg)) {
        std::cerr << "Error: Could not find configuration for selection " << choice << " in config.csv" << std::endl;
        return 1;
    }

    fs::create_directories(cfg.outputPath);
    std::string tExt = cfg.mainExt;
    std::transform(tExt.begin(), tExt.end(), tExt.begin(), ::toupper);

    for (const auto& entry : fs::recursive_directory_iterator(cfg.inputPath)) {
        if (!entry.is_regular_file()) continue;

        std::string fExt = entry.path().extension().string();
        std::transform(fExt.begin(), fExt.end(), fExt.begin(), ::toupper);

        if (fExt == tExt) {
            fs::path xml;
            if (!cfg.sleepExt.empty()) {
                std::string stem = entry.path().stem().string();
                for (const auto& f : fs::directory_iterator(entry.path().parent_path())) {
                    std::string cExt = f.path().extension().string();
                    std::transform(cExt.begin(), cExt.end(), cExt.begin(), ::toupper);
                    if (cExt == cfg.sleepExt && f.path().stem().string().find(stem) != std::string::npos) {
                        xml = f.path();
                        break;
                    }
                }
            }
            processFile(entry.path(), xml, cfg);
        }
    }
    std::cout << "Processing Complete." << std::endl;
    return 0;
}
