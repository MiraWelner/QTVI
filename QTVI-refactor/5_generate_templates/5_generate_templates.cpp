/**
 * @file   5_generate_templates.cpp
 * @brief  Runner for the template generation pipeline.
 *
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstdint>
#include <cmath>
#include <algorithm>

#include "GenerateTemplates.hpp"

 // ============================================================================
 // Read wave_data .bin
 //
 // Must match write_output_binfile layout exactly:
 //   uint64 numBins
 //   per bin:
 //     9 index arrays: ch1(raw,sq,abs), ch2(raw,sq,abs), ch3(raw,sq,abs)
 //     2 PPG index arrays: ppgMaxAmps, ppgMinAmps
 //     4 raw signal arrays: ppg, ecg1, ecg2, ecg3
 //     6 preprocessed signal arrays: ch1(sq,abs), ch2(sq,abs), ch3(sq,abs)
 //     9 noise flag bytes: ch1(raw,sq,abs), ch2(raw,sq,abs), ch3(raw,sq,abs)
 //     pairs: uint64 count + int64[count*2]
 //     ppg_bin_indexs: uint64 count + pair<uint64,uint64>[]
 //     ecg_bin_indexs: uint64 count + pair<uint64,uint64>[]
 // ============================================================================
static std::vector<output_binfile_data> read_wave_data_binfile(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Could not open: " + path);

    constexpr uint64_t MAX_SANE = 50000000;

    uint64_t numBins = 0;
    file.read(reinterpret_cast<char*>(&numBins), 8);
    std::cerr << "Processing " << numBins << " Bins" <<  std::endl;
    if (numBins > 100000) throw std::runtime_error("numBins too large");

    std::vector<output_binfile_data> results;
    results.reserve(numBins);

    auto readIdx = [&](std::vector<std::size_t>& v) {
        uint64_t sz = 0;
        if (!file.read(reinterpret_cast<char*>(&sz), 8)) { v.clear(); return; }
        if (sz > MAX_SANE) throw std::runtime_error("bad idx size: " + std::to_string(sz));
        v.resize(sz);
        if (sz > 0) {
            std::vector<uint64_t> tmp(sz);
            file.read(reinterpret_cast<char*>(tmp.data()), sz * 8);
            for (uint64_t j = 0; j < sz; ++j)
                v[j] = static_cast<std::size_t>(tmp[j] > 0 ? tmp[j] - 1 : 0);
        }
        };

    auto readSignal = [&](std::vector<double>& sig) {
        uint64_t sz = 0;
        if (!file.read(reinterpret_cast<char*>(&sz), 8)) { sig.clear(); return; }
        if (sz > MAX_SANE) throw std::runtime_error("bad signal size: " + std::to_string(sz));
        sig.resize(sz);
        if (sz > 0) file.read(reinterpret_cast<char*>(sig.data()), sz * 8);
        };

    auto readPairVec = [&](std::vector<std::pair<uint64_t, uint64_t>>& v) {
        uint64_t sz = 0;
        if (!file.read(reinterpret_cast<char*>(&sz), 8)) { v.clear(); return; }
        if (sz > MAX_SANE) { v.clear(); return; }
        v.resize(sz);
        if (sz > 0) file.read(reinterpret_cast<char*>(v.data()), sz * 16);
        };

    for (uint64_t i = 0; i < numBins; ++i) {
        results.push_back(output_binfile_data{});
        auto& bin = results.back();

        /* 9 R-peak index arrays: 3 methods x 3 channels */
        readIdx(bin.ch1.raw);
        readIdx(bin.ch1.squared);
        readIdx(bin.ch1.absval);

        readIdx(bin.ch2.raw);
        readIdx(bin.ch2.squared);
        readIdx(bin.ch2.absval);

        readIdx(bin.ch3.raw);
        readIdx(bin.ch3.squared);
        readIdx(bin.ch3.absval);

        /* 2 PPG index arrays */
        readIdx(bin.ppgMaxAmps);
        readIdx(bin.ppgMinAmps);

        /* 4 raw signal arrays */
        readSignal(bin.ppgSignal);
        readSignal(bin.ecgSignal);
        readSignal(bin.ecgSignal2);
        readSignal(bin.ecgSignal3);

        /* 6 preprocessed signal arrays */
        readSignal(bin.ch1.squared_signal);
        readSignal(bin.ch1.absval_signal);
        readSignal(bin.ch2.squared_signal);
        readSignal(bin.ch2.absval_signal);
        readSignal(bin.ch3.squared_signal);
        readSignal(bin.ch3.absval_signal);

        /* 9 noise flags */
        uint8_t flags[9] = {};
        file.read(reinterpret_cast<char*>(flags), 9);
        bin.ch1.raw_noisy = (flags[0] != 0);
        bin.ch1.squared_noisy = (flags[1] != 0);
        bin.ch1.absval_noisy = (flags[2] != 0);
        bin.ch2.raw_noisy = (flags[3] != 0);
        bin.ch2.squared_noisy = (flags[4] != 0);
        bin.ch2.absval_noisy = (flags[5] != 0);
        bin.ch3.raw_noisy = (flags[6] != 0);
        bin.ch3.squared_noisy = (flags[7] != 0);
        bin.ch3.absval_noisy = (flags[8] != 0);

        /* Pairs */
        uint64_t numPairs = 0;
        file.read(reinterpret_cast<char*>(&numPairs), 8);
        if (numPairs > MAX_SANE) throw std::runtime_error("bad pairs count");
        bin.pairs.resize(numPairs);
        if (numPairs > 0) {
            std::vector<int64_t> pairBuf(numPairs * 2);
            file.read(reinterpret_cast<char*>(pairBuf.data()), numPairs * 16);
            for (uint64_t j = 0; j < numPairs; ++j) {
                bin.pairs[j].resize(2);
                bin.pairs[j][0] = (pairBuf[j * 2] == -1)
                    ? -1.0 : static_cast<double>(pairBuf[j * 2] - 1);
                bin.pairs[j][1] = (pairBuf[j * 2 + 1] == -1)
                    ? -1.0 : static_cast<double>(pairBuf[j * 2 + 1] - 1);
            }
        }

        /* ppg/ecg bin index pairs */
        readPairVec(bin.ppg_bin_indexs);
        readPairVec(bin.ecg_bin_indexs);

        bin.index = i;
        bool has_any_ecg = !bin.ecgSignal.empty()
            || !bin.ecgSignal2.empty()
            || !bin.ecgSignal3.empty();
        bool has_any_peaks = !bin.ch1.raw.empty() || !bin.ch1.squared.empty()
            || !bin.ch1.absval.empty();
        bin.bad_segment = !has_any_ecg || (!has_any_peaks && bin.ppgMinAmps.empty());
    }
    return results;
}

// ============================================================================
// Write template_info .bin
// ============================================================================
static void write_template_info_binfile(const std::string& path,
    const std::vector<TemplateInfo>& infos) {
    std::ofstream file(path, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Could not open for writing: " << path << std::endl;
        return;
    }

    char buf[1 << 16];
    file.rdbuf()->pubsetbuf(buf, sizeof(buf));

    uint64_t n = infos.size();
    file.write(reinterpret_cast<const char*>(&n), 8);

    auto writeDouble = [&](double v) {
        file.write(reinterpret_cast<const char*>(&v), 8);
        };

    auto writeVecDouble = [&](const std::vector<double>& v) {
        uint64_t sz = v.size();
        file.write(reinterpret_cast<const char*>(&sz), 8);
        if (sz > 0) file.write(reinterpret_cast<const char*>(v.data()), sz * 8);
        };

    auto writePairVec = [&](const std::vector<std::pair<uint64_t, uint64_t>>& v) {
        uint64_t sz = v.size();
        file.write(reinterpret_cast<const char*>(&sz), 8);
        if (sz > 0) file.write(reinterpret_cast<const char*>(v.data()), sz * 16);
        };

    for (const auto& info : infos) {
        uint64_t idx = info.index;
        file.write(reinterpret_cast<const char*>(&idx), 8);

        writePairVec(info.ppg_bin_indexs);
        writePairVec(info.ecg_bin_indexs);

        uint8_t bad = info.bad_segment ? 1 : 0;
        file.write(reinterpret_cast<const char*>(&bad), 1);

        auto writeChannel = [&](const ChannelTemplates& ch) {
            writeVecDouble(ch.ecgTemplate_raw);
            writeVecDouble(ch.ecgTemplate_squared);
            writeVecDouble(ch.ecgTemplate_absval);
            writeDouble(ch.alignment_point_raw);
            writeDouble(ch.alignment_point_squared);
            writeDouble(ch.alignment_point_absval);
            writeDouble(ch.avg_r_expand_raw);
            writeDouble(ch.avg_r_expand_squared);
            writeDouble(ch.avg_r_expand_absval);
            };

        writeChannel(info.ch1);
        writeChannel(info.ch2);
        writeChannel(info.ch3);

        writeVecDouble(info.ppgTemplate);
    }
}

// ============================================================================
// Config
// ============================================================================
struct TemplateConfig {
    std::string dataType;
    std::string wavePath;
    std::string templatePath;
};

static TemplateConfig parseConfig(const std::string& configPath, int type_row) {
    std::ifstream file(configPath);
    if (!file.is_open()) throw std::runtime_error("Could not open config.csv");

    std::string line;
    std::getline(file, line); // Skip header

    int currentIdx = 0;
    while (std::getline(file, line)) {
        if (currentIdx++ == type_row - 1) {
            std::stringstream ss(line);
            std::vector<std::string> row;
            std::string val;
            while (std::getline(ss, val, ',')) row.push_back(val);

            if (row.size() < 9)
                throw std::runtime_error(
                    "Config row needs at least 9 columns (0-8). "
                    "Col 7=wavePath, 8=templatePath");

            return { row[0], row[7], row[8] };
        }
    }
    throw std::runtime_error("Config type row not found");
}

// ============================================================================
// Main
// ============================================================================
int main() {
    try {
        int choice;
        std::cout << "Enter Data Type:\n1) MESA\n2) Bittium\n3) CHAOS" << std::endl;
        if (!(std::cin >> choice)) return 1;

        TemplateConfig cfg = parseConfig("config.csv", choice);

        if (!std::filesystem::exists(cfg.templatePath)) {
            std::filesystem::create_directories(cfg.templatePath);
        }

        for (const auto& entry : std::filesystem::directory_iterator(cfg.wavePath)) {
            if (entry.path().extension() != ".bin") continue;

            std::string stem = entry.path().stem().string();
            auto pos = stem.find("_wave_markings");
            if (pos == std::string::npos) continue;

            std::string currentID = stem.substr(0, pos);

            std::string outputPath = cfg.templatePath + "/" +
                currentID + "_template_info.bin";
       

            std::cout << "Processing templates: " << currentID << "..." << std::endl;

            auto wave_data = read_wave_data_binfile(entry.path().string());

            auto template_info = GenerateTemplates(wave_data);

            write_template_info_binfile(outputPath, template_info);
        }

        std::cout << "Template Generation Complete." << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "CRITICAL ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}