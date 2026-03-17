#include "load_config.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <map>

static std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

static std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ',')) {
        fields.push_back(trim(field));
    }
    return fields;
}

DataConfig load_config(const std::filesystem::path& config_path) {
    std::ifstream file(config_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open config file: " + config_path.string());
    }

    std::string header_line;
    if (!std::getline(file, header_line)) {
        throw std::runtime_error("Config file is empty");
    }
    std::vector<std::string> headers = split_csv(header_line);

    std::map<std::string, size_t> col_idx;
    for (size_t i = 0; i < headers.size(); ++i) {
        col_idx[headers[i]] = i;
    }

    auto get_field = [&](const std::vector<std::string>& row, const std::string& col) -> std::string {
        auto it = col_idx.find(col);
        if (it == col_idx.end() || it->second >= row.size()) return "";
        return row[it->second];
        };

    auto get_int_field = [&](const std::vector<std::string>& row, const std::string& col) -> int {
        std::string val = get_field(row, col);
        if (val.empty()) return 0;
        try { return std::stoi(val); }
        catch (...) { return 0; }
        };

    // Parse all rows
    std::vector<DataConfig> configs;
    std::string line;
    while (std::getline(file, line)) {
        if (trim(line).empty()) continue;
        std::vector<std::string> fields = split_csv(line);
        DataConfig cfg;
        cfg.data_type = get_field(fields, "DATA_TYPE");
        cfg.wave_bin_path = get_field(fields, "wave_bin_path");
        cfg.generate_template_path = get_field(fields, "generate_template_path");
        cfg.ecg_sampling_rate = get_int_field(fields, "ecg_sampling_rate");
        cfg.ppg_sampling_rate = get_int_field(fields, "ppg_sampling_rate");
        if (!cfg.data_type.empty()) {
            configs.push_back(cfg);
        }
    }

    if (configs.empty()) {
        throw std::runtime_error("No data rows found in config file");
    }

    // Prompt user to select
    std::cout << "\n=== Dataset Selection ===\n";
    for (size_t i = 0; i < configs.size(); ++i) {
        std::cout << "  [" << (i + 1) << "] " << configs[i].data_type << "\n";
    }
    std::cout << "\nSelect a dataset (1-" << configs.size() << "): ";

    int choice = 0;
    while (true) {
        std::string input;
        std::getline(std::cin, input);
        try {
            choice = std::stoi(trim(input));
        }
        catch (...) {
            choice = 0;
        }
        if (choice >= 1 && choice <= static_cast<int>(configs.size())) break;
        std::cout << "Invalid selection. Please enter 1-" << configs.size() << ": ";
    }

    const DataConfig& selected = configs[choice - 1];

    return selected;
}

static WaveFile load_single_bin_file(const std::filesystem::path& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open bin file: " + path.string());
    }

    auto read_u64 = [&]() -> uint64_t {
        uint64_t val;
        file.read(reinterpret_cast<char*>(&val), 8);
        return val;
        };

    auto read_index_vec = [&]() -> std::vector<size_t> {
        uint64_t sz = read_u64();
        std::vector<size_t> v(sz);
        return v;
    };

    auto read_signal_vec = [&]() -> std::vector<double> {
        uint64_t sz = read_u64();
        std::vector<double> v(sz);
        if (sz > 0) {
            file.read(reinterpret_cast<char*>(v.data()), sz * 8);
        }
        return v;
        };

    WaveFile wf;
    wf.filename = path.stem().string();

    uint64_t num_bins = read_u64();
    wf.bins.resize(num_bins);



    for (uint64_t b = 0; b < num_bins; ++b) {
        WaveBin& bin = wf.bins[b];
        bin.ecg_sampling_rate = read_u64();
        bin.ppg_sampling_rate = read_u64();

        auto read_index_vec = [&]() -> std::vector<size_t> {
            uint64_t sz = read_u64();
            std::vector<size_t> v(sz);
            for (uint64_t i = 0; i < sz; ++i) {
                v[i] = static_cast<size_t>(read_u64() - 1);
            }
            return v;
            };

        bin.ecg_r_index = read_index_vec();
        bin.ppg_max_amps = read_index_vec();
        bin.ppg_min_amps = read_index_vec();
        bin.ppg_signal = read_signal_vec();
        bin.ecg_signal = read_signal_vec();

        uint64_t num_pairs = read_u64();
        bin.pairs.resize(num_pairs);
        for (uint64_t p = 0; p < num_pairs; ++p) {
            int64_t ppg_idx, ecg_idx;
            file.read(reinterpret_cast<char*>(&ppg_idx), 8);
            file.read(reinterpret_cast<char*>(&ecg_idx), 8);
            bin.pairs[p][0] = (ppg_idx < 0) ? -1 : ppg_idx - 1;
            bin.pairs[p][1] = (ecg_idx < 0) ? -1 : ecg_idx - 1;
        }
    }

    return wf;
}

std::vector<WaveFile> load_bin_files(const DataConfig& config) {
    std::filesystem::path wave_dir(config.wave_bin_path);
    if (!std::filesystem::exists(wave_dir)) {
        throw std::runtime_error("Wave bin directory does not exist: " + wave_dir.string());
    }

    std::vector<std::filesystem::path> paths;
    for (const auto& entry : std::filesystem::directory_iterator(wave_dir)) {
        if (entry.is_regular_file()) {
            paths.push_back(entry.path());
        }
    }
    std::sort(paths.begin(), paths.end());

    std::cout << "Found " << paths.size() << " files in wave_bin_path" << std::endl;
    std::vector<WaveFile> wave_files;
    wave_files.reserve(paths.size());
    for (const auto& p : paths) {
        if (p.extension() == ".bin") {
            std::cout << "  Loading " << p.filename().string() << "...";
            wave_files.push_back(load_single_bin_file(p));
            std::cout << " " << wave_files.back().bins.size() << " bins\n";
        }
    }

    std::cout << "\nLoaded " << wave_files.size() << " wave files.\n\n";
    return wave_files;
}