#ifndef LOAD_CONFIG_HPP
#define LOAD_CONFIG_HPP

#include <string>
#include <vector>
#include <array>
#include <cstdint>
#include <filesystem>

struct DataConfig {
    std::string data_type;
    std::string wave_bin_path;
    std::string generate_template_path;
    int ecg_sampling_rate = 0;
    int ppg_sampling_rate = 0;
};

// Parses config.csv, prompts the user to select a dataset, and returns the chosen config
DataConfig load_config(const std::filesystem::path& config_path);

struct WaveBin {
    uint64_t ecg_sampling_rate = 0;
    uint64_t ppg_sampling_rate = 0;
    std::vector<size_t> ecg_r_index;
    std::vector<size_t> ppg_max_amps;
    std::vector<size_t> ppg_min_amps;
    std::vector<double> ppg_signal;
    std::vector<double> ecg_signal;
    std::vector<std::array<int64_t, 2>> pairs; // [ppg_idx, ecg_idx], -1 = NaN
};

struct WaveFile {
    std::string filename;
    std::vector<WaveBin> bins;
};

// Loads all wave bin files from the selected config's wave_bin_path directory
std::vector<WaveFile> load_bin_files(const DataConfig& config);

#endif // LOAD_CONFIG_HPP