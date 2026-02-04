// FindWaveBoundRunner.cpp
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <filesystem>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include "FindWaveBounds.h"
#include "ConfigReader.h"

namespace fs = std::filesystem;

// Structure to hold file information
struct FileInfo {
    std::string name;
    std::string folder;
    std::string full_path;
};

// Function to recursively find files matching a pattern
std::vector<FileInfo> findFiles(const std::string& base_path, const std::string& pattern) {
    std::vector<FileInfo> files;

    if (!fs::exists(base_path)) {
        std::cerr << "Error: Path does not exist: " << base_path << std::endl;
        return files;
    }

    for (const auto& entry : fs::recursive_directory_iterator(base_path)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();
            if (filename.find(pattern) != std::string::npos) {
                FileInfo info;
                info.name = filename;
                info.folder = entry.path().parent_path().string();
                info.full_path = entry.path().string();
                files.push_back(info);
            }
        }
    }

    return files;
}

int main() {
    // clear (variables are local in C++, so no need for explicit clear)

    // Read properties from config file
    std::map<std::string, std::string> props = readProps("config.txt");

    std::string Annealed_segments = props["FWB_input_path"];
    std::string output_path = props["FWB_output_path"];
    bool Skip_Existing = (props["Skip_Existing"] == "1" || props["Skip_Existing"] == "true");

    // Find all files matching the pattern
    std::vector<FileInfo> list = findFiles(Annealed_segments, "annealedSegments.mat");

    // Flip the list (reverse it)
    std::reverse(list.begin(), list.end());

    double time = 0.0;

    for (size_t i = 0; i < list.size(); i++) {
        // set name
        std::string name = list[i].name;

        // Find the start index of '_annealedSegments' in the name
        std::regex pattern("_annealedSegments");
        std::smatch match;
        if (std::regex_search(name, match, pattern)) {
            size_t start_idx = match.position(0);
            name = name.substr(0, start_idx);
        }

        // Check if output file exists and Skip_Existing is true
        std::string output_file = (fs::path(output_path) / (name + "_wave_data.mat")).string();
        if (Skip_Existing && fs::exists(output_file)) {
            std::cout << name << "_wave_data.mat exists skipping because Skip_Existing = 1 in config." << std::endl;
            continue;
        }

        std::string file = list[i].full_path;

        // Start timer
        auto tStart = std::chrono::high_resolution_clock::now();

        std::cout << "Beginning analysis of " << file << " | " << (i + 1) << " of " << list.size() << std::endl;

        double avg_time = time / (i + 1);
        std::cout << "Avg Time (s): " << avg_time << std::endl;

        double est_finish = (avg_time * (list.size() - i - 1)) / 60.0;
        std::cout << "Est finish (min): " << est_finish << std::endl;

        std::cout << "Output loc " << output_path << std::endl;

        int r = FindWaveBounds(file, output_path);

        // Commented out code from original:
        // if (r == 1) {
        //     // movefile equivalent would be fs::rename or fs::copy + fs::remove
        //     // fs::rename(file, "/hdd/data/mesa/Manuscript 1/2020 November run/2 Anneal Segs/");
        // }

        std::cout << "____________________________________________________________________________________________________" << std::endl << std::endl;

        // Stop timer and add elapsed time
        auto tEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = tEnd - tStart;
        time += elapsed.count();
    }

    return 0;
}