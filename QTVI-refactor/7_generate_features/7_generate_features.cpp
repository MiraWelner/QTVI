/**
 * @file main_generate_features.hpp
 * @brief Entry point for the PPG feature extraction pipeline. This pipeline takes in the annealed data created in 
 *        step 3, the R peak marking data created in step 4, and the PPG peak and endpoint data created in step 6
 *
 */

#include "runner.hpp"
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    try {
        std::string config_csv = "config.csv";
        ppg::run_pipeline(config_csv);
    }
    catch (const std::exception& e) {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
    return 0;
}