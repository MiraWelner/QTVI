
#include <Eigen/Dense>
#include <iostream>
#include "load_config.hpp"
 
int main() {
    DataConfig config = load_config("config.csv");
    std::vector<WaveFile> wave_files = load_bin_files(config);


    return 0;
}
