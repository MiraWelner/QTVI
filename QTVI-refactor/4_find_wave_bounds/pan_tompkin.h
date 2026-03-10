#ifndef PAN_TOMPKIN_H
#define PAN_TOMPKIN_H

#include <vector>
#include <fstream>

#include <cstddef>

using namespace std;

// Pan-Tompkins QRS detection result
struct PanTompkinResult {
    vector<size_t> qrs_i_raw;   // index of R waves
    vector<double> qrs_amp_raw; // amplitude of R waves
    int delay;                  // delay in samples
};

/**
 * Complete 1-1 translation of Sedghamiz's Pan-Tompkins algorithm.
 */

PanTompkinResult pan_tompkin(const std::vector<double>& ecg_input, double fs, int gr, const std::string& fileID);

#endif // PAN_TOMPKIN_H
