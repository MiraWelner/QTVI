#ifndef RPEAKDETECT_H
#define RPEAKDETECT_H

#include <vector>
#include <cstddef>

using namespace std;

// Result structure to match MATLAB [R_index, hrv, R_t, R_amp, S_t, S_amp]
struct RPeakDetectResult {
    vector<size_t> R_index;  // Indices of R peaks
    vector<double> hrv;      // Heart Rate Variability (RR intervals)
    vector<double> R_t;      // Time of R peaks
    vector<double> R_amp;    // Amplitude of R peaks
    vector<double> S_t;      // Time of S points
    vector<double> S_amp;    // Amplitude of S points
};

/**
 * 1-1 Translation of rpeakdetect.m
 * @param data Input ECG signal vector
 * @param samp_freq Sampling frequency (Hz), default 256
 * @param thresh Threshold for integrated peaks, default 0.2
 * @param testmode Plotting flag (not used in C++), default 0
 */
RPeakDetectResult rpeakdetect(const vector<double>& data, double samp_freq = 256.0, double thresh = 0.2, int testmode = 0);

#endif // RPEAKDETECT_H
