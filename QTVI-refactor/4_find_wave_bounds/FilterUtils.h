// ============================================================================
// File: FilterUtils.h
// Digital filtering utilities including Butterworth filters and filtfilt
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"

// Butterworth lowpass or highpass filter design (3rd order)
void butter(int N, double Wn, const string& type, vector<double>& b, vector<double>& a);

// Butterworth bandpass filter design (Nth order, ZPK pipeline)
void butter(int N, const vector<double>& Wn, vector<double>& b, vector<double>& a);

// Direct Form II Transposed IIR filter
vector<double> filter(const vector<double>& b, const vector<double>& a, const vector<double>& x);

// Zero-phase forward and reverse filtering
vector<double> filtfilt(const vector<double>& b, const vector<double>& a, const vector<double>& x);

// Full convolution: output length = a.size() + b.size() - 1
vector<double> conv(const vector<double>& a, const vector<double>& b);

// 1D median filter
std::vector<double> medfilt1(const std::vector<double>& x, int n);