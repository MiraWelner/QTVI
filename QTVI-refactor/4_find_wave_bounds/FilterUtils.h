// ============================================================================
// File: FilterUtils.h
// Digital filtering utilities including Butterworth filters and filtfilt
// ============================================================================
#ifndef FILTERUTILS_H
#define FILTERUTILS_H

#include "SignalProcessingTypes.h"

// Butterworth filter design
void butter(int N, double Wn, const string& type, vector<double>& b, vector<double>& a);
void butter(int N, const vector<double>& Wn, vector<double>& b, vector<double>& a);

// Apply IIR filter
vector<double> filter(const vector<double>& b, const vector<double>& a, const vector<double>& x);

// Zero-phase forward and reverse filtering
vector<double> filtfilt(const vector<double>& b, const vector<double>& a, const vector<double>& x);

// Convolution
vector<double> conv(const vector<double>& a, const vector<double>& b);

// Median filter
vector<double> medfilt1(const vector<double>& x, int window_size);

#endif // FILTERUTILS_H#pragma once
