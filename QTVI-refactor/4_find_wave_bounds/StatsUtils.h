// ============================================================================
// File: StatsUtils.h
// Statistical utility functions
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"

double mean(const vector<double>& x);
double std_dev(const vector<double>& x);
double median(const vector<double>& x);
double sum(const vector<double>& x);

pair<double, size_t> max_element_index(const vector<double>& x, size_t start, size_t end);
pair<double, size_t> min_element_index(const vector<double>& x, size_t start, size_t end);

vector<double> movmean(const vector<double>& data, size_t window);
vector<double> diff(const vector<double>& x);
vector<double> sort(const vector<double>& x);
vector<size_t> find(const vector<bool>& condition);

void detrend(vector<double>& x);