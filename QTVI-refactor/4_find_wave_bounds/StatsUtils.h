// ============================================================================
// File: StatsUtils.h
// Statistical utility functions
// ============================================================================
#ifndef STATSUTILS_H
#define STATSUTILS_H

#include "SignalProcessingTypes.h"

// Mean calculation
double mean(const vector<double>& x);

// Standard deviation
double std(const vector<double>& x);

// Median
double median(const vector<double>& x);

// Max element and index
pair<double, size_t> max_element_index(const vector<double>& x, size_t start, size_t end);

// Min element and index
pair<double, size_t> min_element_index(const vector<double>& x, size_t start, size_t end);

// Moving mean
vector<double> movmean(const vector<double>& data, size_t window);

// Diff function
vector<double> diff(const vector<double>& x);

// Sum function
double sum(const vector<double>& x);

// Sort function
vector<double> sort(const vector<double>& x);

// Find indices where condition is true
vector<size_t> find(const vector<bool>& condition);

#endif // STATSUTILS_H
#pragma once
