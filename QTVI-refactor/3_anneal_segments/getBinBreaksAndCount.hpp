#pragma once

#include "common.hpp"

BinBreaksResult getBinBreaksAndCount(
    uint64_t total_len,
    uint64_t bin_size_samples,
    double ppgSR,
    double min_bin_size_mins);