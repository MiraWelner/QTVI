#pragma once

#include "common.hpp"

std::vector<Section> markForMovement(
    const std::vector<Exclusion>& exclusions,
    const std::vector<uint64_t>& bin_breaks,
    int bin_count,
    double ppgSR,
    uint64_t bin_size_samples,
    double min_bin_size_mins);
