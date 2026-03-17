#pragma once

#include "common.hpp"

void splitOverlappingBins(
    std::vector<Exclusion>& exclusions,
    const std::vector<uint64_t>& bin_breaks);