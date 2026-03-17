#pragma once

#include "common.hpp"

std::vector<std::pair<uint64_t, uint64_t>> getExclusionIntervals(
    uint64_t a,
    uint64_t b,
    const std::vector<std::pair<uint64_t, uint64_t>>& breaks);
