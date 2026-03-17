#pragma once

#include "common.hpp"

std::vector<int> RoundToClosestBin(
    const std::vector<uint64_t>& bin_breaks,
    const std::vector<uint64_t>& indices);