#pragma once

#include "common.hpp"

std::vector<FinalSegment> AnnealSegments(
    const RawData& data,
    const std::vector<std::pair<double, double>>& noiseSEG,
    double targetLenMins);
