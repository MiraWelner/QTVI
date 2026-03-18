// ============================================================================
// File: FindWaveBounds_EKGandPPG.h
// Top-level processing: detect R-peaks and PPG pulses, pair them
// ============================================================================
#pragma once

#include <vector>
#include <string>
#include "FindWaveBounds.h"

std::vector<WaveData> FindWaveBounds_EKGandPPG(
    const std::vector<AnnealedSegment>& annealedSegments,
    int dbg_plot,
    bool use_R_algorithm,
    std::string fileID);