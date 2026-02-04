#ifndef FINDWAVEBOUNDS_H
#define FINDWAVEBOUNDS_H

#include <string>
#include <vector>
#include <map>

// Structure to hold wave data (equivalent to Matlab struct)
struct WaveData {
    // Add fields as needed based on what FindWaveBounds_EKGandPPG returns
    // This would need to match the actual Matlab structure
    std::vector<double> data;
    // Add other necessary fields
};

// Structure to hold annealed segments
struct AnnealedSegments {
    // Add fields to match the Matlab annealedSegments structure
    std::vector<double> segments;
    // Add other necessary fields
};

// Function declarations
int FindWaveBounds(const std::string& anneal_path, const std::string& output_path);
WaveData FindWaveBounds_EKGandPPG(const AnnealedSegments& annealedSegments, int param1, int param2);

#endif // FINDWAVEBOUNDS_H
