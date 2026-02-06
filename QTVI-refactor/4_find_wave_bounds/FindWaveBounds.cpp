// ============================================================================
// File: FindWaveBounds.cpp
// ============================================================================
#include "FindWaveBounds.h"
#include "FindWaveBounds_EKGandPPG.h"
#include <iostream>

// Helper function to extract name from path
string extractName(const string& path) {
    size_t lastSlash = path.find_last_of("/\\");
    string filename = (lastSlash == string::npos) ? path : path.substr(lastSlash + 1);

    size_t pos = filename.find("_annealedSegments");
    if (pos != string::npos) {
        return filename.substr(0, pos);
    }

    return filename;
}

int FindWaveBounds(const string& anneal_path, const string& output_path) {
    try {
        string name = extractName(anneal_path);

        // Load annealedSegments
        // NOTE: This would need to be implemented based on your file format
        // For now, this is a placeholder
        vector<AnnealedSegment> annealedSegments;
        // annealedSegments = loadAnnealedSegments(anneal_path);

        std::cout << "Finding individual beats..." << std::endl;
        vector<WaveData> wave_data = FindWaveBounds_EKGandPPG(annealedSegments, 0, true);

        std::cout << "Saving..." << std::endl;
        // Save wave_data
        // NOTE: This would need to be implemented based on your file format
        // saveWaveData(output_path + "/" + name + "_wave_data", wave_data);

        return 1;
    }
    catch (const std::exception& e) {
        std::cerr << "Error in FindWaveBounds: " << e.what() << std::endl;
        return 0;
    }
}