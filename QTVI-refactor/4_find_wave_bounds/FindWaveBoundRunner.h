// ============================================================================
// File: FindWaveBoundRunner.h
// Runner for FindWaveBounds processing
// ============================================================================
#ifndef FINDWAVEBOUNDRUNNER_H
#define FINDWAVEBOUNDRUNNER_H

#include "SignalProcessingTypes.h"

// Process all files in directory
void FindWaveBoundRunner(const string& annealedSegmentsPath,
    const string& outputPath,
    bool skipExisting);

#endif // FINDWAVEBOUNDRUNNER_H
