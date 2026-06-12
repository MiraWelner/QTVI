/**
 * @file   user_annotation_handler.h
 * @brief  Manages annotation segments (noise, arrhythmia, artifacts) and
 *         exports them to CSV and binary format.
 */
#pragma once

#include <vector>
#include <string>
#include <cstddef>

struct AnnotationSegment {
    int      startSample;
    int      endSample;
    std::string label;         // "ECG1", "ECG2", "ECG3", "PPG", or "ABP"
    std::string marking_type;  // "Noise/Art.", "AF", "SVT", etc.

    AnnotationSegment(int s, int e, const std::string& l)
        : startSample(s), endSample(e), label(l) {
    }
};

class annotation_handler {
public:
    explicit annotation_handler(double fs);

    void reserve(int n);
    void addSegment(int start, int end,
        const std::string& label, const std::string& marking_type);

    void exportCSV(const std::string& filename)    const;
    void exportBinary(const std::string& filename) const;

    const std::vector<AnnotationSegment>& getSegments() const { return m_segments; }

private:
    std::vector<AnnotationSegment> m_segments;
    double m_sampleRate;
};