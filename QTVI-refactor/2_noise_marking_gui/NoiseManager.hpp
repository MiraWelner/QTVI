#ifndef NOISE_MANAGER_HPP
#define NOISE_MANAGER_HPP

#include <vector>
#include <string>
#include <cstddef>

struct AnnotationSegment {
    /*
        This struct represents a segment that has been selected by the user
        at representing an arrhythmia, noise, or artifact.
    */
    size_t startSample;
    size_t endSample;
    std::string label; // ECG, PPG
    std::string marking_type; // Noise/Artifact, Conduction Delay, AF, SVT, VT, PVC, PAC, Benign Arrhythmia, Significant Arrhythmia

    AnnotationSegment(size_t s, size_t e, const std::string& l)
        : startSample(s), endSample(e), label(l) {
    }
};

class NoiseManager {
public:
    explicit NoiseManager(double fs);
    void reserve(size_t n);
    void addSegment(size_t start, size_t end, const std::string& label, const std::string& marking_type);
    void exportCSV(const std::string& filename) const;
    void exportBinary(const std::string& filename) const;
    const std::vector<AnnotationSegment>& getSegments() const { return m_segments; }

private:
    std::vector<AnnotationSegment> m_segments;
    double m_sampleRate;
};

#endif
