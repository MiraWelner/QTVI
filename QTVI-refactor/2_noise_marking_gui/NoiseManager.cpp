/**
 * @file   NoiseManager.cpp
 * @brief  Implementation of annotation segment storage and export.
 *
 *   labelId: 0=unknown, 1=PPG, 2=ECG1, 3=ECG2, 4=ECG3, 5=ABP
 *   typeId:  0=unknown, 1=Noise/Artifact, 2=AF, 3=SVT, 4=VT,
 *            5=PVC, 6=PAC, 7=Benign Arrhythmia, 8=Significant Arrhythmia
 */
#include "NoiseManager.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <unordered_map>

NoiseManager::NoiseManager(double fs) : m_sampleRate(fs) {}

void NoiseManager::reserve(size_t n) { m_segments.reserve(n); }

void NoiseManager::addSegment(size_t start, size_t end,
    const std::string& label,
    const std::string& marking_type) {
    AnnotationSegment seg(std::min(start, end), std::max(start, end), label);
    seg.marking_type = marking_type;
    m_segments.push_back(seg);
}

void NoiseManager::exportCSV(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) return;

    file << "start_sample,end_sample,start_sec,end_sec,label,marking_type\n";
    for (const auto& seg : m_segments) {
        file << seg.startSample << ","
            << seg.endSample << ","
            << std::fixed << std::setprecision(3)
            << (seg.startSample / m_sampleRate) << ","
            << (seg.endSample / m_sampleRate) << ","
            << seg.label << ","
            << seg.marking_type << "\n";
    }
}

void NoiseManager::exportBinary(const std::string& filename) const {
    static const std::unordered_map<std::string, double> labelMap = {
        {"PPG", 1.0}, {"ECG1", 2.0}, {"ECG2", 3.0}, {"ECG3", 4.0}, {"ABP", 5.0}
    };
    static const std::unordered_map<std::string, double> typeMap = {
        {"Noise/Artifact", 1.0}, {"AF", 2.0}, {"SVT", 3.0}, {"VT", 4.0},
        {"PVC", 5.0}, {"PAC", 6.0},
        {"Benign Arrhythmia", 7.0}, {"Significant Arrhythmia", 8.0}
    };

    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) return;

    uint64_t count = m_segments.size();
    file.write(reinterpret_cast<const char*>(&count), sizeof(count));

    for (const auto& seg : m_segments) {
        auto labelIt = labelMap.find(seg.label);
        auto typeIt = typeMap.find(seg.marking_type);

        const double row[6] = {
            static_cast<double>(seg.startSample),
            static_cast<double>(seg.endSample),
            seg.startSample / m_sampleRate,
            seg.endSample / m_sampleRate,
            (labelIt != labelMap.end()) ? labelIt->second : 0.0,
            (typeIt != typeMap.end()) ? typeIt->second : 0.0,
        };
        file.write(reinterpret_cast<const char*>(row), sizeof(row));
    }
}