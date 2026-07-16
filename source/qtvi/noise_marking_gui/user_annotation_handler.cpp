/**
 * @file   user_annotation_handler.cpp
 * @brief  handles the markings that the user annotates on all the markable channels
 *
 */
#include "user_annotation_handler.h"
#include "annotation_types.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <unordered_map>

void annotation_handler::reserve(int n) { m_segments.reserve(n); }

void annotation_handler::addSegment(int start, int end,
    const std::string& label, const std::string& marking_type, double sampleRate) {
    AnnotationSegment seg(std::min(start, end), std::max(start, end), label, sampleRate);
    seg.marking_type = marking_type;
    m_segments.push_back(seg);
}

void annotation_handler::exportCSV(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) return;

    file << "start_sample,end_sample,start_sec,end_sec,label,marking_type\n";
    for (const auto& seg : m_segments) {
        file << seg.startSample << ","
            << seg.endSample << ","
            << std::fixed << std::setprecision(3)
            << (seg.startSample / seg.sampleRate) << ","
            << (seg.endSample / seg.sampleRate) << ","
            << seg.label << ","
            << seg.marking_type << "\n";
    }
}

void annotation_handler::exportBinary(const std::string& filename) const {
    static const std::unordered_map<std::string, double> labelMap = {
        {"PPG", 1.0}, {"ECG1", 2.0}, {"ECG2", 3.0}, {"ECG3", 4.0}, {"ABP", 5.0}, {"ACCEL", 6.0}, {"ART", 7.0}, {"ART_PULM", 8.0}
    };

    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) return;

    uint64_t count = m_segments.size();
    file.write(reinterpret_cast<const char*>(&count), sizeof(count));

    for (const auto& seg : m_segments) {
        auto labelIt = labelMap.find(seg.label);

        const double row[6] = {
            static_cast<double>(seg.startSample),
            static_cast<double>(seg.endSample),
            seg.startSample / seg.sampleRate,
            seg.endSample / seg.sampleRate,
            (labelIt != labelMap.end()) ? labelIt->second : 0.0,
            annotation_types::codeFor(seg.marking_type),
        };
        file.write(reinterpret_cast<const char*>(row), sizeof(row));
    }
}