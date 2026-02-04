#include "NoiseManager.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>

NoiseManager::NoiseManager(double fs) : m_sampleRate(fs) {}

void NoiseManager::reserve(size_t n) {
    m_segments.reserve(n);
}

void NoiseManager::addSegment(size_t start, size_t end, const std::string& label, const std::string& marking_type) {
    AnnotationSegment seg(std::min(start, end), std::max(start, end), label);
    seg.marking_type = marking_type;
    m_segments.push_back(seg);
}

void NoiseManager::exportCSV(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) return;

    file << "start_sample,end_sample,start_sec,end_sec,label, marking_type\n";
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
    /*
        This is very similar to the exportCSV function, except it writes a structured
		array to a bin and it uses numeric labels instead of string labels for efficiency.
    */
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) return;

    const uint64_t count = m_segments.size();
    file.write(reinterpret_cast<const char*>(&count), sizeof(count));

    for (const auto& seg : m_segments) {

		//convert the signal label to a numeric ID
        double labelId = 0.0;
        if (seg.label == "PPG")      labelId = 1.0;
        else if (seg.label == "ECG") labelId = 2.0;
        else if (seg.label == "BOTH") labelId = 3.0;

		//convert the marking type to a numeric ID
        double typeID = 0.0;
        if (seg.marking_type == "Noise/Artifact")      typeID = 1.0;
        else if (seg.marking_type == "AF") typeID = 2.0;
        else if (seg.marking_type == "SVT") typeID = 3.0;
        else if (seg.marking_type == "VT") typeID = 4.0;
        else if (seg.marking_type == "PVC") typeID = 5.0;
        else if (seg.marking_type == "PAC") typeID = 6.0;
        else if (seg.marking_type == "Benign Arrhythmia") typeID = 7.0;
        else if (seg.marking_type == "Significant Arrhythmia") typeID = 8.0;

        const double row[6] = {
            static_cast<double>(seg.startSample),
            static_cast<double>(seg.endSample),
            seg.startSample / m_sampleRate,
            seg.endSample / m_sampleRate,
            labelId,
            typeID,
        };

        file.write(reinterpret_cast<const char*>(row), sizeof(row));
    }
}

