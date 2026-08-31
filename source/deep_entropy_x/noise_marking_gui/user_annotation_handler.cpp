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
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

void annotation_handler::reserve(int n) { m_segments.reserve(n); }

void annotation_handler::addSegment(int start, int end,
    const std::string& label, const std::string& marking_type, double sampleRate) {
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    addSegment(start, end, label, marking_type, sampleRate, NaN, NaN);
}

void annotation_handler::addSegment(int start, int end,
    const std::string& label, const std::string& marking_type, double sampleRate,
    double threshold, double blanking) {
    AnnotationSegment seg(std::min(start, end), std::max(start, end), label, sampleRate);
    seg.marking_type = marking_type;
    seg.threshold = threshold;
    seg.blanking = blanking;
    m_segments.push_back(seg);
}

void annotation_handler::exportCSV(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) return;

    // Parameters last, after the identity of the span. Reading left to right the
    // row says where, on what, of what kind, and then under what settings --
    // which is also the order in which those fields become relevant.
    file << "start_sample,end_sample,start_sec,end_sec,label,marking_type,"
            "threshold,blanking_ms\n";
    for (const auto& seg : m_segments) {
        file << seg.startSample << ","
            << seg.endSample << ","
            << std::fixed << std::setprecision(3)
            << (seg.startSample / seg.sampleRate) << ","
            << (seg.endSample / seg.sampleRate) << ","
            << seg.label << ","
            << seg.marking_type << ",";
        // EMPTY, not 0 and not "nan". A blank field round-trips through every
        // CSV reader as missing, which is what it is; "nan" parses as a number
        // in some and as a string in others, and 0 is a legitimate threshold.
        if (!std::isnan(seg.threshold)) file << seg.threshold;
        file << ",";
        if (!std::isnan(seg.blanking)) file << seg.blanking;
        file << "\n";
    }
}

void annotation_handler::exportBinary(const std::string& filename) const {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) return;

    namespace nm = noise_markings;

    // Rows are built first, so the count in the header agrees with what follows
    // it. A row whose channel label or marking type is not in its table is
    // dropped and counted: both lookups come from the tables this build defines,
    // so a miss is a programming error rather than operator data, and it is
    // reported rather than silently omitted -- silent omission is what hid the
    // writer/reader channel-map divergence before the two were unified.
    std::vector<std::array<double, nm::kColumns>> rows;
    rows.reserve(m_segments.size());
    std::size_t unknown_channel = 0, unknown_type = 0;

    for (const auto& seg : m_segments) {
        const uint8_t chan = nm::code_for_channel(seg.label);
        if (chan == 0) { ++unknown_channel; continue; }
        const int type = annotation_types::markCode(seg.marking_type);
        if (type == 0) { ++unknown_type; continue; }

        std::array<double, nm::kColumns> row{};
        row[nm::kStartSample] = static_cast<double>(seg.startSample);
        row[nm::kEndSample] = static_cast<double>(seg.endSample);
        row[nm::kStartSec] = seg.startSample / seg.sampleRate;
        row[nm::kEndSec] = seg.endSample / seg.sampleRate;
        row[nm::kChannelCode] = static_cast<double>(chan);
        row[nm::kAnnotationCode] = static_cast<double>(type);
        row[nm::kThreshold] = seg.threshold;
        row[nm::kBlankingMs] = seg.blanking;
        rows.push_back(row);
    }

    if (unknown_channel || unknown_type)
        std::fprintf(stderr,
            "[noise-markings] dropped %zu row(s) with an unknown channel and "
            "%zu with an unknown marking type while writing %s\n",
            unknown_channel, unknown_type, filename.c_str());

    const uint64_t count = rows.size();
    file.write(nm::kMagic, sizeof(nm::kMagic));
    file.write(reinterpret_cast<const char*>(&nm::kVersion), sizeof(nm::kVersion));
    file.write(reinterpret_cast<const char*>(&count), sizeof(count));
    for (const auto& row : rows)
        file.write(reinterpret_cast<const char*>(row.data()),
            sizeof(double) * nm::kColumns);
}
