/**
 * @file   user_annotation_handler.h
 * @brief  Manages annotation segments (noise, arrhythmia, artifacts) and
 *         exports them to CSV and binary format.
 */
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

// ===========================================================================
// Channel codes
// ===========================================================================
//
// ONE TABLE. The writer previously held a static unordered_map<string,double>
// and the reader a switch statement over the same eight pairs, in a different
// file. Two copies of a mapping that has to agree exactly, with nothing
// checking that it did -- and both sides fail quietly: an unknown label writes
// code 0, an unknown code makes the reader skip the row. Add a channel to one
// and not the other and spans vanish on reload with no error.
//
// Same reasoning, and the same fix, as annotation_types.hpp: the table is the
// single definition and both directions are looked up from it.
namespace noise_markings {

    struct ChannelCode {
        const char* label;
        uint8_t     code;      // 0 is reserved for "unknown"
    };

    inline constexpr std::array<ChannelCode, 8> channel_codes = { {
        { "PPG",      1 },
        { "ECG1",     2 },
        { "ECG2",     3 },
        { "ECG3",     4 },
        { "ABP",      5 },
        { "ACCEL",    6 },
        { "ART",      7 },
        { "ART_PULM", 8 },
    } };

    constexpr bool labels_equal(const char* a, const char* b) {
        while (*a && *a == *b) { ++a; ++b; }
        return *a == *b;
    }

    // 0 when the label is not in the table. Callers must treat 0 as a refusal
    // to write the row, not as a channel.
    constexpr uint8_t code_for_channel(const char* label) {
        for (const auto& c : channel_codes)
            if (labels_equal(c.label, label)) return c.code;
        return 0;
    }
    inline uint8_t code_for_channel(const std::string& label) {
        return code_for_channel(label.c_str());
    }

    // nullptr when the code is not in the table.
    inline const char* channel_for_code(uint8_t code) {
        for (const auto& c : channel_codes)
            if (c.code == code) return c.label;
        return nullptr;
    }

    constexpr bool channel_codes_valid() {
        for (std::size_t i = 0; i < channel_codes.size(); ++i) {
            if (channel_codes[i].code == 0) return false;          // 0 = unknown
            if (channel_codes[i].label == nullptr
                || channel_codes[i].label[0] == '\0') return false;
            for (std::size_t j = i + 1; j < channel_codes.size(); ++j) {
                if (channel_codes[i].code == channel_codes[j].code) return false;
                if (labels_equal(channel_codes[i].label, channel_codes[j].label))
                    return false;
            }
        }
        return true;
    }
    static_assert(channel_codes_valid(),
        "noise_markings: channel codes must be unique, non-zero, and uniquely labelled");

    // =======================================================================
    // Binary format
    // =======================================================================
    //
    // ONE FORMAT, NO LEGACY PATH. The pre-versioned layout -- a bare uint64 row
    // count followed by rows of six doubles -- is not read. A file without the
    // magic below is REFUSED, with a message naming it, rather than parsed on a
    // guess: those files carry no threshold or blanking values, so loading one
    // would silently restore every parameter-edit span at the config defaults
    // and move the R peaks inside it. Refusing says so; reading does not.
    // Re-save affected files from the CSV, or re-mark them.
    //
    // EVERY FIELD IS A DOUBLE, including the two codes and the two sample
    // indices. One homogeneous row of kColumns doubles reads and writes as a
    // single contiguous block, so there is no field-by-field traversal to keep
    // in step, no struct padding to declare, and no sizeof assertion to
    // maintain. The column meanings are named by the enum below rather than left
    // to a comment, which is what makes the block self-describing at both ends.
    //
    // Integer exactness is not at risk at these magnitudes: a double holds every
    // integer below 2^53, and a sample index reaches that after roughly 285,000
    // years at 1 kHz.
    inline constexpr char     kMagic[8] = { 'N','M','K','B','0','0','0','1' };
    inline constexpr uint32_t kVersion = 1;

    // Column order. Reader and writer both index by these names, so neither can
    // drift from the other by a position.
    enum Column : int {
        kStartSample = 0,
        kEndSample,
        kStartSec,
        kEndSec,
        kChannelCode,
        kAnnotationCode,
        // NaN where the marking type carries no parameters, which is every type
        // but the paramEdit one. NaN rather than 0 because a threshold of 0 is a
        // meaningful value (detect everything) and has to stay distinguishable
        // from "this row never carried one". NaN round-trips exactly through an
        // IEEE-754 double, so no separate presence flag is needed.
        kThreshold,
        kBlankingMs,
        kColumns
    };

}  // namespace noise_markings

struct AnnotationSegment {
    int startSample;
    int endSample;
    std::string label;
    std::string marking_type;
    double sampleRate;

    // THE PARAMETER VALUES THAT MADE THIS SPAN MEAN SOMETHING.
    //
    // A parameter-edit span is not an observation about the signal, it is an
    // instruction: inside it, use THIS detection threshold and THIS blanking
    // period instead of the config defaults. The span was being saved and the
    // two numbers were not, so a reloaded file restored the highlight and lost
    // what it meant -- every override reverted to cfg.threshold and
    // cfg.blanking_period, and the R peaks inside it moved. Nothing reported
    // that, because a span with the right extent and the wrong parameters looks
    // identical on screen.
    //
    // NaN = not applicable to this marking type. See noise_markings::Row.
    double threshold = std::numeric_limits<double>::quiet_NaN();
    double blanking = std::numeric_limits<double>::quiet_NaN();

    AnnotationSegment(int s, int e, const std::string& l, double sr)
        : startSample(s), endSample(e), label(l), sampleRate(sr) {
    }
};

class annotation_handler {
public:
    annotation_handler() = default;
    void reserve(int n);
    void addSegment(int start, int end, const std::string& label, const std::string& marking_type, double sampleRate);
    // Overload for parameter-edit spans. Pass the threshold and blanking period
    // that apply inside the span; both are written to the CSV and the binary and
    // restored on reload.
    void addSegment(int start, int end, const std::string& label, const std::string& marking_type,
        double sampleRate, double threshold, double blanking);
    void exportCSV(const std::string& filename)    const;
    void exportBinary(const std::string& filename) const;
    const std::vector<AnnotationSegment>& getSegments() const { return m_segments; }

private:
    std::vector<AnnotationSegment> m_segments;
};