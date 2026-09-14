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
 // FOR loadSpans BELOW. Without <fstream> the std::ifstream is an incomplete
 // type, which MSVC reports as "operator '!' cannot be applied to an operand
 // of type std::basic_istream" and "read: function does not take 1 arguments"
 // rather than as a missing include.
#include <fstream>
#include <algorithm>   // std::swap, for a span drawn right-to-left
#include <cstring>     // std::memcpy, for reinterpreting the legacy row count below

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
    // TWO LAYOUTS, BOTH READ (see loadSpans below). The CURRENT layout starts
    // with the 8-byte magic and carries threshold/blanking per row; the LEGACY
    // (pre-versioned) layout is a bare uint64 row count followed by rows of
    // six doubles with no parameter columns at all. Both are accepted for
    // exactly the reason bin_chunk_loader.cpp's readNoiseMarkingsBin() accepts
    // both: refusing the legacy layout here does not protect the two fields it
    // lacks (this Span has no threshold/blanking fields to lose in the first
    // place -- see the comment on Span below), it just throws away every
    // legacy file's exclusion spans entirely, silently. See the comment on
    // loadSpans() for the full story.
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


    // =======================================================================
    // Read
    // =======================================================================
    //
    // The counterpart of annotation_handler::exportBinary, and INSIDE THIS
    // NAMESPACE on purpose: the format is defined once, by the Column enum
    // above, and both directions index by those names. A reader in another
    // file is exactly the writer/reader divergence the channel-code table at
    // the top of this header exists to prevent.
    //
    // WHY TEMPLATE GENERATION NEEDS IT. Section 4.6 partitions the morphology
    // bank by operator class BEFORE clustering, so the classes are an INPUT to
    // generation, not an annotation applied afterwards.
    // jbank::BinBankInput::mark_code is the field that carries them, and
    // nothing ever filled it -- because nothing read this file.
    //
    // NO THRESHOLD/BLANKING HERE. Span carries only what generation and the
    // anneal step's exclusion logic actually key on: the sample range, the
    // channel, and the annotation code. The parameter-edit VALUES (threshold,
    // blanking) live only on AnnotationSegment / the GUI's own override lists,
    // which is a completely separate path -- this struct never held them and
    // never needed to.
    struct Span {
        int64_t start_sample = 0;
        int64_t end_sample = 0;
        uint8_t channel_code = 0;
        uint8_t annotation_code = 0;
    };

    struct LoadResult {
        bool read = false;          // the file opened and a layout was recognized
        std::string path;
        std::string error;
        std::vector<Span> spans;
    };

    // A MISSING FILE IS NOT AN ERROR: the record was never noise-marked, every
    // slice stays unmarked, and the bank has one partition -- the behaviour
    // before partitioning existed.
    //
    // LEGACY (PRE-VERSIONED) FILES ARE NOW READ, NOT REFUSED. This function
    // used to refuse any file whose first 8 bytes did not match kMagic, on the
    // reasoning that the pre-versioned layout carries no threshold/blanking
    // and reading it on a guess would restore a parameter-edit span with those
    // values silently defaulted. That reasoning does not apply to this
    // function at all: Span (above) has no threshold or blanking fields --
    // this reader was refusing legacy files to protect data it never returns
    // in the first place. The actual effect of the refusal was that every
    // record noise-marked before the magic/threshold/blanking format existed
    // loaded here as "wrong magic", came back with an EMPTY span list, and
    // every one of its noise-marked regions silently kept contributing R
    // peaks, beats, and morphology votes to generation and templating --
    // while the GUI's own separate reader (bin_chunk_loader.cpp's
    // readNoiseMarkingsBin, which already tolerates the legacy layout) showed
    // the very same markings on screen as present and applied. Markings
    // visible in the GUI but silently not excluded from processing is exactly
    // the bug this fixes.
    //
    // A genuine parameter-edit span from before parameter storage existed has
    // nothing to recover regardless -- the threshold/blanking values were
    // never captured anywhere, in this format or any other -- and that is a
    // concern for AnnotationSegment/the GUI's override lists, not for this
    // Span-only reader.
    inline LoadResult loadSpans(const std::string& path) {
        LoadResult out;
        out.path = path;

        std::ifstream f(path, std::ios::binary);
        if (!f) { out.error = "not found"; return out; }

        char header[8] = {};
        if (!f.read(header, sizeof(header))) {
            out.error = "truncated header"; return out;
        }

        bool magicMatches = true;
        for (std::size_t i = 0; i < sizeof(kMagic); ++i)
            if (header[i] != kMagic[i]) { magicMatches = false; break; }

        if (!magicMatches) {
            // LEGACY LAYOUT. The same 8 bytes just read are a bare uint64 row
            // count instead of a magic; each row is six raw doubles --
            // start_sample, end_sample, start_sec, end_sec, channel_code,
            // annotation_code, in that order and nothing else. Mirrors
            // bin_chunk_loader.cpp's appendLegacyNoiseMarkingRows exactly
            // (same layout, same guard), except this reader only needs four
            // of those six columns.
            uint64_t legacyCount = 0;
            std::memcpy(&legacyCount, header, sizeof(legacyCount));
            // A COUNT OFF DISK IS NOT A COUNT UNTIL IT IS CHECKED -- same
            // guard as the current-format path below, and for the same
            // reason: an implausible value here is a sign this file is
            // neither layout, not something to read millions of rows against.
            if (legacyCount > (1ull << 22)) {
                out.error = "neither a valid magic header nor a plausible legacy row count";
                return out;
            }

            out.read = true;
            out.spans.reserve(static_cast<std::size_t>(legacyCount));
            enum LegacyColumn : int {
                kLStartSample = 0, kLEndSample, kLStartSec, kLEndSec,
                kLChannelCode, kLAnnotationCode, kLegacyColumns
            };
            for (uint64_t r = 0; r < legacyCount; ++r) {
                double row[kLegacyColumns];
                if (!f.read(reinterpret_cast<char*>(row), sizeof(row))) {
                    out.error = "truncated at row (legacy) " + std::to_string(r);
                    break;      // keep what parsed; the rest is unreadable
                }
                Span sp;
                sp.start_sample = static_cast<int64_t>(row[kLStartSample]);
                sp.end_sample = static_cast<int64_t>(row[kLEndSample]);
                sp.channel_code = static_cast<uint8_t>(row[kLChannelCode]);
                sp.annotation_code = static_cast<uint8_t>(row[kLAnnotationCode]);
                // Drawn right-to-left: the GUI stores the drag as-is.
                if (sp.end_sample < sp.start_sample)
                    std::swap(sp.start_sample, sp.end_sample);
                out.spans.push_back(sp);
            }
            return out;
        }

        uint32_t version = 0;
        uint64_t count = 0;
        if (!f.read(reinterpret_cast<char*>(&version), sizeof(version))
            || !f.read(reinterpret_cast<char*>(&count), sizeof(count))) {
            out.error = "truncated header"; return out;
        }
        if (version != kVersion) {
            out.error = "version " + std::to_string(version)
                + " != " + std::to_string(kVersion);
            return out;
        }
        // A COUNT OFF DISK IS NOT A COUNT UNTIL IT IS CHECKED. reserve() on an
        // unchecked value turns a truncated file into length_error from an
        // allocator instead of an error naming the file.
        if (count > (1ull << 22)) {
            out.error = "implausible row count"; return out;
        }

        out.read = true;
        out.spans.reserve(static_cast<std::size_t>(count));
        for (uint64_t r = 0; r < count; ++r) {
            std::array<double, kColumns> row{};
            if (!f.read(reinterpret_cast<char*>(row.data()),
                sizeof(double) * kColumns)) {
                out.error = "truncated at row " + std::to_string(r);
                break;      // keep what parsed; the rest is unreadable
            }
            Span sp;
            sp.start_sample = static_cast<int64_t>(row[kStartSample]);
            sp.end_sample = static_cast<int64_t>(row[kEndSample]);
            sp.channel_code = static_cast<uint8_t>(row[kChannelCode]);
            sp.annotation_code = static_cast<uint8_t>(row[kAnnotationCode]);
            // Drawn right-to-left: the GUI stores the drag as-is.
            if (sp.end_sample < sp.start_sample)
                std::swap(sp.start_sample, sp.end_sample);
            out.spans.push_back(sp);
        }
        return out;
    }

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