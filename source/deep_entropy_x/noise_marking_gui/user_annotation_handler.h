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
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <cstring>

namespace noise_markings {

    struct ChannelCode {
        //label = humanreadable name, code = number, 0=unknown (it shouldn't ever be this)
        const char* label;
        uint8_t     code;
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


    inline constexpr char     noise_marking_magic[8] = { 'N','M','K','B','0','0','0','1' };
    inline constexpr uint32_t noise_marking_version = 0;
    // VERSION 1 IS THE SAME LAYOUT. Files on disk (e.g. the MESA bl_dmu set)
    // were written by a build whose constant was 1; it was later set to 0 with
    // no change to the header or the eight columns, so a v1 file reads exactly
    // as a v0 one. Accepted, never written. Any OTHER version is still refused.
    inline constexpr uint32_t noise_marking_version_same_layout = 1;

    // Column order. Reader and writer both index by these names, so neither can
    // drift from the other by a position.
    enum Column : int {
        start_location_in_samples = 0,
        end_location_in_samples,
        start_location_in_seconds,
        end_location_in_seconds,
        channel_marked,
        marking_type,
        threshold_value,
        blanking_miliseconds,
        columns
    };

    struct Span {
        int64_t start_sample = 0;
        int64_t end_sample = 0;
        uint8_t channel_code = 0;
        uint8_t annotation_code = 0;
    };

    struct LoadResult {
        bool read = false;          // the file opened, magic and version matched
        std::string path;
        std::string error;
        std::vector<Span> spans;
    };

    struct Row {
        double  start_sample = 0.0;
        double  end_sample = 0.0;
        double  start_sec = 0.0;
        double  end_sec = 0.0;
        uint8_t channel_code = 0;
        uint8_t annotation_code = 0;
        double  threshold = std::numeric_limits<double>::quiet_NaN();
        double  blanking_ms = std::numeric_limits<double>::quiet_NaN();
    };

    struct RowsResult {
        bool read = false;          // the file opened, and its layout was recognised
        // True for a pre-versioned file (see loadRows). Its rows carry no
        // threshold / blanking, so those read NaN.
        bool legacy = false;
        std::string path;
        std::string error;          // set on failure, or on a partial read
        std::vector<Row> rows;
    };

    inline RowsResult loadRows(const std::string& path) {
        RowsResult out;
        out.path = path;

        std::ifstream f(path, std::ios::binary);
        if (!f) { out.error = "not found"; return out; }

        f.seekg(0, std::ios::end);
        const uint64_t fileSize = static_cast<uint64_t>(f.tellg());
        f.seekg(0, std::ios::beg);

        char header[sizeof(noise_marking_magic)] = {};
        if (!f.read(header, sizeof(header))) {
            out.error = "truncated header"; return out;
        }

        // ---- TWO LAYOUTS, TOLD APART BY THE FIRST EIGHT BYTES ----------------
        //
        // CURRENT: the magic, a version, a row count, then rows of `columns`
        // doubles indexed by Column.
        //
        // LEGACY (pre-versioned, no magic): a bare uint64 row count, then rows
        // of SIX doubles -- start_sample, end_sample, start_sec, end_sec,
        // channel_code, annotation_code -- the first six Columns, in the same
        // order. These files predate threshold/blanking storage, so a
        // parameter-edit span in one has no recorded value: both read NaN, the
        // same sentinel the current format uses for "not applicable", so the
        // GUI treats it as an unset override rather than defaulting it.
        //
        // The legacy count and the magic are both 8 bytes, so one read decides.
        // A non-magic header is accepted as a legacy count ONLY if the file is
        // exactly that many six-double rows long -- anything else is a damaged
        // or foreign file, refused rather than read as rows.
        //
        // CHANNEL CODES ARE READ WITH TODAY'S TABLE. A legacy file written
        // before the writer and reader channel maps were unified may use codes
        // that mean a different channel now; that cannot be detected here.
        bool isCurrent = true;
        for (std::size_t i = 0; i < sizeof(noise_marking_magic); ++i)
            if (header[i] != noise_marking_magic[i]) { isCurrent = false; break; }

        if (!isCurrent) {
            constexpr uint64_t kLegacyColumns = 6;
            uint64_t count = 0;
            std::memcpy(&count, header, sizeof(count));
            if (count > (1ull << 22)
                || fileSize != sizeof(count) + count * kLegacyColumns * sizeof(double)) {
                out.error = "bad magic, and not a legacy noise-markings file either "
                    "(size does not match its row count) -- re-export it from the "
                    "marking GUI";
                return out;
            }
            out.read = true;
            out.legacy = true;
            out.rows.reserve(static_cast<std::size_t>(count));
            for (uint64_t r = 0; r < count; ++r) {
                std::array<double, kLegacyColumns> row{};
                if (!f.read(reinterpret_cast<char*>(row.data()),
                    sizeof(double) * kLegacyColumns)) {
                    out.error = "truncated at row " + std::to_string(r);
                    break;
                }
                Row rw;   // threshold / blanking_ms stay NaN
                rw.start_sample = row[start_location_in_samples];
                rw.end_sample = row[end_location_in_samples];
                rw.start_sec = row[start_location_in_seconds];
                rw.end_sec = row[end_location_in_seconds];
                rw.channel_code = static_cast<uint8_t>(row[channel_marked]);
                rw.annotation_code = static_cast<uint8_t>(row[marking_type]);
                out.rows.push_back(rw);
            }
            std::fprintf(stderr, "[noise-markings] %s is a legacy (pre-versioned) "
                "file: %zu row(s), no threshold/blanking. Re-export it from the "
                "marking GUI to upgrade it.\n", path.c_str(), out.rows.size());
            return out;
        }

        uint32_t version = 0;
        uint64_t count = 0;
        if (!f.read(reinterpret_cast<char*>(&version), sizeof(version))
            || !f.read(reinterpret_cast<char*>(&count), sizeof(count))) {
            out.error = "truncated header"; return out;
        }
        // A version mismatch is a hard stop, not a best-effort read: guessing
        // the row stride is how spans land at the wrong times while looking
        // plausible.
        if (version != noise_marking_version
            && version != noise_marking_version_same_layout) {
            out.error = "version " + std::to_string(version) + " is not "
                + std::to_string(noise_marking_version) + " or "
                + std::to_string(noise_marking_version_same_layout);
            return out;
        }
        // A COUNT OFF DISK IS NOT A COUNT UNTIL IT IS CHECKED. reserve() on an
        // unchecked value turns a truncated file into length_error from an
        // allocator instead of an error naming the file.
        if (count > (1ull << 22)) {
            out.error = "implausible row count"; return out;
        }

        out.read = true;
        out.rows.reserve(static_cast<std::size_t>(count));
        for (uint64_t r = 0; r < count; ++r) {
            std::array<double, columns> row{};
            if (!f.read(reinterpret_cast<char*>(row.data()),
                sizeof(double) * columns)) {
                out.error = "truncated at row " + std::to_string(r);
                break;      // keep what parsed; the rest is unreadable
            }
            Row rw;
            rw.start_sample = row[start_location_in_samples];
            rw.end_sample = row[end_location_in_samples];
            rw.start_sec = row[start_location_in_seconds];
            rw.end_sec = row[end_location_in_seconds];
            rw.channel_code = static_cast<uint8_t>(row[channel_marked]);
            rw.annotation_code = static_cast<uint8_t>(row[marking_type]);
            rw.threshold = row[threshold_value];
            rw.blanking_ms = row[blanking_miliseconds];
            out.rows.push_back(rw);
        }
        return out;
    }

    // ADAPTER. Sample indices and raw codes only; make_averaged_templates
    // partitions morphology by annotation_code, so it needs the code rather
    // than a label, and works in samples rather than seconds.
    inline LoadResult loadSpans(const std::string& path) {
        const RowsResult rr = loadRows(path);

        LoadResult out;
        out.path = rr.path;
        out.error = rr.error;
        out.read = rr.read;
        out.spans.reserve(rr.rows.size());
        for (const Row& rw : rr.rows) {
            Span sp;
            sp.start_sample = static_cast<int64_t>(rw.start_sample);
            sp.end_sample = static_cast<int64_t>(rw.end_sample);
            sp.channel_code = rw.channel_code;
            sp.annotation_code = rw.annotation_code;
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
    void export_marking_binfile(const std::string& filename) const;
    const std::vector<AnnotationSegment>& getSegments() const { return m_segments; }

private:
    std::vector<AnnotationSegment> m_segments;
};