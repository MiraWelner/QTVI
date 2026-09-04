#pragma once
//
// boundary_training_log.hpp
//
// Accumulates operator-confirmed landmark boundaries as labeled training data
// for the eventual boundary-detection CNN (Section 9.10).
//
// Each record captures a fixed +/-100 ms segment of signal around a landmark,
// the operator-confirmed sample index within that segment (the label), which
// fit-and-select model the Phase 1 routine chose plus its residual, the
// individual the beat came from, and a little clinical metadata. The confirmed
// index comes from the operator correction in the Phase 1 B2 focus mode.
//
// Records are appended to disk immediately as CSV rows, so they survive
// crashes and accumulate across sessions and individuals toward the target of
// >= 10,000 labeled boundaries per landmark type. One file per landmark type
// (keeps the set trivially splittable); the fixed-length +/-100 ms segment is
// written as trailing sample columns s0..sN.
//
// This component only reuses anchor_fit.hpp (for FitType); it has no other
// dependency on the marking/template pipeline, so it can be called from
// wherever the operator confirmation happens.
//

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include "template_anchoring\anchor_fit.hpp"   // anchor_fit::FitType

namespace boundary_training {

    // The five landmark types this log accumulates (spec goal). Kept separate
    // from AnchorType so this component stays independent of the marking enum;
    // the caller maps whatever it uses onto these.
    enum class Landmark { Q_ONSET, S_END, J_POINT, P_ONSET, T_OFFSET };

    inline const char* landmark_name(Landmark l) {
        switch (l) {
        case Landmark::Q_ONSET:  return "Q_ONSET";
        case Landmark::S_END:    return "S_END";
        case Landmark::J_POINT:  return "J_POINT";
        case Landmark::P_ONSET:  return "P_ONSET";
        case Landmark::T_OFFSET: return "T_OFFSET";
        }
        return "UNKNOWN";
    }

    inline const char* fittype_name(anchor_fit::FitType t) {
        switch (t) {
        case anchor_fit::FitType::LINEAR:     return "LINEAR";
        case anchor_fit::FitType::SIGMOID:    return "SIGMOID";
        case anchor_fit::FitType::FRACTIONAL: return "FRACTIONAL";
        case anchor_fit::FitType::FLAT:       return "FLAT";
        }
        return "UNKNOWN";
    }

    // One labeled boundary. `segment` is +/-100 ms around the landmark
    // (200 samples at 1000 Hz). Two indices INTO `segment` are recorded:
    //   autoDetect -- the Phase 1 fit-and-select algorithm's ORIGINAL guess
    //                 (the seeded position, never a prior manual edit);
    //   expertMark -- the operator-confirmed position from B2 focus mode.
    // expertMark is the training LABEL; autoDetect is kept so the auto
    // detector's error (expertMark - autoDetect) can be measured, and so
    // large-correction examples can be filtered/weighted during training.
    struct BoundaryTrainingRecord {
        std::vector<double> segment; // +/-100 ms (200 samples at 1000 Hz)
        int confirmedIndex; // label
        anchor_fit::FitType fitType; // which fit-and-select model won (incl. FLAT)
        double fitRSS;
        std::string individualID;
        bool bbb; double heartRate, qrsDurationMs; // clinical metadata
    };

    // A log destination: a directory holding ONE boundary log file for the
    // whole run. landmark and anchor are columns in each row (individualID
    // doubles as the source filename), so all landmarks/anchors/templates
    // share the single file.
    struct BoundaryTrainingLog {
        std::string dir;

        explicit BoundaryTrainingLog(std::string directory = std::string())
            : dir(std::move(directory)) {
            if (!dir.empty()) {
                std::error_code ec;
                std::filesystem::create_directories(dir, ec);
            }
        }

        std::string path() const {
            return dir.empty() ? std::string("boundary_training.csv")
                : dir + "/boundary_training.csv";
        }
    };

    namespace detail {
        // CSV field escaping: quote the field and double internal quotes if it
        // contains a comma, quote, or newline (RFC-4180). individualID is the
        // only free-text field; everything else is numeric/enumerated.
        inline std::string csv_escape(const std::string& s) {
            const bool needs = s.find_first_of(",\"\n\r") != std::string::npos;
            if (!needs) return s;
            std::string out = "\"";
            for (char c : s) { if (c == '"') out += "\"\""; else out += c; }
            out += "\"";
            return out;
        }
    }

    // Append one boundary to the single run log as a CSV row. `anchor` (the
    // template's alignment pass) and `landmark` are columns, so one file holds
    // every template/anchor/landmark. (No filename column -- it equals
    // individualID.) confirmedIndex is left BLANK when the operator didn't
    // touch the marker (rec.confirmedIndex < 0); the row is still written.
    // Writes the header row first if the file is new. Columns:
    //   anchor,landmark,individualID,confirmedIndex,fitType,fitRSS,bbb,
    //   heartRate,qrsDurationMs,nSegment,s0,s1,...,s{n-1}
    // Returns false if the file couldn't be opened or the write failed.
    inline bool logBoundary(BoundaryTrainingLog& log,
        const std::string& anchor,
        Landmark landmark,
        const BoundaryTrainingRecord& rec) {
        const std::string path = log.path();

        bool needHeader = true;
        {
            std::error_code ec;
            if (std::filesystem::exists(path, ec) &&
                std::filesystem::file_size(path, ec) > 0) needHeader = false;
        }

        std::ofstream f(path, std::ios::app);
        if (!f) return false;

        if (needHeader) {
            f << "anchor,landmark,individualID,confirmedIndex,fitType,"
                "fitRSS,bbb,heartRate,qrsDurationMs,nSegment";
            for (size_t i = 0; i < rec.segment.size(); ++i) f << ",s" << i;
            f << "\n";
        }

        f << std::setprecision(9);
        f << detail::csv_escape(anchor) << ','
            << landmark_name(landmark) << ','
            << detail::csv_escape(rec.individualID) << ',';
        // confirmedIndex is blank when the operator didn't touch the marker
        // (sentinel < 0 => empty CSV field); the row is still logged.
        if (rec.confirmedIndex >= 0) f << rec.confirmedIndex;
        f << ','
            << fittype_name(rec.fitType) << ','
            << rec.fitRSS << ','
            << (rec.bbb ? 1 : 0) << ','
            << rec.heartRate << ','
            << rec.qrsDurationMs << ','
            << rec.segment.size();
        for (double v : rec.segment) f << ',' << v;
        f << "\n";

        return f.good();
    }

}   // namespace boundary_training