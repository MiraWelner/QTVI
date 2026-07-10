#pragma once
/*
* @file   annotation_types.hpp
* @brief  Single source of truth for the marking/annotation types. Text, color,
*         binary-export code, post-beat tag, and detection-suppression flag all
*         live here so they can't drift apart across the GUI, the exporter, and
*         the peak detector. The `label` strings are the canonical combo-box /
*         marking_type strings -- match them exactly everywhere.
*/

#include <QString>
#include <QColor>
#include <array>
#include <string>

namespace annotation_types {

    struct AnnotationType {
        const char* label;            // exact marking_type string, e.g. "1) Noise/Art."
        int         code;             // numeric code for the binary export
        int r, g, b, a;               // highlight color (RGBA, a = alpha)
        bool        suppressesDetection;  // true => no R peaks detected in the span
        bool        postEligible;     // beat AFTER one of these gets this code in the post column
        bool        paramEdit;        // true => dragging edits threshold/blanking, not a stored annotation
        bool        invertEdit;       // true => dragging flips inversion in the span, not a stored annotation
        bool        includeInThreshold; // true => this span's signal IS kept in the
        //   reference-window threshold / mean-R-R stats
        //   (i.e. NOT excluded), so it has no effect on
        //   detection or the placement of following R
        //   peaks. false => excluded from those stats
        //   (the historical behaviour for every type).
    };

    inline constexpr std::array<AnnotationType, 13> noise_types = { {
            //  label               code  r    g    b    a   suppress  post   param  invert inclThr
        { "1) R Peak Noise",   1,  235, 220, 0,   60, true,  false, false },  // yellow
        { "2) Minor Noise",    2,  150, 80,  30,  60, false, false, false, false, true },  // brown; no effect on detection/thresholds
        { "3) Blank.+Thresh.", 3,  128, 128, 128, 60, false, false, true  },  // gray
        { "4) Cond. Delay",    4,  70,  60,  200, 60, false, false, false },  // indigo
        { "5) AF",             5,  230, 25,  25,  60, false, true,  false },  // red
        { "6) SVT",            6,  0,   160, 160, 60, false, true,  false },  // teal
        { "7) VT",             7,  30,  170, 70,  60, false, true,  false },  // green
        { "8) PVC",            8,  130, 190, 30,  60, false, true,  false },  // lime
        { "9) PAC",            9,  245, 130, 0,   60, false, true,  false },  // orange
        { "Benign Arr.",       10, 200, 40,  170, 60, false, false, false },  // magenta
        { "Sig. Arr.",         11, 0,   120, 210, 60, false, false, false },  // blue
        { "Other",             12, 110, 120, 140, 60, false, false, false },  // slate
        { "Invert/Noninvert",  13, 180, 100, 0,   60, false, false, false, true }  // orange; override, not stored
        } };

    inline const AnnotationType* find(const QString& label) {
        for (const auto& t : noise_types)
            if (label == QLatin1String(t.label)) return &t;
        return nullptr;
    }
    inline const AnnotationType* find(const std::string& label) {
        return find(QString::fromStdString(label));
    }

    // Highlight color for a type; default translucent black if unknown
    // (matches the old updateNoiseHighlights fallback).
    inline QColor colorFor(const QString& label) {
        if (const auto* t = find(label)) return QColor(t->r, t->g, t->b, t->a);
        return QColor(0, 0, 0, 100);
    }
    // Numeric export code; 0 if unknown (matches the old typeMap miss).
    inline double codeFor(const std::string& label) {
        if (const auto* t = find(label)) return static_cast<double>(t->code);
        return 0.0;
    }

    inline int markCode(const std::string& label) {
        const auto* t = find(label);
        return t ? t->code : 0;
    }
    inline int postCode(const std::string& label) {
        const auto* t = find(label);
        return (t && t->postEligible) ? t->code : 0;
    }

    // Whether peaks should be suppressed inside this type's spans (noise/artifact).
    inline bool suppressesDetection(const std::string& label) {
        const auto* t = find(label);
        return t && t->suppressesDetection;
    }

    inline bool isParamEdit(const QString& label) {
        const auto* t = find(label);
        return t && t->paramEdit;
    }

    inline bool isInvertEdit(const QString& label) {
        const auto* t = find(label);
        return t && t->invertEdit;
    }

    // Whether this span's signal is KEPT in the reference-window threshold /
    // mean-R-R statistics. When true, the span has no effect on detection or
    // on the placement of following R peaks (e.g. "Minor Noise"). When false
    // (the default for every other type), the span is excluded from those
    // stats, as it historically always was.
    inline bool includeInThreshold(const std::string& label) {
        const auto* t = find(label);
        return t && t->includeInThreshold;
    }
    inline bool includeInThreshold(const QString& label) {
        const auto* t = find(label);
        return t && t->includeInThreshold;
    }

}  // namespace annotation_types