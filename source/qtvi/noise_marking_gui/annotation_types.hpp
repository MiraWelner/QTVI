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

    };

    inline constexpr std::array<AnnotationType, 10> noise_types = { {
        { "1) R Peak Noise",  1, 255, 255, 0,   60, true,  false },
        { "2) Minor Noise",    2, 157, 60,  0,   60, false, false },
        { "3) Cond. Delay",    3, 128, 0,   128, 60, false, false },
        { "4) AF",             4, 255, 0,   0,   60, false, true  },
        { "5) SVT",            5, 0,   0,   255, 60, false, true  },
        { "6) VT",             6, 0,   255, 0,   60, false, true  },
        { "7) PVC",            7, 128, 255, 0,   60, false, true  },
        { "8) PAC",            8, 255, 128, 0,   60, false, true  },
        { "9) Benign Arr.",    9, 255, 128, 255, 60, false, false },
        { "10) Sig. Arr.",     10, 0,  255, 255, 60, false, false }
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

}  // namespace annotation_types