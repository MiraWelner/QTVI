#pragma once
/*
* @file   annotation_types.hpp
* @brief  Single source of truth for the marking/annotation types. Text, color,
*         binary-export code, post-beat tag, and detection-suppression flag all
*         live here so they can't drift apart across the GUI, the exporter, and
*         the peak detector. The `label` strings are the canonical combo-box /
*         marking_type strings -- match them exactly everywhere.
*
*         QT-FREE, AND THAT IS THE POINT. Every field of AnnotationType is a
*         const char*, an int or a bool -- the table never needed Qt. Only five
*         LOOKUP HELPERS did: colorFor() returns a QColor and four overloads took
*         a QString. Those five now live in chart_utils.hpp, which is already the
*         Qt-side color helper header and is already included by all three
*         translation units that call them.
*
*         What that buys: the template-generation side (template_bank.hpp,
*         bin_pipeline.hpp, seed_pool.hpp, nsvt_detect.hpp) is Qt-free and can
*         now include THIS header directly. It used to mirror the class codes as
*         literals because it couldn't, and annotation_code_check.hpp existed to
*         catch the two copies drifting apart at startup -- a check that was
*         never called from anywhere, so the drift went unguarded for as long as
*         the mirror existed. With one copy there is nothing to check, and
*         code_for_label() below resolves at compile time so a renamed label is a
*         build failure instead of a silent 0.
*/

#include <array>
#include <cstdint>
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

    // ROW ORDER IS THE OPERATOR'S DROPDOWN ORDER. The combo box is populated
    // from this array in order (gui_handler.cpp), so moving a row here moves it
    // in the dropdown. The order below reproduces exactly what the .ui file used
    // to hardcode, Invert/Noninvert included at position 10 -- operators have
    // months of muscle memory in that list and it must not shift because the
    // labels moved house.
    //
    // CODE ORDER IS THEREFORE NOT ROW ORDER, and that is fine. Invert/Noninvert
    // is code 13 sitting tenth. Nothing derives a code from a position: every
    // consumer either looks up by label (code_for_label, find) or scans for a
    // matching code (bin_chunk_loader.cpp:59). The static_asserts below check
    // that codes are unique and non-zero, which is the property that actually
    // matters, rather than that they ascend.
    inline constexpr std::array<AnnotationType, 13> noise_types = { {
            //  label               code  r    g    b    a   suppress  post   param  invert inclThr
        { "1) R Peak Noise",   1,  235, 220, 0,   60, true,  false, false },  // yellow
        { "2) Minor Noise",    2,  150, 80,  30,  60, false, false, false, false, true },  // brown; no effect on detection/thresholds
        { "3) Blank.+Thresh.", 3,  128, 128, 128, 60, false, false, true  },  // gray
        { "4) PVC",            4,  130, 190, 30,  60, false, true,  false },  // lime
        { "5) PAC",            5,  245, 130, 0,   60, false, true,  false },  // orange
        { "6) Cond. Delay",    6,  70,  60,  200, 60, false, false, false },  // indigo
        { "7) AF",             7,  230, 25,  25,  60, false, true,  false },  // red
        { "8) SVT",            8,  0,   160, 160, 60, false, true,  false },  // teal
        { "9) VT",             9,  30,  170, 70,  60, false, true,  false },  // green
        { "Invert/Noninvert",  13, 180, 100, 0,   60, false, false, false, true },  // orange; override, not stored
        { "Benign Arr.",       10, 200, 40,  170, 60, false, false, false },  // magenta
        { "Sig. Arr.",         11, 0,   120, 210, 60, false, false, false },  // blue
        { "Other",             12, 110, 120, 140, 60, false, false, false }   // slate
        } };

    // Char-by-char label comparison. Constexpr, which is what lets everything
    // below happen at compile time -- QString::operator== is a runtime function
    // and was the reason the old checks had to run at startup.
    constexpr bool labels_equal(const char* a, const char* b) {
        while (*a && *a == *b) { ++a; ++b; }
        return *a == *b;
    }

    // ---------------------------------------------------------------------
    // Table-wide invariants, checked at compile time
    // ---------------------------------------------------------------------
    //
    // Every consumer of this table assumes all of the following, and until now
    // none of them was stated anywhere. They are cheap to check and each one
    // fails quietly if violated, which is the combination that produces an
    // archive that looks correct and is not:
    //
    //   DUPLICATE CODE   -> two types export as the same number, and every
    //                       reader that scans for a code (bin_chunk_loader.cpp)
    //                       silently takes whichever comes first. Old exports
    //                       become ambiguous with no way to tell.
    //   ZERO CODE        -> 0 is the unlabeled sentinel downstream
    //                       (tbank::kUnlabeled, morphology_csv's label_code),
    //                       so that type's spans read as "no annotation".
    //   DUPLICATE LABEL  -> find() returns the first match, so one of the two is
    //                       unreachable and its colour and flags never apply.
    //   EMPTY LABEL      -> matches nothing and cannot be selected, but appears
    //                       in the dropdown as a blank row.
    //   ZERO ALPHA       -> the highlight is drawn fully transparent, so marked
    //                       spans are invisible while the marking is recorded.
    //   RGB OUT OF RANGE -> QColor silently clamps, so the swatch and the
    //                       highlight disagree with the table.

    constexpr bool codes_unique() {
        for (std::size_t i = 0; i < noise_types.size(); ++i)
            for (std::size_t j = i + 1; j < noise_types.size(); ++j)
                if (noise_types[i].code == noise_types[j].code) return false;
        return true;
    }
    constexpr bool codes_nonzero() {
        for (const auto& t : noise_types) if (t.code == 0) return false;
        return true;
    }
    constexpr bool labels_unique() {
        for (std::size_t i = 0; i < noise_types.size(); ++i)
            for (std::size_t j = i + 1; j < noise_types.size(); ++j)
                if (labels_equal(noise_types[i].label, noise_types[j].label))
                    return false;
        return true;
    }
    constexpr bool labels_nonempty() {
        for (const auto& t : noise_types)
            if (t.label == nullptr || t.label[0] == '\0') return false;
        return true;
    }
    constexpr bool colors_in_range() {
        for (const auto& t : noise_types) {
            if (t.r < 0 || t.r > 255) return false;
            if (t.g < 0 || t.g > 255) return false;
            if (t.b < 0 || t.b > 255) return false;
            if (t.a <= 0 || t.a > 255) return false;   // 0 alpha = invisible
        }
        return true;
    }

    static_assert(codes_unique(), "annotation_types: two rows share an export code");
    static_assert(codes_nonzero(), "annotation_types: a row has code 0, which is the unlabeled sentinel");
    static_assert(labels_unique(), "annotation_types: two rows share a label; find() would reach only the first");
    static_assert(labels_nonempty(), "annotation_types: a row has an empty label");
    static_assert(colors_in_range(), "annotation_types: a row has an out-of-range or fully transparent color");

    // ---------------------------------------------------------------------
    // Compile-time lookup by label
    // ---------------------------------------------------------------------
    //
    // annotation_code_check.hpp's stated reason for being a RUNTIME check was
    // that "find() returns a pointer into a constexpr array by string
    // comparison, which is not a constant expression under C++17." The
    // comparison was the non-constexpr part, not the array -- QString's
    // operator== is a runtime function. A char-by-char compare over string
    // literals IS a constant expression, so the lookup happens at compile time
    // and the check becomes a static_assert at the point of use.

    // The code for a label, or 0 if no entry carries it. 0 is also the
    // "unlabeled" sentinel downstream, so a miss is indistinguishable from a
    // legitimate no-class at runtime -- which is exactly how the old drift would
    // have stayed invisible. Callers that require the label to exist must
    // static_assert on the result being non-zero; see template_bank.hpp.
    constexpr uint8_t code_for_label(const char* label) {
        for (const auto& t : noise_types)
            if (labels_equal(t.label, label)) return static_cast<uint8_t>(t.code);
        return 0;
    }

    // ---------------------------------------------------------------------
    // The two types identified by BEHAVIOUR rather than by name
    // ---------------------------------------------------------------------
    //
    // paramEdit and invertEdit each mark exactly one row: dragging in such a
    // span edits a threshold or flips inversion instead of storing an
    // annotation. getAllMarkings() has to write those spans out under the right
    // marking_type, and it used to do so by typing the label
    // ("3) Blank.+Thresh." and "Invert/Noninvert") directly into the function.
    //
    // Those were the last two copies of a label string outside this table, and
    // the most likely to be missed in a rename: they sit mid-function rather
    // than in a list next to their siblings, so nothing about editing the table
    // would draw anyone's eye to them. Worse, they were silent on failure --
    // find() returns nullptr for an unknown label and every helper degrades
    // quietly, so a retitled row would have exported those spans as code 0.
    //
    // Looking them up BY FLAG rather than by label removes the string from the
    // call site entirely, instead of relocating it. The flag is also the more
    // honest key: what getAllMarkings() actually needs is "the type that means
    // parameter edit", which is what the flag says and what the label only
    // implies.
    constexpr int count_with_param_edit() {
        int n = 0;
        for (const auto& t : noise_types) if (t.paramEdit) ++n;
        return n;
    }
    constexpr int count_with_invert_edit() {
        int n = 0;
        for (const auto& t : noise_types) if (t.invertEdit) ++n;
        return n;
    }

    // Exactly one of each, and the code assumes it. Two rows flagged paramEdit
    // would make the choice between them arbitrary and order-dependent; zero
    // would make the accessor below return an empty label that exports as code
    // 0. Both are build failures rather than runtime surprises.
    static_assert(count_with_param_edit() == 1,
        "annotation_types: exactly one row must set paramEdit");
    static_assert(count_with_invert_edit() == 1,
        "annotation_types: exactly one row must set invertEdit");

    constexpr const char* param_edit_label() {
        for (const auto& t : noise_types) if (t.paramEdit) return t.label;
        return "";
    }
    constexpr const char* invert_edit_label() {
        for (const auto& t : noise_types) if (t.invertEdit) return t.label;
        return "";
    }

    inline constexpr const char* kParamEditLabel = param_edit_label();
    inline constexpr const char* kInvertEditLabel = invert_edit_label();

    inline constexpr uint8_t kParamEditCode = code_for_label(kParamEditLabel);
    inline constexpr uint8_t kInvertEditCode = code_for_label(kInvertEditLabel);

    // ---------------------------------------------------------------------
    // Runtime lookups
    // ---------------------------------------------------------------------

    inline const AnnotationType* find(const char* label) {
        for (const auto& t : noise_types)
            if (labels_equal(t.label, label)) return &t;
        return nullptr;
    }
    inline const AnnotationType* find(const std::string& label) {
        return find(label.c_str());
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

    // Whether this span's signal is KEPT in the reference-window threshold /
    // mean-R-R statistics. When true, the span has no effect on detection or on
    // the placement of following R peaks (e.g. "Minor Noise"). When false (the
    // default for every other type), the span is excluded from those stats, as
    // it historically always was.
    inline bool includeInThreshold(const std::string& label) {
        const auto* t = find(label);
        return t && t->includeInThreshold;
    }

}  // namespace annotation_types