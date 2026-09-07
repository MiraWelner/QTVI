#pragma once
//
// Adapter layer: reads the template_io binary format and projects it into
// the TemplateBin shape the viewer expects. Also writes/reads the
// separate template_markings.bin produced by the viewer.
//
// Marker set per bin -- THE HUMAN-MOVABLE BARS ONLY:
//   ECG (per lead, per slot, per anchor):  P-onset, Q-onset, S-end, T-end
//   PPG (shared):                          Onset, Dicrotic notch, End
//   Arterial (ABP / ART / ART_PULM):       Onset, Peak, Dicrotic, Peak2, End
//
// Nothing else is persisted. Every reactive glyph (ECG P-peak and T-peak,
// PPG t50 / t80 / t80_rise / pw80 / peak2) and every detector output (the
// *_auto fields) is recomputed from the bars or re-detected on load, because
// a stored copy is a second answer that outlives the bars it came from.
//
// ALL MARKER POSITIONS ARE DOUBLES, in the file and in memory. The detectors
// and the refine/fit stages produce sub-sample positions; rounding them into
// int on the way into a bar threw that away and quantised every millisecond
// column to the sample period. -1 remains the "unmarked / not applicable"
// sentinel; test with `< 0`, not `== -1`.
//
// Structural integers in the file (slot counts, anchor counts, anchor tags,
// bin index, issue flags) stay integers. Only positions are widened.
//
// IQR vectors (ecgTemplate_raw_iqr per channel, ppg_template_iqr) come
// straight from the template file. Empty std => the widget renders the
// trace without a gray band.
//

#include <array>
#include <vector>
#include <string>
#include <map>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <utility>
#include <fstream>
#include <ostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <cstring>

#include "template_generation\template_io.hpp"
#include "template_marking_gui\feature_marks.hpp"
#include "template_morphology_grouping\template_bank.hpp"
#include "template_anchoring\anchor_view.hpp"


// Which half of the markings CSV a call emits. Declared HERE, above every
// function that names it in a parameter list -- it used to sit further down,
// between the two writers, so the first one referenced it before declaration.
enum class MarkingsCsvSection { EcgOnly, PulseOnly };


// ---------------------------------------------------------------------------
// In-memory model used by the viewer.
// ---------------------------------------------------------------------------

struct ChannelTemplateData {
    std::vector<double> ecgTemplate_raw;
    std::vector<double> ecg_template_raw_iqr;
    std::vector<double> ecgTemplate_squared, ecgTemplate_absval;    // unused by viewer
    double alignment_point_raw = 0, alignment_point_squared = 0;
    double alignment_point_absval = 0;
    // True R column in the template (from alignment). Used directly as the R
    // fiducial; replaces the old avg_r_expand positioning constant.
    //
    // STAYS INT: this is a template_io file field, not a marker. It is the
    // column the aligner placed R on, and frameShift differences two of them
    // to get an integer column offset.
    int r_col_raw = -1;
    int r_col_squared = -1;
    int r_col_absval = -1;

    // Median RR of this channel's beats, in samples. THE PLOT WIDTH, and only
    // that: the beat matrix behind the template is framed on the bin's LONGEST
    // RR so that no beat is ever clipped, which leaves the array several times
    // wider than a normal beat whenever the bin contains a pause or a missed R
    // detection. The panel used to take its x-axis from that array length, so a
    // 0.9 s beat was drawn into a sixth of the frame. -1 means unknown, and the
    // panel then falls back to the array length.
    int median_rr_samples = -1;
};

// One bank slot's average for one alignment.
//
// AT NAMESPACE SCOPE, NOT NESTED IN TemplateBin, and that is not a style
// choice. As a nested type, every use had to be spelled
// TemplateBin::AnchoredBankSlot, and MSVC failed to parse those declarations --
// both `map<int, array<vector<...>,3>>` inside a member function and
// `vector<TemplateBin::AnchoredBankSlot> slots(n)` in the loader. Once a
// declarator fails, the following `slots[sl]` is read as an ATTRIBUTE or a
// LAMBDA INTRODUCER, which is where "'sl': attribute not found" and
// "expected a '{' introducing a lambda body" came from, plus a syntax error on
// every `.member` after it. Flat name, no nesting, no lookup.
struct AnchoredBankSlot {
    std::vector<double> tmpl;
    std::vector<double> tmpl_iqr;
    uint32_t n_members = 0;
};

struct TemplateBin {
    // Section 4.6 template bank per ECG channel, index 0..2 == CH1..CH3.
    // Slot 0 of each is the sinus seed and corresponds to that channel's
    // ecgTemplate, except that it excludes ectopy. Channels are allowed to
    // disagree on template count -- a morphology separable on one lead may not
    // be on another -- so these are NOT parallel across the array.
    std::array<tbank::TemplateBank, 3> ecg_bank;

    // ---- AND THE PULSE BANK ---------------------------------------------
    //
    // template_io::BinTemplates has carried ppg_bank since v4 and this struct
    // did not, so `dst.ecg_bank = src.ecg_bank` in fromTemplateFile() loaded the
    // pulse partition off disk and then dropped it. The consequences were all
    // downstream and all silent: the viewer had no pulse bank to label, so a
    // class confirmation could not reach the PPG cohort of the morphology it
    // confirmed, and a right-click had nowhere to record a pulse verdict.
    //
    // Slot i is group i, the same group ECG slot i is -- the projection walks
    // the joint groups in order for every channel, so the four banks are
    // parallel by construction. That is what makes labelling and marking by
    // slot correct across all of them.
    tbank::TemplateBank ppg_bank;

    uint64_t index = 0;
    std::vector<std::pair<uint64_t, uint64_t>> ppg_bin_indexs;
    std::vector<std::pair<uint64_t, uint64_t>> ecg_bin_indexs;
    bool bad_segment = false;
    // Slice counts (post drop-rules) that fed each channel's raw-method
    // median, plus PPG. All driven by ch1.raw R-pairs under Patch B, so
    // they normally read equal; per-channel storage lets a future filter
    // drop them per-channel. 0 = channel absent / no beats.
    uint64_t ch1_n_beats_raw = 0;
    uint64_t ch2_n_beats_raw = 0;
    uint64_t ch3_n_beats_raw = 0;
    uint64_t ppg_n_beats = 0;

    ChannelTemplateData ch1, ch2, ch3;

    // ====================================================================
    // ALL FOUR ALIGNMENTS, LOADED AT ONCE
    //
    // ch1/ch2/ch3 above stay the R-aligned base: the grid draws them, every
    // existing reference to them keeps meaning what it meant, and R is the
    // frame every other alignment's columns are translated into. The other
    // three live here, keyed by (int)AnchorType.
    //
    // readTemplateInfoBin used to project ONE anchor block into ch1..3, which
    // is why seeing a second alignment meant regenerating the templates file
    // and reopening the window. Four sequential screens were a consequence of
    // this one field being a scalar.
    // ====================================================================
    std::map<int, std::array<ChannelTemplateData, 3>> anchored;

    // Glyph positions for ONE alignment. Glyphs are measurements, not
    // judgements, so every alignment detects all of them on its own average
    // and all four are reported -- <landmark>_auto_P through _auto_T are four
    // measurements of one landmark on four waveforms. The flat *_auto_ch
    // fields are R's copy, kept flat so the grid's glyph drawing and every
    // existing consumer are untouched; this map carries the other three.
    struct AnchorAuto {
        double p_begin[3] = { -1, -1, -1 };
        double p_peak[3] = { -1, -1, -1 };
        double q_begin[3] = { -1, -1, -1 };
        double q_peak[3] = { -1, -1, -1 };
        double r_peak[3] = { -1, -1, -1 };
        double s_end[3] = { -1, -1, -1 };
        double t_end[3] = { -1, -1, -1 };
        bool   q_begin_found[3] = { false, false, false };
    };
    std::map<int, AnchorAuto> auto_by_anchor;

    // ---- PER-ANCHOR BANK SLOT AVERAGES --------------------------------
    //
    // anchored above holds one aligned average per (anchor, lead) -- the bin's
    // whole-channel template, which is what slot 0 draws. The BANK slots (the
    // _A / _B columns) had no aligned variant, so their focus panel showed the
    // slot's own R-aligned average under whatever "[X-aligned]" header the
    // clicked bar produced. The switch did nothing and the label was wrong.
    //
    // Keyed [anchor tag][lead][slot], mirroring template_io's bank_anchors.
    // Produced offline in alignTemplatesFromCache as column-wise reductions
    // over the SAME aligned beat matrix the whole-channel average comes from,
    // so a slot and its bin share one frame.
    // FLAT KEY: anchorTag * 4 + lead. Deliberately not
    // map<int, array<vector<...>,3>>: MSVC could not parse a declaration of
    // that nested type inside a member function, and once the declarator
    // failed it read the following `(...)[lead]` as a LAMBDA INTRODUCER --
    // producing "expected a '{' introducing a lambda body" and a cascade of
    // syntax errors on the lines below. One level of nesting, no arrays.
    std::map<int, std::vector<AnchoredBankSlot> > anchored_bank;
    static int bankSlotKey(int lead, AnchorType a) {
        return static_cast<int>(a) * 4 + lead;
    }

    // This slot's aligned average, or nullptr when the file has none -- a
    // pre-v6 templates file, a bin the alignment skipped, or a slot with no
    // members. Callers fall back to the slot's own unaligned average AND say so
    // in the panel header: a header naming an alignment the data is not in is
    // the defect this exists to fix, so silently falling back would reinstate
    // it one level down.
    const AnchoredBankSlot* bankSlotFor(int lead, int slot, AnchorType a) const
    {
        if (lead < 0 || lead > 2) return nullptr;
        if (slot < 0) return nullptr;
        auto it = anchored_bank.find(bankSlotKey(lead, a));
        if (it == anchored_bank.end()) return nullptr;
        if (static_cast<size_t>(slot) >= it->second.size()) return nullptr;
        if (it->second[slot].tmpl.empty()) return nullptr;
        return &it->second[slot];
    }

    // This alignment's template for one lead. Falls back to the R base when
    // the anchor block is absent -- a templates file built before this change
    // has none, and every bar's close-up then shows the R template, which is
    // the old single-pass behaviour rather than a failure.
    // STRICT: this alignment's template, or nullptr. NO FALLBACK.
    //
    // chFor below substitutes the R base when the anchor block is absent, which
    // is right for the CSV writers -- an R-only file must still produce columns
    // -- and WRONG for the focus panel, where it silently showed the R-aligned
    // average under a "[P-aligned]" header. A view that claims an alignment it
    // is not showing is worse than an empty view, because there is no way to
    // tell the two apart by eye.
    //
    // So the panel uses this one and draws nothing when it returns nullptr.
    const ChannelTemplateData* chForStrict(int lead, AnchorType a) const {
        if (lead < 0 || lead > 2) return nullptr;
        if (a == AnchorType::R_PEAK) {
            const ChannelTemplateData* base[3] = { &ch1, &ch2, &ch3 };
            return base[lead]->ecgTemplate_raw.empty() ? nullptr : base[lead];
        }
        auto it = anchored.find(static_cast<int>(a));
        if (it == anchored.end()) return nullptr;
        if (it->second[lead].ecgTemplate_raw.empty()) return nullptr;
        return &it->second[lead];
    }

    const ChannelTemplateData& chFor(int lead, AnchorType a) const {
        const ChannelTemplateData* base[3] = { &ch1, &ch2, &ch3 };
        if (lead < 0 || lead > 2) return ch1;
        if (a == AnchorType::R_PEAK) return *base[lead];
        auto it = anchored.find(static_cast<int>(a));
        if (it == anchored.end()) return *base[lead];
        const ChannelTemplateData& c = it->second[lead];
        return c.ecgTemplate_raw.empty() ? *base[lead] : c;
    }

    // This alignment's glyphs. Absent anchor -> the flat (R) fields, which is
    // also how loadSubject captures each pass's result on the way past.
    AnchorAuto autoFor(AnchorType a) const {
        auto it = auto_by_anchor.find(static_cast<int>(a));
        if (it != auto_by_anchor.end()) return it->second;
        AnchorAuto out;
        for (int c = 0; c < 3; ++c) {
            out.p_begin[c] = p_begin_auto_ch[c];
            out.p_peak[c] = p_peak_auto_ch[c];
            out.q_begin[c] = q_begin_auto_ch[c];
            out.q_peak[c] = q_peak_auto_ch[c];
            out.r_peak[c] = r_peak_auto_ch[c];
            out.s_end[c] = s_end_auto_ch[c];
            out.t_end[c] = t_end_auto_ch[c];
            out.q_begin_found[c] = q_begin_found_auto_ch[c];
        }
        return out;
    }

    // Samples to ADD to a column measured in `from`'s frame to express it in
    // `to`'s. Both alignments put the same beat's R on their own r_col, so the
    // difference of the two r_cols is the whole conversion. Normally 0 or a
    // sample or two: alignTemplatesFromCache aligns each anchor against the
    // median-length snippet, which does not move, so R lands near the base
    // column -- near, not exactly, which is why this is a function and not an
    // assumption. 0 when either r_col is unknown.
    //
    // STAYS INT. r_col is a file field and their difference is a whole-sample
    // column offset; adding it to a double position is exact.
    int frameShift(int lead, AnchorType from, AnchorType to) const {
        const int a = chFor(lead, from).r_col_raw;
        const int b = chFor(lead, to).r_col_raw;
        return (a < 0 || b < 0) ? 0 : (b - a);
    }

    // ---- THE BAR SET THE OPERATOR ACTUALLY PLACED ------------------------
    //
    // Each bar lives in its owning alignment's set, so no single
    // markers_by_anchor entry holds a whole beat any more. Anything needing
    // P-onset AND Q-onset AND S-end AND T-end together -- QRS duration, QT,
    // the reactive P and T peaks, the global interval lines, the VCG rows --
    // assembles them here, each translated into `frame`'s columns on the way.
    //
    // BARS ONLY, and that is now the whole of BankMarkerSet: the reactive
    // P-peak is no longer stored anywhere, so there is nothing else to pull.
    tbank::BankMarkerSet userMarks(int lead, int slot,
        AnchorType frame = AnchorType::R_PEAK) const
    {
        tbank::BankMarkerSet out;   // all -1
        // Pointer-to-member is DOUBLE now, matching BankMarkerSet's fields.
        auto pull = [&](int marker, double tbank::BankMarkerSet::* field) {
            const AnchorType owner = anchor_view::anchorFor(marker);
            const double v = slotMarks(lead, slot, owner).*field;
            if (v < 0.0) return;
            out.*field = v + frameShift(lead, owner, frame);
            };
        pull(anchor_view::kPBegin, &tbank::BankMarkerSet::p_begin);
        pull(anchor_view::kQBegin, &tbank::BankMarkerSet::q_begin);
        pull(anchor_view::kSEnd, &tbank::BankMarkerSet::s_end);
        pull(anchor_view::kTEnd, &tbank::BankMarkerSet::t_end);
        return out;
    }

    // ---- THE GLYPH COUNTERPART OF userMarks() -----------------------------
    //
    // Same idea as userMarks(), but for the frozen auto-detect dot instead of
    // the draggable bar: each landmark is read from ITS OWNING alignment via
    // autoFor(), then frame-shifted into `frame`'s columns exactly the way
    // userMarks() shifts the bar. Without this, a glyph reader that reaches
    // into the flat *_auto_ch fields directly gets whatever alignment ran
    // LAST in the seeding loop (R_PEAK) regardless of which alignment a given
    // landmark actually belongs to -- so a T-end glyph would show R_PEAK's own
    // T-end detection, unshifted, while the T-end bar correctly shows J_POINT's
    // detection translated into R's frame. Two different measurements, not
    // just two different rounding paths, and no additive shift reconciles
    // them after the fact.
    //
    // NO ROUNDING ON THE WAY IN any more. AnchorAuto is double and so is
    // BankMarkerSet, so the detector's sub-sample position survives the
    // frame shift intact; the lround these lines used to do was the last
    // place a glyph got quantised.
    tbank::BankMarkerSet autoMarks(int lead, AnchorType frame = AnchorType::R_PEAK) const
    {
        tbank::BankMarkerSet out;   // all -1
        // AnchorAuto stores each landmark as a double[3] (per lead), so this
        // is spelled out per field rather than one generic lambda over a
        // pointer-to-member -- you can't take a pointer-to-member to an
        // array element the way userMarks() does for BankMarkerSet's scalar
        // fields.
        {
            const AnchorType owner = anchor_view::anchorFor(anchor_view::kPBegin);
            const double v = autoFor(owner).p_begin[lead];
            if (v >= 0.0) out.p_begin = v + frameShift(lead, owner, frame);
        }
        {
            const AnchorType owner = anchor_view::anchorFor(anchor_view::kQBegin);
            const double v = autoFor(owner).q_begin[lead];
            if (v >= 0.0) out.q_begin = v + frameShift(lead, owner, frame);
        }
        {
            const AnchorType owner = anchor_view::anchorFor(anchor_view::kSEnd);
            const double v = autoFor(owner).s_end[lead];
            if (v >= 0.0) out.s_end = v + frameShift(lead, owner, frame);
        }
        {
            const AnchorType owner = anchor_view::anchorFor(anchor_view::kTEnd);
            const double v = autoFor(owner).t_end[lead];
            if (v >= 0.0) out.t_end = v + frameShift(lead, owner, frame);
        }
        return out;
    }

    // ---- NO STORED GLYPH ANY MORE ----------------------------------------
    //
    // syncReactiveGlyphs lived here. It cached the reactive P-peak into
    // BankMarkerSet::p_peak, and that field is gone: the bars fully determine
    // it, so a stored copy was a second answer that could disagree with them,
    // and persisting it made the disagreement survive a reload.
    //
    // Every reader now calls FeatureMarks::reactive_ecg on the bar set it
    // already has. The markings CSV below does exactly that (rxUser / rxAuto),
    // and BinPlotWidget::reactiveGlyphs does it per paint, so the X on screen
    // and the CSV column come from one function with one set of brackets.
    //
    // CALLERS TO UPDATE: TemplateViewerWindow called this after every bar edit
    // and once after seeding. Those calls should simply be deleted -- there is
    // nothing left to sync -- but the marker get/set switches that read and
    // wrote BankMarkerSet::p_peak for BinPlotWidget::EcgPPeak need to return
    // reactive_ecg's result instead, and setMarker on EcgPPeak should be a
    // no-op because a glyph is not draggable.

    std::vector<double> ppgTemplate;
    std::vector<double> ppg_template_iqr;

    // Foot-anchored arterial background traces (empty when absent).
    std::vector<double> abpTemplate;
    std::vector<double> artTemplate;
    std::vector<double> artPulmTemplate;
    // Per-sample std for each arterial template (empty when absent / not
    // computed). Same length as the matching template when present.
    std::vector<double> abpTemplate_iqr;
    std::vector<double> artTemplate_iqr;
    std::vector<double> artPulmTemplate_iqr;

    //error markings made via user right click
    bool    bad_r_ch[3] = { false, false, false };
    uint8_t bad_ppg = 0;   // 0 = ok, 1 = bad, 2 = no ppg

    // Per-anchor USER marker positions. Each alignment anchor (R, Q_ONSET,
    // J_POINT, T_PEAK, ...) has its OWN independent set of draggable ECG
    // markers, because a marker's sample column is only meaningful relative to
    // the alignment it was placed on. Keyed by AnchorType; a missing anchor
    // means "not marked yet for that anchor" (seed fresh on load). Use marks().
    // ---- ONE HOME FOR EVERY LANDMARK: THE BANK TEMPLATE ------------------
    //
    // There used to be two. A bin-level MarkerSet held slot 0's landmarks as
    // per-lead arrays of 3, while slots 1..N held theirs as scalars in
    // BankTemplate::markers_by_anchor. Same landmark, two shapes, two homes --
    // and the consequences were not cosmetic:
    //
    //   * onMarkerMoved (the slot-0 drag path) propagated through the BIN's set
    //     and so had no slot dimension at all. It could not move a sub-template
    //     column no matter what, which is why dragging _A never moved _E.
    //   * Three functions disagreed about what slot 0 was: markingSlotsForBin
    //     treated it as a bank slot, leadsForBinTemplate fell back to chN_raw
    //     for it, onMarkerMoved wrote it to the bin.
    //   * Every operation spanning both had to know which storage applied, and
    //     there were four propagation loops between them.
    //
    // Now slot 0 is a slot like any other. The bank is ALREADY per lead, so the
    // array of 3 was itself a leftover from before banks existed; BankMarkerSet
    // is the correct shape and is the only one.
    //
    // Reach them through slotMarks(). Nothing outside this struct should index
    // ecg_bank[...].templates[...].markers_by_anchor directly.

    // Writable, and CREATES THE SLOT if the bank is shorter -- which is what the
    // markings reader needs, since it builds bins with no banks and has to put
    // the saved landmarks somewhere before the merge. A created slot carries
    // marks and nothing else: no tmpl, no members, so it earns no column and
    // draws nothing until the real bank arrives.
    tbank::BankMarkerSet& slotMarks(int lead, int slot, AnchorType a) {
        static tbank::BankMarkerSet kScratch;
        if (lead < 0 || lead > 2 || slot < 0) { kScratch = {}; return kScratch; }
        tbank::TemplateBank& bk = ecg_bank[lead];
        if (slot >= static_cast<int>(bk.templates.size()))
            bk.templates.resize(static_cast<size_t>(slot) + 1);
        return bk.templates[slot].marks(static_cast<int32_t>(a));
    }

    // Read-only, and NEVER creates. Returns an all -1 set for a slot that does
    // not exist, because BankTemplate::marks() is operator[] on a map and a read
    // through it inserts -- the trap that used to suppress a template's
    // auto-detection permanently.
    const tbank::BankMarkerSet& slotMarks(int lead, int slot, AnchorType a) const {
        static const tbank::BankMarkerSet kEmpty;
        if (lead < 0 || lead > 2 || slot < 0) return kEmpty;
        const tbank::TemplateBank& bk = ecg_bank[lead];
        if (slot >= static_cast<int>(bk.templates.size())) return kEmpty;
        return bk.templates[slot].marks(static_cast<int32_t>(a));
    }

    // Slots that carry a saved marker set, for the writer and the merge. Not
    // the bank's size: a bin read from the markings file has only the slots the
    // operator's session had.
    //
    // NAME IS A MILD LIE and worth fixing separately: it returns
    // templates.size(), i.e. EVERY slot, not only the marked ones. The writer
    // is correct because it walks all of them; a future reader who trusts the
    // name would wrongly assume unmarked slots are skipped.
    int markedSlotCount(int lead) const {
        if (lead < 0 || lead > 2) return 0;
        return static_cast<int>(ecg_bank[lead].templates.size());
    }

    // R peak column: auto-only, re-derived each pass from the template r_col;
    // NOT per-anchor and NOT persisted. Kept flat, and now DOUBLE like every
    // other position -- the refine stage produces a sub-sample R.
    double r_peak_ch[3] = { -1, -1, -1 };
    //sample indicies for each of the 3 ECG channels, auto-detected. -1 = unmarked / not applicable.
    // Auto-detected ECG landmark positions, stored as DOUBLE (sub-sample
    // precision from the fit/refine stages; -1 = unset). NOT serialized --
    // recomputed every loadSubject.
    double p_peak_auto_ch[3] = { -1, -1, -1 };
    double q_peak_auto_ch[3] = { -1, -1, -1 };
    double q_begin_auto_ch[3] = { -1, -1, -1 };
    double r_peak_auto_ch[3] = { -1, -1, -1 };
    double s_end_auto_ch[3] = { -1, -1, -1 };
    // (t_begin_auto_ch removed with the marker -- see BankMarkerSet. The
    //  T-peak bracket is s_end/t_end on both the auto and the bar side.)
    double t_end_auto_ch[3] = { -1, -1, -1 };
    double p_begin_auto_ch[3] = { -1, -1, -1 };
    bool q_begin_found_auto_ch[3] = { false, false, false };

    // ---- PPG -------------------------------------------------------------
    //
    // THE THREE BARS. Sub-sample positions into ppgTemplate, shared across
    // channels. These are the only PPG values in template_markings.bin.
    double ppg_onset = -1;
    double ppg_dicrotic = -1;
    double ppg_end = -1;

    // REACTIVE, NOT PERSISTED. All five come back from
    // FeatureMarks::reactive_ppg bracketed by the bars above plus the
    // auto-detected systolic peak; ppg_peak itself is auto-only (markerAtX
    // refuses to hand it out). They live here as a cache for consumers that
    // want them without re-deriving, are rebuilt on every load, and appear in
    // no binary format.
    double ppg_peak = -1;
    double ppg_peak2 = -1;
    double ppg_t50 = -1;   // 50% up the upslope, foot -> systolic peak
    double ppg_t80 = -1;   // 80% up the upslope, foot -> systolic peak
    // T80_rise: upslope position at the SAME absolute amplitude as the 80%
    // DOWNSLOPE crossing (not an 80%-of-onset->peak level). pw80 = t80 -
    // t80_rise, the pulse width at that level. No glyph.
    double ppg_t80_rise = -1.0;
    double ppg_pw80 = -1.0;

    // Construction-time PPG fiducials (peak = max in [R1,R2], foot = min in
    // [R1,peak]), computed from the real R-pair interval at template build and
    // carried through the template file. seed_all uses these directly instead
    // of re-detecting the peak on the multi-pulse template. -1 if unavailable.
    double ppg_peak_construct = -1;
    double ppg_onset_construct = -1;

    // Arterial channels (ABP / ART / ART_PULM): same marker set as PPG
    // (onset, peak, dicrotic, peak2, end) and same issue flag semantics
    // (0 = ok, 1 = bad, 2 = channel absent). Positions are sub-sample indices
    // into the matching *Template vector above; -1 = unmarked / not applicable.
    //
    // ALL FIVE ARE BARS on every arterial channel, unlike PPG: markerAtX
    // excludes PpgPeak / PpgPeak2 / PpgT50 / PpgT80 but nothing in the
    // Abp*/Art*/ArtPulm* range, because these channels have no detector to
    // freeze a glyph from. So all fifteen are persisted.
    uint8_t abp_issue = 0;
    double abp_onset = -1, abp_peak = -1, abp_dicrotic = -1, abp_peak2 = -1, abp_end = -1;

    uint8_t art_issue = 0;
    double art_onset = -1, art_peak = -1, art_dicrotic = -1, art_peak2 = -1, art_end = -1;

    uint8_t art_pulm_issue = 0;
    double art_pulm_onset = -1, art_pulm_peak = -1, art_pulm_dicrotic = -1,
        art_pulm_peak2 = -1, art_pulm_end = -1;

    /*these initialized autodetect values(none will actually spend their whole lives at - 1) represnt the initial
    positions of movable markers and the permanant positions of nonmovable markers.
    DOUBLES: the pulse detectors (find_foot_pulseox, the dicrotic finder, the
    VPG/APG derivative zero-crossings) all return interpolated positions, and
    these used to round them on the way in.*/
    double ppg_onset_auto = -1, ppg_peak_auto = -1, ppg_dicrotic_auto = -1, ppg_peak2_auto = -1, ppg_end_auto = -1;
    double ppg_t80_auto = -1, ppg_t50_auto = -1;
    double ppg_u_auto = -1, ppg_v_auto = -1, ppg_w_auto = -1;
    double ppg_a_auto = -1, ppg_b_auto = -1, ppg_c_auto = -1, ppg_d_auto = -1, ppg_e_auto = -1, ppg_f_auto = -1;
    bool ppg_dicrotic_found_auto = false, ppg_peak2_found_auto = false;
    double abp_onset_auto = -1, abp_peak_auto = -1, abp_dicrotic_auto = -1, abp_peak2_auto = -1, abp_end_auto = -1;
    double art_onset_auto = -1, art_peak_auto = -1, art_dicrotic_auto = -1, art_peak2_auto = -1, art_end_auto = -1;
    double ppg_p1_auto = -1, ppg_p2_auto = -1;
    double art_pulm_onset_auto = -1, art_pulm_peak_auto = -1, art_pulm_dicrotic_auto = -1, art_pulm_peak2_auto = -1, art_pulm_end_auto = -1;

    // Derived PPG indices (DeepEntropyX Section 6.3). Like every other *_auto
    // field these are recomputed by seed_all each pass and are NOT part of the
    // serialized template_markings.bin format. SI stays NaN until a subject
    // height is threaded into seed_all.
    double ppg_ba_auto = NAN, ppg_ca_auto = NAN, ppg_da_auto = NAN;
    double ppg_ea_auto = NAN, ppg_fa_auto = NAN;
    double ppg_agi_auto = NAN, ppg_ri_auto = NAN, ppg_si_auto = NAN;
    uint16_t ppg_found_mask_auto = 0;

    // Three-tier dicrotic-notch provenance (E-5), recomputed by seed_all.
    // A TIER AND A CONFIDENCE, not positions -- tier stays int.
    int    ppg_dn_tier_auto = 3;          // 1=IEM, 2=Windkessel, 3=absent
    double ppg_dn_confidence_auto = 0.0;

    // ---- THE PPG COUNTERPART OF THE REMOVED syncReactiveGlyphs -----------
    //
    // The bars are the input, every reactive value is the output, and this is
    // the ONLY thing that assigns them. ppg_peak is taken from ppg_peak_auto
    // rather than trusted from wherever it currently sits, because it is
    // auto-only and reactive_ppg needs it as a bracket.
    //
    // Call after seeding and after a markings-bin restore. The bin stores bars
    // only, so on restore these have to be rederived or they keep whatever the
    // auto-seed left while the bars come from the file.
    //
    // NOTE: reactive_ppg takes `dicrotic` but does not read it -- t50, t80,
    // t80_rise, pw80 and peak2 are all derived from onset, peak and end. It is
    // passed anyway so this call matches BinPlotWidget::reactiveGlyphs.
    void syncReactivePpg() {
        ppg_peak = ppg_peak_auto;
        if (ppgTemplate.empty()) return;
        const FeatureMarks::ReactivePpg rx = FeatureMarks::reactive_ppg(
            ppgTemplate, ppg_onset, ppg_peak, ppg_dicrotic, ppg_end);
        ppg_t50 = rx.t50;
        ppg_t80 = rx.t80;
        ppg_peak2 = rx.peak2;
        ppg_t80_rise = rx.t80_rise;
        ppg_pw80 = rx.pw80;
    }
};




// ---------------------------------------------------------------------------
// Read: convert template_io::TemplateFile -> std::vector<TemplateBin>
// ---------------------------------------------------------------------------
inline std::vector<TemplateBin> readTemplateInfoBin(const std::string& path,
    AnchorType /*anchor*/ = AnchorType::R_PEAK) {
    template_io::TemplateFile tf = template_io::read_template_binfile(path);

    // EVERY ANCHOR, NOT ONE. This used to take an anchor and project that one
    // block into chN, which is what made a second alignment cost a template
    // regeneration and a window reload. The parameter is vestigial -- there is
    // nothing left to select -- and is kept only so existing call sites
    // compile; delete it once they are all updated.

    std::vector<TemplateBin> bins(tf.bins.size());
    for (size_t i = 0; i < tf.bins.size(); ++i) {
        const auto& src = tf.bins[i];
        auto& dst = bins[i];

        dst.index = static_cast<uint64_t>(i);
        dst.bad_segment = src.bad_segment;
        dst.ch1_n_beats_raw = src.ch1_n_beats_raw;
        dst.ch2_n_beats_raw = src.ch2_n_beats_raw;
        dst.ch3_n_beats_raw = src.ch3_n_beats_raw;
        dst.ppg_n_beats = src.ppg_n_beats;
        // template_io carries these as whole columns; widening is exact.
        dst.ppg_peak_construct = src.ppg_peak_col;
        dst.ppg_onset_construct = src.ppg_onset_col;
        dst.ppgTemplate = src.ppgTemplate;
        dst.ppg_template_iqr = src.ppg_template_iqr;

        // Section 4.6 banks. Empty on a pre-v3 templates file, which is the
        // correct reading: one template per channel IS a bank of size one, and
        // slot 0 of an absent bank is the chN_raw template already copied
        // above.
        dst.ecg_bank = src.ecg_bank;
        dst.ppg_bank = src.ppg_bank;
        dst.abpTemplate = src.abpTemplate;
        dst.artTemplate = src.artTemplate;
        dst.artPulmTemplate = src.artPulmTemplate;
        dst.abpTemplate_iqr = src.abpTemplate_iqr;
        dst.artTemplate_iqr = src.artTemplate_iqr;
        dst.artPulmTemplate_iqr = src.artPulmTemplate_iqr;

        // chN is ALWAYS the R base: what the grid draws, and the frame every
        // other alignment's columns are translated into.
        auto project = [](const template_io::ChannelMethodTemplate& c,
            ChannelTemplateData& d) {
                d.ecgTemplate_raw = c.ecgTemplate;
                d.ecg_template_raw_iqr = c.ecg_template_iqr;
                d.alignment_point_raw = c.alignment_point;
                d.r_col_raw = c.r_col;
                d.median_rr_samples = c.median_rr_samples;
            };
        project(src.ch1_raw, dst.ch1);
        project(src.ch2_raw, dst.ch2);
        project(src.ch3_raw, dst.ch3);

        // Every alignment block the file carries, R included if it is in
        // there (chFor short-circuits R to the base regardless, so a
        // duplicate costs nothing and an absent one costs nothing either).
        // Guard the per-bin index: an anchor block shorter than the base
        // shouldn't happen, but never index past it.
        for (const auto& kv : tf.raw_anchors) {
            if (i >= kv.second.size()) continue;
            std::array<ChannelTemplateData, 3> trio;
            for (int c = 0; c < 3; ++c) project(kv.second[i][c], trio[c]);
            dst.anchored[kv.first] = std::move(trio);
        }

        // Per-anchor BANK SLOT averages (v6 section). Absent on any file
        // written before that section existed, in which case bankSlotFor
        // returns nullptr and the sub-template panels keep their old,
        // now-honestly-labelled behaviour.
        // NO LOCAL DECLARATIONS IN THIS LOOP, deliberately. Every earlier form
        // declared something whose type mentioned AnchoredBankSlot -- a local
        // vector, a reference to a map element -- and MSVC failed to parse the
        // declarator every time. Once that happens it reads the following
        // `slots[sl]` as an attribute list or a lambda introducer, so the real
        // error is buried under "'sl': attribute not found", "expected a '{'
        // introducing a lambda body", and a syntax error on every `.member`
        // after it. Written long-hand: subscripting and assignment only.
        for (const auto& kv : tf.bank_anchors) {
            if (i >= kv.second.size()) continue;
            for (int c = 0; c < 3; ++c) {
                dst.anchored_bank[kv.first * 4 + c].resize(kv.second[i][c].size());
                for (size_t sl = 0; sl < kv.second[i][c].size(); ++sl) {
                    dst.anchored_bank[kv.first * 4 + c][sl].tmpl =
                        kv.second[i][c][sl].tmpl;
                    dst.anchored_bank[kv.first * 4 + c][sl].tmpl_iqr =
                        kv.second[i][c][sl].tmpl_iqr;
                    dst.anchored_bank[kv.first * 4 + c][sl].n_members =
                        kv.second[i][c][sl].n_members;
                }
            }
        }
    }
    return bins;
}

// ---------------------------------------------------------------------------
// template_markings.bin layout -- BARS ONLY, POSITIONS AS FLOAT64:
//
//   header:
//     uint64  numBins
//
//   per bin:
//     uint64  index
//     uint8   bad_r_ch1, bad_r_ch2, bad_r_ch3
//     uint8   ppg_issue          (0 = ok, 1 = bad, 2 = no ppg)
//
//     -- ECG bars, per lead 0..2: --
//     int32   slotCount
//       per slot:
//         int32   anchorCount
//           per anchor:
//             int32   anchorTag
//             float64 p_begin, q_begin, s_end, t_end
//
//     -- PPG bars: --
//     float64 ppg_onset, ppg_dicrotic, ppg_end
//
//     -- arterial, one block each for ABP, ART, ART_PULM: --
//     uint8   <chan>_issue
//     float64 <chan>_onset, _peak, _dicrotic, _peak2, _end
//
// Positions use -1 as the "unmarked / not applicable" sentinel; test `< 0`.
// Structural fields (counts, tags, index, issue flags) stay integral.
//
// NO VERSION FIELD EXISTS -- the header is a bare bin count -- so a markings
// file written under any earlier layout misparses from its first marker set
// onward. That covers the removal of t_begin, the removal of p_peak, the
// removal of the six auto-only PPG values, and this widening to float64.
// Re-mark rather than migrate.
// ---------------------------------------------------------------------------
inline void writeTemplateMarkingsBin(const std::string& path,
    const std::vector<TemplateBin>& bins) {
    std::ofstream f(path, std::ios::binary);
    if (!f.is_open())
        throw std::runtime_error("cannot open for write: " + path);

    //header: just bin count
    uint64_t n = bins.size();
    f.write(reinterpret_cast<const char*>(&n), 8);

    auto w8 = [&](uint8_t v) { f.write(reinterpret_cast<const char*>(&v), 1); };
    auto w32 = [&](int v) { int32_t i = v; f.write(reinterpret_cast<const char*>(&i), 4); };
    auto w64d = [&](double v) { f.write(reinterpret_cast<const char*>(&v), 8); };

    for (const auto& b : bins) {
        uint64_t idx = b.index;
        f.write(reinterpret_cast<const char*>(&idx), 8);

        w8(b.bad_r_ch[0] ? 1 : 0);
        w8(b.bad_r_ch[1] ? 1 : 0);
        w8(b.bad_r_ch[2] ? 1 : 0);
        w8(b.bad_ppg);

        // ECG user landmarks, PER LEAD, PER SLOT, PER ANCHOR:
        //
        //   per lead 0..2:
        //     int32    slotCount
        //     per slot:
        //       int32    anchorCount
        //       per anchor:
        //         int32    anchorTag
        //         float64  p_begin, q_begin, s_end, t_end
        //
        // FOUR FIELDS: THE BARS ONLY. p_peak left with t_begin -- it is a
        // reactive glyph, fully determined by the bars that bracket it, so
        // storing it only created a second answer that could disagree with
        // them. Readers call FeatureMarks::reactive_ecg instead.
        //
        // THE SLOT DIMENSION IS THE POINT. The old layout was per anchor only,
        // with each field an array of 3 leads -- it could store slot 0 and
        // nothing else, so a sub-template column's landmarks had nowhere to go
        // in this file and were kept in _templates.bin instead, where the
        // pipeline overwrote them on the next anchor pass.
        //
        // LEAD IS THE OUTER LOOP because a bank is per lead and the leads
        // legitimately disagree on slot count: a morphology separable on CH1 may
        // not be on CH2. A single shared slotCount would have to be the max and
        // pad the rest.
        //
        // Counts and tags stay int32: they are structure, not position.
        //
        // R peak is auto-only, re-derived from each template's own r_col every
        // pass, and is not written here.
        for (int lead = 0; lead < 3; ++lead) {
            const int nSlots = b.markedSlotCount(lead);
            w32(nSlots);
            for (int slot = 0; slot < nSlots; ++slot) {
                const auto& byAnchor =
                    b.ecg_bank[lead].templates[slot].markers_by_anchor;
                w32(static_cast<int>(byAnchor.size()));
                for (const auto& kv : byAnchor) {
                    w32(kv.first);                      // AnchorType tag
                    const tbank::BankMarkerSet& m = kv.second;
                    w64d(m.p_begin); w64d(m.q_begin);
                    w64d(m.s_end);   w64d(m.t_end);
                }
            }
        }

        // PPG: THE THREE BARS ONLY. t50 / peak / peak2 / t80 / t80_rise /
        // pw80 are auto-only glyphs -- markerAtX refuses to hand any of them
        // out -- and all five reactive ones come back from reactive_ppg,
        // bracketed by these bars. TemplateBin::syncReactivePpg rederives them
        // after every load.
        w64d(b.ppg_onset);
        w64d(b.ppg_dicrotic);
        w64d(b.ppg_end);

        // Arterial block: ABP, ART, ART_PULM (issue + 5 positions each). All
        // five ARE bars on these channels -- see the struct note.
        w8(b.abp_issue);
        w64d(b.abp_onset); w64d(b.abp_peak); w64d(b.abp_dicrotic); w64d(b.abp_peak2); w64d(b.abp_end);
        w8(b.art_issue);
        w64d(b.art_onset); w64d(b.art_peak); w64d(b.art_dicrotic); w64d(b.art_peak2); w64d(b.art_end);
        w8(b.art_pulm_issue);
        w64d(b.art_pulm_onset); w64d(b.art_pulm_peak); w64d(b.art_pulm_dicrotic);
        w64d(b.art_pulm_peak2); w64d(b.art_pulm_end);
    }
}


struct EcgFeatures {
    // Sub-sample peak positions throughout. q and s come from the refined
    // finders in FeatureMarks; r is now double as well, since r_peak_ch and
    // AnchorAuto::r_peak both carry the refined column.
    double q_idx = -1.0, s_idx = -1.0;
    double r_idx = -1.0;
    double qrs_ms = NAN, qt_ms = NAN;
};

// Every position in and out is a sub-sample double. The two FeatureMarks
// finders below still take an int R column, so the rounding happens HERE, at
// the one call that needs it, instead of at every caller.
inline EcgFeatures computeEcgFeatures(const std::vector<double>& ecg,
    double p_peak, double q_begin, double r_peak, double s_end, double t_end,
    double rateHz)
{
    EcgFeatures f;
    const double N = static_cast<double>(ecg.size());
    const double msPerSamp = (rateHz > 0.0) ? 1000.0 / rateHz : NAN;

    if (q_begin >= 0.0 && s_end >= q_begin) f.qrs_ms = (s_end - q_begin) * msPerSamp;
    if (q_begin >= 0.0 && t_end >= q_begin) f.qt_ms = (t_end - q_begin) * msPerSamp;

    if (r_peak >= 0.0 && r_peak <= N - 1.0) f.r_idx = r_peak;
    const int rInt = (r_peak >= 0.0) ? static_cast<int>(std::lround(r_peak)) : -1;
    // Q for the q_peak column: same canonical finder compute_q_onset uses.
    // Mirrors compute_s_peak's signature below. -1 when there is no Q trough.
    f.q_idx = FeatureMarks::compute_q_peak(ecg, rInt, rateHz);   // sub-sample
    // S for |R|+|S| = first opposite-polarity trough after R (robust; not the
    // max over [R, s_end], which depends on where s_end sits).
    f.s_idx = FeatureMarks::compute_s_peak(ecg, rInt, rateHz);   // sub-sample
    return f;
}

// (EcgColSpec / ecgCols deleted: declared, never referenced, and wrong --
//  it claimed all six ECG markers were user-placed, which was never true and
//  is emphatically not true now that four bars and a set of glyphs are
//  distinguished. The live column order is ecgPointNames, inside
//  writeTemplateMarkingsCsv.)
inline constexpr const char* ppgCols[] = { "ppg_onset","ppg_p50","ppg_peak","ppg_dicr","ppg_peak2","ppg_t80","ppg_t80_rise","ppg_end" };
inline constexpr const char* abpCols[] = { "abp_onset","abp_peak","abp_dicr","abp_peak2","abp_end" };
inline constexpr const char* artCols[] = { "art_onset","art_peak","art_dicr","art_peak2","art_end" };
inline constexpr const char* artPulmCols[] = { "art_pulm_onset","art_pulm_peak","art_pulm_dicr","art_pulm_peak2","art_pulm_end" };


//one unifornm table for pp autodetected pulses
// idx is a POINTER TO DOUBLE now, matching the widened *_auto fields.
struct PulseAutoGlyph {
    const char* name;
    double TemplateBin::* idx;
    bool TemplateBin::* found;
    const char* foundName;
};
inline constexpr PulseAutoGlyph ppg_and_artpulse_automated_markers[] = {
    { "ppg_foot",          &TemplateBin::ppg_onset_auto,    nullptr,                                 nullptr },
    { "ppg_systolic_peak", &TemplateBin::ppg_peak_auto,     nullptr,                                 nullptr },
    { "ppg_dicr",          &TemplateBin::ppg_dicrotic_auto, &TemplateBin::ppg_dicrotic_found_auto,   "ppg_notch_found" },
    { "ppg_end",           &TemplateBin::ppg_end_auto,      nullptr,                                 nullptr },
    { "vpg_u",             &TemplateBin::ppg_u_auto,        nullptr,                                 nullptr },
    { "vpg_v",             &TemplateBin::ppg_v_auto,        nullptr,                                 nullptr },
    { "vpg_w",             &TemplateBin::ppg_w_auto,        nullptr,                                 nullptr },
    { "apg_a",             &TemplateBin::ppg_a_auto,        nullptr,                                 nullptr },
    { "apg_b",             &TemplateBin::ppg_b_auto,        nullptr,                                 nullptr },
    { "apg_c",             &TemplateBin::ppg_c_auto,        nullptr,                                 nullptr },
    { "apg_d",             &TemplateBin::ppg_d_auto,        nullptr,                                 nullptr },
    { "apg_e",             &TemplateBin::ppg_e_auto,        nullptr,                                 nullptr },
    { "apg_f",             &TemplateBin::ppg_f_auto,        nullptr,                                 nullptr },
    { "jpg_p1",            &TemplateBin::ppg_p1_auto,       nullptr,                                 nullptr },
    { "jpg_p2",            &TemplateBin::ppg_p2_auto,       nullptr,                                 nullptr },
};

// ---------------------------------------------------------------------------
// writeTemplateMarkingsCsv
//
// TAKES A STREAM, NOT A PATH. The caller builds each alignment's part in an
// ostringstream, suffixes its value columns and merges the set in one write,
// so nothing is staged through a temp file. A path-taking wrapper follows at
// the end of this header for callers outside the viewer.
//
// For every marker we emit six columns, in this order:
//   {name}_y_normalized_autodetect
//   {name}_y_normalized_user
//   {name}_y_mv_raw_autodetect
//   {name}_y_mv_raw_user
//   {name}_x_ms_autodetect
//   {name}_x_ms_user
//
// Normalization:
//   ECG:   y / Global_Ref_ecg(ch), where Global_Ref_ecg = median across bins
//          of (|R_peak_y| + |S_peak_y|) using that variant's own R/S positions.
//   Pulse: (100 * (y - foot_y) / foot_y) / Global_Ref_pulse(chan), where
//          Global_Ref_pulse = median across bins of 100*(peak - foot)/foot.
//          Autodetect uses foot_auto; user uses foot.
// Both refs computed once per subject inside this function; they use the
// autodetect positions so the "reference" is stable regardless of user edits.
//
// EVERY POSITION IS A DOUBLE END TO END, so every amplitude lookup goes
// through FeatureMarks::sample_at and every _x_ms column carries its
// fraction. The bar side is no longer widened-from-int at the call site --
// BankMarkerSet is double -- so the "_user" millisecond columns are now as
// precise as the "_auto" ones.
//
// Intervals (qrs_ms, qt_ms) emit as pairs: {name}_ms_autodetect, {name}_ms_user.
// Q peak and S peak are computed inside the QRS -- their autodetect variant
// uses (Q_begin_auto, R_peak_auto, S_end_auto); user variant uses the current
// user markers.
// ---------------------------------------------------------------------------
inline void writeTemplateMarkingsCsv(std::ostream& f,
    const std::vector<TemplateBin>& bins,
    const std::string& fileID,
    double sampleRateHz,
    AnchorType anchor,
    MarkingsCsvSection section)
{
    const bool wantEcg = (section == MarkingsCsvSection::EcgOnly);
    const bool wantPulse = (section == MarkingsCsvSection::PulseOnly);


    // ---- helpers -----------------------------------------------------------
    auto medianFinite = [](std::vector<double> v) -> double {
        v.erase(std::remove_if(v.begin(), v.end(),
            [](double x) { return !std::isfinite(x); }), v.end());
        if (v.empty()) return std::nan("");
        const size_t mid = v.size() / 2;
        std::nth_element(v.begin(), v.begin() + mid, v.end());
        return v[mid];
        };

    // ---- global refs (subject-wide, computed from AUTO positions) ----------
    double ecgRef[3] = { std::nan(""), std::nan(""), std::nan("") };
    if (wantEcg) {
        for (int c = 0; c < 3; ++c) {
            std::vector<double> vals;
            for (const auto& b : bins) {
                if (b.bad_segment || b.bad_r_ch[c]) continue;

                const std::vector<double>& ecg = b.chFor(c, AnchorType::R_PEAK).ecgTemplate_raw;
                if (ecg.empty()) continue;
                const auto aaR = b.autoFor(AnchorType::R_PEAK);
                // NO ROUNDING: computeEcgFeatures takes doubles now.
                EcgFeatures ft = computeEcgFeatures(ecg,
                    aaR.p_peak[c], aaR.q_begin[c], aaR.r_peak[c],
                    aaR.s_end[c], aaR.t_end[c],
                    sampleRateHz);
                if (ft.r_idx < 0.0 || ft.s_idx < 0.0) continue;
                const double last = static_cast<double>(ecg.size()) - 1.0;
                if (ft.r_idx > last || ft.s_idx > last) continue;
                // |R| + |S| at the SUB-SAMPLE positions: interpolated, not
                // read from a rounded column. On a steep S limb the two differ
                // by more than the amplitude tolerance this feeds -- and R is
                // now interpolated too, for the same reason.
                const double ry = FeatureMarks::sample_at(ecg, ft.r_idx);
                const double sy = FeatureMarks::sample_at(ecg, ft.s_idx);
                if (std::isnan(ry) || std::isnan(sy)) continue;
                vals.push_back(std::abs(ry) + std::abs(sy));
            }
            ecgRef[c] = medianFinite(std::move(vals));
        }
    }
    double refPpg = std::nan(""), refAbp = std::nan(""), refArt = std::nan(""), refArtPulm = std::nan("");
    if (wantPulse) {
        // Pointer-to-member is double now; the foot/peak amplitudes are
        // interpolated rather than subscripted.
        auto pulseRefAuto = [&](const std::vector<double> TemplateBin::* trace,
            double TemplateBin::* footAuto,
            double TemplateBin::* peakAuto,
            bool checkPpgIssue) -> double
            {
                std::vector<double> vals;
                for (const auto& b : bins) {
                    if (b.bad_segment) continue;
                    if (checkPpgIssue && b.bad_ppg != 0) continue;
                    const auto& v = b.*trace;
                    const double fi = b.*footAuto;
                    const double pi = b.*peakAuto;
                    const double last = static_cast<double>(v.size()) - 1.0;
                    if (fi < 0.0 || pi < 0.0 || fi > last || pi > last) continue;
                    const double fy = FeatureMarks::sample_at(v, fi);
                    const double py = FeatureMarks::sample_at(v, pi);
                    if (std::isnan(fy) || std::isnan(py) || std::abs(fy) < 1e-12) continue;
                    vals.push_back(100.0 * (py - fy) / fy);
                }
                return medianFinite(std::move(vals));
            };
        refPpg = pulseRefAuto(&TemplateBin::ppgTemplate, &TemplateBin::ppg_onset_auto, &TemplateBin::ppg_peak_auto, true);
        refAbp = pulseRefAuto(&TemplateBin::abpTemplate, &TemplateBin::abp_onset_auto, &TemplateBin::abp_peak_auto, false);
        refArt = pulseRefAuto(&TemplateBin::artTemplate, &TemplateBin::art_onset_auto, &TemplateBin::art_peak_auto, false);
        refArtPulm = pulseRefAuto(&TemplateBin::artPulmTemplate, &TemplateBin::art_pulm_onset_auto, &TemplateBin::art_pulm_peak_auto, false);
    }

    // ---- header ------------------------------------------------------------
    // Keys. ppg_issue is PULSE metadata -> only in the pulse section.
    f << "file_id,bin_index,bad_r_ch1,bad_r_ch2,bad_r_ch3";
    if (wantPulse) f << ",ppg_issue";

    // ECG point columns + 2 interval columns per channel.
    //
    // p_begin ADDED. It was absent: the list started at p_peak, so the P-onset
    // bar -- one of the four the operator places, and the entire reason the P
    // alignment exists -- had no column in this file at all. (The aligned-
    // template CSV did emit it, so the gap was here only.)
    //
    // t_begin REMOVED with the marker. Nothing ever set it, so its columns
    // were structurally blank, and the T-peak bracket that referenced it
    // reported nothing while the on-screen X sat correctly between the J-point
    // and T-end bars.
    static const char* ecgPointNames[] = {
        "p_begin", "p_peak", "q_begin", "q_peak", "r_peak", "s_peak",
        "s_end",   "t_peak", "t_end"
    };
    static const char* ecgIntervalNames[] = { "qrs", "qt" };

    auto emitEcgPointHeader = [&](const char* name, int c) {
        // r_peak is the only column with no user variant anywhere: it is the
        // alignment anchor, re-derived from r_col every load, never placed, and
        // it holds no BankMarkerSet field. Every other glyph DOES get a user
        // column -- p_peak and t_peak are reactive, so they have a distinct
        // value under the detector's brackets and under the operator's bars,
        // and the bar-bracketed one is the X on screen.
        //
        // ONE SOURCE for this rule: the row loop asks anchor_view too, so the
        // header and the body cannot disagree about column count.
        const bool userToo = anchor_view::hasUserColumn(name);
        f << ',' << name << "_ch" << c << "_y_norm_auto";
        if (userToo) f << ',' << name << "_ch" << c << "_y_norm_user";
        f << ',' << name << "_ch" << c << "_y_mv_auto";
        if (userToo) f << ',' << name << "_ch" << c << "_y_mv_user";
        f << ',' << name << "_ch" << c << "_x_ms_auto";
        if (userToo) f << ',' << name << "_ch" << c << "_x_ms_user";
        };
    auto emitIntervalHeader = [&](const char* name, int c) {
        f << ',' << name << "_ch" << c << "_ms_auto"
            << ',' << name << "_ch" << c << "_ms_user";
        };
    // Pulse: 6 cols per marker.
    auto emitPulsePointHeader = [&](const char* name) {
        f << ',' << name << "_y_norm_auto"
            << ',' << name << "_y_norm_user"
            << ',' << name << "_y_mv_auto"
            << ',' << name << "_y_mv_user"
            << ',' << name << "_x_ms_auto"
            << ',' << name << "_x_ms_user";
        };
    // Autodetected computed features (no user bar; derived from AUTO markers).
    auto emitAutoFeatHeader = [&](const char* name) {
        f << ',' << name << "_x_ms" << ',' << name << "_y_mv";
        };

    if (wantEcg) {
        for (int c = 1; c <= 3; ++c) {
            for (const char* n : ecgPointNames)     emitEcgPointHeader(n, c);
            // SAME LOOP as the points, so the two interval columns land
            // immediately after this channel's nine point groups and the row
            // loop below can mirror the order without a second pass.
            for (const char* n : ecgIntervalNames)  emitIntervalHeader(n, c);
        }
        for (int c = 1; c <= 3; ++c) {
            char b[64];
            for (const char* g : { "p_wave_auto", "q_onset_auto",
                                   "r_wave_auto", "t_peak_auto" }) {
                std::snprintf(b, sizeof b, "%s_ch%d", g, c);
                emitAutoFeatHeader(b);
            }
        }
    }

    if (wantPulse) {
        for (const char* n : ppgCols)     emitPulsePointHeader(n);
        // PW80 width (t80 - t80_rise), in ms. Single value per bar-set, not a
        // pulse point -- a width has no y/position, so it gets its own two
        // columns rather than the 6-subcolumn pulse layout.
        f << ",ppg_pw80_ms_auto";
        f << ",abp_issue";
        for (const char* n : abpCols)     emitPulsePointHeader(n);
        f << ",art_issue";
        for (const char* n : artCols)     emitPulsePointHeader(n);
        f << ",art_pulm_issue";
        for (const char* n : artPulmCols) emitPulsePointHeader(n);
        for (const auto& gl : ppg_and_artpulse_automated_markers) {
            char nb[64];
            std::snprintf(nb, sizeof nb, "%s_auto", gl.name);
            emitAutoFeatHeader(nb);
            if (gl.foundName) f << ',' << gl.foundName;
        }
        for (const char* g : { "vpg_u_auto", "vpg_v_auto", "vpg_w_auto" })
            emitAutoFeatHeader(g);
    }
    f << '\n';

    // ---- row loop ----------------------------------------------------------
    const double toMs = (sampleRateHz > 0.0) ? 1000.0 / sampleRateHz : 1.0;
    // Default 6 significant figures would drop the sub-sample fraction on
    // 4-digit millisecond values.
    f << std::setprecision(10);

    // Emit one 6-column ECG point group: normalized (auto/user), raw
    // (auto/user), x_ms (auto/user). Any missing piece leaves that field blank.
    // BOTH SIDES ARE SUB-SAMPLE now, so the amplitude is interpolated rather
    // than read from a rounded column -- the local lround this used to do was
    // the last quantisation left in the ECG path.
    auto emitEcgPoint = [&](const std::vector<double>& ecg,
        double idx_auto, double idx_user, double ref, bool userToo)
        {
            auto y_of = [&](double idx) -> double {
                if (idx < 0.0) return std::nan("");
                if (idx > static_cast<double>(ecg.size()) - 1.0) return std::nan("");
                const double y = FeatureMarks::sample_at(ecg, idx);
                return std::isnan(y) ? std::nan("") : y;
                };
            const double y_a = y_of(idx_auto);
            const double y_u = y_of(idx_user);
            const bool refOk = std::isfinite(ref) && ref != 0.0;

            f << ',';   if (std::isfinite(y_a) && refOk) f << (y_a / ref);
            if (userToo) { f << ','; if (std::isfinite(y_u) && refOk) f << (y_u / ref); }
            f << ',';   if (std::isfinite(y_a)) f << y_a;
            if (userToo) { f << ','; if (std::isfinite(y_u)) f << y_u; }
            f << ',';   if (idx_auto >= 0.0) f << (idx_auto * toMs);
            if (userToo) { f << ','; if (idx_user >= 0.0) f << (idx_user * toMs); }
        };

    // Emit one 6-column pulse point group. footIdx_auto / footIdx_user are
    // the "onset" indices for their respective variants (used to compute the
    // local ratio (y - foot)/foot). ref is the channel's global PI median.
    // The amplitude at a fractional position is INTERPOLATED (FeatureMarks::
    // sample_at) rather than read from a rounded column, and the millisecond
    // column carries the fraction, so T80 is no longer quantised to 3.9 ms at
    // 256 Hz -- which matters because Section 6.3's T80 entropy result turns on
    // small differences in exactly that interval.
    auto emitPulsePoint = [&](const std::vector<double>& v,
        double idx_auto, double idx_user,
        double foot_auto, double foot_user, double ref)
        {
            auto y_of = [&](double idx) -> double {
                if (idx < 0.0 || idx >(double)v.size() - 1.0) return std::nan("");
                return FeatureMarks::sample_at(v, idx);
                };
            auto normOf = [&](double y, double fIdx) -> double {
                if (!std::isfinite(y)) return std::nan("");
                if (fIdx < 0.0 || fIdx >(double)v.size() - 1.0) return std::nan("");
                const double fy = FeatureMarks::sample_at(v, fIdx);
                if (std::isnan(fy) || std::abs(fy) < 1e-12) return std::nan("");
                if (!std::isfinite(ref) || ref == 0.0)      return std::nan("");
                const double local = 100.0 * (y - fy) / fy;
                return local / ref;
                };
            const double y_a = y_of(idx_auto);
            const double y_u = y_of(idx_user);
            const double n_a = normOf(y_a, foot_auto);
            const double n_u = normOf(y_u, foot_user);

            f << ',';   if (std::isfinite(n_a)) f << n_a;
            f << ',';   if (std::isfinite(n_u)) f << n_u;
            f << ',';   if (std::isfinite(y_a)) f << y_a;
            f << ',';   if (std::isfinite(y_u)) f << y_u;
            f << ',';   if (idx_auto >= 0.0)    f << (idx_auto * toMs);
            f << ',';   if (idx_user >= 0.0)    f << (idx_user * toMs);
        };

    auto emitIntervalPair = [&](double auto_ms, double user_ms) {
        f << ',';   if (std::isfinite(auto_ms)) f << auto_ms;
        f << ',';   if (std::isfinite(user_ms)) f << user_ms;
        };

    // Autodetected computed feature point (used by both ECG and PPG glyph
    // blocks -- defined once here so it's in scope for both).
    // Sub-sample position in, interpolated amplitude out.
    auto emitAutoFeatPt = [&](const std::vector<double>& sig, double idx) {
        f << ',';   if (idx >= 0.0) f << (idx * toMs);
        const double y = (idx >= 0.0) ? FeatureMarks::sample_at(sig, idx)
            : std::nan("");
        f << ',';   if (std::isfinite(y)) f << y;
        };

    for (const auto& b : bins) {
        f << fileID << ',' << b.index << ','
            << (b.bad_r_ch[0] ? 1 : 0) << ','
            << (b.bad_r_ch[1] ? 1 : 0) << ','
            << (b.bad_r_ch[2] ? 1 : 0);
        if (wantPulse) f << ',' << static_cast<int>(b.bad_ppg);

        // (the chs[3] local is gone: every ECG lookup below goes through
        //  b.chFor(lead, anchor), which selects this block's alignment.)

        if (wantEcg)
            for (int c = 0; c < 3; ++c) {
                // EVERYTHING ECG IN THIS BLOCK IS MEASURED ON THIS ALIGNMENT.
                // Its own average supplies the amplitudes, its own detections
                // supply the glyph columns, and only the one bar it owns gets a
                // user column. Four blocks, four frames.
                const std::vector<double>& ecg = b.chFor(c, anchor).ecgTemplate_raw;
                const double ref = ecgRef[c];
                const auto aa = b.autoFor(anchor);

                // ---- BARS: ONE ALIGNMENT EACH ------------------------------
                // slotMarks(c, 0, anchor) already holds only the bar this
                // alignment owns (maskFor seeds nothing else into it), but the
                // filter is explicit because the column's MEANING depends on
                // it: a P onset measured against the P-aligned average is a
                // different number from one measured against R, so reporting
                // it under another suffix would attribute it to a waveform it
                // was never compared to.
                tbank::BankMarkerSet umk;
                {
                    const tbank::BankMarkerSet& own = b.slotMarks(c, 0, anchor);
                    if (anchor_view::owns(anchor, anchor_view::kPBegin)) umk.p_begin = own.p_begin;
                    if (anchor_view::owns(anchor, anchor_view::kQBegin)) umk.q_begin = own.q_begin;
                    if (anchor_view::owns(anchor, anchor_view::kSEnd))   umk.s_end = own.s_end;
                    if (anchor_view::owns(anchor, anchor_view::kTEnd))   umk.t_end = own.t_end;
                }

                // ---- GLYPHS: EVERY ALIGNMENT -------------------------------
                // A glyph is a measurement, not a judgement, so all four
                // alignments report all of them. The reactive ones need the
                // WHOLE bar set, which spans three alignments, so they are
                // assembled and translated into THIS alignment's frame -- the
                // same bars the operator placed, re-measured on this waveform.
                //
                // reactive_ecg is the function BinPlotWidget::reactiveGlyphs
                // calls, with the same bracket bars (P peak between P-onset and
                // Q-onset; T peak between S-END and T-end). It is also now the
                // ONLY source of p_peak anywhere -- the stored copy in
                // BankMarkerSet is gone, so screen and file cannot diverge.
                const tbank::BankMarkerSet whole = b.userMarks(c, 0, anchor);
                const FeatureMarks::ReactiveEcg rxUser = FeatureMarks::reactive_ecg(
                    ecg, whole.p_begin, whole.q_begin, whole.s_end, whole.t_end, sampleRateHz);
                const FeatureMarks::ReactiveEcg rxAuto = FeatureMarks::reactive_ecg(
                    ecg, aa.p_begin[c], aa.q_begin[c],
                    aa.s_end[c], aa.t_end[c], sampleRateHz);

                EcgFeatures ftAuto = computeEcgFeatures(ecg,
                    aa.p_peak[c], aa.q_begin[c], aa.r_peak[c],
                    aa.s_end[c], aa.t_end[c],
                    sampleRateHz);
                // Derived from the assembled bars, not from `umk`: q_peak,
                // s_peak and the QRS/QT intervals need a whole beat's
                // brackets, and no single alignment's marker set holds one any
                // more. The two intervals come out identical in all four
                // blocks (a duration is frame-free); the two positions differ
                // between blocks by the frame shift.
                EcgFeatures ftUser = computeEcgFeatures(ecg,
                    rxUser.p_peak, whole.q_begin, b.r_peak_ch[c],
                    whole.s_end, whole.t_end, sampleRateHz);

                // Order MUST match ecgPointNames:
                //   p_begin(bar), p_peak(glyph), q_begin(bar), q_peak(computed),
                //   r_peak(glyph), s_peak(computed), s_end(bar),
                //   t_peak(glyph), t_end(bar)
                // BOTH SIDES ARE SUB-SAMPLE now -- the (double) casts on the
                // bar side are gone with BankMarkerSet's int fields.
                struct P { const char* name; double a; double u; };
                const P pts[] = {
                    { "p_begin", aa.p_begin[c],  umk.p_begin  },
                    { "p_peak",  rxAuto.p_peak,  rxUser.p_peak },   // reactive glyph, both sides
                    { "q_begin", aa.q_begin[c],  umk.q_begin  },
                    { "q_peak",  ftAuto.q_idx,   ftUser.q_idx },
                    { "r_peak",  aa.r_peak[c],   -1.0         },   // no user column at all
                    { "s_peak",  ftAuto.s_idx,   ftUser.s_idx },
                    { "s_end",   aa.s_end[c],    umk.s_end    },
                    { "t_peak",  rxAuto.t_peak,  rxUser.t_peak },   // reactive glyph, both sides
                    { "t_end",   aa.t_end[c],    umk.t_end    }
                };
                for (const P& pt : pts) {
                    // SAME FUNCTION the header emitter called, so the two
                    // cannot disagree about how many columns this point has.
                    // This was a bare `k != 3` index test, which meant adding a
                    // point column silently shifted which one lost its user
                    // variant.
                    emitEcgPoint(ecg, pt.a, pt.u, ref,
                        anchor_view::hasUserColumn(pt.name));
                }

                // ---- INTERVALS: ecgIntervalNames ORDER, qrs THEN qt -------
                // A duration is frame-free -- (s_end - q_begin) is the same
                // number whichever alignment's columns it was measured in --
                // so these two come out identical in all four blocks, unlike
                // the positions above. They are emitted per block anyway
                // because the column set has to be the same shape in every
                // part for the sidecar merge to line up.
                //
                // Both come from computeEcgFeatures, which now takes doubles,
                // so the user side is a sub-sample duration rather than a
                // whole number of samples times the sample period.
                emitIntervalPair(ftAuto.qrs_ms, ftUser.qrs_ms);
                emitIntervalPair(ftAuto.qt_ms, ftUser.qt_ms);
            }

        if (wantPulse) {
            // PPG: onset, p50, peak, dicrotic, peak2, t80, t80_rise, end
            // (matches ppgCols). p50/t80/t80_rise/peak2 are reactive --
            // bracketed by onset/peak/end -- so BOTH SIDES of all four come
            // from the shared FeatureMarks::reactive_ppg, the same call the
            // on-screen glyph makes: rxAuto under the detector's brackets,
            // rxUser under the operator's. peak2's auto side used to read
            // ppg_peak2_auto instead, which made it the only reactive column
            // whose two halves came from different functions. The cached
            // ppg_t50 / ppg_t80 / ppg_peak2 fields are deliberately not used
            // for these columns: they are a convenience copy, and reading them
            // here would let a stale cache disagree with the screen.
            const FeatureMarks::ReactivePpg rxAuto = FeatureMarks::reactive_ppg(
                b.ppgTemplate, b.ppg_onset_auto, b.ppg_peak_auto, b.ppg_dicrotic_auto, b.ppg_end_auto);
            const FeatureMarks::ReactivePpg rxUser = FeatureMarks::reactive_ppg(
                b.ppgTemplate, b.ppg_onset, b.ppg_peak, b.ppg_dicrotic, b.ppg_end);
            emitPulsePoint(b.ppgTemplate, b.ppg_onset_auto, b.ppg_onset,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            emitPulsePoint(b.ppgTemplate, rxAuto.t50, rxUser.t50,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            emitPulsePoint(b.ppgTemplate, b.ppg_peak_auto, b.ppg_peak,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            emitPulsePoint(b.ppgTemplate, b.ppg_dicrotic_auto, b.ppg_dicrotic,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            emitPulsePoint(b.ppgTemplate, rxAuto.peak2, rxUser.peak2,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            emitPulsePoint(b.ppgTemplate, rxAuto.t80, rxUser.t80,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            // T80_rise: upslope point at t80's level (a position, like t80).
            emitPulsePoint(b.ppgTemplate, rxAuto.t80_rise, rxUser.t80_rise,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            emitPulsePoint(b.ppgTemplate, b.ppg_end_auto, b.ppg_end,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            // PW80 width, ms only, autodetect bracketing (t80 - t80_rise).
            // Single value; blank when unavailable. Matches the one
            // ppg_pw80_ms_auto header column.
            f << ',';  if (rxAuto.pw80 >= 0.0) f << (rxAuto.pw80 * toMs);

            f << ',' << static_cast<int>(b.abp_issue);
            emitPulsePoint(b.abpTemplate, b.abp_onset_auto, b.abp_onset,
                b.abp_onset_auto, b.abp_onset, refAbp);
            emitPulsePoint(b.abpTemplate, b.abp_peak_auto, b.abp_peak,
                b.abp_onset_auto, b.abp_onset, refAbp);
            emitPulsePoint(b.abpTemplate, b.abp_dicrotic_auto, b.abp_dicrotic,
                b.abp_onset_auto, b.abp_onset, refAbp);
            emitPulsePoint(b.abpTemplate, b.abp_peak2_auto, b.abp_peak2,
                b.abp_onset_auto, b.abp_onset, refAbp);
            emitPulsePoint(b.abpTemplate, b.abp_end_auto, b.abp_end,
                b.abp_onset_auto, b.abp_onset, refAbp);

            f << ',' << static_cast<int>(b.art_issue);
            emitPulsePoint(b.artTemplate, b.art_onset_auto, b.art_onset,
                b.art_onset_auto, b.art_onset, refArt);
            emitPulsePoint(b.artTemplate, b.art_peak_auto, b.art_peak,
                b.art_onset_auto, b.art_onset, refArt);
            emitPulsePoint(b.artTemplate, b.art_dicrotic_auto, b.art_dicrotic,
                b.art_onset_auto, b.art_onset, refArt);
            emitPulsePoint(b.artTemplate, b.art_peak2_auto, b.art_peak2,
                b.art_onset_auto, b.art_onset, refArt);
            emitPulsePoint(b.artTemplate, b.art_end_auto, b.art_end,
                b.art_onset_auto, b.art_onset, refArt);

            f << ',' << static_cast<int>(b.art_pulm_issue);
            emitPulsePoint(b.artPulmTemplate, b.art_pulm_onset_auto, b.art_pulm_onset,
                b.art_pulm_onset_auto, b.art_pulm_onset, refArtPulm);
            emitPulsePoint(b.artPulmTemplate, b.art_pulm_peak_auto, b.art_pulm_peak,
                b.art_pulm_onset_auto, b.art_pulm_onset, refArtPulm);
            emitPulsePoint(b.artPulmTemplate, b.art_pulm_dicrotic_auto, b.art_pulm_dicrotic,
                b.art_pulm_onset_auto, b.art_pulm_onset, refArtPulm);
            emitPulsePoint(b.artPulmTemplate, b.art_pulm_peak2_auto, b.art_pulm_peak2,
                b.art_pulm_onset_auto, b.art_pulm_onset, refArtPulm);
            emitPulsePoint(b.artPulmTemplate, b.art_pulm_end_auto, b.art_pulm_end,
                b.art_pulm_onset_auto, b.art_pulm_onset, refArtPulm);
        } // end if (wantPulse) pulse point groups

        // Autodetected ECG glyph columns (p_wave, q_onset, r_wave, t_peak).
        // p_wave/q_onset/r_wave ARE this alignment's stored detections --
        // emitted straight, with no parallel recompute that could disagree
        // with them. (r_wave in particular: R has exactly one definition, the
        // alignment's own r_peak. It is never re-derived by an argmax here.)
        // p_wave and t_peak are reactive, bracketed here by the DETECTOR's
        // bars to match this group's "_auto" name -- t_peak by s_end/t_end,
        // which is the bracket reactiveGlyphs uses, not the t_begin this line
        // used to pass and nothing ever set.
        if (wantEcg)
            for (int c = 0; c < 3; ++c) {
                const std::vector<double>& ecg = b.chFor(c, anchor).ecgTemplate_raw;
                const auto aa = b.autoFor(anchor);
                const FeatureMarks::ReactiveEcg rx = FeatureMarks::reactive_ecg(
                    ecg, aa.p_begin[c], aa.q_begin[c],
                    aa.s_end[c], aa.t_end[c], sampleRateHz);
                emitAutoFeatPt(ecg, rx.p_peak);
                emitAutoFeatPt(ecg, aa.q_begin[c]);
                emitAutoFeatPt(ecg, aa.r_peak[c]);
                emitAutoFeatPt(ecg, rx.t_peak);
            }
        if (wantPulse) {
            for (const auto& gl : ppg_and_artpulse_automated_markers) {
                emitAutoFeatPt(b.ppgTemplate, b.*gl.idx);
                if (gl.found) f << ',' << (b.*gl.found ? 1 : 0);
            }
        }

        f << '\n';
    }
}

// Path-taking wrapper, for callers outside the viewer's one-save merge. The
// ostream overload above is what the viewer uses, so a save can build each
// alignment's part in memory instead of staging it through a temp file.
inline void writeTemplateMarkingsCsv(const std::string& path,
    const std::vector<TemplateBin>& bins,
    const std::string& fileID,
    double sampleRateHz,
    AnchorType anchor,
    MarkingsCsvSection section)
{
    std::ofstream f(path);
    if (!f.is_open())
        throw std::runtime_error("cannot open for write: " + path);
    writeTemplateMarkingsCsv(f, bins, fileID, sampleRateHz, anchor, section);
}

inline std::vector<TemplateBin> readTemplateMarkingsBin(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open())
        throw std::runtime_error("cannot open for read: " + path);

    auto r8 = [&]() -> uint8_t {
        uint8_t v = 0; f.read(reinterpret_cast<char*>(&v), 1); return v;
        };
    auto r32 = [&]() -> int {
        int32_t v = 0; f.read(reinterpret_cast<char*>(&v), 4); return v;
        };
    auto r64d = [&]() -> double {
        double v = 0.0; f.read(reinterpret_cast<char*>(&v), 8); return v;
        };

    //read header - it tells you how many bins to expect
    uint64_t n = 0;
    f.read(reinterpret_cast<char*>(&n), 8);

    std::vector<TemplateBin> bins(n);
    for (uint64_t i = 0; i < n; ++i) {
        auto& b = bins[i];
        f.read(reinterpret_cast<char*>(&b.index), 8);

        b.bad_r_ch[0] = (r8() != 0);
        b.bad_r_ch[1] = (r8() != 0);
        b.bad_r_ch[2] = (r8() != 0);
        b.bad_ppg = r8();

        // Mirrors the write order exactly: lead, slot, anchor. Counts and tags
        // are int32; positions are float64.
        //
        // The bins this function returns have NO BANKS -- it reads a marking
        // file, not a template file. slotMarks() therefore creates each slot as
        // it goes, carrying marks and nothing else: no tmpl, no members, so a
        // created slot earns no column and draws nothing. The caller merges
        // these onto the real bins built from _templates.bin, bounds-checking
        // every position against that template's own length.
        for (int lead = 0; lead < 3; ++lead) {
            const int nSlots = r32();
            for (int slot = 0; slot < nSlots; ++slot) {
                const int nAnchors = r32();
                for (int a = 0; a < nAnchors; ++a) {
                    const int tag = r32();
                    tbank::BankMarkerSet& m =
                        b.slotMarks(lead, slot, static_cast<AnchorType>(tag));
                    m.p_begin = r64d(); m.q_begin = r64d();
                    m.s_end = r64d();   m.t_end = r64d();
                    // p_peak deliberately absent: not in the record, not on
                    // the struct. Readers call reactive_ecg on these bars.
                }
            }
        }

        // PPG bars only. The caller calls syncReactivePpg() after the merge to
        // rederive t50 / t80 / t80_rise / pw80 / peak2 from these three plus
        // the auto-detected systolic peak.
        b.ppg_onset = r64d();
        b.ppg_dicrotic = r64d();
        b.ppg_end = r64d();

        b.abp_issue = r8();
        b.abp_onset = r64d(); b.abp_peak = r64d(); b.abp_dicrotic = r64d();
        b.abp_peak2 = r64d(); b.abp_end = r64d();
        b.art_issue = r8();
        b.art_onset = r64d(); b.art_peak = r64d(); b.art_dicrotic = r64d();
        b.art_peak2 = r64d(); b.art_end = r64d();
        b.art_pulm_issue = r8();
        b.art_pulm_onset = r64d(); b.art_pulm_peak = r64d(); b.art_pulm_dicrotic = r64d();
        b.art_pulm_peak2 = r64d(); b.art_pulm_end = r64d();
    }
    return bins;
}