#pragma once
/*The output of <id>_template_markings.csv
*
* ---------------------------------------------------------------------------
* A _user COLUMN EXISTS ONLY WHERE AN OPERATOR VALUE EXISTS
* ---------------------------------------------------------------------------
*
* There is no single "six columns per marker" rule. There are five shapes, and
* which one a column group gets is decided by two questions: is there an
* operator value at all, and if so, does it belong to this block?
*
* 1. ECG POINT, GLYPH -- 3 columns:
*      {name}_ch{c}_y_norm_auto
*      {name}_ch{c}_y_mv_auto
*      {name}_ch{c}_x_ms_auto
*    p_peak, q_peak, r_peak, s_peak, t_peak. markerAtX never hands a glyph out
*    for a drag, so there is no placement to report. p_peak and t_peak ARE
*    re-measured between the operator's bars, and that value is the X on
*    screen -- but it comes from userMarks(), the bar set assembled across all
*    four alignments, so emitting it here put one placed mark into all four
*    blocks, three of them under a waveform it was never compared to.
*
* 2. ECG POINT, BAR IN ITS OWNING BLOCK -- 6 columns: the three above plus
*    _y_norm_user / _y_mv_user / _x_ms_user, interleaved auto-then-user per
*    quantity. P onset in the _P part, Q onset in _Q, J point in _R, T end in
*    _J (anchor_view::owns). In the other three parts that same bar emits
*    shape 1: a P onset measured against the P-aligned average is a different
*    number from one measured against R, and reporting it under another suffix
*    would attribute it to a waveform it was never compared to.
*
* 3. ECG INTERVAL -- 2 columns in the R part, 1 elsewhere:
*      {name}_ch{c}_ms_auto  [, {name}_ch{c}_ms_user]
*    qrs needs q_onset and s_end; qt needs q_onset and t_end -- bars from
*    three different alignments -- so `owns` cannot answer for an interval the
*    way it does for a point. A duration is frame-free, so the user half is
*    reported once, under R, rather than four identical times. No y columns at
*    all: a duration has no amplitude.
*
* 4. PULSE POINT -- 6 columns for a BAR, 3 for a glyph:
*      {name}_y_norm_auto [, _user], _y_mv_auto [, _user], _x_ms_auto [, _user]
*    The PPG channel has THREE bars among its eight columns -- onset,
*    dicrotic, end -- because markerAtX skips PpgPeak, PpgT50, PpgT80 and
*    PpgPeak2. t50, peak, peak2, t80 and t80_rise are auto-only, and their
*    _user half was reactive_ppg re-measured between the three bars they do
*    own, labelled as a placement. The ARTERIAL channels are not like this:
*    all five of their markers are draggable, so all five keep both halves.
*    There is no alignment dimension here -- pulse marks are per bin and the
*    PULSE part is written once -- so problem (2) cannot arise.
*
* 5. AUTODETECTED COMPUTED FEATURE -- 2 columns:
*      {name}_x_ms, {name}_y_mv
*    No auto/user suffix (there is only an auto side) and no normalized
*    column. The derivative landmarks (vpg_u/v/w, apg_a..f, jpg_p1/p2), the
*    p_wave/q_onset/r_wave/t_peak auto features, and the pulse glyph block.
*    Some carry a trailing found-flag column. ppg_pw80_ms_auto is a lone
*    column on the same footing: a width has no position either.
*
* ---------------------------------------------------------------------------
*/

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
#include "fiducial_marker_finding\feature_marks.hpp"
#include "template_morphology_grouping\template_bank.hpp"
#include "fiducial_marker_finding\anchor_view.hpp"

enum class MarkingsCsvSection { EcgOnly, PulseOnly, EcgAndPulse };


// ---------------------------------------------------------------------------
// In-memory model used by the viewer.
// ---------------------------------------------------------------------------

struct ChannelTemplateData {
    std::vector<double> ecgTemplate_raw;
    std::vector<double> ecg_template_raw_iqr;
    std::vector<double> ecgTemplate_squared, ecgTemplate_absval;    // unused by viewer
    double alignment_point_raw = 0, alignment_point_squared = 0;
    double alignment_point_absval = 0;
    int r_col_raw = -1;
    int r_col_squared = -1;
    int r_col_absval = -1;
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
    std::array<tbank::TemplateBank, 3> ecg_bank;
    tbank::TemplateBank ppg_bank;

    uint64_t index = 0;
    std::vector<std::pair<uint64_t, uint64_t>> ppg_bin_indexs;
    std::vector<std::pair<uint64_t, uint64_t>> ecg_bin_indexs;
    bool bad_segment = false;
    uint64_t ch1_n_beats_raw = 0;
    uint64_t ch2_n_beats_raw = 0;
    uint64_t ch3_n_beats_raw = 0;
    uint64_t ppg_n_beats = 0;
    ChannelTemplateData ch1, ch2, ch3;

    // THE OPERATOR'S PER-CHANNEL "Lead Reversed" ANSWER, stamped on every bin
    // by the viewer at load (TemplateViewerWindow::setLeadPolarity) and by
    // template generation when it builds the bins. Carried HERE, on the bin,
    // rather than threaded through alignedLandmarks / ecgDetect / seedSlotBars
    // and the CSV emitters: every one of those already has a TemplateBin and a
    // lead index in hand, so the sign and the channel cannot get out of step.
    // The two functions that take a BARE TRACE with no lead -- computeEcgFeatures
    // and ecgDetectOn -- take an explicit sgn instead, because there is nothing
    // on their inputs to read it from.
    //
    // NOT SERIALIZED. This is a per-RECORD fact duplicated across bins, not
    // per-bin state, and writeTemplateMarkingsBin does not persist it -- the
    // checkbox is re-read from the noise-marking stage on every run, so a copy
    // on disk would be a second source that could disagree with it.
    LeadPolarity polarity;
    std::map<int, std::array<ChannelTemplateData, 3>> anchored;
    struct AnchorAuto {
        double p_begin[3] = { -1, -1, -1 };
        double p_peak[3] = { -1, -1, -1 };
        double q_onset[3] = { -1, -1, -1 };
        double q_peak[3] = { -1, -1, -1 };
        double r_peak[3] = { -1, -1, -1 };
        double s_end[3] = { -1, -1, -1 };
        double t_end[3] = { -1, -1, -1 };
        bool   q_onset_found[3] = { false, false, false };
    };
    std::map<int, std::vector<AnchoredBankSlot> > anchored_bank;
    static int bankSlotKey(int lead, AnchorType a) {
        return static_cast<int>(a) * 4 + lead;
    }
    const AnchoredBankSlot* bankSlotFor(int lead, int slot, AnchorType a) const
    {
        if (lead < 0 || lead > 2) return nullptr;
        if (slot < 0) return nullptr;
        auto it = anchored_bank.find(bankSlotKey(lead, a));
        if (it == anchored_bank.end()) return nullptr;
        if (static_cast<size_t>(slot) >= it->second.size()) return nullptr;
        if (it->second[slot].tmpl.empty()) return nullptr;
        // ALL-NaN IS ABSENT TOO -- see the note above.
        for (double v : it->second[slot].tmpl)
            if (!std::isnan(v)) return &it->second[slot];
        return nullptr;
    }
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
        static const ChannelTemplateData kEmpty{};
        const ChannelTemplateData* base[3] = { &ch1, &ch2, &ch3 };
        if (lead < 0 || lead > 2) return kEmpty;
        if (a == AnchorType::R_PEAK) return *base[lead];
        auto it = anchored.find(static_cast<int>(a));
        if (it == anchored.end()) return kEmpty;
        return it->second[lead];
    }
    int frameShift(int lead, AnchorType from, AnchorType to) const {
        const int a = chFor(lead, from).r_col_raw;
        const int b = chFor(lead, to).r_col_raw;
        return (a < 0 || b < 0) ? 0 : (b - a);
    }
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
        pull(anchor_view::p_begin, &tbank::BankMarkerSet::p_begin);
        pull(anchor_view::q_begin, &tbank::BankMarkerSet::q_onset);
        pull(anchor_view::j_point, &tbank::BankMarkerSet::s_end);
        pull(anchor_view::t_end, &tbank::BankMarkerSet::t_end);
        return out;
    }

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
    // The array bound: EVERY slot, including empty ones. Ask hasVisiblePanel
    // whether a slot is real.
    int slotCount(int lead) const {
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
    double q_onset_auto_ch[3] = { -1, -1, -1 };
    double r_peak_auto_ch[3] = { -1, -1, -1 };
    double s_end_auto_ch[3] = { -1, -1, -1 };
    // (t_begin_auto_ch removed with the marker -- see BankMarkerSet. The
    //  T-peak bracket is s_end/t_end on both the auto and the bar side.)
    double t_end_auto_ch[3] = { -1, -1, -1 };
    double p_begin_auto_ch[3] = { -1, -1, -1 };
    bool q_onset_found_auto_ch[3] = { false, false, false };

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
    // freeze a glyph from. So all fifteen are persisted -- and all fifteen get
    // a _user CSV column, which is why pulseHasUserColumn below names only PPG
    // points.
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
inline std::vector<TemplateBin> binsFromTemplateFile(const template_io::TemplateFile& tf) {
    // EVERY ANCHOR, NOT ONE. This used to take an anchor and project that one
    // block into chN, which is what made a second alignment cost a template
    // regeneration and a window reload.
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

        // Section 4.6 banks. Empty when no bank reached this bin, which is the
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

        // Per-anchor BANK SLOT averages. Absent for a bin the alignment
        // skipped or a slot with no members, in which case bankSlotFor returns
        // nullptr and the caller reports the gap rather than substituting the
        // slot's R-aligned average.
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

// Path-taking wrapper, for callers that genuinely have a file: the standalone
// template_marking tool, and a reload of a finished record. The anchor
// parameter is vestigial -- nothing left to select -- and is kept only so
// existing call sites compile.
inline std::vector<TemplateBin> readTemplateInfoBin(const std::string& path,
    AnchorType /*anchor*/ = AnchorType::R_PEAK) {
    return binsFromTemplateFile(template_io::read_template_binfile(path));
}

// ---------------------------------------------------------------------------
// template_markings.bin layout -- BARS ONLY, POSITIONS AS FLOAT64:
//
//   header:
//     uint32  magic   = kMarkMagic ("TMRK")
//     uint32  version = kMarkVersion (1)
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
//             float64 p_begin, q_onset, s_end, t_end
//
//     -- PPG bars, PER SLOT (same shape as the ECG block above): --
//     int32   slotCount
//       per slot:
//         float64 onset, dicrotic, end
//
//     -- arterial, one block each for ABP, ART, ART_PULM: --
//     uint8   <chan>_issue
//     float64 <chan>_onset, _peak, _dicrotic, _peak2, _end
//
// Positions use -1 as the "unmarked / not applicable" sentinel; test `< 0`.
// Structural fields (counts, tags, index, issue flags) stay integral.
//
// THE SET WRITTEN HERE IS THE SET THAT EARNS A _user CSV COLUMN. Four ECG
// bars per (lead, slot, anchor), three PPG bars, five per arterial channel --
// and nothing else, because nothing else is placed by hand. If a landmark is
// absent from this record it is auto-only, and the CSV writer below must not
// emit a user half for it.
//
// THERE IS NOW A MAGIC AND A VERSION, which is what this comment used to say
// was missing. The header was a bare bin count, so a file written under any
// earlier layout misparsed from its first marker set onward with no error --
// it simply read a bin count out of whatever the old first field was. The
// magic catches that: a pre-version file fails the check and is reported
// rather than silently misread.
//
// VERSION 1 IS THE FIRST AND ONLY VERSION. Nothing has shipped, so there is no
// migration path and none is wanted: a file whose version is not 1 is
// rejected. When a second version arrives, branch in readTemplateMarkingsBin
// on the field rather than adding another unversioned section.
//
// PPG BARS ARE PER SLOT, as of version 1. They used to be three bin-level
// values, which meant a bin's morphology columns all shared one set of pulse
// bars while each column DREW its own pulse (ppg_bank.templates[slot].tmpl) --
// so the foot could sit nowhere near a minimum on every column but one. The
// bars now live where the waveform does, in
// ppg_bank.templates[slot].pulse_marks, for slot 0 exactly as for the rest.
// ---------------------------------------------------------------------------
// "TMRK", little-endian. Any 4 bytes would do; these are readable in a hex
// dump, which is worth something the first time a file will not open.
inline constexpr uint32_t kMarkMagic = 0x4B524D54u;
inline constexpr uint32_t kMarkVersion = 1u;

inline void writeTemplateMarkingsBin(const std::string& path,
    const std::vector<TemplateBin>& bins) {
    std::ofstream f(path, std::ios::binary);
    if (!f.is_open())
        throw std::runtime_error("cannot open for write: " + path);

    // header: magic, version, bin count
    uint32_t magic = kMarkMagic;
    uint32_t ver = kMarkVersion;
    f.write(reinterpret_cast<const char*>(&magic), 4);
    f.write(reinterpret_cast<const char*>(&ver), 4);
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
        //         float64  p_begin, q_onset, s_end, t_end
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
            const int nSlots = b.slotCount(lead);
            w32(nSlots);
            for (int slot = 0; slot < nSlots; ++slot) {
                const auto& byAnchor =
                    b.ecg_bank[lead].templates[slot].markers_by_anchor;
                w32(static_cast<int>(byAnchor.size()));
                for (const auto& kv : byAnchor) {
                    w32(kv.first);                      // AnchorType tag
                    const tbank::BankMarkerSet& m = kv.second;
                    w64d(m.p_begin); w64d(m.q_onset);
                    w64d(m.s_end);   w64d(m.t_end);
                }
            }
        }

        // PPG: THE THREE BARS ONLY, PER SLOT. t50 / peak / peak2 / t80 /
        // t80_rise / pw80 are auto-only glyphs -- markerAtX refuses to hand any
        // of them out -- and all five reactive ones come back from reactive_ppg
        // bracketed by these bars, so storing them would be storing a cache.
        //
        // PER SLOT, because each morphology column draws its OWN pulse
        // (ppg_bank.templates[slot].tmpl). Three bin-level values meant every
        // column shared one set of bars against a different waveform each.
        // Slot 0 is a slot like any other here.
        {
            const int nPulse = static_cast<int>(b.ppg_bank.templates.size());
            w32(nPulse);
            for (int slot = 0; slot < nPulse; ++slot) {
                const tbank::BankPulseMarkerSet& pm =
                    b.ppg_bank.templates[slot].pulse_marks;
                w64d(pm.onset);
                w64d(pm.dicrotic);
                w64d(pm.end);
            }
        }

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
inline EcgFeatures computeEcgFeatures(const std::vector<double>& ecg, double p_peak, double q_onset, double r_peak, double s_end, double t_end,
    double rateHz, double sgn, curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto)
{
    EcgFeatures f;
    const double N = static_cast<double>(ecg.size());
    const double msPerSamp = (rateHz > 0.0) ? 1000.0 / rateHz : NAN;

    if (q_onset >= 0.0 && s_end >= q_onset) f.qrs_ms = (s_end - q_onset) * msPerSamp;
    if (q_onset >= 0.0 && t_end >= q_onset) f.qt_ms = (t_end - q_onset) * msPerSamp;

    if (r_peak >= 0.0 && r_peak <= N - 1.0) f.r_idx = r_peak;
    const int rInt = (r_peak >= 0.0) ? static_cast<int>(std::lround(r_peak)) : -1;
    // Q for the q_peak column: same canonical finder compute_q_onset uses.
    // Mirrors compute_s_peak's signature below. -1 when there is no Q trough.
    f.q_idx = FeatureMarks::find_q_peak(ecg, rInt, rateHz, 1.0, peakMode);    // sub-sample
    // S for |R|+|S| = first opposite-polarity trough after R (robust; not the
    // max over [R, s_end], which depends on where s_end sits).
    f.s_idx = FeatureMarks::find_s_peak(ecg, rInt, rateHz, 1.0, peakMode);   // sub-sample
    return f;
}

// (EcgColSpec / ecgCols deleted: declared, never referenced, and wrong --
//  it claimed all six ECG markers were user-placed, which was never true and
//  is emphatically not true now that four bars and a set of glyphs are
//  distinguished. The live column order is ecgPointNames, inside
//  writeTemplateMarkingsCsv.)
inline constexpr const char* ppgCols[] = { "ppg_onset","ppg_t50","ppg_peak","ppg_dicr","ppg_peak2","ppg_t80","ppg_t80_rise","ppg_end" };
inline constexpr const char* abpCols[] = { "abp_onset","abp_peak","abp_dicr","abp_peak2","abp_end" };
inline constexpr const char* artCols[] = { "art_onset","art_peak","art_dicr","art_peak2","art_end" };
inline constexpr const char* artPulmCols[] = { "art_pulm_onset","art_pulm_peak","art_pulm_dicr","art_pulm_peak2","art_pulm_end" };


// ---- WHICH PULSE COLUMNS HAVE AN OPERATOR VALUE AT ALL --------------------
//
// markerAtX's skip list is the authority: it refuses to hand out PpgPeak,
// PpgT50, PpgT80 and PpgPeak2, so FIVE of the eight PPG columns are auto-only
// and their _user half was never a placement. It was reactive_ppg re-measured
// between the three bars the operator does own (onset, dicrotic, end),
// labelled as though someone had put it there -- the same mislabelling the
// ECG glyph columns had.
//
// THE ARTERIAL CHANNELS ARE GENUINELY DIFFERENT. Nothing in the Abp* / Art* /
// ArtPulm* range is skipped, because those channels have no detector to freeze
// a glyph from, so all five of each are draggable and all five keep both
// halves. Hence this predicate names PPG points only.
//
// KEEP IN STEP WITH markerAtX AND WITH writeTemplateMarkingsBin. A name listed
// here that the widget still hands out loses a real column; one omitted that
// the widget skips fabricates one. The .bin is the cross-check: exactly the
// positions written there are the ones that earn a user column.
inline bool pulseHasUserColumn(const char* pointName) {
    static const char* kAutoOnly[] = {
        "ppg_t50", "ppg_peak", "ppg_peak2", "ppg_t80", "ppg_t80_rise"
    };
    for (const char* a : kAutoOnly)
        if (std::strcmp(pointName, a) == 0) return false;
    return true;
}


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

// ONE SOURCE FOR "DOES THIS TEMPLATE HAVE A PANEL ON SCREEN", shared by
// leadsForBinTemplate and the markings CSV so the row set IS the column set.

// THIS ALIGNMENT'S LANDMARKS FOR ONE CHANNEL, detected on the spot.
//
// Replaces the auto_by_anchor cache: it held the result of running
// detect_template_landmarks over a stored template -- a pure function of data
// already in memory -- and autoFor() returned the R-aligned field set on a
// miss, so a caller asking for P silently got R. Recomputing removes both the
// staleness and the fallback. An absent alignment yields valid=false.
inline FeatureMarks::TemplateLandmarks alignedLandmarks(
    const TemplateBin& b, int lead, AnchorType a, double sampleRate,
    curve_fit::FitMode fitMode = curve_fit::FitMode::Auto,
    curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto)
{
    FeatureMarks::TemplateLandmarks lm;
    if (lead < 0 || lead > 2) return lm;
    const ChannelTemplateData& cd = b.chFor(lead, a);
    if (cd.ecgTemplate_raw.empty() || cd.r_col_raw < 0) return lm;
    return FeatureMarks::detect_template_landmarks(cd.ecgTemplate_raw,
        static_cast<int>(std::lround(cd.r_col_raw)), sampleRate,
        b.polarity.sign(lead), fitMode, peakMode);
}

// ---------------------------------------------------------------------------
// PEAK PLACEMENT -- the one conversion from a detector seed to a drawn,
// focused and written position. The glyph, the focus mark and the CSV all go
// through this, or they are three numbers for one landmark.
//
// Returns the fit WHOLE: position and provenance together, so a caller cannot
// report the model from here and the location from elsewhere.
//
// seedPos is the detector's answer (detect_template_landmarks / ecgDetect): it
// decides where to look. The returned position is the placement.
//
// This makes the consumers agree, not correct -- the seed's bracket search
// cannot distinguish an inverted monophasic R from an upright biphasic one, so
// where the seed is wrong all three are wrong together.
// ---------------------------------------------------------------------------
enum class EcgPeak { P, Q, R, S, T };

// Per-landmark fit windows, resolved in one place so no call site carries its
// own copy of the mapping.
inline double peakSigmaFor(EcgPeak w) {
    switch (w) {
    case EcgPeak::P: return upsample_for_fit::peak_sigma::P;
    case EcgPeak::Q: return upsample_for_fit::peak_sigma::Q;
    case EcgPeak::S: return upsample_for_fit::peak_sigma::S;
    case EcgPeak::T: return upsample_for_fit::peak_sigma::T;
    default:         return upsample_for_fit::peak_sigma::R;
    }
}
inline int peakHalfWidthFor(EcgPeak w) {
    switch (w) {
    case EcgPeak::P: return upsample_for_fit::peak_halfwidth::P;
    case EcgPeak::Q: return upsample_for_fit::peak_halfwidth::Q;
    case EcgPeak::S: return upsample_for_fit::peak_halfwidth::S;
    case EcgPeak::T: return upsample_for_fit::peak_halfwidth::T;
    default:         return upsample_for_fit::peak_halfwidth::R;
    }
}

// The contest, once. Callers that only need the position use placeEcgPeak
// below; the focus panel needs the curves too and takes the whole thing.
inline upsample_for_fit::PeakCandidates ecgPeakCandidates(
    const std::vector<double>& tmpl, EcgPeak which, double seedPos,
    curve_fit::PeakFitMode mode)
{
    upsample_for_fit::PeakCandidates none;
    const int n = static_cast<int>(tmpl.size());
    // Absent stays absent: -1 means absent everywhere downstream, and clamping
    // to column 0 would make a missing landmark look like one found at an edge.
    if (n < 5 || !(seedPos >= 0.0) || seedPos > static_cast<double>(n - 1))
        return none;
    const int seed = std::clamp(static_cast<int>(std::lround(seedPos)), 0, n - 1);
    // max(3, ...) so a bad table entry cannot give a degenerate fit window.
    return upsample_for_fit::peakCandidates(
        tmpl, seed, peakSigmaFor(which),
        std::max(3, peakHalfWidthFor(which)), mode);
}

inline upsample_for_fit::peak_fit placeEcgPeak(
    const std::vector<double>& tmpl, EcgPeak which, double seedPos,
    curve_fit::PeakFitMode mode)
{
    const upsample_for_fit::PeakCandidates pc =
        ecgPeakCandidates(tmpl, which, seedPos, mode);
    upsample_for_fit::peak_fit out;
    if (!pc.valid || pc.winner < 0) return out;   // position stays -1 = absent
    out = pc.draw[pc.winner];
    out.position = pc.placement;   // the guarded contest's answer, not the
    // unguarded draw fit's own vertex
    return out;
}


// ONE ACCESSOR FOR "THIS SLOT, THIS ALIGNMENT". The waveform and its R column
// travel together: handing out one without the other is how a bar came to hold
// a number measured in a different alignment's frame.
struct SlotView {
    const std::vector<double>* tmpl = nullptr;   // this alignment's average
    const std::vector<double>* iqr = nullptr;    // its spread, same source
    int  r_col = -1;                             // this alignment's R column
    bool valid = false;
};

inline SlotView slotView(const TemplateBin& b, int lead, int slot, AnchorType a)
{
    SlotView v;
    if (lead < 0 || lead > 2 || slot < 0) return v;

    // NO SPECIAL CASE FOR R. It used to return the slot's BankTemplate::tmpl
    // -- that slot averaged over ALL its members -- while leadsForBinTemplate
    // and the focus panel draw bankSlotFor(.., R)->tmpl, averaged over
    // members_clean. Two different populations, so the detector measured one
    // waveform and the glyph was painted over another: the X sat at its own
    // array's apex, several samples off the apex of the trace beneath it.
    //
    // prepareViewerJob aligns all four anchors INCLUDING R, so the R entry
    // exists and no fallback is needed. A null here is a writer gap, which the
    // caller reports -- not a reason to substitute a different average.

    // Row subsets of the bin's aligned matrix, so the bin's per-anchor r_col
    // is their R column.
    const AnchoredBankSlot* asl = b.bankSlotFor(lead, slot, a);
    if (!asl) return v;
    v.tmpl = &asl->tmpl;
    v.iqr = &asl->tmpl_iqr;
    v.r_col = static_cast<int>(std::lround(b.chFor(lead, a).r_col_raw));
    v.valid = (v.r_col >= 0);
    return v;
}

// Does this ECG channel exist: the raw per-channel template, the same test
// leadsForBin uses for max_leads. NOT chFor(lead, anchor), which returns a
// sized entry for all three channels under P/Q/J whether they hold data or not.
inline bool ecgChannelPresent(const TemplateBin& b, int lead) {
    if (lead < 0 || lead > 2) return false;
    return !b.chFor(lead, AnchorType::R_PEAK).ecgTemplate_raw.empty();
}

// TRUE when the grid draws a panel for this (channel, slot) in this alignment.
// Mirrors leadsForBinTemplate's two accepting branches: the bank branch, and
// the pre-bank slot-0 fallback for a file with no bank at all.
//
// pulseThin is part of it on purpose -- a thin PPG cohort suppresses the ECG
// panel, and it gates the fallback too, which is the defect that let slot 0
// keep appearing however tightly the callers were gated.
inline bool hasVisiblePanel(const TemplateBin& b, int lead, int slot,
    AnchorType gridAnchor)
{
    // No fallbacks: a slot has a bank template with an average for this
    // alignment, or it has no panel.
    if (lead < 0 || lead > 2 || slot < 0) return false;
    if (!ecgChannelPresent(b, lead)) return false;

    const tbank::TemplateBank& bank = b.ecg_bank[lead];
    if (slot >= bank.size()) return false;              // ragged: shorter bank
    const tbank::BankTemplate& tp = bank.templates[slot];
    if (tp.tmpl.empty()) return false;
    if (tp.tooFewBeats(/*is_ppg=*/false)) return false;
    // EVERY SLOT, INCLUDING 0. This was `slot != 0 && ...`, exempting slot 0
    // from the predicate that governed every other column -- and
    // leadsForBinTemplate carried the same exemption, so the two had to be
    // removed together or the page would count a column it could not draw.
    //
    // REPORTED, because a rejected slot 0 is a bin losing its dominant
    // morphology from the grid AND from the markings CSV, and that must not be
    // silent. seed_pool seeds slot 0 from the clean pool, so a non-REGULAR
    // slot 0 means the seed pool was contaminated -- a fact worth seeing
    // rather than a column worth faking.
    if (!tp.wantsLandmarkMarking()) {
        if (slot == 0)
            fprintf(stderr, "[visible] bin %llu lead %d slot 0 is not REGULAR"
                " (presumed category %d) -- NO _A COLUMN for this bin\n",
                (unsigned long long)b.index, lead,
                static_cast<int>(tp.presumedCategory()));
        return false;
    }
    // A thin pulse cohort suppresses the ECG panel for that slot.
    if (slot < b.ppg_bank.size()
        && b.ppg_bank.templates[slot].tooFewBeats(/*is_ppg=*/true))
        return false;

    return slotView(b, lead, slot, gridAnchor).valid;
}

// The slots any channel draws, in slot order -- the column set for one bin.
inline std::vector<int> visibleSlots(const TemplateBin& b, AnchorType gridAnchor)
{
    std::vector<int> out;
    for (int slot = 0; slot < tbank::max_templates_per_bin * 4; ++slot)
        for (int lead = 0; lead < 3; ++lead)
            if (hasVisiblePanel(b, lead, slot, gridAnchor)) {
                out.push_back(slot);
                break;
            }
    return out;
}

// ===========================================================================
// THE ONLY PLACE ECG LANDMARKS ARE COMPUTED.
//
// Three callers used to assemble the detector's arguments themselves -- the bar
// seeding, the focus panel, and BinPlotWidget::reactiveGlyphs -- and every P
// onset bug today was one of those assemblies differing from another: the
// bin-wide average instead of the slot's, R's column instead of the alignment's,
// a hardcoded Auto instead of the live radio, the /ref-normalised array instead
// of the raw one. Matching them by hand does not hold, because nothing stops
// the next one drifting.
//
// So the arguments are assembled ONCE, here. Callers pass what identifies the
// template -- (bin, lead, slot, alignment) -- and the bars, and get positions
// back. slotView supplies the trace and its R column together, which is what
// makes the wrong pairing unrepresentable.
// ===========================================================================
struct EcgFiducials {
    // Every ECG landmark, in the alignment's own columns. -1 = absent.
    double p_begin = -1.0, p_peak = -1.0;
    double q_onset = -1.0, q_peak = -1.0;
    double r_peak = -1.0;
    double s_peak = -1.0, s_end = -1.0;
    double t_peak = -1.0, t_end = -1.0;
    bool   q_onset_found = false;
    bool   valid = false;
};

// SPLIT ALONG THE LINE THE BARS DRAW. ecgDetect is a pure function of the
// trace, so a caller may cache it; ecgFiducialsFrom is the two bar-bracketed
// peaks and is cheap enough to call per repaint. ecgFiducials is the two back
// to back, with the original signature.
struct EcgDetection {
    FeatureMarks::TemplateLandmarks lm;
    double s_peak = -1.0;
    // The trace the landmarks were measured on; the reactive half re-brackets
    // on the SAME array. Non-owning and interior to the bin, so a holder must
    // drop the cache when the bin or the slot changes.
    const std::vector<double>* tmpl = nullptr;
    bool valid = false;
};

// TRACE GIVEN EXPLICITLY. A caller that already holds the waveform on screen --
// BinPlotWidget does, in m_ecg, notch-filtered and amplitude-scaled -- must
// detect on THAT array, or its glyphs describe a different signal from the one
// under them. `tmpl` has to outlive the returned EcgDetection, which holds a
// pointer to it for the reactive half.
inline EcgDetection ecgDetectOn(const std::vector<double>& tmpl, int r_col,
    double sampleRate, double sgn,
    curve_fit::FitMode onOffsetMode = curve_fit::FitMode::Auto,
    curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto)
{
    EcgDetection d;
    if (tmpl.size() < 3 || r_col < 0) return d;

    d.lm = FeatureMarks::detect_template_landmarks(tmpl, r_col, sampleRate, sgn,
        onOffsetMode, peakMode);
    if (!d.lm.valid) return d;

    d.s_peak = FeatureMarks::find_s_peak(tmpl, r_col, sampleRate, sgn, peakMode);
    d.tmpl = &tmpl;
    d.valid = true;
    return d;
}

inline EcgDetection ecgDetect(const TemplateBin& b, int lead, int slot,
    AnchorType a, double sampleRate,
    curve_fit::FitMode onOffsetMode = curve_fit::FitMode::Auto,
    curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto)
{
    EcgDetection d;
    const SlotView sv = slotView(b, lead, slot, a);
    if (!sv.valid) return d;

    d.lm = FeatureMarks::detect_template_landmarks(*sv.tmpl, sv.r_col,
        sampleRate, b.polarity.sign(lead), onOffsetMode, peakMode);
    if (!d.lm.valid) return d;

    // S peak has no field on TemplateLandmarks; it is the same finder the
    // interval code uses, on this alignment's trace and R column.
    d.s_peak = FeatureMarks::find_s_peak(*sv.tmpl, sv.r_col, sampleRate,
        b.polarity.sign(lead), peakMode);
    d.tmpl = sv.tmpl;
    d.valid = true;
    return d;
}

// `bars` is the REACTIVE input: p_peak and t_peak are re-measured between the
// operator's P-onset/Q-onset and S-end/T-end bars, so they follow a drag. Pass
// a default-constructed set to get the detector's own brackets instead.
//
// EACH PEAK TESTS ITS OWN TWO BRACKETS. A per-anchor marker set holds only the
// bars its alignment shows (anchor_view::showsBar), so a half-populated set
// must degrade one peak rather than poison both -- asking bars.isUnset() once
// did the latter, and compute_p_peak clamps a -1 bracket to the trace edge
// instead of reporting absence.
inline EcgFiducials ecgFiducialsFrom(const EcgDetection& d, double sampleRate,
    curve_fit::PeakFitMode peakMode, const tbank::BankMarkerSet& bars)
{
    EcgFiducials out;
    if (!d.valid || !d.tmpl) return out;
    const FeatureMarks::TemplateLandmarks& lm = d.lm;

    out.p_begin = lm.p_begin;
    out.q_onset = lm.q_onset;
    out.q_onset_found = lm.q_onset_found;
    out.q_peak = lm.q_peak;
    out.r_peak = lm.r_peak;
    out.s_end = lm.s_end;
    out.t_end = lm.t_end;
    out.s_peak = d.s_peak;

    // P AND T PEAK ARE BRACKET-DERIVED, not detected on their own.
    // TemplateLandmarks has no t_peak field at all -- reactive_ecg is the one
    // function that measures both, between the bars that bracket them.
    const bool haveP = (bars.p_begin >= 0.0 && bars.q_onset >= 0.0);
    const bool haveT = (bars.s_end >= 0.0 && bars.t_end >= 0.0);
    const FeatureMarks::ReactiveEcg rx = FeatureMarks::reactive_ecg(*d.tmpl,
        haveP ? bars.p_begin : lm.p_begin,
        haveP ? bars.q_onset : lm.q_onset,
        haveT ? bars.s_end : lm.s_end,
        haveT ? bars.t_end : lm.t_end,
        sampleRate, peakMode);
    out.p_peak = rx.p_peak;
    out.t_peak = rx.t_peak;

    out.valid = true;
    return out;
}

inline EcgFiducials ecgFiducials(const TemplateBin& b, int lead, int slot,
    AnchorType a, double sampleRate,
    curve_fit::FitMode onOffsetMode = curve_fit::FitMode::Auto,
    curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto,
    const tbank::BankMarkerSet& bars = tbank::BankMarkerSet{})
{
    return ecgFiducialsFrom(
        ecgDetect(b, lead, slot, a, sampleRate, onOffsetMode, peakMode),
        sampleRate, peakMode, bars);
}

// SEEDING, with the same one assembly. seed_bank_template owns the
// admissibility mask (it includes landmark_admissibility.hpp, which the viewer
// cannot -- the graph would cycle), so the bars go through it; slotView still
// supplies the trace and R column, so the inputs match ecgFiducials exactly.
inline bool seedSlotBars(TemplateBin& b, int lead, int slot, AnchorType a,
    double sampleRate,
    curve_fit::FitMode onOffsetMode = curve_fit::FitMode::Auto,
    curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto)
{
    const SlotView sv = slotView(b, lead, slot, a);
    if (!sv.valid) return false;
    FeatureMarks::seed_bank_template(*sv.tmpl, sv.r_col, sampleRate,
        b.polarity.sign(lead), a,
        b.slotMarks(lead, slot, a), onOffsetMode, peakMode);
    return true;
}

// PQRST_A and friends: the name the grid shows, minus its "Ch1 " prefix.
// Letter from tbank::letterRanks, the same function leadsForBinTemplate uses --
// the raw slot index skips letters when a lower slot is empty.
inline std::string bankSlotName(const tbank::TemplateBank& bank, int slot) {
    if (slot < 0 || static_cast<size_t>(slot) >= bank.templates.size())
        return "none";
    const uint8_t code = bank.templates[slot].label_code;
    std::string cls = "PQRST";
    if (code != tbank::kUnlabeled) {
        switch (code) {
        case tbank::kCodePvc:        cls = "PVC";   break;
        case tbank::kCodePac:        cls = "PAC";   break;
        case tbank::kCodeVt:         cls = "VT";    break;
        case tbank::kCodeMinorNoise: cls = "NOISE"; break;
        default:                     cls = "CODE" + std::to_string(code); break;
        }
    }
    int letterIdx = 0;
    const std::vector<uint8_t> letters = tbank::letterRanks(bank);
    if (static_cast<size_t>(slot) < letters.size())
        letterIdx = letters[static_cast<size_t>(slot)];
    return cls + "_" + std::string(1, static_cast<char>('A' + (letterIdx % 26)));
}

inline void writeTemplateMarkingsCsv(std::ostream& f,
    const std::vector<TemplateBin>& bins,
    const std::string& fileID,
    double sampleRateHz,
    AnchorType anchor,
    MarkingsCsvSection section,
    curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto,
    curve_fit::FitMode fitMode = curve_fit::FitMode::Auto)
{
    const bool wantEcg = (section == MarkingsCsvSection::EcgOnly
        || section == MarkingsCsvSection::EcgAndPulse);
    const bool wantPulse = (section == MarkingsCsvSection::PulseOnly
        || section == MarkingsCsvSection::EcgAndPulse);

    // WHICH ALIGNMENTS THIS CALL EMITS. EcgAndPulse walks all four and
    // suffixes each block itself; the single-alignment modes emit just the one
    // asked for, unsuffixed, exactly as before. `anchor` is unused for
    // EcgAndPulse.
    std::vector<AnchorType> anchors;
    if (section == MarkingsCsvSection::EcgAndPulse)
        anchors.assign(anchor_view::anchor_array.begin(),
            anchor_view::anchor_array.end());
    else
        anchors.push_back(anchor);
    const bool suffixed = (section == MarkingsCsvSection::EcgAndPulse);
    auto sfxFor = [&](AnchorType a) {
        return suffixed ? (std::string("_") + anchor_view::label(a))
            : std::string();
        };

    // The R block is the one that carries the interval user columns. See
    // emitIntervalHeader. Reassigned per alignment in both loops below;
    // emitIntervalPair captures it by reference.
    auto intervalUserIn = [](AnchorType a) { return a == AnchorType::R_PEAK; };
    bool intervalUserHere = intervalUserIn(anchors.front());


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
                const FeatureMarks::TemplateLandmarks aaR =
                    alignedLandmarks(b, c, AnchorType::R_PEAK, sampleRateHz,
                        fitMode, peakMode);
                if (!aaR.valid) continue;
                EcgFeatures ft = computeEcgFeatures(ecg, aaR.p_peak, aaR.q_onset, aaR.r_peak, aaR.s_end, aaR.t_end, sampleRateHz, b.polarity.sign(c), peakMode);
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
    // ONE ROW PER (bin, channel, template slot) -- one row per panel on
    // screen. channel is a row key, not a column suffix: the template name is
    // only unique within (bin, channel), so CH1 and CH2 can each hold a
    // PQRST_A and they are different templates. n_members is memberCount;
    // n_clean is cleanCount, the one the visibility predicates test.
    // bad_ecg / bad_ppg, 1 = yes and 0 = no, INDEPENDENT. The right-click
    // cycle can set either or both (Good -> BAD ECG -> BAD PPG -> both), so a
    // panel marked both bad reads 1,1. They were bad_r and ppg_issue, which
    // named the mechanism rather than the verdict and read as alternatives.
    //
    // bad_ppg is 1 when the pulse is absent as well as when it was marked bad:
    // the file says whether there is a usable pulse, not why there isn't.
    // rr_interval_ms is the template's OWN mean R-R, not the bin's: it is
    // averaged over this template's member slices, so a PVC column and its
    // sinus sibling in the same bin report different intervals, which is the
    // point of having it per row. Blank when no interval was measured --
    // a reloaded pre-v5 bank, or a template whose slices carried no RR -- and
    // never 0, which would read as a measured zero-length cycle.
    f << "file_id,bin_index,channel,template,bad_ecg,bad_ppg,n_members,n_clean"
        ",rr_interval_ms";

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
        "p_begin", "p_peak", "q_onset", "q_peak", "r_peak", "s_peak",
        "s_end",   "t_peak", "t_end"
    };
    static const char* ecgIntervalNames[] = { "qrs", "qt" };

    auto emitEcgPointHeader = [&](const char* name,
        AnchorType anchor, const std::string& sfx) {
            // A _user triple exists only for a BAR, and only in the block of the
            // alignment that owns it: P onset under _P, Q onset under _Q, J point
            // under _R, T end under _J. Glyphs get the auto triple alone -- there
            // is no operator value to report, because markerAtX never hands a
            // glyph out for a drag. p_peak and t_peak ARE re-measured between the
            // operator's bars, and that value is the X on screen, but it is
            // derived from userMarks() -- the bar set assembled across all four
            // alignments -- so emitting it here put one placed mark into all four
            // blocks, three of them under a waveform it was never compared to.
            //
            // ONE SOURCE for this rule: the row loop asks anchor_view too, with
            // the same two arguments, so the header and the body cannot disagree
            // about column count.
            const bool userToo = anchor_view::hasUserColumn(name, anchor);
            f << ',' << name << "_y_norm_auto" << sfx;
            if (userToo) f << ',' << name << "_y_norm_user" << sfx;
            f << ',' << name << "_y_mv_auto" << sfx;
            if (userToo) f << ',' << name << "_y_mv_user" << sfx;
            f << ',' << name << "_x_ms_auto" << sfx;
            if (userToo) f << ',' << name << "_x_ms_user" << sfx;
        };
    auto emitIntervalHeader = [&](const char* name,
        AnchorType anchor, const std::string& sfx) {
            f << ',' << name << "_ms_auto" << sfx;
            // NO OWNING BLOCK. qrs needs q_onset and s_end, qt needs q_onset and
            // t_end -- bars from three different alignments -- so `owns` cannot
            // answer for an interval the way it does for a point. A duration is
            // frame-free, so the one copy under R is the whole answer rather than
            // an arbitrary pick, and the other three blocks carried identical
            // numbers.
            if (intervalUserIn(anchor))
                f << ',' << name << "_ms_user" << sfx;
        };
    // Pulse: 6 cols for a BAR, 3 for a glyph. Five of the eight PPG points are
    // auto-only (see pulseHasUserColumn); all five on each arterial channel
    // are bars. No alignment dimension here.
    auto emitPulsePointHeader = [&](const char* name) {
        const bool userToo = pulseHasUserColumn(name);
        f << ',' << name << "_y_norm_auto";
        if (userToo) f << ',' << name << "_y_norm_user";
        f << ',' << name << "_y_mv_auto";
        if (userToo) f << ',' << name << "_y_mv_user";
        f << ',' << name << "_x_ms_auto";
        if (userToo) f << ',' << name << "_x_ms_user";
        };
    // Autodetected computed features (no user bar; derived from AUTO markers).
    auto emitAutoFeatHeader = [&](const char* name, const std::string& sfx = {}) {
        f << ',' << name << "_x_ms" << sfx << ',' << name << "_y_mv" << sfx;
        };

    if (wantEcg) {
        for (AnchorType a : anchors) {
            const std::string sfx = sfxFor(a);
            for (const char* n : ecgPointNames)     emitEcgPointHeader(n, a, sfx);
            // SAME LOOP as the points, so the two interval columns land
            // immediately after the nine point groups and the row loop below can
            // mirror the order without a second pass.
            for (const char* n : ecgIntervalNames)  emitIntervalHeader(n, a, sfx);
            for (const char* g : { "p_wave_auto", "q_onset_auto",
                                   "r_wave_auto", "t_peak_auto" })
                emitAutoFeatHeader(g, sfx);
        }
    }

    if (wantPulse) {
        // No ppg_issue header here any more: the pulse verdict is the bad_ppg
        // ROW KEY above. Emitting it here too would have left a column with no
        // value and shifted every pulse column by one.
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
        // (vpg_u/v/w ARE NOT EMITTED AGAIN HERE. They are already in
        //  ppg_and_artpulse_automated_markers, which the loop above walks and
        //  which the row loop walks exactly once -- so this second pass added
        //  six header names that no row ever filled, leaving every data row
        //  six fields short of the header it was written under.)
    }
    f << '\n';

    // ---- row loop ----------------------------------------------------------
    const double toMs = (sampleRateHz > 0.0) ? 1000.0 / sampleRateHz : 1.0;
    // Default 6 significant figures would drop the sub-sample fraction on
    // 4-digit millisecond values.
    f << std::setprecision(10);

    // Emit one ECG point group: normalized, raw, x_ms -- auto always, user
    // only when `userToo`. Any missing piece leaves that field blank.
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

    // Emit one pulse point group. TAKES ITS COLUMN NAME, so the user half is
    // gated by the SAME pulseHasUserColumn call the header emitter made --
    // there is no second rule to keep in step, exactly as on the ECG side.
    // footIdx_auto / footIdx_user are the "onset" indices for their respective
    // variants (used to compute the local ratio (y - foot)/foot). ref is the
    // channel's global PI median.
    //
    // The amplitude at a fractional position is INTERPOLATED (FeatureMarks::
    // sample_at) rather than read from a rounded column, and the millisecond
    // column carries the fraction, so T80 is no longer quantised to 3.9 ms at
    // 256 Hz -- which matters because Section 6.3's T80 entropy result turns on
    // small differences in exactly that interval.
    auto emitPulsePoint = [&](const char* name, const std::vector<double>& v,
        double idx_auto, double idx_user,
        double foot_auto, double foot_user, double ref)
        {
            const bool userToo = pulseHasUserColumn(name);
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
            const double y_u = userToo ? y_of(idx_user) : std::nan("");
            const double n_a = normOf(y_a, foot_auto);
            const double n_u = userToo ? normOf(y_u, foot_user) : std::nan("");

            f << ',';   if (std::isfinite(n_a)) f << n_a;
            if (userToo) { f << ','; if (std::isfinite(n_u)) f << n_u; }
            f << ',';   if (std::isfinite(y_a)) f << y_a;
            if (userToo) { f << ','; if (std::isfinite(y_u)) f << y_u; }
            f << ',';   if (idx_auto >= 0.0)    f << (idx_auto * toMs);
            if (userToo) { f << ','; if (idx_user >= 0.0) f << (idx_user * toMs); }
        };

    // Interval pair. The user half only in the R block -- see
    // emitIntervalHeader for why an interval has no owning alignment.
    auto emitIntervalPair = [&](double auto_ms, double user_ms) {
        f << ',';   if (std::isfinite(auto_ms)) f << auto_ms;
        if (intervalUserHere) { f << ','; if (std::isfinite(user_ms)) f << user_ms; }
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
        for (int c = 0; c < 3; ++c) {
            static const char* kChan[3] = { "CH1", "CH2", "CH3" };
            const tbank::TemplateBank& bank = b.ecg_bank[c];
            // SAME RANGE the grid scans (markingSlotsForBin walks
            // max_templates_per_bin * 4), not bank.size(): slot indices are not
            // dense, and hasVisiblePanel is what decides, not the bound.
            for (int slot = 0; slot < tbank::max_templates_per_bin * 4; ++slot) {
                // The loop bound used to be the slot count itself, which is every slot
                // including the empty ones -- six identical rows per channel, all named
                // PQRST_A with n_members 0, because letterRanks ranks an unlabelled
                // empty slot the same as the first real one. hasVisiblePanel decides
                // now; the bound is only a bound.
                // A ROW EXACTLY WHEN THE GRID DRAWS A PANEL. One call, shared with
                // leadsForBinTemplate, so this cannot drift from what is on screen.
                // Tested per alignment below; the row exists if ANY alignment draws it,
                // and each block blanks individually when its own average is absent.
                bool anyVisible = false;
                for (AnchorType a : anchors)
                    if (hasVisiblePanel(b, c, slot, a)) { anyVisible = true; break; }
                if (!anyVisible) continue;
                // NOT bank.templates[slot] unguarded: hasVisiblePanel also accepts the
                // pre-bank slot-0 fallback, where the bank is empty and that index does
                // not exist. Both counts are 0 there, which is the honest answer -- the
                // trace is the bin's own chN_raw, not a template with members.
                const bool haveBankTmpl = (slot < bank.size());
                const int nMembers = haveBankTmpl ? bank.templates[slot].memberCount() : 0;
                const int nClean = haveBankTmpl ? bank.templates[slot].cleanCount() : 0;
                // 0.0 is the unset value of mean_rr_ms, so it is tested rather
                // than written: the pre-bank slot-0 fallback has no template to
                // ask at all, and both cases mean "not measured".
                const double rrMs = haveBankTmpl
                    ? bank.templates[slot].mean_rr_ms : 0.0;
                // NO FALLBACK TO THE R BASE for a missing alignment. bankSlotFor is
                // explicit that substituting the R-aligned average under a non-R tag is
                // worse than omitting it, so a slot with no average for an alignment
                // gets an EMPTY trace and blanked marks for that block. The emitters
                // already blank on a negative index or a non-finite amplitude, so the
                // column count is unaffected -- which is why this is done by nulling
                // the inputs rather than by counting commas.
                static const std::vector<double> kNoTrace;
                f << fileID << ',' << b.index << ',' << kChan[c]
                    << ',' << bankSlotName(b.ecg_bank[c], slot)
                    << ',' << (b.bad_r_ch[c] ? 1 : 0)
                    // ANY NON-ZERO IS BAD. bad_ppg is 0 = ok, 1 = marked bad, 2 = no
                    // pulse present -- and no pulse is the same verdict as marked bad for
                    // a consumer of this file: there is no usable pulse either way.
                    << ',' << ((b.bad_ppg != 0) ? 1 : 0)
                    << ',' << nMembers
                    << ',' << nClean;
                f << ',';
                if (rrMs > 0.0) f << rrMs;

                // (every ECG lookup below goes through b.bankSlotFor(c, slot, anchor),
                //  which selects this row's template in this block's alignment.)

                // ONE ROW, ALL FOUR ALIGNMENTS. The alignment loop lives here
                // rather than in four calls stitched afterwards, so a row is owned
                // by one function. The ECG auto-feature groups moved INSIDE this
                // loop: the header emits them per block, and with both sections live
                // in one row the old order (all ECG points, then pulse points, then
                // ECG auto-features) no longer lines up with it.
                if (wantEcg)
                    for (AnchorType anchor : anchors) {
                        intervalUserHere = intervalUserIn(anchor);
                        {
                            // EVERYTHING ECG IN THIS BLOCK IS MEASURED ON THIS ALIGNMENT.
                            // Its own average supplies the amplitudes, its own detections
                            // supply the glyph columns, and only the one bar it owns gets a
                            // user column. Four blocks, four frames.
                            // This slot's own average, not the bin-wide one.
                            const AnchoredBankSlot* asl =
                                hasVisiblePanel(b, c, slot, anchor)
                                ? b.bankSlotFor(c, slot, anchor) : nullptr;
                            const std::vector<double>& ecg = asl ? asl->tmpl : kNoTrace;
                            const double ref = ecgRef[c];
                            // NO CACHE READ. Every field this row uses is overwritten from
                            // lm below, so autoFor's value was never consulted -- it only
                            // supplied a struct, and with it the flat-field fallback that
                            // hands back R's positions for another alignment.
                            TemplateBin::AnchorAuto aa{};
                            if (asl) {
                                const int rSeed = static_cast<int>(
                                    std::lround(b.chFor(c, anchor).r_col_raw));
                                const FeatureMarks::TemplateLandmarks lm =
                                    FeatureMarks::detect_template_landmarks(
                                        asl->tmpl, rSeed, sampleRateHz,
                                        b.polarity.sign(c), fitMode, peakMode);
                                aa.p_begin[c] = lm.p_begin;  aa.p_peak[c] = lm.p_peak;
                                aa.q_onset[c] = lm.q_onset;  aa.r_peak[c] = lm.r_peak;
                                aa.s_end[c] = lm.s_end;    aa.t_end[c] = lm.t_end;
                            }
                            else {
                                aa.p_begin[c] = aa.p_peak[c] = aa.q_onset[c] =
                                    aa.r_peak[c] = aa.s_end[c] = aa.t_end[c] = -1.0;
                            }

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
                                const tbank::BankMarkerSet& own = b.slotMarks(c, slot, anchor);
                                if (anchor_view::owns(anchor, anchor_view::p_begin)) umk.p_begin = own.p_begin;
                                if (anchor_view::owns(anchor, anchor_view::q_begin)) umk.q_onset = own.q_onset;
                                if (anchor_view::owns(anchor, anchor_view::j_point))   umk.s_end = own.s_end;
                                if (anchor_view::owns(anchor, anchor_view::t_end))   umk.t_end = own.t_end;
                            }

                            //bracket t peak by send tbegin
                            const tbank::BankMarkerSet whole = asl ? b.userMarks(c, slot, anchor)
                                : tbank::BankMarkerSet{};
                            const FeatureMarks::ReactiveEcg user_placed_s_and_t_bars_for_bracketing_tpeak = FeatureMarks::reactive_ecg(ecg, whole.p_begin, whole.q_onset, whole.s_end, whole.t_end, sampleRateHz, peakMode);
                            const FeatureMarks::ReactiveEcg auto_s_and_t_bars_for_bracketing_tpeak = FeatureMarks::reactive_ecg(ecg, aa.p_begin[c], aa.q_onset[c], aa.s_end[c], aa.t_end[c], sampleRateHz, peakMode);
                            // Empty, not computed, with no trace: computeEcgFeatures on an
                            // empty vector still returns an s_idx.
                            EcgFeatures ftAuto = asl
                                ? computeEcgFeatures(ecg, aa.p_peak[c], aa.q_onset[c], aa.r_peak[c], aa.s_end[c], aa.t_end[c], sampleRateHz, b.polarity.sign(c), peakMode)
                                : EcgFeatures{};
                            // Derived from the assembled bars, not from `umk`: q_peak,
                            // s_peak and the QRS/QT intervals need a whole beat's
                            // brackets, and no single alignment's marker set holds one any
                            // more. The two intervals come out identical in all four
                            // blocks (a duration is frame-free), which is exactly why only
                            // the R block emits their user half.
                            EcgFeatures ftUser = asl
                                ? computeEcgFeatures(ecg, user_placed_s_and_t_bars_for_bracketing_tpeak.p_peak, whole.q_onset, b.r_peak_ch[c], whole.s_end, whole.t_end, sampleRateHz, b.polarity.sign(c), peakMode)
                                : EcgFeatures{};   // see ftAuto

                            // Order MUST match ecgPointNames:
                            //   p_begin(bar), p_peak(glyph), q_onset(bar), q_peak(computed),
                            //   r_peak(glyph), s_peak(computed), s_end(bar),
                            //   t_peak(glyph), t_end(bar)
                            //
                            // The `u` field of a GLYPH is -1 and unused: hasUserColumn
                            // returns false for it, so emitEcgPoint never reads it. It is
                            // left as -1 rather than as the reactive value to make the
                            // absence explicit at the table rather than only at the emit.
                            struct P { const char* name; double a; double u; };
                            const P pts[] = {
                                { "p_begin", aa.p_begin[c],  umk.p_begin  },
                                { "p_peak",  auto_s_and_t_bars_for_bracketing_tpeak.p_peak,  -1.0         },   // glyph: auto only
                                { "q_onset", aa.q_onset[c],  umk.q_onset  },
                                { "q_peak",  ftAuto.q_idx,   -1.0         },   // glyph: auto only
                                { "r_peak",  aa.r_peak[c],   -1.0         },   // glyph: auto only
                                { "s_peak",  ftAuto.s_idx,   -1.0         },   // glyph: auto only
                                { "s_end",   aa.s_end[c],    umk.s_end    },
                                { "t_peak",  auto_s_and_t_bars_for_bracketing_tpeak.t_peak,  -1.0         },   // glyph: auto only
                                { "t_end",   aa.t_end[c],    umk.t_end    }
                            };
                            for (const P& pt : pts) {
                                // SAME FUNCTION the header emitter called, so the two
                                // cannot disagree about how many columns this point has.
                                // This was a bare `k != 3` index test, which meant adding a
                                // point column silently shifted which one lost its user
                                // variant.
                                emitEcgPoint(ecg, pt.a, pt.u, ref,
                                    anchor_view::hasUserColumn(pt.name, anchor));
                            }

                            // ---- INTERVALS: ecgIntervalNames ORDER, qrs THEN qt -------
                            // A duration is frame-free -- (s_end - q_onset) is the same
                            // number whichever alignment's columns it was measured in --
                            // so the auto side is emitted per block (the column set has to
                            // be the same shape within a part) while the user side appears
                            // once, under R. emitIntervalPair applies that gate.
                            //
                            // Both come from computeEcgFeatures, which now takes doubles,
                            // so the user side is a sub-sample duration rather than a
                            // whole number of samples times the sample period.
                            emitIntervalPair(ftAuto.qrs_ms, ftUser.qrs_ms);
                            emitIntervalPair(ftAuto.qt_ms, ftUser.qt_ms);
                        }

                        // Autodetected ECG glyph columns (p_wave, q_onset, r_wave, t_peak).
                        // p_wave/q_onset/r_wave ARE this alignment's stored detections --
                        // emitted straight, with no parallel recompute that could disagree
                        // with them. (r_wave in particular: R has exactly one definition, the
                        // alignment's own r_peak. It is never re-derived by an argmax here.)
                        // p_wave and t_peak are reactive, bracketed here by the DETECTOR's
                        // bars to match this group's "_auto" name -- t_peak by s_end/t_end,
                        // which is the bracket reactiveGlyphs uses, not the t_begin this line
                        // used to pass and nothing ever set.
                        {
                            // Same slot trace and same re-derived auto set as the point
                            // groups above, so these glyph columns cannot disagree with
                            // them.
                            const AnchoredBankSlot* asl2 =
                                hasVisiblePanel(b, c, slot, anchor)
                                ? b.bankSlotFor(c, slot, anchor) : nullptr;
                            const std::vector<double>& ecg = asl2 ? asl2->tmpl : kNoTrace;
                            TemplateBin::AnchorAuto aa{};   // see the point block above
                            if (asl2) {
                                const int rSeed = static_cast<int>(
                                    std::lround(b.chFor(c, anchor).r_col_raw));
                                const FeatureMarks::TemplateLandmarks lm =
                                    FeatureMarks::detect_template_landmarks(
                                        asl2->tmpl, rSeed, sampleRateHz,
                                        b.polarity.sign(c), fitMode, peakMode);
                                aa.p_begin[c] = lm.p_begin;  aa.p_peak[c] = lm.p_peak;
                                aa.q_onset[c] = lm.q_onset;  aa.r_peak[c] = lm.r_peak;
                                aa.s_end[c] = lm.s_end;    aa.t_end[c] = lm.t_end;
                            }
                            else {
                                aa.p_begin[c] = aa.p_peak[c] = aa.q_onset[c] =
                                    aa.r_peak[c] = aa.s_end[c] = aa.t_end[c] = -1.0;
                            }
                            const FeatureMarks::ReactiveEcg rx = FeatureMarks::reactive_ecg(
                                ecg, aa.p_begin[c], aa.q_onset[c],
                                aa.s_end[c], aa.t_end[c], sampleRateHz, peakMode);
                            emitAutoFeatPt(ecg, rx.p_peak);
                            emitAutoFeatPt(ecg, aa.q_onset[c]);
                            emitAutoFeatPt(ecg, aa.r_peak[c]);
                            emitAutoFeatPt(ecg, rx.t_peak);
                        }
                    }
                // (bad_ecg and bad_ppg are row keys -- emitted once, above.)
                if (wantPulse) {
                    // PPG: onset, t50, peak, dicrotic, peak2, t80, t80_rise, end
                    //onset, dicrotic, and end are the only user movable bars
                    // THIS SLOT'S PULSE AND THIS SLOT'S BARS. The rows here
                    // are already per (lead, slot); only the pulse columns were
                    // still bin-level, so every column of a bin reported one
                    // set of pulse marks measured against a waveform that was
                    // only correct for one of them. pulse_marks is where the
                    // bars live now (and what the .bin serializes per slot).
                    //
                    // The *_auto columns come from the same place: they are the
                    // detector's own positions on THIS slot's pulse, where
                    // b.ppg_*_auto describe b.ppgTemplate -- a different
                    // waveform for every column but one.
                    // BY VALUE, not by reference: the ternary below has a
                    // temporary in one branch, so it yields a prvalue and a
                    // const& would bind to a copy in both cases anyway. Saying
                    // so beats relying on lifetime extension to read right.
                    const tbank::BankPulseMarkerSet pmU =
                        (slot >= 0 && slot < static_cast<int>(b.ppg_bank.templates.size()))
                        ? b.ppg_bank.templates[slot].pulse_marks
                        : tbank::BankPulseMarkerSet{};
                    const std::vector<double>& pulseU =
                        (slot >= 0 && slot < static_cast<int>(b.ppg_bank.templates.size()))
                        ? b.ppg_bank.templates[slot].tmpl
                        : b.ppgTemplate;
                    const FeatureMarks::ReactivePpg rxUser =
                        FeatureMarks::reactive_ppg(pulseU, pmU.onset,
                            pmU.peak_auto, pmU.dicrotic, pmU.end);
                    const FeatureMarks::ReactivePpg rxAuto =
                        FeatureMarks::reactive_ppg(pulseU, pmU.onset_auto,
                            pmU.peak_auto, pmU.dicrotic_auto, pmU.end_auto);
                    // EVERY ARGUMENT FROM THIS SLOT: its pulse, its detector
                    // columns, its bars, and its own foot as the perfusion
                    // baseline. ppg_peak is auto-only (markerAtX hands it out
                    // to nobody), so its user half is the detector's column --
                    // not a placement, and the same answer the bin fields gave.
                    emitPulsePoint("ppg_onset", pulseU, pmU.onset_auto, pmU.onset, pmU.onset_auto, pmU.onset, refPpg);
                    emitPulsePoint("ppg_t50", pulseU, rxAuto.t50, rxUser.t50, pmU.onset_auto, pmU.onset, refPpg);
                    emitPulsePoint("ppg_peak", pulseU, pmU.peak_auto, pmU.peak_auto, pmU.onset_auto, pmU.onset, refPpg);
                    emitPulsePoint("ppg_dicr", pulseU, pmU.dicrotic_auto, pmU.dicrotic, pmU.onset_auto, pmU.onset, refPpg);
                    emitPulsePoint("ppg_peak2", pulseU, rxAuto.peak2, rxUser.peak2, pmU.onset_auto, pmU.onset, refPpg);
                    emitPulsePoint("ppg_t80", pulseU, rxAuto.t80, rxUser.t80, pmU.onset_auto, pmU.onset, refPpg);
                    // T80_rise: upslope point at t80's level (a position, like t80).
                    emitPulsePoint("ppg_t80_rise", pulseU, rxAuto.t80_rise, rxUser.t80_rise, pmU.onset_auto, pmU.onset, refPpg);
                    emitPulsePoint("ppg_end", pulseU, pmU.end_auto, pmU.end, pmU.onset_auto, pmU.onset, refPpg);
                    // PW80 width, ms only, autodetect bracketing (t80 - t80_rise).
                    // Single value; blank when unavailable. Matches the one
                    // ppg_pw80_ms_auto header column.
                    f << ',';  if (rxAuto.pw80 >= 0.0) f << (rxAuto.pw80 * toMs);

                    // ARTERIAL: all five are bars on every one of these channels, so
                    // pulseHasUserColumn returns true for all fifteen names and each
                    // group emits its full six columns.
                    f << ',' << static_cast<int>(b.abp_issue);
                    emitPulsePoint("abp_onset", b.abpTemplate, b.abp_onset_auto, b.abp_onset,
                        b.abp_onset_auto, b.abp_onset, refAbp);
                    emitPulsePoint("abp_peak", b.abpTemplate, b.abp_peak_auto, b.abp_peak,
                        b.abp_onset_auto, b.abp_onset, refAbp);
                    emitPulsePoint("abp_dicr", b.abpTemplate, b.abp_dicrotic_auto, b.abp_dicrotic,
                        b.abp_onset_auto, b.abp_onset, refAbp);
                    emitPulsePoint("abp_peak2", b.abpTemplate, b.abp_peak2_auto, b.abp_peak2,
                        b.abp_onset_auto, b.abp_onset, refAbp);
                    emitPulsePoint("abp_end", b.abpTemplate, b.abp_end_auto, b.abp_end,
                        b.abp_onset_auto, b.abp_onset, refAbp);

                    f << ',' << static_cast<int>(b.art_issue);
                    emitPulsePoint("art_onset", b.artTemplate, b.art_onset_auto, b.art_onset,
                        b.art_onset_auto, b.art_onset, refArt);
                    emitPulsePoint("art_peak", b.artTemplate, b.art_peak_auto, b.art_peak,
                        b.art_onset_auto, b.art_onset, refArt);
                    emitPulsePoint("art_dicr", b.artTemplate, b.art_dicrotic_auto, b.art_dicrotic,
                        b.art_onset_auto, b.art_onset, refArt);
                    emitPulsePoint("art_peak2", b.artTemplate, b.art_peak2_auto, b.art_peak2,
                        b.art_onset_auto, b.art_onset, refArt);
                    emitPulsePoint("art_end", b.artTemplate, b.art_end_auto, b.art_end,
                        b.art_onset_auto, b.art_onset, refArt);

                    f << ',' << static_cast<int>(b.art_pulm_issue);
                    emitPulsePoint("art_pulm_onset", b.artPulmTemplate, b.art_pulm_onset_auto, b.art_pulm_onset,
                        b.art_pulm_onset_auto, b.art_pulm_onset, refArtPulm);
                    emitPulsePoint("art_pulm_peak", b.artPulmTemplate, b.art_pulm_peak_auto, b.art_pulm_peak,
                        b.art_pulm_onset_auto, b.art_pulm_onset, refArtPulm);
                    emitPulsePoint("art_pulm_dicr", b.artPulmTemplate, b.art_pulm_dicrotic_auto, b.art_pulm_dicrotic,
                        b.art_pulm_onset_auto, b.art_pulm_onset, refArtPulm);
                    emitPulsePoint("art_pulm_peak2", b.artPulmTemplate, b.art_pulm_peak2_auto, b.art_pulm_peak2,
                        b.art_pulm_onset_auto, b.art_pulm_onset, refArtPulm);
                    emitPulsePoint("art_pulm_end", b.artPulmTemplate, b.art_pulm_end_auto, b.art_pulm_end,
                        b.art_pulm_onset_auto, b.art_pulm_onset, refArtPulm);
                } // end if (wantPulse) pulse point groups
                if (wantPulse) {
                    for (const auto& gl : ppg_and_artpulse_automated_markers) {
                        emitAutoFeatPt(b.ppgTemplate, b.*gl.idx);
                        if (gl.found) f << ',' << (b.*gl.found ? 1 : 0);
                    }
                }

                f << '\n';
            }   // slot
        }     // channel
    }       // bin
}

// Path-taking wrapper, for callers outside the viewer's one-save merge. The
// ostream overload above is what the viewer uses, so a save can build each
// alignment's part in memory instead of staging it through a temp file.
inline void writeTemplateMarkingsCsv(const std::string& path,
    const std::vector<TemplateBin>& bins,
    const std::string& fileID,
    double sampleRateHz,
    AnchorType anchor,
    MarkingsCsvSection section,
    curve_fit::PeakFitMode peakMode = curve_fit::PeakFitMode::Auto,
    curve_fit::FitMode fitMode = curve_fit::FitMode::Auto)
{
    std::ofstream f(path);
    if (!f.is_open())
        throw std::runtime_error("cannot open for write: " + path);
    writeTemplateMarkingsCsv(f, bins, fileID, sampleRateHz, anchor, section,
        peakMode, fitMode);
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

    // Header: magic, version, bin count. The magic is what makes a
    // pre-version file an ERROR rather than a silent misparse -- without it the
    // old first field (a bin count) was read as a bin count and everything
    // after it shifted by eight bytes with no complaint.
    uint32_t magic = 0, ver = 0;
    f.read(reinterpret_cast<char*>(&magic), 4);
    f.read(reinterpret_cast<char*>(&ver), 4);
    if (magic != kMarkMagic)
        throw std::runtime_error(
            "not a template-markings file (bad magic), or written before the "
            "version field existed -- re-mark: " + path);
    if (ver != kMarkVersion)
        throw std::runtime_error(
            "template-markings version " + std::to_string(ver)
            + " is not supported (this build writes and reads version "
            + std::to_string(kMarkVersion) + "): " + path);

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
        for (int lead = 0; lead < 3; ++lead) {
            const int nSlots = r32();
            for (int slot = 0; slot < nSlots; ++slot) {
                const int nAnchors = r32();
                for (int a = 0; a < nAnchors; ++a) {
                    const int tag = r32();
                    tbank::BankMarkerSet& m =
                        b.slotMarks(lead, slot, static_cast<AnchorType>(tag));
                    m.p_begin = r64d(); m.q_onset = r64d();
                    m.s_end = r64d();   m.t_end = r64d();
                }
            }
        }
        // PPG bars, per slot. slotMarks' pulse counterpart: the bank may be
        // shorter than the record (a bin read from the markings file has no
        // templates yet), so the slots are created to hold the marks and the
        // real bank merges onto them later -- the same thing the ECG block
        // above relies on slotMarks for.
        {
            const int nPulse = r32();
            if (nPulse > 0) {
                if (static_cast<int>(b.ppg_bank.templates.size()) < nPulse)
                    b.ppg_bank.templates.resize(static_cast<size_t>(nPulse));
                for (int slot = 0; slot < nPulse; ++slot) {
                    tbank::BankPulseMarkerSet& pm =
                        b.ppg_bank.templates[slot].pulse_marks;
                    pm.onset = r64d();
                    pm.dicrotic = r64d();
                    pm.end = r64d();
                }
            }
        }

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