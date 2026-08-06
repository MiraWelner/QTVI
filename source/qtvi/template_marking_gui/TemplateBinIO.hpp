#pragma once
//
// Adapter layer: reads the template_io binary format and projects it into
// the TemplateBin shape the viewer expects. Also writes/reads the
// separate template_markings.bin produced by the viewer.
//
// Marker set per bin:
//   ECG (per channel):  P-onset, Q-begin, T-begin, T-end
//   PPG (shared):       Onset, Peak, Dicrotic notch, 50% recovery, End
//
// All marker sample indices use -1 as the "unmarked / not applicable"
// sentinel.
//
// IQR vectors (ecgTemplate_raw_iqr per channel, ppg_template_iqr) come
// straight from the template file. Empty std => the widget renders the
// trace without a gray band.
//

#include <vector>
#include <string>
#include <map>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <utility>
#include <fstream>
#include <stdexcept>
#include <algorithm> 
#include <cstring>

#include "template_generation\template_io.hpp"
#include "template_marking_gui\feature_marks.hpp"


// ---------------------------------------------------------------------------
// In-memory model used by the viewer.
// ---------------------------------------------------------------------------

struct ChannelTemplateData {
    std::vector<double> ecgTemplate_raw;
    std::vector<double> ecg_template_raw_iqr;
    std::vector<double> ecgTemplate_squared;   // unused by viewer
    std::vector<double> ecgTemplate_absval;    // unused by viewer
    double alignment_point_raw = 0;
    double alignment_point_squared = 0;
    double alignment_point_absval = 0;
    // True R column in the template (from alignment). Used directly as the R
    // fiducial; replaces the old avg_r_expand positioning constant.
    int r_col_raw = -1;
    int r_col_squared = -1;
    int r_col_absval = -1;
};

struct TemplateBin {
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
    struct MarkerSet {
        int p_begin_ch[3] = { -1, -1, -1 };
        int p_peak_ch[3] = { -1, -1, -1 };
        int q_begin_ch[3] = { -1, -1, -1 };
        int s_end_ch[3] = { -1, -1, -1 };
        int t_begin_ch[3] = { -1, -1, -1 };
        int t_end_ch[3] = { -1, -1, -1 };
    };
    std::map<AnchorType, MarkerSet> markers_by_anchor;

    // Accessor: returns the marker set for `a`, creating an all-unmarked one
    // on first access. Every read/write of a per-anchor ECG user marker goes
    // through here.
    MarkerSet& marks(AnchorType a) { return markers_by_anchor[a]; }
    const MarkerSet& marks(AnchorType a) const {
        static const MarkerSet kEmpty;
        auto it = markers_by_anchor.find(a);
        return (it == markers_by_anchor.end()) ? kEmpty : it->second;
    }

    // R peak column: auto-only, re-derived each pass from the template r_col;
    // NOT per-anchor and NOT persisted. Kept flat.
    int r_peak_ch[3] = { -1, -1, -1 };
    //sample indicies for each of the 3 ECG channels, auto-detected. -1 = unmarked / not applicable.
    // Auto-detected ECG landmark positions, stored as DOUBLE (sub-sample
    // precision from the fit/refine stages; -1 = unset). NOT serialized --
    // recomputed every loadSubject -- so widening to double is not a .bin
    // format change. Consumers that need an integer sample index round at use.
    double p_peak_auto_ch[3] = { -1, -1, -1 };
    double q_begin_auto_ch[3] = { -1, -1, -1 };
    double r_peak_auto_ch[3] = { -1, -1, -1 };
    double s_end_auto_ch[3] = { -1, -1, -1 };
    double t_begin_auto_ch[3] = { -1, -1, -1 };
    double t_end_auto_ch[3] = { -1, -1, -1 };
    double p_begin_auto_ch[3] = { -1, -1, -1 };

    // PPG: sample indices into ppgTemplate. Shared across channels.
    int ppg_onset = -1;
    int ppg_t50 = -1;   // 50% up the upslope, foot -> systolic peak
    int ppg_t80 = -1; // 80% up the upslope, foot -> systolic peak
    int ppg_peak = -1;
    int ppg_dicrotic = -1;
    int ppg_peak2 = -1;
    int ppg_end = -1;

    // Construction-time PPG fiducials (peak = max in [R1,R2], foot = min in
    // [R1,peak]), computed from the real R-pair interval at template build and
    // carried through the template file. seed_all uses these directly instead
    // of re-detecting the peak on the multi-pulse template. -1 if unavailable.
    int ppg_peak_construct = -1;
    int ppg_onset_construct = -1;

    // Arterial channels (ABP / ART / ART_PULM): same marker set as PPG
    // (onset, peak, dicrotic, 50%, end) and same issue flag semantics
    // (0 = ok, 1 = bad, 2 = channel absent). Indices are into the matching
    // *Template vector above; -1 = unmarked / not applicable.
    uint8_t abp_issue = 0;
    int abp_onset = -1, abp_peak = -1, abp_dicrotic = -1, abp_peak2 = -1, abp_end = -1;

    uint8_t art_issue = 0;
    int art_onset = -1, art_peak = -1, art_dicrotic = -1, art_peak2 = -1, art_end = -1;

    uint8_t art_pulm_issue = 0;
    int art_pulm_onset = -1, art_pulm_peak = -1, art_pulm_dicrotic = -1,
        art_pulm_peak2 = -1, art_pulm_end = -1;

    // Auto-detect mirror fields for pulse channels (same rule as ECG auto
    // mirrors above: filled every loadSubject, NOT serialized). Preserve
    // the original auto positions so the CSV can emit both variants.
    int ppg_onset_auto = -1, ppg_t50_auto = -1, ppg_peak_auto = -1,
        ppg_dicrotic_auto = -1, ppg_peak2_auto = -1, ppg_end_auto = -1;
    int ppg_t80_auto = -1;
    // Whether the dicrotic notch / diastolic peak / pulse end were
    // genuinely found by the detector vs. fell back to a placeholder
    // position. Drive X-vs-O glyph rendering; set once by seed_all's
    // single PPG detection pass. Not written to CSV -- the CSV only
    // needs the x/y values for each landmark, since a fallback landmark
    // still has valid coordinates.
    bool ppg_dicrotic_found_auto = false;
    bool ppg_peak2_found_auto = false;
    bool ppg_end_found_auto = false;
    int abp_onset_auto = -1, abp_peak_auto = -1, abp_dicrotic_auto = -1,
        abp_peak2_auto = -1, abp_end_auto = -1;
    int art_onset_auto = -1, art_peak_auto = -1, art_dicrotic_auto = -1,
        art_peak2_auto = -1, art_end_auto = -1;
    int art_pulm_onset_auto = -1, art_pulm_peak_auto = -1, art_pulm_dicrotic_auto = -1,
        art_pulm_peak2_auto = -1, art_pulm_end_auto = -1;

};

// ---------------------------------------------------------------------------
// Read: convert template_io::TemplateFile -> std::vector<TemplateBin>
// ---------------------------------------------------------------------------
inline std::vector<TemplateBin> readTemplateInfoBin(const std::string& path,
    AnchorType anchor = AnchorType::R_PEAK) {
    template_io::TemplateFile tf = template_io::read_template_binfile(path);

    // Non-R anchors live in tf.raw_anchors (per-bin 3-channel raw). R_PEAK
    // (and any anchor not present in the file) uses the scalar base in
    // bins[i].chN_raw. Look the requested anchor up once.
    const int anchorTag = static_cast<int>(anchor);
    const auto anchorIt = tf.raw_anchors.find(anchorTag);
    const bool useAnchor = (anchor != AnchorType::R_PEAK)
        && (anchorIt != tf.raw_anchors.end());

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
        dst.ppg_peak_construct = src.ppg_peak_col;
        dst.ppg_onset_construct = src.ppg_onset_col;
        dst.ppgTemplate = src.ppgTemplate;
        dst.ppg_template_iqr = src.ppg_template_iqr;
        dst.abpTemplate = src.abpTemplate;
        dst.artTemplate = src.artTemplate;
        dst.artPulmTemplate = src.artPulmTemplate;
        dst.abpTemplate_iqr = src.abpTemplate_iqr;
        dst.artTemplate_iqr = src.artTemplate_iqr;
        dst.artPulmTemplate_iqr = src.artPulmTemplate_iqr;

        // Pick the raw ECG source: the requested anchor's block if present
        // and non-R, else the scalar base (R_PEAK). Guard the per-bin index
        // in case an anchor block is shorter than the base (shouldn't happen,
        // but never index past it).
        const bool anchorHasBin = useAnchor && (i < anchorIt->second.size());
        const template_io::ChannelMethodTemplate& c1 =
            anchorHasBin ? anchorIt->second[i][0] : src.ch1_raw;
        const template_io::ChannelMethodTemplate& c2 =
            anchorHasBin ? anchorIt->second[i][1] : src.ch2_raw;
        const template_io::ChannelMethodTemplate& c3 =
            anchorHasBin ? anchorIt->second[i][2] : src.ch3_raw;

        dst.ch1.ecgTemplate_raw = c1.ecgTemplate;
        dst.ch1.ecg_template_raw_iqr = c1.ecg_template_iqr;
        dst.ch1.alignment_point_raw = c1.alignment_point;
        dst.ch1.r_col_raw = c1.r_col;

        dst.ch2.ecgTemplate_raw = c2.ecgTemplate;
        dst.ch2.ecg_template_raw_iqr = c2.ecg_template_iqr;
        dst.ch2.alignment_point_raw = c2.alignment_point;
        dst.ch2.r_col_raw = c2.r_col;

        dst.ch3.ecgTemplate_raw = c3.ecgTemplate;
        dst.ch3.ecg_template_raw_iqr = c3.ecg_template_iqr;
        dst.ch3.alignment_point_raw = c3.alignment_point;
        dst.ch3.r_col_raw = c3.r_col;
    }
    return bins;
}

// ---------------------------------------------------------------------------
// template_markings.bin layout:
//
//   header:
//     uint64  numBins
//
//   per bin:
//     uint64  index
//     uint8   bad_r_ch1, bad_r_ch2, bad_r_ch3
//     uint8   ppg_issue          (0 = ok, 1 = bad, 2 = no ppg)
//     int32   p_peak_ch1..3, q_begin_ch1..3, r_peak_ch1..3,
//             s_end_ch1..3, t_begin_ch1..3, t_end_ch1..3
//     int32   ppg_onset, ppg_p50, ppg_peak, ppg_dicrotic, ppg_peak2, ppg_end, ppg_t80
//     -- arterial, one block each for ABP, ART, ART_PULM: --
//     uint8   <chan>_issue
//     int32   <chan>_onset, _peak, _dicrotic, _peak2, _end
//
// All int32 fields use -1 as the "unmarked / not applicable" sentinel.
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

    for (const auto& b : bins) {
        uint64_t idx = b.index;
        f.write(reinterpret_cast<const char*>(&idx), 8);

        w8(b.bad_r_ch[0] ? 1 : 0);
        w8(b.bad_r_ch[1] ? 1 : 0);
        w8(b.bad_r_ch[2] ? 1 : 0);
        w8(b.bad_ppg);

        // Per-anchor ECG user markers: [uint32 count] then, per anchor,
        // [int32 anchorTag][6 marker fields x 3 channels]. R peak is auto-only
        // and re-derived each pass, so it is NOT written here.
        w32(static_cast<int>(b.markers_by_anchor.size()));
        for (const auto& kv : b.markers_by_anchor) {
            w32(static_cast<int>(kv.first));           // AnchorType tag
            const TemplateBin::MarkerSet& m = kv.second;
            for (int c = 0; c < 3; ++c) w32(m.p_begin_ch[c]);
            for (int c = 0; c < 3; ++c) w32(m.p_peak_ch[c]);
            for (int c = 0; c < 3; ++c) w32(m.q_begin_ch[c]);
            for (int c = 0; c < 3; ++c) w32(m.s_end_ch[c]);
            for (int c = 0; c < 3; ++c) w32(m.t_begin_ch[c]);
            for (int c = 0; c < 3; ++c) w32(m.t_end_ch[c]);
        }

        w32(b.ppg_onset);
        w32(b.ppg_t50);
        w32(b.ppg_peak);
        w32(b.ppg_dicrotic);
        w32(b.ppg_peak2);
        w32(b.ppg_end);
        w32(b.ppg_t80);

        // Arterial block: ABP, ART, ART_PULM (issue + 5 indices each).
        w8(b.abp_issue);
        w32(b.abp_onset); w32(b.abp_peak); w32(b.abp_dicrotic); w32(b.abp_peak2); w32(b.abp_end);
        w8(b.art_issue);
        w32(b.art_onset); w32(b.art_peak); w32(b.art_dicrotic); w32(b.art_peak2); w32(b.art_end);
        w8(b.art_pulm_issue);
        w32(b.art_pulm_onset); w32(b.art_pulm_peak); w32(b.art_pulm_dicrotic);
        w32(b.art_pulm_peak2); w32(b.art_pulm_end);
    }
}


struct EcgFeatures {
    int q_idx = -1, r_idx = -1, s_idx = -1;   // peak sample positions
    double qrs_ms = NAN, qt_ms = NAN;
};

inline EcgFeatures computeEcgFeatures(const std::vector<double>& ecg, int p_peak, int q_begin, int r_peak, int s_end, int t_end, double rateHz)
{
    EcgFeatures f;
    const int N = static_cast<int>(ecg.size());
    const double msPerSamp = (rateHz > 0.0) ? 1000.0 / rateHz : NAN;
    auto inRange = [&](int i) { return i >= 0 && i < N; };

    if (q_begin >= 0 && s_end >= q_begin) f.qrs_ms = (s_end - q_begin) * msPerSamp;
    if (q_begin >= 0 && t_end >= q_begin) f.qt_ms = (t_end - q_begin) * msPerSamp;

    if (inRange(r_peak)) f.r_idx = r_peak;
    f.q_idx = FeatureMarks::compute_q_peak(ecg, q_begin, r_peak);
    // S for |R|+|S| = first opposite-polarity trough after R (robust; not the
    // max over [R, s_end], which depends on where s_end sits).
    f.s_idx = FeatureMarks::compute_s_peak(ecg, r_peak, rateHz);
    return f;
}

// Per-channel ECG column order. isUser => straight off a marker (gets _user).
// isInterval => qrs/qt, ms only (no y, no x/user).
// Under v4, the six ECG markers are ALL user-placed:
//   p_peak, q_begin, r_peak, s_end, t_peak, t_end.
// Q peak and S peak inside the QRS complex are still computed
// (argmin/max between q_begin and s_end); QRS and QT intervals derive
// from q_begin and t_end as before.
struct EcgColSpec { const char* name; bool isUser; bool isInterval; };
static const EcgColSpec ecgCols[] = {
    {"p_peak",  true,  false}, {"q_begin", true,  false},
    {"q_peak",  false, false}, {"r_peak",  true,  false}, {"s_peak", false, false},
    {"s_end",   true,  false}, {"t_begin", true,  false}, {"t_end",  true,  false},
    {"qrs",     false, true},  {"qt",      false, true},
};
inline constexpr const char* ppgCols[] = { "ppg_onset","ppg_p50","ppg_peak","ppg_dicrotic","ppg_peak2","ppg_t80","ppg_end" };
inline constexpr const char* abpCols[] = { "abp_onset","abp_peak","abp_dicrotic","abp_peak2","abp_end" };
inline constexpr const char* artCols[] = { "art_onset","art_peak","art_dicrotic","art_peak2","art_end" };
inline constexpr const char* artPulmCols[] = { "art_pulm_onset","art_pulm_peak","art_pulm_dicrotic","art_pulm_peak2","art_pulm_end" };

// ---------------------------------------------------------------------------
// writeTemplateMarkingsCsv
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
// Intervals (qrs_ms, qt_ms) emit as pairs: {name}_ms_autodetect, {name}_ms_user.
// Q peak and S peak are computed inside the QRS -- their autodetect variant
// uses (Q_begin_auto, R_peak_auto, S_end_auto); user variant uses the current
// user markers.
// ---------------------------------------------------------------------------
enum class MarkingsCsvSection { EcgOnly, PulseOnly };

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
                const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
                const auto& ecg = chs[c]->ecgTemplate_raw;
                if (ecg.empty()) continue;
                EcgFeatures ft = computeEcgFeatures(ecg,
                    (int)std::lround(b.p_peak_auto_ch[c]), (int)std::lround(b.q_begin_auto_ch[c]), (int)std::lround(b.r_peak_auto_ch[c]),
                    (int)std::lround(b.s_end_auto_ch[c]), (int)std::lround(b.t_end_auto_ch[c]),
                    sampleRateHz);
                if (ft.r_idx < 0 || ft.s_idx < 0) continue;
                if (ft.r_idx >= (int)ecg.size() || ft.s_idx >= (int)ecg.size()) continue;
                const double ry = ecg[ft.r_idx], sy = ecg[ft.s_idx];
                if (std::isnan(ry) || std::isnan(sy)) continue;
                vals.push_back(std::abs(ry) + std::abs(sy));
            }
            ecgRef[c] = medianFinite(std::move(vals));
        }
    }
    double refPpg = std::nan(""), refAbp = std::nan(""), refArt = std::nan(""), refArtPulm = std::nan("");
    if (wantPulse) {
        auto pulseRefAuto = [&](const std::vector<double> TemplateBin::* trace,
            const int TemplateBin::* footAuto,
            const int TemplateBin::* peakAuto,
            bool checkPpgIssue) -> double
            {
                std::vector<double> vals;
                for (const auto& b : bins) {
                    if (b.bad_segment) continue;
                    if (checkPpgIssue && b.bad_ppg != 0) continue;
                    const auto& v = b.*trace;
                    const int fi = b.*footAuto;
                    const int pi = b.*peakAuto;
                    if (fi < 0 || pi < 0 || fi >= (int)v.size() || pi >= (int)v.size()) continue;
                    const double fy = v[fi], py = v[pi];
                    if (std::isnan(fy) || std::isnan(py) || std::abs(fy) < 1e-12) continue;
                    vals.push_back(100.0 * (py - fy) / fy);
                }
                return medianFinite(std::move(vals));
            };
        refPpg = pulseRefAuto(&TemplateBin::ppgTemplate,
            &TemplateBin::ppg_onset_auto, &TemplateBin::ppg_peak_auto, true);
        refAbp = pulseRefAuto(&TemplateBin::abpTemplate,
            &TemplateBin::abp_onset_auto, &TemplateBin::abp_peak_auto, false);
        refArt = pulseRefAuto(&TemplateBin::artTemplate,
            &TemplateBin::art_onset_auto, &TemplateBin::art_peak_auto, false);
        refArtPulm = pulseRefAuto(&TemplateBin::artPulmTemplate,
            &TemplateBin::art_pulm_onset_auto, &TemplateBin::art_pulm_peak_auto, false);
    }

    // ---- header ------------------------------------------------------------
    // Keys. ppg_issue is PULSE metadata -> only in the pulse section.
    f << "file_id,bin_index,bad_r_ch1,bad_r_ch2,bad_r_ch3";
    if (wantPulse) f << ",ppg_issue";

    // ECG 8 point columns + 2 interval columns per channel.
    static const char* ecgPointNames[] = {
        "p_peak", "q_begin", "q_peak", "r_peak", "s_peak",
        "s_end",  "t_peak", "t_begin",  "t_end"
    };
    static const char* ecgIntervalNames[] = { "qrs", "qt" };

    auto emitEcgPointHeader = [&](const char* name, int c) {
        // r_peak is the only autodetect-only column (R is never user-placed).
        // t_peak DOES get a user column: it's reactive, so it has a distinct
        // value under the auto brackets and under the user brackets, and the
        // user one is what the operator sees on screen. Keep this rule and the
        // row loop's `k != 3` test in step -- they are the same rule.
        const bool userToo = (std::strcmp(name, "r_peak") != 0);
        f << ',' << name << "_ch" << c << "_y_normalized_autodetect";
        if (userToo) f << ',' << name << "_ch" << c << "_y_normalized_user";
        f << ',' << name << "_ch" << c << "_y_mv_raw_autodetect";
        if (userToo) f << ',' << name << "_ch" << c << "_y_mv_raw_user";
        f << ',' << name << "_ch" << c << "_x_ms_autodetect";
        if (userToo) f << ',' << name << "_ch" << c << "_x_ms_user";
        };
    auto emitIntervalHeader = [&](const char* name, int c) {
        f << ',' << name << "_ch" << c << "_ms_autodetect"
            << ',' << name << "_ch" << c << "_ms_user";
        };
    // Pulse: 6 cols per marker.
    auto emitPulsePointHeader = [&](const char* name) {
        f << ',' << name << "_y_normalized_autodetect"
            << ',' << name << "_y_normalized_user"
            << ',' << name << "_y_mv_raw_autodetect"
            << ',' << name << "_y_mv_raw_user"
            << ',' << name << "_x_ms_autodetect"
            << ',' << name << "_x_ms_user";
        };
    // Autodetected computed features (no user bar; derived from AUTODETECT markers).
    auto emitAutoFeatHeader = [&](const char* name) {
        f << ',' << name << "_x_ms" << ',' << name << "_y_mv_raw";
        };

    if (wantEcg) {
        for (int c = 1; c <= 3; ++c) {
            for (const char* n : ecgPointNames)     emitEcgPointHeader(n, c);
        }
        for (int c = 1; c <= 3; ++c) {
            char b[64];
            for (const char* g : { "p_wave_autodetect", "q_onset_autodetect",
                                   "r_wave_autodetect", "t_peak_autodetect" }) {
                std::snprintf(b, sizeof b, "%s_ch%d", g, c);
                emitAutoFeatHeader(b);
            }
        }
    }

    if (wantPulse) {
        for (const char* n : ppgCols)     emitPulsePointHeader(n);
        f << ",abp_issue";
        for (const char* n : abpCols)     emitPulsePointHeader(n);
        f << ",art_issue";
        for (const char* n : artCols)     emitPulsePointHeader(n);
        f << ",art_pulm_issue";
        for (const char* n : artPulmCols) emitPulsePointHeader(n);
        for (const char* g : { "ppg_foot_autodetect", "ppg_systolic_peak_autodetect",
                               "ppg_dicrotic_autodetect", "ppg_end_autodetect" })
            emitAutoFeatHeader(g);
        f << ",ppg_notch_found";
        f << ",ppg_end_found";
    }
    f << '\n';

    // ---- row loop ----------------------------------------------------------
    const double toMs = (sampleRateHz > 0.0) ? 1000.0 / sampleRateHz : 1.0;

    // Emit one 6-column ECG point group: normalized (auto/user), raw
    // (auto/user), x_ms (auto/user). Any missing piece leaves that field blank.
    auto emitEcgPoint = [&](const std::vector<double>& ecg,
        int idx_auto, int idx_user, double ref, bool userToo)
        {
            auto y_of = [&](int idx) -> double {
                if (idx < 0 || idx >= (int)ecg.size()) return std::nan("");
                const double y = ecg[idx];
                return std::isnan(y) ? std::nan("") : y;
                };
            const double y_a = y_of(idx_auto);
            const double y_u = y_of(idx_user);
            const bool refOk = std::isfinite(ref) && ref != 0.0;

            f << ',';   if (std::isfinite(y_a) && refOk) f << (y_a / ref);
            if (userToo) { f << ','; if (std::isfinite(y_u) && refOk) f << (y_u / ref); }
            f << ',';   if (std::isfinite(y_a)) f << y_a;
            if (userToo) { f << ','; if (std::isfinite(y_u)) f << y_u; }
            f << ',';   if (idx_auto >= 0) f << (idx_auto * toMs);
            if (userToo) { f << ','; if (idx_user >= 0) f << (idx_user * toMs); }
        };

    // Emit one 6-column pulse point group. footIdx_auto / footIdx_user are
    // the "onset" indices for their respective variants (used to compute the
    // local ratio (y - foot)/foot). ref is the channel's global PI median.
    auto emitPulsePoint = [&](const std::vector<double>& v,
        int idx_auto, int idx_user,
        int foot_auto, int foot_user, double ref)
        {
            auto y_of = [&](int idx) -> double {
                if (idx < 0 || idx >= (int)v.size()) return std::nan("");
                const double y = v[idx];
                return std::isnan(y) ? std::nan("") : y;
                };
            auto normOf = [&](double y, int fIdx) -> double {
                if (!std::isfinite(y)) return std::nan("");
                if (fIdx < 0 || fIdx >= (int)v.size()) return std::nan("");
                const double fy = v[fIdx];
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
            f << ',';   if (idx_auto >= 0)      f << (idx_auto * toMs);
            f << ',';   if (idx_user >= 0)      f << (idx_user * toMs);
        };

    auto emitIntervalPair = [&](double auto_ms, double user_ms) {
        f << ',';   if (std::isfinite(auto_ms)) f << auto_ms;
        f << ',';   if (std::isfinite(user_ms)) f << user_ms;
        };

    // Autodetected computed feature point (used by both ECG and PPG glyph
    // blocks -- defined once here so it's in scope for both).
    auto emitAutoFeatPt = [&](const std::vector<double>& sig, int idx) {
        f << ',';   if (idx >= 0) f << (idx * toMs);
        f << ',';   if (idx >= 0 && idx < (int)sig.size() && !std::isnan(sig[idx])) f << sig[idx];
        };

    for (const auto& b : bins) {
        f << fileID << ',' << b.index << ','
            << (b.bad_r_ch[0] ? 1 : 0) << ','
            << (b.bad_r_ch[1] ? 1 : 0) << ','
            << (b.bad_r_ch[2] ? 1 : 0);
        if (wantPulse) f << ',' << static_cast<int>(b.bad_ppg);

        const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };

        if (wantEcg)
            for (int c = 0; c < 3; ++c) {
                const auto& ecg = chs[c]->ecgTemplate_raw;
                const double ref = ecgRef[c];

                EcgFeatures ftAuto = computeEcgFeatures(ecg,
                    (int)std::lround(b.p_peak_auto_ch[c]), (int)std::lround(b.q_begin_auto_ch[c]), (int)std::lround(b.r_peak_auto_ch[c]),
                    (int)std::lround(b.s_end_auto_ch[c]), (int)std::lround(b.t_end_auto_ch[c]),
                    sampleRateHz);
                const TemplateBin::MarkerSet& umk = b.marks(anchor);
                EcgFeatures ftUser = computeEcgFeatures(ecg,
                    umk.p_peak_ch[c], umk.q_begin_ch[c], b.r_peak_ch[c],
                    umk.s_end_ch[c], umk.t_end_ch[c],
                    sampleRateHz);

                // T-peak is reactive: there's no T-peak bar, it's the extremum
                // between the T-begin and T-end bars. So it has two honest
                // values, and both are emitted through the one shared
                // FeatureMarks::reactive_ecg -- the autodetect column bracketed
                // by the *_auto_ch columns, the user column bracketed by the
                // operator's own bars (i.e. exactly the X drawn on screen).
                const int tPeakAuto = FeatureMarks::reactive_ecg(ecg,
                    (int)std::lround(b.t_begin_auto_ch[c]),
                    (int)std::lround(b.t_end_auto_ch[c])).t_peak;
                const int tPeakUser = FeatureMarks::reactive_ecg(ecg,
                    umk.t_begin_ch[c], umk.t_end_ch[c]).t_peak;

                // Order MUST match ecgPointNames:
                //   p_peak, q_begin, q_peak(computed), r_peak, s_peak(computed),
                //   s_end,  t_peak(reactive), t_begin, t_end
                struct P { int a; int u; };
                const P pts[] = {
                    { (int)std::lround(b.p_peak_auto_ch[c]),  umk.p_peak_ch[c]  },
                    { (int)std::lround(b.q_begin_auto_ch[c]), umk.q_begin_ch[c] },
                    { ftAuto.q_idx,         ftUser.q_idx    },
                    { (int)std::lround(b.r_peak_auto_ch[c]),  b.r_peak_ch[c]  },
                    { ftAuto.s_idx,         ftUser.s_idx    },
                    { (int)std::lround(b.s_end_auto_ch[c]),   umk.s_end_ch[c]   },
                    { tPeakAuto,            tPeakUser       },   // reactive, both sides
                    { (int)std::lround(b.t_begin_auto_ch[c]), umk.t_begin_ch[c] },
                    { (int)std::lround(b.t_end_auto_ch[c]),   umk.t_end_ch[c]   }
                };
                for (int k = 0; k < 9; ++k) {
                    // r_peak (index 3) is the only autodetect-only column --
                    // same rule as emitEcgPointHeader's userToo above.
                    const bool userToo = (k != 3);
                    emitEcgPoint(ecg, pts[k].a, pts[k].u, ref, userToo);
                }
            }

        if (wantPulse) {
            // PPG: onset, p50, peak, dicrotic, peak2, t80, end (matches ppgCols).
            // p50/t80 are reactive -- amplitude crossings bracketed by
            // onset/peak/end -- so each side is derived from its own bar set via
            // the shared FeatureMarks::reactive_ppg, the same call the on-screen
            // glyph makes. The stored ppg_t50/ppg_t80 fields are the detector's
            // originals and are deliberately not used for these columns.
            const FeatureMarks::ReactivePpg rxAuto = FeatureMarks::reactive_ppg(
                b.ppgTemplate, b.ppg_onset_auto, b.ppg_peak_auto, b.ppg_end_auto);
            const FeatureMarks::ReactivePpg rxUser = FeatureMarks::reactive_ppg(
                b.ppgTemplate, b.ppg_onset, b.ppg_peak, b.ppg_end);
            emitPulsePoint(b.ppgTemplate, b.ppg_onset_auto, b.ppg_onset,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            emitPulsePoint(b.ppgTemplate, rxAuto.t50, rxUser.t50,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            emitPulsePoint(b.ppgTemplate, b.ppg_peak_auto, b.ppg_peak,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            emitPulsePoint(b.ppgTemplate, b.ppg_dicrotic_auto, b.ppg_dicrotic,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            emitPulsePoint(b.ppgTemplate, b.ppg_peak2_auto, b.ppg_peak2,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            emitPulsePoint(b.ppgTemplate, rxAuto.t80, rxUser.t80,
                b.ppg_onset_auto, b.ppg_onset, refPpg);
            emitPulsePoint(b.ppgTemplate, b.ppg_end_auto, b.ppg_end,
                b.ppg_onset_auto, b.ppg_onset, refPpg);

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
        // p_wave/q_onset/r_wave ARE the stored autodetect columns -- emitted
        // straight, with no parallel recompute that could disagree with them.
        // (r_wave in particular: R has exactly one definition, r_peak_auto_ch
        // from alignment. It is never re-derived by an argmax here.) Only
        // t_peak is reactive, bracketed by the auto T-begin/T-end to match this
        // group's "_autodetect" name.
        if (wantEcg)
            for (int c = 0; c < 3; ++c) {
                const auto& ecg = chs[c]->ecgTemplate_raw;
                emitAutoFeatPt(ecg, (int)std::lround(b.p_peak_auto_ch[c]));
                emitAutoFeatPt(ecg, (int)std::lround(b.q_begin_auto_ch[c]));
                emitAutoFeatPt(ecg, (int)std::lround(b.r_peak_auto_ch[c]));
                emitAutoFeatPt(ecg, FeatureMarks::reactive_ecg(ecg,
                    (int)std::lround(b.t_begin_auto_ch[c]),
                    (int)std::lround(b.t_end_auto_ch[c])).t_peak);
            }
        if (wantPulse) {
            // PPG glyph values are just the bin's own auto fields now --
            // no separate recompute (see FeatureMarks::detect_ppg_fiducials).
            emitAutoFeatPt(b.ppgTemplate, b.ppg_onset_auto);
            emitAutoFeatPt(b.ppgTemplate, b.ppg_peak_auto);
            emitAutoFeatPt(b.ppgTemplate, b.ppg_dicrotic_auto);
            f << ',' << (b.ppg_dicrotic_found_auto ? 1 : 0);
            emitAutoFeatPt(b.ppgTemplate, b.ppg_end_auto);
            f << ',' << (b.ppg_end_found_auto ? 1 : 0);
        }

        f << '\n';
    }
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

        // Per-anchor ECG user markers: [count] then (tag, 6x3 fields) each.
        // Mirrors the write order. R peak is not stored (auto-only per pass).
        {
            const int nAnchors = r32();
            for (int a = 0; a < nAnchors; ++a) {
                const int tag = r32();
                TemplateBin::MarkerSet m;
                for (int c = 0; c < 3; ++c) m.p_begin_ch[c] = r32();
                for (int c = 0; c < 3; ++c) m.p_peak_ch[c] = r32();
                for (int c = 0; c < 3; ++c) m.q_begin_ch[c] = r32();
                for (int c = 0; c < 3; ++c) m.s_end_ch[c] = r32();
                for (int c = 0; c < 3; ++c) m.t_begin_ch[c] = r32();
                for (int c = 0; c < 3; ++c) m.t_end_ch[c] = r32();
                b.markers_by_anchor[static_cast<AnchorType>(tag)] = m;
            }
        }

        b.ppg_onset = r32();
        b.ppg_t50 = r32();
        b.ppg_peak = r32();
        b.ppg_dicrotic = r32();
        b.ppg_peak2 = r32();
        b.ppg_end = r32();
        b.ppg_t80 = r32();

        b.abp_issue = r8();
        b.abp_onset = r32(); b.abp_peak = r32(); b.abp_dicrotic = r32();
        b.abp_peak2 = r32(); b.abp_end = r32();
        b.art_issue = r8();
        b.art_onset = r32(); b.art_peak = r32(); b.art_dicrotic = r32();
        b.art_peak2 = r32(); b.art_end = r32();
        b.art_pulm_issue = r8();
        b.art_pulm_onset = r32(); b.art_pulm_peak = r32(); b.art_pulm_dicrotic = r32();
        b.art_pulm_peak2 = r32(); b.art_pulm_end = r32();
    }
    return bins;
}