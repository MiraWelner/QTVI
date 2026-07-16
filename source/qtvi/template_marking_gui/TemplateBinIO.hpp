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
// Std vectors (ecgTemplate_raw_std per channel, ppgTemplate_std) come
// straight from the template file. Empty std => the widget renders the
// trace without a gray band.
//

#include <vector>
#include <string>
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
    std::vector<double> ecgTemplate_raw_std;
    std::vector<double> ecgTemplate_squared;   // unused by viewer
    std::vector<double> ecgTemplate_absval;    // unused by viewer
    double alignment_point_raw = 0;
    double alignment_point_squared = 0;
    double alignment_point_absval = 0;
    double avg_r_expand_raw = 0;
    double avg_r_expand_squared = 0;
    double avg_r_expand_absval = 0;
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
    std::vector<double> ppgTemplate_std;

    // Foot-anchored arterial background traces (empty when absent).
    std::vector<double> abpTemplate;
    std::vector<double> artTemplate;
    std::vector<double> artPulmTemplate;
    // Per-sample std for each arterial template (empty when absent / not
    // computed). Same length as the matching template when present.
    std::vector<double> abpTemplate_std;
    std::vector<double> artTemplate_std;
    std::vector<double> artPulmTemplate_std;

    // Markings produced by the viewer. ECG markings are per-channel
    // (the morphology differs by lead). PPG markings are shared across
    // channels of the same bin (the PPG trace is identical in all three
    // channel views).
    bool    bad_r_ch[3] = { false, false, false };
    uint8_t ppg_issue = 0;   // 0 = ok, 1 = bad, 2 = no ppg

    // ECG: per-channel sample indices into the channel's ecgTemplate_raw.
    // Six user-placed landmarks per channel, in temporal order:
    //   P peak, Q onset, R peak, S end, T peak, T end.
    // All -1 = unmarked / not applicable.
    int p_peak_ch[3] = { -1, -1, -1 };
    int q_begin_ch[3] = { -1, -1, -1 };
    int r_peak_ch[3] = { -1, -1, -1 };
    int s_end_ch[3] = { -1, -1, -1 };
    int t_peak_ch[3] = { -1, -1, -1 };
    int t_end_ch[3] = { -1, -1, -1 };

    // Auto-detect mirror fields (populated at every loadSubject in
    // seedBinMarkers, NOT serialized to the .bin). These preserve the
    // original auto-detected positions so CSVs can emit both the
    // autodetect and user versions of every column even after the user
    // has dragged markers around.
    int p_peak_auto_ch[3] = { -1, -1, -1 };
    int q_begin_auto_ch[3] = { -1, -1, -1 };
    int r_peak_auto_ch[3] = { -1, -1, -1 };
    int s_end_auto_ch[3] = { -1, -1, -1 };
    int t_peak_auto_ch[3] = { -1, -1, -1 };
    int t_end_auto_ch[3] = { -1, -1, -1 };

    // PPG: sample indices into ppgTemplate. Shared across channels.
    int ppg_onset = -1;
    int ppg_p50 = -1;   // 50% up the upslope, foot -> systolic peak
    int ppg_tac80 = -1; // 80% up the upslope, foot -> systolic peak
    int ppg_peak = -1;
    int ppg_dicrotic = -1;
    int ppg_peak2 = -1;
    int ppg_end = -1;

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
    int ppg_onset_auto = -1, ppg_p50_auto = -1, ppg_peak_auto = -1,
        ppg_dicrotic_auto = -1, ppg_peak2_auto = -1, ppg_end_auto = -1;
    int ppg_tac80_auto = -1;
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
inline std::vector<TemplateBin> readTemplateInfoBin(const std::string& path) {
    template_io::TemplateFile tf = template_io::read_template_binfile(path);

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
        dst.ppgTemplate = src.ppgTemplate;
        dst.ppgTemplate_std = src.ppgTemplate_std;
        dst.abpTemplate = src.abpTemplate;
        dst.artTemplate = src.artTemplate;
        dst.artPulmTemplate = src.artPulmTemplate;
        dst.abpTemplate_std = src.abpTemplate_std;
        dst.artTemplate_std = src.artTemplate_std;
        dst.artPulmTemplate_std = src.artPulmTemplate_std;

        dst.ch1.ecgTemplate_raw = src.ch1_raw.ecgTemplate;
        dst.ch1.ecgTemplate_raw_std = src.ch1_raw.ecgTemplate_std;
        dst.ch1.alignment_point_raw = src.ch1_raw.alignment_point;
        dst.ch1.avg_r_expand_raw = src.ch1_raw.avg_r_expand;

        dst.ch2.ecgTemplate_raw = src.ch2_raw.ecgTemplate;
        dst.ch2.ecgTemplate_raw_std = src.ch2_raw.ecgTemplate_std;
        dst.ch2.alignment_point_raw = src.ch2_raw.alignment_point;
        dst.ch2.avg_r_expand_raw = src.ch2_raw.avg_r_expand;

        dst.ch3.ecgTemplate_raw = src.ch3_raw.ecgTemplate;
        dst.ch3.ecgTemplate_raw_std = src.ch3_raw.ecgTemplate_std;
        dst.ch3.alignment_point_raw = src.ch3_raw.alignment_point;
        dst.ch3.avg_r_expand_raw = src.ch3_raw.avg_r_expand;
    }
    return bins;
}

// ---------------------------------------------------------------------------
// template_markings.bin layout:
//
//   VERSIONING: v1 files begin directly with uint64 numBins. v2 files begin
//   with a sentinel (uint64 = 0xFFFFFFFFFFFFFFFF, an impossible bin count)
//   followed by uint32 version (=2), then uint64 numBins. readTemplateMarkingsBin
//   detects which by peeking the first 8 bytes, so old files still load.
//
//   per bin (v4 layout; v1..v3 handled as legacy reads):
//     uint64  index
//     uint8   bad_r_ch1, bad_r_ch2, bad_r_ch3
//     uint8   ppg_issue          (0 = ok, 1 = bad, 2 = no ppg)
//     int32   p_peak_ch1..3, q_begin_ch1..3, r_peak_ch1..3,
//             s_end_ch1..3, t_peak_ch1..3, t_end_ch1..3
//     int32   ppg_onset, ppg_p50, ppg_peak, ppg_dicrotic, ppg_peak2, ppg_end
//     -- arterial (v2+), one block each for ABP, ART, ART_PULM: --
//     uint8   <chan>_issue
//     int32   <chan>_onset, _peak, _dicrotic, _peak2, _end
//
// All int32 fields use -1 as the "unmarked / not applicable" sentinel.
// ---------------------------------------------------------------------------
static constexpr uint64_t kTemplateMarkingsV2Sentinel = ~uint64_t(0);
static constexpr uint32_t kTemplateMarkingsVersion = 5;   // v5 adds ppg_tac80 (P80)

inline void writeTemplateMarkingsBin(const std::string& path,
    const std::vector<TemplateBin>& bins) {
    std::ofstream f(path, std::ios::binary);
    if (!f.is_open())
        throw std::runtime_error("cannot open for write: " + path);

    // v2 header: sentinel + version, then the real bin count.
    uint64_t sentinel = kTemplateMarkingsV2Sentinel;
    f.write(reinterpret_cast<const char*>(&sentinel), 8);
    uint32_t ver = kTemplateMarkingsVersion;
    f.write(reinterpret_cast<const char*>(&ver), 4);

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
        w8(b.ppg_issue);

        // v4 ECG markers: six per channel in temporal order.
        for (int c = 0; c < 3; ++c) w32(b.p_peak_ch[c]);
        for (int c = 0; c < 3; ++c) w32(b.q_begin_ch[c]);
        for (int c = 0; c < 3; ++c) w32(b.r_peak_ch[c]);
        for (int c = 0; c < 3; ++c) w32(b.s_end_ch[c]);
        for (int c = 0; c < 3; ++c) w32(b.t_peak_ch[c]);
        for (int c = 0; c < 3; ++c) w32(b.t_end_ch[c]);

        w32(b.ppg_onset);
        w32(b.ppg_p50);
        w32(b.ppg_peak);
        w32(b.ppg_dicrotic);
        w32(b.ppg_peak2);
        w32(b.ppg_end);
        w32(b.ppg_tac80);   // v5: P80

        // v2 arterial block: ABP, ART, ART_PULM (issue + 5 indices each).
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
    int q_idx = -1, r_idx = -1, s_idx = -1, t_idx = -1;   // peak sample positions
    double qrs_ms = NAN, qt_ms = NAN;
};

inline EcgFeatures computeEcgFeatures(const std::vector<double>& ecg, int p_peak, int q_begin, int r_peak, int s_end, int t_peak, int t_end, double rateHz)
{
    EcgFeatures f;
    const int N = static_cast<int>(ecg.size());
    const double msPerSamp = (rateHz > 0.0) ? 1000.0 / rateHz : NAN;
    auto inRange = [&](int i) { return i >= 0 && i < N; };

    if (q_begin >= 0 && s_end >= q_begin) f.qrs_ms = (s_end - q_begin) * msPerSamp;
    if (q_begin >= 0 && t_end >= q_begin) f.qt_ms = (t_end - q_begin) * msPerSamp;

    if (inRange(r_peak)) f.r_idx = r_peak;
    if (inRange(t_peak)) f.t_idx = t_peak;
    f.q_idx = FeatureMarks::compute_q_peak(ecg, q_begin, r_peak);
    f.s_idx = FeatureMarks::compute_s_peak(ecg, r_peak, s_end);
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
    {"s_end",   true,  false}, {"t_peak",  true,  false}, {"t_end",  true,  false},
    {"qrs",     false, true},  {"qt",      false, true},
};
inline constexpr const char* ppgCols[] = { "ppg_onset","ppg_p50","ppg_peak","ppg_dicrotic","ppg_peak2","ppg_p80","ppg_end" };
inline constexpr const char* abpCols[] = { "abp_onset","abp_peak","abp_dicrotic","abp_peak2","abp_end" };
inline constexpr const char* artCols[] = { "art_onset","art_peak","art_dicrotic","art_peak2","art_end" };
inline constexpr const char* artPulmCols[] = { "art_pulm_onset","art_pulm_peak","art_pulm_dicrotic","art_pulm_peak2","art_pulm_end" };

// ---------------------------------------------------------------------------
// writeTemplateMarkingsCsv
//
// For every marker we emit six columns, in this order:
//   {name}_y_mv_normalized_autodetect
//   {name}_y_mv_normalized_user
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
inline void writeTemplateMarkingsCsv(const std::string& path,
    const std::vector<TemplateBin>& bins,
    const std::string& fileID,
    double sampleRateHz)
{
    std::ofstream f(path);
    if (!f.is_open())
        throw std::runtime_error("cannot open for write: " + path);

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
    for (int c = 0; c < 3; ++c) {
        std::vector<double> vals;
        for (const auto& b : bins) {
            if (b.bad_segment || b.bad_r_ch[c]) continue;
            const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
            const auto& ecg = chs[c]->ecgTemplate_raw;
            if (ecg.empty()) continue;
            EcgFeatures ft = computeEcgFeatures(ecg,
                b.p_peak_auto_ch[c], b.q_begin_auto_ch[c], b.r_peak_auto_ch[c],
                b.s_end_auto_ch[c], b.t_peak_auto_ch[c], b.t_end_auto_ch[c],
                sampleRateHz);
            if (ft.r_idx < 0 || ft.s_idx < 0) continue;
            if (ft.r_idx >= (int)ecg.size() || ft.s_idx >= (int)ecg.size()) continue;
            const double ry = ecg[ft.r_idx], sy = ecg[ft.s_idx];
            if (std::isnan(ry) || std::isnan(sy)) continue;
            vals.push_back(std::abs(ry) + std::abs(sy));
        }
        ecgRef[c] = medianFinite(std::move(vals));
    }

    auto pulseRefAuto = [&](const std::vector<double> TemplateBin::* trace,
        const int TemplateBin::* footAuto,
        const int TemplateBin::* peakAuto,
        bool checkPpgIssue) -> double
        {
            std::vector<double> vals;
            for (const auto& b : bins) {
                if (b.bad_segment) continue;
                if (checkPpgIssue && b.ppg_issue != 0) continue;
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
    const double refPpg = pulseRefAuto(&TemplateBin::ppgTemplate,
        &TemplateBin::ppg_onset_auto, &TemplateBin::ppg_peak_auto, true);
    const double refAbp = pulseRefAuto(&TemplateBin::abpTemplate,
        &TemplateBin::abp_onset_auto, &TemplateBin::abp_peak_auto, false);
    const double refArt = pulseRefAuto(&TemplateBin::artTemplate,
        &TemplateBin::art_onset_auto, &TemplateBin::art_peak_auto, false);
    const double refArtPulm = pulseRefAuto(&TemplateBin::artPulmTemplate,
        &TemplateBin::art_pulm_onset_auto, &TemplateBin::art_pulm_peak_auto, false);

    // ---- header ------------------------------------------------------------
    f << "file_id,bin_index,bad_r_ch1,bad_r_ch2,bad_r_ch3,ppg_issue";

    // ECG 8 point columns + 2 interval columns per channel.
    static const char* ecgPointNames[] = {
        "p_peak", "q_begin", "q_peak", "r_peak", "s_peak",
        "s_end",  "t_peak",  "t_end"
    };
    static const char* ecgIntervalNames[] = { "qrs", "qt" };

    auto emitEcgPointHeader = [&](const char* name, int c) {
        const bool userToo = (std::strcmp(name, "r_peak") != 0);
        f << ',' << name << "_ch" << c << "_y_mv_normalized_autodetect";
        if (userToo) f << ',' << name << "_ch" << c << "_y_mv_normalized_user";
        f << ',' << name << "_ch" << c << "_y_mv_raw_autodetect";
        if (userToo) f << ',' << name << "_ch" << c << "_y_mv_raw_user";
        f << ',' << name << "_ch" << c << "_x_ms_autodetect";
        if (userToo) f << ',' << name << "_ch" << c << "_x_ms_user";
        };
    auto emitIntervalHeader = [&](const char* name, int c) {
        f << ',' << name << "_ch" << c << "_ms_autodetect"
            << ',' << name << "_ch" << c << "_ms_user";
        };
    for (int c = 1; c <= 3; ++c) {
        for (const char* n : ecgPointNames)     emitEcgPointHeader(n, c);
        for (const char* n : ecgIntervalNames)  emitIntervalHeader(n, c);
    }

    // Pulse: 6 cols per marker.
    auto emitPulsePointHeader = [&](const char* name) {
        f << ',' << name << "_y_mv_normalized_autodetect"
            << ',' << name << "_y_mv_normalized_user"
            << ',' << name << "_y_mv_raw_autodetect"
            << ',' << name << "_y_mv_raw_user"
            << ',' << name << "_x_ms_autodetect"
            << ',' << name << "_x_ms_user";
        };
    for (const char* n : ppgCols)     emitPulsePointHeader(n);
    f << ",abp_issue";
    for (const char* n : abpCols)     emitPulsePointHeader(n);
    f << ",art_issue";
    for (const char* n : artCols)     emitPulsePointHeader(n);
    f << ",art_pulm_issue";
    for (const char* n : artPulmCols) emitPulsePointHeader(n);
    // Computed glyph positions (derived from the USER markers).
    auto emitGlyphHeader = [&](const char* name) {
        f << ',' << name << "_x_ms" << ',' << name << "_y_mv_raw";
        };
    for (int c = 1; c <= 3; ++c) {
        char b[64];
        for (const char* g : { "p_wave_glyph", "q_onset_glyph", "r_wave_glyph",
                               "s_end_glyph", "t_peak_glyph", "t_end_glyph" }) {
            std::snprintf(b, sizeof b, "%s_ch%d", g, c);
            emitGlyphHeader(b);
        }
    }
    for (const char* g : { "ppg_foot_glyph", "ppg_p1_glyph",
                           "ppg_dic_glyph", "ppg_p2_glyph" })
        emitGlyphHeader(g);
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

    for (const auto& b : bins) {
        f << fileID << ',' << b.index << ','
            << (b.bad_r_ch[0] ? 1 : 0) << ','
            << (b.bad_r_ch[1] ? 1 : 0) << ','
            << (b.bad_r_ch[2] ? 1 : 0) << ','
            << static_cast<int>(b.ppg_issue);

        const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
        for (int c = 0; c < 3; ++c) {
            const auto& ecg = chs[c]->ecgTemplate_raw;
            const double ref = ecgRef[c];

            EcgFeatures ftAuto = computeEcgFeatures(ecg,
                b.p_peak_auto_ch[c], b.q_begin_auto_ch[c], b.r_peak_auto_ch[c],
                b.s_end_auto_ch[c], b.t_peak_auto_ch[c], b.t_end_auto_ch[c],
                sampleRateHz);
            EcgFeatures ftUser = computeEcgFeatures(ecg,
                b.p_peak_ch[c], b.q_begin_ch[c], b.r_peak_ch[c],
                b.s_end_ch[c], b.t_peak_ch[c], b.t_end_ch[c],
                sampleRateHz);

            // Order MUST match ecgPointNames:
            //   p_peak, q_begin, q_peak(computed), r_peak, s_peak(computed),
            //   s_end,  t_peak, t_end
            struct P { int a; int u; };
            const P pts[] = {
                { b.p_peak_auto_ch[c],  b.p_peak_ch[c]  },
                { b.q_begin_auto_ch[c], b.q_begin_ch[c] },
                { ftAuto.q_idx,         ftUser.q_idx    },
                { b.r_peak_auto_ch[c],  b.r_peak_ch[c]  },
                { ftAuto.s_idx,         ftUser.s_idx    },
                { b.s_end_auto_ch[c],   b.s_end_ch[c]   },
                { b.t_peak_auto_ch[c],  b.t_peak_ch[c]  },
                { b.t_end_auto_ch[c],   b.t_end_ch[c]   }
            };
            for (int k = 0; k < 8; ++k) {
                const bool userToo = (k != 3);   // ecgPointNames[3] == "r_peak"
                emitEcgPoint(ecg, pts[k].a, pts[k].u, ref, userToo);
            }
            // Intervals: qrs, qt (order matches ecgIntervalNames).
            emitIntervalPair(ftAuto.qrs_ms, ftUser.qrs_ms);
            emitIntervalPair(ftAuto.qt_ms, ftUser.qt_ms);
        }

        // PPG: onset, p50, peak, dicrotic, peak2, end (matches ppgCols).
        emitPulsePoint(b.ppgTemplate, b.ppg_onset_auto, b.ppg_onset,
            b.ppg_onset_auto, b.ppg_onset, refPpg);
        emitPulsePoint(b.ppgTemplate, b.ppg_p50_auto, b.ppg_p50,
            b.ppg_onset_auto, b.ppg_onset, refPpg);
        emitPulsePoint(b.ppgTemplate, b.ppg_peak_auto, b.ppg_peak,
            b.ppg_onset_auto, b.ppg_onset, refPpg);
        emitPulsePoint(b.ppgTemplate, b.ppg_dicrotic_auto, b.ppg_dicrotic,
            b.ppg_onset_auto, b.ppg_onset, refPpg);
        emitPulsePoint(b.ppgTemplate, b.ppg_peak2_auto, b.ppg_peak2,
            b.ppg_onset_auto, b.ppg_onset, refPpg);
        emitPulsePoint(b.ppgTemplate, b.ppg_tac80_auto, b.ppg_tac80,
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

        // ---- Computed glyph positions (from the USER markers) -------------
        auto emitGlyphPt = [&](const std::vector<double>& sig, int idx) {
            f << ',';   if (idx >= 0) f << (idx * toMs);
            f << ',';   if (idx >= 0 && idx < (int)sig.size() && !std::isnan(sig[idx])) f << sig[idx];
            };
        for (int c = 0; c < 3; ++c) {
            const auto& ecg = chs[c]->ecgTemplate_raw;
            const FeatureMarks::EcgGlyphs gl = FeatureMarks::compute_ecg_glyphs(
                ecg, b.p_peak_ch[c], b.q_begin_ch[c], b.s_end_ch[c],
                b.t_peak_ch[c], b.t_end_ch[c], sampleRateHz);
            emitGlyphPt(ecg, gl.p_wave);
            emitGlyphPt(ecg, gl.q_onset);
            emitGlyphPt(ecg, gl.r_wave);
            emitGlyphPt(ecg, gl.s_end);
            emitGlyphPt(ecg, gl.t_peak);
            emitGlyphPt(ecg, gl.t_end);
        }
        {
            const FeatureMarks::PpgGlyphs pgl = FeatureMarks::compute_ppg_glyphs(
                b.ppgTemplate, b.ppg_onset, b.ppg_dicrotic, b.ppg_peak2);
            emitGlyphPt(b.ppgTemplate, pgl.foot);
            emitGlyphPt(b.ppgTemplate, pgl.p1);
            emitGlyphPt(b.ppgTemplate, pgl.dic);
            emitGlyphPt(b.ppgTemplate, pgl.p2);
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

    // Peek the first 8 bytes. v2+ begin with the sentinel; anything else
    // is a v1 file whose first 8 bytes are the real bin count.
    uint64_t first = 0;
    f.read(reinterpret_cast<char*>(&first), 8);

    uint32_t version = 1;
    uint64_t n = 0;
    if (first == kTemplateMarkingsV2Sentinel) {
        f.read(reinterpret_cast<char*>(&version), 4);   // 2 or 3
        f.read(reinterpret_cast<char*>(&n), 8);
    }
    else {
        n = first;   // v1: first field WAS the bin count
    }
    const bool has_arterial = (version >= 2);
    const bool has_p50 = (version >= 3);
    const bool has_ecg_peaks = (version >= 4);   // v4: six groups (peaks), not five
    const bool has_tac80 = (version >= 5);   // v5: ppg_tac80 (P80)

    std::vector<TemplateBin> bins(n);
    for (uint64_t i = 0; i < n; ++i) {
        auto& b = bins[i];
        f.read(reinterpret_cast<char*>(&b.index), 8);

        b.bad_r_ch[0] = (r8() != 0);
        b.bad_r_ch[1] = (r8() != 0);
        b.bad_r_ch[2] = (r8() != 0);
        b.ppg_issue = r8();

        if (has_ecg_peaks) {
            // v4: P peak, Q begin, R peak, S end, T peak, T end.
            for (int c = 0; c < 3; ++c) b.p_peak_ch[c] = r32();
            for (int c = 0; c < 3; ++c) b.q_begin_ch[c] = r32();
            for (int c = 0; c < 3; ++c) b.r_peak_ch[c] = r32();
            for (int c = 0; c < 3; ++c) b.s_end_ch[c] = r32();
            for (int c = 0; c < 3; ++c) b.t_peak_ch[c] = r32();
            for (int c = 0; c < 3; ++c) b.t_end_ch[c] = r32();
        }
        else {
            // v1-v3 legacy: P onset, Q begin, S end, T begin, T end. We
            // drop P onset and T begin (they don't map to any current
            // marker) and leave the new peak fields at -1 for re-seeding.
            int scratch[3];
            for (int c = 0; c < 3; ++c) scratch[c] = r32();   // was p_begin
            for (int c = 0; c < 3; ++c) b.q_begin_ch[c] = r32();
            for (int c = 0; c < 3; ++c) b.s_end_ch[c] = r32();
            for (int c = 0; c < 3; ++c) scratch[c] = r32();   // was t_begin
            for (int c = 0; c < 3; ++c) b.t_end_ch[c] = r32();
            (void)scratch;
        }

        b.ppg_onset = r32();
        if (has_p50) b.ppg_p50 = r32();
        b.ppg_peak = r32();
        b.ppg_dicrotic = r32();
        b.ppg_peak2 = r32();
        b.ppg_end = r32();
        if (has_tac80) b.ppg_tac80 = r32();

        if (has_arterial) {
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
        // v1: arterial fields keep their default (-1 / issue 0).
    }
    return bins;
}