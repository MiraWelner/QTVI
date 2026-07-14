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

#include "template_generation\template_io.hpp"

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
    int p_begin_ch[3] = { -1, -1, -1 };
    int q_begin_ch[3] = { -1, -1, -1 };
    int s_end_ch[3] = { -1, -1, -1 };
    int t_begin_ch[3] = { -1, -1, -1 };
    int t_end_ch[3] = { -1, -1, -1 };

    // PPG: sample indices into ppgTemplate. Shared across channels.
    int ppg_onset = -1;
    int ppg_p50 = -1;   // 50% up the upslope, foot -> systolic peak
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
//   per bin (v2 adds the arterial block; v1 stops after the ppg fields):
//     uint64  index
//     uint8   bad_r_ch1, bad_r_ch2, bad_r_ch3
//     uint8   ppg_issue          (0 = ok, 1 = bad, 2 = no ppg)
//     int32   p_begin_ch1..3, q_begin_ch1..3, t_begin_ch1..3, t_end_ch1..3
//     int32   ppg_onset, ppg_peak, ppg_dicrotic, ppg_peak2, ppg_end
//     -- v2 only, one block each for ABP, ART, ART_PULM: --
//     uint8   <chan>_issue
//     int32   <chan>_onset, _peak, _dicrotic, _50, _end
//
// All int32 fields use -1 as the "unmarked / not applicable" sentinel.
// ---------------------------------------------------------------------------
static constexpr uint64_t kTemplateMarkingsV2Sentinel = ~uint64_t(0);
static constexpr uint32_t kTemplateMarkingsVersion = 3;   // v3 adds ppg_p50

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

        for (int c = 0; c < 3; ++c) w32(b.p_begin_ch[c]);
        for (int c = 0; c < 3; ++c) w32(b.q_begin_ch[c]);
        for (int c = 0; c < 3; ++c) w32(b.s_end_ch[c]);
        for (int c = 0; c < 3; ++c) w32(b.t_begin_ch[c]);
        for (int c = 0; c < 3; ++c) w32(b.t_end_ch[c]);

        w32(b.ppg_onset);
        w32(b.ppg_p50);
        w32(b.ppg_peak);
        w32(b.ppg_dicrotic);
        w32(b.ppg_peak2);
        w32(b.ppg_end);

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

inline EcgFeatures computeEcgFeatures(const std::vector<double>& ecg,
    int p_begin, int q_begin, int s_end, int t_begin, int t_end, double rateHz)
{
    EcgFeatures f;
    const int N = static_cast<int>(ecg.size());
    const double msPerSamp = (rateHz > 0.0) ? 1000.0 / rateHz : NAN;
    auto inRange = [&](int i) { return i >= 0 && i < N; };

    if (q_begin >= 0 && s_end >= q_begin) f.qrs_ms = (s_end - q_begin) * msPerSamp;
    if (q_begin >= 0 && t_end >= q_begin) f.qt_ms = (t_end - q_begin) * msPerSamp;
    if (N == 0) return f;

    // isoelectric baseline (median of pre-P segment) � used only to decide QRS polarity
    double baseline = 0.0;
    {
        int hi = std::min(p_begin > 0 ? p_begin : q_begin, N);
        if (hi > 0) {
            std::vector<double> w(ecg.begin(), ecg.begin() + hi);
            std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end());
            baseline = w[w.size() / 2];
        }
    }

    if (q_begin >= 0 && s_end > q_begin && inRange(q_begin) && inRange(s_end)) {
        int rmax = q_begin, rmin = q_begin;
        for (int i = q_begin; i <= s_end; ++i) {
            if (ecg[i] > ecg[rmax]) rmax = i;
            if (ecg[i] < ecg[rmin]) rmin = i;
        }
        const bool up = std::abs(ecg[rmax] - baseline) >= std::abs(ecg[rmin] - baseline);
        f.r_idx = up ? rmax : rmin;

        int q = q_begin;
        for (int i = q_begin; i <= f.r_idx; ++i)
            if (up ? ecg[i] < ecg[q] : ecg[i] > ecg[q]) q = i;
        f.q_idx = q;

        int s = f.r_idx;
        for (int i = f.r_idx; i <= s_end; ++i)
            if (up ? ecg[i] < ecg[s] : ecg[i] > ecg[s]) s = i;
        f.s_idx = s;
    }

    if (t_begin >= 0 && t_end > t_begin && inRange(t_begin)) {
        int t = t_begin;
        const int hi = std::min(t_end, N - 1);
        for (int i = t_begin; i <= hi; ++i)
            if (std::abs(ecg[i] - baseline) > std::abs(ecg[t] - baseline)) t = i;
        f.t_idx = t;
    }
    return f;
}

// Per-channel ECG column order. isUser => straight off a marker (gets _user).
// isInterval => qrs/qt, ms only (no y, no x/user).
struct EcgColSpec { const char* name; bool isUser; bool isInterval; };
static const EcgColSpec ecgCols[] = {
    {"p_begin", true,  false}, {"q_begin", true,  false},
    {"q_peak",  false, false}, {"r_peak",  false, false}, {"s_peak", false, false},
    {"s_end",   true,  false}, {"t_begin", true,  false},
    {"t_peak",  false, false}, {"t_end",   true,  false},
    {"qrs",     false, true},  {"qt",      false, true},
};
inline constexpr const char* ppgCols[] = { "ppg_onset","ppg_p50","ppg_peak","ppg_dicrotic","ppg_peak2","ppg_end" };
inline constexpr const char* abpCols[] = { "abp_onset","abp_peak","abp_dicrotic","abp_peak2","abp_end" };
inline constexpr const char* artCols[] = { "art_onset","art_peak","art_dicrotic","art_peak2","art_end" };
inline constexpr const char* artPulmCols[] = { "art_pulm_onset","art_pulm_peak","art_pulm_dicrotic","art_pulm_peak2","art_pulm_end" };

inline void writeTemplateMarkingsCsv(const std::string& path,
    const std::vector<TemplateBin>& bins,
    const std::string& fileID,
    double sampleRateHz)
{
    std::ofstream f(path);
    if (!f.is_open())
        throw std::runtime_error("cannot open for write: " + path);

    // header (generated so it can't drift from the row loop)
    f << "file_id,bin_index,bad_r_ch1,bad_r_ch2,bad_r_ch3,ppg_issue";
    for (int c = 1; c <= 3; ++c) for (const auto& col : ecgCols) {
        const char* u = col.isUser ? "_user" : "";
        if (col.isInterval) f << ',' << col.name << "_ch" << c << "_ms";
        else {
            f << ',' << col.name << "_ch" << c << "_x_ms" << u;
            f << ',' << col.name << "_ch" << c << "_y_mv" << u;
        }
    }
    // Per pulse marker: x location, raw (pre-normalization) sample value,
    // and PI-normalized value (100*(v-foot)/foot / global_ref).
    for (const char* n : ppgCols)
        f << ',' << n << "_x_ms_user"
        << ',' << n << "_y_mv_user"
        << ',' << n << "_y_norm_user";
    f << ",abp_issue";
    for (const char* n : abpCols)
        f << ',' << n << "_x_ms_user"
        << ',' << n << "_y_mv_user"
        << ',' << n << "_y_norm_user";
    f << ",art_issue";
    for (const char* n : artCols)
        f << ',' << n << "_x_ms_user"
        << ',' << n << "_y_mv_user"
        << ',' << n << "_y_norm_user";
    f << ",art_pulm_issue";
    for (const char* n : artPulmCols)
        f << ',' << n << "_x_ms_user"
        << ',' << n << "_y_mv_user"
        << ',' << n << "_y_norm_user";
    f << '\n';

    const double toMs = 1000.0 / sampleRateHz;

    // Per-channel PI global reference = median across bins of
    //   100 * (peak_y - foot_y) / foot_y
    // Skips bad/absent bins for each channel. Matches the formula in
    // NormalizeFeatures.hpp exactly.
    auto medianFinite = [](std::vector<double> v) -> double {
        v.erase(std::remove_if(v.begin(), v.end(),
            [](double x) { return !std::isfinite(x); }), v.end());
        if (v.empty()) return std::nan("");
        const size_t mid = v.size() / 2;
        std::nth_element(v.begin(), v.begin() + mid, v.end());
        return v[mid];
        };
    auto computeRef = [&](const std::vector<double> TemplateBin::* trace,
        const int TemplateBin::* footIdx, const int TemplateBin::* peakIdx,
        bool includePpgIssue) -> double
        {
            std::vector<double> vals;
            for (const auto& b : bins) {
                if (b.bad_segment) continue;
                if (includePpgIssue && b.ppg_issue != 0) continue;
                const auto& v = b.*trace;
                const int fi = b.*footIdx;
                const int pi = b.*peakIdx;
                if (fi < 0 || pi < 0 || fi >= (int)v.size() || pi >= (int)v.size()) continue;
                const double fy = v[fi], py = v[pi];
                if (std::isnan(fy) || std::isnan(py) || std::abs(fy) < 1e-12) continue;
                vals.push_back(100.0 * (py - fy) / fy);
            }
            return medianFinite(std::move(vals));
        };
    const double refPpg = computeRef(&TemplateBin::ppgTemplate,
        &TemplateBin::ppg_onset, &TemplateBin::ppg_peak, true);
    const double refAbp = computeRef(&TemplateBin::abpTemplate,
        &TemplateBin::abp_onset, &TemplateBin::abp_peak, false);
    const double refArt = computeRef(&TemplateBin::artTemplate,
        &TemplateBin::art_onset, &TemplateBin::art_peak, false);
    const double refArtPulm = computeRef(&TemplateBin::artPulmTemplate,
        &TemplateBin::art_pulm_onset, &TemplateBin::art_pulm_peak, false);

    for (const auto& b : bins) {
        f << fileID << ',' << b.index << ','
            << (b.bad_r_ch[0] ? 1 : 0) << ','
            << (b.bad_r_ch[1] ? 1 : 0) << ','
            << (b.bad_r_ch[2] ? 1 : 0) << ','
            << static_cast<int>(b.ppg_issue);

        // emit x_ms,y_mv for a template point; blank both if idx unmarked/out of range
        auto xy = [&](const std::vector<double>& v, int idx) {
            const bool ok = (idx >= 0 && idx < static_cast<int>(v.size()));
            f << ','; if (ok) f << (idx * toMs);
            f << ','; if (ok && !std::isnan(v[idx])) f << v[idx];
            };
        // emit x_ms, y_mv (pre-norm), y_norm for one pulse marker.
        // footIdx is the same "onset" index the channel uses; ref is the
        // channel's global PI median. Any missing piece leaves y_norm blank.
        auto xyn = [&](const std::vector<double>& v, int idx, int footIdx, double ref) {
            const bool ok = (idx >= 0 && idx < (int)v.size());
            f << ','; if (ok) f << (idx * toMs);
            f << ','; if (ok && !std::isnan(v[idx])) f << v[idx];
            f << ',';
            if (!ok || footIdx < 0 || footIdx >= (int)v.size()) return;
            const double fy = v[footIdx];
            if (std::isnan(fy) || std::abs(fy) < 1e-12) return;
            if (!std::isfinite(ref) || ref == 0.0) return;
            if (std::isnan(v[idx])) return;
            const double local = 100.0 * (v[idx] - fy) / fy;
            f << (local / ref);
            };
        auto ms1 = [&](double v) { f << ','; if (!std::isnan(v)) f << v; };

        const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
        for (int c = 0; c < 3; ++c) {
            const auto& ecg = chs[c]->ecgTemplate_raw;
            EcgFeatures ft = computeEcgFeatures(ecg,
                b.p_begin_ch[c], b.q_begin_ch[c], b.s_end_ch[c],
                b.t_begin_ch[c], b.t_end_ch[c], sampleRateHz);
            const int pts[] = {
                b.p_begin_ch[c], b.q_begin_ch[c], ft.q_idx, ft.r_idx, ft.s_idx,
                b.s_end_ch[c],   b.t_begin_ch[c], ft.t_idx, b.t_end_ch[c]
            };
            for (int idx : pts) xy(ecg, idx);
            ms1(ft.qrs_ms);
            ms1(ft.qt_ms);
        }
        // PPG: onset, P50, peak, dicrotic, peak2, end — matches ppgCols order.
        xyn(b.ppgTemplate, b.ppg_onset, b.ppg_onset, refPpg);
        xyn(b.ppgTemplate, b.ppg_p50, b.ppg_onset, refPpg);
        xyn(b.ppgTemplate, b.ppg_peak, b.ppg_onset, refPpg);
        xyn(b.ppgTemplate, b.ppg_dicrotic, b.ppg_onset, refPpg);
        xyn(b.ppgTemplate, b.ppg_peak2, b.ppg_onset, refPpg);
        xyn(b.ppgTemplate, b.ppg_end, b.ppg_onset, refPpg);

        f << ',' << static_cast<int>(b.abp_issue);
        xyn(b.abpTemplate, b.abp_onset, b.abp_onset, refAbp);
        xyn(b.abpTemplate, b.abp_peak, b.abp_onset, refAbp);
        xyn(b.abpTemplate, b.abp_dicrotic, b.abp_onset, refAbp);
        xyn(b.abpTemplate, b.abp_peak2, b.abp_onset, refAbp);
        xyn(b.abpTemplate, b.abp_end, b.abp_onset, refAbp);

        f << ',' << static_cast<int>(b.art_issue);
        xyn(b.artTemplate, b.art_onset, b.art_onset, refArt);
        xyn(b.artTemplate, b.art_peak, b.art_onset, refArt);
        xyn(b.artTemplate, b.art_dicrotic, b.art_onset, refArt);
        xyn(b.artTemplate, b.art_peak2, b.art_onset, refArt);
        xyn(b.artTemplate, b.art_end, b.art_onset, refArt);

        f << ',' << static_cast<int>(b.art_pulm_issue);
        xyn(b.artPulmTemplate, b.art_pulm_onset, b.art_pulm_onset, refArtPulm);
        xyn(b.artPulmTemplate, b.art_pulm_peak, b.art_pulm_onset, refArtPulm);
        xyn(b.artPulmTemplate, b.art_pulm_dicrotic, b.art_pulm_onset, refArtPulm);
        xyn(b.artPulmTemplate, b.art_pulm_peak2, b.art_pulm_onset, refArtPulm);
        xyn(b.artPulmTemplate, b.art_pulm_end, b.art_pulm_onset, refArtPulm);
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

    std::vector<TemplateBin> bins(n);
    for (uint64_t i = 0; i < n; ++i) {
        auto& b = bins[i];
        f.read(reinterpret_cast<char*>(&b.index), 8);

        b.bad_r_ch[0] = (r8() != 0);
        b.bad_r_ch[1] = (r8() != 0);
        b.bad_r_ch[2] = (r8() != 0);
        b.ppg_issue = r8();

        for (int c = 0; c < 3; ++c) b.p_begin_ch[c] = r32();
        for (int c = 0; c < 3; ++c) b.q_begin_ch[c] = r32();
        for (int c = 0; c < 3; ++c) b.s_end_ch[c] = r32();
        for (int c = 0; c < 3; ++c) b.t_begin_ch[c] = r32();
        for (int c = 0; c < 3; ++c) b.t_end_ch[c] = r32();

        b.ppg_onset = r32();
        if (has_p50) b.ppg_p50 = r32();
        b.ppg_peak = r32();
        b.ppg_dicrotic = r32();
        b.ppg_peak2 = r32();
        b.ppg_end = r32();

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