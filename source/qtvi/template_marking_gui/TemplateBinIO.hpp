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

    ChannelTemplateData ch1, ch2, ch3;
    std::vector<double> ppgTemplate;
    std::vector<double> ppgTemplate_std;

    // Markings produced by the viewer. ECG markings are per-channel
    // (the morphology differs by lead). PPG markings are shared across
    // channels of the same bin (the PPG trace is identical in all three
    // channel views).
    bool    bad_r_ch[3] = { false, false, false };
    uint8_t ppg_issue = 0;   // 0 = ok, 1 = bad, 2 = no ppg

    // ECG: per-channel sample indices into the channel's ecgTemplate_raw.
    int p_begin_ch[3] = { -1, -1, -1 };
    int q_begin_ch[3] = { -1, -1, -1 };
    int t_begin_ch[3] = { -1, -1, -1 };
    int t_end_ch[3] = { -1, -1, -1 };

    // PPG: sample indices into ppgTemplate. Shared across channels.
    int ppg_onset = -1;
    int ppg_peak = -1;
    int ppg_dicrotic = -1;
    int ppg_50 = -1;
    int ppg_end = -1;
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
        dst.ppgTemplate = src.ppgTemplate;
        dst.ppgTemplate_std = src.ppgTemplate_std;

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
//   uint64   numBins
//   per bin:
//     uint64  index
//     uint8   bad_r_ch1, bad_r_ch2, bad_r_ch3
//     uint8   ppg_issue          (0 = ok, 1 = bad, 2 = no ppg)
//     int32   p_begin_ch1, p_begin_ch2, p_begin_ch3
//     int32   q_begin_ch1, q_begin_ch2, q_begin_ch3
//     int32   t_begin_ch1, t_begin_ch2, t_begin_ch3
//     int32   t_end_ch1,   t_end_ch2,   t_end_ch3
//     int32   ppg_onset, ppg_peak, ppg_dicrotic, ppg_50, ppg_end
//
// All int32 fields use -1 as the "unmarked / not applicable" sentinel.
// ---------------------------------------------------------------------------
inline void writeTemplateMarkingsBin(const std::string& path,
    const std::vector<TemplateBin>& bins) {
    std::ofstream f(path, std::ios::binary);
    if (!f.is_open())
        throw std::runtime_error("cannot open for write: " + path);

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
        for (int c = 0; c < 3; ++c) w32(b.t_begin_ch[c]);
        for (int c = 0; c < 3; ++c) w32(b.t_end_ch[c]);

        w32(b.ppg_onset);
        w32(b.ppg_peak);
        w32(b.ppg_dicrotic);
        w32(b.ppg_50);
        w32(b.ppg_end);
    }
}

// CSV mirror of writeTemplateMarkingsBin: one row per bin, same fields in the
// same order. Marker sample indices are written as-is; the -1 sentinel
// (unmarked / n/a) becomes an empty cell. If sampleRateHz > 0, a matching
// *_sec column is appended for every marker (sample / rate).
//
// NOTE on _sec: this uses ONE rate for both ECG and PPG markers, matching how
// the viewer itself labels both traces (BinPlotWidget::m_sampleRate). If ECG
// and PPG ever run at different rates, the PPG _sec columns would be off --
// leave sampleRateHz at 0 in that case and work in sample indices.
inline void writeTemplateMarkingsCsv(const std::string& path,
    const std::vector<TemplateBin>& bins,
    double sampleRateHz = 0.0)
{
    std::ofstream f(path);
    if (!f.is_open())
        throw std::runtime_error("cannot open for write: " + path);

    const bool secs = sampleRateHz > 0.0;

    f << "bin_index,bad_r_ch1,bad_r_ch2,bad_r_ch3,ppg_issue,ppg_issue_text,"
        "p_begin_ch1,p_begin_ch2,p_begin_ch3,"
        "q_begin_ch1,q_begin_ch2,q_begin_ch3,"
        "t_begin_ch1,t_begin_ch2,t_begin_ch3,"
        "t_end_ch1,t_end_ch2,t_end_ch3,"
        "ppg_onset,ppg_peak,ppg_dicrotic,ppg_50,ppg_end";
    if (secs)
        f << ",p_begin_ch1_sec,p_begin_ch2_sec,p_begin_ch3_sec,"
        "q_begin_ch1_sec,q_begin_ch2_sec,q_begin_ch3_sec,"
        "t_begin_ch1_sec,t_begin_ch2_sec,t_begin_ch3_sec,"
        "t_end_ch1_sec,t_end_ch2_sec,t_end_ch3_sec,"
        "ppg_onset_sec,ppg_peak_sec,ppg_dicrotic_sec,ppg_50_sec,ppg_end_sec";
    f << '\n';

    auto issueText = [](uint8_t v) -> const char* {
        switch (v) { case 0: return "ok"; case 1: return "bad"; case 2: return "no_ppg"; default: return "unknown"; }
        };

    for (const auto& b : bins) {
        f << b.index << ','
            << (b.bad_r_ch[0] ? 1 : 0) << ','
            << (b.bad_r_ch[1] ? 1 : 0) << ','
            << (b.bad_r_ch[2] ? 1 : 0) << ','
            << static_cast<int>(b.ppg_issue) << ','
            << issueText(b.ppg_issue);

        const int marks[] = {
            b.p_begin_ch[0], b.p_begin_ch[1], b.p_begin_ch[2],
            b.q_begin_ch[0], b.q_begin_ch[1], b.q_begin_ch[2],
            b.t_begin_ch[0], b.t_begin_ch[1], b.t_begin_ch[2],
            b.t_end_ch[0],   b.t_end_ch[1],   b.t_end_ch[2],
            b.ppg_onset, b.ppg_peak, b.ppg_dicrotic, b.ppg_50, b.ppg_end
        };
        for (int m : marks) { f << ','; if (m >= 0) f << m; }          // -1 -> empty
        if (secs) for (int m : marks) { f << ','; if (m >= 0) f << (m / sampleRateHz); }
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

    uint64_t n = 0;
    f.read(reinterpret_cast<char*>(&n), 8);

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
        for (int c = 0; c < 3; ++c) b.t_begin_ch[c] = r32();
        for (int c = 0; c < 3; ++c) b.t_end_ch[c] = r32();

        b.ppg_onset = r32();
        b.ppg_peak = r32();
        b.ppg_dicrotic = r32();
        b.ppg_50 = r32();
        b.ppg_end = r32();
    }
    return bins;
}