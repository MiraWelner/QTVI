#pragma once
//
// Adapter layer: reads the new template_io binary format and projects it
// into the TemplateBin shape the viewer expects. Writing markings stays
// in this file -- markings are a separate file from the template file,
// so write_template_binfile (in template_io) is irrelevant here.
//
// Marker set per bin:
//   ECG (per channel):  P-onset, Q-begin, T-begin, T-end
//   PPG (shared):       Onset, Peak, Dicrotic notch, 50% recovery, End
//
// All marker sample indices use -1 as the "unmarked / not applicable"
// sentinel.
//

#include <vector>
#include <string>
#include <cstdint>
#include <cmath>
#include <utility>
#include <fstream>
#include <stdexcept>

#include "template_io.hpp"

// ---------------------------------------------------------------------------
// Markings file format version.
//
//   v1 (legacy, no magic): 3 ECG markers (Q, Tb, Te) + 2 PPG markers (On, Pk)
//   v2 (current):          4 ECG markers (P, Q, Tb, Te) + 5 PPG markers
//                          (On, Pk, Dc, 50, End)
//
// v2 files start with the 8-byte magic "TMARK\0\0\0" followed by a
// uint32 version, before the bin count. v1 files start straight with
// the uint64 bin count -- they're detected by the missing magic.
// ---------------------------------------------------------------------------
inline constexpr char     kTMarkMagic[8] = { 'T','M','A','R','K',0,0,0 };
inline constexpr uint32_t kTMarkVersion = 2;

// ---------------------------------------------------------------------------
// In-memory model used by the viewer.
// ---------------------------------------------------------------------------

struct ChannelTemplateData {
    std::vector<double> ecgTemplate_raw;
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

    // Markings produced by the viewer. ECG markings are per-channel
    // (the morphology differs by lead). PPG markings are shared across
    // channels of the same bin (the PPG trace is identical in all three
    // channel views).
    bool    bad_r_ch[3] = { false, false, false };
    uint8_t ppg_issue = 0;   // 0 = ok, 1 = bad, 2 = no ppg

    // ECG: per-channel sample indices into the channel's ecgTemplate_raw.
    int p_begin_ch[3] = { -1, -1, -1 };   // NEW
    int q_begin_ch[3] = { -1, -1, -1 };
    int t_begin_ch[3] = { -1, -1, -1 };
    int t_end_ch[3] = { -1, -1, -1 };

    // PPG: sample indices into ppgTemplate. Shared across channels.
    int ppg_onset = -1;
    int ppg_peak = -1;
    int ppg_dicrotic = -1;   // NEW
    int ppg_50 = -1;   // NEW
    int ppg_end = -1;   // NEW
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

        dst.ch1.ecgTemplate_raw = src.ch1_raw.ecgTemplate;
        dst.ch1.alignment_point_raw = src.ch1_raw.alignment_point;
        dst.ch1.avg_r_expand_raw = src.ch1_raw.avg_r_expand;

        dst.ch2.ecgTemplate_raw = src.ch2_raw.ecgTemplate;
        dst.ch2.alignment_point_raw = src.ch2_raw.alignment_point;
        dst.ch2.avg_r_expand_raw = src.ch2_raw.avg_r_expand;

        dst.ch3.ecgTemplate_raw = src.ch3_raw.ecgTemplate;
        dst.ch3.alignment_point_raw = src.ch3_raw.alignment_point;
        dst.ch3.avg_r_expand_raw = src.ch3_raw.avg_r_expand;
    }
    return bins;
}

// ---------------------------------------------------------------------------
// Write template_markings.bin. Format v2:
//
//   char[8]  magic              "TMARK\0\0\0"
//   uint32   version            (= 2)
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

    // Magic + version header (v2).
    f.write(kTMarkMagic, 8);
    f.write(reinterpret_cast<const char*>(&kTMarkVersion), 4);

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

        for (int c = 0; c < 3; ++c) w32(b.p_begin_ch[c]);   // NEW
        for (int c = 0; c < 3; ++c) w32(b.q_begin_ch[c]);
        for (int c = 0; c < 3; ++c) w32(b.t_begin_ch[c]);
        for (int c = 0; c < 3; ++c) w32(b.t_end_ch[c]);

        w32(b.ppg_onset);
        w32(b.ppg_peak);
        w32(b.ppg_dicrotic);   // NEW
        w32(b.ppg_50);         // NEW
        w32(b.ppg_end);        // NEW
    }
}

// ---------------------------------------------------------------------------
// Read template_markings.bin. Handles both v1 (legacy, no magic) and v2.
// v1 files are detected by the absence of the magic header; their
// missing marker fields default to -1.
// ---------------------------------------------------------------------------
inline std::vector<TemplateBin> readTemplateMarkingsBin(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open())
        throw std::runtime_error("cannot open for read: " + path);

    // Peek at the first 8 bytes to detect the magic.
    char head[8] = { 0 };
    f.read(head, 8);
    if (!f) throw std::runtime_error("markings file too short: " + path);

    bool v2 = (std::memcmp(head, kTMarkMagic, 8) == 0);

    auto r8 = [&]() -> uint8_t {
        uint8_t v = 0; f.read(reinterpret_cast<char*>(&v), 1); return v;
        };
    auto r32 = [&]() -> int {
        int32_t v = 0; f.read(reinterpret_cast<char*>(&v), 4); return v;
        };

    uint64_t n = 0;
    uint32_t version = 1;

    if (v2) {
        f.read(reinterpret_cast<char*>(&version), 4);
        f.read(reinterpret_cast<char*>(&n), 8);
    }
    else {
        // v1: the 8 bytes we just read ARE the bin count.
        std::memcpy(&n, head, 8);
    }

    std::vector<TemplateBin> bins(n);
    for (uint64_t i = 0; i < n; ++i) {
        auto& b = bins[i];
        f.read(reinterpret_cast<char*>(&b.index), 8);

        b.bad_r_ch[0] = (r8() != 0);
        b.bad_r_ch[1] = (r8() != 0);
        b.bad_r_ch[2] = (r8() != 0);
        b.ppg_issue = r8();

        if (v2) {
            for (int c = 0; c < 3; ++c) b.p_begin_ch[c] = r32();
        }
        for (int c = 0; c < 3; ++c) b.q_begin_ch[c] = r32();
        for (int c = 0; c < 3; ++c) b.t_begin_ch[c] = r32();
        for (int c = 0; c < 3; ++c) b.t_end_ch[c] = r32();

        b.ppg_onset = r32();
        b.ppg_peak = r32();

        if (v2) {
            b.ppg_dicrotic = r32();
            b.ppg_50 = r32();
            b.ppg_end = r32();
        }
    }
    return bins;
}