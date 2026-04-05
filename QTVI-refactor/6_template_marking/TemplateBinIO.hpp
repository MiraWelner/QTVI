#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include <cmath>
#include <utility>
#include <fstream>

struct ChannelTemplateData {
    std::vector<double> ecgTemplate_raw;
    std::vector<double> ecgTemplate_squared;
    std::vector<double> ecgTemplate_absval;
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

    // Markings (set during review)
    bool bad_r_ch[3] = { false, false, false };
    uint8_t ppg_issue = 0;   // 0 = ok, 1 = bad, 2 = no ppg
    int dicrotic = -1;
    int onset = -1;
    int peak = -1;
    int end_idx = -1;
};

// ---------------------------------------------------------------------------
// Read helpers
// ---------------------------------------------------------------------------
namespace detail {

    inline uint64_t readU64(std::ifstream& f) {
        uint64_t v = 0; f.read(reinterpret_cast<char*>(&v), 8); return v;
    }
    inline double readF64(std::ifstream& f) {
        double v = 0; f.read(reinterpret_cast<char*>(&v), 8); return v;
    }
    inline uint8_t readU8(std::ifstream& f) {
        uint8_t v = 0; f.read(reinterpret_cast<char*>(&v), 1); return v;
    }
    inline std::vector<double> readDoubleVec(std::ifstream& f) {
        uint64_t sz = readU64(f);
        if (sz == 0) return {};
        std::vector<double> v(sz);
        f.read(reinterpret_cast<char*>(v.data()), sz * 8);
        return v;
    }
    inline std::vector<std::pair<uint64_t, uint64_t>> readPairVec(std::ifstream& f) {
        uint64_t sz = readU64(f);
        if (sz == 0) return {};
        std::vector<std::pair<uint64_t, uint64_t>> v(sz);
        f.read(reinterpret_cast<char*>(v.data()), sz * 16);
        return v;
    }
    inline ChannelTemplateData readChannel(std::ifstream& f) {
        ChannelTemplateData ch;
        ch.ecgTemplate_raw = readDoubleVec(f);
        ch.ecgTemplate_squared = readDoubleVec(f);
        ch.ecgTemplate_absval = readDoubleVec(f);
        ch.alignment_point_raw = readF64(f);
        ch.alignment_point_squared = readF64(f);
        ch.alignment_point_absval = readF64(f);
        ch.avg_r_expand_raw = readF64(f);
        ch.avg_r_expand_squared = readF64(f);
        ch.avg_r_expand_absval = readF64(f);
        return ch;
    }

} // namespace detail

// ---------------------------------------------------------------------------
// Read template_info.bin
// ---------------------------------------------------------------------------
inline std::vector<TemplateBin> readTemplateInfoBin(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) return {};

    uint64_t n = detail::readU64(f);
    std::vector<TemplateBin> bins(n);

    for (uint64_t i = 0; i < n; ++i) {
        auto& b = bins[i];
        b.index = detail::readU64(f);
        b.ppg_bin_indexs = detail::readPairVec(f);
        b.ecg_bin_indexs = detail::readPairVec(f);
        b.bad_segment = detail::readU8(f) != 0;
        b.ch1 = detail::readChannel(f);
        b.ch2 = detail::readChannel(f);
        b.ch3 = detail::readChannel(f);
        b.ppgTemplate = detail::readDoubleVec(f);
    }
    return bins;
}

// ---------------------------------------------------------------------------
// Write template_markings.bin
//
//   uint64  numBins
//   per bin:
//     uint64  index
//     uint8   bad_r_ch1
//     uint8   bad_r_ch2
//     uint8   bad_r_ch3
//     uint8   ppg_issue    (0 = ok, 1 = bad, 2 = no ppg)
//     int32   dicrotic     (-1 = NaN)
//     int32   onset
//     int32   peak
//     int32   end_idx
// ---------------------------------------------------------------------------
inline void writeTemplateMarkingsBin(const std::string& path,
    const std::vector<TemplateBin>& bins) {
    std::ofstream f(path, std::ios::binary);
    if (!f.is_open()) return;

    uint64_t n = bins.size();
    f.write(reinterpret_cast<const char*>(&n), 8);

    for (const auto& b : bins) {
        uint64_t idx = b.index;
        f.write(reinterpret_cast<const char*>(&idx), 8);

        auto w8 = [&](uint8_t v) { f.write(reinterpret_cast<const char*>(&v), 1); };
        w8(b.bad_r_ch[0] ? 1 : 0);
        w8(b.bad_r_ch[1] ? 1 : 0);
        w8(b.bad_r_ch[2] ? 1 : 0);
        w8(b.ppg_issue);

        auto w32 = [&](int v) { int32_t i = v; f.write(reinterpret_cast<const char*>(&i), 4); };
        w32(b.dicrotic);
        w32(b.onset);
        w32(b.peak);
        w32(b.end_idx);
    }
}