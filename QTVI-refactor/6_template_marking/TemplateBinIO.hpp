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
    bool templateBad = false;
    bool bad_r = false;
    bool bad_ppg = false;
    int dicrotic = -1;   // sample index, -1 = not set / NaN
    int onset = 0;
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
//     uint8   templateBad
//     uint8   bad_r
//     uint8   bad_ppg
//     int32   dicrotic   (-1 = NaN)
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

        uint8_t tb = b.templateBad ? 1 : 0;
        uint8_t br = b.bad_r ? 1 : 0;
        uint8_t bp = b.bad_ppg ? 1 : 0;
        f.write(reinterpret_cast<const char*>(&tb), 1);
        f.write(reinterpret_cast<const char*>(&br), 1);
        f.write(reinterpret_cast<const char*>(&bp), 1);

        int32_t d = b.dicrotic;
        int32_t o = b.onset;
        int32_t p = b.peak;
        int32_t e = b.end_idx;
        f.write(reinterpret_cast<const char*>(&d), 4);
        f.write(reinterpret_cast<const char*>(&o), 4);
        f.write(reinterpret_cast<const char*>(&p), 4);
        f.write(reinterpret_cast<const char*>(&e), 4);
    }
}