#pragma once
//
// cgm_bin_io.hpp
//
// <patientID>_<CGM file name>.bin: one processed CGM file.
// Written to <output_folder>/cgm_output/. No Qt, so it is unit-testable.
//
// LAYOUT, little-endian, packed.
//
//   HEADER (72 bytes)
//     uint32   magic     0x424D4743 ("CGMB" in a hex dump)
//     uint32   version   0
//     uint64   n         number of 15-min grid slots (length of each column)
//     float64  mean      mg/dL                 } the spec's cgmBin over every
//     float64  sd        mg/dL                 } used value in the recording
//     float64  cv        fraction (sd / mean)  } (after exclusion and
//     float64  roc       mg/dL per reading     } single-gap fill).
//     float64  tar       % above 180           } sd and cv are NaN with fewer
//     float64  tbr       % below 70            } than 8 values; all are NaN
//     float64  tir       % 70 to 180           } with fewer than 2.
//   BODY (n entries per column, each column written whole)
//     float64  time[n]         Unix-epoch ms of the slot, CGM device clock
//                              (local time written as if UTC)
//     float64  raw[n]          reading as recorded, mg/dL; NaN if no reading
//     float64  normalized[n]   (value - p2) / (p98 - p2) of the value used for
//                              metrics; NaN where that value is missing
//     uint8    flags[n]        one code per reading:
//                                0 observed, 1 interpolated, 2 missing,
//                                3 out of range, 4 compression
//
// File size is 72 + 25 n bytes. tir + tar + tbr = 100.
//
// A reader that sees a bad magic, another version, or a file whose size does
// not match n rejects it.
//

#include "cgm_pipeline.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace cgm_bin_io {

    using namespace cgm_pipeline;

    inline constexpr std::uint32_t kCgmMagic = 0x424D4743u;
    inline constexpr std::uint32_t kCgmVersion = 0;

    // Per-reading flag codes in the file. The pipeline keeps a bit mask per
    // slot; each slot is in exactly one of these states (a reading cannot be
    // both out of range and compression, since compression ignores
    // out-of-range readings).
    inline constexpr std::uint8_t kFlagObserved = 0;
    inline constexpr std::uint8_t kFlagInterpolated = 1;
    inline constexpr std::uint8_t kFlagMissing = 2;
    inline constexpr std::uint8_t kFlagOutOfRange = 3;
    inline constexpr std::uint8_t kFlagCompression = 4;

    inline std::uint8_t flagCode(std::uint8_t slotFlags) {
        if (slotFlags & compression) return kFlagCompression;
        if (slotFlags & out_of_range_code) return kFlagOutOfRange;
        if (slotFlags & interpolated_code) return kFlagInterpolated;
        if (slotFlags & observed_code) return kFlagObserved;
        return kFlagMissing;
    }

    struct CgmBinFile {
        // header summary, whole recording
        double mean = kNaN, sd = kNaN, cv = kNaN, roc = kNaN;
        double tar = kNaN, tbr = kNaN, tir = kNaN;
        // body columns
        std::vector<double>       time;        // epoch ms
        std::vector<double>       raw;
        std::vector<double>       normalized;
        std::vector<std::uint8_t> flags;       // 0..4, see kFlag*
        std::size_t size() const { return raw.size(); }
    };

    /// Build the file contents from a pipeline result.
    inline CgmBinFile makeCgmBinFile(const CgmPipelineResult& R,
        const CgmParams& p = CgmParams{}) {
        CgmBinFile F;
        const CgmGrid& G = R.grid;
        const std::size_t n = G.size();
        F.time.resize(n);
        F.raw = G.raw;
        F.flags.resize(n);
        for (std::size_t k = 0; k < n; ++k) F.flags[k] = flagCode(G.flags[k]);
        F.normalized.assign(n, kNaN);

        const double span = R.p98 - R.p2;
        const bool canScale = std::isfinite(span) && span > 0.0;
        std::vector<double> used;
        used.reserve(n);
        for (std::size_t k = 0; k < n; ++k) {
            F.time[k] = (G.t0Min + static_cast<double>(k) * G.stepMin) * 60000.0;
            if (!std::isfinite(G.value[k])) continue;
            used.push_back(G.value[k]);
            if (canScale) F.normalized[k] = (G.value[k] - R.p2) / span;
        }

        const CgmBinFeatures f = cgmBin(used);
        if (f.valid) {
            F.mean = f.mean; F.roc = f.roc;
            F.tar = f.tar; F.tbr = f.tbr; F.tir = f.tir;
            if (static_cast<int>(used.size()) >= p.minPointsSd) { F.sd = f.sd; F.cv = f.cv; }
        }
        return F;
    }

    namespace detail {
        template <class T> void put(std::ofstream& f, const T& v) {
            f.write(reinterpret_cast<const char*>(&v), sizeof(T));
        }
        template <class T> bool get(std::ifstream& f, T& v) {
            return static_cast<bool>(f.read(reinterpret_cast<char*>(&v), sizeof(T)));
        }
    }

    namespace detail {
        inline void putCol(std::ofstream& f, const std::vector<double>& v) {
            f.write(reinterpret_cast<const char*>(v.data()),
                static_cast<std::streamsize>(v.size() * sizeof(double)));
        }
        inline bool getCol(std::ifstream& f, std::vector<double>& v, std::uint64_t n) {
            v.resize(n);
            return static_cast<bool>(f.read(reinterpret_cast<char*>(v.data()),
                static_cast<std::streamsize>(n * sizeof(double))));
        }
    }

    inline bool writeCgmBin(const std::string& path, const CgmBinFile& F) {
        using namespace detail;
        const std::size_t n = F.raw.size();
        if (F.time.size() != n || F.normalized.size() != n || F.flags.size() != n)
            return false;
        std::ofstream f(path, std::ios::binary | std::ios::trunc);
        if (!f) return false;
        put<std::uint32_t>(f, kCgmMagic);
        put<std::uint32_t>(f, kCgmVersion);
        put<std::uint64_t>(f, static_cast<std::uint64_t>(n));
        for (double v : { F.mean, F.sd, F.cv, F.roc, F.tar, F.tbr, F.tir }) put<double>(f, v);
        putCol(f, F.time);
        putCol(f, F.raw);
        putCol(f, F.normalized);
        f.write(reinterpret_cast<const char*>(F.flags.data()), static_cast<std::streamsize>(n));
        return static_cast<bool>(f);
    }

    inline bool readCgmBin(const std::string& path, CgmBinFile& F,
        std::string* error = nullptr) {
        using namespace detail;
        auto fail = [&](const char* why) { if (error) *error = why; return false; };
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (!f) return fail("cannot open");
        const auto size = static_cast<std::uint64_t>(f.tellg());
        f.seekg(0);

        std::uint32_t magic = 0, version = 0;
        std::uint64_t n = 0;
        if (!get(f, magic) || magic != kCgmMagic) return fail("not a CGMB file");
        if (!get(f, version) || version != kCgmVersion) return fail("unsupported version");
        if (!get(f, n)) return fail("truncated header");
        if (size < 72 || n >(size - 72) / 25 || size != 72 + 25 * n) return fail("size does not match n");
        if (!(get(f, F.mean) && get(f, F.sd) && get(f, F.cv) && get(f, F.roc) &&
            get(f, F.tar) && get(f, F.tbr) && get(f, F.tir))) return fail("truncated header");

        F.flags.resize(n);
        const bool ok = getCol(f, F.time, n) && getCol(f, F.raw, n) && getCol(f, F.normalized, n) &&
            static_cast<bool>(f.read(reinterpret_cast<char*>(F.flags.data()), static_cast<std::streamsize>(n)));
        if (!ok) return fail("truncated body");
        return true;
    }

}  // namespace cgm_bin_io