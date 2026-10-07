#pragma once
//
// cgm_pipeline.hpp
//
// Continuous glucose monitor (CGM) binning with sensor warm-up exclusion, gap
// handling, artifact flagging and standard glycemic metrics. Header-only.
//
// Input is the FreeStyle Libre (Pro) patient-report CSV: two preamble rows, a
// column-header row starting "Device,Serial Number,Device Timestamp,...", then
// one row per record. Only Record Type 0 ("Historic Glucose mg/dL", the
// 15-minute stored trace) is used; scans, strips, notes and insulin rows are
// skipped. Both timestamp styles seen in our exports parse:
//     "09-10-2026 01:21 PM"   (MM-DD-YYYY hh:mm AM/PM)
//     "7/9/2026 11:03"        (M/D/YYYY HH:MM, 24 h)
//
// STAGES (runCgmPipeline), in this order -- the order matters:
//
//   1. WARM-UP.     Drop every reading up to 12 h after the FIRST reading
//                   of its sensor (keyed by Serial Number, so a mid-recording
//                   sensor change gets its own warm-up). A reading at exactly
//                   +12:00 is dropped too; the first kept reading is after it.
//   2. GRID.        Snap readings onto a regular 15-min grid anchored at the
//                   first post-warm-up reading. Empty slots are NaN.
//   3. FLAGS.       Out-of-range: < 40 or > 500 mg/dL.
//                   Compression: a DROP below 55 mg/dL (previous reading >= 55)
//                   lasting MORE than 15 min (>= 2 consecutive readings),
//                   followed within 30 min by a RAPID RETURN
//                   (>= 2 mg/dL/min from the last low reading, to >= 70).
//                   Flagged readings are EXCLUDED (set missing) by default.
//   4. GAPS.        A single missing slot with valid neighbours on both sides
//                   is linearly interpolated. Runs of 2+ stay missing.
//                   Excluded (flagged) readings are NEVER interpolated over: a
//                   lone sub-40 reading may be a real severe hypo and must not be
//                   replaced by a plausible-looking midpoint.
//   5. SCALING.     Individual 2nd / 98th percentile of the recording's own
//                   OBSERVED, unflagged readings (interpolated points are not
//                   used to set the scale). gn = (g - p2) / (p98 - p2).
//                   Not clipped, so ~4% of points fall outside [0, 1].
//   6. BINS.        Fixed bins of binMinutes from grid slot 0. Per bin:
//                     valid   needs >= 2 points (30 min) AND <= 20% missing
//                     sd, cv  need  >= 8 points (2 h), else NaN (sdValid=false)
//                   The spec's cgmBin computes the metrics; scoreBin feeds it
//                   the present readings and applies the rules around it.
//                   Metrics on RAW mg/dL: mean, sd, cv, roc, tir/tar/tbr.
//                   Metrics on NORMALIZED glucose: nmean, nsd, nroc.
//                   Clinical thresholds (70, 180) ALWAYS use raw glucose.
//
// DEFINITIONS
//   tir  % of valid points with 70 <= g <= 180   (both edges are in range)
//   tar  % with g > 180
//   tbr  % with g < 70
//   tir + tar + tbr == 100 for every valid bin (same denominator, disjoint).
//   cv   sd / mean as a FRACTION (x100 for the usual clinical %).
//   roc  spec formula on the present readings: (last - first) / (n - 1),
//        mg/dL per reading. With no empty slots inside the bin that is per
//        15 min; an empty slot inside the bin stretches each step to 30 min.
//   "missing" = slots still NaN AFTER single-gap interpolation, i.e. truly
//   absent runs of 2+ plus excluded flagged readings. An interpolated point
//   counts as present.
//
// BIN SIZE. Use a CGM-specific bin length, NOT config_entry::bin_size_minutes
// (that one is the ECG/PPG bin, typically a few minutes, which would make every
// CGM bin shorter than the 30-min floor). SD/CV can only ever be valid when
// binMinutes >= 120.
//

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <istream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace cgm_pipeline {

    inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // ---- Parameters ------------------------------------------------------
    struct CgmParams {
        double sampleMinutes = 15.0;     // Libre historic cadence
        double warmupHours = 12.0;
        double binMinutes = 60.0;

        int    minPointsSd = 8;          // 2 h at 15-min sampling
        double maxMissingFrac = 0.20;    // bin invalid if missing > 20%

        double oorLow = 40.0;            // flag g < oorLow
        double oorHigh = 500.0;          // flag g > oorHigh

        double compLow = 55.0;           // compression: below this...
        double compMinMinutes = 15.0;    // ...for MORE than this...
        double compReturnWindowMinutes = 30.0;  // ...then back within this
        double compReturnRate = 2.0;     // mg/dL/min, from the last low reading
        double compReturnLevel = 70.0;   // and back to at least this

        bool   excludeOutOfRange = true;
        bool   excludeCompression = true;

        double pLow = 2.0;               // normalization percentiles
        double pHigh = 98.0;
    };

    // ---- Spec per-bin feature struct and function -----------------------
    // Exactly as given in the spec. cgmBin takes PRESENT readings only (no
    // NaN); scoreBin() below strips empty slots before calling it, and applies
    // the 8-point SD/CV rule and the 20% missing rule afterwards.
    struct CgmBinFeatures { double mean, sd, cv, roc, tir, tar, tbr; bool valid; };
    inline CgmBinFeatures cgmBin(const std::vector<double>& g /*mg/dL, 15-min*/) {
        CgmBinFeatures f{}; if (g.size() < 2) { f.valid = false; return f; }
        double m = 0; for (double x : g) m += x; m /= g.size();
        double s = 0; for (double x : g) s += (x - m) * (x - m); s = std::sqrt(s / (g.size() - 1));
        int in = 0, ab = 0, be = 0; for (double x : g) { if (x > 180)++ab; else if (x < 70)++be; else ++in; }
        f.mean = m; f.sd = s; f.cv = s / m; f.roc = (g.back() - g.front()) / (g.size() - 1);
        f.tir = 100.0 * in / g.size(); f.tar = 100.0 * ab / g.size(); f.tbr = 100.0 * be / g.size();
        f.valid = true;
        return f;
    }

    // ---- Readings --------------------------------------------------------
    struct CgmReading {
        double      tMin = 0.0;       // minutes since 1970-01-01 (device clock)
        double      mgdl = kNaN;      // mg/dL as exported
        std::string serial;           // sensor serial, keys warm-up
    };

    // ---- Grid ------------------------------------------------------------
    enum flag : std::uint8_t {
        observed_code = 1 << 0,  // a reading landed in this slot
        interpolated_code = 1 << 1,  // single-gap fill
        out_of_range_code = 1 << 2,  // < 40 or > 500
        compression = 1 << 3,  // probable compression artifact
        excluded = 1 << 4,  // flagged and removed from metrics
    };

    struct CgmGrid {
        double t0Min = 0.0;                // time of slot 0
        double stepMin = 15.0;
        std::vector<double>       raw;     // as read (NaN if empty / token)
        std::vector<double>       value;   // after exclusion + interpolation
        std::vector<std::uint8_t> flags;
        std::size_t size() const { return value.size(); }
    };

    struct CgmBin {
        int    index = 0;
        double startMin = 0.0;        // grid time of the bin's first slot
        int    expected = 0;          // slots the bin spans
        int    nObserved = 0;
        int    nInterpolated = 0;
        int    nFlagged = 0;          // out-of-range or compression
        int    nMissing = 0;          // expected - finite after interpolation
        double missingFrac = 0.0;
        int    n = 0;                 // present readings passed to cgmBin
        bool   sdValid = false;       // valid and n >= 8; else raw.sd/raw.cv are NaN
        CgmBinFeatures raw{};         // spec metrics on mg/dL; raw.valid is THE flag
        double nmean = kNaN, nsd = kNaN, nroc = kNaN;   // on normalized glucose
    };

    struct CgmPipelineResult {
        CgmGrid             grid;
        std::vector<CgmBin> bins;
        double p2 = kNaN, p98 = kNaN;
        int nInput = 0;
        int nWarmupDropped = 0;
        int nDuplicateSlot = 0;       // second reading in an occupied slot
        int nOutOfRange = 0;
        int nCompression = 0;
        int nInterpolated = 0;
    };

    // ======================================================================
    //  Stage 1: warm-up
    // ======================================================================
    inline std::vector<CgmReading> dropWarmup(const std::vector<CgmReading>& in,
        const CgmParams& p, int* nDropped = nullptr) {
        std::map<std::string, double> firstBySerial;
        for (const auto& r : in) {
            auto it = firstBySerial.find(r.serial);
            if (it == firstBySerial.end() || r.tMin < it->second)
                firstBySerial[r.serial] = r.tMin;
        }
        const double warm = p.warmupHours * 60.0;
        std::vector<CgmReading> out;
        out.reserve(in.size());
        int dropped = 0;
        for (const auto& r : in) {
            if (warm > 0.0 && r.tMin - firstBySerial[r.serial] <= warm) { ++dropped; continue; }
            out.push_back(r);
        }
        if (nDropped) *nDropped = dropped;
        return out;
    }

    // ======================================================================
    //  Stage 2: grid
    // ======================================================================
    inline CgmGrid buildGrid(std::vector<CgmReading> rs, const CgmParams& p,
        int* nDuplicate = nullptr) {
        CgmGrid G;
        G.stepMin = p.sampleMinutes;
        int dup = 0;
        if (rs.empty()) { if (nDuplicate) *nDuplicate = 0; return G; }
        std::stable_sort(rs.begin(), rs.end(),
            [](const CgmReading& a, const CgmReading& b) { return a.tMin < b.tMin; });
        G.t0Min = rs.front().tMin;
        const long long nSlots =
            std::llround((rs.back().tMin - G.t0Min) / G.stepMin) + 1;
        G.raw.assign(static_cast<std::size_t>(nSlots), kNaN);
        G.flags.assign(static_cast<std::size_t>(nSlots), 0);
        for (const auto& r : rs) {
            const auto k = static_cast<std::size_t>(
                std::llround((r.tMin - G.t0Min) / G.stepMin));
            if (G.flags[k] & observed_code) { ++dup; continue; }  // keep first
            if (!std::isfinite(r.mgdl)) continue;
            G.raw[k] = r.mgdl;
            G.flags[k] |= observed_code;
        }
        G.value = G.raw;
        if (nDuplicate) *nDuplicate = dup;
        return G;
    }

    // ======================================================================
    //  Stage 3: flags
    // ======================================================================
    inline int flagOutOfRange(CgmGrid& G, const CgmParams& p) {
        int n = 0;
        for (std::size_t i = 0; i < G.size(); ++i) {
            const double x = G.raw[i];
            if (std::isfinite(x) && (x < p.oorLow || x > p.oorHigh))
                G.flags[i] |= out_of_range_code;
            if (G.flags[i] & out_of_range_code) ++n;
        }
        return n;
    }

    /// Run on raw values, after flagOutOfRange. Missing or out-of-range slots
    /// break a low run.
    inline int flagCompression(CgmGrid& G, const CgmParams& p) {
        const int N = static_cast<int>(G.size());
        auto usable = [&](int i) {
            return i >= 0 && i < N && std::isfinite(G.raw[i]) &&
                !(G.flags[i] & out_of_range_code);
            };
        auto isLow = [&](int i) { return usable(i) && G.raw[i] < p.compLow; };
        const int retSlots = static_cast<int>(
            std::floor(p.compReturnWindowMinutes / G.stepMin + 1e-9));

        int n = 0;
        for (int s = 0; s < N; ++s) {
            if (!isLow(s)) continue;
            int e = s;
            while (isLow(e + 1)) ++e;
            const bool dropped = usable(s - 1) && G.raw[s - 1] >= p.compLow;
            const double dur = (e - s + 1) * G.stepMin;
            bool hit = false;
            if (dropped && dur > p.compMinMinutes) {
                for (int j = e + 1; j <= e + retSlots && j < N; ++j) {
                    if (!usable(j)) continue;
                    const double rate = (G.raw[j] - G.raw[e]) / ((j - e) * G.stepMin);
                    if (G.raw[j] >= p.compReturnLevel && rate >= p.compReturnRate) {
                        hit = true; break;
                    }
                }
            }
            if (hit) for (int k = s; k <= e; ++k) { G.flags[k] |= compression; ++n; }
            s = e;  // resume after the run
        }
        return n;
    }

    inline void applyExclusions(CgmGrid& G, const CgmParams& p) {
        for (std::size_t i = 0; i < G.size(); ++i) {
            const bool ex = ((G.flags[i] & out_of_range_code) && p.excludeOutOfRange) ||
                ((G.flags[i] & compression) && p.excludeCompression);
            G.value[i] = ex ? kNaN : G.raw[i];
            if (ex) G.flags[i] |= excluded;
        }
    }

    // ======================================================================
    //  Stage 4: gaps
    // ======================================================================
    inline int interpolateSingleGaps(CgmGrid& G) {
        int n = 0;
        for (std::size_t i = 1; i + 1 < G.size(); ++i) {
            const bool absent = !(G.flags[i] & observed_code);
            if (absent && !std::isfinite(G.value[i]) &&
                std::isfinite(G.value[i - 1]) && std::isfinite(G.value[i + 1])) {
                G.value[i] = 0.5 * (G.value[i - 1] + G.value[i + 1]);
                G.flags[i] |= interpolated_code;
                ++n;
            }
        }
        return n;
    }

    // ======================================================================
    //  Stage 5: scaling
    // ======================================================================
    /// Linear-interpolated percentile (Hyndman-Fan type 7, MATLAB/NumPy
    /// default), pct in [0, 100]. NaNs ignored.
    inline double percentile(std::vector<double> x, double pct) {
        x.erase(std::remove_if(x.begin(), x.end(),
            [](double v) { return !std::isfinite(v); }), x.end());
        if (x.empty()) return kNaN;
        std::sort(x.begin(), x.end());
        const double h = (x.size() - 1) * pct / 100.0;
        const auto lo = static_cast<std::size_t>(std::floor(h));
        const std::size_t hi = std::min(lo + 1, x.size() - 1);
        return x[lo] + (h - lo) * (x[hi] - x[lo]);
    }

    // ======================================================================
    //  Stage 6: bins
    // ======================================================================
    /// Score one bin or window from its slot values g (NaN = missing):
    ///   1. keep the present readings and run the spec's cgmBin on them
    ///      (it marks fewer than 2 readings -- under 30 min -- invalid);
    ///   2. more than 20% of the slots missing -> invalid;
    ///   3. fewer than 8 readings -> sd and cv blanked (sdValid = false);
    ///   4. normalized mean / sd / roc from cgmBin on the rescaled readings.
    inline void scoreBin(CgmBin& B, const std::vector<double>& g,
        const CgmParams& p, double p2, double p98) {
        std::vector<double> x;
        x.reserve(g.size());
        for (double v : g) if (std::isfinite(v)) x.push_back(v);
        B.n = static_cast<int>(x.size());
        B.raw = cgmBin(x);

        B.nMissing = B.expected - B.n;
        B.missingFrac = B.expected > 0 ? static_cast<double>(B.nMissing) / B.expected : 1.0;
        if (B.missingFrac > p.maxMissingFrac + 1e-12) B.raw.valid = false;

        B.sdValid = B.raw.valid && B.n >= p.minPointsSd;
        if (!B.sdValid) { B.raw.sd = kNaN; B.raw.cv = kNaN; }

        const double span = p98 - p2;
        if (B.raw.valid && std::isfinite(span) && span > 0.0) {
            std::vector<double> xn(x.size());
            for (std::size_t i = 0; i < x.size(); ++i) xn[i] = (x[i] - p2) / span;
            const CgmBinFeatures fn = cgmBin(xn);
            B.nmean = fn.mean;
            B.nroc = fn.roc;
            B.nsd = B.sdValid ? fn.sd : kNaN;
        }
    }

    inline std::vector<CgmBin> binGrid(const CgmGrid& G, const CgmParams& p,
        double p2, double p98) {
        std::vector<CgmBin> bins;
        const int N = static_cast<int>(G.size());
        const int binSlots = std::max(1,
            static_cast<int>(std::llround(p.binMinutes / G.stepMin)));

        for (int b = 0, start = 0; start < N; ++b, start += binSlots) {
            const int end = std::min(N, start + binSlots);
            CgmBin B;
            B.index = b;
            B.startMin = G.t0Min + start * G.stepMin;
            B.expected = end - start;

            std::vector<double> g(G.value.begin() + start, G.value.begin() + end);
            for (int i = start; i < end; ++i) {
                if (G.flags[i] & observed_code) ++B.nObserved;
                if (G.flags[i] & interpolated_code) ++B.nInterpolated;
                if (G.flags[i] & (out_of_range_code | compression)) ++B.nFlagged;
            }
            scoreBin(B, g, p, p2, p98);
            bins.push_back(B);
        }
        return bins;
    }

    // ======================================================================
    //  Whole pipeline
    // ======================================================================
    inline CgmPipelineResult runCgmPipeline(const std::vector<CgmReading>& readings,
        const CgmParams& p = CgmParams{}) {
        CgmPipelineResult R;
        R.nInput = static_cast<int>(readings.size());
        const auto kept = dropWarmup(readings, p, &R.nWarmupDropped);
        R.grid = buildGrid(kept, p, &R.nDuplicateSlot);
        R.nOutOfRange = flagOutOfRange(R.grid, p);
        R.nCompression = flagCompression(R.grid, p);
        applyExclusions(R.grid, p);
        R.nInterpolated = interpolateSingleGaps(R.grid);

        std::vector<double> ref;
        ref.reserve(R.grid.size());
        for (std::size_t i = 0; i < R.grid.size(); ++i)
            if ((R.grid.flags[i] & observed_code) && !(R.grid.flags[i] & excluded))
                ref.push_back(R.grid.value[i]);
        R.p2 = percentile(ref, p.pLow);
        R.p98 = percentile(ref, p.pHigh);

        R.bins = binGrid(R.grid, p, R.p2, R.p98);
        return R;
    }

    // ======================================================================
    //  Glucose features for ECG/PPG template bins (windows, not bins)
    // ======================================================================
    //
    // Template bins are typically 15 min. A Libre reads every 15 min, so a
    // template bin holds ONE reading: below the 30-min floor, and nowhere near
    // the 8 points SD/CV need. Re-binning glucose at 15 min therefore makes
    // every bin invalid.
    //
    // Instead each template bin gets two WINDOWS on the CGM grid:
    //   level        max(levelMinutes, template bin) = 30 min by default
    //                -> mean, roc, tir / tar / tbr, nmean, nroc
    //   variability  variabilityMinutes = 120 min by default
    //                -> sd, cv, nsd
    // Each window applies the same rules as a CGM bin (>= 2 points, <= 20%
    // missing, >= 8 points for SD/CV), with slots outside the CGM recording
    // counted as missing.
    //
    // WINDOWS OVERLAP. Neighbouring 15-min templates share 7 of the 8 points
    // in their 120-min window, so these features are a rolling series, not
    // independent samples. Account for that in any per-bin statistics.
    //
    // CLOCKS. recordStartMin must be on the CGM device clock, in the same
    // epoch-minute units as CgmReading::tMin. The ECG .bin header carries no
    // wall-clock start, so this comes from the caller; any offset between the
    // CGM and the ECG recorder clocks has to be corrected before the call.
    struct CgmTemplateWindow {
        enum class Align { Centered, Trailing };
        double levelMinutes = 30.0;
        double variabilityMinutes = 120.0;
        Align  align = Align::Centered;   // Trailing: window ends at bin end
    };

    struct CgmTemplateFeatures {
        int    bin = 0;                   // template bin index
        double startMin = 0.0;            // template bin start, CGM clock
        double glucoseAtMid = kNaN;       // linear between the 2 nearest slots
        CgmBin level{};                   // use level.raw.{mean,roc,tir,tar,tbr,valid}
        CgmBin variability{};             // use variability.raw.{sd,cv}, variability.sdValid
    };

    /// Features of the CGM grid slots whose times fall in [a, b).
    inline CgmBin featuresInWindow(const CgmGrid& G, const CgmParams& p,
        double a, double b, double p2, double p98) {
        CgmBin B;
        B.index = -1;
        B.startMin = a;
        const int N = static_cast<int>(G.size());
        const long long kLo = static_cast<long long>(std::ceil((a - G.t0Min) / G.stepMin - 1e-9));
        const long long kHi = static_cast<long long>(std::ceil((b - G.t0Min) / G.stepMin - 1e-9)) - 1;
        if (kHi < kLo) { B.raw = cgmBin({}); return B; }
        B.expected = static_cast<int>(kHi - kLo + 1);

        std::vector<double> g(static_cast<std::size_t>(B.expected), kNaN);
        for (long long k = kLo; k <= kHi; ++k) {
            if (k < 0 || k >= N) continue;   // outside the recording: missing
            const auto i = static_cast<std::size_t>(k);
            g[static_cast<std::size_t>(k - kLo)] = G.value[i];
            if (G.flags[i] & observed_code) ++B.nObserved;
            if (G.flags[i] & interpolated_code) ++B.nInterpolated;
            if (G.flags[i] & (out_of_range_code | compression)) ++B.nFlagged;
        }
        scoreBin(B, g, p, p2, p98);
        return B;
    }

    /// Template bin i covers [recordStartMin + i*L, recordStartMin + (i+1)*L).
    inline std::vector<CgmTemplateFeatures> cgmForTemplateBins(
        const CgmPipelineResult& R, double recordStartMin,
        double templateBinMinutes, int nBins,
        const CgmParams& p = CgmParams{},
        const CgmTemplateWindow& w = CgmTemplateWindow{}) {
        std::vector<CgmTemplateFeatures> out;
        if (nBins <= 0 || !(templateBinMinutes > 0.0)) return out;
        out.reserve(static_cast<std::size_t>(nBins));
        const CgmGrid& G = R.grid;
        const int N = static_cast<int>(G.size());
        const double Wl = std::max(w.levelMinutes, templateBinMinutes);
        const double Wv = std::max(w.variabilityMinutes, templateBinMinutes);

        auto window = [&](double s, double W, double& a, double& b) {
            if (w.align == CgmTemplateWindow::Align::Trailing) {
                b = s + templateBinMinutes; a = b - W;
            }
            else {
                const double mid = s + 0.5 * templateBinMinutes;
                a = mid - 0.5 * W; b = mid + 0.5 * W;
            }
            };

        for (int i = 0; i < nBins; ++i) {
            CgmTemplateFeatures T;
            T.bin = i;
            T.startMin = recordStartMin + i * templateBinMinutes;
            double a, b;
            window(T.startMin, Wl, a, b);
            T.level = featuresInWindow(G, p, a, b, R.p2, R.p98);
            window(T.startMin, Wv, a, b);
            T.variability = featuresInWindow(G, p, a, b, R.p2, R.p98);

            if (N > 0) {
                const double x = (T.startMin + 0.5 * templateBinMinutes - G.t0Min) / G.stepMin;
                const long long k0 = static_cast<long long>(std::floor(x));
                if (k0 >= 0 && k0 < N) {
                    const double v0 = G.value[static_cast<std::size_t>(k0)];
                    const double fr = x - static_cast<double>(k0);
                    if (fr < 1e-9) T.glucoseAtMid = v0;
                    else if (k0 + 1 < N) {
                        const double v1 = G.value[static_cast<std::size_t>(k0 + 1)];
                        T.glucoseAtMid = v0 + fr * (v1 - v0);   // NaN if either is
                    }
                }
            }
            out.push_back(T);
        }
        return out;
    }

    // ======================================================================
    //  Libre CSV reader
    // ======================================================================
    namespace detail {
        inline std::vector<std::string> splitCsvLine(const std::string& line) {
            std::vector<std::string> out;
            std::string cur;
            bool q = false;
            for (std::size_t i = 0; i < line.size(); ++i) {
                const char c = line[i];
                if (q) {
                    if (c == '"' && i + 1 < line.size() && line[i + 1] == '"') { cur += '"'; ++i; }
                    else if (c == '"') q = false;
                    else cur += c;
                }
                else if (c == '"') q = true;
                else if (c == ',') { out.push_back(cur); cur.clear(); }
                else if (c != '\r' && c != '\n') cur += c;
            }
            out.push_back(cur);
            return out;
        }

        inline std::string trim(const std::string& s) {
            std::size_t a = 0, b = s.size();
            while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
            while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
            return s.substr(a, b - a);
        }

        // Howard Hinnant's days_from_civil.
        inline long long daysFromCivil(long long y, unsigned m, unsigned d) {
            y -= m <= 2;
            const long long era = (y >= 0 ? y : y - 399) / 400;
            const unsigned yoe = static_cast<unsigned>(y - era * 400);
            const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
            const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
            return era * 146097 + static_cast<long long>(doe) - 719468;
        }
    }  // namespace detail

    /// "MM-DD-YYYY hh:mm AM/PM", "M/D/YYYY HH:MM[:SS]", or ISO
    /// "YYYY-MM-DD HH:MM[:SS]". Minutes since 1970-01-01, device-local clock.
    inline bool parseTimestampMinutes(const std::string& s, double& outMin) {
        std::vector<long long> nums;
        std::string ampm;
        long long cur = -1;
        for (char c : s) {
            if (std::isdigit(static_cast<unsigned char>(c))) {
                cur = (cur < 0 ? 0 : cur * 10) + (c - '0');
            }
            else {
                if (cur >= 0) { nums.push_back(cur); cur = -1; }
                if (std::isalpha(static_cast<unsigned char>(c)))
                    ampm += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            }
        }
        if (cur >= 0) nums.push_back(cur);
        if (nums.size() < 5) return false;

        long long y, mo, d;
        if (nums[0] > 31) { y = nums[0]; mo = nums[1]; d = nums[2]; }   // ISO
        else { mo = nums[0]; d = nums[1]; y = nums[2]; }                 // US
        long long h = nums[3], mi = nums[4];
        const double sec = nums.size() > 5 ? static_cast<double>(nums[5]) : 0.0;
        if (ampm == "PM" && h < 12) h += 12;
        if (ampm == "AM" && h == 12) h = 0;
        if (y < 100) y += 2000;
        if (mo < 1 || mo > 12 || d < 1 || d > 31 || h > 23 || mi > 59) return false;
        outMin = static_cast<double>(detail::daysFromCivil(y,
            static_cast<unsigned>(mo), static_cast<unsigned>(d))) * 1440.0 +
            h * 60.0 + mi + sec / 60.0;
        return true;
    }

    inline std::vector<CgmReading> parseLibreCsv(std::istream& is) {
        std::vector<CgmReading> out;
        std::string line;
        int cTime = -1, cType = -1, cHist = -1, cSerial = -1;
        bool header = false;
        while (std::getline(is, line)) {
            const auto f = detail::splitCsvLine(line);
            if (!header) {
                for (int i = 0; i < static_cast<int>(f.size()); ++i) {
                    const std::string t = detail::trim(f[i]);
                    if (t == "Device Timestamp") cTime = i;
                    else if (t == "Record Type") cType = i;
                    else if (t == "Historic Glucose mg/dL") cHist = i;
                    else if (t == "Serial Number") cSerial = i;
                }
                header = (cTime >= 0 && cHist >= 0);
                continue;
            }
            const int need = std::max(cTime, cHist);
            if (static_cast<int>(f.size()) <= need) continue;
            if (cType >= 0 && cType < static_cast<int>(f.size()) &&
                detail::trim(f[cType]) != "0") continue;   // historic only

            const std::string gs = detail::trim(f[cHist]);
            if (gs.empty()) continue;
            CgmReading r;
            if (!parseTimestampMinutes(f[cTime], r.tMin)) continue;
            if (cSerial >= 0 && cSerial < static_cast<int>(f.size()))
                r.serial = detail::trim(f[cSerial]);

            try {
                std::size_t used = 0;
                r.mgdl = std::stod(gs, &used);
                if (used != gs.size()) continue;
            }
            catch (...) { continue; }
            out.push_back(r);
        }
        return out;
    }

    inline std::vector<CgmReading> readLibreCsv(const std::string& path) {
        std::ifstream f(path);
        if (!f) return {};
        return parseLibreCsv(f);
    }

    // ======================================================================
    //  Output
    // ======================================================================
    namespace detail {
        // Howard Hinnant's civil_from_days, the inverse of daysFromCivil.
        inline void civilFromDays(long long z, long long& y, unsigned& m, unsigned& d) {
            z += 719468;
            const long long era = (z >= 0 ? z : z - 146096) / 146097;
            const unsigned doe = static_cast<unsigned>(z - era * 146097);
            const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
            const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
            const unsigned mp = (5 * doy + 2) / 153;
            d = doy - (153 * mp + 2) / 5 + 1;
            m = mp < 10 ? mp + 3 : mp - 9;
            y = static_cast<long long>(yoe) + era * 400 + (m <= 2);
        }
    }  // namespace detail

    /// "YYYY-MM-DD HH:MM", on the same clock the minutes came from (device
    /// local for Libre). Empty for a non-finite time.
    inline std::string formatEpochMinutes(double tMin) {
        if (!std::isfinite(tMin)) return {};
        const long long total = std::llround(tMin);
        long long days = total / 1440, rem = total % 1440;
        if (rem < 0) { rem += 1440; --days; }
        long long y; unsigned m, d;
        detail::civilFromDays(days, y, m, d);
        char buf[64];
        std::snprintf(buf, sizeof buf, "%04lld-%02u-%02u %02lld:%02lld",
            y, m, d, rem / 60, rem % 60);
        return buf;
    }

    /// One row per bin. start_time is readable (device clock); start_epoch_ms
    /// is the same instant in the unit the .bin files use.
    inline bool writeBinsCsv(const std::string& path, const CgmPipelineResult& R) {
        std::ofstream o(path);
        if (!o) return false;
        o << "bin,start_time,start_epoch_ms,expected,observed,interpolated,flagged,"
            "missing,missing_frac,valid,sd_valid,n,mean,sd,cv,roc,tir,tar,tbr,"
            "nmean,nsd,nroc\n";
        auto num = [&](double v) { if (std::isfinite(v)) o << v; };
        for (const auto& b : R.bins) {
            o << b.index << ',' << formatEpochMinutes(b.startMin) << ','
                << static_cast<long long>(std::llround(b.startMin * 60000.0)) << ','
                << b.expected << ',' << b.nObserved << ',' << b.nInterpolated << ','
                << b.nFlagged << ',' << b.nMissing << ',' << b.missingFrac << ','
                << (b.raw.valid ? 1 : 0) << ',' << (b.sdValid ? 1 : 0) << ','
                << b.n << ',';
            num(b.raw.mean); o << ','; num(b.raw.sd); o << ','; num(b.raw.cv); o << ',';
            num(b.raw.roc); o << ','; num(b.raw.tir); o << ','; num(b.raw.tar); o << ',';
            num(b.raw.tbr); o << ','; num(b.nmean); o << ','; num(b.nsd); o << ',';
            num(b.nroc); o << '\n';
        }
        return static_cast<bool>(o);
    }

    /// One row per 15-min grid slot, flags spelled out, so the processed trace
    /// can be checked against the source CSV by eye.
    inline bool writeGridCsv(const std::string& path, const CgmPipelineResult& R) {
        std::ofstream o(path);
        if (!o) return false;
        o << "slot,time,epoch_ms,raw_mgdl,value_mgdl,observed,interpolated,"
            "out_of_range,compression,excluded\n";
        auto num = [&](double v) { if (std::isfinite(v)) o << v; };
        const CgmGrid& G = R.grid;
        for (std::size_t k = 0; k < G.size(); ++k) {
            const double t = G.t0Min + static_cast<double>(k) * G.stepMin;
            const std::uint8_t f = G.flags[k];
            o << k << ',' << formatEpochMinutes(t) << ','
                << static_cast<long long>(std::llround(t * 60000.0)) << ',';
            num(G.raw[k]); o << ','; num(G.value[k]); o << ','
                << ((f & observed_code) ? 1 : 0) << ',' << ((f & interpolated_code) ? 1 : 0) << ','
                << ((f & out_of_range_code) ? 1 : 0) << ',' << ((f & compression) ? 1 : 0) << ','
                << ((f & excluded) ? 1 : 0) << '\n';
        }
        return static_cast<bool>(o);
    }

}  // namespace cgm_pipeline