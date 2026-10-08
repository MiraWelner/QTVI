#pragma once
/**
 * @file   sqi_ppg.hpp
 * @brief  Per-pulse PPG Signal Quality Index (Task A, Section 4.4).
 *
 *         The PPG needs its own index: motion corrupts it more readily than
 *         the ECG, and the ADC clips during strong pulses. Four components:
 *
 *           templateCorr  Pearson r of the pulse against the bin's PPG
 *                         template, both resampled to 200 samples (cubic) so
 *                         pulses of different duration are comparable.
 *           chiSq         chi-squared distance against the same template,
 *                         0 = identical shape, 1 = maximally different.
 *           clipping      1 - fraction of the pulse's samples within one
 *                         quantization step of either ADC rail. 1 = no
 *                         clipping. Rails are found empirically as the 0.1st
 *                         and 99.9th percentiles of the WHOLE record.
 *           motion        1 - motionFlag, where motionFlag is the share of the
 *                         pulse window the accelerometer flags as motion
 *                         (Task G): 1 = no motion.
 *
 *         composite = 0.40 corr + 0.25 (1 - chiSq) + 0.20 clipping + 0.15 motion,
 *         included when >= 0.70. An excluded pulse's reason -- the component
 *         that cost it the most -- is logged.
 *
 *         Wiring: writePpgSQICsv() is called from analysis_job::finalize()
 *         next to writeEcgSQICsv(). One CSV per input file, in
 *         cfg.training_log.
 *
 *         CHI-SQUARED, MADE BOUNDED. The composite uses (1 - chiSq), which
 *         only makes sense for a distance in [0, 1]. So both curves are first
 *         min-max scaled to [0, 1] and normalized to unit sum (p, q), and
 *             chiSq = 0.5 * sum (p - q)^2 / (p + q),
 *         the standard chi-squared distance between distributions, which is
 *         0 for identical shapes and at most 1.
 *
 *         MOTION CONVENTION. Here motionFlag is 1 for motion (the spec's
 *         `motion = 1 - motionFlag`), the opposite of the ECG SQI's 1 = clean.
 *         ppgMotionFlag() converts the accelerometer's per-interval answer.
 *         No accelerometer = no penalty (motionFlag 0).
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "accel/accel_pipeline.hpp"
#include "config_file_handling/config.hpp"
#include "peak_finding/peakfinding_structs.hpp"         // output_binfile_data
#include "prep_for_peakfinding/beat_times.hpp"          // Splice
#include "prep_for_peakfinding/data_bin_header.hpp"
#include "stats_utils.hpp"                              // ::pearson
#include "template_generation/template_structs.hpp"

namespace sqi_ppg {

    inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    inline constexpr int kPulseSamples = 200;

    // ---- Helpers the spec code calls --------------------------------------

    /// Natural cubic spline resample of v onto n evenly spaced points spanning
    /// the same interval. Interior NaNs are filled linearly first, edge NaNs
    /// take the nearest value. Empty if v has fewer than 2 finite samples.
    inline std::vector<double> resampleTo(const std::vector<double>& in, int n) {
        std::vector<double> v = in;
        const int m = static_cast<int>(v.size());
        int first = -1, last = -1;
        for (int i = 0; i < m; ++i) if (std::isfinite(v[i])) { if (first < 0) first = i; last = i; }
        if (first < 0 || first == last || n < 2) return {};
        for (int i = 0; i < first; ++i) v[i] = v[first];
        for (int i = last + 1; i < m; ++i) v[i] = v[last];
        for (int i = first + 1; i < last; ++i) {
            if (std::isfinite(v[i])) continue;
            int j = i; while (!std::isfinite(v[j])) ++j;
            for (int k = i; k < j; ++k)
                v[k] = v[i - 1] + (v[j] - v[i - 1]) * double(k - (i - 1)) / double(j - (i - 1));
            i = j;
        }
        // Natural spline second derivatives (unit knot spacing).
        std::vector<double> M(m, 0.0), c(m, 0.0), d(m, 0.0);
        for (int i = 1; i < m - 1; ++i) {
            const double rhs = 6.0 * (v[i + 1] - 2.0 * v[i] + v[i - 1]);
            const double den = 4.0 - c[i - 1];
            c[i] = 1.0 / den;
            d[i] = (rhs - d[i - 1]) / den;
        }
        for (int i = m - 2; i >= 1; --i) M[i] = d[i] - c[i] * M[i + 1];
        std::vector<double> out(n);
        for (int k = 0; k < n; ++k) {
            const double x = double(k) * (m - 1) / (n - 1);
            int i = std::min(m - 2, static_cast<int>(std::floor(x)));
            const double t = x - i, u = 1.0 - t;
            out[k] = u * v[i] + t * v[i + 1] + ((u * u * u - u) * M[i] + (t * t * t - t) * M[i + 1]) / 6.0;
        }
        return out;
    }

    /// Pearson r, 0 when fewer than 4 pairs or undefined (the ECG SQI's policy).
    inline double pearson(const std::vector<double>& a, const std::vector<double>& b) {
        const PearsonResult pr = ::pearson(a, b);
        return (pr.n_overlap < 4 || !pr.defined()) ? 0.0 : pr.r;
    }

    /// Chi-squared distance in [0, 1] between two curves' shapes (see header).
    inline double chiSquaredDistance(const std::vector<double>& a, const std::vector<double>& b) {
        const std::size_t n = std::min(a.size(), b.size());
        if (n == 0) return 1.0;
        auto unit = [n](const std::vector<double>& v) {
            double lo = std::numeric_limits<double>::infinity(), hi = -lo;
            for (std::size_t i = 0; i < n; ++i) if (std::isfinite(v[i])) { lo = std::min(lo, v[i]); hi = std::max(hi, v[i]); }
            std::vector<double> p(n, 0.0);
            double s = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                p[i] = (std::isfinite(v[i]) && hi > lo) ? (v[i] - lo) / (hi - lo) : 0.0;
                s += p[i];
            }
            if (s > 0.0) for (double& x : p) x /= s;
            return p;
            };
        const std::vector<double> p = unit(a), q = unit(b);
        double chi = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            if (p[i] + q[i] > 0.0) chi += (p[i] - q[i]) * (p[i] - q[i]) / (p[i] + q[i]);
        return std::min(1.0, 0.5 * chi);
    }

    // ---- ADC rails ----------------------------------------------------------
    struct AdcRails {
        double lo = kNaN, hi = kNaN;   // 0.1st / 99.9th percentile of the record
        double step = 0.0;             // quantization step
        bool known() const { return std::isfinite(lo) && std::isfinite(hi); }
        // A sample at or beyond these counts as clipped.
        double loThreshold() const { return known() ? lo + step : -std::numeric_limits<double>::infinity(); }
        double hiThreshold() const { return known() ? hi - step : std::numeric_limits<double>::infinity(); }
    };

    /// Rails from every sample of the record. The quantization step is the
    /// smallest gap between two distinct sample values.
    inline AdcRails detectAdcRails(std::vector<double> v) {
        v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
        AdcRails r;
        if (v.size() < 2) return r;
        std::sort(v.begin(), v.end());
        auto pct = [&](double p) {
            const double h = (v.size() - 1) * p / 100.0;
            const auto i = static_cast<std::size_t>(std::floor(h));
            const std::size_t j = std::min(i + 1, v.size() - 1);
            return v[i] + (h - i) * (v[j] - v[i]);
            };
        r.lo = pct(0.1);
        r.hi = pct(99.9);
        double step = std::numeric_limits<double>::infinity();
        for (std::size_t i = 1; i < v.size(); ++i)
            if (v[i] > v[i - 1]) step = std::min(step, v[i] - v[i - 1]);
        r.step = std::isfinite(step) ? step : 0.0;
        return r;
    }

    /// The record's PPG (channel slot 4) from the original data .bin: the
    /// native samples when stored, else the upsampled ones. Native samples are
    /// what the ADC produced, so they carry its true quantization step.
    inline std::vector<double> readPpgRecordFromBin(const std::string& path) {
        constexpr int kSlotPpg = 4;
        const data_bin_header::Header h = data_bin_header::read(path);
        if (!h.ok) return {};
        uint64_t at = data_bin_header::kHeaderBytes, upOff = 0;
        for (int i = 0; i < kSlotPpg; ++i) at += 8ull * h.sizes_up[i] + 16ull * h.sizes_raw[i];
        upOff = at;
        const uint64_t rawOff = upOff + 8ull * h.sizes_up[kSlotPpg];
        std::ifstream f(path, std::ios::binary);
        if (!f) return {};
        std::vector<double> out;
        const uint64_t nRaw = h.sizes_raw[kSlotPpg];
        if (nRaw > 1 && rawOff + 16ull * nRaw <= h.file_size) {
            std::vector<double> pairs(2 * nRaw);
            f.seekg(static_cast<std::streamoff>(rawOff));
            f.read(reinterpret_cast<char*>(pairs.data()), static_cast<std::streamsize>(16ull * nRaw));
            if (f) { out.resize(nRaw); for (uint64_t i = 0; i < nRaw; ++i) out[i] = pairs[2 * i + 1]; return out; }
            f.clear();
        }
        const uint64_t nUp = h.sizes_up[kSlotPpg];
        if (nUp > 1 && upOff + 8ull * nUp <= h.file_size) {
            out.resize(nUp);
            f.seekg(static_cast<std::streamoff>(upOff));
            f.read(reinterpret_cast<char*>(out.data()), static_cast<std::streamsize>(8ull * nUp));
            if (!f) out.clear();
        }
        return out;
    }

    /// Accelerometer interval answer (1 clean, 0 motion, -1 unavailable) to
    /// this file's motionFlag (1 = motion, 0 = none; unavailable -> 0).
    inline double ppgMotionFlag(int accelFlag) { return accelFlag == 0 ? 1.0 : 0.0; }

    // ---- Spec struct and function, as given -------------------------------
    struct PulseSQI {
        double templateCorr, chiSq, clipping, motion, composite; bool include;
    };
    inline PulseSQI computePpgSQI(const std::vector<double>& pulse,
        const std::vector<double>& medianTemplate200,
        double adcLo, double adcHi, double motionFlag) {
        std::vector<double> r = resampleTo(pulse, 200); // cubic, Phase 1 Task N
        PulseSQI q{};
        q.templateCorr = pearson(r, medianTemplate200);
        q.chiSq = chiSquaredDistance(r, medianTemplate200);
        int clipped = 0;
        for (double v : pulse) if (v <= adcLo || v >= adcHi) ++clipped;
        q.clipping = 1.0 - (double)clipped / pulse.size(); // 1.0 means no clipping
        q.motion = 1.0 - motionFlag;
        q.composite = 0.40 * q.templateCorr + 0.25 * (1.0 - q.chiSq)
            + 0.20 * q.clipping + 0.15 * q.motion;
        q.include = q.composite >= 0.70;
        return q;
    }

    /// Why a pulse was excluded: the component that cost the composite the
    /// most. "" for an included pulse.
    inline const char* exclusionReason(const PulseSQI& q) {
        if (q.include) return "";
        const double loss[4] = { 0.40 * (1.0 - q.templateCorr), 0.25 * q.chiSq,
                                 0.20 * (1.0 - q.clipping), 0.15 * (1.0 - q.motion) };
        static const char* const names[4] = { "template", "chi_squared", "clipping", "motion" };
        int worst = 0;
        for (int i = 1; i < 4; ++i) if (loss[i] > loss[worst]) worst = i;
        return names[worst];
    }

    // ---- File-level driver ------------------------------------------------
    //
    // Scores every kept pulse in every bin against that bin's PPG template and
    // writes one row per pulse to <cfg.training_log>/<stem>_ppg_quality.csv.
    //
    // A pulse's window is its R-pair: [R(s) + lag, R(s+1) + lag], s being the
    // pulse row's R-pair ordinal (ppg_slice_of_row) -- the same window the
    // pulse slicer searches. motionFlag is the share of that window the
    // accelerometer flags; with no accelerometer it is 0 (no penalty).
    inline void writePpgSQICsv(const config_entry& cfg,
        const std::string& stem,
        const template_structs::TemplateFile& tmpl,
        const template_structs::BeatsFile& beats,
        const std::vector<std::vector<uint32_t>>& ppgSliceOfRow,   // [bin][row]
        const std::vector<output_binfile_data>& peakResults,
        double ecgRecRateHz,                 // rate R positions are in
        double ppgLagMs,
        const AdcRails& rails,
        const accel_pipeline::AccelResult* accel) {
        const auto it = beats.per_channel_beats.find("PPG");
        if (it == beats.per_channel_beats.end()) return;   // no PPG in this record

        const std::string outPath = cfg.logs + "/" + stem + "_ppg_quality.csv";
        std::ofstream f(outPath);
        if (!f.is_open()) {
            std::cerr << "  WARNING: could not open " << outPath << " for PPG SQI output\n";
            return;
        }
        f << "bin,beat,template_corr,chi_sq,clipping,motion,composite,is_included,exclusion_reason\n";
        if (!rails.known())
            std::cerr << "  [ppg sqi] no ADC rails for " << stem << "; clipping scored as 1\n";

        const double lagS = ppgLagMs / 1000.0;
        for (std::size_t bin = 0; bin < tmpl.bins.size() && bin < it->second.size(); ++bin) {
            const auto& bt = tmpl.bins[bin];
            if (bt.ppgTemplate.empty()) continue;
            const std::vector<double> tmpl200 = resampleTo(bt.ppgTemplate, kPulseSamples);
            if (tmpl200.empty()) continue;

            // Recording seconds of R-pair s, NaN if unknown.
            const output_binfile_data* pr = bin < peakResults.size() ? &peakResults[bin] : nullptr;
            const beat_times::Splice sp(pr ? pr->ecg_bin_indexs : std::vector<std::pair<uint64_t, uint64_t>>{},
                std::numeric_limits<uint64_t>::max());
            auto rSec = [&](std::size_t s) {
                if (!pr || !(ecgRecRateHz > 0.0) || s >= pr->ch1.raw.size()) return kNaN;
                uint64_t orig = 0;
                return sp.original(static_cast<uint64_t>(pr->ch1.raw[s]), orig)
                    ? static_cast<double>(orig) / ecgRecRateHz : kNaN;
                };
            auto motionOf = [&](std::size_t row) {
                if (!accel || bin >= ppgSliceOfRow.size() || row >= ppgSliceOfRow[bin].size()) return 0.0;
                const std::size_t s = ppgSliceOfRow[bin][row];
                const double t0 = rSec(s) + lagS;
                double t1 = rSec(s + 1) + lagS;
                if (!std::isfinite(t0)) return 0.0;
                if (!std::isfinite(t1) || t1 <= t0) t1 = t0 + 1.0;
                const double frac = accel->motionFraction(t0, t1);   // Task G
                return std::isfinite(frac) ? frac : 0.0;           // no coverage: no penalty
                };

            const auto& pulses = it->second[bin];   // [beat][sample]
            for (std::size_t bi = 0; bi < pulses.size(); ++bi) {
                if (pulses[bi].empty()) continue;
                const PulseSQI q = computePpgSQI(pulses[bi], tmpl200,
                    rails.loThreshold(), rails.hiThreshold(), motionOf(bi));
                f << bin << ',' << bi << ',' << q.templateCorr << ',' << q.chiSq << ','
                    << q.clipping << ',' << q.motion << ',' << q.composite << ','
                    << (q.include ? "INCLUDE" : "EXCLUDE") << ',' << exclusionReason(q) << '\n';
            }
        }
    }

}  // namespace sqi_ppg