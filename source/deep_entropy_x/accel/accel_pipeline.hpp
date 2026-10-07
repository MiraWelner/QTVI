#pragma once
//
// accel_pipeline.hpp  (Task G, Section 8)
//
// Vector magnitude, four-state activity classification, percentile
// normalization, and a motion flag that feeds ECG and PPG SQI. Header-only.
//
// STAGES (runAccelPipeline)
//   1. VM.        VM(t) = sqrt(X^2 + Y^2 + Z^2) per sample, in g. NaN if any
//                 axis is missing.
//   2. MOTION.    A sample is flagged when VM > 1.05 g.
//   3. EPOCHS.    Fixed epochs of epochSec (30 s, the sleep-epoch length) from
//                 the start of the record. Per epoch: meanVM, sdVM (n - 1),
//                 zOverVM = |mean Z| / meanVM and yOverVM = |mean Y| / meanVM
//                 (the share of gravity on each axis, i.e. posture), and the
//                 fraction of samples flagged as motion.
//   4. CLASSIFY.  The spec's classify() on those four numbers, in raw g.
//   5. NORMALIZE. Each epoch's meanVM and sdVM are also rescaled by the
//                 recording's own 2nd and 98th percentiles of that quantity
//                 (valid epochs only): n = (v - p2) / (p98 - p2), not clipped.
//                 Classification and the motion flag use raw g, never these.
//
// EPOCH VALIDITY. An epoch is valid when at least 80% of its samples are
// present. Invalid epochs still carry their numbers but should not be used.
//
// FEEDING SQI. sqiMotionFlag(t0, t1) answers for any time interval (e.g. one
// beat) in the convention computeEcgSQI already takes as motionFlag:
//   1 clean, 0 motion (some sample in the interval has VM > 1.05 g),
//  -1 unavailable (no accelerometer samples in the interval).
// computeEcgSQI multiplies its composite by this (with -1 counting as 1), so a
// motion beat scores 0 and is EXCLUDEd. sqiMotionTerm() is that same rule as a
// function, for the PPG SQI to use the same way.
//
// DOWNSTREAM. activityAt(t) gives the state at a recording time, for dynamic
// template selection (9.5), the local baseline (10) and sleep-onset alignment
// (10.3).
//
// UNITS. The thresholds are in g. If the recorder stores milli-g or m/s^2, set
// AccelParams::unitsPerG (1000, or 9.80665) so samples are converted to g.
//
// TIME. Seconds from the start of the record: sample i is at i / fs, the same
// axis beat_times.hpp and record_sleep.hpp use.
//

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "prep_for_peakfinding/data_bin_header.hpp"

namespace accel_pipeline {

    inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // ---- Spec classification, as given ----------------------------------
    enum class Activity { SLEEP_SUPINE, SEATED_REST, LIGHT, MODERATE_VIGOROUS };
    inline Activity classify(double meanVM, double sdVM, double zOverVM, double yOverVM) {
        if (meanVM > 1.10 || sdVM > 0.10) return Activity::MODERATE_VIGOROUS;
        if (meanVM >= 1.02 || sdVM >= 0.02) return Activity::LIGHT;
        if (zOverVM > 0.9) return Activity::SLEEP_SUPINE; // posture from gravity axis
        if (yOverVM > 0.7) return Activity::SEATED_REST;
        return Activity::SEATED_REST;
    }
    // VM(t) = sqrt(X^2+Y^2+Z^2); motion artifact flag when VM > 1.05 g

    inline const char* activityName(Activity a) {
        switch (a) {
        case Activity::SLEEP_SUPINE:      return "SLEEP_SUPINE";
        case Activity::SEATED_REST:       return "SEATED_REST";
        case Activity::LIGHT:             return "LIGHT";
        case Activity::MODERATE_VIGOROUS: return "MODERATE_VIGOROUS";
        }
        return "?";
    }

    // ---- Parameters -------------------------------------------------------
    struct AccelParams {
        double epochSec = 30.0;          // classification epoch
        double motionThresholdG = 1.05;  // motion flag when VM > this
        double unitsPerG = 1.0;          // stored units per g (1000 for mg)
        double minEpochCoverage = 0.80;  // valid epoch: >= 80% samples present
        double pLow = 2.0;               // normalization percentiles
        double pHigh = 98.0;
    };

    // ---- Input ------------------------------------------------------------
    struct AccelData {
        double fs = 0.0;                 // samples per second
        std::vector<double> x, y, z;     // stored units
        bool present() const {
            return fs > 0.0 && !x.empty() && x.size() == y.size() && x.size() == z.size();
        }
    };

    // Accelerometer slots in the data .bin (gui_handler.h ChannelIdx).
    inline constexpr int kSlotAccelX = 5, kSlotAccelY = 6, kSlotAccelZ = 7;

    /// Upsampled X/Y/Z and their rate from a file_to_bin v1 data .bin. Empty
    /// (present() == false) when the file is unreadable or has no
    /// accelerometer -- a channel of size 0, or the single -1 sentinel.
    inline AccelData readAccelFromBin(const std::string& path) {
        AccelData a;
        const data_bin_header::Header h = data_bin_header::read(path);
        if (!h.ok) return a;

        // Per channel: sizes_up doubles, then 2 * sizes_raw doubles.
        uint64_t off[data_bin_header::kNumChannels];
        uint64_t at = data_bin_header::kHeaderBytes;
        for (int i = 0; i < data_bin_header::kNumChannels; ++i) {
            off[i] = at;
            at += 8ull * h.sizes_up[i] + 16ull * h.sizes_raw[i];
        }
        std::ifstream f(path, std::ios::binary);
        if (!f) return a;
        auto load = [&](int slot, std::vector<double>& v) {
            const uint64_t n = h.sizes_up[slot];
            if (n == 0 || off[slot] + 8ull * n > h.file_size) return false;
            v.resize(n);
            f.seekg(static_cast<std::streamoff>(off[slot]));
            f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(8ull * n));
            if (!f) return false;
            return !(n == 1 && v[0] == -1.0);
            };
        if (!load(kSlotAccelX, a.x) || !load(kSlotAccelY, a.y) || !load(kSlotAccelZ, a.z)) return {};
        const std::size_t n = std::min({ a.x.size(), a.y.size(), a.z.size() });
        a.x.resize(n); a.y.resize(n); a.z.resize(n);
        a.fs = static_cast<double>(h.up_rates[kSlotAccelX]);
        return a.present() ? a : AccelData{};
    }

    // ---- Stage 1: vector magnitude ----------------------------------------
    inline std::vector<double> vectorMagnitude(const AccelData& a, const AccelParams& p = AccelParams{}) {
        std::vector<double> vm(a.x.size(), kNaN);
        const double s = (p.unitsPerG > 0.0) ? 1.0 / p.unitsPerG : 1.0;
        for (std::size_t i = 0; i < vm.size(); ++i) {
            const double x = a.x[i] * s, y = a.y[i] * s, z = a.z[i] * s;
            if (std::isfinite(x) && std::isfinite(y) && std::isfinite(z))
                vm[i] = std::sqrt(x * x + y * y + z * z);
        }
        return vm;
    }

    // ---- Results ----------------------------------------------------------
    struct Epoch {
        int    index = 0;
        double startSec = 0.0;
        int    expected = 0;          // samples the epoch spans
        int    n = 0;                 // samples present
        bool   valid = false;         // n >= minEpochCoverage * expected
        double meanVM = kNaN, sdVM = kNaN;          // g
        double zOverVM = kNaN, yOverVM = kNaN;      // posture
        Activity activity = Activity::SEATED_REST;
        double motionFrac = 0.0;      // share of present samples with VM > 1.05
        bool   motion = false;        // any sample with VM > 1.05
        double nMeanVM = kNaN, nSdVM = kNaN;        // percentile-normalized
    };

    struct AccelResult {
        double fs = 0.0;
        std::vector<double>  vm;         // per sample, g (NaN = missing)
        std::vector<uint8_t> motion;     // per sample, 1 = VM > 1.05 g
        std::vector<Epoch>   epochs;
        double meanP2 = kNaN, meanP98 = kNaN;   // percentiles of epoch meanVM
        double sdP2 = kNaN, sdP98 = kNaN;       // percentiles of epoch sdVM

        /// 1 clean, 0 motion, -1 unavailable -- for samples with t0 <= t <= t1.
        int sqiMotionFlag(double t0, double t1) const {
            if (!(fs > 0.0) || vm.empty() || !(t1 >= t0)) return -1;
            const long long n = static_cast<long long>(vm.size());
            const long long i0 = std::max(0LL, static_cast<long long>(std::ceil(t0 * fs - 1e-9)));
            const long long i1 = std::min(n - 1, static_cast<long long>(std::floor(t1 * fs + 1e-9)));
            bool any = false;
            for (long long i = i0; i <= i1; ++i) {
                if (!std::isfinite(vm[static_cast<std::size_t>(i)])) continue;
                any = true;
                if (motion[static_cast<std::size_t>(i)]) return 0;
            }
            return any ? 1 : -1;
        }

        /// Share of the present samples with t0 <= t <= t1 that are flagged as
        /// motion, 0..1; NaN when there are none (no accelerometer coverage).
        double motionFraction(double t0, double t1) const {
            if (!(fs > 0.0) || vm.empty() || !(t1 >= t0)) return kNaN;
            const long long n = static_cast<long long>(vm.size());
            const long long i0 = std::max(0LL, static_cast<long long>(std::ceil(t0 * fs - 1e-9)));
            const long long i1 = std::min(n - 1, static_cast<long long>(std::floor(t1 * fs + 1e-9)));
            int present = 0, moving = 0;
            for (long long i = i0; i <= i1; ++i) {
                if (!std::isfinite(vm[static_cast<std::size_t>(i)])) continue;
                ++present;
                moving += motion[static_cast<std::size_t>(i)];
            }
            return present ? static_cast<double>(moving) / present : kNaN;
        }

        /// The epoch containing time t, or nullptr.
        const Epoch* epochAt(double t) const {
            if (epochs.empty() || !(t >= 0.0)) return nullptr;
            const double len = epochs.size() > 1 ? epochs[1].startSec - epochs[0].startSec
                : static_cast<double>(epochs[0].expected) / fs;
            const auto e = static_cast<std::size_t>(std::floor(t / len));
            return e < epochs.size() ? &epochs[e] : nullptr;
        }

        /// Activity at time t; false when t has no valid epoch.
        bool activityAt(double t, Activity& out) const {
            const Epoch* e = epochAt(t);
            if (!e || !e->valid) return false;
            out = e->activity;
            return true;
        }
    };

    /// The SQI composite's motion factor, the rule computeEcgSQI applies:
    /// -1 (unavailable) leaves the score alone, 0 (motion) zeroes it.
    inline double sqiMotionTerm(int motionFlag) {
        return (motionFlag < 0) ? 1.0 : static_cast<double>(motionFlag);
    }

    /// Linear-interpolated percentile (type 7), NaNs ignored.
    inline double percentile(std::vector<double> v, double pct) {
        v.erase(std::remove_if(v.begin(), v.end(), [](double d) { return !std::isfinite(d); }), v.end());
        if (v.empty()) return kNaN;
        std::sort(v.begin(), v.end());
        const double h = (v.size() - 1) * pct / 100.0;
        const auto lo = static_cast<std::size_t>(std::floor(h));
        const std::size_t hi = std::min(lo + 1, v.size() - 1);
        return v[lo] + (h - lo) * (v[hi] - v[lo]);
    }

    // ---- Whole pipeline -----------------------------------------------------
    inline AccelResult runAccelPipeline(const AccelData& a, const AccelParams& p = AccelParams{}) {
        AccelResult R;
        if (!a.present()) return R;
        R.fs = a.fs;

        // Stages 1-2
        R.vm = vectorMagnitude(a, p);
        R.motion.assign(R.vm.size(), 0);
        for (std::size_t i = 0; i < R.vm.size(); ++i)
            if (std::isfinite(R.vm[i]) && R.vm[i] > p.motionThresholdG) R.motion[i] = 1;

        // Stages 3-4
        const double s = (p.unitsPerG > 0.0) ? 1.0 / p.unitsPerG : 1.0;
        const std::size_t N = R.vm.size();
        const auto len = static_cast<std::size_t>(std::max(1.0, std::round(p.epochSec * a.fs)));
        for (std::size_t start = 0, e = 0; start < N; start += len, ++e) {
            const std::size_t end = std::min(N, start + len);
            Epoch E;
            E.index = static_cast<int>(e);
            E.startSec = static_cast<double>(start) / a.fs;
            E.expected = static_cast<int>(len);

            double sum = 0.0, sy = 0.0, sz = 0.0;
            int nMotion = 0;
            for (std::size_t i = start; i < end; ++i) {
                if (!std::isfinite(R.vm[i])) continue;
                sum += R.vm[i]; sy += a.y[i] * s; sz += a.z[i] * s;
                nMotion += R.motion[i];
                ++E.n;
            }
            if (E.n > 0) {
                E.meanVM = sum / E.n;
                double ss = 0.0;
                for (std::size_t i = start; i < end; ++i)
                    if (std::isfinite(R.vm[i])) ss += (R.vm[i] - E.meanVM) * (R.vm[i] - E.meanVM);
                E.sdVM = E.n > 1 ? std::sqrt(ss / (E.n - 1)) : 0.0;
                E.zOverVM = E.meanVM > 0.0 ? std::fabs(sz / E.n) / E.meanVM : kNaN;
                E.yOverVM = E.meanVM > 0.0 ? std::fabs(sy / E.n) / E.meanVM : kNaN;
                E.activity = classify(E.meanVM, E.sdVM, E.zOverVM, E.yOverVM);
                E.motionFrac = static_cast<double>(nMotion) / E.n;
                E.motion = nMotion > 0;
            }
            E.valid = E.n >= p.minEpochCoverage * static_cast<double>(E.expected);
            R.epochs.push_back(E);
        }

        // Stage 5
        std::vector<double> means, sds;
        for (const Epoch& E : R.epochs)
            if (E.valid) { means.push_back(E.meanVM); sds.push_back(E.sdVM); }
        R.meanP2 = percentile(means, p.pLow);  R.meanP98 = percentile(means, p.pHigh);
        R.sdP2 = percentile(sds, p.pLow);      R.sdP98 = percentile(sds, p.pHigh);
        const double mSpan = R.meanP98 - R.meanP2, sSpan = R.sdP98 - R.sdP2;
        for (Epoch& E : R.epochs) {
            if (!E.valid) continue;
            if (mSpan > 0.0) E.nMeanVM = (E.meanVM - R.meanP2) / mSpan;
            if (sSpan > 0.0) E.nSdVM = (E.sdVM - R.sdP2) / sSpan;
        }
        return R;
    }

}  // namespace accel_pipeline#pragma once
