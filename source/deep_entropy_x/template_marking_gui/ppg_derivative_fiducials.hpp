#pragma once
/*
ppg_derivative_fiducials.hpp

Given a PPG pulse, get the first and second derivatives to calculate the derivative-based fiducial points from pyPPG: 
a Python toolbox for comprehensive photoplethysmography signal analysis by Márton Á Goda, Peter H Charlton and Joachim A Behar

VPG (velocity of ppg - first derivative) markers
- u: The highest amplitude between the pulse onset and systolic peak on PPG′
- v: The lowest amplitude between the u-point and diastolic peak on PPG′
- w: The first local maximum or inflection point after the dicrotic notch on PPG

APG (acceleration of ppg - second derivative) markers
- a: The highest amplitude between pulse onset and systolic peak on PPG″
- b: The first local minimum after the a-point on PPG″
- c: The local maximum with the highest amplitude between the b-point and e-point, or if no local maximum is present, 
     then the inflection point on PPG″
- d: The local minimum with the lowest amplitude between the c-point and e-point, or if no local minimum is present, 
     then the inflection point on PPG″
- e: The local maximum with the highest amplitude after the b-point and before the diastolic peak on PPG″

- f: The first local minimum after the e-point on PPG″


JPG (jerk of ppg - second derivative) markers
- p1 The first local maximum after the b-point on PPG‴
- p1 The first local maximum after the b-point on PPG‴

@author Mira Welner
@date 08-21-2026
*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace ppg_deriv {

    /* from the paper: "10 ms MAF for the PPG derivatives: to eliminate the high-frequency content in the PPG derivatives, 
    a 10 ms standard flat (boxcar or top-hat) MAF with 45 Hz cut-off frequency was applied.*/
    inline constexpr double ten_ms_maf = 10.0;

    inline std::vector<double> compute_vpg(const std::vector<double>& ppg, double fs,
        double smoothMs = ten_ms_maf) {
        const int N = static_cast<int>(ppg.size());
        std::vector<double> d(N, std::numeric_limits<double>::quiet_NaN());
        if (fs <= 0.0) return d;
        // At least a 5-point operator, else this degenerates back to a
        // 3-point central difference and the smoothing is notional.
        const int m = std::max(2, static_cast<int>(std::lround(0.5 * (smoothMs / 1000.0) * fs)));
        if (N < 2 * m + 1) return d;
        double norm = 0.0;
        for (int k = -m; k <= m; ++k) norm += static_cast<double>(k) * k;
        for (int i = m; i < N - m; ++i) {
            double s = 0.0;
            bool ok = true;
            for (int k = -m; k <= m && ok; ++k) {
                const double x = ppg[i + k];
                if (std::isnan(x)) ok = false; else s += k * x;
            }
            // NaN propagates rather than being interpolated across. The
            // template tails carry NaN by construction, and inventing
            // samples there would put phantom extrema in w's range.
            if (ok) d[i] = s / norm;
        }
        return d;
    }

    // -1 = not detected, matching every other landmark field in the tree.
    // No midpoint placeholders: a missing w is a real morphological fact
    // in a damped pulse, and a fabricated one would feed Aw/Au a silently
    // wrong number.
    struct VpgFiducials {
        int u = -1;
        int v = -1;
        int w = -1;
    };

    namespace detail {
        inline int argmax_in(const std::vector<double>& s, int lo, int hi) {
            lo = std::max(0, lo); hi = std::min(hi, static_cast<int>(s.size()) - 1);
            int best = -1; double bv = -std::numeric_limits<double>::infinity();
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(s[i]) && s[i] > bv) { bv = s[i]; best = i; }
            return best;
        }
        inline int argmin_in(const std::vector<double>& s, int lo, int hi) {
            lo = std::max(0, lo); hi = std::min(hi, static_cast<int>(s.size()) - 1);
            int best = -1; double bv = std::numeric_limits<double>::infinity();
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(s[i]) && s[i] < bv) { bv = s[i]; best = i; }
            return best;
        }
        inline int firstLocalMax(const std::vector<double>& s, int lo, int hi) {
            lo = std::max(1, lo); hi = std::min(hi, static_cast<int>(s.size()) - 2);
            for (int i = lo; i <= hi; ++i) {
                if (std::isnan(s[i - 1]) || std::isnan(s[i]) || std::isnan(s[i + 1])) continue;
                if (s[i] >= s[i - 1] && s[i] > s[i + 1]) return i;
            }
            return -1;
        }
    } // namespace detail

    // on / sp / off are FeatureMarks::detect_ppg_fiducials' already-resolved
    // g.onset, g.peak, g.end. Wc is that function's visible-window clamp:
    // nothing is placed past Wc-1, same contract as its cl() lambda.
    //
    // dp is the diastolic peak, table 3's upper bound for v. Pass g.peak2
    // when g.peak2_found, and -1 otherwise -- NOT the unconditional
    // g.peak2, because when peak2_found is false that field holds a
    // midpoint placeholder, and a bound derived from a fabricated landmark
    // is worse than the coarser duration cap this falls back to.
    //
    // Reads nothing else. In particular it never touches g.dicrotic: table
    // 3 defines w off the notch, but table 4's correction restates it in
    // terms of v (and, once APG lands, e and f), and that is the form
    // implemented below.
    inline VpgFiducials detect_vpg(const std::vector<double>& vpg,
        int on, int sp, int off, int Wc, int dp = -1)
    {
        VpgFiducials r;
        if (vpg.empty() || on < 0 || sp <= on) return r;

        const int hiW = std::min(Wc, static_cast<int>(vpg.size())) - 1;
        if (hiW <= sp) return r;
        auto cl = [&](int x) { return std::clamp(x, 0, hiW); };

        const int pulseEnd = (off > sp) ? std::min(off, hiW) : hiW;
        const double dur = static_cast<double>(pulseEnd - on);
        if (dur < 4.0) return r;
        const int cap60 = cl(on + static_cast<int>(std::lround(0.60 * dur)));
        const int cap80 = cl(on + static_cast<int>(std::lround(0.80 * dur)));

        // ---- u: max VPG on [on, sp] ------------------------------------
        r.u = detail::argmax_in(vpg, on, sp);
        if (r.u < 0) return r;             // no finite VPG in systole

        // ---- v: min VPG on (u, vHi] ------------------------------------
        // vHi is the diastolic peak when it was genuinely detected (table
        // 3's bound), else 60% of pulse duration -- the same fraction
        // table 4 uses to cap e, and a safe proxy since the systolic
        // downslope always precedes it.
        //
        // argmin, not lowest-local-min: on a heavily damped pulse the
        // downslope can be monotone through the whole window, and the
        // steepest point is then an endpoint, which is not a local
        // minimum. Table 3 says "lowest amplitude", so argmin is the
        // literal reading anyway.
        const int vHi = (dp > r.u) ? std::min(dp, hiW) : cap60;
        r.v = detail::argmin_in(vpg, r.u + 1, vHi);
        if (r.v < 0) return r;

        // ---- w: first local max after v, capped at 80% -----------------
        r.w = detail::firstLocalMax(vpg, r.v + 1, cap80);

        return r;
    }


    // =====================================================================
    // APG (PPG'') fiducials: a, b, c, d, e, f
    // =====================================================================
    //
    // DEPENDENCY ORDER IS NOT OPTIONAL. Table 3 defines these in terms of
    // each other -- b bounds c and e, e bounds c and d, c bounds d, e
    // bounds f -- so detect_apg runs them in the only order that resolves:
    //
    //     a -> b -> e -> c -> d -> f
    //
    // and NOT the alphabetical order the table lists them in.

    // -1 = not detected, matching every other landmark field in the tree.
    //
    // NO INFLECTION FALLBACK, deliberately. Table 3 says c and d fall back to
    // "the inflection point" when no local extremum exists between their
    // bounds. Implementing that means locating an extremum of the JPG and
    // reporting it as an APG landmark -- a materially weaker measurement that
    // a consumer cannot distinguish from a real extremum unless it travels
    // with a flag. Rather than carry the flag, c and d are reported only when
    // a genuine APG extremum exists; -1 otherwise, which is self-describing.
    //
    // Consequence: on Dawber class 3/4 (shoulder) waves pyPPG will report a
    // c/d where this returns -1. Those are the pulses where the c/d position
    // is least trustworthy anyway.
    struct ApgFiducials {
        int a = -1, b = -1, c = -1, d = -1, e = -1, f = -1;
    };

    // Chained SG differences. Each stage smooths, matching pyPPG's
    // difference-then-moving-average per derivative order.
    inline std::vector<double> compute_apg(const std::vector<double>& vpg, double fs,
        double smoothMs = ten_ms_maf) {
        return compute_vpg(vpg, fs, smoothMs);
    }
    // Unused until p1/p2 land (they are the only points defined on the JPG).
    // Kept here because it is the same operator and the alternative is
    // rediscovering that fact later.
    inline std::vector<double> compute_jpg(const std::vector<double>& apg, double fs,
        double smoothMs = ten_ms_maf) {
        return compute_vpg(apg, fs, smoothMs);
    }

    namespace detail {

        inline bool finite3(const std::vector<double>& s, int i) {
            return !std::isnan(s[i - 1]) && !std::isnan(s[i]) && !std::isnan(s[i + 1]);
        }
        inline bool isLocalMax(const std::vector<double>& s, int i) {
            return finite3(s, i) && s[i] >= s[i - 1] && s[i] > s[i + 1];
        }
        inline bool isLocalMin(const std::vector<double>& s, int i) {
            return finite3(s, i) && s[i] <= s[i - 1] && s[i] < s[i + 1];
        }
        inline int firstLocalMin(const std::vector<double>& s, int lo, int hi) {
            lo = std::max(1, lo); hi = std::min(hi, static_cast<int>(s.size()) - 2);
            for (int i = lo; i <= hi; ++i) if (isLocalMin(s, i)) return i;
            return -1;
        }
        // Table 3's "local maximum with the highest amplitude", which is NOT
        // argmax_in: argmax over a monotone stretch returns an endpoint, and
        // an endpoint is not a local maximum. That distinction is the whole
        // reason c and d can legitimately come back -1.
        inline int highestLocalMax(const std::vector<double>& s, int lo, int hi) {
            lo = std::max(1, lo); hi = std::min(hi, static_cast<int>(s.size()) - 2);
            int best = -1; double bv = -std::numeric_limits<double>::infinity();
            for (int i = lo; i <= hi; ++i)
                if (isLocalMax(s, i) && s[i] > bv) { bv = s[i]; best = i; }
            return best;
        }
        inline int lowestLocalMin(const std::vector<double>& s, int lo, int hi) {
            lo = std::max(1, lo); hi = std::min(hi, static_cast<int>(s.size()) - 2);
            int best = -1; double bv = std::numeric_limits<double>::infinity();
            for (int i = lo; i <= hi; ++i)
                if (isLocalMin(s, i) && s[i] < bv) { bv = s[i]; best = i; }
            return best;
        }

    } // namespace detail

    // apg from compute_apg.
    //
    // on / sp / off are detect_ppg_fiducials' g.onset / g.peak / g.end. dp is
    // the diastolic peak -- pass g.peak2 only when g.peak2_found, and -1
    // otherwise, since that field holds a midpoint placeholder when the
    // detection fell through and a bound taken from a fabricated landmark is
    // worse than the duration cap this falls back to.
    inline ApgFiducials detect_apg(const std::vector<double>& apg,
        int on, int sp, int off, int Wc, int dp = -1)
    {
        ApgFiducials r;
        if (apg.empty() || on < 0 || sp <= on) return r;

        const int hiW = std::min(Wc, static_cast<int>(apg.size())) - 1;
        if (hiW <= sp) return r;
        auto cl = [&](int x) { return std::clamp(x, 0, hiW); };

        const int pulseEnd = (off > sp) ? std::min(off, hiW) : hiW;
        const double dur = static_cast<double>(pulseEnd - on);
        if (dur < 4.0) return r;
        // Table 4's upper bounds, as fractions of pulse duration: 60% for e,
        // 80% for f.
        const int cap60 = cl(on + static_cast<int>(std::lround(0.60 * dur)));
        const int cap80 = cl(on + static_cast<int>(std::lround(0.80 * dur)));

        // ---- a: highest amplitude on [on, sp] ---------------------------
        r.a = detail::argmax_in(apg, on, sp);
        if (r.a < 0) return r;          // no finite APG in systole

        // ---- b: first local minimum after a ----------------------------
        r.b = detail::firstLocalMin(apg, r.a + 1, cap60);
        if (r.b < 0) return r;          // b gates c, d and e

        // ---- e: tallest local max after b, before dp, capped at 60% -----
        {
            const int hi = (dp > r.b) ? std::min(dp, cap60) : cap60;
            r.e = detail::highestLocalMax(apg, r.b + 1, hi);
        }

        // ---- c: tallest local max on (b, e) ----------------------------
        if (r.e > r.b + 1)
            r.c = detail::highestLocalMax(apg, r.b + 1, r.e - 1);

        // ---- d: lowest local min on (c, e) -----------------------------
        if (r.c >= 0 && r.e > r.c + 1)
            r.d = detail::lowestLocalMin(apg, r.c + 1, r.e - 1);

        // ---- f: first local minimum after e, capped at 80% --------------
        if (r.e >= 0) r.f = detail::firstLocalMin(apg, r.e + 1, cap80);

        return r;
    }


    // =====================================================================
    // JPG (third derivative) fiducials: p1, p2
    // =====================================================================
    //
    // p1 is the early systolic component of the pulse wave and p2 the late
    // systolic component (Takazawa 1998); their amplitude ratio is the
    // augmentation index.
    //
    // BOTH ARE BOUNDED BY APG POINTS, so detect_apg must run first: p1 keys
    // off b, p2 off b and d. That makes p2 the most fragile point in the
    // whole set -- d is itself the least reliable APG point (it needs a true
    // local minimum between c and e), and when d is -1 there is no defensible
    // upper bound for p2, so p2 is -1 too. Substituting e for the missing d
    // would widen the window past the point table 3 specifies and admit
    // diastolic minima that are not p2.
    struct JpgFiducials {
        int p1 = -1;
        int p2 = -1;
    };

    namespace detail {
        inline int lastLocalMin(const std::vector<double>& s, int lo, int hi) {
            lo = std::max(1, lo); hi = std::min(hi, static_cast<int>(s.size()) - 2);
            for (int i = hi; i >= lo; --i) if (isLocalMin(s, i)) return i;
            return -1;
        }
    } // namespace detail

    // jpg from compute_jpg. b and d come from detect_apg; pass them straight
    // through, -1 and all -- a -1 bound disables the point that needs it
    // rather than falling back to a fabricated window.
    //
    // on / off / Wc as elsewhere; the 60%-of-duration cap bounds p1, matching
    // the cap on e (p1 is a systolic component, so a diastolic ripple is
    // never the answer).
    inline JpgFiducials detect_jpg(const std::vector<double>& jpg,
        int on, int off, int Wc, int b, int d)
    {
        JpgFiducials r;
        if (jpg.empty() || on < 0 || b < 0) return r;

        const int hiW = std::min(Wc, static_cast<int>(jpg.size())) - 1;
        if (hiW <= b) return r;
        auto cl = [&](int x) { return std::clamp(x, 0, hiW); };

        const int pulseEnd = (off > b) ? std::min(off, hiW) : hiW;
        const double dur = static_cast<double>(pulseEnd - on);
        if (dur < 4.0) return r;
        const int cap60 = cl(on + static_cast<int>(std::lround(0.60 * dur)));

        // ---- p1: first local maximum after b ----------------------------
        r.p1 = detail::firstLocalMax(jpg, b + 1, cap60);

        // ---- p2: last local minimum after b and before d ----------------
        if (d > b + 1) r.p2 = detail::lastLocalMin(jpg, b + 1, std::min(d - 1, hiW));

        return r;
    }

} // namespace ppg_deriv