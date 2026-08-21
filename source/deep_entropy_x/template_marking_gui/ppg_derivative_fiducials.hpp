#pragma once
//
// ppg_derivative_fiducials.hpp
//
// Fiducial points on the PPG derivatives, implemented STRICTLY from table 3 of
// Goda, Charlton & Behar, "pyPPG", Physiol. Meas. 45 (2024) 045001. The
// definitions, verbatim:
//
//   PPG' (VPG)
//     u   The highest amplitude between the pulse onset and systolic peak
//     v   The lowest amplitude between the u-point and diastolic peak
//     w   The first local maximum or inflection point after the dicrotic notch
//
//   PPG'' (APG)
//     a   The highest amplitude between pulse onset and systolic peak
//     b   The first local minimum after the a-point
//     c   The local maximum with the highest amplitude between the b-point and
//         e-point, or if no local maximum is present, then the inflection point
//     d   The local minimum with the lowest amplitude between the c-point and
//         e-point, or if no local minimum is present, then the inflection point
//     e   The local maximum with the highest amplitude after the b-point and
//         before the diastolic peak
//     f   The first local minimum after the e-point
//
//   PPG''' (JPG)
//     p1  The first local maximum after the b-point
//     p2  The last local minimum after the b-point and before the d-point

//

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace ppg_deriv {

    // ---------------------------------------------------------------------
    // Derivatives
    // ---------------------------------------------------------------------
    // Savitzky-Golay first difference, chained. For an SG fit the centre
    // first-derivative coefficients of the quadratic and linear fits coincide
    // (even-order terms cancel by symmetry), so this reduces to the
    // least-squares slope, sum_k k*v[i+k] / sum_k k^2 -- no matrix solve.
    //
    // The window is a time constant, so the operator is sample-rate
    // independent. pyPPG applies a 10 ms moving average to each derivative
    // stage (~45 Hz cutoff); kSmoothMs matches that, with a flatter passband
    // than the boxcar.
    inline constexpr double kSmoothMs = 10.0;

    inline std::vector<double> derivative(const std::vector<double>& x, double fs,
        double smoothMs = kSmoothMs) {
        const int N = static_cast<int>(x.size());
        std::vector<double> d(N, std::numeric_limits<double>::quiet_NaN());
        if (fs <= 0.0) return d;
        const int m = std::max(2, static_cast<int>(std::lround(0.5 * (smoothMs / 1000.0) * fs)));
        if (N < 2 * m + 1) return d;
        double norm = 0.0;
        for (int k = -m; k <= m; ++k) norm += static_cast<double>(k) * k;
        for (int i = m; i < N - m; ++i) {
            double s = 0.0;
            bool ok = true;
            for (int k = -m; k <= m && ok; ++k) {
                const double val = x[i + k];
                if (std::isnan(val)) ok = false; else s += k * val;
            }
            // NaN propagates rather than being interpolated across: template
            // tails carry NaN by construction, and inventing samples there
            // would put phantom extrema inside the search windows.
            if (ok) d[i] = s / norm;
        }
        return d;
    }

    struct Derivatives {
        std::vector<double> vpg;   // PPG'
        std::vector<double> apg;   // PPG''
        std::vector<double> jpg;   // PPG'''
    };

    // fs is the PPG channel rate in Hz (ppgRate at the call site, NOT the ECG
    // sample rate -- the two differ in every dataset this tree handles).
    inline Derivatives compute(const std::vector<double>& ppg, double fs,
        double smoothMs = kSmoothMs) {
        Derivatives d;
        d.vpg = derivative(ppg, fs, smoothMs);
        d.apg = derivative(d.vpg, fs, smoothMs);
        d.jpg = derivative(d.apg, fs, smoothMs);
        return d;
    }

    // ---------------------------------------------------------------------
    // Results
    // ---------------------------------------------------------------------
    struct Fiducials {
        int u = -1, v = -1, w = -1;              // PPG'
        int a = -1, b = -1, c = -1, d = -1, e = -1, f = -1;   // PPG''
        int p1 = -1, p2 = -1;                    // PPG'''
    };

    // ---------------------------------------------------------------------
    // Extremum helpers, NaN-skipping
    // ---------------------------------------------------------------------
    namespace detail {

        // "The highest / lowest amplitude" -- a plain argmax/argmin over the
        // range, endpoints included. Table 3 uses this phrasing for u, v and a.
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

        inline bool finite3(const std::vector<double>& s, int i) {
            return !std::isnan(s[i - 1]) && !std::isnan(s[i]) && !std::isnan(s[i + 1]);
        }
        inline bool isLocalMax(const std::vector<double>& s, int i) {
            return finite3(s, i) && s[i] >= s[i - 1] && s[i] > s[i + 1];
        }
        inline bool isLocalMin(const std::vector<double>& s, int i) {
            return finite3(s, i) && s[i] <= s[i - 1] && s[i] < s[i + 1];
        }

        inline int firstLocalMax(const std::vector<double>& s, int lo, int hi) {
            lo = std::max(1, lo); hi = std::min(hi, static_cast<int>(s.size()) - 2);
            for (int i = lo; i <= hi; ++i) if (isLocalMax(s, i)) return i;
            return -1;
        }
        inline int firstLocalMin(const std::vector<double>& s, int lo, int hi) {
            lo = std::max(1, lo); hi = std::min(hi, static_cast<int>(s.size()) - 2);
            for (int i = lo; i <= hi; ++i) if (isLocalMin(s, i)) return i;
            return -1;
        }
        inline int lastLocalMin(const std::vector<double>& s, int lo, int hi) {
            lo = std::max(1, lo); hi = std::min(hi, static_cast<int>(s.size()) - 2);
            for (int i = hi; i >= lo; --i) if (isLocalMin(s, i)) return i;
            return -1;
        }

        // "The local maximum with the highest amplitude" (c, e) and "the local
        // minimum with the lowest amplitude" (d). NOT argmax/argmin: over a
        // monotone stretch those return an endpoint, and an endpoint is not a
        // local extremum. That distinction is precisely why table 3 gives c and
        // d an inflection fallback.
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

        // "The inflection point" of a signal, on [lo, hi]. An inflection is a
        // stationary point of the derivative, so it is an extremum one order
        // up: for a flattening on a descending curve -- the shoulder case that
        // makes c and d fall through -- that is a MAXIMUM of the next
        // derivative (least-negative slope), whichever way the absent feature
        // would have bent. Returns -1 if the next derivative is unavailable.
        inline int inflection(const std::vector<double>& nextDeriv, int lo, int hi) {
            if (nextDeriv.empty()) return -1;
            return highestLocalMax(nextDeriv, lo, hi);
        }

    } // namespace detail

    // ---------------------------------------------------------------------
    // Detection
    // ---------------------------------------------------------------------
    //
    // d      the three derivatives, from compute()
    // on     pulse onset          (detect_ppg_fiducials' g.onset)
    // sp     systolic peak        (g.peak)
    // dn     dicrotic notch       (g.dicrotic)   -- bounds w
    // dp     diastolic peak       (g.peak2)      -- bounds v and e
    // Wc     visible-window clamp; nothing is placed past Wc-1
    //
    // dn and dp are passed as-is. They are the paper's bounds and this function
    // uses them; it does not inspect any *_found flag, and it does not
    // substitute anything when one is -1 -- it simply searches to the end of
    // the window, since table 3 supplies no alternative bound.
    inline Fiducials detect(const Derivatives& d, int on, int sp, int dn, int dp, int Wc)
    {
        Fiducials r;
        if (d.apg.empty() || d.vpg.empty()) return r;
        if (on < 0 || sp <= on) return r;

        const int hiW = std::min(Wc, static_cast<int>(d.apg.size())) - 1;
        if (hiW <= sp) return r;

        const std::vector<double>& V = d.vpg;
        const std::vector<double>& A = d.apg;
        const std::vector<double>& J = d.jpg;

        // Upper bound for the definitions that name the diastolic peak (v, e).
        // No dp => no bound, so search the window.
        const int dpHi = (dp > sp) ? std::min(dp, hiW) : hiW;

        // ===== PPG'' ======================================================

        // a: highest amplitude between pulse onset and systolic peak.
        r.a = detail::argmax_in(A, on, sp);
        if (r.a < 0) return r;

        // b: first local minimum after the a-point.
        r.b = detail::firstLocalMin(A, r.a + 1, hiW);
        if (r.b < 0) return r;          // b bounds c, d, e, p1 and p2

        // e: local maximum with the highest amplitude after the b-point and
        //    before the diastolic peak.
        r.e = detail::highestLocalMax(A, r.b + 1, dpHi);

        // c: local maximum with the highest amplitude between b and e; if no
        //    local maximum is present, the inflection point.
        if (r.e > r.b + 1) {
            r.c = detail::highestLocalMax(A, r.b + 1, r.e - 1);
            if (r.c < 0) r.c = detail::inflection(J, r.b + 1, r.e - 1);
        }

        // d: local minimum with the lowest amplitude between c and e; if no
        //    local minimum is present, the inflection point.
        if (r.c >= 0 && r.e > r.c + 1) {
            r.d = detail::lowestLocalMin(A, r.c + 1, r.e - 1);
            if (r.d < 0) r.d = detail::inflection(J, r.c + 1, r.e - 1);
        }

        // f: first local minimum after the e-point.
        if (r.e >= 0) r.f = detail::firstLocalMin(A, r.e + 1, hiW);

        // ===== PPG' =======================================================

        // u: highest amplitude between the pulse onset and systolic peak.
        r.u = detail::argmax_in(V, on, sp);

        // v: lowest amplitude between the u-point and diastolic peak.
        if (r.u >= 0 && dpHi > r.u) r.v = detail::argmin_in(V, r.u + 1, dpHi);

        // w: first local maximum or inflection point after the dicrotic notch.
        if (dn >= 0 && dn < hiW) {
            r.w = detail::firstLocalMax(V, dn + 1, hiW);
            if (r.w < 0) r.w = detail::inflection(A, dn + 1, hiW);
        }

        // ===== PPG''' =====================================================

        // p1: first local maximum after the b-point.
        r.p1 = detail::firstLocalMax(J, r.b + 1, hiW);

        // p2: last local minimum after the b-point and before the d-point.
        if (r.d > r.b + 1) r.p2 = detail::lastLocalMin(J, r.b + 1, r.d - 1);

        return r;
    }

} // namespace ppg_deriv