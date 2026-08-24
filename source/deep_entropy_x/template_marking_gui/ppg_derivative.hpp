#pragma once
//
// ppg_derivative.hpp
//
// E-3 (Savitzky-Golay derivative bank) + E-4 (derivative fiducial points),
// DeepEntropyX Phase 2, Sections 6.2-6.3.
//
// E-3: all derivatives come from a Savitzky-Golay differentiator, not finite
// differences. Naive differencing amplifies high-frequency noise in proportion
// to frequency, and the 4th derivative that p1/p2 need is unusable without it.
// SG fits a local polynomial by least squares and differentiates it
// analytically, suppressing noise while preserving peak position.
//
// E-4: eleven fiducials on those derivatives, resolved as an ordered extremum
// walk (each point sought only after its predecessor):
//
//   PPG' (VPG, d1)   u  global maximum before the systolic peak
//                    v  global minimum after u
//                    w  first local maximum after v
//   PPG'' (APG, d2)  a  first positive peak (global max before systolic peak)
//                    b  first negative trough after a
//                    c  first local maximum after b
//                    d  first local minimum after c
//                    e  first local maximum after d   (dicrotic notch)
//                    f  first local minimum after e
//   4th derivative   p1 first zero-crossing of d4 after b  (forward wave)
//                    p2 first zero-crossing of d4 after e  (reflected wave)
//
// Derived indices (b/a..f/a, AGI, RI, SI, foundMask) follow in computeIndices().
//

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace ppg_deriv {

    // =====================================================================
    // E-3: Savitzky-Golay derivative bank
    // =====================================================================

    // Savitzky-Golay convolution kernel: fits a polynomial of order `order`
    // over [-h, +h] by least squares and returns the kernel for the `deriv`-th
    // derivative, already scaled (x deriv!) so a convolution yields the
    // derivative in units per sample.
    inline std::vector<double> sgCoeffs(int h, int order, int deriv) {
        const int n = 2 * h + 1, m = order + 1;
        // Vandermonde A[i][j] = (i - h)^j
        std::vector<std::vector<double>> A(n, std::vector<double>(m, 1.0));
        for (int i = 0; i < n; ++i)
            for (int j = 1; j < m; ++j) A[i][j] = A[i][j - 1] * (i - h);
        // Normal equations (A^T A) c = e_deriv; the deriv-th row of the
        // pseudoinverse is the kernel.
        std::vector<std::vector<double>> ATA(m, std::vector<double>(m + 1, 0.0));
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < m; ++k)
                for (int i = 0; i < n; ++i) ATA[j][k] += A[i][j] * A[i][k];
            ATA[j][m] = (j == deriv) ? 1.0 : 0.0;   // RHS selects the derivative
        }
        // Gauss-Jordan elimination with partial pivoting.
        for (int c = 0; c < m; ++c) {
            int piv = c;
            for (int r = c + 1; r < m; ++r)
                if (std::fabs(ATA[r][c]) > std::fabs(ATA[piv][c])) piv = r;
            std::swap(ATA[c], ATA[piv]);
            for (int r = 0; r < m; ++r) {
                if (r == c) continue;
                const double f = ATA[r][c] / ATA[c][c];
                for (int k = c; k <= m; ++k) ATA[r][k] -= f * ATA[c][k];
            }
        }
        std::vector<double> b(m);
        for (int j = 0; j < m; ++j) b[j] = ATA[j][m] / ATA[j][j];
        std::vector<double> kern(n, 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) kern[i] += b[j] * A[i][j];
        double fact = 1.0; for (int dd = 2; dd <= deriv; ++dd) fact *= dd;   // deriv!
        for (double& v : kern) v *= fact;
        return kern;
    }

    struct DerivBank {
        std::vector<double> d0, d1, d2, d3, d4;   // d0 is the SG-smoothed signal
    };

    // hLow / hHigh are the half-widths at a 500 Hz REFERENCE rate (config fields
    // sg_halfwidth_low = 12 for d1/d2, sg_halfwidth_high = 20 for d3/d4), and
    // order is sg_poly_order. The half-width is scaled by fs/500 so the effective
    // time window -- and thus the smoothing -- is identical across datasets.
    inline DerivBank buildDerivatives(const std::vector<double>& x, double fs,
        int hLow = 12, int hHigh = 20, int order = 4)
    {
        constexpr double kRefHz = 500.0;
        const int hMin = std::max(2, (order + 1) / 2);   // need 2h+1 >= order+1
        auto scaleH = [&](int hRef) {
            const int h = (fs > 0.0)
                ? static_cast<int>(std::lround(hRef * fs / kRefHz)) : hRef;
            return std::max(h, hMin);
            };
        const int hL = scaleH(hLow);
        const int hH = scaleH(hHigh);
        const int N = static_cast<int>(x.size());

        auto conv = [&](const std::vector<double>& k, int h, int d) {
            std::vector<double> y(N, 0.0);
            const double scale = std::pow(fs, d);   // per-sample -> per-second
            for (int i = 0; i < N; ++i) {
                double s = 0.0;
                for (int j = -h; j <= h; ++j) {
                    const int idx = std::clamp(i + j, 0, N - 1);   // edge replication
                    s += k[j + h] * x[idx];
                }
                y[i] = s * scale;
            }
            return y;
            };

        DerivBank b;
        if (N == 0) return b;
        b.d0 = conv(sgCoeffs(hL, order, 0), hL, 0);
        b.d1 = conv(sgCoeffs(hL, order, 1), hL, 1);
        b.d2 = conv(sgCoeffs(hL, order, 2), hL, 2);
        b.d3 = conv(sgCoeffs(hH, order, 3), hH, 3);
        b.d4 = conv(sgCoeffs(hH, order, 4), hH, 4);
        return b;
    }

    // =====================================================================
    // E-4: fiducial points
    // =====================================================================

    struct Fiducials {
        int u = -1, v = -1, w = -1;              // d1 (VPG)
        int a = -1, b = -1, c = -1, d = -1, e = -1, f = -1;   // d2 (APG)
        int p1 = -1, p2 = -1;                    // d4
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
        // First sign change of the value between two finite adjacent samples,
        // strictly after `from`, up to `hi`. NaN gaps are not crossings.
        inline int firstZeroCross(const std::vector<double>& s, int from, int hi) {
            const int lo = std::max(from, 1);
            hi = std::min(hi, static_cast<int>(s.size()) - 1);
            for (int i = lo; i <= hi; ++i) {
                if (std::isnan(s[i - 1]) || std::isnan(s[i])) continue;
                if ((s[i - 1] > 0.0) != (s[i] > 0.0)) return i;
            }
            return -1;
        }

    } // namespace detail

    // on     pulse onset  -- UNUSED: u and a are bounded by [0, systolic peak).
    // sp     systolic peak (upper bound of the u/a search)
    // dn     dicrotic notch -- UNUSED by these definitions.
    // dp     diastolic peak -- UNUSED by these definitions.
    // Wc     visible-window clamp; nothing is placed past Wc-1, and the spec's
    //        "end of pulse" maps to this window's end.
    inline Fiducials detect(const DerivBank& d, int on, int sp, int dn, int dp, int Wc)
    {
        (void)on; (void)dn; (void)dp;

        Fiducials r;
        if (d.d2.empty() || d.d1.empty() || d.d4.empty()) return r;
        if (sp <= 0) return r;

        const int nA = static_cast<int>(d.d2.size());
        const int hiW = std::min(Wc, nA) - 1;
        if (hiW <= sp) return r;

        const std::vector<double>& V = d.d1;   // VPG
        const std::vector<double>& A = d.d2;   // APG
        const std::vector<double>& D4 = d.d4;  // 4th derivative

        const int prePeakHi = std::min(sp - 1, hiW);   // [0, sp)

        // ---- d1 (VPG): u -> v -> w ----
        r.u = detail::argmax_in(V, 0, prePeakHi);
        if (r.u >= 0 && hiW > r.u) r.v = detail::argmin_in(V, r.u + 1, hiW);
        if (r.v >= 0) r.w = detail::firstLocalMax(V, r.v + 1, hiW);

        // ---- d2 (APG): a -> b -> c -> d -> e -> f ----
        r.a = detail::argmax_in(A, 0, prePeakHi);
        if (r.a >= 0) r.b = detail::firstLocalMin(A, r.a + 1, hiW);
        if (r.b >= 0) r.c = detail::firstLocalMax(A, r.b + 1, hiW);
        if (r.c >= 0) r.d = detail::firstLocalMin(A, r.c + 1, hiW);
        if (r.d >= 0) r.e = detail::firstLocalMax(A, r.d + 1, hiW);
        if (r.e >= 0) r.f = detail::firstLocalMin(A, r.e + 1, hiW);

        // ---- d4: p1, p2 ----
        if (r.b >= 0) r.p1 = detail::firstZeroCross(D4, r.b + 1, hiW);
        if (r.e >= 0) r.p2 = detail::firstZeroCross(D4, r.e + 1, hiW);

        return r;
    }

    // =====================================================================
    // Derived indices (Section 6.3)
    // =====================================================================

    struct Indices {
        double ba = std::numeric_limits<double>::quiet_NaN();
        double ca = std::numeric_limits<double>::quiet_NaN();
        double da = std::numeric_limits<double>::quiet_NaN();
        double ea = std::numeric_limits<double>::quiet_NaN();
        double fa = std::numeric_limits<double>::quiet_NaN();
        double agi = std::numeric_limits<double>::quiet_NaN();
        double ri = std::numeric_limits<double>::quiet_NaN();
        double si = std::numeric_limits<double>::quiet_NaN();
        uint16_t foundMask = 0;   // bits 0..10: u,v,w,a,b,c,d,e,f,p1,p2
    };

    inline uint16_t buildFoundMask(const Fiducials& F) {
        const int pts[11] = { F.u, F.v, F.w, F.a, F.b, F.c, F.d, F.e, F.f, F.p1, F.p2 };
        uint16_t m = 0;
        for (int k = 0; k < 11; ++k)
            if (pts[k] >= 0) m |= static_cast<uint16_t>(1u << k);
        return m;
    }

    // pulse         original PPG pulse (RI reads its amplitude at p1/p2)
    // d             the derivative bank (ratios/AGI use APG amplitude, d2)
    // heightMeters  subject height; NaN or <= 0 leaves SI as NaN (no population
    //               substitute), per spec.
    inline Indices computeIndices(const std::vector<double>& pulse,
        const DerivBank& d, const Fiducials& F,
        double fs, double heightMeters)
    {
        Indices I;
        I.foundMask = buildFoundMask(F);

        const std::vector<double>& A = d.d2;
        const int nA = static_cast<int>(A.size());
        auto apgOk = [&](int i) { return i >= 0 && i < nA && !std::isnan(A[i]); };

        if (apgOk(F.a) && std::fabs(A[F.a]) > 1e-9) {
            const double a = A[F.a];
            auto r = [&](int i) { return apgOk(i) ? A[i] / a : std::numeric_limits<double>::quiet_NaN(); };
            I.ba = r(F.b); I.ca = r(F.c); I.da = r(F.d); I.ea = r(F.e); I.fa = r(F.f);
            if (apgOk(F.b) && apgOk(F.c) && apgOk(F.d) && apgOk(F.e))
                I.agi = (A[F.b] - A[F.c] - A[F.d] - A[F.e]) / a;
        }

        const int nP = static_cast<int>(pulse.size());
        auto pulseOk = [&](int i) { return i >= 0 && i < nP && !std::isnan(pulse[i]); };
        if (fs > 0.0 && pulseOk(F.p1) && pulseOk(F.p2) && std::fabs(pulse[F.p1]) > 1e-9) {
            I.ri = pulse[F.p2] / pulse[F.p1];
            const double dt = (F.p2 - F.p1) / fs;
            I.si = (dt > 1e-6 && heightMeters > 0.0 && !std::isnan(heightMeters))
                ? heightMeters / dt
                : std::numeric_limits<double>::quiet_NaN();
        }
        return I;
    }

} // namespace ppg_deriv