#pragma once
// ═══════════════════════════════════════════════════════════════════════════════
// ppg_utils.hpp — Utility functions faithfully matching MATLAB behavior.
// Includes workspace-based variants of hot-path functions to avoid
// per-beat heap allocations in the SQI inner loop.
// ═══════════════════════════════════════════════════════════════════════════════
#include "ppg_features.hpp"

namespace ppg {

    // ─── 1D Linear Interpolation ────────────────────────────────────────────────
    inline std::vector<double> interp1_linear(const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& xq) {
        std::vector<double> result(xq.size());
        for (size_t i = 0; i < xq.size(); i++) {
            double xi = xq[i];
            if (xi <= x.front()) { result[i] = y.front(); continue; }
            if (xi >= x.back()) { result[i] = y.back();  continue; }
            auto it = std::lower_bound(x.begin(), x.end(), xi);
            int j = (int)(it - x.begin());
            if (j == 0) j = 1;
            double t = (xi - x[j - 1]) / (x[j] - x[j - 1]);
            result[i] = y[j - 1] + t * (y[j] - y[j - 1]);
        }
        return result;
    }

    // Workspace version: writes into pre-allocated output buffer
    inline void interp1_linear_ws(const double* x, const double* y, int n,
        const double* xq, double* out, int nq) {
        for (int i = 0; i < nq; i++) {
            double xi = xq[i];
            if (xi <= x[0]) { out[i] = y[0];     continue; }
            if (xi >= x[n - 1]) { out[i] = y[n - 1];   continue; }
            int lo = 0, hi = n - 1;
            while (lo < hi - 1) {
                int mid = (lo + hi) / 2;
                if (x[mid] <= xi) lo = mid; else hi = mid;
            }
            double t = (xi - x[lo]) / (x[hi] - x[lo]);
            out[i] = y[lo] + t * (y[hi] - y[lo]);
        }
    }

    // ─── InterX-style crossing finder ──────────────────────────────────────────
    inline std::pair<double, double> interx_last_crossing(
        const double* y_data, int n, double level) {
        double last_x = NaN;
        for (int k = 0; k < n - 1; k++) {
            double y0 = y_data[k] - level;
            double y1 = y_data[k + 1] - level;
            if ((y0 <= 0 && y1 >= 0) || (y0 >= 0 && y1 <= 0)) {
                if (y0 == y1) last_x = (double)(k + 1);
                else { double t = -y0 / (y1 - y0); last_x = (double)(k + 1) + t; }
            }
        }
        return { last_x, std::isnan(last_x) ? NaN : level };
    }

    inline std::pair<double, double> interx_last_crossing_range(
        const std::vector<double>& beat, int from, int to, double level) {
        double last_x = NaN;
        for (int k = from; k < to - 1; k++) {
            double y0 = beat[k] - level;
            double y1 = beat[k + 1] - level;
            if ((y0 <= 0 && y1 >= 0) || (y0 >= 0 && y1 <= 0)) {
                if (y0 == y1) last_x = (double)(k - from + 1);
                else { double t = -y0 / (y1 - y0); last_x = (double)(k - from + 1) + t; }
            }
        }
        return { last_x, std::isnan(last_x) ? NaN : level };
    }

    // ─── PLA (exact MATLAB PLA.m) ───────────────────────────────────────────────
    inline std::pair<std::vector<double>, std::vector<int>> PLA(
        const std::vector<double>& input, int s, double th) {
        int n = (int)input.size();
        if (n == 0) return { input, {0} };
        std::vector<int> pla_points;
        pla_points.reserve(n / s + 10);
        pla_points.push_back(0);
        int s1 = s, i = 0;
        while (i < n) {
            int i_plus_s = std::min(i + s1, n - 1);
            bool interrupt = false;
            while (!interrupt) {
                int j = i + 1;
                while (j <= i_plus_s) {
                    double distance = input[i_plus_s] - input[i];
                    double denom = (double)(i_plus_s - i);
                    if (denom == 0) { j++; continue; }
                    double dcur = input[j] - input[i] - (distance * (j - i) / denom);
                    if (std::abs(dcur) > th) {
                        s1 = j - i; i_plus_s = i + s1; j = i + 1; interrupt = true; continue;
                    }
                    j++;
                }
                if (interrupt) { pla_points.push_back(j - 2); i = j - 3; s1 = s; }
                else { if (i_plus_s >= n - 1) break; else i_plus_s = std::min(i_plus_s + s1, n - 1); }
            }
            i++;
        }
        pla_points.push_back(n - 1);
        std::sort(pla_points.begin(), pla_points.end());
        pla_points.erase(std::unique(pla_points.begin(), pla_points.end()), pla_points.end());
        return { input, pla_points };
    }

    // Workspace PLA: writes indices into pre-allocated buffer, returns count
    inline int PLA_ws(const double* input, int n, int s, double th, int* pla_out, int max_pla) {
        if (n == 0) { pla_out[0] = 0; return 1; }
        int pp = 0;
        pla_out[pp++] = 0;
        int s1 = s, i = 0;
        while (i < n) {
            int i_plus_s = std::min(i + s1, n - 1);
            bool interrupt = false;
            while (!interrupt) {
                int j = i + 1;
                while (j <= i_plus_s) {
                    double distance = input[i_plus_s] - input[i];
                    double denom = (double)(i_plus_s - i);
                    if (denom == 0) { j++; continue; }
                    double dcur = input[j] - input[i] - (distance * (j - i) / denom);
                    if (std::abs(dcur) > th) {
                        s1 = j - i; i_plus_s = i + s1; j = i + 1; interrupt = true; continue;
                    }
                    j++;
                }
                if (interrupt) { if (pp < max_pla) pla_out[pp++] = j - 2; i = j - 3; s1 = s; }
                else { if (i_plus_s >= n - 1) break; else i_plus_s = std::min(i_plus_s + s1, n - 1); }
            }
            i++;
        }
        if (pp < max_pla) pla_out[pp++] = n - 1;
        std::sort(pla_out, pla_out + pp);
        pp = (int)(std::unique(pla_out, pla_out + pp) - pla_out);
        return pp;
    }

    // ─── Smoothing ──────────────────────────────────────────────────────────────
    inline std::vector<double> nanfastsmooth(const std::vector<double>& data, int width) {
        int n = (int)data.size();
        if (n == 0) return data;
        std::vector<double> out(n);
        int half = width / 2;
        for (int i = 0; i < n; i++) {
            double sum = 0; int count = 0;
            for (int j = std::max(0, i - half); j <= std::min(n - 1, i + half); j++)
                if (!std::isnan(data[j])) { sum += data[j]; count++; }
            out[i] = (count > 0) ? sum / count : NaN;
        }
        return out;
    }

    // ─── PPG Median Filter ─────────────────────────────────────────────────────
    // The actual MATLAB PPGmedianfilter source is not available in the uploaded files.
    // The original C++ implementation was a no-op. We keep it as a no-op to match.
    // If you have the MATLAB source for PPGmedianfilter, implement it here.
    inline std::vector<double> PPGmedianfilter(const std::vector<double>& wave, int /*order*/, double /*Fs*/) {
        return wave;
    }

    // ─── Per-thread SQI workspace ───────────────────────────────────────────────
    // Allocated ONCE per thread, reused across all beats. No heap allocs in hot loop.
    struct SqiWorkspace {
        static constexpr int MAX_BEAT = 2048;  // 3*Fs at 512Hz + margin
        static constexpr int MAX_PLA = 512;   // PLA typically produces 100-300 points
        static constexpr int MAX_DTW = MAX_PLA * MAX_PLA; // 262144

        // Beat buffers
        double bx[MAX_BEAT], by[MAX_BEAT], ix_buf[MAX_BEAT], yi[MAX_BEAT];
        double d2[MAX_BEAT], d2n[MAX_BEAT];
        double y1n[MAX_BEAT], y2n[MAX_BEAT];
        double interp_out[MAX_BEAT];
        double xo[MAX_BEAT], xf[MAX_BEAT];

        // PLA
        int pla2[MAX_PLA];

        // DTW matrices
        double w[MAX_DTW];
        double D[MAX_DTW];
        int tbi[MAX_DTW], tbj[MAX_DTW];
        double ym2[MAX_BEAT];
        int path_p[MAX_PLA * 2], path_q[MAX_PLA * 2];

        // Frechet — stack buffer for small cases, heap fallback for large
        double ca_stack[MAX_BEAT * MAX_BEAT > 0 ? 1 : 1]; // placeholder
        std::vector<double> ca_heap;

        double* get_frechet_buf(int n, int m) {
            size_t needed = (size_t)n * m;
            // Always use heap for Frechet since n*m can be large
            if (ca_heap.size() < needed) ca_heap.resize(needed);
            return ca_heap.data();
        }
    };

    // ─── Workspace-based DTW pipeline ───────────────────────────────────────────

    inline void simmx_dtw_ws(const double* y1, const int* pla1, int n1,
        const double* y2, const int* pla2, int n2,
        double* w_out) {
        for (int i = 0; i < n1 - 1; i++) {
            double slope1 = (pla1[i + 1] != pla1[i]) ?
                (y1[pla1[i + 1]] - y1[pla1[i]]) / (double)(pla1[i + 1] - pla1[i]) : 0.0;
            for (int j = 0; j < n2 - 1; j++) {
                double slope2 = (pla2[j + 1] != pla2[j]) ?
                    (y2[pla2[j + 1]] - y2[pla2[j]]) / (double)(pla2[j + 1] - pla2[j]) : 0.0;
                w_out[i * n2 + j] = std::abs(slope1 - slope2);
            }
        }
    }

    inline int dp_dtw2_ws(const double* w, int n, int m,
        double* D_buf, int* tbi_buf, int* tbj_buf,
        int* p_out, int* q_out) {
        for (int k = 0; k < n * m; k++) { D_buf[k] = 1e18; tbi_buf[k] = -1; tbj_buf[k] = -1; }
        auto d = [&](int i, int j) -> double& { return D_buf[i * m + j]; };
        auto ti = [&](int i, int j) -> int& { return tbi_buf[i * m + j]; };
        auto tj = [&](int i, int j) -> int& { return tbj_buf[i * m + j]; };

        d(0, 0) = w[0];
        for (int i = 1; i < n; i++) { d(i, 0) = d(i - 1, 0) + w[i * m]; ti(i, 0) = i - 1; tj(i, 0) = 0; }
        for (int j = 1; j < m; j++) { d(0, j) = d(0, j - 1) + w[j]; ti(0, j) = 0; tj(0, j) = j - 1; }
        for (int i = 1; i < n; i++)
            for (int j = 1; j < m; j++) {
                double c[3] = { d(i - 1,j), d(i,j - 1), d(i - 1,j - 1) };
                int ci[3] = { i - 1,i,i - 1 }, cj[3] = { j,j - 1,j - 1 };
                int best = 0;
                for (int k = 1; k < 3; k++) if (c[k] < c[best]) best = k;
                d(i, j) = c[best] + w[i * m + j]; ti(i, j) = ci[best]; tj(i, j) = cj[best];
            }
        int len = 0, ci2 = n - 1, cj2 = m - 1;
        while (ci2 >= 0 && cj2 >= 0) {
            p_out[len] = ci2; q_out[len] = cj2; len++;
            int ni = ti(ci2, cj2), nj = tj(ci2, cj2);
            if (ni < 0) break; ci2 = ni; cj2 = nj;
        }
        for (int i = 0; i < len / 2; i++) { std::swap(p_out[i], p_out[len - 1 - i]); std::swap(q_out[i], q_out[len - 1 - i]); }
        return len;
    }

    inline void draw_dtw_ws(const double* y1, const int* pla1, int n1_pla,
        const int* p, const double* y2, const int* pla2, int n2_pla,
        const int* q, int path_len,
        double* ym2, int n1_sig) {
        std::memset(ym2, 0, n1_sig * sizeof(double));
        int y2_max = 0; // we need y2 size — infer from pla2
        for (int k = 0; k < n2_pla; k++) if (pla2[k] > y2_max) y2_max = pla2[k];
        for (int k = 0; k < path_len - 1; k++) {
            int i1s = (p[k] < n1_pla) ? pla1[p[k]] : pla1[n1_pla - 1];
            int i1e = (p[k + 1] < n1_pla) ? pla1[p[k + 1]] : pla1[n1_pla - 1];
            int i2s = (q[k] < n2_pla) ? pla2[q[k]] : pla2[n2_pla - 1];
            int i2e = (q[k + 1] < n2_pla) ? pla2[q[k + 1]] : pla2[n2_pla - 1];
            for (int idx = i1s; idx <= std::min(i1e, n1_sig - 1); idx++) {
                double t = (i1e != i1s) ? (double)(idx - i1s) / (i1e - i1s) : 0.0;
                int src = i2s + (int)(t * (i2e - i2s));
                src = std::max(0, std::min(src, y2_max));
                ym2[idx] = y2[src];
            }
        }
    }

    inline double DiscreteFrechetDist_ws(const double* P, int n, const double* Q, int m, double* ca) {
        if (n == 0 || m == 0) return 0;
        auto at = [&](int i, int j) -> double& { return ca[i * m + j]; };
        at(0, 0) = std::abs(P[0] - Q[0]);
        for (int i = 1; i < n; i++) at(i, 0) = std::max(at(i - 1, 0), std::abs(P[i] - Q[0]));
        for (int j = 1; j < m; j++) at(0, j) = std::max(at(0, j - 1), std::abs(P[0] - Q[j]));
        for (int i = 1; i < n; i++)
            for (int j = 1; j < m; j++)
                at(i, j) = std::max(std::min({ at(i - 1,j), at(i - 1,j - 1), at(i,j - 1) }), std::abs(P[i] - Q[j]));
        return at(n - 1, m - 1);
    }

    // ─── Non-workspace versions (used outside hot path) ─────────────────────────

    inline DtwResult simmx_dtw(const std::vector<double>& y1, const std::vector<int>& pla1,
        const std::vector<double>& y2, const std::vector<int>& pla2) {
        int n1 = (int)pla1.size(), n2 = (int)pla2.size();
        DtwResult r; r.rows = n1; r.cols = n2;
        r.w.assign((size_t)n1 * n2, 0.0); r.ta = pla1; r.tb = pla2;
        simmx_dtw_ws(y1.data(), pla1.data(), n1, y2.data(), pla2.data(), n2, r.w.data());
        return r;
    }
    inline DpResult dp_dtw2(const DtwResult& dr) {
        int n = dr.rows, m = dr.cols;
        std::vector<double> D(n * m); std::vector<int> tbi(n * m), tbj(n * m), p(n + m), q(n + m);
        int len = dp_dtw2_ws(dr.w.data(), n, m, D.data(), tbi.data(), tbj.data(), p.data(), q.data());
        p.resize(len); q.resize(len);
        return { p, q, D[n * m - 1] };
    }
    inline std::vector<double> draw_dtw(const std::vector<double>& y1, const std::vector<int>& pla1,
        const std::vector<int>& p,
        const std::vector<double>& y2, const std::vector<int>& pla2,
        const std::vector<int>& q) {
        int n1 = (int)y1.size();
        std::vector<double> ym2(n1, 0.0);
        draw_dtw_ws(y1.data(), pla1.data(), (int)pla1.size(), p.data(), y2.data(), pla2.data(), (int)pla2.size(), q.data(), (int)p.size(), ym2.data(), n1);
        return ym2;
    }
    inline double DiscreteFrechetDist(const std::vector<double>& P, const std::vector<double>& Q) {
        int n = (int)P.size(), m = (int)Q.size();
        if (n == 0 || m == 0) return 0;
        std::vector<double> ca((size_t)n * m);
        return DiscreteFrechetDist_ws(P.data(), n, Q.data(), m, ca.data());
    }

    // ─── Shear Transform ────────────────────────────────────────────────────────
    inline std::vector<double> shear_transform(const std::vector<int>&, const std::vector<double>& y_vals) {
        int n = (int)y_vals.size();
        if (n < 2) return y_vals;
        double slope = (y_vals.back() - y_vals.front()) / (double)(n - 1);
        std::vector<double> transformed(n);
        for (int i = 0; i < n; i++) transformed[i] = y_vals[i] - (y_vals.front() + slope * i);
        return transformed;
    }
    inline bool orthoginal_dist_thresh(const std::vector<double>&, const std::vector<double>& norm_line,
        const std::vector<double>& norm_press, double thresh) {
        for (size_t i = 0; i < norm_press.size() && i < norm_line.size(); i++)
            if (std::abs(norm_press[i] - norm_line[i]) > thresh) return true;
        return false;
    }

    // ─── Run-Length Encoding ────────────────────────────────────────────────────
    inline RunLengthResult RunLength(const std::vector<double>& data) {
        RunLengthResult result;
        if (data.empty()) return result;
        double cur = data[0]; int count = 1;
        for (size_t i = 1; i < data.size(); i++) {
            bool same = (std::isnan(cur) && std::isnan(data[i])) || (cur == data[i]);
            if (same) count++; else { result.values.push_back(cur); result.lengths.push_back(count); cur = data[i]; count = 1; }
        }
        result.values.push_back(cur); result.lengths.push_back(count);
        return result;
    }

    // ─── Correlation Coefficient ────────────────────────────────────────────────
    inline double corrcoef(const double* a, const double* b, int len) {
        if (len < 2) return 0;
        double ma = 0, mb = 0;
        for (int i = 0; i < len; i++) { ma += a[i]; mb += b[i]; }
        ma /= len; mb /= len;
        double num = 0, da = 0, db = 0;
        for (int i = 0; i < len; i++) { double ai = a[i] - ma, bi = b[i] - mb; num += ai * bi; da += ai * ai; db += bi * bi; }
        double denom = std::sqrt(da * db);
        return (denom > 0) ? num / denom : 0;
    }

} // namespace ppg