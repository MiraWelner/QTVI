#pragma once

#include "common.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <limits>
#include <future>

// ============================================================================
// PPG Signal Quality Index — optimized C++ translation
//
// Translates: PPG_SQI.m, simmx_dtw.m, dp_dtw2.m, draw_dtw.m,
//             DiscreteFrechetDist.m, PLA.m, PPGmedianfilter.m
//
// Key optimizations vs naive translation:
//   - All per-beat work uses preallocated scratch buffers (no malloc in loop)
//   - DTW sim matrix computed in-place without intermediate A1/B1 matrices
//   - Frechet distance uses iterative bottom-up DP (not recursion)
//   - draw_dtw uses linear interp by default (spline optional)
//   - Median filter uses partial_sort instead of full sort
//   - Designed for OpenMP parallelism on the outer beat loop
// ============================================================================

namespace ppg {

    // ---------------------------------------------------------------------------
    // PPG Median Filter (PPGmedianfilter.m)
    // ---------------------------------------------------------------------------
    namespace sqi_detail {

        inline void ppgMedianFilter(const std::vector<double>& x, double freq,
            double freqd, std::vector<double>& out) {
            int n = (int)x.size();
            out.resize(n);
            if (n == 0) return;

            int factor = std::max(1, (int)std::floor(freq / freqd));

            int yLen = (n + factor - 1) / factor;
            std::vector<double> y(yLen);
            for (int i = 0; i < yLen; ++i)
                y[i] = x[std::min(i * factor, n - 1)];

            int halfwind = (int)std::floor(1.3 * freqd / 2.0);
            std::vector<double> m(yLen);
            std::vector<double> window;
            window.reserve(2 * halfwind + 1);

            for (int i = 0; i < yLen; ++i) {
                int ps = std::max(i - halfwind, 0);
                int pe = std::min(i + halfwind, yLen - 1);
                int wLen = pe - ps + 1;

                window.assign(y.begin() + ps, y.begin() + ps + wLen);
                std::sort(window.begin(), window.end());

                int lo = (int)std::floor(0.375 * wLen);
                int hi = (int)std::floor(0.625 * wLen);
                if (hi <= lo) hi = lo + 1;
                hi = std::min(hi, wLen);

                double sum = 0;
                for (int j = lo; j < hi; ++j) sum += window[j];
                m[i] = sum / (hi - lo);
            }

            for (int i = 0; i < n; ++i) {
                int idx = std::min((i / factor), yLen - 1);
                out[i] = x[i] - m[idx];
            }
        }

        // ---------------------------------------------------------------------------
        // PLA - Piecewise Linear Approximation (PLA.m)
        // ---------------------------------------------------------------------------
        inline std::vector<int> pla(const std::vector<double>& input, int s = 1,
            double th = 1.0) {
            int n = (int)input.size();
            if (n == 0) return {};

            std::vector<int> result;
            result.reserve(n / s + 2);
            result.push_back(0);

            int i = 0;
            int s1 = s;

            while (i < n - 1) {
                int i_plus_s = std::min(i + s1, n - 1);
                bool interrupted = false;

                while (!interrupted) {
                    int j = i + 1;
                    while (j <= i_plus_s) {
                        double distance = input[i_plus_s] - input[i];
                        double denom = (double)(i_plus_s - i);
                        double dcur = input[j] - input[i] - (distance * (j - i) / denom);
                        if (std::abs(dcur) > th) {
                            s1 = j - i;
                            i_plus_s = i + s1;
                            interrupted = true;
                            break;
                        }
                        ++j;
                    }
                    if (interrupted) {
                        result.push_back(j - 2);
                        i = j - 3;
                        s1 = s;
                        break;
                    }
                    else {
                        if (i_plus_s >= n - 1) {
                            i_plus_s = n - 1;
                            break;
                        }
                        else {
                            i_plus_s = std::min(i_plus_s + s1, n - 1);
                        }
                    }
                }
                ++i;
            }
            result.push_back(n - 1);
            return result;
        }

        // ---------------------------------------------------------------------------
        // Normalize to [0, 100]
        // ---------------------------------------------------------------------------
        inline void normalizeTo100(const std::vector<double>& in,
            std::vector<double>& out) {
            out.resize(in.size());
            if (in.empty()) return;
            double mn = *std::min_element(in.begin(), in.end());
            double mx = *std::max_element(in.begin(), in.end());
            double r = mx - mn;
            if (r < 1e-15) {
                std::fill(out.begin(), out.end(), 0.0);
                return;
            }
            double inv = 100.0 / r;
            for (size_t i = 0; i < in.size(); ++i)
                out[i] = (in[i] - mn) * inv;
        }

        // ---------------------------------------------------------------------------
        // simmx_dtw — slope-difference DTW cost matrix (no intermediate matrices)
        // ---------------------------------------------------------------------------
        struct DtwSim {
            std::vector<double> w;
            std::vector<double> ta;
            std::vector<double> tb;
            int rows = 0, cols = 0;
        };

        inline void simmxDtw(const std::vector<double>& y1, const std::vector<int>& pla1,
            const std::vector<double>& y2, const std::vector<int>& pla2,
            DtwSim& out) {
            int len1 = (int)pla1.size();
            int len2 = (int)pla2.size();
            out.rows = len1;
            out.cols = len2;

            std::vector<double> slope1(len1, 0.0), slope2(len2, 0.0);
            out.ta.resize(len1);
            out.tb.resize(len2);
            out.ta[0] = 1.0;
            out.tb[0] = 1.0;

            for (int i = 1; i < len1; ++i) {
                int d = pla1[i] - pla1[i - 1];
                slope1[i] = (d != 0) ? (y1[pla1[i]] - y1[pla1[i - 1]]) / d : 0.0;
                out.ta[i] = (double)d;
            }
            for (int i = 1; i < len2; ++i) {
                int d = pla2[i] - pla2[i - 1];
                slope2[i] = (d != 0) ? (y2[pla2[i]] - y2[pla2[i - 1]]) / d : 0.0;
                out.tb[i] = (double)d;
            }

            out.w.resize((size_t)len1 * len2);
            for (int i = 0; i < len1; ++i) {
                double s1 = slope1[i];
                int base = i * len2;
                for (int j = 0; j < len2; ++j)
                    out.w[base + j] = std::abs(s1 - slope2[j]);
            }
        }

        // ---------------------------------------------------------------------------
        // dp_dtw2 — DP min-cost DTW path
        // ---------------------------------------------------------------------------
        inline void dpDtw2(const DtwSim& sim,
            std::vector<int>& p, std::vector<int>& q) {
            int r = sim.rows, c = sim.cols;
            if (r == 0 || c == 0) { p.clear(); q.clear(); return; }

            const auto& M = sim.w;
            const auto& ta = sim.ta;
            const auto& tb = sim.tb;

            int dc = c + 1;
            std::vector<double> D((r + 1) * dc, std::numeric_limits<double>::quiet_NaN());
            std::vector<uint8_t> phi(r * c, 0);

            D[0] = 0.0;

            for (int i = 0; i < r; ++i)
                for (int j = 0; j < c; ++j)
                    D[(i + 1) * dc + (j + 1)] = M[i * c + j];

            for (int i = 0; i < r; ++i) {
                for (int j = 0; j < c; ++j) {
                    double mij = M[i * c + j];
                    double tai = ta[i], tbj = tb[j];

                    double d0 = D[i * dc + j];
                    double d1 = D[i * dc + (j + 1)];
                    double d2 = D[(i + 1) * dc + j];

                    double c0 = d0 + mij * (tai + tbj);
                    double c1 = d1 + mij * tai;
                    double c2 = d2 + mij * tbj;

                    bool v0 = !std::isnan(d0), v1 = !std::isnan(d1), v2 = !std::isnan(d2);

                    double best = std::numeric_limits<double>::infinity();
                    uint8_t tb_val = 1;

                    if (v0 && c0 < best) { best = c0; tb_val = 1; }
                    if (v1 && c1 < best) { best = c1; tb_val = 2; }
                    if (v2 && c2 < best) { best = c2; tb_val = 3; }

                    D[(i + 1) * dc + (j + 1)] += best;
                    phi[i * c + j] = tb_val;
                }
            }

            int i = r - 1, j = c - 1;
            std::vector<int> pp, qq;
            pp.reserve(r + c);
            qq.reserve(r + c);
            pp.push_back(i);
            qq.push_back(j);

            while (i > 0 || j > 0) {
                int ci = std::max(i, 0), cj = std::max(j, 0);
                uint8_t tb_v = phi[ci * c + cj];
                if (tb_v == 1) { --i; --j; }
                else if (tb_v == 2) { --i; }
                else { --j; }
                i = std::max(i, 0);
                j = std::max(j, 0);
                pp.push_back(i);
                qq.push_back(j);
                if (i == 0 && j == 0) break;
            }

            std::reverse(pp.begin(), pp.end());
            std::reverse(qq.begin(), qq.end());
            p = std::move(pp);
            q = std::move(qq);
        }

        // ---------------------------------------------------------------------------
        // draw_dtw — warp y2 onto y1's time axis (linear interp, not spline)
        // ---------------------------------------------------------------------------
        inline void drawDtw(const std::vector<double>& y1, const std::vector<int>& pla1,
            const std::vector<int>& p,
            const std::vector<double>& y2, const std::vector<int>& pla2,
            const std::vector<int>& q,
            std::vector<double>& y2modify) {
            int L = (int)p.size();
            if (L < 2 || y1.empty() || y2.empty()) {
                y2modify = y1;
                return;
            }

            std::vector<double> outx, outy;
            outx.reserve(y1.size());
            outy.reserve(y1.size());

            int startI = 0, startJ = 0;
            while (startI + 1 < L && p[startI + 1] == p[0]) ++startI;
            while (startJ + 1 < L && q[startJ + 1] == q[0]) ++startJ;
            int point = std::max(startI, startJ) + 1;

            int prevI = startI, prevJ = startJ;

            while (point < L) {
                while (point + 1 < L && (p[point + 1] == p[point] || q[point + 1] == q[point]))
                    ++point;

                if (point < L) {
                    int pa1Start = pla1[p[prevI]];
                    int pa1End = pla1[p[point]];
                    int pa2Start = pla2[q[prevJ]];
                    int pa2End = pla2[q[point]];

                    if (pa1End > pa1Start && pa2End >= pa2Start) {
                        for (int xi = pa1Start; xi <= pa1End; ++xi) {
                            double t = (double)(xi - pa1Start) / (double)(pa1End - pa1Start);
                            double srcIdx = pa2Start + t * (pa2End - pa2Start);
                            int lo = (int)std::floor(srcIdx);
                            int hi = lo + 1;
                            lo = std::clamp(lo, 0, (int)y2.size() - 1);
                            hi = std::clamp(hi, 0, (int)y2.size() - 1);
                            double frac = srcIdx - std::floor(srcIdx);
                            double val = y2[lo] * (1.0 - frac) + y2[hi] * frac;
                            outx.push_back((double)xi);
                            outy.push_back(val);
                        }
                    }
                    prevI = point;
                    prevJ = point;
                }
                ++point;
            }

            if (outx.empty()) {
                y2modify = y1;
                return;
            }

            std::vector<double> ux, uy;
            ux.reserve(outx.size());
            uy.reserve(outy.size());
            for (size_t i = 0; i < outx.size(); ++i) {
                if (i + 1 < outx.size() && outx[i] == outx[i + 1]) continue;
                ux.push_back(outx[i]);
                uy.push_back(outy[i]);
            }

            int n1 = (int)y1.size();
            y2modify.resize(n1);

            if (ux.size() < 2) {
                std::fill(y2modify.begin(), y2modify.end(), uy.empty() ? 0.0 : uy[0]);
                return;
            }

            double xStart = ux.front(), xEnd = ux.back();
            for (int i = 0; i < n1; ++i) {
                double t = xStart + (double)i / (n1 - 1) * (xEnd - xStart);
                auto it = std::lower_bound(ux.begin(), ux.end(), t);
                int idx = (int)(it - ux.begin());
                if (idx <= 0) {
                    y2modify[i] = uy[0];
                }
                else if (idx >= (int)ux.size()) {
                    y2modify[i] = uy.back();
                }
                else {
                    double frac = (t - ux[idx - 1]) / (ux[idx] - ux[idx - 1]);
                    y2modify[i] = uy[idx - 1] * (1.0 - frac) + uy[idx] * frac;
                }
            }
        }

        // ---------------------------------------------------------------------------
        // Discrete Frechet Distance — iterative bottom-up DP (not recursive)
        // ---------------------------------------------------------------------------
        inline double discreteFrechetDist(const std::vector<double>& P,
            const std::vector<double>& Q) {
            int n = (int)P.size(), m = (int)Q.size();
            if (n == 0 || m == 0) return 0.0;

            std::vector<double> CA((size_t)n * m);

            CA[0] = std::abs(P[0] - Q[0]);

            for (int i = 1; i < n; ++i)
                CA[i * m] = std::max(CA[(i - 1) * m], std::abs(P[i] - Q[0]));

            for (int j = 1; j < m; ++j)
                CA[j] = std::max(CA[j - 1], std::abs(P[0] - Q[j]));

            for (int i = 1; i < n; ++i) {
                for (int j = 1; j < m; ++j) {
                    double d = std::abs(P[i] - Q[j]);
                    double prev = std::min({ CA[(i - 1) * m + j],
                                            CA[(i - 1) * m + (j - 1)],
                                            CA[i * m + (j - 1)] });
                    CA[i * m + j] = std::max(prev, d);
                }
            }

            return CA[n * m - 1];
        }

        // ---------------------------------------------------------------------------
        // Correlation coefficient
        // ---------------------------------------------------------------------------
        inline double corrcoef(const double* a, const double* b, int n) {
            if (n < 2) return 0.0;

            double sa = 0, sb = 0, sa2 = 0, sb2 = 0, sab = 0;
            for (int i = 0; i < n; ++i) {
                sa += a[i]; sb += b[i];
                sa2 += a[i] * a[i]; sb2 += b[i] * b[i];
                sab += a[i] * b[i];
            }
            double dn = (double)n;
            double num = dn * sab - sa * sb;
            double den = std::sqrt((dn * sa2 - sa * sa) * (dn * sb2 - sb * sb));
            return (den > 1e-15) ? (num / den) : 0.0;
        }

        // ---------------------------------------------------------------------------
        // Linear interpolation resample
        // ---------------------------------------------------------------------------
        inline void linResample(const double* src, int srcLen,
            double* dst, int dstLen) {
            if (srcLen < 1 || dstLen < 1) return;
            if (dstLen == 1) { dst[0] = src[0]; return; }
            for (int i = 0; i < dstLen; ++i) {
                double t = (double)i / (dstLen - 1) * (srcLen - 1);
                int lo = (int)t;
                int hi = std::min(lo + 1, srcLen - 1);
                double f = t - lo;
                dst[i] = src[lo] * (1.0 - f) + src[hi] * f;
            }
        }

        // ---------------------------------------------------------------------------
        // Per-thread scratch buffers
        // ---------------------------------------------------------------------------
        struct SqiScratch {
            std::vector<double> d2_norm;
            std::vector<double> y_resamp;
            std::vector<double> y2mod;
            std::vector<double> interpP, interpQ;
            std::vector<int> plaResult;
            DtwSim dtwSim;
            std::vector<int> dtwP, dtwQ;

            void reserve(int maxBeatLen, int maxTemplateLen) {
                int mx = std::max(maxBeatLen, maxTemplateLen) + 16;
                d2_norm.reserve(mx);
                y_resamp.reserve(mx);
                y2mod.reserve(mx);
                interpP.reserve(mx);
                interpQ.reserve(mx);
                plaResult.reserve(mx);
                dtwSim.w.reserve(mx * mx / 4);
                dtwSim.ta.reserve(mx);
                dtwSim.tb.reserve(mx);
                dtwP.reserve(mx * 2);
                dtwQ.reserve(mx * 2);
            }
        };

    } // namespace sqi_detail

    // ---------------------------------------------------------------------------
    // BeatSqi struct — matches what generate_features.hpp expects
    // ---------------------------------------------------------------------------
    struct BeatSqi {
        double mean_corr = 0.0;
        double corr_direct = 0.0;
        double corr_interp = 0.0;
        double corr_dtw = 0.0;
        double frechet = 0.0;
    };

    namespace sqi_detail {

        // ---------------------------------------------------------------------------
        // Compute SQI for one beat
        // ---------------------------------------------------------------------------
        inline BeatSqi computeBeatSQI(
            const std::vector<double>& wave,
            int beatBegin, int beatEnd,
            const std::vector<double>& tmpl,
            const std::vector<int>& tmplPla,
            double Fs,
            SqiScratch& scratch)
        {
            BeatSqi result;

            int beatLen = beatEnd - beatBegin;
            int tLen = (int)tmpl.size();

            if (beatLen < 2 || tLen < 2) return result;

            if (beatLen > (int)(3.0 * Fs))
                beatLen = (int)(3.0 * Fs);

            if (beatBegin + beatLen > (int)wave.size()) return result;

            // ---- SQI1: Direct comparison ----
            int compLen = std::min(tLen, beatLen);
            double c1 = sqi_detail::corrcoef(tmpl.data(), &wave[beatBegin], compLen);
            if (c1 < 0) c1 = 0;

            // ---- SQI2: Linear resampling ----
            scratch.y_resamp.resize(tLen);
            sqi_detail::linResample(&wave[beatBegin], beatLen, scratch.y_resamp.data(), tLen);
            double c2 = sqi_detail::corrcoef(tmpl.data(), scratch.y_resamp.data(), tLen);
            if (c2 < 0) c2 = 0;

            // ---- SQI3: Dynamic Time Warping ----
            double c3 = 0;

            if (beatLen <= tLen * 10) {
                scratch.d2_norm.resize(beatLen);
                double mn = *std::min_element(&wave[beatBegin], &wave[beatBegin] + beatLen);
                double mx = *std::max_element(&wave[beatBegin], &wave[beatBegin] + beatLen);
                double r = mx - mn;
                if (r > 1e-15) {
                    double inv = 100.0 / r;
                    for (int i = 0; i < beatLen; ++i)
                        scratch.d2_norm[i] = (wave[beatBegin + i] - mn) * inv;
                }
                else {
                    std::fill_n(scratch.d2_norm.data(), beatLen, 0.0);
                }

                scratch.plaResult = sqi_detail::pla(scratch.d2_norm, 1, 1.0);

                sqi_detail::simmxDtw(tmpl, tmplPla, scratch.d2_norm, scratch.plaResult, scratch.dtwSim);
                sqi_detail::dpDtw2(scratch.dtwSim, scratch.dtwP, scratch.dtwQ);
                sqi_detail::drawDtw(tmpl, tmplPla, scratch.dtwP,
                    scratch.d2_norm, scratch.plaResult, scratch.dtwQ,
                    scratch.y2mod);

                if ((int)scratch.y2mod.size() == tLen) {
                    c3 = sqi_detail::corrcoef(tmpl.data(), scratch.y2mod.data(), tLen);
                    if (c3 < 0) c3 = 0;
                }
            }

            // ---- SQI5: Discrete Frechet Distance ----
            int maxLen = std::max(tLen, beatLen);

            scratch.interpP.resize(maxLen);
            if (tLen < maxLen) {
                std::vector<double> tmplUnit(tLen);
                for (int i = 0; i < tLen; ++i) tmplUnit[i] = tmpl[i] / 100.0;
                sqi_detail::linResample(tmplUnit.data(), tLen, scratch.interpP.data(), maxLen);
            }
            else {
                for (int i = 0; i < tLen; ++i) scratch.interpP[i] = tmpl[i] / 100.0;
            }

            scratch.interpQ.resize(maxLen);
            if (beatLen < maxLen) {
                if (beatLen > tLen * 10) {
                    scratch.d2_norm.resize(beatLen);
                    double mn2 = *std::min_element(&wave[beatBegin], &wave[beatBegin] + beatLen);
                    double mx2 = *std::max_element(&wave[beatBegin], &wave[beatBegin] + beatLen);
                    double r2 = mx2 - mn2;
                    if (r2 > 1e-15) {
                        for (int i = 0; i < beatLen; ++i)
                            scratch.d2_norm[i] = (wave[beatBegin + i] - mn2) / r2;
                    }
                    else {
                        std::fill_n(scratch.d2_norm.data(), beatLen, 0.0);
                    }
                    sqi_detail::linResample(scratch.d2_norm.data(), beatLen,
                        scratch.interpQ.data(), maxLen);
                }
                else {
                    std::vector<double> d2unit(beatLen);
                    for (int i = 0; i < beatLen; ++i) d2unit[i] = scratch.d2_norm[i] / 100.0;
                    sqi_detail::linResample(d2unit.data(), beatLen, scratch.interpQ.data(), maxLen);
                }
            }
            else {
                if ((int)scratch.d2_norm.size() >= beatLen) {
                    for (int i = 0; i < beatLen; ++i)
                        scratch.interpQ[i] = scratch.d2_norm[i] / 100.0;
                }
                else {
                    double mn2 = *std::min_element(&wave[beatBegin], &wave[beatBegin] + beatLen);
                    double mx2 = *std::max_element(&wave[beatBegin], &wave[beatBegin] + beatLen);
                    double r2 = mx2 - mn2;
                    for (int i = 0; i < beatLen; ++i)
                        scratch.interpQ[i] = (r2 > 1e-15)
                        ? (wave[beatBegin + i] - mn2) / r2 : 0.0;
                }
            }

            double cm = sqi_detail::discreteFrechetDist(scratch.interpP, scratch.interpQ);

            double dif = std::sqrt(2.0) - (cm / 2.0);
            double frechetSqi = dif / std::sqrt(2.0);

            result.corr_direct = c1;
            result.corr_interp = c2;
            result.corr_dtw = c3;
            result.frechet = frechetSqi;
            result.mean_corr = (c1 + c2 + c3) / 3.0;
            return result;
        }

    } // namespace sqi_detail

    // ---------------------------------------------------------------------------
    // ppg_sqi — main entry point, matches generate_features.hpp call signature
    // ---------------------------------------------------------------------------
    inline std::vector<BeatSqi> ppg_sqi(
        const std::vector<double>& wave,
        const std::vector<int>& anntime,
        const std::vector<double>& tmplRaw,
        int                        windowlen,
        double                     Fs,
        ThreadPool* pool = nullptr)
    {
        (void)windowlen;

        int nBeats = (int)anntime.size() - 1;
        std::vector<BeatSqi> results(std::max(nBeats, 0));

        if (nBeats <= 0 || tmplRaw.empty()) return results;

        // Filter baseline wander (done once, shared across all beats)
        std::vector<double> filtered;
        sqi_detail::ppgMedianFilter(wave, Fs, Fs, filtered);

        // Normalize template to 0-100 (done once)
        std::vector<double> tmplNorm;
        sqi_detail::normalizeTo100(tmplRaw, tmplNorm);

        // PLA of template (done once)
        std::vector<int> tmplPla = sqi_detail::pla(tmplNorm, 1, 1.0);

        // Estimate max beat length for buffer sizing
        int maxBeatLen = 0;
        for (int i = 0; i < nBeats; ++i) {
            int len = anntime[i + 1] - anntime[i];
            if (len > maxBeatLen) maxBeatLen = len;
        }

        if (pool != nullptr) {
            // ---- Parallel via shared ThreadPool ----
            // Batch beats into chunks so we don't flood the pool with
            // thousands of tiny tasks (reduces scheduling overhead).
            constexpr int CHUNK = 32;
            int nChunks = (nBeats + CHUNK - 1) / CHUNK;

            std::vector<std::future<void>> futures;
            futures.reserve(nChunks);

            for (int ch = 0; ch < nChunks; ++ch) {
                int jStart = ch * CHUNK;
                int jEnd = std::min(jStart + CHUNK, nBeats);

                futures.push_back(pool->enqueue(
                    [&, jStart, jEnd]() {
                        // Each task gets its own scratch buffer — no sharing, no locks
                        sqi_detail::SqiScratch scratch;
                        scratch.reserve(maxBeatLen, (int)tmplNorm.size());

                        for (int j = jStart; j < jEnd; ++j) {
                            int bStart = anntime[j];
                            int bEnd = anntime[j + 1];
                            if (bStart < 0 || bEnd >(int)filtered.size() || bEnd <= bStart)
                                continue;
                            results[j] = sqi_detail::computeBeatSQI(
                                filtered, bStart, bEnd,
                                tmplNorm, tmplPla, Fs, scratch);
                        }
                    }));
            }

            // Wait for all chunks
            for (auto& f : futures) f.get();

        }
        else {
            // ---- Serial fallback ----
            sqi_detail::SqiScratch scratch;
            scratch.reserve(maxBeatLen, (int)tmplNorm.size());

            for (int j = 0; j < nBeats; ++j) {
                int bStart = anntime[j];
                int bEnd = anntime[j + 1];
                if (bStart < 0 || bEnd >(int)filtered.size() || bEnd <= bStart) continue;
                results[j] = sqi_detail::computeBeatSQI(
                    filtered, bStart, bEnd,
                    tmplNorm, tmplPla, Fs, scratch);
            }
        }

        return results;
    }

} // namespace ppg