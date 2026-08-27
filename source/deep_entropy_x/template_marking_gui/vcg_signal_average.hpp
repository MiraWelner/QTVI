#pragma once
//
// vcg_signal_average.hpp
//
// Signal-averages the VCG per bin and extracts the loop features, then writes
// them to <dir>/<id>_vcg.csv beside the template markings CSV.
//
// ---------------------------------------------------------------------------
// ORDER OF OPERATIONS: transform per beat, THEN aggregate
// ---------------------------------------------------------------------------
// The per-lead templates are a column-wise NaN-skipping MEDIAN over the beats
// that survived alignment's rejections. The transform is linear, so it
// commutes with a mean -- but not with a median. transform(median of beats) !=
// median of transform(beats), and the two differ most exactly where it
// matters: on the asymmetric, high-slew QRS peak that dominates every feature
// below.
//
// So each beat is transformed to XYZ first, and the median is taken over the
// per-beat VCGs, column by column, exactly as create_ecg_templates does for a
// single lead. That also leaves the per-beat loops in hand, which is the only
// way to compute beat-to-beat variability at all.
//
// Same R-peak alignment and outlier rejection as the per-lead templates: this
// consumes the KEPT beat sets (create_ecg_templates' out_kept_beats), already
// on the shared axis with every beat's R at r_col. Two different things are
// going on there and the difference matters:
//
//   Tukey RR-length passes, wave-score rejections, and baseline_source==NONE
//   are applied UPSTREAM of the capture. Those beats are simply not in the set
//   and nothing here re-applies the rules -- re-deriving them would let the two
//   definitions drift.
//
//   Ectopic beats ARE in the set. create_ecg_templates captures out_kept_beats
//   BEFORE the ectopic mask on purpose, so the per-beat record keeps the
//   ectopy and flags it in the parallel kept_rhythm vector, even though the
//   median the operator marks excludes it. So supplying BinBeats::rhythm is
//   what makes the VCG match the templates; without it the loops would include
//   PVCs that the per-lead templates do not.
//
// ---------------------------------------------------------------------------
// AXIS
// ---------------------------------------------------------------------------
// Each channel's kept beats sit on that channel's own axis with R at its own
// r_col, so beats are combined on the shared R-RELATIVE axis: offset k is read
// at (that channel's r_col + k). Interval boundaries from
// global_intervals::GlobalIntervals are already R-relative and index straight
// into it.
//
// Units: amplitudes are whatever the beats carry (mV if raw). Areas are
// mV^2 per projection plane; velocity is mV/s; angles are degrees; the
// planarity and variability ratios are dimensionless.
//

#include "vcg.hpp"
#include "global_intervals.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace vcg_avg {

    inline constexpr int kNumEcgCh = 3;
    inline const double kNaN = std::numeric_limits<double>::quiet_NaN();

    /// Rhythm verdict values, matching create_ecg_templates' kept_rhythm.
    enum : std::uint8_t { RHYTHM_NORMAL = 0, RHYTHM_PVC = 1, RHYTHM_VOTED_PVC = 2 };

    /// The kept, aligned beats for one bin: three channels, each a set of
    /// equal-length beats on that channel's shared axis.
    ///
    /// EXCLUSIONS -- what is already applied and what this struct still has to
    /// do:
    ///   Tukey / wave-score / baseline-NONE  -- ALREADY GONE. out_kept_beats is
    ///       captured after those rejections, so nothing here re-applies them.
    ///   Ectopic (PVC / voted PVC)           -- STILL PRESENT. The capture in
    ///       create_ecg_templates happens BEFORE the ectopic mask, precisely so
    ///       the per-beat output keeps the ectopy and flags it. The median the
    ///       operator marks excludes it; this beat set does not. Supply
    ///       `rhythm` and the loops will.
    struct BinBeats {
        const std::vector<std::vector<double>>* beats[kNumEcgCh] = { nullptr,nullptr,nullptr };

        /// Per-kept-beat verdict, parallel to beats[c] (create_ecg_templates'
        /// kept_rhythm_raw / kept_beats_by_channel's rhythm vector). Null means
        /// "no verdicts supplied": see excludeEctopic.
        const std::vector<std::uint8_t>* rhythm[kNumEcgCh] = { nullptr,nullptr,nullptr };

        int r_col[kNumEcgCh] = { -1,-1,-1 };

        /// Drop any beat flagged non-NORMAL in ANY channel. A PVC's loop
        /// differs morphologically by definition, so leaving even a few in
        /// dominates loop variability and drags the averaged loop off the
        /// sinus morphology every other feature is measured against.
        bool excludeEctopic = true;
    };

    /// One beat's loop, on the R-relative axis.
    struct Loop {
        std::vector<vcg::VCGSample> pts;
        int firstOffset = 0;   ///< R-relative offset of pts[0]
        int indexForOffset(double off) const {
            const int i = static_cast<int>(std::lround(off)) - firstOffset;
            return (i >= 0 && i < static_cast<int>(pts.size())) ? i : -1;
        }
    };

    // -----------------------------------------------------------------------
    // Build per-beat loops and the signal-averaged loop
    // -----------------------------------------------------------------------
    /**
     * @brief Transform every kept beat to XYZ on the shared R-relative axis.
     *
     *        Beat j of channel 0 is paired with beat j of channels 1 and 2:
     *        under this pipeline all three channels are sliced from the same
     *        R-pair list, so index j is the same cardiac cycle in each. If the
     *        three sets differ in size that assumption is broken, and this
     *        returns empty rather than pairing beats that are not the same
     *        beat.
     */
    inline std::vector<Loop> perBeatLoops(const BinBeats& in, int pre, int post,
        const vcg::VcgMatrix& mat = vcg::kIdentity,
        int* outNumEctopicExcluded = nullptr) {
        std::vector<Loop> out;
        if (outNumEctopicExcluded) *outNumEctopicExcluded = 0;
        for (int c = 0; c < kNumEcgCh; ++c)
            if (!in.beats[c] || in.beats[c]->empty() || in.r_col[c] < 0) return out;

        const std::size_t nBeats = in.beats[0]->size();
        for (int c = 1; c < kNumEcgCh; ++c)
            if (in.beats[c]->size() != nBeats) return out;   // not the same beats
        if (pre < 0 || post < 0 || pre + post < 1) return out;

        // Ectopic exclusion, taken as the UNION across channels: beat j is one
        // cardiac cycle seen in three leads, so if any lead called it ectopic
        // the cycle is ectopic and all three of its lanes go. Excluding it in
        // one channel only would leave the remaining lanes contributing a
        // vector with a missing component.
        std::vector<bool> drop(nBeats, false);
        std::size_t nDropped = 0;
        if (in.excludeEctopic) {
            for (int c = 0; c < kNumEcgCh; ++c) {
                if (!in.rhythm[c]) continue;
                // A verdict vector that does not line up with the beats cannot
                // be trusted to name the right beats, so it is ignored rather
                // than applied to whichever beats happen to be in range.
                if (in.rhythm[c]->size() != nBeats) continue;
                for (std::size_t j = 0; j < nBeats; ++j)
                    if ((*in.rhythm[c])[j] != RHYTHM_NORMAL) drop[j] = true;
            }
            for (std::size_t j = 0; j < nBeats; ++j) if (drop[j]) ++nDropped;
        }
        if (outNumEctopicExcluded) *outNumEctopicExcluded = static_cast<int>(nDropped);

        const int n = pre + post + 1;
        out.reserve(nBeats);
        for (std::size_t j = 0; j < nBeats; ++j) {
            if (drop[j]) continue;
            std::vector<double> lane[kNumEcgCh];
            for (int c = 0; c < kNumEcgCh; ++c) {
                const std::vector<double>& b = (*in.beats[c])[j];
                lane[c].assign(n, kNaN);
                for (int i = 0; i < n; ++i) {
                    const int col = in.r_col[c] + (i - pre);
                    if (col >= 0 && col < static_cast<int>(b.size()))
                        lane[c][i] = b[col];
                }
            }
            const std::vector<const std::vector<double>*> lanes{ &lane[0], &lane[1], &lane[2] };
            const vcg::VcgResult r = vcg::reconstructVCG(lanes, mat);
            if (!r.valid) continue;
            Loop L; L.pts = r.samples; L.firstOffset = -pre;
            out.push_back(std::move(L));
        }
        return out;
    }

    /// Column-wise NaN-skipping median of the per-beat loops -- the same
    /// aggregation create_ecg_templates uses for one lead, applied to X, Y and
    /// Z independently.
    inline Loop medianLoop(const std::vector<Loop>& loops) {
        Loop out;
        if (loops.empty()) return out;
        const std::size_t n = loops[0].pts.size();
        out.firstOffset = loops[0].firstOffset;
        out.pts.assign(n, vcg::VCGSample{ kNaN, kNaN, kNaN });

        std::vector<double> col;
        col.reserve(loops.size());
        auto med = [&col](void) -> double {
            if (col.empty()) return kNaN;
            std::sort(col.begin(), col.end());
            const std::size_t m = col.size();
            return (m % 2 == 0) ? 0.5 * (col[m / 2 - 1] + col[m / 2]) : col[m / 2];
            };

        for (std::size_t i = 0; i < n; ++i) {
            for (int axis = 0; axis < 3; ++axis) {
                col.clear();
                for (const Loop& L : loops) {
                    if (i >= L.pts.size()) continue;
                    const vcg::VCGSample& s = L.pts[i];
                    const double v = (axis == 0) ? s.x : (axis == 1) ? s.y : s.z;
                    if (!std::isnan(v)) col.push_back(v);
                }
                const double m = med();
                if (axis == 0) out.pts[i].x = m;
                else if (axis == 1) out.pts[i].y = m;
                else               out.pts[i].z = m;
            }
        }
        return out;
    }

    // -----------------------------------------------------------------------
    // Features
    // -----------------------------------------------------------------------
    enum class Plane { XY, XZ, YZ };

    /**
     * @brief Signed area enclosed by the loop's projection, by the shoelace
     *        formula over [lo, hi].
     *
     *        The polygon is closed from the last point back to the first, as
     *        the shoelace formula requires; a QRS loop does not return exactly
     *        to its origin, so leaving it open would omit that closing
     *        triangle. Sign carries the rotation direction (positive =
     *        counter-clockwise in the plane's axis order), which is
     *        diagnostically meaningful, so it is NOT taken as absolute here.
     *
     * @return NaN if fewer than 3 usable points.
     */
    inline double shoelaceArea(const Loop& L, int lo, int hi, Plane pl) {
        if (L.pts.empty()) return kNaN;
        lo = std::max(0, lo);
        hi = std::min(hi, static_cast<int>(L.pts.size()) - 1);
        if (hi - lo + 1 < 3) return kNaN;

        auto proj = [pl](const vcg::VCGSample& s, double& a, double& b) {
            switch (pl) {
            case Plane::XY: a = s.x; b = s.y; break;
            case Plane::XZ: a = s.x; b = s.z; break;
            case Plane::YZ: a = s.y; b = s.z; break;
            }
            };

        std::vector<std::pair<double, double>> poly;
        poly.reserve(hi - lo + 1);
        for (int i = lo; i <= hi; ++i) {
            const vcg::VCGSample& s = L.pts[i];
            if (std::isnan(s.x) || std::isnan(s.y) || std::isnan(s.z)) continue;
            double a, b; proj(s, a, b);
            poly.emplace_back(a, b);
        }
        if (poly.size() < 3) return kNaN;

        double acc = 0.0;
        for (std::size_t i = 0; i < poly.size(); ++i) {
            const auto& p = poly[i];
            const auto& q = poly[(i + 1) % poly.size()];   // closes the polygon
            acc += p.first * q.second - q.first * p.second;
        }
        return 0.5 * acc;
    }

    /**
     * @brief Eigenvalues of the loop points' scatter matrix, descending.
     *        This is the SVD of the mean-centred loop expressed through its
     *        normal matrix; for a 3x3 a Jacobi rotation is exact enough and
     *        avoids a linear-algebra dependency.
     */
    struct LoopSpread { double lambda[3] = { kNaN,kNaN,kNaN }; int nPts = 0; bool valid = false; };

    inline LoopSpread loopSpread(const Loop& L, int lo, int hi) {
        LoopSpread out;
        if (L.pts.empty()) return out;
        lo = std::max(0, lo);
        hi = std::min(hi, static_cast<int>(L.pts.size()) - 1);

        double mean[3] = { 0,0,0 }; int cnt = 0;
        for (int i = lo; i <= hi; ++i) {
            const vcg::VCGSample& s = L.pts[i];
            if (std::isnan(s.x) || std::isnan(s.y) || std::isnan(s.z)) continue;
            mean[0] += s.x; mean[1] += s.y; mean[2] += s.z; ++cnt;
        }
        if (cnt < 3) return out;
        for (int k = 0; k < 3; ++k) mean[k] /= cnt;

        double C[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
        for (int i = lo; i <= hi; ++i) {
            const vcg::VCGSample& s = L.pts[i];
            if (std::isnan(s.x) || std::isnan(s.y) || std::isnan(s.z)) continue;
            const double d[3] = { s.x - mean[0], s.y - mean[1], s.z - mean[2] };
            for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) C[r][c] += d[r] * d[c];
        }
        for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) C[r][c] /= (cnt - 1);

        for (int sweep = 0; sweep < 32; ++sweep) {
            const double off = std::fabs(C[0][1]) + std::fabs(C[0][2]) + std::fabs(C[1][2]);
            if (off < 1e-15) break;
            for (int p = 0; p < 2; ++p) for (int q = p + 1; q < 3; ++q) {
                if (std::fabs(C[p][q]) < 1e-18) continue;
                const double th = 0.5 * std::atan2(2.0 * C[p][q], C[q][q] - C[p][p]);
                const double cs = std::cos(th), sn = std::sin(th);
                for (int k = 0; k < 3; ++k) {
                    const double akp = C[k][p], akq = C[k][q];
                    C[k][p] = cs * akp - sn * akq; C[k][q] = sn * akp + cs * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    const double apk = C[p][k], aqk = C[q][k];
                    C[p][k] = cs * apk - sn * aqk; C[q][k] = sn * apk + cs * aqk;
                }
            }
        }
        // A covariance matrix is positive semi-definite, so a negative
        // eigenvalue here is Jacobi round-off on a degenerate (flat or
        // collinear) loop -- around -1e-18 in practice. Clamped at zero:
        // otherwise the planarity ratio comes out slightly negative, which
        // reads as a nonsensical measurement rather than "perfectly planar".
        double ev[3] = { std::max(0.0, C[0][0]),
                         std::max(0.0, C[1][1]),
                         std::max(0.0, C[2][2]) };
        std::sort(ev, ev + 3, std::greater<double>());
        for (int k = 0; k < 3; ++k) out.lambda[k] = ev[k];
        out.nPts = cnt;
        out.valid = true;
        return out;
    }

    /**
     * @brief Peak spatial velocity: max ||d(XYZ)/dt|| over [lo, hi], in
     *        amplitude units per second.
     *
     *        Central differences where both neighbours exist, so the estimate
     *        is not biased half a sample forward the way a forward difference
     *        is -- on a QRS upstroke that shift lands the peak on the wrong
     *        sample. Any triple touching a NaN is skipped rather than treated
     *        as a step, which would otherwise register as a huge false peak.
     */
    inline double peakSpatialVelocity(const Loop& L, int lo, int hi, double fs) {
        if (L.pts.empty() || !(fs > 0.0)) return kNaN;
        lo = std::max(1, lo);
        hi = std::min(hi, static_cast<int>(L.pts.size()) - 2);
        double best = kNaN;
        for (int i = lo; i <= hi; ++i) {
            const vcg::VCGSample& a = L.pts[i - 1];
            const vcg::VCGSample& b = L.pts[i + 1];
            if (std::isnan(a.x) || std::isnan(a.y) || std::isnan(a.z) ||
                std::isnan(b.x) || std::isnan(b.y) || std::isnan(b.z)) continue;
            const double dx = (b.x - a.x) * 0.5 * fs;
            const double dy = (b.y - a.y) * 0.5 * fs;
            const double dz = (b.z - a.z) * 0.5 * fs;
            const double v = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (std::isnan(best) || v > best) best = v;
        }
        return best;
    }

    /**
     * @brief Beat-to-beat loop variability: RMS 3-D distance from each beat's
     *        loop to the signal-averaged loop over [lo, hi], averaged across
     *        beats, then divided by the averaged loop's peak spatial magnitude.
     *
     *        Normalized so it is comparable between bins and subjects: the raw
     *        RMS scales with signal amplitude, which would make a
     *        high-amplitude lead set look unstable. Dimensionless; 0 = every
     *        beat identical to the average.
     *
     * @return NaN with fewer than 2 beats (variability is undefined for one).
     */
    inline double loopVariability(const std::vector<Loop>& loops, const Loop& avg,
        int lo, int hi) {
        if (loops.size() < 2 || avg.pts.empty()) return kNaN;
        lo = std::max(0, lo);
        hi = std::min(hi, static_cast<int>(avg.pts.size()) - 1);
        if (hi < lo) return kNaN;

        double scale = 0.0;
        for (int i = lo; i <= hi; ++i) {
            const vcg::VCGSample& s = avg.pts[i];
            if (std::isnan(s.x) || std::isnan(s.y) || std::isnan(s.z)) continue;
            scale = std::max(scale, std::sqrt(s.x * s.x + s.y * s.y + s.z * s.z));
        }
        if (!(scale > 0.0)) return kNaN;

        double sumRms = 0.0; int nUsed = 0;
        for (const Loop& L : loops) {
            double acc = 0.0; int cnt = 0;
            for (int i = lo; i <= hi; ++i) {
                if (i >= static_cast<int>(L.pts.size())) break;
                const vcg::VCGSample& p = L.pts[i];
                const vcg::VCGSample& q = avg.pts[i];
                if (std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.z) ||
                    std::isnan(q.x) || std::isnan(q.y) || std::isnan(q.z)) continue;
                const double dx = p.x - q.x, dy = p.y - q.y, dz = p.z - q.z;
                acc += dx * dx + dy * dy + dz * dz;
                ++cnt;
            }
            if (cnt == 0) continue;
            sumRms += std::sqrt(acc / cnt);
            ++nUsed;
        }
        if (nUsed == 0) return kNaN;
        return (sumRms / nUsed) / scale;
    }

    // -----------------------------------------------------------------------
    // One bin's full feature row
    // -----------------------------------------------------------------------
    struct BinFeatures {
        int    binIndex = -1;
        int    nBeats = 0;         ///< beats that produced a loop (post-exclusion)
        int    nEctopicExcluded = 0;  ///< beats dropped as PVC / voted PVC
        double qrstAngle_deg = kNaN;
        double qrsArea_xy = kNaN, qrsArea_xz = kNaN, qrsArea_yz = kNaN;
        double tArea_xy = kNaN, tArea_xz = kNaN, tArea_yz = kNaN;
        double qrsLambda1 = kNaN, qrsLambda2 = kNaN, qrsLambda3 = kNaN;
        double planarity = kNaN;   ///< lambda3 / (l1+l2+l3): OUT-of-plane fraction.
        ///< 0 = perfectly planar loop. Stated as a
        ///< fraction of total spread so it does not
        ///< blow up when lambda1 is small.
        double planarity_l2_l1 = kNaN;  ///< lambda2/lambda1: loop roundness vs. a line
        double peakSpatialVelocity = kNaN;
        double loopVariability = kNaN;
        double qrsDuration_ms = kNaN;
        bool   valid = false;
        std::string note;          ///< why not, when !valid
    };

    /**
     * @brief Everything for one bin.
     *
     * @param g    Global intervals for this bin: supplies the QRS window and,
     *             via qtInterval, the T end. All R-relative, so they index the
     *             loops directly.
     * @param fs   ECG sample rate, for velocity and the ms columns.
     */
    inline BinFeatures analyzeBin(int binIndex, const BinBeats& in,
        const global_intervals::GlobalIntervals& g,
        double fs, int marginSamples = 10,
        const vcg::VcgMatrix& mat = vcg::kIdentity) {
        BinFeatures f;
        f.binIndex = binIndex;

        if (!g.valid) { f.note = "global intervals not established"; return f; }

        // Window: the QRS plus the T wave, both R-relative. tEnd comes from the
        // QT interval, which is measured from the global QRS onset.
        const double onset = g.qrsOnset, offset = g.qrsOffset;
        double tEndOff = kNaN;
        if (!std::isnan(g.qtInterval_ms) && fs > 0.0)
            tEndOff = onset + g.qtInterval_ms * fs / 1000.0;

        const int pre = static_cast<int>(std::ceil(-onset)) + marginSamples;
        const int post = static_cast<int>(std::ceil(std::isnan(tEndOff) ? offset : tEndOff))
            + marginSamples;

        int nEctopic = 0;
        const std::vector<Loop> loops = perBeatLoops(in, pre, post, mat, &nEctopic);
        f.nEctopicExcluded = nEctopic;
        if (loops.empty()) {
            // Distinguish "the beats were all ectopic" from "there were no
            // beats": the first is a rhythm finding about this bin, the second
            // is a plumbing fault, and they need different follow-up.
            f.note = (nEctopic > 0)
                ? ("all " + std::to_string(nEctopic) + " beats excluded as ectopic")
                : "no usable beats (channel missing or beat counts differ)";
            return f;
        }
        const Loop avg = medianLoop(loops);
        f.nBeats = static_cast<int>(loops.size());

        const int qLo = avg.indexForOffset(onset);
        const int qHi = avg.indexForOffset(offset);
        if (qLo < 0 || qHi <= qLo) { f.note = "QRS window outside the loop"; return f; }

        f.qrsArea_xy = shoelaceArea(avg, qLo, qHi, Plane::XY);
        f.qrsArea_xz = shoelaceArea(avg, qLo, qHi, Plane::XZ);
        f.qrsArea_yz = shoelaceArea(avg, qLo, qHi, Plane::YZ);

        const LoopSpread sp = loopSpread(avg, qLo, qHi);
        if (sp.valid) {
            f.qrsLambda1 = sp.lambda[0]; f.qrsLambda2 = sp.lambda[1]; f.qrsLambda3 = sp.lambda[2];
            const double tot = sp.lambda[0] + sp.lambda[1] + sp.lambda[2];
            if (tot > 0.0) f.planarity = sp.lambda[2] / tot;
            if (sp.lambda[0] > 0.0) f.planarity_l2_l1 = sp.lambda[1] / sp.lambda[0];
        }

        f.peakSpatialVelocity = peakSpatialVelocity(avg, qLo, qHi, fs);
        f.loopVariability = loopVariability(loops, avg, qLo, qHi);
        if (fs > 0.0) f.qrsDuration_ms = (offset - onset) * 1000.0 / fs;

        // T loop and the QRS-T angle need the T window. Left NaN rather than
        // guessed when QT was not measured -- a T loop over an assumed window
        // is not a T loop.
        if (!std::isnan(tEndOff)) {
            const int tHi = avg.indexForOffset(tEndOff);
            if (tHi > qHi) {
                f.tArea_xy = shoelaceArea(avg, qHi, tHi, Plane::XY);
                f.tArea_xz = shoelaceArea(avg, qHi, tHi, Plane::XZ);
                f.tArea_yz = shoelaceArea(avg, qHi, tHi, Plane::YZ);
                f.qrstAngle_deg = vcg::spatialQRSTAngle(avg.pts, qLo, qHi, tHi);
            }
        }

        f.valid = true;
        return f;
    }


    // -----------------------------------------------------------------------
    // 6. Save-time path: build the loop from the stored TEMPLATES
    // -----------------------------------------------------------------------
    //
    // The per-beat path above needs the kept beats, which exist only inside
    // make_averaged_templates. At save time all that survives on a TemplateBin
    // is the three channel TEMPLATES -- so this path reconstructs the loop
    // from those instead, and needs no change to the .bin format.
    //
    // The cost of that, stated plainly: a template is a column-wise MEDIAN of
    // beats, and the transform is linear, so
    //   transform(median of beats) != median of transform(beats).
    // This path computes the former. The difference is largest on the
    // high-slew QRS peak. For angle, areas, planarity and velocity it is a
    // small bias, and it is IDENTICAL for every bin and every subject, so
    // comparisons across bins stay valid.
    //
    // loopVariability is the exception -- it is a statement about the scatter
    // BETWEEN beats and cannot be recovered from their median at all. It stays
    // NaN on this path. To get it, compute it at build time (where the beats
    // are) and carry one double per bin.
    //

    /// Build the loop from a bin's three channel templates, on the shared
    /// R-relative axis (each channel read at its OWN r_col + offset).
    inline Loop loopFromTemplates(const TemplateBin& b, int pre, int post,
        const vcg::VcgMatrix& mat = vcg::kIdentity) {
        Loop out;
        out.firstOffset = -pre;
        if (pre < 0 || post < 0 || pre + post < 1) return out;

        const std::vector<double>* tpl[kNumEcgCh] = {
            &b.ch1.ecgTemplate_raw, &b.ch2.ecgTemplate_raw, &b.ch3.ecgTemplate_raw };
        double rc[kNumEcgCh];
        for (int c = 0; c < kNumEcgCh; ++c) {
            rc[c] = (b.r_peak_auto_ch[c] >= 0.0) ? b.r_peak_auto_ch[c]
                : static_cast<double>(b.r_peak_ch[c]);
            // All three channels are required: a loop from two axes is a
            // projection that looks plausible and is wrong.
            if (rc[c] < 0.0 || tpl[c]->size() < 3) return out;
        }

        const int n = pre + post + 1;
        std::vector<double> lane[kNumEcgCh];
        for (int c = 0; c < kNumEcgCh; ++c) {
            lane[c].assign(n, kNaN);
            for (int i = 0; i < n; ++i)
                lane[c][i] = FeatureMarks::sample_at(*tpl[c], rc[c] + (i - pre));
        }
        const std::vector<const std::vector<double>*> lanes{ &lane[0], &lane[1], &lane[2] };
        const vcg::VcgResult r = vcg::reconstructVCG(lanes, mat);
        if (!r.valid) return out;
        out.pts = r.samples;
        return out;
    }


    /**
     * @brief The derived VCG scalar laid out on channel `refCh`'s COLUMN axis,
     *        so it plots against exactly the same x axis as that channel's
     *        panel and stacks under it correctly.
     *
     *        Index j of the result is column j of refCh's template. Every
     *        channel is still read at ITS OWN r_col plus the same R-relative
     *        offset (j - refCh's r_col), so the combination is per-instant even
     *        though the output is indexed in one channel's columns. Columns
     *        where any channel has no data come back NaN, which the plot skips.
     *
     * @param outRCol  Receives refCh's R column, for BinPlotWidget's
     *                 rPeakSample argument. Optional.
     * @return Empty when a channel or an R column is missing.
     */
    inline std::vector<double> derivedTraceOnChannelAxis(const TemplateBin& b,
        int refCh, vcg::DerivedLead which,
        const vcg::VcgMatrix& mat = vcg::kIdentity,
        double* outRCol = nullptr) {
        std::vector<double> out;
        if (refCh < 0 || refCh >= kNumEcgCh) return out;

        const std::vector<double>* tpl[kNumEcgCh] = {
            &b.ch1.ecgTemplate_raw, &b.ch2.ecgTemplate_raw, &b.ch3.ecgTemplate_raw };
        double rc[kNumEcgCh];
        for (int c = 0; c < kNumEcgCh; ++c) {
            rc[c] = (b.r_peak_auto_ch[c] >= 0.0) ? b.r_peak_auto_ch[c]
                : static_cast<double>(b.r_peak_ch[c]);
            if (rc[c] < 0.0 || tpl[c]->size() < 3) return out;
        }
        if (outRCol) *outRCol = rc[refCh];

        const int n = static_cast<int>(tpl[refCh]->size());
        std::vector<double> lane[kNumEcgCh];
        for (int c = 0; c < kNumEcgCh; ++c) {
            lane[c].assign(n, kNaN);
            for (int j = 0; j < n; ++j)
                lane[c][j] = FeatureMarks::sample_at(*tpl[c], rc[c] + (j - rc[refCh]));
        }
        const std::vector<const std::vector<double>*> lanes{ &lane[0], &lane[1], &lane[2] };
        const vcg::VcgResult r = vcg::reconstructVCG(lanes, mat);
        if (!r.valid) return out;
        return vcg::derivedLeadTrace(r, which);
    }

    /// Extract every interval-dependent feature from an already-built loop.
    /// loopVariability is left NaN; pass it separately if it was computed at
    /// build time.
    inline BinFeatures extractFeatures(int binIndex, const Loop& loop,
        const global_intervals::GlobalIntervals& g,
        double fs) {
        BinFeatures f;
        f.binIndex = binIndex;
        if (loop.pts.empty()) { f.note = "no VCG loop (channel or R column missing)"; return f; }
        if (!g.valid) { f.note = "global intervals not established"; return f; }

        const int qLo = loop.indexForOffset(g.qrsOnset);
        const int qHi = loop.indexForOffset(g.qrsOffset);
        if (qLo < 0 || qHi <= qLo) { f.note = "QRS window outside the loop"; return f; }

        f.qrsArea_xy = shoelaceArea(loop, qLo, qHi, Plane::XY);
        f.qrsArea_xz = shoelaceArea(loop, qLo, qHi, Plane::XZ);
        f.qrsArea_yz = shoelaceArea(loop, qLo, qHi, Plane::YZ);

        const LoopSpread sp = loopSpread(loop, qLo, qHi);
        if (sp.valid) {
            f.qrsLambda1 = sp.lambda[0]; f.qrsLambda2 = sp.lambda[1]; f.qrsLambda3 = sp.lambda[2];
            const double tot = sp.lambda[0] + sp.lambda[1] + sp.lambda[2];
            if (tot > 0.0)          f.planarity = sp.lambda[2] / tot;
            if (sp.lambda[0] > 0.0) f.planarity_l2_l1 = sp.lambda[1] / sp.lambda[0];
        }

        f.peakSpatialVelocity = peakSpatialVelocity(loop, qLo, qHi, fs);
        if (fs > 0.0) f.qrsDuration_ms = (g.qrsOffset - g.qrsOnset) * 1000.0 / fs;

        // T loop and QRS-T angle need a measured T end. Left NaN rather than
        // guessed: a T loop over an assumed window is not a T loop.
        if (!std::isnan(g.qtInterval_ms) && fs > 0.0) {
            const double tEndOff = g.qrsOnset + g.qtInterval_ms * fs / 1000.0;
            const int tHi = loop.indexForOffset(tEndOff);
            if (tHi > qHi) {
                f.tArea_xy = shoelaceArea(loop, qHi, tHi, Plane::XY);
                f.tArea_xz = shoelaceArea(loop, qHi, tHi, Plane::XZ);
                f.tArea_yz = shoelaceArea(loop, qHi, tHi, Plane::YZ);
                f.qrstAngle_deg = vcg::spatialQRSTAngle(loop.pts, qLo, qHi, tHi);
            }
        }
        f.valid = true;
        return f;
    }

    /// One call for the save path: template loop -> features. `marginSamples`
    /// widens the window past the measured QT so the T loop cannot be clipped.
    inline BinFeatures analyzeBinFromTemplates(int binIndex, const TemplateBin& b,
        const global_intervals::GlobalIntervals& g,
        double fs, int marginSamples = 10,
        const vcg::VcgMatrix& mat = vcg::kIdentity) {
        int pre = 0, post = 0;
        if (g.valid) {
            pre = static_cast<int>(std::ceil(-g.qrsOnset)) + marginSamples;
            post = static_cast<int>(std::ceil(g.qrsOffset)) + marginSamples;
            if (!std::isnan(g.qtInterval_ms) && fs > 0.0) {
                const double tEndOff = g.qrsOnset + g.qtInterval_ms * fs / 1000.0;
                post = std::max(post, static_cast<int>(std::ceil(tEndOff)) + marginSamples);
            }
        }
        if (pre <= 0)  pre = 40 + marginSamples;
        if (post <= 0) post = 60 + marginSamples;
        return extractFeatures(binIndex, loopFromTemplates(b, pre, post, mat), g, fs);
    }

    // -----------------------------------------------------------------------
    // CSV
    // -----------------------------------------------------------------------
    inline const char* csvHeader() {
        return "subject_id,bin_index,n_beats,n_ectopic_excluded,"
            "qrst_angle_deg,"
            "qrs_area_xy,qrs_area_xz,qrs_area_yz,"
            "t_area_xy,t_area_xz,t_area_yz,"
            "qrs_lambda1,qrs_lambda2,qrs_lambda3,"
            "planarity_out_of_plane_fraction,planarity_l2_over_l1,"
            "peak_spatial_velocity,loop_variability,qrs_duration_ms,note";
    }

    /// Empty cell for a NaN, so a missing value is never read as 0 by whatever
    /// loads this next.
    inline std::string num(double v) {
        if (std::isnan(v)) return std::string();
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.6g", v);
        return std::string(buf);
    }

    /**
     * @brief Write one row per bin to <dir>/<id>_vcg.csv.
     *
     * @param dir  Same directory the template markings CSV goes to.
     * @return false if the file could not be opened.
     */
    inline bool writeVcgCsv(const std::string& dir,
        const std::string& subjectId,
        const std::vector<BinFeatures>& rows) {
        const std::string path = dir + "/" + subjectId + "_vcg.csv";
        std::ofstream f(path);
        if (!f) return false;
        f << csvHeader() << "\n";
        for (const BinFeatures& r : rows) {
            f << subjectId << ',' << r.binIndex << ',' << r.nBeats << ','
                << r.nEctopicExcluded << ','
                << num(r.qrstAngle_deg) << ','
                << num(r.qrsArea_xy) << ',' << num(r.qrsArea_xz) << ',' << num(r.qrsArea_yz) << ','
                << num(r.tArea_xy) << ',' << num(r.tArea_xz) << ',' << num(r.tArea_yz) << ','
                << num(r.qrsLambda1) << ',' << num(r.qrsLambda2) << ',' << num(r.qrsLambda3) << ','
                << num(r.planarity) << ',' << num(r.planarity_l2_l1) << ','
                << num(r.peakSpatialVelocity) << ',' << num(r.loopVariability) << ','
                << num(r.qrsDuration_ms) << ',' << r.note << "\n";
        }
        return static_cast<bool>(f);
    }

}  // namespace vcg_avg