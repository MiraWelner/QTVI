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
// ---------------------------------------------------------------------------
// THE BASIS IS RECORD-LOCAL, AND THE CSV SAYS SO
// ---------------------------------------------------------------------------
// The save path (analyzeBinWithBasis) projects every bin onto ONE SVD basis
// built from this record's own templates (buildRecordBasis). Its axes are
// PC1/PC2/PC3 of THIS record: not Frank X/Y/Z, and not the same directions in
// the next record. So:
//
//   comparable across bins of one record   every plane area, the QRS-T angle
//   comparable across records              only the rotation-invariant
//                                          features: the QRS-T angle, the
//                                          lambdas and planarity, peak spatial
//                                          velocity and loop variability
//   NOT comparable across records          the per-plane areas (xy/xz/yz are
//                                          PC1-PC2, PC1-PC3, PC2-PC3 here)
//
// Every row carries basis = record_local_svd and the basis's rank and lead
// signs, and the basis itself is written to <id>_vcg_template_basis.csv.
//
// TWO CHECKS DECIDE WHAT THE BASIS CAN SUPPORT, both run before it is built:
//
//   INVERSION. vcg::checkLimbLeadPolarity tests Einthoven's law (II = I + III)
//   under all eight sign combinations. When the leads are consistent with
//   limb leads, the winning signs are applied to the leads before the basis
//   sees them, exactly as the noise-marking VCG does -- one inverted lead
//   otherwise rotates the whole loop.
//
//   RANK. Limb leads obey II = I + III exactly, so three limb leads span only
//   TWO dimensions: the frontal plane. A third "axis" is then noise, and every
//   quantity measured along it -- the two out-of-plane areas, lambda3, the
//   out-of-plane planarity -- is a measurement of noise that looks like a
//   measurement of the heart. So the rank is 2 when the leads pass the limb-
//   lead test, or when sigma3/sigma1 is below kRankFloor; the loop is then
//   held in the PC1-PC2 plane (PC3 set to 0) and those quantities are left
//   empty. Rank 1 (a line) supports no loop at all.
//

#include "noise_marking_gui/vcg.hpp"
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

        // Shared with the basis derivation in vcg.hpp rather than swept here:
        // this file used to carry its own copy of the same Jacobi iteration,
        // whose convergence test was absolute (off < 1e-15). Loop covariances
        // are in mV^2 and routinely land near that threshold, so the sweep
        // could stop before converging on exactly the small-amplitude loops
        // whose planarity is most in question. symmetricEigen3 tests relative
        // to the diagonal, which is scale-free. Eigenvectors are returned but
        // unused: only the spread magnitudes matter here.
        double ev[3], evec[3][3];
        vcg::symmetricEigen3(C, ev, evec);
        (void)evec;
        // A covariance matrix is positive semi-definite, so a negative
        // eigenvalue here is Jacobi round-off on a degenerate (flat or
        // collinear) loop -- around -1e-18 in practice. Clamped at zero:
        // otherwise the planarity ratio comes out slightly negative, which
        // reads as a nonsensical measurement rather than "perfectly planar".
        for (int k = 0; k < 3; ++k) ev[k] = std::max(0.0, ev[k]);
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
    // The record-local basis (see the header)
    // -----------------------------------------------------------------------
    /// sigma3/sigma1 below this is treated as rank 2 when the limb-lead test
    /// did not already say so. A starting value, not a validated threshold.
    inline constexpr double kRankFloor = 0.05;
    /// Same for sigma2/sigma1 and rank 1.
    inline constexpr double kRank1Floor = 0.05;

    struct RecordBasis {
        vcg::OrthoBasis basis;
        int    leadSign[kNumEcgCh] = { 1, 1, 1 };  ///< applied to ECG1/2/3 first
        bool   polarityChecked = false;      ///< the limb-lead test could run
        bool   limbLeads = false;            ///< ...and the leads passed it
        double limbResidual = kNaN;          ///< its normalized residual
        int    rank = 0;                     ///< 3, 2, or 0 when unusable
        double s3_over_s1 = kNaN;
        long long nInstants = 0;
        bool   valid = false;
        std::string why;                     ///< when !valid
        static constexpr const char* kTag = "record_local_svd";

        std::string signsText() const {
            std::string t;
            for (int c = 0; c < kNumEcgCh; ++c) t += (leadSign[c] < 0) ? '-' : '+';
            return t;
        }
    };

    /// The three lanes of one bin on the R-relative axis [-pre, post], each
    /// channel read at its own R column. False when a channel or R is missing.
    inline bool binLanes(const time_bin& b, int pre, int post,
        std::vector<double>(&lane)[kNumEcgCh]) {
        if (pre < 0 || post < 0 || pre + post < 1) return false;
        const std::vector<double>* tpl[kNumEcgCh] = {
            &b.ch1.ecgTemplate_raw, &b.ch2.ecgTemplate_raw, &b.ch3.ecgTemplate_raw };
        double rc[kNumEcgCh];
        for (int c = 0; c < kNumEcgCh; ++c) {
            rc[c] = (b.r_peak_auto_ch[c] >= 0.0) ? b.r_peak_auto_ch[c]
                : static_cast<double>(b.r_peak_ch[c]);
            if (rc[c] < 0.0 || tpl[c]->size() < 3) return false;
        }
        const int n = pre + post + 1;
        for (int c = 0; c < kNumEcgCh; ++c) {
            lane[c].assign(n, kNaN);
            for (int i = 0; i < n; ++i)
                lane[c][i] = FeatureMarks::sample_at(*tpl[c], rc[c] + (i - pre));
        }
        return true;
    }

    /**
     * @brief ONE basis for the whole record, from its own templates.
     *
     *        Accumulated over every usable bin's three-lead template on
     *        [R - preMs, R + postMs] -- the stretch the loops are cut from --
     *        so the basis describes the same waveforms it is applied to.
     *        Templates rather than raw signal because they are what the save
     *        path has, and because the median has already removed the noise a
     *        raw-signal covariance would fold into PC3.
     *
     *        Order: inversion check, sign correction, accumulate, finish, sign
     *        convention (fixSigns), rank. See the header for why each.
     */
    inline RecordBasis buildRecordBasis(const std::vector<time_bin>& bins,
        double fs, double preMs = 200.0, double postMs = 600.0) {
        RecordBasis rb;
        if (!(fs > 0.0)) { rb.why = "no ECG sample rate"; return rb; }
        const int pre = static_cast<int>(std::lround(preMs * 0.001 * fs));
        const int post = static_cast<int>(std::lround(postMs * 0.001 * fs));

        // Every usable bin's lanes, kept: the sign correction is decided from
        // all of them and then applied to all of them.
        std::vector<std::vector<double>> all[kNumEcgCh];
        for (const time_bin& b : bins) {
            if (b.bad_segment) continue;
            std::vector<double> lane[kNumEcgCh];
            if (!binLanes(b, pre, post, lane)) continue;
            for (int c = 0; c < kNumEcgCh; ++c) all[c].push_back(std::move(lane[c]));
        }
        if (all[0].empty()) { rb.why = "no bin has all three ECG templates"; return rb; }

        auto accumulate = [&](vcg::OrthoAccumulator& acc) {
            for (std::size_t k = 0; k < all[0].size(); ++k)
                acc.addLanes({ &all[0][k], &all[1][k], &all[2][k] });
            };

        // ---- INVERSION --------------------------------------------------
        {
            vcg::OrthoAccumulator probe;
            accumulate(probe);
            const vcg::PolarityCheckResult pc = vcg::checkLimbLeadPolarity(probe);
            rb.polarityChecked = pc.why.empty();
            rb.limbResidual = pc.normalizedResidual;
            rb.limbLeads = rb.polarityChecked && pc.consistentWithLimbLeads;
            // Only a PASSING test is allowed to flip a lead: when no sign
            // combination makes the leads obey Einthoven's law they are not
            // limb leads, and the best of eight is an arbitrary answer.
            if (rb.limbLeads)
                for (int c = 0; c < kNumEcgCh; ++c) rb.leadSign[c] = pc.sign[c];
        }
        for (int c = 0; c < kNumEcgCh; ++c)
            if (rb.leadSign[c] < 0)
                for (auto& lane : all[c]) for (double& v : lane) v = -v;

        // ---- BASIS --------------------------------------------------------
        vcg::OrthoAccumulator acc;
        accumulate(acc);
        rb.nInstants = acc.n;
        rb.basis = acc.finish(/*minSamples=*/500);
        if (!rb.basis.valid) { rb.why = rb.basis.why; return rb; }
        {
            // The sign convention needs the traces; one concatenated lane set.
            std::vector<double> cat[kNumEcgCh];
            for (int c = 0; c < kNumEcgCh; ++c)
                for (const auto& lane : all[c]) cat[c].insert(cat[c].end(), lane.begin(), lane.end());
            vcg::fixSigns(rb.basis, { &cat[0], &cat[1], &cat[2] });
        }
        rb.basis.label = "record-local SVD (PC1/PC2/PC3)";
        rb.basis.mat.name = "record-local SVD (PC1/PC2/PC3)";

        // ---- RANK ---------------------------------------------------------
        const double s1 = rb.basis.sigma[0];
        rb.s3_over_s1 = (s1 > 0.0) ? rb.basis.sigma[2] / s1 : kNaN;
        const double s2_over_s1 = (s1 > 0.0) ? rb.basis.sigma[1] / s1 : kNaN;
        if (!(s2_over_s1 >= kRank1Floor)) {
            rb.why = "the three leads span one direction (rank 1): no loop";
            rb.rank = 1;
            return rb;
        }
        rb.rank = (rb.limbLeads || !(rb.s3_over_s1 >= kRankFloor)) ? 2 : 3;
        rb.valid = true;
        return rb;
    }

    // -----------------------------------------------------------------------
    // One bin's full feature row
    // -----------------------------------------------------------------------
    struct BinFeatures {
        int    binIndex = -1;
        int    nBeats = 0;         ///< beats that produced a loop (post-exclusion)
        // ALWAYS 0 NOW, and the CSV column with it. It was written only by
        // the per-beat path (perBeatLoops excluded ectopics as it built the
        // loops); analyzeBinFromTemplates works from the per-lead templates,
        // which the ectopic mask has already been applied to upstream, so
        // there is nothing for this to count. Kept so the CSV layout does not
        // move.
        int    nEctopicExcluded = 0;  ///< always 0: see above
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
    inline Loop loopFromTemplates(const time_bin& b, int pre, int post,
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
    inline std::vector<double> derivedTraceOnChannelAxis(const time_bin& b,
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

    /// The loop on the record basis: lead signs applied, projected onto
    /// PC1/PC2/PC3 about the signal's own zero (NOT the record mean -- the
    /// QRS-T angle is measured from the origin, and a rotation preserves it),
    /// and held in the PC1-PC2 plane when the basis is rank 2.
    inline Loop loopOnRecordBasis(const time_bin& b, int pre, int post,
        const RecordBasis& rb) {
        Loop out;
        out.firstOffset = -pre;
        if (!rb.valid) return out;
        std::vector<double> lane[kNumEcgCh];
        if (!binLanes(b, pre, post, lane)) return out;
        for (int c = 0; c < kNumEcgCh; ++c)
            if (rb.leadSign[c] < 0) for (double& v : lane[c]) v = -v;
        const vcg::VcgResult r = vcg::reconstructVCG(
            std::vector<const std::vector<double>*>{ &lane[0], & lane[1], & lane[2] },
            rb.basis.mat);
        if (!r.valid) return out;
        out.pts = r.samples;
        if (rb.rank == 2)
            for (vcg::VCGSample& s : out.pts) if (!std::isnan(s.z)) s.z = 0.0;
        return out;
    }

    /// Extract every interval-dependent feature from an already-built loop.
    /// loopVariability is left NaN; pass it separately if it was computed at
    /// build time.
    inline BinFeatures extractFeatures(int binIndex, const Loop& loop,
        const global_intervals::GlobalIntervals& g,
        double fs, int rank = 3) {
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
        // RANK 2: the third axis is not signal (see the header), so nothing
        // measured along it is reported. The in-plane quantities stand.
        if (rank == 2) {
            f.qrsArea_xz = f.qrsArea_yz = kNaN;
            f.tArea_xz = f.tArea_yz = kNaN;
            f.qrsLambda3 = kNaN;
            f.planarity = kNaN;
        }
        f.valid = true;
        return f;
    }

    /// One call for the save path: template loop -> features. `marginSamples`
    /// widens the window past the measured QT so the T loop cannot be clipped.
    inline BinFeatures analyzeBinFromTemplates(int binIndex, const time_bin& b,
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

    /// The save path on the record basis. Same window rule as above; the loop
    /// is projected onto `rb`, and the rank decides which features exist.
    inline BinFeatures analyzeBinWithBasis(int binIndex, const time_bin& b,
        const global_intervals::GlobalIntervals& g, double fs,
        const RecordBasis& rb, int marginSamples = 10) {
        if (!rb.valid) {
            BinFeatures f;
            f.binIndex = binIndex;
            f.note = "no record basis: " + rb.why;
            return f;
        }
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
        return extractFeatures(binIndex, loopOnRecordBasis(b, pre, post, rb), g, fs, rb.rank);
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
            "peak_spatial_velocity,loop_variability,qrs_duration_ms,note,"
            // Appended, so every existing column keeps its position. One
            // basis per record, repeated on every row.
            "basis,basis_rank,lead_signs,limb_lead_residual,basis_s3_over_s1";
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
        const std::vector<BinFeatures>& rows,
        // The basis every row was measured on. Null = the old identity
        // transform, reported as such so the file is never ambiguous.
        const RecordBasis* rb = nullptr) {
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
                << num(r.qrsDuration_ms) << ',' << r.note << ',';
            if (rb && rb->valid)
                f << RecordBasis::kTag << ',' << rb->rank << ',' << rb->signsText()
                << ',' << num(rb->limbResidual) << ',' << num(rb->s3_over_s1) << "\n";
            else if (rb)
                f << RecordBasis::kTag << ",0,," << num(rb->limbResidual) << ",\n";
            else
                f << "identity_not_orthogonal,,,,\n";
        }
        return static_cast<bool>(f);
    }

}  // namespace vcg_avg