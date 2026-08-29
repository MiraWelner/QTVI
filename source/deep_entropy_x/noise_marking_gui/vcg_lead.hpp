#pragma once
//
// vcg_lead.hpp
//
// Turns the recorded ECG channels into ONE derived lead that the rest of the
// GUI treats as an ordinary channel: it is markable, it scrolls, it gets a
// peak finder and annotations like any other. It is rendered in the
// derived-lead chart view (ui->kors_matrix).
//
// This is the GUI-side adapter. The transform itself lives in vcg.hpp; this
// file's only job is the part specific to continuous scrolling data:
// getting three channels onto ONE time grid before combining them, and
// producing output in the exact shape channelRefs() hands out
// (a QVector<double> of upsampled samples plus a QVector<QPointF> of raw
// (t, v) pairs in chunk-local seconds).
//
// ---------------------------------------------------------------------------
// WHY THE RAW (t, v) PAIRS AND NOT THE UPSAMPLED VECTORS
// ---------------------------------------------------------------------------
// A VCG is a per-instant linear combination, so all inputs must be at the
// same instant. The upsampled vectors are indexed by sample, and the three ECG
// channels can sit at different upsampled rates
// (channel_upsampled_rates[CH_ECG1..3]) -- index i is then a different time in
// each, and combining by index would silently mix instants that are
// milliseconds apart. The raw pairs carry an explicit timestamp per sample, so
// this builds the derived lead on ECG1's timestamps and INTERPOLATES the other
// channels at those times. If all three rates happen to be equal the
// interpolation is an identity and costs a comparison per sample.
//
// A gap (NaN) in any contributing channel makes that instant NaN in the
// output rather than 0: a partial instant is not a vector, and a zero there
// would read downstream as a real low-voltage segment.
//

#include "vcg.hpp"

#include <QPointF>
#include <QString>
#include <QVector>

#include <cmath>
#include <cstdio>
#include <string>
#include <limits>

namespace vcg_lead {

    /// Label the derived channel is known by everywhere in the GUI.
    inline const char* kLabel = "VCG";

    struct Config {
        /// Built ONCE per file and reused for every chunk. A basis rebuilt per
        /// chunk gives a markable channel whose meaning drifts as the operator
        /// scrolls: an annotation placed in chunk 3 would refer to a different
        /// linear combination of the leads than the same trace in chunk 4.
        /// build() populates this on first use; the FILE-open path must clear
        /// it (see rebuild's note), and the chunk-load path must not.
        vcg::OrthoBasis       ortho;
        vcg::OrthoAccumulator orthoAcc;

        /// Per-lead polarity correction (ECG1/ECG2/ECG3), decided ONCE from
        /// the first chunk via vcg::checkLimbLeadPolarity and applied to
        /// every chunk from then on -- same lifecycle as `ortho`: it exists
        /// specifically so a sign decided on chunk 1 stays in force for
        /// chunk 40, rather than the derived lead's polarity drifting as the
        /// operator scrolls. -1 = that recorded channel is inverted and gets
        /// negated before use. Left at {1,1,1} (a no-op) when the first
        /// chunk's leads were not consistent with limb leads at all -- see
        /// build()'s note; this does not force a "best of 8" answer onto
        /// channels the check cannot actually adjudicate.
        int leadSign[3] = { 1, 1, 1 };

        /// Current "operator has flagged this lead inverted" state for
        /// ECG1/ECG2/ECG3 -- the ecg_N_reverse checkbox, whether checked by
        /// hand or by autoDetectLeadPolarity. Unlike leadSign/ortho, this is
        /// NOT cached: the caller must set it fresh from the live checkbox
        /// state (ui->ecg_N_reverse->isChecked()) immediately before every
        /// call to build()/rebuild(), because the operator can check or
        /// uncheck the box at any time and the very next rebuild has to see
        /// that change. true = withhold the VCG entirely for this call (see
        /// build()'s note) rather than compute one from a lead that is
        /// currently flagged untrustworthy.
        bool leadFlaggedInverted[3] = { false, false, false };

        /// Where the basis is written when one is built, and the name it is
        /// filed under. Empty path = do not write. Set at FILE open, next to
        /// the `ortho` reset: a basis is built exactly once per file, so the
        /// write fires exactly once per file with no bookkeeping flag.
        std::string basisCsvPath;
        std::string basisCsvSubject;

        /// PC1, the max-variance axis. NOT VectorMagnitude: that is invariant
        /// under any orthonormal basis, so selecting it would make the whole
        /// SVD step a no-op (vcg.hpp, SVD section).
        vcg::DerivedLead derived = vcg::DerivedLead::X;
    };

    struct Built {
        QVector<double>  upsampled;   ///< same length/rate as ECG1's upsampled vector
        QVector<QPointF> raw;         ///< (chunk-local seconds, value)
        bool     valid = false;
        QString  why;                 ///< why not, for the status bar / panel title
        QString  basisName;           ///< shown in the panel title; names the method
    };

    // -----------------------------------------------------------------------
    // Interpolation onto the reference channel's timestamps
    // -----------------------------------------------------------------------
    /**
     * @brief Value of `src` at time t, linearly interpolated.
     *
     * @param cursor  Carried between calls. Both series are time-ordered, so
     *                the search resumes where the last one ended, making the
     *                whole resample O(n) instead of O(n log n).
     *
     * @return NaN before the first or after the last sample of `src` (the
     *         derived lead does not extrapolate past a channel's coverage) and
     *         NaN if either bracketing sample is NaN.
     */
    inline double sampleAtTime(const QVector<QPointF>& src, double t, int& cursor) {
        const int n = src.size();
        const double kNaN = std::numeric_limits<double>::quiet_NaN();
        if (n == 0) return kNaN;
        if (cursor < 0) cursor = 0;
        if (cursor > n - 1) cursor = n - 1;

        while (cursor > 0 && src[cursor].x() > t) --cursor;
        while (cursor + 1 < n && src[cursor + 1].x() <= t) ++cursor;

        if (t < src[0].x() || t > src[n - 1].x()) return kNaN;
        if (cursor + 1 >= n) {
            return (t == src[n - 1].x()) ? src[n - 1].y() : kNaN;
        }

        const double t0 = src[cursor].x(), t1 = src[cursor + 1].x();
        const double v0 = src[cursor].y(), v1 = src[cursor + 1].y();
        if (std::isnan(v0) || std::isnan(v1)) return kNaN;
        if (!(t1 > t0)) return v0;
        const double f = (t - t0) / (t1 - t0);
        return v0 * (1.0 - f) + v1 * f;
    }

    /// True when a channel has no usable content (empty, or all NaN).
    inline bool isEmptyChannel(const QVector<QPointF>& raw) {
        for (const QPointF& p : raw) if (!std::isnan(p.y())) return false;
        return true;
    }

    // -----------------------------------------------------------------------
    // Build
    // -----------------------------------------------------------------------
    /**
     * @brief Build the derived lead for the currently loaded chunk.
     *
     * @param ecgRaw       Raw (t, v) pairs for ECG1/2/3, in chunk-local
     *                     seconds. ecgRaw[0] is the REFERENCE grid: the output
     *                     has one sample per ECG1 sample, at ECG1's times.
     * @param upsampledLen Length of ECG1's upsampled vector, so the derived
     *                     upsampled vector matches what every consumer of a
     *                     channel expects. 0 => derive from raw only.
     *
     * @return `valid == false` with `why` set when the inputs cannot support a
     *         reconstruction. Callers should then leave the derived-lead panel
     *         empty and inactive rather than plotting a partial result.
     */
    inline Built build(const QVector<QPointF>* ecgRaw[3],
        int upsampledLen,
        Config& cfg) {
        Built out;

        // Refuse outright if any contributing lead is currently flagged
        // inverted -- whether the operator checked the box by hand or
        // autoDetectLeadPolarity set it. Checked FIRST, unconditionally,
        // before even leadSign's own internal correction below gets a
        // chance to run: "flagged inverted" means the recorded lead is not
        // trusted AS RECORDED, and this GUI's live signal (unlike the
        // offline pipeline, which actually applies the flag to correct the
        // signal) never gets corrected in place -- the checkbox here is
        // only consulted, later, downstream. A VCG computed from a lead in
        // that state -- even one this file's own polarity measurement could
        // internally compensate for -- would be silently untrustworthy, so
        // it is withheld and flagged instead. Clearing the flag (the
        // operator unchecks the box) is what brings the VCG back.
        for (int c = 0; c < 3; ++c) {
            if (cfg.leadFlaggedInverted[c]) {
                out.why = QString("ECG%1 is flagged inverted; VCG withheld until corrected")
                    .arg(c + 1);
                return out;
            }
        }

        for (int c = 0; c < 3; ++c) {
            if (!ecgRaw[c] || isEmptyChannel(*ecgRaw[c])) {
                // The single-channel case lands here and is normal, not a
                // fault: one projection cannot give three axes. The caller
                // shows nothing and collapses the panel.
                out.why = QString("VCG needs three ECG channels; ECG%1 is absent")
                    .arg(c + 1);
                return out;
            }
        }
        const QVector<QPointF>& ref = *ecgRaw[0];
        const int n = ref.size();
        if (n < 2) { out.why = "ECG1 has fewer than 2 samples"; return out; }

        // Resample every channel onto ECG1's timestamps. ch0 is already on
        // that grid by definition and is copied, not interpolated, so the
        // reference channel is never altered by its own resample.
        std::vector<double> lane[3];
        for (int c = 0; c < 3; ++c) lane[c].resize(n);
        int cursor[3] = { 0, 0, 0 };
        for (int i = 0; i < n; ++i) {
            const double t = ref[i].x();
            lane[0][i] = ref[i].y();
            for (int c = 1; c < 3; ++c)
                lane[c][i] = sampleAtTime(*ecgRaw[c], t, cursor[c]);
        }

        // ------------------------------------------------------------------
        // Basis
        // ------------------------------------------------------------------
        // The lanes are on one grid by now, which is the precondition the basis
        // computation shares with the transform: a covariance accumulated
        // across mismatched instants describes a mixture, not the lead
        // geometry.
        const std::vector<const std::vector<double>*> chans{ &lane[0], &lane[1], &lane[2] };

        vcg::VcgMatrix mat;
        if (!cfg.ortho.valid) {          // first chunk of this file only
            // Measured, not asked: decide which (if any) recorded lead is
            // polarity-inverted relative to the other two, per Einthoven's
            // law -- BEFORE this chunk contributes to the whole-file basis,
            // so the basis itself is built from already-corrected leads. A
            // throwaway accumulator: cfg.orthoAcc must only ever receive
            // sign-corrected data, so it cannot be probed with directly.
            vcg::OrthoAccumulator probeAcc;
            probeAcc.addLanes(chans);
            const vcg::PolarityCheckResult pc = vcg::checkLimbLeadPolarity(probeAcc);
            if (!pc.why.empty()) {
                fprintf(stderr, "[vcg] polarity check skipped: %s -- leaving "
                    "lead signs as recorded\n", pc.why.c_str());
            }
            else if (!pc.consistentWithLimbLeads) {
                fprintf(stderr, "[vcg] polarity check: best sign combo "
                    "ECG1=%+d ECG2=%+d ECG3=%+d still has normalized residual "
                    "%.4f (threshold 0.10) -- not explainable by a sign flip "
                    "alone (wrong lead set, or a per-channel gain/calibration "
                    "mismatch); leaving lead signs as recorded\n",
                    pc.sign[0], pc.sign[1], pc.sign[2], pc.normalizedResidual);
            }
            else {
                for (int c = 0; c < 3; ++c) cfg.leadSign[c] = pc.sign[c];
                if (pc.sign[0] < 0 || pc.sign[1] < 0 || pc.sign[2] < 0) {
                    fprintf(stderr, "[vcg] auto-corrected polarity: ECG%d is "
                        "inverted relative to the other two (normalized "
                        "residual %.4f) -- flipped for the rest of this file\n",
                        (pc.sign[0] < 0) ? 1 : (pc.sign[1] < 0) ? 2 : 3,
                        pc.normalizedResidual);
                }
            }

            // Apply the (possibly identity) correction before it contributes
            // to the basis. Cached signs are applied to every chunk further
            // below, so this line is this chunk's copy of that same step.
            for (int c = 0; c < 3; ++c)
                if (cfg.leadSign[c] < 0)
                    for (int i = 0; i < n; ++i) lane[c][i] = -lane[c][i];

            cfg.orthoAcc.addLanes(chans);
            vcg::OrthoBasis nb = cfg.orthoAcc.finish();
            vcg::fixSigns(nb, chans);
            if (!nb.valid) {
                // Not enough usable data to define a basis. Refusing beats
                // silently falling back to kIdentity, which would relabel
                // three correlated channels as orthogonal axes without
                // saying so anywhere the operator can see.
                out.why = QString::fromStdString(nb.why);
                return out;
            }
            cfg.ortho = nb;
            // The one moment a basis comes into existence. Writing here
            // rather than from the caller is what makes "once per file"
            // structural instead of something a flag has to remember.
            if (!cfg.basisCsvPath.empty()) {
                std::string bwhy;
                if (!vcg::writeBasisCsv(cfg.ortho, cfg.basisCsvPath,
                    cfg.basisCsvSubject, bwhy))
                    fprintf(stderr, "[vcg] basis CSV failed: %s\n", bwhy.c_str());
                else
                    fprintf(stderr, "[vcg] wrote %s\n", cfg.basisCsvPath.c_str());
            }
        }
        else {
            // Every later chunk: apply the SAME cached correction decided on
            // chunk 1, so this chunk is on the same polarity convention the
            // basis (and its mean, subtracted below) was built from.
            for (int c = 0; c < 3; ++c)
                if (cfg.leadSign[c] < 0)
                    for (int i = 0; i < n; ++i) lane[c][i] = -lane[c][i];
        }
        mat = cfg.ortho.mat;

        // The basis was derived from mean-removed data, so project the same
        // quantity: leaving the mean in adds a DC offset along whichever
        // axes it happens to project onto, which the basis never saw and
        // the sign convention was not measured against.
        for (int c = 0; c < 3; ++c)
            for (int i = 0; i < n; ++i) lane[c][i] -= cfg.ortho.mean[c];

        // ------------------------------------------------------------------
        // Transform
        // ------------------------------------------------------------------
        const vcg::VcgResult xyz = vcg::reconstructVCG(chans, mat);
        if (!xyz.valid) {
            out.why = QString::fromStdString(xyz.why);
            return out;
        }
        out.basisName = QString("%1  %2  (planarity %3)")
            .arg(QString::fromStdString(cfg.ortho.label))
            .arg(QString::fromUtf8(vcg::orthoAxisName(cfg.derived)))
            .arg(cfg.ortho.planarity(), 0, 'f', 3);

        // ------------------------------------------------------------------
        // Flatten to the single markable trace
        // ------------------------------------------------------------------
        const std::vector<double> scalar = vcg::derivedLeadTrace(xyz, cfg.derived);
        if (static_cast<int>(scalar.size()) != n) {
            out.why = "internal: derived trace length mismatch";
            return out;
        }

        out.raw.resize(n);
        for (int i = 0; i < n; ++i)
            out.raw[i] = QPointF(ref[i].x(), scalar[i]);

        // The upsampled vector is what the renderer draws and what the marker
        // hit-testing indexes, so it has to be the same length as the other
        // channels' upsampled vectors. Resample the derived trace onto that
        // uniform grid spanning the same time range.
        if (upsampledLen > 1) {
            const double t0 = ref.front().x(), t1 = ref.back().x();
            const double dt = (t1 - t0) / static_cast<double>(upsampledLen - 1);
            out.upsampled.resize(upsampledLen);
            int cur = 0;
            for (int i = 0; i < upsampledLen; ++i)
                out.upsampled[i] = sampleAtTime(out.raw, t0 + i * dt, cur);
        }
        else {
            out.upsampled.resize(n);
            for (int i = 0; i < n; ++i) out.upsampled[i] = scalar[i];
        }

        out.valid = true;
        return out;
    }

    /**
     * @brief Convenience wrapper matching how the GUI holds its channels.
     *        Returns false and leaves the destinations untouched on failure,
     *        so a stale-but-valid lead is never half-overwritten by a failed
     *        rebuild.
     *
     * @param cfg  NON-CONST: the SVD basis is cached here on the first chunk
     *             and reused for the rest of the file. The caller must clear
     *             it when a NEW FILE is opened --
     *                 cfg.ortho = vcg::OrthoBasis{};
     *                 cfg.orthoAcc.reset();
     *                 cfg.leadSign[0] = cfg.leadSign[1] = cfg.leadSign[2] = 1;
     *             -- and must NOT clear it on a chunk load. Clearing per chunk
     *             is the failure mode this cache exists to prevent: the derived
     *             lead would silently change meaning under the operator's
     *             existing annotations. Carrying it across files is the other
     *             direction of the same bug.
     *
     *             Separately, and on EVERY call (not just file-open): the
     *             caller must set cfg.leadFlaggedInverted[0..2] fresh from
     *             the live ecg_N_reverse checkbox state right before calling
     *             this --
     *                 cfg.leadFlaggedInverted[0] = ui->ecg_1_reverse->isChecked();
     *                 cfg.leadFlaggedInverted[1] = ui->ecg_2_reverse->isChecked();
     *                 cfg.leadFlaggedInverted[2] = ui->ecg_3_reverse->isChecked();
     *             -- so a box the operator just checked (or unchecked)
     *             withholds (or restores) the VCG on the very next rebuild,
     *             not on the next file open.
     */
    inline bool rebuild(const QVector<QPointF>& e1,
        const QVector<QPointF>& e2,
        const QVector<QPointF>& e3,
        int upsampledLen,
        Config& cfg,
        QVector<double>& dstUpsampled,
        QVector<QPointF>& dstRaw,
        QString& status) {
        const QVector<QPointF>* in[3] = { &e1, &e2, &e3 };
        Built b = build(in, upsampledLen, cfg);
        status = b.valid ? b.basisName : b.why;
        if (!b.valid) return false;
        dstUpsampled = std::move(b.upsampled);
        dstRaw = std::move(b.raw);
        return true;
    }

}  // namespace vcg_lead