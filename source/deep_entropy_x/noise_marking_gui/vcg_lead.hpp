#pragma once
//
// vcg_lead.hpp
//
// Turns the recorded ECG channels into ONE derived lead that the rest of the
// GUI treats as an ordinary channel: it is markable, it scrolls, it gets a
// peak finder and annotations like any other. It is rendered in the
// kors_matrix chart view.
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
#include <limits>

namespace vcg_lead {

    /// Label the derived channel is known by everywhere in the GUI.
    inline const char* kLabel = "VCG";

    struct Config {
        /// The 3x3 transform: rows X/Y/Z, columns ECG1/ECG2/ECG3.
        vcg::VcgMatrix   matrix = vcg::kIdentity;
        vcg::DerivedLead derived = vcg::DerivedLead::VectorMagnitude;
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
     *         reconstruction. Callers should then leave the kors_matrix panel
     *         empty and inactive rather than plotting a partial result.
     */
    inline Built build(const QVector<QPointF>* ecgRaw[3],
        int upsampledLen,
        const Config& cfg) {
        Built out;

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
        // Transform
        // ------------------------------------------------------------------
        const std::vector<const std::vector<double>*> chans{ &lane[0], &lane[1], &lane[2] };
        const vcg::VcgResult xyz = vcg::reconstructVCG(chans, cfg.matrix);
        if (!xyz.valid) {
            out.why = QString::fromStdString(xyz.why);
            return out;
        }
        out.basisName = QString("%1 / %2")
            .arg(QString::fromUtf8(xyz.basisName))
            .arg(QString::fromUtf8(vcg::derivedLeadName(cfg.derived)));

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
     */
    inline bool rebuild(const QVector<QPointF>& e1,
        const QVector<QPointF>& e2,
        const QVector<QPointF>& e3,
        int upsampledLen,
        const Config& cfg,
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