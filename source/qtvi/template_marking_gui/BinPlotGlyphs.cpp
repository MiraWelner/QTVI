// ============================================================================
// BinPlotGlyphs.cpp
//
// Feature-glyph QC layer for BinPlotWidget. Landmark positions are captured
// ONCE at setData() (captureGlyphSnapshot) from the initial auto-placed
// markers and frozen -- they do NOT follow marker drags, so they stay a fixed
// reference for where the algorithm first placed things. drawFeatureGlyphs()
// just renders the stored snapshot.
//
// Marks (all small opaque black):
//   ECG - P peak, R peak, T peak (computed) + Q onset, S end, T end (markers)
//   PPG - foot (global argmin), systolic peak (#1), dicrotic notch OR an 'o'
//         at the [peak,end] midpoint when no notch is present, diastolic peak.
// ECG glyphs use the left-axis (yLo,yHi) scale + ECG x-geometry; PPG glyphs
// use the shared right-axis (pLo,pHi) scale + foot-anchored x-geometry.
// ============================================================================
#include "BinPlotWidget.hpp"
#include "template_marking_gui\feature_marks.hpp"
#include "feature_marks.hpp"
#include <QPainter>
#include <algorithm>
#include <cmath>
#include <vector>

void BinPlotWidget::captureGlyphSnapshot() {
    m_glyphs = GlyphSnapshot{};   // reset to all -1 / false

    auto median_of = [](std::vector<double> w) -> double {
        if (w.empty()) return 0.0;
        std::sort(w.begin(), w.end());
        return w[w.size() / 2];
        };
    auto inTrace = [](const std::vector<double>& v, int i) {
        return i >= 0 && i < (int)v.size() && !std::isnan(v[i]);
        };

    // ---- ECG landmarks ---------------------------------------------------
    // Six reactive X glyphs, auto-computed but tracking the movable markers
    // live: P peak, Q begin, R peak, S end, T peak, T end. Q-peak / S-peak
    // are no longer part of the set.
    if ((int)m_ecg.size() >= 3) {
        auto e = FeatureMarks::compute_ecg_glyphs(m_ecg,
            m_markers[EcgPPeak], m_markers[EcgQBegin],
            m_markers[EcgSEnd], m_markers[EcgTPeak], m_markers[EcgTEnd],
            m_sampleRate);
        m_glyphs.ecgPPeak = e.p_peak_glyph;
        m_glyphs.ecgQ = e.q_begin_glyph;
        m_glyphs.ecgRPeak = e.r_peak_glyph;
        m_glyphs.ecgS = e.s_end_glyph;
        m_glyphs.ecgTPeak = e.t_peak_glyph;
        m_glyphs.ecgTend = e.t_end_glyph;
        m_glyphs.ecgQPeak = -1;   // dropped from the reactive set
        m_glyphs.ecgSPeak = -1;
    }

    // ---- PPG landmarks ---------------------------------------------------
    //   foot     = the PpgOnset marker (already seeded to argmin after R)
    //   P1       = argmax over (foot, foot + 0.75*(N-foot))  ["first 75%"]
    //              X if a real local max exists in that window; else O at
    //              the window midpoint.
    //   P50      = 50% point on the upslope foot -> P1 (value-based, i.e.
    //              the sample nearest 0.5*(P1-foot) + foot amplitude).
    //              X if the trace actually crosses that value on the
    //              rising limb; else O at the temporal midpoint.
    //   end      = argmin over (P1, N)
    //   dicrotic = FIRST interior local min over (P1, end); X if found,
    //              O at midpoint of (P1, end) if not
    //   P2       = local max over (dicrotic, end); X if found, else O
    //              at the temporal midpoint of that window.
    if (m_hasPPG && (int)m_ppg.size() >= 3) {
        const std::vector<double>& v = m_ppg;
        const int N = (int)v.size();
        const int foot = m_markers[PpgOnset];

        // P1: max in first 75% after foot.
        int p1 = -1, p1O = -1;
        if (foot >= 0) {
            const int p1hi = foot + (3 * (N - foot)) / 4;   // first 75%
            double best = -1e300;
            for (int i = foot + 1; i < p1hi && i < N - 1; ++i) {
                if (std::isnan(v[i]) || std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                if (v[i] >= v[i - 1] && v[i] >= v[i + 1] && v[i] > best) {
                    best = v[i]; p1 = i;
                }
            }
            if (p1 < 0) p1O = (foot + p1hi) / 2;
        }

        // Pick the "real" P1 sample position for downstream searches:
        // fall back to the O position if no real max was found, so end/
        // notch/P2 still have a reference.
        const int p1Ref = (p1 >= 0) ? p1 : p1O;

        // 50%-upslope: sample nearest half amplitude between foot and P1.
        int p50 = -1, p50O = -1;
        if (foot >= 0 && p1Ref > foot &&
            !std::isnan(v[foot]) && !std::isnan(v[p1Ref])) {
            const double target = 0.5 * (v[foot] + v[p1Ref]);
            double bestDiff = 1e300;
            for (int i = foot; i <= p1Ref; ++i) {
                if (std::isnan(v[i])) continue;
                const double d = std::abs(v[i] - target);
                if (d < bestDiff) { bestDiff = d; p50 = i; }
            }
            // Sanity: if we couldn't get within a reasonable band of the
            // target (extreme jitter or all-NaN), mark absent.
            const double band = 0.25 * std::abs(v[p1Ref] - v[foot]);
            if (p50 < 0 || bestDiff > band) {
                p50 = -1;
                p50O = (foot + p1Ref) / 2;
            }
        }
        else if (foot >= 0 && p1Ref > foot) {
            p50O = (foot + p1Ref) / 2;
        }

        int end = -1;
        if (p1Ref >= 0) {
            double best = 1e300;
            for (int i = p1Ref + 1; i < N; ++i)
                if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; end = i; }
        }

        // Dicrotic notch detection unchanged: least-negative slope in the
        // middle 50% of the downslope.
        int notchIdx = -1;
        int notchOIdx = -1;
        if (p1Ref >= 0 && end > p1Ref + 4) {
            int firstMin = -1;
            for (int i = p1Ref + 1; i < end && i < N - 1; ++i) {
                if (std::isnan(v[i]) || std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                if (v[i] <= v[i - 1] && v[i] <= v[i + 1]) { firstMin = i; break; }
            }
            const int slopeStart = (firstMin > 0) ? firstMin : p1Ref;

            const int span = end - slopeStart;
            if (span >= 4) {
                const int mid_lo = slopeStart + span / 4;
                const int mid_hi = slopeStart + (3 * span) / 4;

                int bestIdx = -1;
                double bestSlope = -1e300;
                for (int i = std::max(1, mid_lo); i <= std::min(N - 2, mid_hi); ++i) {
                    if (std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                    const double slope = 0.5 * (v[i + 1] - v[i - 1]);
                    if (slope > bestSlope) { bestSlope = slope; bestIdx = i; }
                }

                if (bestIdx >= 0 && bestSlope >= 0.0) notchIdx = bestIdx;
                else                                  notchOIdx = (slopeStart + end) / 2;
            }
        }

        // P2: local max between notch and end. O at temporal midpoint of
        // that window if no local max exists (or notch is absent).
        int p2 = -1, p2O = -1;
        const int p2Lo = (notchIdx >= 0) ? notchIdx : (notchOIdx >= 0 ? notchOIdx : p1Ref);
        if (p2Lo >= 0 && end > p2Lo + 1) {
            double best = -1e300;
            for (int i = p2Lo + 1; i < end; ++i) {
                if (std::isnan(v[i]) || std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                if (v[i] >= v[i - 1] && v[i] >= v[i + 1] && v[i] > best) {
                    best = v[i]; p2 = i;
                }
            }
            if (p2 < 0) p2O = (p2Lo + end) / 2;
        }

        m_glyphs.ppgFoot = foot;
        m_glyphs.ppgP50 = p50;   m_glyphs.ppgP50OFallback = p50O;
        m_glyphs.ppgP1 = p1;    m_glyphs.ppgP1OFallback = p1O;
        m_glyphs.ppgDic = notchIdx;
        m_glyphs.ppgP2 = p2;    m_glyphs.ppgP2OFallback = p2O;
        m_glyphs.ppgEnd = end;
        m_glyphs.ppgNotch = (notchIdx >= 0);
        m_glyphs.ppgNoNotchO = notchOIdx;
    }

    m_glyphs.valid = true;
}

void BinPlotWidget::drawFeatureGlyphs(QPainter& p,
    double yLo, double yHi, double pLo, double pHi, int ph) const
{
    if (!m_glyphs.valid) return;

    auto plotY = [&](double val, double lo, double hi) {
        const double r = (hi - lo > 1e-10) ? (hi - lo) : 1.0;
        return margin_top + ph - (val - lo) / r * ph;
        };
    auto glyph = [&](double x, double y) {   // opaque black "X"
        p.setBrush(Qt::NoBrush);
        p.setPen(QPen(Qt::black, 1.8));
        const double s = 4.0;
        p.drawLine(QPointF(x - s, y - s), QPointF(x + s, y + s));
        p.drawLine(QPointF(x - s, y + s), QPointF(x + s, y - s));
        };

    if (m_showEcgTrace && (int)m_ecg.size() >= 3) {
        const std::vector<double>& v = m_ecg;
        const int N = (int)v.size();
        double baseline = 0.0;
        {
            const int bhi = std::min(10, N);
            std::vector<double> bw;
            for (int i = 0; i < bhi; ++i) if (!std::isnan(v[i])) bw.push_back(v[i]);
            if (!bw.empty()) { std::sort(bw.begin(), bw.end()); baseline = bw[bw.size() / 2]; }
        }
        auto g = [&](int idx) {
            if (idx < 0 || idx >= N) return;
            const double val = std::isnan(v[idx]) ? baseline : v[idx];
            glyph(xFromSample(idx, /*isEcg=*/true), plotY(val, yLo, yHi));
            };
        // Six ECG X glyphs, from the reactive computes.
        g(m_glyphs.ecgPPeak);   // P wave
        g(m_glyphs.ecgQ);       // Q onset
        g(m_glyphs.ecgRPeak);   // R wave
        g(m_glyphs.ecgS);       // S end
        g(m_glyphs.ecgTPeak);   // T peak
        g(m_glyphs.ecgTend);    // T end
    }

    if (m_showPpgTrace && m_hasPPG && (int)m_ppg.size() >= 3) {
        const std::vector<double>& v = m_ppg;
        const int N = (int)v.size();
        // Clip PPG glyphs against the PPG's own visible extent, not the
        // ECG length. Marker lines already draw against the PPG length --
        // matching that here avoids "marker line visible but no glyph"
        // when ECG and PPG templates aren't exactly the same length.
        auto g = [&](int idx) {
            if (idx < 0 || idx >= N) return;
            const double val = std::isnan(v[idx]) ? pLo : v[idx];
            glyph(xFromSample(idx, /*isEcg=*/false), plotY(val, pLo, pHi));
            };
        auto circ = [&](int idx) {
            if (idx < 0 || idx >= N) return;
            const double val = std::isnan(v[idx]) ? pLo : v[idx];
            const double x = xFromSample(idx, /*isEcg=*/false);
            const double y = plotY(val, pLo, pHi);
            p.setBrush(Qt::NoBrush);
            p.setPen(QPen(Qt::black, 1.8));
            p.drawEllipse(QPointF(x, y), 4.0, 4.0);
            };
        // Five glyphs: foot, P50 (X/O), P1 (X/O), dicrotic (X/O), P2 (X/O),
        // end (X). Foot and end always draw as X at the seeded positions.
        // Everything else falls back to O at a sensible midpoint.
        g(m_glyphs.ppgFoot);
        if (m_glyphs.ppgP50 >= 0) g(m_glyphs.ppgP50);
        else                      circ(m_glyphs.ppgP50OFallback);
        if (m_glyphs.ppgP1 >= 0)  g(m_glyphs.ppgP1);
        else                      circ(m_glyphs.ppgP1OFallback);
        if (m_glyphs.ppgNotch)    g(m_glyphs.ppgDic);
        else                      circ(m_glyphs.ppgNoNotchO);
        if (m_glyphs.ppgP2 >= 0)  g(m_glyphs.ppgP2);
        else                      circ(m_glyphs.ppgP2OFallback);
        g(m_glyphs.ppgEnd);
    }
}