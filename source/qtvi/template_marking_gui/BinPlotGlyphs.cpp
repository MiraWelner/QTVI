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

    // ---- ECG landmarks (from the initial markers) ----
    //
    // Peaks (P, R, T) get an X if a genuine local extremum (matching R
    // polarity for P/T; max deviation from baseline for R) exists in the
    // bracketing window, otherwise an O at the window midpoint marks the
    // "expected here but not found" location. Marker landmarks (Q, S end,
    // T end) are drawn as X wherever the marker sits; if the marker slot
    // is unset (-1) they simply don't draw.
    if ((int)m_ecg.size() >= 3) {
        const std::vector<double>& v = m_ecg;
        const int N = (int)v.size();
        const int pOn = m_markers[EcgP];        // P onset
        const int qOn = m_markers[EcgQBegin];   // Q onset
        const int sEn = m_markers[EcgSEnd];     // S end
        const int tBe = m_markers[EcgTBegin];   // T begin
        const int tEn = m_markers[EcgTEnd];     // T end

        // Baseline from the pre-P region (or the first few samples).
        const int bhi = (pOn > 1 && pOn < N) ? pOn : std::min(10, N);
        std::vector<double> bw;
        for (int i = 0; i < bhi; ++i) if (!std::isnan(v[i])) bw.push_back(v[i]);
        const double baseline = median_of(bw);

        // Utility: interior local extremum matching a polarity, over [lo, hi].
        // "up" => local max; "!up" => local min. Skips NaN neighbours. Returns
        // -1 if none found.
        auto interiorExtremum = [&](int lo, int hi, bool up, int& idx)->bool {
            idx = -1;
            lo = std::max(1, lo); hi = std::min(N - 2, hi);
            double best = up ? -1e300 : 1e300;
            for (int i = lo; i <= hi; ++i) {
                if (std::isnan(v[i]) || std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                const bool isMax = (v[i] >= v[i - 1] && v[i] >= v[i + 1]);
                const bool isMin = (v[i] <= v[i - 1] && v[i] <= v[i + 1]);
                if ((up && isMax && v[i] > best) || (!up && isMin && v[i] < best)) {
                    best = v[i]; idx = i;
                }
            }
            return idx >= 0;
            };

        // R peak: greatest deviation from baseline within [Q, S end]. R must
        // exist for the beat to exist, so this window "should" always yield
        // something; if it doesn't (missing markers), fall back to O at the
        // window midpoint.
        int rIdx = -1, rOIdx = -1;
        bool rFound = false;
        if (qOn >= 0 && sEn > qOn) {
            int lo = std::max(0, qOn), hi = std::min(N - 1, sEn);
            double best = -1.0;
            for (int i = lo; i <= hi; ++i) {
                if (std::isnan(v[i])) continue;
                const double d = std::abs(v[i] - baseline);
                if (d > best) { best = d; rIdx = i; }
            }
            rFound = (rIdx >= 0);
            if (!rFound) rOIdx = (qOn + sEn) / 2;
        }
        // Polarity for P/T comes from R's deviation direction (or "up" if R
        // itself is missing -- an assumption, but the O will show that anyway).
        const bool up = (rFound && inTrace(v, rIdx))
            ? (v[rIdx] - baseline) >= 0 : true;

        // P peak: interior local extremum (same polarity as R) in [P on, Q on].
        int pIdx = -1, pOIdx = -1;
        bool pFound = false;
        if (pOn >= 0 && qOn > pOn) {
            pFound = interiorExtremum(pOn, qOn, up, pIdx);
            if (!pFound) pOIdx = (pOn + qOn) / 2;
        }

        // T peak: interior local extremum (same polarity as R) in [T beg, T end].
        int tIdx = -1, tOIdx = -1;
        bool tFound = false;
        if (tBe >= 0 && tEn > tBe) {
            tFound = interiorExtremum(tBe, tEn, up, tIdx);
            if (!tFound) tOIdx = (tBe + tEn) / 2;
        }

        m_glyphs.ecgP = pIdx;   m_glyphs.ecgPOFallback = pOIdx;
        m_glyphs.ecgR = rIdx;   m_glyphs.ecgROFallback = rOIdx;
        m_glyphs.ecgT = tIdx;   m_glyphs.ecgTOFallback = tOIdx;
        m_glyphs.ecgQ = qOn;    m_glyphs.ecgS = sEn;   m_glyphs.ecgTend = tEn;
    }

    // ---- PPG landmarks ---------------------------------------------------
    //   foot     = the PpgOnset marker (already seeded to argmin after R)
    //   P1       = argmax over (foot, N)
    //   end      = argmin over (P1, N)
    //   dicrotic = FIRST interior local min over (P1, end); X if found,
    //              O at midpoint of (P1, end) if not
    //   P2       = local max over (dicrotic, end); no glyph if none
    if (m_hasPPG && (int)m_ppg.size() >= 3) {
        const std::vector<double>& v = m_ppg;
        const int N = (int)v.size();
        const int foot = m_markers[PpgOnset];

        int p1 = -1;
        if (foot >= 0) {
            double best = -1e300;
            for (int i = foot + 1; i < N; ++i)
                if (!std::isnan(v[i]) && v[i] > best) { best = v[i]; p1 = i; }
        }

        int end = -1;
        if (p1 >= 0) {
            double best = 1e300;
            for (int i = p1 + 1; i < N; ++i)
                if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; end = i; }
        }

        // Dicrotic notch detection (rule ii: least-negative slope).
        //
        //   1. firstMin = first local min after P1 (the "shoulder" that
        //      starts the downslope proper). If none, use P1 itself.
        //   2. downslope = [firstMin, end]
        //   3. Search the MIDDLE 50% of the downslope for the point where
        //      the first derivative is closest to zero (least-negative
        //      slope). On a healthy PPG this is where the dicrotic notch's
        //      subtle upward kink lives, even when no full local min forms.
        //   4. If the least-negative slope is genuinely non-negative there
        //      -> X on that point (real notch).
        //   5. If the slope stays negative through the whole middle 50%
        //      -> O at the middle of the downslope (visual "notch would be
        //      here" placeholder).
        int notchIdx = -1;
        int notchOIdx = -1;
        if (p1 >= 0 && end > p1 + 4) {
            // Step 1: first local min after p1.
            int firstMin = -1;
            for (int i = p1 + 1; i < end && i < N - 1; ++i) {
                if (std::isnan(v[i]) || std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                if (v[i] <= v[i - 1] && v[i] <= v[i + 1]) { firstMin = i; break; }
            }
            const int slopeStart = (firstMin > 0) ? firstMin : p1;

            // Step 2/3: middle 50% of [slopeStart, end].
            const int span = end - slopeStart;
            if (span >= 4) {
                const int mid_lo = slopeStart + span / 4;
                const int mid_hi = slopeStart + (3 * span) / 4;

                int bestIdx = -1;
                double bestSlope = -1e300;   // want least-negative (highest)
                for (int i = std::max(1, mid_lo); i <= std::min(N - 2, mid_hi); ++i) {
                    if (std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                    const double slope = 0.5 * (v[i + 1] - v[i - 1]);   // central diff
                    if (slope > bestSlope) { bestSlope = slope; bestIdx = i; }
                }

                if (bestIdx >= 0 && bestSlope >= 0.0) {
                    // Genuine notch (slope actually turned non-negative here).
                    notchIdx = bestIdx;
                }
                else {
                    // No real inflection -- placeholder O at the middle of
                    // the downslope.
                    notchOIdx = (slopeStart + end) / 2;
                }
            }
        }

        int p2 = -1;
        if (notchIdx >= 0 && end > notchIdx + 1) {
            double best = -1e300;
            for (int i = notchIdx + 1; i < end; ++i) {
                if (std::isnan(v[i]) || std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                if (v[i] >= v[i - 1] && v[i] >= v[i + 1] && v[i] > best) {
                    best = v[i]; p2 = i;
                }
            }
        }

        m_glyphs.ppgFoot = foot;
        m_glyphs.ppgP1 = p1;
        m_glyphs.ppgDic = notchIdx;
        m_glyphs.ppgP2 = p2;
        m_glyphs.ppgNotch = (notchIdx >= 0);
        m_glyphs.ppgNoNotchO = notchOIdx;

        fprintf(stderr,
            "[glyph] N=%d foot=%d p1=%d dic=%d p2=%d noNotchO=%d end=%d\n",
            N, foot, p1, notchIdx, p2, m_glyphs.ppgNoNotchO, end);
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
        // O-glyph at expected-but-not-found positions (ECG version, on the
        // ECG axis). Same visual style as the PPG side's O.
        auto ecgCirc = [&](int idx) {
            if (idx < 0 || idx >= N) return;
            const double val = std::isnan(v[idx]) ? baseline : v[idx];
            const double x = xFromSample(idx, /*isEcg=*/true);
            const double y = plotY(val, yLo, yHi);
            p.setBrush(Qt::NoBrush);
            p.setPen(QPen(Qt::black, 1.8));
            p.drawEllipse(QPointF(x, y), 4.0, 4.0);
            };
        // P/R/T peaks: X when a real extremum was found, else O at the
        // fallback (window midpoint).
        if (m_glyphs.ecgP >= 0) g(m_glyphs.ecgP);
        else                    ecgCirc(m_glyphs.ecgPOFallback);
        if (m_glyphs.ecgR >= 0) g(m_glyphs.ecgR);
        else                    ecgCirc(m_glyphs.ecgROFallback);
        if (m_glyphs.ecgT >= 0) g(m_glyphs.ecgT);
        else                    ecgCirc(m_glyphs.ecgTOFallback);
        // Q onset, S end, T end: always X (movable-marker positions;
        // unset markers naturally hit the idx<0 guard in g()).
        g(m_glyphs.ecgQ); g(m_glyphs.ecgS); g(m_glyphs.ecgTend);
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
        // Four glyphs: foot, P1, dicrotic (X if found, O at descent midpoint
        // if not), P2. Each only when its own value is available (idx >= 0).
        g(m_glyphs.ppgFoot);
        g(m_glyphs.ppgP1);
        if (m_glyphs.ppgNotch) g(m_glyphs.ppgDic);
        else                   circ(m_glyphs.ppgNoNotchO);
        g(m_glyphs.ppgP2);
    }
}