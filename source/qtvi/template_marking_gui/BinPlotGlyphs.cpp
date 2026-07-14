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

    // ---- PPG landmarks (all driven by the anchor-pulse markers set by
    // seedBinMarkers, which restrict search to the R-window; the global
    // argmax used previously would land on the SECOND pulse under Patch B
    // slicing, then cascade every other landmark off the anchor pulse). ----
    if (m_hasPPG && (int)m_ppg.size() >= 3) {
        const std::vector<double>& v = m_ppg;
        const int N = (int)v.size();
        const int foot = m_markers[PpgOnset];      // anchor foot
        const int p1 = m_markers[PpgPeak];       // anchor systolic peak
        const int end = m_markers[PpgEnd];        // anchor end (trough after)

        // Notch test: any interior local minimum on the downstroke [p1, end].
        bool hasNotch = false;
        int  notchIdx = -1;
        if (p1 >= 0 && end > p1) {
            const int lo = std::max(0, p1 + 1), hi = std::min(N - 1, end);
            double best = 1e300;
            for (int i = lo; i < hi; ++i) {
                if (std::isnan(v[i]) || std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                if (v[i] <= v[i - 1] && v[i] <= v[i + 1] && v[i] < best) {
                    best = v[i]; notchIdx = i; hasNotch = true;
                }
            }
        }

        // Diastolic peak: local max between the notch and the end. Only if
        // a real local max exists; no fallback (see previous edit).
        int p2 = -1;
        if (hasNotch && end > notchIdx) {
            const int lo = std::max(0, notchIdx + 1), hi = std::min(N - 1, end);
            double best = -1e300;
            for (int i = lo + 1; i < hi; ++i) {
                if (std::isnan(v[i])) continue;
                if (v[i] >= v[i - 1] && v[i] >= v[i + 1] && v[i] > best) {
                    best = v[i]; p2 = i;
                }
            }
        }

        m_glyphs.ppgFoot = foot;
        m_glyphs.ppgP1 = p1;
        m_glyphs.ppgDic = notchIdx;        // -1 if none; draw branch handles that
        m_glyphs.ppgP2 = p2;
        m_glyphs.ppgNotch = hasNotch;
        m_glyphs.ppgNoNotchO = (!hasNotch && p1 >= 0 && end > p1)
            ? (p1 + end) / 2 : -1;

        fprintf(stderr,
            "[glyph] ppg: N=%d foot=%d p1=%d end=%d hasNotch=%d notch=%d "
            "p2=%d noNotchO=%d\n",
            N, foot, p1, end, (int)hasNotch, notchIdx, p2,
            m_glyphs.ppgNoNotchO);
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
        const int ppgClip = pulseClipN();
        auto g = [&](int idx) {
            if (idx < 0 || idx >= N || idx >= ppgClip) return;
            const double val = std::isnan(v[idx]) ? pLo : v[idx];
            glyph(xFromSample(idx, /*isEcg=*/false), plotY(val, pLo, pHi));
            };
        auto circ = [&](int idx) {
            if (idx < 0 || idx >= N || idx >= ppgClip) return;
            const double val = std::isnan(v[idx]) ? pLo : v[idx];
            const double x = xFromSample(idx, /*isEcg=*/false);
            const double y = plotY(val, pLo, pHi);
            p.setBrush(Qt::NoBrush);
            p.setPen(QPen(Qt::black, 1.8));
            p.drawEllipse(QPointF(x, y), 4.0, 4.0);
            };
        g(m_glyphs.ppgFoot);
        g(m_glyphs.ppgP1);
        if (m_glyphs.ppgNotch) g(m_glyphs.ppgDic);
        else                   circ(m_glyphs.ppgNoNotchO);
        g(m_glyphs.ppgP2);
    }
}