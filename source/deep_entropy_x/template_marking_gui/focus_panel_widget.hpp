#pragma once
//
// focus_panel_widget.hpp
//
// The landmark close-up. When the operator selects a landmark in a
// BinPlotWidget, this panel draws that landmark's anchored average zoomed in
// around its column: the mean trace, a +/- 1 sd band, the candidate curves the
// DETECTOR fit, and a dotted line at the detector's own placement.
//
// THE PANEL FITS NOTHING AND DETECTS NOTHING. Everything it draws arrives
// through the setters below, from whoever ran the detector. That is the only
// invariant in this file that matters: a curve this panel produced itself would
// answer a different question from the one that placed the mark (different
// window, different weighting, different vertex) and would be drawn beside it.
//
// It is driven entirely by (marker, mean, sd, n, column), so which alignment
// the mean came from, and which channel, are the caller's business.
//
#include <QWidget>
#include <QString>
#include <cstdint>
#include <limits>
#include <vector>
#include <algorithm>   // std::max in setPeakFitKind
#include "subsample_refine.hpp"   // PeakCandidates / TransitionCandidates

class QPainter;

class FocusPanelWidget : public QWidget {
    Q_OBJECT
public:
    explicit FocusPanelWidget(QWidget* parent = nullptr);

    // Point the panel at one landmark of one anchored average. Clears every
    // derived array and every supplied fit, so a caller that supplies none
    // cannot inherit the previous landmark's.
    //   mean    : per-column center line.
    //   sd      : per-column spread, SAME LENGTH AS MEAN or longer; a shorter
    //             one disables the band.
    //   nBeats  : beats contributing, for the header readout.
    //   landmarkCol : zoom window centre.
    //   label   : panel header, e.g. "R peak [P-aligned]".
    //   framingBias : -1 frames the landmark toward the RIGHT edge (it ENDS
    //                 this segment), +1 toward the LEFT (it STARTS it),
    //                 0 centered.
    void setFocus(const std::vector<double>& mean,
        const std::vector<double>& sd,
        int nBeats,
        int landmarkCol,
        const QString& label,
        int halfWindowSamples = 100,
        int framingBias = 0);

    // Per-sample SD in MSEC: the amplitude SD divided by the template's local
    // |dV/dt|. floorMask marks columns where the slope was clamped at the floor
    // -- flat regions, peak tops -- where the value is a lower bound rather
    // than a measurement, and which are shaded. Display only; the band above is
    // the amplitude spread and is unaffected. Call after setFocus.
    void setSdMs(const std::vector<double>& sdMs,
        const std::vector<uint8_t>& floorMask,
        const std::vector<double>& deriv,
        double slopeFloor);

    // Clear the panel (no landmark selected).
    void clearFocus();

    // The detector's position for the focused landmark, in this trace's
    // columns. The dotted fiducial is drawn here and nowhere else. Call after
    // setFocus; -1 = not supplied, and the line falls back to the bar column.
    void setDetectorFiducial(double col);

    // Which family of curve placed this landmark, so the overlay matches the
    // model the detector actually used. None = NOT PLACED BY A FIT at all: an
    // ECG Q onset with no trough (compute_q_onset's R-upstroke fallback, drawn
    // as a circle rather than an X on the grid), or any pulse landmark that is
    // a bracketed search or an interpolated crossing. Drawing candidates there
    // would claim a contest that never happened.
    enum class FitKind { None, Transition, PeakQuadratic, PeakCubic };

    // TWO OVERLOADS, NO DEFAULTS. Transition and None have no peak parameters
    // to supply and requiring them would force the caller to invent values;
    // m_peakHalfWidth goes to -1 here so a width left from a previous selection
    // can never be read.
    void setFitKind(FitKind k) {
        m_fitKind = k;
        m_peakHalfWidth = -1;
        update();
    }

    // Peaks: both required, and both must be the values the detector fitted
    // with (subsample_refine::peak_sigma / peak_halfwidth for ECG,
    // pulse_sigma / pulse_halfwidth for pulse). A different span here draws a
    // curve over a fit nothing else in the system used.
    void setPeakFitKind(FitKind k, double peakSigma, int peakHalfWidth) {
        m_fitKind = k;
        m_peakSigma = peakSigma;
        m_peakHalfWidth = std::max(3, peakHalfWidth);
        update();
    }

    // The exact fits the detector ran. Required for anything to be drawn: the
    // panel cannot reproduce them from sigma and half-width alone, because it
    // does not know the integer seed the detector fitted around. Ignored for
    // the other family; cleared by setFocus and clearFocus.
    void setPeakCandidates(const subsample_refine::PeakCandidates& c) {
        m_peakCands = c; update();
    }
    void setTransitionCandidates(const subsample_refine::TransitionCandidates& c) {
        m_transCands = c; update();
    }

    // One tested model curve. Green when it is the one that placed the mark,
    // red otherwise.
    struct Candidate {
        std::vector<double> curve;   // per-column, NaN outside the fit
        bool selected = false;
        QString label;               // model name (shown for the winner)
        // Where THIS model puts the fiducial, in sub-sample columns; NaN if it
        // has no placement to report.
        double position = std::numeric_limits<double>::quiet_NaN();
    };

protected:
    void paintEvent(QPaintEvent*) override;

private:
    FitKind m_fitKind = FitKind::Transition;
    double  m_peakSigma = 4.0;
    int     m_peakHalfWidth = -1;   // -1 = no peak in focus
    subsample_refine::PeakCandidates       m_peakCands;
    subsample_refine::TransitionCandidates m_transCands;
    std::vector<double>  m_mean;
    std::vector<double>  m_sd;
    std::vector<double>  m_sdMs;       // per-column SD in msec
    std::vector<double>  m_deriv;      // per-column |dV/dt|, amp/sample
    double  m_slopeFloor = 0.0;        // denominator used where the slope is below it
    std::vector<uint8_t> m_floorMask;  // 1 where the slope floor engaged
    int    m_nBeats = 0;
    int    m_landmarkCol = -1;
    // Sub-sample column the fiducial was last DRAWN at, recorded in paintEvent
    // for the footer readout. Quadratic and cubic vertices can differ by well
    // under one sample -- a sub-pixel shift at this zoom -- so the number is
    // the only reliable way to see the mark move with the fit-model radio.
    double m_lastFidCol = -1.0;
    double m_detectorFid = -1.0;
    int    m_half = 30;
    int    m_framingBias = 0;   // -1 right-edge, +1 left-edge, 0 centered
    QString m_label;
    bool   m_active = false;

    // Every model the DETECTOR tested for this landmark over [lo, hi], winner
    // flagged. Assembled from the supplied candidates; fits nothing.
    std::vector<Candidate> candidateCurves(int lo, int hi) const;
};