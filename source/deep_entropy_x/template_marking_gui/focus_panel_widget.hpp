#pragma once
//
// FocusPanelWidget.hpp  (B2 focus mode)
//
// When the operator selects a landmark in a BinPlotWidget, this panel
// renders that landmark's anchored average, zoomed in around the landmark
// column, with:
//   - the mean trace (center line),
//   - the fitted curve (from curve_fit::selectBestFit over the zoom
//     window), and
//   - a 95% confidence band: mean +/- 1.96 * se per column, where
//     se = sd / sqrt(nBeats).
//
// The panel reads the template's per-sample mean + sd + beat count directly
// from the data the viewer already holds (TemplateBin.ecgTemplate_raw /
// _raw_iqr [which holds STD, ddof=1] / n_beats_raw, and the PPG analogues),
// so it needs no re-anchoring or beat matrix of its own.
//
// ONE PANEL PER SELECTION. There used to be two -- a QRS view and a JT view --
// with the J-point drawn in both, framed to opposite edges, because with a
// single alignment on screen neither panel alone could show both of its
// neighbourhoods. TemplateViewerWindow now switches this panel's WAVEFORM to
// the clicked landmark's own alignment (anchor_view::anchorFor), so there is
// one view per landmark and nothing to keep in step.
//
// This class needed no change for that: it is driven entirely by
// (marker, mean, sd, n, column), so which alignment the mean came from is the
// caller's business.
//
#include <QWidget>
#include <QString>
#include <cstdint>
#include <limits>
#include <vector>
#include "subsample_refine.hpp"   // subsample_refine::TransitionCandidates

class QPainter;

class FocusPanelWidget : public QWidget {
    Q_OBJECT
public:
    explicit FocusPanelWidget(QWidget* parent = nullptr);

    // Point the panel at one landmark of one anchored average.
    //   mean    : per-column center line (the template mean).
    //   sd      : per-column standard deviation (ddof=1) across contributing
    //             beats -- same length as mean.
    //   nBeats  : number of beats contributing (the n in se = sd/sqrt(n)).
    //   landmarkCol : the column the selected landmark sits on; the zoom
    //                 window is centered here.
    //   label   : shown in the panel header (e.g. "R peak", "J-point (QRS)").
    //   halfWindowSamples : zoom half-width in samples around landmarkCol.
    //   framingBias : -1 frames the landmark toward the RIGHT edge (it ENDS
    //                 this segment, e.g. J-point as QRS-end), +1 toward the
    //                 LEFT edge (it STARTS this segment, e.g. J-point as
    //                 JT-start), 0 = centered (default).
    void setFocus(const std::vector<double>& mean,
        const std::vector<double>& sd,
        int nBeats,
        int landmarkCol,
        const QString& label,
        int halfWindowSamples = 100,
        int framingBias = 0);

    // Per-sample SD in MSEC: each column's amplitude SD divided by the
    // template's local |dV/dt| there (localAbsSlope / slopeFloor in
    // TemplateViewerWindow.cpp). floorMask marks the columns where the slope
    // was clamped at the floor -- flat regions, peak tops -- where the value
    // is a lower bound rather than a measurement, and which are shaded.
    //
    // Display only: the band above is unchanged (it is the amplitude CI).
    // Empty vectors leave the panel exactly as it was.
    void setSdMs(const std::vector<double>& sdMs,
        const std::vector<uint8_t>& floorMask,
        const std::vector<double>& deriv,
        double slopeFloor);

    // Clear the panel (no landmark selected).
    void clearFocus();

    // Which curve the red dashed overlay draws, so it matches the model that
    // actually PLACED this landmark: a peak's weighted quadratic/cubic vs an
    // onset/offset's curve_fit transition model. Set per landmark by the owner
    // right after setFocus; defaults to Transition (the onset/offset case).
    enum class FitKind { Transition, PeakQuadratic, PeakCubic };
    // peakSigma is the SAME weighting the detector used for this landmark (so
    // the drawn fit and the placement fit can't diverge); ignored for
    // transitions. Set per landmark by the owner right after setFocus.
    void setFitKind(FitKind k, double peakSigma = 4.0) {
        m_fitKind = k; m_peakSigma = peakSigma; update();
    }
    // Forced peak model (Fit-Peaks radio); Auto = BIC quad-vs-cubic.
    void setPeakFitMode(curve_fit::PeakFitMode m) { m_panelPeakMode = m; update(); }

    // Supply the EXACT candidate curves the detector fit for this transition
    // landmark (sample-indexed closures + winner), so the panel draws the fits
    // that placed the mark rather than a re-fit over the visible window. Cleared
    // by clearFocus; ignored for peaks.
    void setTransitionCandidates(const subsample_refine::TransitionCandidates& c) {
        m_transCands = c; update();
    }

    // One tested model curve, plus whether the selector chose it. Drawn green
    // when selected, red otherwise.
    struct Candidate {
        std::vector<double> curve;   // per-column, NaN outside the fit
        bool selected = false;
        QString label;               // model name (shown for the winner)
        // WHERE THIS MODEL PUTS THE FIDUCIAL, in sub-sample columns; NaN if it
        // has no placement to report. The dotted line is drawn at the SELECTED
        // candidate's position, which is what makes it follow the fit-model
        // radios -- it used to be pinned to the integer bar column and so never
        // responded to them at all.
        double position = std::numeric_limits<double>::quiet_NaN();
    };

protected:
    void paintEvent(QPaintEvent*) override;

private:
    FitKind m_fitKind = FitKind::Transition;
    double  m_peakSigma = 4.0;
    curve_fit::PeakFitMode m_panelPeakMode = curve_fit::PeakFitMode::Auto;
    subsample_refine::TransitionCandidates m_transCands;   // supplied exact transition fits
    std::vector<double> m_mean;
    std::vector<double> m_sd;
    std::vector<double>  m_sdMs;       // per-column SD in msec
    std::vector<double>  m_deriv;      // per-column |dV/dt| (SG derivative), amp/sample
    double  m_slopeFloor = 0.0;        // denominator used where the slope is below it
    std::vector<uint8_t> m_floorMask;  // 1 where the slope floor engaged
    int    m_nBeats = 0;
    int    m_landmarkCol = -1;
    // Sub-sample column the fiducial was last DRAWN at. Recorded in paintEvent
    // for the footer readout: quadratic and cubic vertices often differ by well
    // under one sample, which at this zoom is a sub-pixel shift, so the number
    // is the only reliable way to see that the mark moved with the radio.
    double m_lastFidCol = -1.0;
    int    m_half = 30;
    int    m_framingBias = 0;   // -1 right-edge, +1 left-edge, 0 centered
    QString m_label;
    bool   m_active = false;

    // Peak candidate curve over [lo, hi]: a weighted quadratic (cubic=false) or
    // cubic (cubic=true). NaN if the fit degenerated. (Retained; candidateCurves
    // is the live path.)

    // Every model the DETECTOR tested for this landmark, over [lo, hi], with
    // the selector's winner flagged. Peaks: quadratic + cubic (BIC). Onsets/
    // offsets: piecewise-linear + sigmoid + fractional (selectBestFit).
    // Uses the same fitting functions the detector uses, so the drawn curves
    // are the tested curves.
    std::vector<Candidate> candidateCurves(int lo, int hi) const;
};