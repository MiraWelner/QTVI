#pragma once
//
// focus_panel_widget.hpp
//
// The landmark close-up. When the operator selects a landmark in a
// BinPlotWidget, this panel draws that landmark's anchored average zoomed in
// around its column: the mean trace, a +/- 1 sd band, the candidate curves the
// DETECTOR fit, and two vertical fiducial lines.

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
    //   landmarkCol : the landmark's own column, and the zoom window centre.
    //             NOT where the thick bar line is drawn -- that is
    //             setUserFiducial, which a glyph never gets.
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


    void setLandmarkCol(int col) {
        if (col == m_landmarkCol) return;   // a repaint would show the same thing
        m_landmarkCol = col;
        if (m_userFid >= 0.0) m_userFid = static_cast<double>(col);
        m_active = (col >= 0 && !m_mean.empty());
        update();
    }
    int landmarkCol() const { return m_landmarkCol; }

    // The detector's position for the focused landmark, in this trace's
    // columns. The dotted AUTO line is drawn here and nowhere else. Call after
    // setFocus.
    //
    // -1 = NOT SUPPLIED, and no auto line is drawn. It does NOT mean "use the
    // landmark column": a landmark with no detector placement has no auto
    // position to report, and drawing one at the bar would attribute the
    // operator's choice to the detector.
    void setDetectorFiducial(double col);

    // The OPERATOR's position for the focused landmark, in this trace's
    // columns. The THICK BAR LINE is drawn here and nowhere else, and it is
    // drawn for as long as the focus is on a MOVABLE bar -- not just during
    // the drag. The line is what says "this landmark is yours to move and it
    // is currently here"; the dotted auto line beside it is what the detector
    // said. Both readings belong on screen the whole time the bar is selected,
    // and the gap between them is the edit.
    //
    // -1 = THIS LANDMARK CANNOT HAVE ONE, and no bar line is drawn. Every
    // auto-only glyph (R/P/Q/T peak, the pulse T50/peak/peak2/T80) lands here:
    // its column is the detector's answer, so a bar drawn at it would
    // attribute a detection to the operator, and the footer's "bar =" would
    // report an operator placement that was never made.
    //
    // NOTHING WAS CALLING THIS, which is why the line never appeared at all:
    // m_userFid stayed at -1 for the life of the panel and the paint block
    // guarded by it was dead code. setLandmarkCol could not bootstrap one
    // either -- it only TRACKS a fiducial that is already set.
    void setUserFiducial(double col);
    double userFiducial() const { return m_userFid; }
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
    // pulse_window::{peak,foot}Halfwidth + ::sigma for pulse, resolved
    // against that channel's rate). A different span here draws a curve over
    // a fit nothing else in the system used.
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
    void setPeakCandidates(const upsample_for_fit::PeakCandidates& c) {
        m_peakCands = c; update();
    }
    void setTransitionCandidates(const upsample_for_fit::TransitionCandidates& c) {
        m_transCands = c; update();
    }

    // One tested model curve. Solid green when it is the one that placed the
    // mark, solid gray otherwise.
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
    upsample_for_fit::PeakCandidates       m_peakCands;
    upsample_for_fit::TransitionCandidates m_transCands;
    std::vector<double>  m_mean;
    std::vector<double>  m_sd;
    std::vector<double>  m_sdMs;       // per-column SD in msec
    std::vector<double>  m_deriv;      // per-column |dV/dt|, amp/sample
    double  m_slopeFloor = 0.0;        // denominator used where the slope is below it
    std::vector<uint8_t> m_floorMask;  // 1 where the slope floor engaged
    int    m_nBeats = 0;
    // The landmark's own column: the zoom window centre and the index the
    // footer's sd lookups use. NOT the thick line's position.
    int    m_landmarkCol = -1;
    // The operator's own position, or -1 for a landmark that cannot have one.
    // The thick line is drawn here and nowhere else.
    double m_userFid = -1.0;
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