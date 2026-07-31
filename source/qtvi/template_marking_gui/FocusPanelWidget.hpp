#pragma once
//
// FocusPanelWidget.hpp  (B2 focus mode)
//
// When the operator selects a landmark in a BinPlotWidget, this panel
// renders that landmark's anchored average, zoomed in around the landmark
// column, with:
//   - the mean trace (center line),
//   - the fitted curve (from anchor_fit::selectAnchorModel over the zoom
//     window), and
//   - a 95% confidence band: mean +/- 1.96 * se per column, where
//     se = sd / sqrt(nBeats).
//
// The panel reads the template's per-sample mean + sd + beat count directly
// from the data the viewer already holds (TemplateBin.ecgTemplate_raw /
// _raw_iqr [which holds STD, ddof=1] / n_beats_raw, and the PPG analogues),
// so it needs no re-anchoring or beat matrix of its own.
//
// The J-point (S-end) is shared between the QRS and JT anchored views (spec
// I-4 already models J-point as a single AnchorType). Because this panel is
// driven by (marker, mean, sd, n), the owner can refresh BOTH the QRS and JT
// focus panels from the same J-point edit simply by calling setFocus() on
// each with the appropriate window -- see TemplateViewerWindow wiring.
//
#include <QWidget>
#include <QString>
#include <vector>

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
        int halfWindowSamples = 30,
        int framingBias = 0);

    // Clear the panel (no landmark selected).
    void clearFocus();

protected:
    void paintEvent(QPaintEvent*) override;

private:
    std::vector<double> m_mean;
    std::vector<double> m_sd;
    int    m_nBeats = 0;
    int    m_landmarkCol = -1;
    int    m_half = 30;
    int    m_framingBias = 0;   // -1 right-edge, +1 left-edge, 0 centered
    QString m_label;
    bool   m_active = false;

    // Build the fitted curve over [lo, hi] using anchor_fit; returns a
    // per-column vector (NaN outside the fit window). Defined in the .cpp.
    std::vector<double> fittedCurve(int lo, int hi) const;
};