#pragma once
//
// global_interval_lines.hpp
//
// Draws the GLOBAL interval boundaries as vertical reference lines across the
// per-lead template panels, so the earliest onset and latest offset found in
// ANY lead are visible in EVERY lead.
//
// Self-contained on purpose: BinPlotWidget only has to store a
// std::vector<Line> and make one paint() call, so nothing about the boundary
// logic lives inside the widget's already-large paintEvent.
//
// ---------------------------------------------------------------------------
// THE CONVERSION THAT MAKES THIS CORRECT
// ---------------------------------------------------------------------------
// GlobalIntervals reports its boundaries as R-RELATIVE sample offsets, because
// each channel is aligned independently (its own r_col_raw) and therefore has
// its own column space. Drawing the same column number in all three panels
// would put the line at a different physical instant in each one, which is
// exactly the error the global measurement exists to avoid.
//
// forChannel() does the one conversion that fixes it: column = offset + THIS
// channel's own R column. Same instant, three different column numbers, lines
// that visually align across the stacked panels.
//
// The lines are drawn dashed and thin, behind the draggable marker bars. They
// are read-only context imported from the other leads -- if they looked like
// this panel's own markers, a user would try to drag one.
//

#include "global_intervals.hpp"
#include "template_marking_bin_io.hpp"   // TemplateBin

#include <QColor>
#include <QPainter>
#include <QPointF>
#include <QString>

#include <cmath>
#include <vector>

namespace global_interval_lines {

    /// One vertical line, positioned in the SAMPLE COLUMNS of the panel that
    /// will draw it (already converted out of the R-relative axis).
    struct Line {
        double  column = -1.0;
        QColor  color{ 120, 120, 200 };
        QString label;          ///< short tag drawn at the top ("QRS on" / "QRS off")
    };

    inline QColor onsetColor() { return QColor(90, 110, 210); }
    inline QColor offsetColor() { return QColor(150, 90, 190); }

    /// That channel's own R column, sub-sample where available. -1 if absent.
    inline double rColumnFor(const TemplateBin& b, int ch) {
        if (ch < 0 || ch >= 3) return -1.0;
        return (b.r_peak_auto_ch[ch] >= 0.0) ? b.r_peak_auto_ch[ch]
            : static_cast<double>(b.r_peak_ch[ch]);
    }

    /**
     * @brief The boundary lines for ONE channel's panel.
     *
     * @return Empty when the intervals are not established or this channel has
     *         no R column -- with no R there is no way to place a shared-axis
     *         offset in this channel's columns, and a line drawn at the raw
     *         offset would land near sample 0 and look like a real boundary.
     */
    inline std::vector<Line> forChannel(const TemplateBin& b,  const global_intervals::GlobalIntervals& g, int ch) {
        //creates the dotted line boundries representing the min onset and max offset of QRS for once plot box
        std::vector<Line> out;
        if (!g.valid) return out;
        const double rc = rColumnFor(b, ch);
        if (rc < 0.0) return out;

        if (!std::isnan(g.qrsOnset)) {
            Line l;
            l.column = g.qrsOnset + rc;
            l.color = onsetColor();
            out.push_back(l);
        }
        if (!std::isnan(g.qrsOffset)) {
            Line l;
            l.column = g.qrsOffset + rc;
            l.color = offsetColor();
            out.push_back(l);
        }
        return out;
    }

    /**
     * @brief Paint the lines.
     *
     * @param xForColumn  Maps a sample column to a pixel x. Pass the widget's
     *                    own mapping (BinPlotWidget::xFromSample(s, 1.0)) so
     *                    the lines cannot drift from the trace when the panel
     *                    resizes -- one mapping, not two.
     * @param nEcgSamples Length of this panel's ECG trace. A boundary set by
     *                    another lead can fall outside this template's extent;
     *                    such lines are skipped rather than clamped to the
     *                    frame edge, where they would imply a boundary that
     *                    was measured here.
     */
    template <typename XMap>
    inline void paint(QPainter& p,
        const std::vector<Line>& lines,
        XMap xForColumn,
        int topPx, int bottomPx,
        int nEcgSamples) {
        if (lines.empty() || nEcgSamples < 2) return;

        p.save();
        for (const Line& l : lines) {
            if (l.column < 0.0 || l.column > nEcgSamples - 1) continue;

            const double x = xForColumn(l.column);

            QPen pen(l.color);
            pen.setWidthF(0.7);
            pen.setStyle(Qt::DotLine);
            p.setPen(pen);
            p.drawLine(QPointF(x, topPx), QPointF(x, bottomPx));

            if (!l.label.isEmpty()) {
                QFont f = p.font();
                f.setPointSize(5);
                p.setFont(f);
                p.setPen(l.color);
                p.drawText(QPointF(x + 2.0, topPx + 9.0), l.label);
            }
        }
        p.restore();
    }

}  // namespace global_interval_lines#pragma once
