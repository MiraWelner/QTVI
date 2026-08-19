/**
 * @file   grid_overlay.hpp
 * @brief  Static red grid overlay for the markable signal charts.
 *
 *         Draws thin translucent red grid lines across the visible
 *         window of every markable channel (ECG1/2/3, PPG, ABP). The
 *         x-axis is seconds; the y-axis is whatever the channel is
 *         (mV, ADC counts, etc) so the two spacings cannot be equal
 *         in data units. Instead, the y-spacing is chosen at refresh
 *         time so that one grid cell renders as a visual square -- the
 *         number of y-units per pixel and the number of seconds per
 *         pixel are computed from the chart's plot area, and the
 *         y-spacing is scaled to match the x-spacing in pixels.
 *
 *         Lines drawn:
 *           - vertical:   minor every 0.04 s, major every 0.2 s
 *           - horizontal: spaced to match those in pixel terms (i.e.
 *                         squares; spacing in y-units recomputed on
 *                         every refresh from plotArea() + axis ranges)
 *
 *         Owns four QLineSeries per channel (vertical minor, vertical
 *         major, horizontal minor, horizontal major). Each series packs
 *         many segments into one polyline using NaN breakpoints so we
 *         don't need one series per line. All four are detached and
 *         re-attached on every refresh() because the host charts
 *         rebuild their axes on every redraw.
 *
 *         Because the squares depend on the chart's pixel dimensions
 *         and on the y-axis range, the overlay also refreshes itself
 *         on chart resize and on y-axis range-change. Both happen
 *         frequently (window resize, channel autoscale) and would
 *         otherwise leave the grid stretched into rectangles until
 *         the next handle_data_plot().
 */
#pragma once

#include <QObject>
#include <QList>
#include <QMap>
#include <QPen>
#include <QPointer>
#include <QString>
#include <QStringList>

class QLineSeries;
class QChartView;
class QValueAxis;

class noise_marking_gui;

class pulse_overlay : public QObject {
    Q_OBJECT
public:
    /**
     * @param gui     Owning GUI; used to look up chart views, axes, and
     *                window bounds via its public/friend interface.
     * @param labels  Channels to overlay (matches markableChannelLabels()).
     */
    pulse_overlay(noise_marking_gui* gui, const QStringList& labels);
    ~pulse_overlay() override;

    /**
     * @brief Re-attach the overlay series to each channel's current chart
     *        and rebuild the grid geometry (vertical AND horizontal) for
     *        the current window. Must be called at the END of
     *        handle_data_plot(), after axes have been rebuilt. No-op
     *        while disabled.
     */
    void refresh();

    /**
     * @brief Return the overlay's series for `label` (vertical and
     *        horizontal, minor and major -- 4 entries), or an empty
     *        list if no such label exists. The caller does NOT own
     *        the returned pointers.
     *
     *        Used by handle_data_plot()'s wipeChartContent loop to keep
     *        the overlay's series alive across a chart wipe -- the
     *        overlay owns them and deletes them in its own destructor.
     */
    QList<QLineSeries*> seriesForLabel(const QString& label) const;

    /**
     * @brief Enable or disable the overlay. When disabled, all series
     *        are detached from their charts; refresh() is a no-op.
     *        Re-enabling re-attaches on the next refresh().
     */
    void setEnabled(bool on);

    bool isEnabled() const { return m_enabled; }

private:
    // Per-channel grid: two vertical series (minor + major time grid)
    // and two horizontal series (minor + major y-grid). All four are
    // built fresh on every refresh() from the chart's pixel-space
    // dimensions so cells stay square through resizes and rescales.
    struct ChannelSeries {
        QLineSeries* vMinor = nullptr;
        QLineSeries* vMajor = nullptr;
        QLineSeries* hMinor = nullptr;
        QLineSeries* hMajor = nullptr;
    };

    // Wire up the chart-resize and y-axis range-change signals for
    // `label` so the overlay re-squares itself outside the normal
    // handle_data_plot() path. Safe to call repeatedly: connections
    // are guarded by m_hookedCharts / m_hookedAxes so we don't double-
    // subscribe when refresh() rediscovers the same chart.
    void hookChartIfNeeded(const QString& label, QChartView* cv,
        QValueAxis* yAxis);

    noise_marking_gui* m_gui = nullptr;
    QStringList                    m_labels;
    QMap<QString, ChannelSeries>   m_series;
    QPen                           m_minorPen;
    QPen                           m_majorPen;
    bool                           m_enabled = true;

    // Re-entry guard. refresh() calls addSeries() / attachAxis() /
    // setSeries replace, all of which can synchronously fire
    // plotAreaChanged or rangeChanged on the chart we just hooked --
    // and our hook posts another refresh(). Without this flag, every
    // toggle of the grid would recurse until the stack blew. The
    // legitimate case (chart actually resized after refresh() returned)
    // still works because the flag clears before this call's return.
    bool                           m_inRefresh = false;

    // Tracks which charts / axes we've already hooked so the per-frame
    // refresh() doesn't connect the same signal twice every time it
    // walks the channel list. QPointer lets us notice if a chart got
    // deleted out from under us (host wipes charts on data reload).
    QList<QPointer<QObject>> m_hookedCharts;
    QList<QPointer<QObject>> m_hookedAxes;
};