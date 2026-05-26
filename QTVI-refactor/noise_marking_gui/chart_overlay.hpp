/**
 * @file   chart_overlay.hpp
 * @brief  Static red grid overlay for the markable signal charts.
 *
 *         Draws thin translucent red vertical lines across the visible
 *         window of every markable channel (ECG1/2/3, PPG, ABP):
 *           - a minor line every 0.1 s
 *           - a thicker major line every 0.5 s
 *
 *         Owns two QLineSeries per channel (one for minor lines, one
 *         for major lines). Each series packs N vertical segments into
 *         one polyline using NaN breakpoints so we don't need N series.
 *         The series are detached and re-attached on every refresh()
 *         because the host charts rebuild their axes on every redraw.
 */
#pragma once

#include <QObject>
#include <QList>
#include <QMap>
#include <QPen>
#include <QString>
#include <QStringList>

class QLineSeries;
class QChartView;

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
     *        and rebuild the vertical-line geometry for the current
     *        window. Must be called at the END of handle_data_plot(),
     *        after axes have been rebuilt. No-op while disabled.
     */
    void refresh();

    /**
     * @brief Return the overlay's series for `label` (both minor and
     *        major), or an empty list if no such label exists. The
     *        caller does NOT own the returned pointers.
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
    // Per-channel pair: one minor (0.1 s) series, one major (0.5 s) series.
    struct ChannelSeries {
        QLineSeries* minor = nullptr;
        QLineSeries* major = nullptr;
    };

    noise_marking_gui* m_gui = nullptr;
    QStringList                    m_labels;
    QMap<QString, ChannelSeries>   m_series;
    QPen                           m_minorPen;
    QPen                           m_majorPen;
    bool                           m_enabled = true;
};