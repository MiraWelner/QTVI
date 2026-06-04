/**
 * @file   gap_indicator.hpp
 * @brief  Detects and renders time-gap markers on the markable signal
 *         charts (ECG1/2/3, PPG, ABP).
 *
 *         A "gap" is any jump between consecutive raw timestamps in the
 *         currently loaded chunk longer than the threshold defined in
 *         gap_indicator.cpp (kGapThresholdSec). Because the GUI plots
 *         the raw scatter at uniform index-time spacing (one sample =
 *         1/native_rate of chart x, regardless of real-time stamp), the
 *         gap itself is visually invisible -- two samples bracketing a
 *         multi-hour gap sit at consecutive x positions like any other
 *         pair. This class is what tells the user "by the way, those
 *         two adjacent dots are actually 8 hours apart on the clock."
 *
 *         For each detected gap, draws a thin vertical line at the
 *         join, plus a "<dur> s gap" text label. The bracket sits
 *         between two adjacent samples (sub-pixel width).
 *
 *         Lifetime mirrors pulse_overlay:
 *           - Constructed once with the owning GUI.
 *           - rescan() once per chunk load, BEFORE the GUI rewrites the
 *             raw vectors to index-time -- detection needs the real
 *             wall-clock timestamps still in place.
 *           - clearSeries() at the TOP of handle_data_plot(), before
 *             wipeChartContent runs.
 *           - refresh() at the END of handle_data_plot(), after the
 *             host charts have rebuilt their axes.
 *
 *         The gap_indicator is a friend of noise_marking_gui so it can
 *         read the raw (t, v) vectors and per-channel native rates
 *         directly. Window state comes through public getters.
 *
 * @author Mira Welner
 */
#pragma once

#include <QList>
#include <QVector>

class QAbstractSeries;
class noise_marking_gui;

class gap_indicator {
public:
    explicit gap_indicator(noise_marking_gui* gui);
    ~gap_indicator();

    gap_indicator(const gap_indicator&) = delete;
    gap_indicator& operator=(const gap_indicator&) = delete;

    /// Re-scan the GUI's cached raw vectors for timestamp jumps that
    /// exceed the gap threshold. MUST be called BEFORE the raw vectors
    /// are rewritten to index-time -- detection looks at the real
    /// timestamps. Cheap; one pass per chunk load.
    void rescan();

    /// Draw bracket line + duration label for every cached gap that
    /// intersects the current visible window. Call at the END of
    /// handle_data_plot().
    void refresh();

    /// Detach and delete every series this owns. Call at the TOP of
    /// handle_data_plot() before wipeChartContent runs.
    void clearSeries();

    int gapCount() const { return m_gaps.size(); }

private:
    // One gap = one chart-x position (in the same index-time coords the
    // raw vectors use after rewrite) plus the real-time duration that
    // would otherwise be hidden there.
    struct GapEntry {
        double xPos;          ///< chart x where the bracket/label lives
        double durationSec;   ///< real wall-clock gap duration
    };

    noise_marking_gui* m_gui = nullptr;
    QVector<GapEntry>        m_gaps;
    QList<QAbstractSeries*>  m_series;     ///< owned, recreated per refresh
};