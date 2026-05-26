/**
 * @file   gap_indicator.hpp
 * @brief  Detects and renders time-gap markers on the markable signal
 *         charts (ECG1/2/3, PPG, ABP).
 *
 *         A "gap" is any jump between consecutive raw timestamps in the
 *         currently loaded chunk longer than the threshold defined in
 *         gap_indicator.cpp (kGapThresholdSec). For each gap that
 *         intersects the visible window, the indicator draws two dashed
 *         gray vertical lines bracketing the gap plus a "<dur> s gap"
 *         text label.
 *
 *         Lifetime mirrors pulse_overlay:
 *           - Constructed once with the owning GUI.
 *           - rescan() once per chunk load (after raw vectors populated).
 *           - clearSeries() at the TOP of handle_data_plot(), before
 *             wipeChartContent, so the wipe doesn't see (and double-
 *             free) series that we already own.
 *           - refresh() at the END of handle_data_plot(), after the host
 *             charts have rebuilt their axes.
 *
 *         The gap_indicator is a friend of noise_marking_gui so it can
 *         read the raw (t, v) vectors directly. The view/window state
 *         comes through public getters on the GUI to keep the coupling
 *         narrow.
 *
 * @author Mira Welner
 */
#pragma once

#include <QList>
#include <QVector>
#include <QPair>

class QAbstractSeries;
class noise_marking_gui;

class gap_indicator {
public:
    explicit gap_indicator(noise_marking_gui* gui);
    ~gap_indicator();

    gap_indicator(const gap_indicator&) = delete;
    gap_indicator& operator=(const gap_indicator&) = delete;

    /// Re-scan the GUI's cached raw vectors for timestamp jumps that
    /// exceed the gap threshold. Result cached in `m_gaps`. Cheap; call
    /// once per chunk load after the raw blocks are populated.
    void rescan();

    /// Draw bracket lines + duration label for every cached gap that
    /// intersects the current visible window. Call at the END of
    /// handle_data_plot(), after axes are rebuilt by renderWindowedChart.
    void refresh();

    /// Detach and delete every series this owns. Call at the TOP of
    /// handle_data_plot() before wipeChartContent runs.
    void clearSeries();

    int gapCount() const { return m_gaps.size(); }

private:
    noise_marking_gui* m_gui = nullptr;
    QVector<QPair<double, double>> m_gaps;       ///< chunk-local seconds
    QList<QAbstractSeries*>        m_series;     ///< owned, recreated per refresh
};