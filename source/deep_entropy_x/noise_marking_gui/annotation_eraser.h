/**
 * @file   annotation_eraser.h
 * @brief  Right-click-to-delete for annotated regions on the markable charts.
 *
 *         A right-click inside a highlighted annotation removes that annotation
 *         (the most-recently-added one when several overlap, so repeated clicks
 *         peel a stack off one at a time, topmost first). The marking list
 *         (m_genExc) is the source of truth, exactly as undo/clear treat it:
 *         the matching entry is dropped, m_noiseManager is rebuilt from the
 *         survivors, the now-stale logged beats in the erased span are removed
 *         (the redraw re-detects/re-logs whatever is actually there now), and
 *         the dialog is redrawn.
 *
 *         Friend of noise_marking_gui so it can edit m_genExc / m_noiseManager
 *         and trigger a redraw, mirroring gap_indicator's access.
 */
#pragma once

class QChartView;
class QPoint;
class noise_marking_gui;

class annotation_eraser {
public:
    explicit annotation_eraser(noise_marking_gui* gui);

    annotation_eraser(const annotation_eraser&) = delete;
    annotation_eraser& operator=(const annotation_eraser&) = delete;

    /// Delete the annotation under a right-click. @p viewportPos is the click
    /// position in the chart view's viewport coordinates (as delivered to
    /// eventFilter). Returns true if an annotation was erased (the caller
    /// should then treat the mouse event as handled).
    bool handleRightClick(QChartView* cv, const QPoint& viewportPos);

private:
    noise_marking_gui* m_gui = nullptr;
};