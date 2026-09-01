#pragma once
/**
* @file   chart_utils.hpp
* @brief  Sets the colors of all the series and contains various features required both markable and nonmarkable plots.
*/

#include <QtCharts/QChart>
#include <QtCharts/QValueAxis>
#include "annotation_types.hpp"

#include <QString>
#include <QVector>

inline const QColor COLOR_ECG1{ 101, 67, 33 };      // dark brown
inline const QColor COLOR_ECG2{ 0, 110, 0 };        // dark green
inline const QColor COLOR_ECG3{ 0, 0, 139 };        // navy blue
inline const QColor COLOR_PPG{ 120, 20, 110 };      // dark magenta
inline const QColor COLOR_ABP{ 0, 95, 105 };        // dark teal
inline const QColor COLOR_ACCEL_X{ 170, 90, 0 };    // dark orange
inline const QColor COLOR_ACCEL_Y{ 0, 105, 80 };    // dark sea green
inline const QColor COLOR_ACCEL_Z{ 80, 20, 140 };   // dark violet
inline const QColor COLOR_RESP{ 90, 100, 0 };       // dark olive
inline const QColor COLOR_CVP{ 140, 0, 0 };         // dark red / maroon
inline const QColor COLOR_TEMP{ 140, 50, 20 };      // dark sienna
inline const QColor COLOR_MARKER{ 80, 80, 80 };     // dark gray
inline const QColor COLOR_RAW_SCATTER{ 0, 0, 0 };   // black
inline const QColor COLOR_ART{ 150, 40, 40 };       // dark red
inline const QColor COLOR_ART_PULM{ 40, 60, 150 };  // dark blue
inline const QColor COLOR_VCG{ 75, 0, 130 };        // dark indigo
// SHHS respiratory / oximetry traces. AIRFLOW, THOR and ABDO share one chart,
// so these only need to be distinguishable from each other -- they are never
// drawn alongside the ECG/PPG palette above.
inline const QColor COLOR_FLOW{ 0, 100, 130 };      // dark teal-blue
inline const QColor COLOR_THOR{ 120, 70, 0 };       // dark amber
inline const QColor COLOR_ABDO{ 90, 0, 120 };       // dark purple
inline const QColor COLOR_SPO2{ 0, 90, 60 };        // dark green

// ---------------------------------------------------------------------------
// Qt-side lookups into the annotation table
// ---------------------------------------------------------------------------
//
// The table itself is in annotation_types.hpp, which is Qt-free so the Qt-free
// template-generation side can read it directly. These five helpers are the only
// parts that ever needed QString or QColor, and they live here because this is
// already the Qt color-helper header and is already included by every
// translation unit that calls them (gui_handler.cpp, signal_renderer.cpp,
// user_marking_handler.cpp). Reopening the namespace is deliberate: call sites
// keep saying annotation_types::colorFor(...) and nothing had to change.
namespace annotation_types {

    inline const AnnotationType* find(const QString& label) {
        return find(label.toStdString());
    }

    // Highlight color for a type; default translucent black if unknown
    // (matches the old updateNoiseHighlights fallback).
    inline QColor colorFor(const QString& label) {
        if (const auto* t = find(label)) return QColor(t->r, t->g, t->b, t->a);
        return QColor(0, 0, 0, 100);
    }

    inline bool isParamEdit(const QString& label) {
        const auto* t = find(label);
        return t && t->paramEdit;
    }
    inline bool isInvertEdit(const QString& label) {
        const auto* t = find(label);
        return t && t->invertEdit;
    }
    inline bool includeInThreshold(const QString& label) {
        const auto* t = find(label);
        return t && t->includeInThreshold;
    }

}  // namespace annotation_types

inline void wipe_chart(QChart* chart, const QList<QAbstractSeries*>& keep = {}) {
    /*
    * Only ever called in handle_data_plot in signal_renderer. It is to remove things from charts when the signal is changed.
    * The QList 'keep' contains everything that is supposed to be persistent and
    * calculated once: the raw sample, the upsampled sample, and the red grid lines
    * (if they are there). This is for efficiency so they aren't redrawn. However
    * the annotations, etc are removed.
    */
    for (auto* s : chart->series()) {
        if (keep.contains(s)) { chart->removeSeries(s); continue; }
        chart->removeSeries(s);
        delete s;
    }
    for (auto* a : chart->axes()) { chart->removeAxis(a); delete a; }
}


inline void set_padded_y_range(QValueAxis* yAxis, double yMin, double yMax) {
    /*
    * The y range of a plot isn't just the max min diff, there is a pad. It can be tuned here!
    * Handles range in event of flatline.
    */
    if (!(yMin <= yMax)) { yMin = -1.0; yMax = 1.0; }   // NaN/inverted (no valid data): safe default
    double span = yMax - yMin;
    double pad = 0.05 * span;
    if (!span) { pad = 0.05; }//flatline
    yAxis->setRange(yMin - pad, yMax + pad);
}

inline bool is_missing_signal(const QVector<double>& data) {
    /*
    * Returns true if the file doesn't have a the signal. Either the dataset doesn't have this type of data,
    * or the experimenter forgot to attach the lead.
    */
    return data.isEmpty() || (data.size() == 1 && data[0] == -1.0);
}

inline bool sleep_data_present(const QVector<double>& sleepStages) {
    /*
    * Sleep staging comes from the dataset's annotation file: MESA and SHHS1/SHHS2 both ship
    * one, Bittium and CHAOS do not. This is the check for a file whose staging is absent or
    * was written as the missing-channel placeholder.
    */
    return  !sleepStages.isEmpty() && !(sleepStages.size() == 1 && sleepStages[0] == -1.0);
}

inline QString get_timestamp(double seconds) {
    /*
    * The timestamps along the X axis of the charts are formatted as: HH:MM:SS.t, the tenth is
    * because if you don't, the display is weird for 1 second windows
    */
    if (seconds < 0.0) seconds = 0.0;
    int tenths = static_cast<int>(seconds * 10.0 + 0.5);   // round to 0.1 s
    int h = tenths / 36000;
    int m = (tenths % 36000) / 600;
    int s = (tenths % 600) / 10;
    int d = tenths % 10;
    return QString("%1:%2:%3.%4")
        .arg(h, 2, 10, QChar('0'))
        .arg(m, 2, 10, QChar('0'))
        .arg(s, 2, 10, QChar('0'))
        .arg(d);
}

inline QString get_chart_title(const QString& signalName, double nativeHz, double pxPerSample, double bpm = -1.0, double upHz = -1.0) {
    /*On the chart list the native rate, the upsampled rate, the pixel resolution, and the bpm if it is a heart signal.
    The bpm is optional because not all signals have a bpm.*/
    QString space = QString(QChar(0x00A0)).repeated(4);
    QString base = signalName + space
        + QString("Original Frequency: %1 Hz").arg(nativeHz, 0, 'f', 1);
    if (upHz >= 0.0)
        base += space + QString("Upsampled Frequency: %1 Hz").arg(upHz, 0, 'f', 1);
    base += space + QString("Pixel Resolution: %1 px/sample").arg(pxPerSample, 0, 'f', 3);
    if (bpm >= 0.0)
        base += space + QString("%1 bpm").arg(bpm, 0, 'f', 0);
    return base;
}