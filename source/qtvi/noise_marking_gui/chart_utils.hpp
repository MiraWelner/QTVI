#pragma once
/**
 * @file   chart_utils.hpp
 * @brief  Internal helpers shared across gui_handler_*.cpp translation units.
 *         Not part of the public API -- do not include from headers.
 */

#include "annotation_types.hpp"

#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QValueAxis>
#include <QtCharts/QCategoryAxis>
#include <QFont>
#include <QString>
#include <QVector>
#include <cmath>
#include <QtWidgets/QGraphicsLayout>
#include <QMap>

inline const QColor COLOR_ECG1{ 101, 67, 33 };
inline const QColor COLOR_ECG2{ 0,  128, 0 };
inline const QColor COLOR_ECG3{ 0, 0, 139 };
inline const QColor COLOR_PPG{ 48, 25, 52 };
inline const QColor COLOR_ABP{ 52, 25, 48 };
inline const QColor COLOR_ACCEL_X{ "#F39C12" };
inline const QColor COLOR_ACCEL_Y{ "#27AE60" };
inline const QColor COLOR_ACCEL_Z{ "#8E44AD" };
inline const QColor COLOR_RESP{ 0,0,0 };
inline const QColor COLOR_CVP{ 0,0,0 };
inline const QColor COLOR_RAW_SCATTER{ 0, 0, 0, 255 };

// ---------------------------------------------------------------------------
// Chart helpers
// ---------------------------------------------------------------------------

inline void wipeChartContent(QChart* chart,
    const QList<QAbstractSeries*>& keep = {})
{
    if (!chart) return;
    for (auto* s : chart->series()) {
        if (keep.contains(s)) { chart->removeSeries(s); continue; }
        chart->removeSeries(s);
        delete s;
    }
    for (auto* a : chart->axes()) { chart->removeAxis(a); delete a; }
}


inline void setPaddedYRange(QValueAxis* yAxis, double yMin, double yMax) {
    if (yMin > yMax) { yMin = -1.0; yMax = 1.0; }
    double span = yMax - yMin;
    double pad = (span > 1e-9) ? 0.05 * span : 0.5;
    yAxis->setRange(yMin - pad, yMax + pad);
}

// ---------------------------------------------------------------------------
// Signal helpers
// ---------------------------------------------------------------------------

inline bool isMissingSignal(const QVector<double>& data) {
    return data.isEmpty() || (data.size() == 1 && data[0] == -1.0);
}

inline bool isRawUsable(const QVector<QPointF>& v) {
    return v.size() >= 2 && !(v.size() == 1 && v[0].x() == -1.0);
}

inline bool sleepDataPresent(const QVector<double>& sleepStages) {
    return !sleepStages.isEmpty()
        && !(sleepStages.size() == 1 && sleepStages[0] == -1.0);
}

// ---------------------------------------------------------------------------
// Formatting
// ---------------------------------------------------------------------------

inline QString formatHMS(double seconds) {
    int t = static_cast<int>(seconds + 0.5);
    int h = t / 3600;
    int m = (t % 3600) / 60;
    int s = t % 60;
    return QString("%1:%2:%3")
        .arg(h, 2, 10, QChar('0'))
        .arg(m, 2, 10, QChar('0'))
        .arg(s, 2, 10, QChar('0'));
}

inline QString formatChartTitle(const QString& signalName,
    double nativeHz, double pxPerSample,
    double bpm = -1.0)
{
    QString base = QString("%1  -- Original Frequency: %2 Hz -- Pixel Resolution: %3 px/sample")
        .arg(signalName)
        .arg(nativeHz, 0, 'f', 1)
        .arg(pxPerSample, 0, 'f', 3);
    if (bpm >= 0.0)
        base += QString("  --  %1 bpm").arg(bpm, 0, 'f', 0);
    return base;
}