#pragma once
/**
* @file   chart_utils.hpp
* @brief  Sets the colors of all the series and contains various
*         features required by the plotting functionality - both markable and
*         nonmarkable
*/

#include <QtCharts/QChart>
#include <QtCharts/QValueAxis>
#include <QString>
#include <QVector>

inline const QColor COLOR_ECG1{ 101, 67, 33 };
inline const QColor COLOR_ECG2{ 0,  128, 0 };
inline const QColor COLOR_ECG3{ 0, 0, 139 };
inline const QColor COLOR_PPG{ 48, 25, 52 };
inline const QColor COLOR_ABP{ 52, 25, 48 };
inline const QColor COLOR_ACCEL_X{ 243, 156, 18 };
inline const QColor COLOR_ACCEL_Y{ 39, 174, 96 };
inline const QColor COLOR_ACCEL_Z{ 142, 68, 173 };
inline const QColor COLOR_RESP{ 0,130,0};
inline const QColor COLOR_CVP{ 130,0,0};
inline const QColor COLOR_RAW_SCATTER{ 0, 0, 0};


inline void wipe_chart(QChart* chart, const QList<QAbstractSeries*>& keep = {}){
    /*
        The QList 'keep' contains everything that is supposed to be persistant and
        calculated once: the raw sample, the upsampled sample, and the red grid lines
        (if they are there). This is for efficiency so they aren't redrawn. However
        the annotations, etc are removed.
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
        The y range isn't just the max min diff, there is a pad. It can be tuned here!
        Handles range in event of flatline.
    */
    double span = yMax - yMin;
    double pad = 0.05 * span;
    if (!span) {pad = 0.05;}//flatline
    yAxis->setRange(yMin - pad, yMax + pad);
}

inline bool isMissingSignal(const QVector<double>& data) {
    /*
        You can't use dataset type to determine if signal is empty because sometimes the
        experimentor will have forgotten to put a lead on, or something.
    */
    return data.isEmpty() || (data.size() == 1 && data[0] == -1.0);
}

inline bool sleepDataPresent(const QVector<double>& sleepStages) {
    /*
        Currently, sleep data is only in the MESA, and in all MESA files. This is just an extra check
        in case there is a broken MESA file.
    */
    return  !sleepStages.isEmpty()
            && !(sleepStages.size() == 1 
            && sleepStages[0] == -1.0);
}

inline QString get_timestamp(double seconds) {
    /*
        The timestamps along the X axis of the charts are formatted as: HH:MM:SS
        They are set by this function
    */
    int t = static_cast<int>(seconds + 0.5);
    int h = t / 3600;
    int m = (t % 3600) / 60;
    int s = t % 60;
    return QString("%1:%2:%3")
        .arg(h, 2, 10, QChar('0'))
        .arg(m, 2, 10, QChar('0'))
        .arg(s, 2, 10, QChar('0'));
}

inline QString get_chart_title(const QString& signalName, double nativeHz, double pxPerSample, double bpm){
    //Make title for each markable chart, set sig figs in each printed value
    QString space = QString(QChar(0x00A0)).repeated(4); //QString doesn't have tabs and doesn't respsect whitespace for some reason
    QString base = QString("%1" + space + "Original Frequency: %2 Hz" + space + "Pixel Resolution: %3 px/sample" + space + "%4 bpm")
        .arg(signalName)
        .arg(nativeHz, 0, 'f', 1)
        .arg(pxPerSample, 0, 'f', 3)
        .arg(bpm, 0, 'f', 0);
    return base;
}