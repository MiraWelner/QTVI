#pragma once

#include <QtWidgets/QDialog>
#include <QVector>
#include <QtWidgets/QRubberBand>
#include <QPair>
#include "ui_noise_marking_gui.h"
#include "NoiseManager.hpp"
#include "lower_row_buttons.hpp"

struct GenExcStruct {
    QVector<QPair<double, double>> noiseExc;
    QStringList data_type; //ECG or PPG
    QStringList marking_type;  // Noise/Artifact, AF, SVT, VT, PVC, PAC, Benign Arrhythmia, Significant Arrhythmia

};

class QLineSeries;
class QAbstractSeries;
class QAreaSeries;
class QChartView;

class noise_marking_gui : public QDialog {
    Q_OBJECT

        friend class lower_row_buttons;

public:
    explicit noise_marking_gui(QWidget* parent = nullptr);
    QString m_currentMarkingType;
    ~noise_marking_gui() override;

    GenExcStruct getMarkings() const { return m_genExc; }
    void setFileSource(const QString& filePath);

protected:
    bool eventFilter(QObject* watched, QEvent* event) override;
    void finalizeMarking(QChartView* cv, double endX, bool isECG);

private slots:
    void on_skip_interval_box_returnPressed();
    void on_skip_interval_box_editingFinished();
    void on_marking_type_currentTextChanged(const QString& text);
    void start_marking_button_clicked(bool isECG);
    void on_next8hours_clicked();
    void on_prev8hours_clicked();

private:
    std::unique_ptr<Ui::noise_marking_gui> ui;
    std::unique_ptr<NoiseManager> m_noiseManager;
    std::unique_ptr<QRubberBand> m_rubberBand;
    std::unique_ptr<lower_row_buttons> m_buttonHandler;
    QWidget* m_draggedViewport = nullptr;


    bool m_isWaitingForECGStart = false;
    bool m_isWaitingForPPGStart = false;
    bool m_isWaitingForECGEnd = false;
    bool m_isWaitingForPPGEnd = false;
    bool m_isDragging = false;


    double m_ecgStartTimeValue = 0.0;
    double m_ppgStartTimeValue = 0.0;
    double m_ppgSR = 0.0, m_ecgSR = 0.0, m_sleepSR = 0.0;
    double m_currentStartTime = 0.0;
    double m_windowDuration = 10.0;
    double m_skipInterval = 5.0;

    GenExcStruct m_genExc;
    QLineSeries* m_ecgStartMarkerLine = nullptr;
    QLineSeries* m_ppgStartMarkerLine = nullptr;
    QLineSeries* m_ecgAmpSeries = nullptr;
    QLineSeries* m_ppgAmpSeries = nullptr;
    QLineSeries* m_ecgCursorBar = nullptr;
    QLineSeries* m_ppgCursorBar = nullptr;
    QLineSeries* m_hypnoCursorBar = nullptr;
    QList<QAreaSeries*> m_highlights;
    QPoint m_dragStartPos;
    QVector<double> m_ppg, m_ecg, m_ecg2, m_ecg3, m_sleepStages;
    QList<QAbstractSeries*> m_hypnoStageSeries;

    void handle_data_plot();
    void showStartMarker(QChartView* cv, double xValue, bool isECG);
    void handle_ampogram_plot(double sampling_length = 60);
    void updateAmpogramCursor();
    void setupHypnogram();
    void updateNoiseHighlights();
    void clearECGStartMarker();
    void clearPPGStartMarker();
    bool loadChunkFromFile(uint64_t chunkIndex);
    QString m_binFilePath;
    QString formatTimeLabel(double seconds);

    uint64_t m_totalEcgSamples = 0, m_totalPpgSamples = 0, m_totalSleepSamples = 0, m_totalSignal2Samples = 0, m_totalSignal3Samples = 0;

    uint64_t m_currentChunkIndex = 0;
    qint64 m_fileHeaderSize = 0;
    const double CHUNK_DURATION_SEC = 28800.0; // 8 hours (8 * 3600)

};