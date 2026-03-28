/**
 * @file   noise_marking_gui.h
 * @brief  Main dialog for viewing ECG/PPG signals and annotating noise,
 *         arrhythmias, and artifacts. Loads data in 8-hour chunks from
 *         a binary file produced by file_to_bin.
 */
#pragma once

#include <QtWidgets/QDialog>
#include <QVector>
#include <QtWidgets/QRubberBand>
#include <QPair>
#include "ui_noise_marking_gui.h"
#include "NoiseManager.h"
#include "lower_row_buttons.h"

struct GenExcStruct {
    QVector<QPair<double, double>> noiseExc;
    QStringList data_type;     // "ECG1", "ECG2", "ECG3", or "PPG"
    QStringList marking_type;  // "Noise/Artifact", "AF", "SVT", etc.
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
    ~noise_marking_gui() override;

    GenExcStruct getMarkings() const { return m_genExc; }
    void setFileSource(const QString& filePath);

protected:
    bool eventFilter(QObject* watched, QEvent* event) override;
    void finalizeMarking(QChartView* cv, double endX, const QString& signalLabel);

private slots:
    void on_skip_interval_box_returnPressed();
    void on_skip_interval_box_editingFinished();
    void on_marking_type_currentTextChanged(const QString& text);
    void start_marking_button_clicked(const QString& signalLabel);
    void on_next8hours_clicked();
    void on_prev8hours_clicked();

private:
    // --- Core components ---
    std::unique_ptr<Ui::noise_marking_gui> ui;
    std::unique_ptr<NoiseManager> m_noiseManager;
    std::unique_ptr<QRubberBand> m_rubberBand;
    std::unique_ptr<lower_row_buttons> m_buttonHandler;

    // --- Marking state (per signal channel) ---
    struct ChannelMarkingState {
        bool isWaitingForStart = false;
        bool isWaitingForEnd = false;
        double startTimeValue = 0.0;
        QLineSeries* startMarkerLine = nullptr;
    };

    ChannelMarkingState m_markState_ecg1;
    ChannelMarkingState m_markState_ecg2;
    ChannelMarkingState m_markState_ecg3;
    ChannelMarkingState m_markState_ppg;

    ChannelMarkingState& markStateFor(const QString& label);

    QString m_currentMarkingType;
    bool m_isDragging = false;
    GenExcStruct m_genExc;

    // --- Sampling rates ---
    double m_ecgSR = 0.0;
    double m_ppgSR = 0.0;
    double m_sleepSR = 0.0;

    // --- View state ---
    double m_currentStartTime = 0.0;
    double m_windowDuration = 10.0;
    double m_skipInterval = 5.0;

    // --- Chart series ---
    QLineSeries* m_ecgAmpSeries = nullptr;
    QLineSeries* m_ppgAmpSeries = nullptr;
    QLineSeries* m_ecgCursorBar = nullptr;
    QLineSeries* m_ppgCursorBar = nullptr;
    QLineSeries* m_hypnoCursorBar = nullptr;
    QList<QAreaSeries*> m_highlights;
    QList<QAbstractSeries*> m_hypnoStageSeries;

    // --- Drag state ---
    QWidget* m_draggedViewport = nullptr;
    QPoint m_dragStartPos;
    QString m_dragSignalLabel;

    // --- Signal data ---
    QVector<double> m_ecg1, m_ecg2, m_ecg3, m_ppg, m_sleepStages;

    // --- File layout (88-byte header) ---
    QString m_binFilePath;
    qint64 m_fileHeaderSize = 0;
    uint64_t m_totalEcg1Samples = 0;
    uint64_t m_totalEcg2Samples = 0;
    uint64_t m_totalEcg3Samples = 0;
    uint64_t m_totalPpgSamples = 0;
    uint64_t m_totalSleepSamples = 0;
    uint64_t m_totalAbs1Samples = 0;  // read but ignored
    uint64_t m_totalAbs2Samples = 0;  // read but ignored
    uint64_t m_totalAbs3Samples = 0;  // read but ignored
    uint64_t m_currentChunkIndex = 0;
    static constexpr double CHUNK_DURATION_SEC = 28800.0; // 8 hours

    // --- Private methods ---
    bool loadChunkFromFile(uint64_t chunkIndex);
    void handle_data_plot();
    void handle_ampogram_plot(double sampling_length = 60);
    void updateAmpogramCursor();
    void setupHypnogram();
    void updateNoiseHighlights();
    void showStartMarker(QChartView* cv, double xValue, ChannelMarkingState& state, const QColor& color, QPushButton* stopBtn);
    void clearStartMarker(ChannelMarkingState& state);
    QString formatTimeLabel(double seconds);

    // Helpers to map between signal labels and chart views / properties
    QString selectedEcgLabel() const;                        // reads ecg_channel_selector
    QString signalLabelForChartView(QChartView* cv) const;
    QChartView* chartViewForSignalLabel(const QString& label) const;
    double sampleRateForSignal(const QString& label) const;
    QColor colorForSignal(const QString& label) const;
};