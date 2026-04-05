/**
 * @file   gui_handler.h
 * @brief  Main dialog for viewing ECG/PPG signals and annotating noise,
 *         arrhythmias, and artifacts. Loads data in 8-hour chunks from
 *         a binary file produced by file_to_bin.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 */
#pragma once

#include <QtWidgets/QDialog>
#include <QtWidgets/QRubberBand>
#include <QtWidgets/QPushButton>
#include <QtCharts/QLineSeries>
#include <QtCharts/QAreaSeries>
#include <QtCharts/QChartView>
#include <QVector>
#include <QSet>
#include <QColor>
#include <memory>

#include "ui_noise_marking_gui.h"
#include "NoiseManager.h"
#include "lower_row_buttons.h"

 // Collected annotations returned to the caller after dialog closes
struct GenExcStruct {
    QVector<QPair<double, double>> noiseExc;
    QStringList data_type;     // "ECG1", "ECG2", "ECG3", or "PPG"
    QStringList marking_type;  // "Noise/Artifact", "AF", "SVT", etc.
};

class noise_marking_gui : public QDialog {
    Q_OBJECT
        friend class lower_row_buttons;

public:
    explicit noise_marking_gui(QWidget* parent = nullptr);
    ~noise_marking_gui() override;

    GenExcStruct getMarkings() const { return m_genExc; }
    void         setFileSource(const QString& filePath);

protected:
    bool eventFilter(QObject* watched, QEvent* event) override;

private slots:
    void on_skip_interval_box_returnPressed();
    void on_skip_interval_box_editingFinished();
    void on_marking_type_currentTextChanged(const QString& text);
    void on_next8hours_clicked();
    void on_prev8hours_clicked();

private:
    // --- Per-channel marking state machine ---
    //
    //  Idle
    //    ──(Start button)──▸ WaitingForStart
    //    ──(mouse press)───▸ Dragging         (handled outside this struct)
    //
    //  WaitingForStart
    //    ──(click)─────────▸ WaitingForEnd    (start marker placed)
    //
    //  WaitingForEnd
    //    ──(click)─────────▸ WaitingForEnd    (start marker moved)
    //    ──(Stop button)───▸ WaitingForStop
    //    ──(Start button)──▸ Idle             (cancelled)
    //
    //  WaitingForStop
    //    ──(click)─────────▸ Idle             (finalizeMarking)
    //
public:
    enum class MarkPhase { Idle, WaitingForStart, WaitingForEnd, WaitingForStop };

    struct ChannelMarkingState {
        MarkPhase    phase = MarkPhase::Idle;
        double       globalStartTime = 0.0;   // always in global seconds
        QLineSeries* startMarkerLine = nullptr;
    };

    // Called by lower_row_buttons to enter/transition marking phases
    void beginMarking(const QString& signalLabel);
    void beginStopPhase(const QString& signalLabel);
    void beginMarkingAll();
    void beginStopPhaseAll();

private:

    ChannelMarkingState  m_markState_ecg1;
    ChannelMarkingState  m_markState_ecg2;
    ChannelMarkingState  m_markState_ecg3;
    ChannelMarkingState  m_markState_ppg;

    // --- Core components ---
    std::unique_ptr<Ui::noise_marking_gui> ui;
    std::unique_ptr<NoiseManager>          m_noiseManager;
    std::unique_ptr<lower_row_buttons>     m_buttonHandler;

    QSet<QString> m_activeChannels;
    QString       m_currentMarkingType;
    GenExcStruct  m_genExc;

    // --- Sampling rates ---
    double m_ecgSR = 0.0;
    double m_ppgSR = 0.0;
    double m_sleepSR = 0.0;

    // --- View state ---
    double m_currentStartTime = 0.0;
    double m_windowDuration = 10.0;
    double m_skipInterval = 5.0;

    // --- Ampogram / hypnogram series ---
    QLineSeries* m_ecgAmpSeries = nullptr;
    QLineSeries* m_ppgAmpSeries = nullptr;
    QLineSeries* m_ecgCursorBar = nullptr;
    QLineSeries* m_ppgCursorBar = nullptr;
    QLineSeries* m_hypnoCursorBar = nullptr;
    QList<QAreaSeries*>     m_highlights;
    QList<QAbstractSeries*> m_hypnoStageSeries;

    // --- Mark-all state ---
    bool m_markAllActive = false;   // true while "Mark All Signals" mode is in effect

    // --- Drag state ---
    bool     m_isDragging = false;
    QWidget* m_draggedViewport = nullptr;
    QPoint   m_dragStartPos;
    QString  m_dragSignalLabel;

    // --- Signal data ---
    QVector<double> m_ecg1, m_ecg2, m_ecg3, m_ppg, m_sleepStages;

    // --- File layout ---
    QString  m_binFilePath;
    qint64   m_fileHeaderSize = 0;
    uint64_t m_totalEcg1Samples = 0;
    uint64_t m_totalEcg2Samples = 0;
    uint64_t m_totalEcg3Samples = 0;
    uint64_t m_totalPpgSamples = 0;
    uint64_t m_totalSleepSamples = 0;
    uint64_t m_currentChunkIndex = 0;

    static constexpr double CHUNK_DURATION_SEC = 28800.0;  // 8 hours

    // --- Lookup helpers ---
    ChannelMarkingState& markStateFor(const QString& label);
    QString              signalLabelForChartView(QChartView* cv) const;
    QChartView* chartViewForSignalLabel(const QString& label) const;
    double               sampleRateForSignal(const QString& label) const;
    QColor               colorForSignal(const QString& label) const;
    bool                 isChannelActive(const QString& label) const;
    double               totalChunkDuration() const;

    // --- Per-channel button helpers ---
    QPushButton* startButtonForSignal(const QString& label) const;
    QPushButton* stopButtonForSignal(const QString& label) const;
    void         updateButtonStatesForChannel(const QString& label);
    void         updateAllChannelButtonStates();

    // --- Plotting ---
    bool loadChunkFromFile(uint64_t chunkIndex);
    void handle_data_plot();
    void handle_ampogram_plot(double sampling_length = 60);
    void updateAmpogramCursor();
    void setupHypnogram();
    void updateNoiseHighlights();

    // --- Marking ---
    void finalizeMarking(QChartView* cv, double endX, const QString& signalLabel);
    void cancelMarking(const QString& signalLabel);
    void restoreMarkingMarkers();
    void showStartMarker(QChartView* cv, double xValue, ChannelMarkingState& state,
        const QColor& color, QPushButton* stopBtn);
    void clearStartMarker(ChannelMarkingState& state);

    // --- Formatting ---
    QString formatTimeLabel(double seconds);
};