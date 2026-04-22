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
#include <QMap>
#include <QColor>
#include <memory>

#include "ui_noise_marking_gui.h"
#include "NoiseManager.h"
#include "lower_row_buttons.h"

 // Collected annotations returned to the caller after dialog closes
struct GenExcStruct {
    QString filePath;          // source file these markings belong to
    QVector<QPair<double, double>> noiseExc;
    QStringList data_type;     // "ECG1", "ECG2", "ECG3", "PPG", or "ABP"
    QStringList marking_type;  // "Noise/Artifact", "AF", "SVT", etc.
};

class noise_marking_gui : public QDialog {
    Q_OBJECT
        friend class lower_row_buttons;

public:
    explicit noise_marking_gui(QWidget* parent = nullptr);
    ~noise_marking_gui() override;

    GenExcStruct              getMarkings() const;
    QVector<GenExcStruct>     getAllMarkings() const;
    QString                   getFilePath() const { return m_binFilePath; }
    void                      setFileSource(const QString& filePath);

protected:
    bool eventFilter(QObject* watched, QEvent* event) override;

private slots:
    void on_skip_interval_box_returnPressed();
    void on_skip_interval_box_editingFinished();
    void on_marking_type_currentTextChanged(const QString& text);
    void on_next8hours_clicked();
    void on_prev8hours_clicked();
    void handleBrowseFile();

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
    enum class PlotMode { Line, Scatter };

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
    ChannelMarkingState  m_markState_abp;

    // --- Core components ---
    std::unique_ptr<Ui::noise_marking_gui> ui;
    std::unique_ptr<NoiseManager>          m_noiseManager;
    std::unique_ptr<lower_row_buttons>     m_buttonHandler;

    QSet<QString> m_activeChannels;
    QString       m_currentMarkingType;
    GenExcStruct  m_genExc;
    QMap<QString, GenExcStruct> m_fileMarkings;  // stashed markings per file path

    // --- Sampling rates ---
    double m_ecgSR = 0.0;
    double m_ppgSR = 0.0;
    double m_boolSR = 0.0;
    double m_sleepSR = 0.0;

    // --- View state ---
    double m_currentStartTime = 0.0;
    double m_windowDuration = 10.0;
    double m_skipInterval = 5.0;

    // --- Ampogram / hypnogram series ---
    QLineSeries* m_ecgAmpSeries = nullptr;
    QLineSeries* m_ppgAmpSeries = nullptr;
    QLineSeries* m_respAmpSeries = nullptr;   // RESP 8h overview (min/max-decimated raw)
    QLineSeries* m_cvpAmpSeries = nullptr;   // CVP  8h overview (min/max-decimated raw)
    QLineSeries* m_ecgCursorBar = nullptr;
    QLineSeries* m_ppgCursorBar = nullptr;
    QLineSeries* m_respCursorBar = nullptr;
    QLineSeries* m_cvpCursorBar = nullptr;
    QLineSeries* m_hypnoCursorBar = nullptr;
    QList<QAreaSeries*>     m_highlights;
    QList<QAbstractSeries*> m_hypnoStageSeries;

    // --- Mark-all state ---
    bool m_markAllActive = false;   // true while "Mark All Signals" mode is in effect

    // --- Plot style (global, applies to all signal charts) ---
    PlotMode m_plotMode = PlotMode::Line;

    // --- Drag state ---
    bool     m_isDragging = false;
    QWidget* m_draggedViewport = nullptr;
    QPoint   m_dragStartPos;
    QString  m_dragSignalLabel;

    // --- Signal data ---
    QVector<double> m_ecg1, m_ecg2, m_ecg3, m_ppg;
    QVector<double> m_accelX, m_accelY, m_accelZ;
    QVector<double> m_resp;        // respiration (slot 34)
    QVector<double> m_cvp;         // central venous pressure (slot 16, Bittium)
    QVector<double> m_abp;         // arterial blood pressure (slot 35, Bittium)
    QVector<double> m_sleepStages;

    // --- File layout (180-byte header: 45 x uint32 = 4 rates + 40 chan-sizes + 1 sleep) ---
    QString  m_binFilePath;
    static constexpr qint64 FILE_HEADER_SIZE = 180;

    uint32_t m_chanSizes[40] = {};
    uint32_t m_totalSleepSamples = 0;

    static constexpr int CH_ECG1 = 0, CH_ECG2 = 1, CH_ECG3 = 2, CH_PPG = 3;
    static constexpr int CH_ACCEL_X = 4, CH_ACCEL_Y = 5, CH_ACCEL_Z = 6;
    static constexpr int CH_MARKER = 7, CH_TEMP = 8, CH_PACEMAKER = 9;
    static constexpr int CH_EOG_L = 10, CH_EOG_R = 11, CH_EMG = 12;
    static constexpr int CH_EEG1 = 13, CH_EEG2 = 14, CH_EEG3 = 15;
    static constexpr int CH_PRES = 16, CH_FLOW = 17, CH_THOR = 18;
    static constexpr int CH_ABDO = 19, CH_LEG = 20, CH_THERM = 21;
    static constexpr int CH_POS = 22;
    static constexpr int CH_EKG_OFF = 23, CH_EOGL_OFF = 24, CH_EOGR_OFF = 25;
    static constexpr int CH_EMG_OFF = 26, CH_EEG1_OFF = 27, CH_EEG2_OFF = 28;
    static constexpr int CH_EEG3_OFF = 29;
    static constexpr int CH_OXSTATUS = 30, CH_SPO2 = 31;
    static constexpr int CH_HR = 32, CH_DHR = 33, CH_RESP = 34;
    static constexpr int CH_ABP = 35, CH_EEG4 = 36;
    static constexpr int CH_ART = 37, CH_ART_PULM = 38;

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

    // --- File selection ---
    void loadSelectedFile(const QString& filePath);

    // --- Formatting ---
    QString formatTimeLabel(double seconds);
};