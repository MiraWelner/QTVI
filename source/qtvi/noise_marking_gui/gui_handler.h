/**
 * @file    gui_handler.hpp
 * @brief   Main dialog for viewing ECG/PPG signals and annotating noise,
 *          arrhythmias, and artifacts. Loads data in 8-hour chunks from a
 *          binary file produced by file_to_bin.
 *
 * @author  Mira Welner
 * @email   MEW386@pitt.edu
 */
#pragma once

#include <QtWidgets/QDialog>
#include <QtWidgets/QPushButton>
#include <QtCharts/QLineSeries>
#include <QtCharts/QAreaSeries>
#include <QtCharts/QChartView>
#include <QtCharts/QScatterSeries>
#include <QVector>
#include <QSet>
#include <QMap>
#include <QHash>
#include <QColor>
#include <QPointF>
#include <QTimer>
#include <memory>
#include <vector>

#include "ui_noise_marking_gui.h"
#include "user_annotation_handler.h"
#include "user_control_handler.h"
#include "config_entry.hpp"
#include "grid_overlay.hpp"
#include "gap_indicator.hpp"

struct GenExcStruct {
    QString filePath;                          ///< Source file these markings belong to.
    QVector<QPair<double, double>> noiseExc;   ///< [start, end] in global seconds.
    QStringList data_type;                     ///< "ECG1", "ECG2", "ECG3", "PPG", or "ABP".
    QStringList marking_type;                  ///< "Noise/Artifact", "AF", "SVT", ...
};

struct markable_data_series {
    const QVector<double>* data;
    QColor color;
    const QVector<QPointF>* rawData;
};

class beat_log;   // defined in beat_log.hpp; only a pointer is held here
class annotation_eraser;   // defined in annotation_eraser.h; pointer held below

class noise_marking_gui : public QDialog {
    Q_OBJECT
    friend class user_control_handler;
    friend class gap_indicator;
    friend class annotation_eraser;

public:
    /**
     * @param parent Optional Qt parent widget.
     */
    explicit noise_marking_gui(QWidget* parent = nullptr);
    ~noise_marking_gui() override;
    GenExcStruct getMarkings() const;
    QVector<GenExcStruct> getAllMarkings() const;
    QString getFilePath() const { return m_binFilePath; }
    void setFileSource(const QString& filePath);
    double currentStartTime() const { return current_start_time; }
    double windowDuration()   const { return visible_window_size; }
    QChartView* chartViewForSignalLabel(const QString& label) const;
    bool        isChannelActive(const QString& label) const;
    enum class MarkPhase { Idle, WaitingForStart, WaitingForEnd, WaitingForStop };
    enum class MarkScope { One, Ecg, All };
    enum class PlotMode { Line, Scatter };
    void setBeatLog(beat_log* log) { m_beatLog = log; }
    void set_params_to_config_defaults(const config_entry& cfg) {
        // Set the default values for the threshold and blanking period spinboxes based on the config entry.
        m_cfg = cfg;
    }

protected:
    //these override native QT event handlers which is why they are protected
    bool eventFilter(QObject* watched, QEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;

private slots:
    void on_skip_interval_box_editingFinished();
    void on_marking_type_currentTextChanged(const QString& text);
    void on_next8hours_clicked();
    void on_prev8hours_clicked();
    void handleBrowseFile();

private:
    double scrollStepSeconds() const;
    config_entry m_cfg;
    beat_log* m_beatLog = nullptr;
    QTimer* m_logFlushTimer = nullptr;   // flushes the beat log to disk every 30 s

    struct ChannelMarkingState {
        MarkPhase    phase = MarkPhase::Idle;
        double       globalStartTime = 0.0;   ///< Always in global seconds (not chunk-local).
        QLineSeries* startMarkerLine = nullptr;
    };

    struct data_channel_features {
        QChartView* chartView = nullptr;
        QPushButton* startButton = nullptr;
        QPushButton* stopButton = nullptr;
        ChannelMarkingState* state = nullptr;
        const QVector<double>* upsampled_data = nullptr;
        const QVector<QPointF>* dataRaw = nullptr;   ///< raw (t, v) pairs, chunk-local seconds
        const double* sampleRate = nullptr;  //the original sampling rate of the raw data
        QColor  color;
    };
    struct ParamOverride {
        QString channel;
        double  start = 0.0;
        double  end = 0.0;
        double  value = 0.0;
    };

    void setMarkScope(MarkScope scope);
    void toggleAnnotationArm();
    void updateMarkingButtons();
    QStringList scopeChannels(const QString& clickedLabel) const;

    data_channel_features channelRefs(const QString& label) const;
    static const QStringList& markableChannelLabels();

    // --- Per-channel marking state ---
    ChannelMarkingState  m_markState_ecg1;
    ChannelMarkingState  m_markState_ecg2;
    ChannelMarkingState  m_markState_ecg3;
    ChannelMarkingState  m_markState_ppg;
    ChannelMarkingState  m_markState_abp;

    std::unique_ptr<Ui::noise_marking_gui> ui;
    std::unique_ptr<annotation_handler>       m_noiseManager;
    std::unique_ptr<user_control_handler>     m_buttonHandler;
    std::unique_ptr<pulse_overlay>            m_pulseOverlay;
    std::unique_ptr<gap_indicator>            m_gapIndicator;
    std::unique_ptr<annotation_eraser>        m_annotationEraser;

    QSet<QString>               m_activeChannels;
    QString                     m_currentMarkingType;
    GenExcStruct                m_genExc;
    QMap<QString, GenExcStruct> m_fileMarkings;  ///< stashed markings per file path

    // --- Sampling rates ---
    double m_ecgSR = 0.0;     ///< common upsampled grid rate (raw32[0], e.g. 1000 Hz)
    double m_boolSR = 0.0;    ///< low/bool rate for <=1 Hz channels (raw32[1])
    double m_sleepSR = 0.0;   ///< sleep epochs per second (= 1 / sleep_epoch_length)

    // Channel rate helpers -- the single place that maps a channel to its rate.
    static constexpr double kBoolRateThreshold = 1.0;
    double nativeRateFor(int ch) const { return channel_native_rates[ch]; }
    double upsampledRateFor(int ch) const {
        return (channel_native_rates[ch] <= kBoolRateThreshold) ? m_boolSR : m_ecgSR;
    }
    // --- View state ---
    double current_start_time;
    double visible_window_size;
    double percent_of_window_to_shift_on_click;

    // --- Ampogram / hypnogram series ---
    QLineSeries* ecg1_ampogram_series = nullptr;
    QLineSeries* ecg2_ampogram_series = nullptr;
    QLineSeries* ecg3_ampogram_series = nullptr;
    QLineSeries* ppg_ampogram_series = nullptr;
    QLineSeries* m_ecgCursorBar = nullptr;
    QLineSeries* m_ppgCursorBar = nullptr;
    QLineSeries* m_respCursorBar = nullptr;
    QLineSeries* m_cvpCursorBar = nullptr;
    QLineSeries* m_hypnoCursorBar = nullptr;
    QList<QAreaSeries*>     m_highlights;
    QList<QAbstractSeries*> m_hypnoStageSeries;
    QHash<QChartView*, QList<QLineSeries*>> m_persistentLines;
    QHash<QChartView*, QList<QScatterSeries*>> m_persistentRawScatter;
    MarkScope m_markScope = MarkScope::One;
    bool      m_markArmed = false;;

    // --- Plot style (global, applies to all signal charts) ---
    PlotMode m_plotMode = PlotMode::Line;

    // checked = per-window autoscale (drift hidden)
    bool m_filterBaselineDrift = false;

    // --- Drag state ---
    bool     m_isDragging = false;
    QWidget* m_draggedViewport = nullptr;
    QPoint   m_dragStartPos;
    QString  m_dragSignalLabel;
    QMap<QChartView*, QAreaSeries*> m_dragPreviews;

    // --- Signal data (upsampled, 1 kHz unless otherwise noted) ---
    QVector<double> m_ecg1, m_ecg2, m_ecg3, m_ppg, m_accelX, m_accelY, m_accelZ, m_resp;
    QVector<double> m_cvp, m_abp, m_temp, m_marker,  m_pacemaker,  m_sleepStages;

    QVector<QPointF> m_ecg1Raw, m_ecg2Raw, m_ecg3Raw, m_ppgRaw, m_abpRaw;
    QVector<QPointF> m_accelXRaw, m_accelYRaw, m_accelZRaw, m_respRaw, m_cvpRaw;
    QVector<QPointF> m_tempRaw, m_markerRaw,m_pacemakerRaw;
    QString m_binFilePath;
    static constexpr qint64 FILE_HEADER_SIZE = 496;
    static constexpr int NUM_CHANNELS = 40;

    uint32_t upsampled_channel_sizes[NUM_CHANNELS] = {};
    uint32_t raw_channel_sizes[NUM_CHANNELS] = {};
    float    channel_native_rates[NUM_CHANNELS] = {};
    uint32_t total_sleep_samples = 0;

    // Channel indices. Slot 0 is the new Timestamp channel; every other
    // channel has shifted +1 from the v1 layout. CH_RESERVED_40 is a
    // reserved tail slot (always absent on disk).
    enum ChannelIdx {
        CH_TIMESTAMP = 0,
        CH_ECG1, CH_ECG2, CH_ECG3, CH_PPG,
        CH_ACCEL_X, CH_ACCEL_Y, CH_ACCEL_Z,
        CH_MARKER, CH_TEMP, CH_PACEMAKER,
        CH_EOG_L, CH_EOG_R, CH_EMG,
        CH_EEG1, CH_EEG2, CH_EEG3, CH_EEG4,
        CH_CVP, CH_FLOW, CH_THOR, CH_ABDO,
        CH_LEG, CH_THERM, CH_POS,
        CH_EKG_OFF, CH_EOG_L_OFF, CH_EOG_R_OFF, CH_EMG_OFF,
        CH_EEG1_OFF, CH_EEG2_OFF, CH_EEG3_OFF,
        CH_OXSTATUS, CH_SPO2, CH_HR, CH_DHR,
        CH_RESP, CH_ABP,
        CH_ART, CH_ART_PULM
    };

    uint64_t current_chunk_index = 0;

    static constexpr int seconds_in_memory_at_once = 28800; //8 hours 


    // --- Lookup helpers (thin wrappers over channelRefs) ---
    ChannelMarkingState& markStateFor(const QString& label);
    QString              signalLabelForChartView(QChartView* cv) const;
    double               sampleRateForSignal(const QString& label) const;
    QColor               colorForSignal(const QString& label) const;

    /** @return Duration in seconds of the currently-loaded 8-hour chunk. */
    double               totalChunkDuration() const;

    // --- Per-channel button helpers ---
    QPushButton* startButtonForSignal(const QString& label) const;
    QPushButton* stopButtonForSignal(const QString& label) const;
    bool loadChunkFromFile(uint64_t chunkIndex);
    void handle_data_plot();
	void determine_which_nonmarkable_charts_to_plot();
    void plot_nonmarkable(QChartView* view, const QString& title, const QList<markable_data_series>& serieses, double sampling_rate);
    QVector<QPointF> display_peaks_in_window(const QString& label, std::vector<int>* outPostTags = nullptr) const;

    QVector<QPointF> detectPeaks(const QString& label,
        double detStart, double detEnd,
        std::vector<int>* outPostTags = nullptr) const;

    QVector<QPointF> get_bpm(const QString& label, double& outDuration) const;

    std::pair<double, double> statsWindow(double detStart, double detEnd) const;
    void ampogram(double sampling_length = 15);
    void updateAmpogramCursor();
    void syncChunkScrollBar();
    void setupHypnogram();
    void updateNoiseHighlights();
    void finalizeMarking(QChartView* cv, double endX, const QString& signalLabel);
    void cancelMarking(const QString& signalLabel);
    void restoreMarkingMarkers();
    void showStartMarker(QChartView* cv, double xValue, ChannelMarkingState& state, const QColor& color, QPushButton* stopBtn);
    void clearStartMarker(ChannelMarkingState& state);
    void updateDragPreview(QChartView* cv, double x0, double x1, const QColor& color);
    void clearDragPreview();
    void loadSelectedFile(const QString& filePath);

    bool handleMousePress(QChartView* cv, QWidget* viewport, QMouseEvent* event);
    double yScaleForSignal(const QString& label) const;
    bool   invertedForSignal(const QString& label) const;

    void resetUnpinnedGains();

    QVector<ParamOverride> m_thresholdOverrides;
    QVector<ParamOverride> m_blankingOverrides;

    void   finalizeParamEdit(const QStringList& channels, double globalStart, double globalEnd);
    void   commitMarkingSpan(const QString& clickedLabel, double globalStart, double globalEnd);
    double thresholdAt(const QString& label, double globalTime) const;
    double blankingAt(const QString& label, double globalTime) const;
};