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
    bool invertedForSignal(const QString& label) const;
    void set_params_to_config_defaults(const config_entry& cfg) {
        m_cfg = cfg;        // Set the default values for the threshold and blanking period spinboxes based on the config entry.
		ui->notch_filter->setEnabled(m_cfg.notch_filter_hz != 0); //enable notch filter checkbox if the config entry has a non-zero notch filter frequency
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
    void onFixScaleToggled(const QString& label, bool on);
    double scrollStepSeconds() const;
    config_entry m_cfg;
    beat_log* m_beatLog = nullptr;
    QTimer* m_logFlushTimer = nullptr;   // flushes the beat log to disk every 30 s
    QHash<QString, QPair<double, double>> m_fixedYRange; // fixed y range for scaling

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
        double sampleRate = 0.0;
        double nativeRate = 0.0;   ///< config native rate (Hz), for the 'Original Frequency' label
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
    ChannelMarkingState  mark_state_ecg1;
    ChannelMarkingState  mark_state_ecg2;
    ChannelMarkingState  mark_state_ecg3;
    ChannelMarkingState  mark_state_ppg;
    ChannelMarkingState  mark_state_abp;
    ChannelMarkingState  mark_state_accel;
    ChannelMarkingState  mark_state_art;
    ChannelMarkingState  mark_state_art_pulm;

    std::unique_ptr<Ui::noise_marking_gui> ui;
    std::unique_ptr<annotation_handler>       m_noiseManager;
    std::unique_ptr<user_control_handler>     m_buttonHandler;
    std::unique_ptr<pulse_overlay>            m_pulseOverlay;
    std::unique_ptr<annotation_eraser>        m_annotationEraser;

    QSet<QString>               m_activeChannels;
    QString                     m_currentMarkingType;
    GenExcStruct                m_genExc;
    QMap<QString, GenExcStruct> m_fileMarkings;  ///< stashed markings per file path

    // --- Sampling rates ---
    double m_sleepSR = 0.0;   ///< sleep epochs per second (= 1 / sleep_epoch_length)

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
    bool m_notchFilterEnabled = false;

    // --- Drag state ---
    bool     m_isDragging = false;
    QWidget* m_draggedViewport = nullptr;
    QPoint   m_dragStartPos;
    QString  m_dragSignalLabel;
    QMap<QChartView*, QAreaSeries*> m_dragPreviews;

    // --- Signal data (upsampled, 1 kHz unless otherwise noted) ---
    QVector<double> m_ecg1, m_ecg2, m_ecg3, m_ppg, m_accelX, m_accelY, m_accelZ, m_resp;
    QVector<double> m_cvp, m_abp, m_temp, m_marker, m_pacemaker, m_sleepStages, m_art, m_artPulm;


    QVector<QPointF> m_ecg1Raw, m_ecg2Raw, m_ecg3Raw, m_ppgRaw, m_abpRaw;
    QVector<QPointF> m_accelXRaw, m_accelYRaw, m_accelZRaw, m_respRaw, m_cvpRaw;
    QVector<QPointF> m_tempRaw, m_markerRaw, m_pacemakerRaw, m_artRaw, m_artPulmRaw;
    QString m_binFilePath;

    static constexpr uint32_t BIN_HEADER_VERSION = 1;   // header format version (offset 0)
    static constexpr qint64 FILE_HEADER_SIZE = 592;     // version+n_channels+sleep + 4*36 arrays + sleep_count
    static constexpr int NUM_CHANNELS = 36;

    uint32_t upsampled_channel_sizes[NUM_CHANNELS] = {};
    uint32_t raw_channel_sizes[NUM_CHANNELS] = {};
    float    channel_native_rates[NUM_CHANNELS] = {};
    float channel_upsampled_rates[NUM_CHANNELS] = {};
    uint32_t total_sleep_samples = 0;

    // Channel indices. Slot 0 is the new Timestamp channel; every other
    // channel has shifted +1 from the v1 layout. CH_RESERVED_40 is a
    // reserved tail slot (always absent on disk).
    enum ChannelIdx {
        CH_TIMESTAMP = 0,
        CH_ECG1, CH_ECG2, CH_ECG3, CH_PPG,
        CH_ACCEL_X, CH_ACCEL_Y, CH_ACCEL_Z,
        CH_MARKER, CH_TEMP, CH_PACEMAKER_EVENT,
        CH_EOG_L, CH_EOG_R, CH_EMG,
        CH_EEG1, CH_EEG2, CH_EEG3, CH_EEG4,
        CH_CVP, CH_PRES, CH_FLOW, CH_SNORE, CH_THOR, CH_ABDO,
        CH_LEG, CH_AUXAC, CH_THERM, CH_POS,
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
    bool loadChunkFromFile(uint64_t chunkIndex, bool resetScroll = true);
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
    void restoreMarkingFeatureMarkers();
    void showStartMarker(QChartView* cv, double xValue, ChannelMarkingState& state, const QColor& color, QPushButton* stopBtn);
    void clearStartMarker(ChannelMarkingState& state);
    void updateDragPreview(QChartView* cv, double x0, double x1, const QColor& color);
    void clearDragPreview();
    void loadSelectedFile(const QString& filePath);

    bool handleMousePress(QChartView* cv, QWidget* viewport, QMouseEvent* event);
    double yScaleForSignal(const QString& label) const;

    void resetUnpinnedGains();

    QVector<ParamOverride> m_thresholdOverrides;
    QVector<ParamOverride> m_blankingOverrides;
    QVector<ParamOverride> m_invertOverrides;

    bool invertedAt(const QString& label, double globalTime) const;
    void applyInvertOverride(const QStringList& channels, double globalStart, double globalEnd);
    void   finalizeParamEdit(const QStringList& channels, double globalStart, double globalEnd);
    bool   editParamOverrideAt(QChartView* cv, const QPoint& pos);
    bool   promptThresholdBlanking(const QString& header, double& thr, double& blk);
    void   applyParamOverrides(const QStringList& channels, double lo, double hi, double thrVal, double blkVal);
    void   commitMarkingSpan(const QString& clickedLabel, double globalStart, double globalEnd);
    double thresholdAt(const QString& label, double globalTime) const;
    double blankingAt(const QString& label, double globalTime) const;
};