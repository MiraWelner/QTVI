/**
 * @file    gui_handler.hpp
 * @brief   Main dialog for viewing ECG/PPG signals and annotating noise,
 *          arrhythmias, and artifacts. Loads data in 8-hour chunks from a
 *          binary file produced by file_to_bin.
 *
 * @details Expects the paired-channel .bin format (512-byte header):
 *            - 4 rate fields (signal, boolean, pacemaker, sleep-epoch-length)
 *            - 40 upsampled-block sizes
 *            - 40 raw-block sizes (each counts (timestamp, value) PAIRS,
 *              so the on-disk byte length is 2 * size * sizeof(double))
 *            - 40 native sampling rates (float32; 0 = absent)
 *            - 1 sleep-sample count
 *          Channel 0 is a Timestamp channel (seconds from recording start).
 *          It's loaded alongside the rest but not plotted -- the GUI only
 *          plots channels it already plotted, and the timestamp slot is
 *          treated like any other unplotted channel.
 *          For each of the 40 channels, the upsampled block is written
 *          first, followed by the raw block of (t, v) pairs in seconds
 *          from start of recording. For markable channels (ECG1/2/3, PPG,
 *          ABP) the raw block is overlaid on the plot as a grayscale
 *          scatter at the true sample timestamps. Storing timestamps
 *          (instead of assuming uniform spacing) lets irregular and
 *          gap-filled data (CHAOS .dat) render correctly.
 *
 * @author  Mira Welner
 * @email   MEW386@pitt.edu
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
#include <QHash>
#include <QColor>
#include <QPointF>
#include <QTimer>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "ui_noise_marking_gui.h"
#include "user_annotation_handler.h"
#include "user_control_handler.h"
#include "config_entry.hpp"
#include "grid_overlay.hpp"
#include "gap_indicator.hpp"

 /**
  * @brief Collected annotations returned to the caller after the dialog closes.
  *
  *        All three parallel lists (noiseExc, data_type, marking_type) share
  *        the same index space -- entry `i` describes one marked segment:
  *        its [start, end] global-time range, which channel it was marked on,
  *        and what kind of marking it is.
  */
struct GenExcStruct {
    QString filePath;                          ///< Source file these markings belong to.
    QVector<QPair<double, double>> noiseExc;   ///< [start, end] in global seconds.
    QStringList data_type;                     ///< "ECG1", "ECG2", "ECG3", "PPG", or "ABP".
    QStringList marking_type;                  ///< "Noise/Artifact", "AF", "SVT", ...
};

class beat_log;   // defined in beat_log.hpp; only a pointer is held here

class noise_marking_gui : public QDialog {
    Q_OBJECT
        friend class user_control_handler;
    friend class gap_indicator;

public:
    /**
     * @param parent Optional Qt parent widget.
     */
    explicit noise_marking_gui(QWidget* parent = nullptr);
    ~noise_marking_gui() override;

    /**
     * @brief  Markings recorded for the *currently loaded* file only.
     *         For cross-file markings use getAllMarkings().
     * @return GenExcStruct populated with this file's markings (filePath set).
     */
    GenExcStruct getMarkings() const;

    /**
     * @brief  Markings for every file touched during this dialog session.
     * @return One GenExcStruct per file that has at least one marking;
     *         stashed markings from previously-viewed files are included
     *         alongside the currently-loaded file's live markings.
     */
    QVector<GenExcStruct> getAllMarkings() const;

    /**
     * @return Path of the .bin currently loaded (empty if none).
     */
    QString getFilePath() const { return m_binFilePath; }

    /**
     * @brief  Load a .bin file produced by file_to_bin.
     * @param  filePath Absolute or relative path to the .bin.
     *
     * @note   Safe to call more than once; markings for the previous file
     *         are stashed and can be retrieved via getAllMarkings().
     */
    void setFileSource(const QString& filePath);

    double currentStartTime() const { return m_currentStartTime; }
    double windowDuration()   const { return m_windowDuration; }

    QChartView* chartViewForSignalLabel(const QString& label) const;
    bool        isChannelActive(const QString& label) const;
    enum class MarkPhase { Idle, WaitingForStart, WaitingForEnd, WaitingForStop };
    enum class PlotMode { Line, Scatter };
    enum class ParamEdit { None, Threshold, Blanking };
    void setBeatLog(beat_log* log) { m_beatLog = log; }


    void set_params_to_config_defaults(const config_entry& cfg) {
        // Set the default values for the threshold and blanking period spinboxes based on the config entry.
        m_cfg = cfg;
    }

protected:
    bool eventFilter(QObject* watched, QEvent* event) override;

    // Clicking anywhere on the dialog (outside a focused spinbox) steals
    // focus from that spinbox so its editingFinished fires and the
    // clearFocus lambda runs. Without this, ClickFocus-policy spinboxes
    // keep focus indefinitely because clicking on a chart -- which has
    // NoFocus policy -- doesn't transfer focus.
    void mousePressEvent(QMouseEvent* event) override;

private slots:
    void on_skip_interval_box_returnPressed();
    void on_skip_interval_box_editingFinished();
    void on_marking_type_currentTextChanged(const QString& text);
    void on_next8hours_clicked();
    void on_prev8hours_clicked();
    void handleBrowseFile();

private:
    config_entry m_cfg;
    beat_log* m_beatLog = nullptr;
    QTimer* m_logFlushTimer = nullptr;   // flushes the beat log to disk every 30 s


    /**
    * @brief Per-channel state for the marking state machine.
    */
    struct ChannelMarkingState {
        MarkPhase    phase = MarkPhase::Idle;
        double       globalStartTime = 0.0;   ///< Always in global seconds (not chunk-local).
        QLineSeries* startMarkerLine = nullptr;
    };

    // --- Called by lower_row_buttons to enter/transition marking phases ---

    /** @brief Enter WaitingForStart on `signalLabel`'s state machine. */
    void beginMarking(const QString& signalLabel);

    /** @brief Transition `signalLabel` from WaitingForEnd to WaitingForStop. */
    void beginStopPhase(const QString& signalLabel);

    /** @brief Enter WaitingForStart on every active markable channel. */
    void beginMarkingAll();

    /** @brief Transition every active channel from WaitingForEnd to WaitingForStop. */
    void beginStopPhaseAll();

    /** @brief Enter WaitingForStart on ECG1/2/3 only (the "All ECG" buttons). */
    void beginMarkingEcgAll();

    /** @brief Transition ECG1/2/3 from WaitingForEnd to WaitingForStop. */
    void beginStopPhaseEcgAll();

    /*
        The data_channel_features has all the features that
        a single one of the markable channels might have, such as the upsampled
        data, raw data, and R peak threshold
    */
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

    /**
     * @brief  Lookup bundle of all per-channel references.
     * @param  label Channel label ("ECG1", "ECG2", "ECG3", "PPG", or "ABP").
     * @return Populated ChannelRefs (default-constructed fields if label unknown).
     */
    data_channel_features channelRefs(const QString& label) const;

    /// Labels of all markable channels, in display order.
    static const QStringList& markableChannelLabels();

    // --- Per-channel marking state ---
    ChannelMarkingState  m_markState_ecg1;
    ChannelMarkingState  m_markState_ecg2;
    ChannelMarkingState  m_markState_ecg3;
    ChannelMarkingState  m_markState_ppg;
    ChannelMarkingState  m_markState_abp;

    // --- Core components ---
    std::unique_ptr<Ui::noise_marking_gui> ui;
    std::unique_ptr<annotation_handler>          m_noiseManager;
    std::unique_ptr<user_control_handler>     m_buttonHandler;
    std::unique_ptr<pulse_overlay>            m_pulseOverlay;
    std::unique_ptr<gap_indicator>            m_gapIndicator;

    QSet<QString>               m_activeChannels;
    QString                     m_currentMarkingType;
    GenExcStruct                m_genExc;
    QMap<QString, GenExcStruct> m_fileMarkings;  ///< stashed markings per file path

    // --- Sampling rates (upsampled signals) ---
    double m_ecgSR = 0.0;   ///< TARGET_RATE for upsampled ECG/PPG/ABP (1 kHz per file_to_bin)
    double m_ppgSR = 0.0;
    double m_boolSR = 0.0;   ///< 1 Hz channels
    double m_sleepSR = 0.0;   ///< sleep epochs per second (= 1 / sleep_epoch_length)

    // (Native raw sample rates are no longer needed: each raw sample carries
    // its own timestamp, so the plotter doesn't have to infer spacing.)

    // --- View state ---
    double m_currentStartTime = 0.0;
    double m_windowDuration = 10.0;
    double m_skipInterval = 5.0;

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

    // --- Persistent line series, keyed by chart view ---
    //
    // Each chart has one QLineSeries per logical series slot, allocated
    // lazily on first render and kept alive for the lifetime of the dialog.
    // renderWindowedChart() reuses these instead of allocating+deleting per
    // frame -- the OpenGL teardown in ~QLineSeries was costing tens of
    // seconds on a Line->Scatter mode switch, because every chart had a
    // dense GL-backed line that had to release GPU buffers.
    //
    // In Scatter mode the lines are hidden via setVisible(false) but not
    // freed; switching back to Line mode just toggles visibility and calls
    // replace() with the new window's data.
    QHash<QChartView*, QList<QLineSeries*>> m_persistentLines;

    // --- Mark-all state ---
    bool m_markAllActive = false;   ///< true while "Mark All Signals" mode is in effect

    // --- Plot style (global, applies to all signal charts) ---
    PlotMode m_plotMode = PlotMode::Line;

    // checked = per-window autoscale (drift hidden)
    bool m_filterBaselineDrift = false;

    // --- Drag state ---
    bool     m_isDragging = false;
    QWidget* m_draggedViewport = nullptr;
    QPoint   m_dragStartPos;
    QString  m_dragSignalLabel;

    // --- Signal data (upsampled, 1 kHz unless otherwise noted) ---
    QVector<double> m_ecg1, m_ecg2, m_ecg3, m_ppg;
    QVector<double> m_accelX, m_accelY, m_accelZ;
    QVector<double> m_resp;        ///< respiration (slot 34)
    QVector<double> m_cvp;         ///< central venous pressure (slot 16, Bittium)
    QVector<double> m_abp;         ///< arterial blood pressure (slot 35, Bittium)
    QVector<double> m_sleepStages;

    // --- Raw data, plotted as scatter overlay on markable signals
    //     (ECG/PPG/ABP) and on the accel chart. Each element is (t, v)
    //     where t is chunk-local seconds. The loader filters the on-disk
    //     pair stream to this chunk. ---
    QVector<QPointF> m_ecg1Raw, m_ecg2Raw, m_ecg3Raw, m_ppgRaw, m_abpRaw;
    QVector<QPointF> m_accelXRaw, m_accelYRaw, m_accelZRaw, m_respRaw, m_cvpRaw;

    // ------------------------------------------------------------------------
    //  File layout (v2):
    //    512-byte header (128 x uint32-sized fields) =
    //      4 rates + 40 upsampled-sizes + 40 raw-sizes + 40 native-rate
    //      floats + 1 sleep-size
    //    Then, for each of the 40 channels, in order:
    //      upsampled block (m_chanSizes[i] doubles)
    //      raw block       (m_chanSizesRaw[i] * 2 doubles, interleaved (t, v))
    //    Then sleep stages (m_totalSleepSamples doubles).
    //
    //    m_chanSizesRaw[i] counts PAIRS, not individual doubles -- so the
    //    raw block's byte length is m_chanSizesRaw[i] * 2 * sizeof(double).
    //    m_chanNativeRates[i] holds the channel's original sampling rate
    //    in Hz, or 0.0 if the channel is absent. The GUI doesn't consume
    //    these (all plotting uses m_ecgSR/m_ppgSR/m_boolSR); they're just
    //    parsed so the member stays in sync with the file format.
    // ------------------------------------------------------------------------
    QString m_binFilePath;
    static constexpr qint64 FILE_HEADER_SIZE = 500;
    static constexpr int NUM_CHANNELS = 40;

    uint32_t m_chanSizes[NUM_CHANNELS] = {};      ///< upsampled-block sizes
    uint32_t m_chanSizesRaw[NUM_CHANNELS] = {};   ///< raw-block sizes (PAIRS)
    float    m_chanNativeRates[NUM_CHANNELS] = {};///< native sample rate (Hz), 0 = absent
    uint32_t m_totalSleepSamples = 0;

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
        CH_PRES, CH_FLOW, CH_THOR, CH_ABDO,
        CH_LEG, CH_THERM, CH_POS,
        CH_EKG_OFF, CH_EOG_L_OFF, CH_EOG_R_OFF, CH_EMG_OFF,
        CH_EEG1_OFF, CH_EEG2_OFF, CH_EEG3_OFF,
        CH_OXSTATUS, CH_SPO2, CH_HR, CH_DHR,
        CH_RESP, CH_ABP,
        CH_ART, CH_ART_PULM
    };

    uint64_t m_currentChunkIndex = 0;

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

    /** @brief Refresh Start/Stop button styling for a single channel. */
    void updateButtonStatesForChannel(const QString& label);

    /** @brief Refresh Start/Stop button styling for all channels + the mark-all buttons. */
    void updateAllChannelButtonStates();

    /** @brief Highlight whichever parameter-edit button (threshold/blanking) is
 *         currently armed; clear both when no edit mode is active. */
    void updateParamButtonStyles();

    // --- Plotting ---
    /**
     * @brief  Load the requested 8-hour chunk into the per-channel vectors.
     * @param  chunkIndex 0-based chunk index (0 = hours 0-8, 1 = hours 8-16, ...).
     * @return true on success, false if the file could not be opened.
     */
    bool loadChunkFromFile(uint64_t chunkIndex);

    /** @brief Redraw the 10-second signal plots for the current view window. */
    void handle_data_plot();

    /** @brief Run the simple peak finder on `label`'s raw data restricted
    *      to the current view window. Returns chunk-local (t, v) of
    *      detected peaks. Empty when m_showPeaks is false. */
    QVector<QPointF> display_peaks_in_window(const QString& label,
        std::vector<std::string>* outPostTags = nullptr) const;

    QVector<QPointF> detectPeaks(const QString& label,
        double detStart, double detEnd,
        std::vector<std::string>* outPostTags = nullptr) const;

    /** @brief Run the simple peak finder over a BPM-estimation window: the
    *      visible window, extended backwards to a minimum of 10 s when
    *      the visible window is shorter than that. Clamped to the start
    *      of the chunk. Returns chunk-local (t, v) of detected peaks and,
    *      via @p outDuration, the actual duration of the window used
    *      (so the caller can divide peak count by it). Empty when
    *      m_showPeaks is false. */
    QVector<QPointF> get_bpm(const QString& label,
        double& outDuration) const;

    std::pair<double, double> statsWindow(double detStart, double detEnd) const;

    /**
     * @brief  Redraw the 8-hour amplitude overviews and raw-signal overviews.
     * @param  sampling_length Bucket width in seconds for amplitude decimation.
     */
    void ampogram(double sampling_length = 15);

    /** @brief Redraw the black vertical cursor bar on each overview chart. */
    void updateAmpogramCursor();

    /** @brief Build / rebuild the hypnogram (sleep-stage) chart from m_sleepStages. */
    void setupHypnogram();

    /** @brief Overlay translucent colored rectangles for recorded marking segments. */
    void updateNoiseHighlights();

    // --- Marking ---
    /**
     * @brief  Commit a marking that is currently in WaitingForStop.
     * @param  cv         Chart view the marking was placed on.
     * @param  endX       Chunk-local x (seconds) of the end click.
     * @param  signalLabel Channel the marking is for.
     */
    void finalizeMarking(QChartView* cv, double endX, const QString& signalLabel);

    /** @brief Abort an in-progress marking and clear its visual state. */
    void cancelMarking(const QString& signalLabel);

    /** @brief Restore the dashed start markers after a chunk/file change. */
    void restoreMarkingMarkers();

    /**
     * @brief Draw the dashed vertical start-marker line on a chart.
     * @param cv      Target chart view.
     * @param xValue  Chunk-local x (seconds).
     * @param state   Per-channel marking state to update.
     * @param color   Line color (matches the channel's signal color).
     * @param stopBtn Stop button to enable once the marker is shown.
     */
    void showStartMarker(QChartView* cv, double xValue, ChannelMarkingState& state,
        const QColor& color, QPushButton* stopBtn);

    /** @brief Remove the dashed start-marker line for a channel (if any). */
    void clearStartMarker(ChannelMarkingState& state);

    // --- File selection ---
    /**
     * @brief  Load a .bin from disk, replacing the current view.
     * @param  filePath Path to the .bin.
     *
     * @note   Any in-progress markings on the previous file are stashed in
     *         m_fileMarkings and can be retrieved with getAllMarkings().
     */
    void loadSelectedFile(const QString& filePath);

    // --- Formatting ---
    /**
     * @param  seconds Offset from chunk start.
     * @return "HH:MM:SS.ss"-formatted string.
     */
    QString formatTimeLabel(double seconds);


    bool handleMousePress(QChartView* cv, QWidget* viewport, QMouseEvent* event);

    /** @return Y-axis multiplier for a signal: when the chart's "Fix Scale"
 *          checkbox is checked, returns the spinbox value; otherwise 1.
 *          Multiplier > 1 widens the visible y range (signal looks shorter). */
    double yScaleForSignal(const QString& label) const;

    void resetUnpinnedGains();

    ParamEdit m_paramEditMode = ParamEdit::None;
    double    m_paramDragStartGlobal = 0.0;

    struct ParamOverride {
        QString channel;          // "ECG1".."ABP"
        double  start = 0.0;      // global seconds
        double  end = 0.0;
        double  value = 0.0;
    };
    QVector<ParamOverride> m_thresholdOverrides;
    QVector<ParamOverride> m_blankingOverrides;

    void   enterParamEdit(ParamEdit which);
    void   finalizeParamEdit(const QString& label, double globalStart, double globalEnd);
    double thresholdAt(const QString& label, double globalTime) const;
    double blankingAt(const QString& label, double globalTime) const;
};