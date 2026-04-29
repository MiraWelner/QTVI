/**
 * @file    gui_handler.hpp
 * @brief   Main dialog for viewing ECG/PPG signals and annotating noise,
 *          arrhythmias, and artifacts. Loads data in 8-hour chunks from a
 *          binary file produced by file_to_bin.
 *
 * @details Expects the paired-channel .bin format (512-byte header):
 *            - 4 rate fields (signal, boolean, pacemaker, sleep-epoch-length)
 *            - 41 upsampled-block sizes
 *            - 41 raw-block sizes (each counts (timestamp, value) PAIRS,
 *              so the on-disk byte length is 2 * size * sizeof(double))
 *            - 41 native sampling rates (float32; 0 = absent)
 *            - 1 sleep-sample count
 *          Channel 0 is a Timestamp channel (seconds from recording start).
 *          It's loaded alongside the rest but not plotted -- the GUI only
 *          plots channels it already plotted, and the timestamp slot is
 *          treated like any other unplotted channel.
 *          For each of the 41 channels, the upsampled block is written
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
#include <QColor>
#include <QPointF>
#include <memory>

#include "ui_noise_marking_gui.h"
#include "NoiseManager.h"
#include "lower_row_buttons.h"

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

class noise_marking_gui : public QDialog {
    Q_OBJECT
        friend class lower_row_buttons;

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

protected:
    bool eventFilter(QObject* watched, QEvent* event) override;

private slots:
    void on_skip_interval_box_returnPressed();
    void on_skip_interval_box_editingFinished();
    void on_marking_type_currentTextChanged(const QString& text);
    void on_next8hours_clicked();
    void on_prev8hours_clicked();
    void handleBrowseFile();

public:
    enum class MarkPhase { Idle, WaitingForStart, WaitingForEnd, WaitingForStop };
    enum class PlotMode { Line, Scatter };

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

private:
    // ------------------------------------------------------------------------
    //  Channel lookup table
    //
    //  The five markable channels (ECG1/2/3, PPG, ABP) have parallel:
    //      * signal-data vector (upsampled)
    //      * raw-data vector (native rate, overlaid as grayscale scatter)
    //      * chart view
    //      * start / stop buttons
    //      * marking-state machine
    //      * color
    //
    //  Previously these were 5 independent 5-way switches. `ChannelRefs`
    //  collapses them into one lookup, accessed via `channelRefs(label)`.
    // ------------------------------------------------------------------------
    struct ChannelRefs {
        QChartView* chartView = nullptr;
        QPushButton* startButton = nullptr;
        QPushButton* stopButton = nullptr;
        ChannelMarkingState* state = nullptr;
        const QVector<double>* data = nullptr;   ///< upsampled samples (1 kHz)
        const QVector<QPointF>* dataRaw = nullptr;   ///< raw (t, v) pairs, chunk-local seconds
        const double* sampleRate = nullptr;   ///< rate for `data`
        QColor                    color;
    };

    /**
     * @brief  Lookup bundle of all per-channel references.
     * @param  label Channel label ("ECG1", "ECG2", "ECG3", "PPG", or "ABP").
     * @return Populated ChannelRefs (default-constructed fields if label unknown).
     */
    ChannelRefs channelRefs(const QString& label) const;

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
    std::unique_ptr<NoiseManager>          m_noiseManager;
    std::unique_ptr<lower_row_buttons>     m_buttonHandler;

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

    // --- Mark-all state ---
    bool m_markAllActive = false;   ///< true while "Mark All Signals" mode is in effect

    // --- Plot style (global, applies to all signal charts) ---
    PlotMode m_plotMode = PlotMode::Line;

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
    QVector<QPointF> m_accelXRaw, m_accelYRaw, m_accelZRaw;

    // ------------------------------------------------------------------------
    //  File layout (v2):
    //    512-byte header (128 x uint32-sized fields) =
    //      4 rates + 41 upsampled-sizes + 41 raw-sizes + 41 native-rate
    //      floats + 1 sleep-size
    //    Then, for each of the 41 channels, in order:
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
    static constexpr qint64 FILE_HEADER_SIZE = 512;
    static constexpr int NUM_CHANNELS = 41;

    uint32_t m_chanSizes[NUM_CHANNELS] = {};      ///< upsampled-block sizes
    uint32_t m_chanSizesRaw[NUM_CHANNELS] = {};   ///< raw-block sizes (PAIRS)
    float    m_chanNativeRates[NUM_CHANNELS] = {};///< native sample rate (Hz), 0 = absent
    uint32_t m_totalSleepSamples = 0;

    // Channel indices. Slot 0 is the new Timestamp channel; every other
    // channel has shifted +1 from the v1 layout. CH_RESERVED_40 is a
    // reserved tail slot (always absent on disk).
    static constexpr int CH_TIMESTAMP = 0;
    static constexpr int CH_ECG1 = 1, CH_ECG2 = 2, CH_ECG3 = 3, CH_PPG = 4;
    static constexpr int CH_ACCEL_X = 5, CH_ACCEL_Y = 6, CH_ACCEL_Z = 7;
    static constexpr int CH_MARKER = 8, CH_TEMP = 9, CH_PACEMAKER = 10;
    static constexpr int CH_EOG_L = 11, CH_EOG_R = 12, CH_EMG = 13;
    static constexpr int CH_EEG1 = 14, CH_EEG2 = 15, CH_EEG3 = 16;
    static constexpr int CH_PRES = 17, CH_FLOW = 18, CH_THOR = 19;
    static constexpr int CH_ABDO = 20, CH_LEG = 21, CH_THERM = 22;
    static constexpr int CH_POS = 23;
    static constexpr int CH_EKG_OFF = 24, CH_EOGL_OFF = 25, CH_EOGR_OFF = 26;
    static constexpr int CH_EMG_OFF = 27, CH_EEG1_OFF = 28, CH_EEG2_OFF = 29;
    static constexpr int CH_EEG3_OFF = 30;
    static constexpr int CH_OXSTATUS = 31, CH_SPO2 = 32;
    static constexpr int CH_HR = 33, CH_DHR = 34, CH_RESP = 35;
    static constexpr int CH_ABP = 36, CH_EEG4 = 37;
    static constexpr int CH_ART = 38, CH_ART_PULM = 39;
    static constexpr int CH_RESERVED_40 = 40;

    uint64_t m_currentChunkIndex = 0;

    static constexpr double CHUNK_DURATION_SEC = 28800.0;  ///< 8 hours

    // --- Lookup helpers (thin wrappers over channelRefs) ---
    ChannelMarkingState& markStateFor(const QString& label);
    QString              signalLabelForChartView(QChartView* cv) const;
    QChartView* chartViewForSignalLabel(const QString& label) const;
    double               sampleRateForSignal(const QString& label) const;
    QColor               colorForSignal(const QString& label) const;
    bool                 isChannelActive(const QString& label) const;

    /** @return Duration in seconds of the currently-loaded 8-hour chunk. */
    double               totalChunkDuration() const;

    // --- Per-channel button helpers ---
    QPushButton* startButtonForSignal(const QString& label) const;
    QPushButton* stopButtonForSignal(const QString& label) const;

    /** @brief Refresh Start/Stop button styling for a single channel. */
    void updateButtonStatesForChannel(const QString& label);

    /** @brief Refresh Start/Stop button styling for all channels + the mark-all buttons. */
    void updateAllChannelButtonStates();

    // --- Plotting ---
    /**
     * @brief  Load the requested 8-hour chunk into the per-channel vectors.
     * @param  chunkIndex 0-based chunk index (0 = hours 0-8, 1 = hours 8-16, ...).
     * @return true on success, false if the file could not be opened.
     */
    bool loadChunkFromFile(uint64_t chunkIndex);

    /** @brief Redraw the 10-second signal plots for the current view window. */
    void handle_data_plot();

    /**
     * @brief  Redraw the 8-hour amplitude overviews and raw-signal overviews.
     * @param  sampling_length Bucket width in seconds for amplitude decimation.
     */
    void handle_ampogram_plot(double sampling_length = 60);

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
};