/**
 * @file    gui_handler.cpp
 * @brief   Implementation of the noise marking GUI. Displays ECG1/ECG2/ECG3,
 *          PPG, ABP, accel, RESP, CVP, and sleep-stage data and lets the user
 *          mark noise / arrhythmia segments on the markable channels; outputs
 *          are a CSV and a .bin written by the caller (see noise_marking_gui.cpp).
 *
 * @details The signal charts overlay two views of each markable channel:
 *            - The *raw* native-rate samples drawn as a solid black scatter
 *              overlay, always. In Scatter mode, this is the ONLY thing
 *              shown (no upsampled foreground).
 *            - In Line mode, the *upsampled* samples (1 kHz) are drawn as
 *              a purple line underneath the raw black scatter.
 *
 *          Input .bin format: 512-byte header (128 x uint32-sized fields):
 *            Offset  0:   signal_rate    (uint32) -- upsampled rate (1 kHz)
 *            Offset  4:   boolean_rate   (uint32) -- 1 Hz channels
 *            Offset  8:   pacemaker_rate (uint32)
 *            Offset 12:   sleep_rate     (uint32) -- epoch length in seconds
 *            Offset 16:   41 x upsampled-block sizes     (uint32)
 *            Offset 180:  41 x raw-block sizes           (uint32)
 *            Offset 344:  41 x native sampling rates     (float32, Hz; 0 = absent)
 *            Offset 508:  sleep_sample_count             (uint32)
 *          Signal data: for each of the 41 channels, upsampled block then
 *          raw block (both doubles), followed by sleep-stage doubles.
 *          Slot 0 is a Timestamp channel; this GUI parses it like any other
 *          slot but doesn't plot it (matches "ignored channels" behavior).
 *          Missing channels are stored as a single -1.0 (size = 1).
 *
 * @author  Mira Welner
 * @email   MEW386@pitt.edu
 * @date    2026-03-22
 */
#include "gui_handler.hpp"

#include <QtCharts/QAreaSeries>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QXYSeries>
#include <QtCharts/QValueAxis>
#include <QtCharts/QCategoryAxis>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QShortcut>
#include <QButtonGroup>
#include <QRadioButton>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QtCharts/QLegendMarker>
#include <algorithm>
#include <cstring>

 // ============================================================================
 // Constants
 // ============================================================================

static const QColor COLOR_ECG1 = QColor("#BF00FF");
static const QColor COLOR_ECG2 = QColor("#BF00FF");
static const QColor COLOR_ECG3 = QColor("#BF00FF");
static const QColor COLOR_PPG = QColor("#BF00FF");
static const QColor COLOR_ABP = QColor("#BF00FF");
static const QColor COLOR_ACCEL_X = QColor("#F39C12");  // orange
static const QColor COLOR_ACCEL_Y = QColor("#27AE60");  // green
static const QColor COLOR_ACCEL_Z = QColor("#8E44AD");  // purple
static const QColor COLOR_RESP = QColor("#16A085");  // teal
static const QColor COLOR_CVP = QColor("#2980B9");  // blue

/// Mid-gray @ 50% alpha for the raw-sample scatter overlay on markable signals.
static const QColor COLOR_RAW_SCATTER = QColor(0, 0, 0, 255);

static const QMap<QString, QColor> MARKING_COLORS = {
    {"Noise/Artifact",         QColor(255, 255, 0,   30)},
    {"Conduction Delay",       QColor(128, 0,   128, 30)},
    {"AF",                     QColor(255, 0,   0,   30)},
    {"SVT",                    QColor(0,   255, 0,   60)},
    {"VT",                     QColor(0,   0,   255, 60)},
    {"PVC",                    QColor(128, 255, 0,   60)},
    {"PAC",                    QColor(255, 128, 0,   60)},
    {"Benign Arrhythmia",      QColor(255, 128, 255, 60)},
    {"Significant Arrhythmia", QColor(0,   255, 255, 60)}
};

// ============================================================================
// Anonymous helpers
// ============================================================================

namespace {

    void clearAxes(QChart* chart) {
        if (!chart) return;
        for (QAbstractAxis* axis : chart->axes()) chart->removeAxis(axis);
    }

    // Format the chart's title text with the channel name plus its native
    // (un-upsampled) sample rate and the current pixel/second resolution.
    // The title lives inside the chart itself, so it reflows automatically
    // as the chart resizes; we just reset the string on plot-area changes.
    QString formatChartTitle(const QString& signalName,
        double nativeHz, double pxPerSec) {
        return QString("%1  -- Original Frequency: %2 Hz -- Pixel Resolution: %3 px/s")
            .arg(signalName)
            .arg(nativeHz, 0, 'f', 1)
            .arg(pxPerSec, 0, 'f', 2);
    }

    void setupChartDefaults(QChartView* view) {
        auto* chart = new QChart();
        chart->legend()->hide();
        chart->setMargins(QMargins(0, 0, 0, 0));
        chart->setBackgroundRoundness(0);
        view->setChart(chart);
    }

    bool isMissingSignal(const QVector<double>& data) {
        return data.isEmpty() || (data.size() == 1 && data[0] == -1.0);
    }

}  // namespace

// ============================================================================
// Channel lookup table
//
// One central point that answers every per-channel question: what chart is
// this channel displayed on, which Start/Stop buttons drive it, which state
// machine tracks its marking phase, where are its samples and native rate,
// and what color is it drawn in. All other lookup helpers are thin wrappers
// around this.
// ============================================================================

/**
 * @return Ordered list of all markable channel labels.
 */
const QStringList& noise_marking_gui::markableChannelLabels() {
    static const QStringList kLabels{ "ECG1", "ECG2", "ECG3", "PPG", "ABP" };
    return kLabels;
}

noise_marking_gui::ChannelRefs noise_marking_gui::channelRefs(const QString& label) const {
    // We cast away const for the state pointer so callers like markStateFor()
    // can get a mutable reference; the ChannelRefs itself is a pure lookup
    // bundle (no ownership, no lifetime implications).
    auto* self = const_cast<noise_marking_gui*>(this);

    ChannelRefs r;
    if (label == "ECG1") {
        r.chartView = ui->ecg_axis_1;
        r.startButton = ui->start_ecg1_mark;
        r.stopButton = ui->stop_ecg1_mark;
        r.state = &self->m_markState_ecg1;
        r.data = &m_ecg1;
        r.dataRaw = &m_ecg1Raw;
        r.sampleRate = &m_ecgSR;
        r.color = COLOR_ECG1;
    }
    else if (label == "ECG2") {
        r.chartView = ui->ecg_axis_2;
        r.startButton = ui->start_ecg2_mark;
        r.stopButton = ui->stop_ecg2_mark;
        r.state = &self->m_markState_ecg2;
        r.data = &m_ecg2;
        r.dataRaw = &m_ecg2Raw;
        r.sampleRate = &m_ecgSR;
        r.color = COLOR_ECG2;
    }
    else if (label == "ECG3") {
        r.chartView = ui->ecg_axis_3;
        r.startButton = ui->start_ecg3_mark;
        r.stopButton = ui->stop_ecg3_mark;
        r.state = &self->m_markState_ecg3;
        r.data = &m_ecg3;
        r.dataRaw = &m_ecg3Raw;
        r.sampleRate = &m_ecgSR;
        r.color = COLOR_ECG3;
    }
    else if (label == "PPG") {
        r.chartView = ui->ppg_axis;
        r.startButton = ui->startNoisePPG;
        r.stopButton = ui->stopNoisePPG;
        r.state = &self->m_markState_ppg;
        r.data = &m_ppg;
        r.dataRaw = &m_ppgRaw;
        r.sampleRate = &m_ppgSR;
        r.color = COLOR_PPG;
    }
    else if (label == "ABP") {
        r.chartView = ui->accel_or_abg_axis;
        r.startButton = ui->startNoiseABP;
        r.stopButton = ui->stopNoiseABP;
        r.state = &self->m_markState_abp;
        r.data = &m_abp;
        r.dataRaw = &m_abpRaw;
        r.sampleRate = &m_ecgSR;
        r.color = COLOR_ABP;
    }
    return r;
}

// ----------------------------------------------------------------------------
//  Thin wrappers around channelRefs()
// ----------------------------------------------------------------------------

noise_marking_gui::ChannelMarkingState& noise_marking_gui::markStateFor(const QString& label) {
    ChannelRefs r = channelRefs(label);
    // Default to PPG's state if label is unknown (preserves original behavior).
    return r.state ? *r.state : m_markState_ppg;
}

QString noise_marking_gui::signalLabelForChartView(QChartView* cv) const {
    if (cv == ui->ecg_axis_1) return "ECG1";
    if (cv == ui->ecg_axis_2) return "ECG2";
    if (cv == ui->ecg_axis_3) return "ECG3";
    if (cv == ui->ppg_axis)   return "PPG";

    // accel_or_abg_axis is shared: it's "ABP" only when no accel channels
    // are present on the currently loaded file.
    if (cv == ui->accel_or_abg_axis) {
        bool anyAccel = !isMissingSignal(m_accelX)
            || !isMissingSignal(m_accelY)
            || !isMissingSignal(m_accelZ);
        if (!anyAccel && !isMissingSignal(m_abp)) return "ABP";
    }

    return {};
}

QChartView* noise_marking_gui::chartViewForSignalLabel(const QString& label) const {
    return channelRefs(label).chartView;
}

double noise_marking_gui::sampleRateForSignal(const QString& label) const {
    ChannelRefs r = channelRefs(label);
    return r.sampleRate ? *r.sampleRate : m_ecgSR;
}

QColor noise_marking_gui::colorForSignal(const QString& label) const {
    ChannelRefs r = channelRefs(label);
    return r.chartView ? r.color : COLOR_PPG;
}

bool noise_marking_gui::isChannelActive(const QString& label) const {
    return m_activeChannels.contains(label);
}

double noise_marking_gui::totalChunkDuration() const {
    if (m_ecg1.size() > 1 && m_ecgSR > 0) return m_ecg1.size() / m_ecgSR;
    if (m_ppg.size() > 1 && m_ppgSR > 0)  return m_ppg.size() / m_ppgSR;
    return 0.0;
}

// ============================================================================
// Per-channel button helpers
// ============================================================================

QPushButton* noise_marking_gui::startButtonForSignal(const QString& label) const {
    return channelRefs(label).startButton;
}

QPushButton* noise_marking_gui::stopButtonForSignal(const QString& label) const {
    return channelRefs(label).stopButton;
}

void noise_marking_gui::updateButtonStatesForChannel(const QString& label) {
    QPushButton* startBtn = startButtonForSignal(label);
    QPushButton* stopBtn = stopButtonForSignal(label);
    if (!startBtn || !stopBtn) return;

    bool active = isChannelActive(label);

    if (!active) {
        startBtn->setEnabled(false);
        stopBtn->setEnabled(false);
        startBtn->setStyleSheet("color: gray;");
        stopBtn->setStyleSheet("color: gray;");
        return;
    }

    // Channel is active — restore button appearance based on marking phase
    MarkPhase phase = markStateFor(label).phase;

    switch (phase) {
    case MarkPhase::WaitingForStart:
        startBtn->setEnabled(true);
        startBtn->setStyleSheet("background-color: #f39c12; color: white;");
        stopBtn->setEnabled(false);
        stopBtn->setStyleSheet("");
        break;
    case MarkPhase::WaitingForEnd:
        startBtn->setEnabled(true);
        startBtn->setStyleSheet("background-color: #f39c12; color: white;");
        stopBtn->setEnabled(true);
        stopBtn->setStyleSheet("");
        break;
    case MarkPhase::WaitingForStop:
        startBtn->setEnabled(true);
        startBtn->setStyleSheet("");
        stopBtn->setEnabled(true);
        stopBtn->setStyleSheet("background-color: #e74c3c; color: white;");
        break;
    case MarkPhase::Idle:
    default:
        startBtn->setEnabled(true);
        startBtn->setStyleSheet("");
        stopBtn->setEnabled(false);
        stopBtn->setStyleSheet("");
        break;
    }
}

void noise_marking_gui::updateAllChannelButtonStates() {
    for (const QString& label : markableChannelLabels())
        updateButtonStatesForChannel(label);

    bool anyActive = !m_activeChannels.isEmpty();
    ui->start_all_mark->setEnabled(anyActive);
    ui->start_all_mark->setStyleSheet(anyActive ? "" : "color: gray;");

    if (!anyActive || !m_markAllActive) {
        ui->stop_all_mark->setEnabled(false);
        ui->stop_all_mark->setStyleSheet(anyActive ? "" : "color: gray;");
        return;
    }

    // Mark-all is in progress — restore button styling based on channel phases
    bool anyWaitingEnd = false;
    bool anyWaitingStop = false;
    for (const QString& lbl : markableChannelLabels()) {
        if (!isChannelActive(lbl)) continue;
        MarkPhase p = markStateFor(lbl).phase;
        if (p == MarkPhase::WaitingForEnd) anyWaitingEnd = true;
        if (p == MarkPhase::WaitingForStop) anyWaitingStop = true;
    }

    if (anyWaitingStop) {
        ui->start_all_mark->setStyleSheet("");
        ui->stop_all_mark->setEnabled(true);
        ui->stop_all_mark->setStyleSheet("background-color: #e74c3c; color: white;");
    }
    else if (anyWaitingEnd) {
        ui->start_all_mark->setStyleSheet("background-color: #f39c12; color: white;");
        ui->stop_all_mark->setEnabled(true);
        ui->stop_all_mark->setStyleSheet("");
    }
    else {
        ui->start_all_mark->setStyleSheet("background-color: #f39c12; color: white;");
        ui->stop_all_mark->setEnabled(false);
        ui->stop_all_mark->setStyleSheet("");
    }
}

// ============================================================================
// Construction / Destruction
// ============================================================================

noise_marking_gui::noise_marking_gui(QWidget* parent)
    : QDialog(parent)
    , ui(std::make_unique<Ui::noise_marking_gui>())
    , m_noiseManager(std::make_unique<NoiseManager>(256.0))
    , m_buttonHandler(std::make_unique<lower_row_buttons>(this))
{
    ui->setupUi(this);

    // Size to 75% × 90% of screen
    if (auto* screen = QGuiApplication::primaryScreen()) {
        const QRect avail = screen->availableGeometry();
        resize(avail.width() * 3 / 4, static_cast<int>(avail.height() * 0.9));
        move(avail.center() - rect().center());
    }

    m_buttonHandler->setupConnections();

    // Prevent buttons from stealing keyboard focus (arrow keys navigate)
    for (auto* btn : findChildren<QPushButton*>())
        btn->setFocusPolicy(Qt::NoFocus);

    // --- Navigation shortcuts ---
    new QShortcut(QKeySequence(Qt::Key_Left), this, [this]() {
        m_currentStartTime = std::max(0.0, m_currentStartTime - m_skipInterval);
        handle_data_plot();
        updateAmpogramCursor();
        });

    new QShortcut(QKeySequence(Qt::Key_Right), this, [this]() {
        double maxStart = std::max(0.0, totalChunkDuration() - m_windowDuration);
        m_currentStartTime = std::min(m_currentStartTime + m_skipInterval, maxStart);
        handle_data_plot();
        updateAmpogramCursor();
        });

    // --- Chart setup ---
    const QList<QChartView*> charts = {
        ui->resp_cvp_axis, ui->ecg_ampogram_axis, ui->amp_ppg_axis,
        ui->resp_cvp_axis,
        ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3, ui->ppg_axis,
        ui->accel_or_abg_axis
    };
    for (auto* view : charts) {
        if (!view) continue;
        view->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        view->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        view->setRenderHint(QPainter::Antialiasing);
        view->setFocusPolicy(Qt::NoFocus);
        view->viewport()->installEventFilter(this);
    }

    m_skipInterval = ui->skip_interval_box->text().toDouble();
    if (m_skipInterval <= 0.0) m_skipInterval = 5.0;

    ui->scatter_line->setCurrentIndex(0);  // Line by default
    ui->scatter_line->setFocusPolicy(Qt::NoFocus);

    connect(ui->scatter_line, QOverload<int>::of(&QComboBox::currentIndexChanged),
        this, [this](int index) {
            m_plotMode = (index == 1) ? PlotMode::Scatter : PlotMode::Line;
            handle_data_plot();
        });

    QList<QChartView*> defaultCharts = { ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3,
                                         ui->ppg_axis, ui->ecg_ampogram_axis, ui->amp_ppg_axis,
                                         ui->resp_cvp_axis, ui->accel_or_abg_axis };
    for (int i = 0; i < defaultCharts.size(); ++i) {
        QChartView* v = defaultCharts[i];
        if (!v) continue;
        setupChartDefaults(v);
        // Reflow the chart's title (which carries native Hz + px/s) whenever
        // the plot area settles. plotAreaChanged fires AFTER QChart has
        // updated its internal layout, so plotArea() is fresh here -- unlike
        // the widget's Resize event, which fires first and reads stale coords.
        // This is what keeps the px/s value accurate after the dialog's
        // initial layout pass and on subsequent window resizes.
        connect(v->chart(), &QChart::plotAreaChanged, v,
            [v, this](const QRectF& /*pa*/) {
                const QString sigName = v->property("signalName").toString();
                if (sigName.isEmpty()) return;
                const double nativeHz = v->property("nativeHz").toDouble();
                const double pxPerSec = (m_windowDuration > 0.0)
                    ? v->chart()->plotArea().width() / m_windowDuration : 0.0;
                v->chart()->setTitle(
                    formatChartTitle(sigName, nativeHz, pxPerSec));
            });
    }

    // Ampogram / 8h-overview series (ECG and PPG only; RESP and CVP are
    // now windowed signal plots, not 8h overviews, so they don't get
    // ampogram series or cursor bars).
    ecg1_ampogram_series = new QLineSeries();
    ecg1_ampogram_series->setName("ECG1");
    ecg2_ampogram_series = new QLineSeries();
    ecg2_ampogram_series->setName("ECG2");
    ecg3_ampogram_series = new QLineSeries();
    ecg3_ampogram_series->setName("ECG3");
    ppg_ampogram_series = new QLineSeries();
    ui->ecg_ampogram_axis->chart()->addSeries(ecg1_ampogram_series);
    ui->ecg_ampogram_axis->chart()->addSeries(ecg2_ampogram_series);
    ui->ecg_ampogram_axis->chart()->addSeries(ecg3_ampogram_series);
    ui->amp_ppg_axis->chart()->addSeries(ppg_ampogram_series);

    // Initially disable all stop buttons
    ui->stop_ecg1_mark->setEnabled(false);
    ui->stop_ecg2_mark->setEnabled(false);
    ui->stop_ecg3_mark->setEnabled(false);
    ui->stopNoisePPG->setEnabled(false);
    ui->stop_all_mark->setEnabled(false);

    // Cursor bars (black vertical line tracking the currently-viewed window)
    auto addCursor = [](QChartView* view, QLineSeries*& series) {
        series = new QLineSeries();
        series->setPen(QPen(Qt::black, 2));
        view->chart()->addSeries(series);
        };
    addCursor(ui->ecg_ampogram_axis, m_ecgCursorBar);
    addCursor(ui->amp_ppg_axis, m_ppgCursorBar);

    // Hypnogram chart
    auto* hypnoChart = new QChart();
    hypnoChart->legend()->hide();
    hypnoChart->setMargins(QMargins(0, 0, 0, 0));
    ui->resp_cvp_axis->setChart(hypnoChart);

    m_hypnoCursorBar = new QLineSeries();
    m_hypnoCursorBar->setPen(QPen(Qt::black, 2));
    hypnoChart->addSeries(m_hypnoCursorBar);

    m_currentMarkingType = ui->marking_type->currentText();

    // Explicit connects for the file selector
    connect(ui->browse_file_button, &QPushButton::clicked,
        this, &noise_marking_gui::handleBrowseFile);

}

noise_marking_gui::~noise_marking_gui() = default;

GenExcStruct noise_marking_gui::getMarkings() const {
    GenExcStruct result = m_genExc;
    result.filePath = m_binFilePath;
    return result;
}

QVector<GenExcStruct> noise_marking_gui::getAllMarkings() const {
    // Start with stashed markings for other files
    QMap<QString, GenExcStruct> all = m_fileMarkings;
    // Overwrite / add the current file's live markings
    GenExcStruct current = m_genExc;
    current.filePath = m_binFilePath;
    all[m_binFilePath] = current;

    QVector<GenExcStruct> result;
    for (auto it = all.cbegin(); it != all.cend(); ++it) {
        if (!it->noiseExc.isEmpty())
            result.append(it.value());
    }
    return result;
}

// ============================================================================
// File Loading
// ============================================================================

void noise_marking_gui::setFileSource(const QString& filePath) {
    loadSelectedFile(filePath);
}

void noise_marking_gui::loadSelectedFile(const QString& filePath) {
    // Stash markings for the file we're leaving so they survive a file switch.
    if (!m_binFilePath.isEmpty()) {
        m_genExc.filePath = m_binFilePath;
        m_fileMarkings[m_binFilePath] = m_genExc;
    }

    m_binFilePath = filePath;
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly)) return;

    // 512-byte header: 4 rates + 41 upsampled-sizes + 41 raw-sizes +
    //                  41 native-rate floats + 1 sleep-size.
    // Read as raw bytes; sizes are uint32 but the native-rate block is
    // float32, so we reinterpret each slice after the read.
    constexpr int kNumChannels = NUM_CHANNELS;   // pulled from the header
    constexpr int kNumFields = 4 + 3 * kNumChannels + 1;   // = 128
    static_assert(FILE_HEADER_SIZE == kNumFields * 4,
        "FILE_HEADER_SIZE must match 4 + 3*N + 1 slot count");

    uint32_t raw32[kNumFields] = {};
    file.read(reinterpret_cast<char*>(raw32), sizeof(raw32));
    file.close();

    m_ecgSR = static_cast<double>(raw32[0]);
    m_boolSR = static_cast<double>(raw32[1]);
    m_ppgSR = m_ecgSR;
    double sleepEpoch = static_cast<double>(raw32[3]);
    m_sleepSR = (sleepEpoch > 0) ? (1.0 / sleepEpoch) : 0;

    constexpr int kSizesUpBase = 4;                                   // = 4
    constexpr int kSizesRawBase = kSizesUpBase + kNumChannels;         // = 45
    constexpr int kNativeRatesBase = kSizesRawBase + kNumChannels;        // = 86
    constexpr int kSleepCountIdx = kNativeRatesBase + kNumChannels;     // = 127

    for (int i = 0; i < kNumChannels; ++i) {
        m_chanSizes[i] = raw32[kSizesUpBase + i];
        m_chanSizesRaw[i] = raw32[kSizesRawBase + i];
        // Native rates are float32 on disk; reinterpret the slot.
        std::memcpy(&m_chanNativeRates[i],
            &raw32[kNativeRatesBase + i], sizeof(float));
    }
    m_totalSleepSamples = raw32[kSleepCountIdx];

    // Raw samples carry their own (t, v) timestamps on disk, so we don't
    // need to infer per-channel raw rates from m_chanNativeRates; it's
    // read to keep the in-memory struct in sync with the on-disk format.
    // `m_chanSizesRaw[i]` counts PAIRS; the byte length of the raw block
    // is 2 * size * sizeof(double).

    // Restore markings if we've seen this file before, otherwise start fresh.
    if (m_fileMarkings.contains(filePath)) {
        m_genExc = m_fileMarkings[filePath];
        m_noiseManager = std::make_unique<NoiseManager>(m_ecgSR);
        for (int i = 0; i < m_genExc.noiseExc.size(); ++i) {
            double sr = sampleRateForSignal(m_genExc.data_type[i]);
            m_noiseManager->addSegment(
                static_cast<size_t>(m_genExc.noiseExc[i].first * sr),
                static_cast<size_t>(m_genExc.noiseExc[i].second * sr),
                m_genExc.data_type[i].toStdString(),
                m_genExc.marking_type[i].toStdString()
            );
        }
    }
    else {
        m_genExc = GenExcStruct();
        m_genExc.filePath = filePath;
        m_noiseManager = std::make_unique<NoiseManager>(m_ecgSR);
    }

    m_currentStartTime = 0.0;
    m_markAllActive = false;

    for (const QString& lbl : markableChannelLabels())
        cancelMarking(lbl);

    setWindowTitle("Marking: " + QFileInfo(filePath).fileName());

    loadChunkFromFile(0);
}

void noise_marking_gui::handleBrowseFile() {
    QString startDir;
    if (!m_binFilePath.isEmpty())
        startDir = QFileInfo(m_binFilePath).absolutePath();

    QString path = QFileDialog::getOpenFileName(
        this, "Select Binary Signal File", startDir,
        "Binary files (*.bin);;All files (*)");

    if (path.isEmpty() || path == m_binFilePath) return;

    loadSelectedFile(path);
}

bool noise_marking_gui::loadChunkFromFile(uint64_t chunkIndex) {
    QFile file(m_binFilePath);
    if (!file.open(QIODevice::ReadOnly)) return false;

    m_currentChunkIndex = chunkIndex;

    // On-disk layout for each channel slot is (upsampled block, raw block).
    // The raw block is a stream of (t, v) PAIRS, so its byte length is
    // 2 * m_chanSizesRaw[i] * sizeof(double). Walk the channels in order,
    // tracking the byte offset (in doubles, from end-of-header) to the
    // START of each slot's upsampled block. The raw block for slot i starts
    // at (upsampled offset) + (upsampled size); the next slot's upsampled
    // block starts at (raw offset) + 2 * (raw pair count).
    uint64_t chanUpOffset[NUM_CHANNELS];
    uint64_t chanRawOffset[NUM_CHANNELS];
    uint64_t running = 0;
    for (int i = 0; i < NUM_CHANNELS; ++i) {
        chanUpOffset[i] = running;
        running += m_chanSizes[i];
        chanRawOffset[i] = running;
        running += m_chanSizesRaw[i] * 2;   // 2 doubles per (t, v) pair
    }
    uint64_t sleepByteOffset = running;

    // Sampling rate of the UPSAMPLED block for a given channel index.
    // 1 Hz channels are tagged by slot number; everything else is at m_ecgSR.
    auto rateForChannel = [this](int chIdx) -> double {
        if (chIdx == CH_MARKER || chIdx == CH_TEMP || chIdx == CH_PACEMAKER)
            return m_boolSR;
        if (chIdx >= CH_EKG_OFF && chIdx <= CH_EEG3_OFF)
            return m_boolSR;
        if (chIdx == CH_OXSTATUS || chIdx == CH_SPO2 || chIdx == CH_HR)
            return m_boolSR;
        return m_ecgSR;
        };

    // Load 8-hour chunk of the upsampled block for channel `chIdx` into `dest`.
    auto loadSignal = [&](QVector<double>& dest, int chIdx) {
        double sr = rateForChannel(chIdx);
        uint64_t totalSamples = m_chanSizes[chIdx];
        uint64_t perChunk = static_cast<uint64_t>(CHUNK_DURATION_SEC * sr);
        uint64_t start = chunkIndex * perChunk;
        uint64_t count = (totalSamples > start)
            ? std::min(perChunk, totalSamples - start) : 0;
        dest.resize(static_cast<int>(count));
        file.seek(FILE_HEADER_SIZE + (chanUpOffset[chIdx] + start) * sizeof(double));
        file.read(reinterpret_cast<char*>(dest.data()), count * sizeof(double));
        };

    // Load the (t, v) pair block for a markable channel, keeping only pairs
    // whose global timestamp falls inside the current 8-hour chunk. Stores
    // chunk-local time in the x-component of each QPointF, so the plotter
    // can use it directly against the chart's chunk-local x axis.
    //
    // We stream the pair block in blocks of ~32k pairs to avoid loading an
    // entire night's raw samples at once. Since timestamps are monotonic,
    // we can bail out early once t >= chunkEnd.
    auto loadRaw = [&](QVector<QPointF>& dest, int chIdx) {
        dest.clear();
        uint64_t totalPairs = m_chanSizesRaw[chIdx];

        // Missing / sentinel: file stores a single (-1.0, -1.0) pair.
        // Keep one sentinel QPointF so isMissingSignal() still flags it.
        if (totalPairs <= 1) {
            if (totalPairs == 1) {
                double pair[2] = { -1.0, -1.0 };
                file.seek(FILE_HEADER_SIZE + chanRawOffset[chIdx] * sizeof(double));
                file.read(reinterpret_cast<char*>(pair), 2 * sizeof(double));
                dest.append(QPointF(pair[0], pair[1]));
            }
            return;
        }

        const double chunkStart = chunkIndex * CHUNK_DURATION_SEC;
        const double chunkEnd = chunkStart + CHUNK_DURATION_SEC;

        constexpr uint64_t BLOCK_PAIRS = 1 << 15;   // 32768 pairs per read
        std::vector<double> buf(BLOCK_PAIRS * 2);

        file.seek(FILE_HEADER_SIZE + chanRawOffset[chIdx] * sizeof(double));
        uint64_t pairsRead = 0;
        bool done = false;
        while (!done && pairsRead < totalPairs) {
            uint64_t thisBlock = std::min<uint64_t>(BLOCK_PAIRS, totalPairs - pairsRead);
            qint64 got = file.read(reinterpret_cast<char*>(buf.data()),
                thisBlock * 2 * sizeof(double));
            if (got <= 0) break;
            uint64_t gotPairs = static_cast<uint64_t>(got) / (2 * sizeof(double));

            for (uint64_t k = 0; k < gotPairs; ++k) {
                double t = buf[k * 2];
                double v = buf[k * 2 + 1];
                if (t < chunkStart) continue;
                if (t >= chunkEnd) { done = true; break; }
                dest.append(QPointF(t - chunkStart, v));
            }
            pairsRead += gotPairs;
        }
        };

    // Upsampled signals (the line/scatter foreground plots).
    loadSignal(m_ecg1, CH_ECG1);
    loadSignal(m_ecg2, CH_ECG2);
    loadSignal(m_ecg3, CH_ECG3);
    loadSignal(m_ppg, CH_PPG);
    loadSignal(m_accelX, CH_ACCEL_X);
    loadSignal(m_accelY, CH_ACCEL_Y);
    loadSignal(m_accelZ, CH_ACCEL_Z);
    loadSignal(m_cvp, CH_PRES);    // slot 16 -- MESA nasal pressure / Bittium CVP
    loadSignal(m_resp, CH_RESP);    // slot 34 -- respiration
    loadSignal(m_abp, CH_ABP);     // slot 35 -- Bittium arterial blood pressure

    // Raw overlays. Markable channels (ECG/PPG/ABP) get scatter overlays
    // for click-to-mark. Accel channels also get raw overlays so their
    // native ~25 Hz sample positions are visible against the upsampled line.
    loadRaw(m_ecg1Raw, CH_ECG1);
    loadRaw(m_ecg2Raw, CH_ECG2);
    loadRaw(m_ecg3Raw, CH_ECG3);
    loadRaw(m_ppgRaw, CH_PPG);
    loadRaw(m_abpRaw, CH_ABP);
    loadRaw(m_accelXRaw, CH_ACCEL_X);
    loadRaw(m_accelYRaw, CH_ACCEL_Y);
    loadRaw(m_accelZRaw, CH_ACCEL_Z);

    // Sleep stages live just past the last channel slot.
    {
        uint64_t perChunk = static_cast<uint64_t>(CHUNK_DURATION_SEC * m_sleepSR);
        uint64_t start = chunkIndex * perChunk;
        uint64_t count = (m_totalSleepSamples > start)
            ? std::min(perChunk, m_totalSleepSamples - start) : 0;
        m_sleepStages.resize(static_cast<int>(count));
        file.seek(FILE_HEADER_SIZE + (sleepByteOffset + start) * sizeof(double));
        file.read(reinterpret_cast<char*>(m_sleepStages.data()), count * sizeof(double));
    }

    file.close();

    int startHr = chunkIndex * 8;
    ui->topLabel->setText(
        QString("     Data Range: Hour %1 to Hour %2").arg(startHr).arg(startHr + 8));

    m_currentStartTime = 0;

    m_activeChannels.clear();
    auto markActive = [this](const QString& label, const QVector<double>& data) {
        bool missing = isMissingSignal(data);
        if (auto* cv = chartViewForSignalLabel(label))
            cv->setVisible(!missing);
        if (!missing)
            m_activeChannels.insert(label);
        };
    markActive("ECG1", m_ecg1);
    markActive("ECG2", m_ecg2);
    markActive("ECG3", m_ecg3);
    markActive("PPG", m_ppg);

    bool anyAccel = !isMissingSignal(m_accelX)
        || !isMissingSignal(m_accelY)
        || !isMissingSignal(m_accelZ);
    if (!anyAccel) markActive("ABP", m_abp);

    // accel_or_abg_axis: show if accel OR ABP has data
    if (ui->accel_or_abg_axis) {
        bool accelPresent = !isMissingSignal(m_accelX)
            || !isMissingSignal(m_accelY)
            || !isMissingSignal(m_accelZ);
        bool abpPresent = !isMissingSignal(m_abp);
        ui->accel_or_abg_axis->setVisible(accelPresent || abpPresent);
    }

    // resp_cvp_axis: show if respiration has data
    if (ui->resp_cvp_axis) {
        ui->resp_cvp_axis->setVisible(!isMissingSignal(m_resp));
    }

    // resp_cvp_axis: show if sleep stages OR CVP has data
    if (ui->resp_cvp_axis) {
        bool sleepPresent = !m_sleepStages.isEmpty()
            && !(m_sleepStages.size() == 1 && m_sleepStages[0] == -1.0);
        bool cvpPresent = !isMissingSignal(m_cvp);
        ui->resp_cvp_axis->setVisible(sleepPresent || cvpPresent);
    }

    updateAllChannelButtonStates();

    handle_ampogram_plot();
    handle_data_plot();
    setupHypnogram();        // may replace resp_cvp_axis contents; must run before cursor update
    updateAmpogramCursor();
    restoreMarkingMarkers();

    uint64_t ecgPerChunk = static_cast<uint64_t>(CHUNK_DURATION_SEC * m_ecgSR);
    ui->prev8hours->setEnabled(chunkIndex > 0);
    ui->next8hours->setEnabled((chunkIndex * ecgPerChunk + m_ecg1.size()) < m_chanSizes[CH_ECG1]);

    return true;
}

void noise_marking_gui::on_next8hours_clicked() { loadChunkFromFile(m_currentChunkIndex + 1); }
void noise_marking_gui::on_prev8hours_clicked() { if (m_currentChunkIndex > 0) loadChunkFromFile(m_currentChunkIndex - 1); }

// ============================================================================
// Hypnogram
// ============================================================================

void noise_marking_gui::setupHypnogram() {
    if (m_sleepSR <= 0.0) return;

    // No sleep stages in this file � the resp_cvp_axis is being
    // used for CVP (handled in handle_data_plot). Do nothing here.
    bool sleepPresent = !m_sleepStages.isEmpty()
        && !(m_sleepStages.size() == 1 && m_sleepStages[0] == -1.0);
    if (!sleepPresent) return;

    auto* chart = ui->resp_cvp_axis->chart();


    if (m_cvpCursorBar && m_cvpCursorBar->chart() == chart) {
        chart->removeSeries(m_cvpCursorBar);
    }

    for (auto* s : m_hypnoStageSeries) { chart->removeSeries(s); delete s; }
    m_hypnoStageSeries.clear();

    struct Stage { int value; QColor color; };
    const QList<Stage> stages = {
        {0, Qt::black}, {1, Qt::darkGreen}, {2, Qt::blue}, {3, Qt::cyan}, {4, Qt::red}
    };

    double dt = 1.0 / m_sleepSR;
    double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;

    if (!m_sleepStages.isEmpty()) {
        for (const auto& st : stages) {
            auto* s = new QScatterSeries();
            s->setColor(st.color);
            s->setMarkerSize(3.0);
            s->setPen(Qt::NoPen);
            s->setMarkerShape(QScatterSeries::MarkerShapeRectangle);
            for (int i = 0; i < m_sleepStages.size(); ++i) {
                if (static_cast<int>(m_sleepStages[i]) == st.value)
                    s->append(globalOffset + i * dt + dt / 2.0, st.value);
            }
            chart->addSeries(s);
            m_hypnoStageSeries.append(s);
        }
    }

    clearAxes(chart);

    // Match the labelling style used for the ampograms.
    chart->setTitle("Sleep stages");
    chart->setTitleFont(QFont("Arial", 8, QFont::Bold));
    chart->setTitleBrush(Qt::black);

    auto* xAxis = new QCategoryAxis();
    xAxis->setRange(globalOffset, globalOffset + CHUNK_DURATION_SEC);
    for (int h = 0; h <= 8; ++h)
        xAxis->append(QString::number(h), globalOffset + h * 3600.0);
    xAxis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
    xAxis->setGridLineVisible(false);
    xAxis->setLabelsFont(QFont("Arial", 6));

    auto* yAxis = new QCategoryAxis();
    for (const auto& st : stages)
        yAxis->append("", st.value + 0.4);
    yAxis->setRange(-0.5, 4.5);
    yAxis->setReverse(true);
    yAxis->setVisible(false);
    yAxis->setGridLineVisible(false);

    chart->addAxis(xAxis, Qt::AlignBottom);
    chart->addAxis(yAxis, Qt::AlignLeft);
    for (auto* s : chart->series()) {
        s->attachAxis(xAxis);
        s->attachAxis(yAxis);
    }

    if (m_hypnoCursorBar) {
        chart->removeSeries(m_hypnoCursorBar);
        chart->addSeries(m_hypnoCursorBar);
        m_hypnoCursorBar->attachAxis(xAxis);
        m_hypnoCursorBar->attachAxis(yAxis);
    }
}

void noise_marking_gui::handle_ampogram_plot(double range) {
    /*
        The ampogram displays the variability in amplitude across the range param for a given channel, which can suggest
        where there may be noise.
    */
    double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;

    auto calculate_amplitude = [range, globalOffset](
        const QVector<double>& data, double sr) {
            QList<QPointF> pts;
            if (data.isEmpty() || sr <= 0.0) return pts;
            double duration = data.size() / sr;
            for (double t = 0; t <= duration - range; t += range) {
                int s = static_cast<int>(t * sr);
                int e = static_cast<int>((t + range) * sr);
                auto [mi, ma] = std::minmax_element(data.begin() + s, data.begin() + e);
                pts.append({ globalOffset + t, *ma - *mi });
            }
            return pts;
        };

    auto create_plot = [this, globalOffset](
        QChartView* view, QLineSeries* series, const QList<QPointF>& pts,
        QLineSeries* cursor, const QColor& color, const QString& title) {
            series->replace(pts);
            series->setPen(QPen(color, 1));

            auto* chart = view->chart();
            clearAxes(chart);
            chart->legend()->hide();

            // Small colored title in the top-left of the chart (matches the
            // labelling style used by plotDisplayChart).
            if (!title.isEmpty()) {
                chart->setTitle(title);
                chart->setTitleFont(QFont("Arial", 8, QFont::Bold));
                chart->setTitleBrush(color);
            }
            else {
                chart->setTitle(QString());
            }

            auto* xAxis = new QCategoryAxis();
            xAxis->setRange(globalOffset, globalOffset + CHUNK_DURATION_SEC);
            for (int h = 0; h <= 8; ++h)
                xAxis->append(QString::number(h), globalOffset + h * 3600.0);
            xAxis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
            xAxis->setGridLineVisible(false);
            xAxis->setLabelsVisible(false);
            chart->addAxis(xAxis, Qt::AlignBottom);
            series->attachAxis(xAxis);
            if (cursor) cursor->attachAxis(xAxis);

            double yMin = 0, yMax = 1.0;
            if (!pts.isEmpty()) {
                auto [mi, ma] = std::minmax_element(pts.begin(), pts.end(),
                    [](const QPointF& a, const QPointF& b) { return a.y() < b.y(); });
                double pad = std::max(0.5, (ma->y() - mi->y()) * 0.05);
                yMin = mi->y() - pad;
                yMax = ma->y() + pad;
            }

            auto* yAxis = new QValueAxis();
            yAxis->setRange(yMin, yMax);
            yAxis->setVisible(false);
            yAxis->setGridLineVisible(false);
            chart->addAxis(yAxis, Qt::AlignLeft);
            series->attachAxis(yAxis);
            if (cursor) cursor->attachAxis(yAxis);
        };

    auto ecg1Pts = calculate_amplitude(m_ecg1, m_ecgSR);
    auto ecg2Pts = calculate_amplitude(m_ecg2, m_ecgSR);
    auto ecg3Pts = calculate_amplitude(m_ecg3, m_ecgSR);


    auto* chart = ui->ecg_ampogram_axis->chart();
    chart->removeSeries(ecg2_ampogram_series);
    chart->removeSeries(ecg3_ampogram_series);;

    // 1. Plot ECG1 (sets up axes + cursor)
    create_plot(ui->ecg_ampogram_axis, ecg1_ampogram_series,
        ecg1Pts, m_ecgCursorBar, QColor("#BF00FF"), "ECG amplitude");

    // 2. Add ECG2 and ECG3
    auto* xAxis = chart->axes(Qt::Horizontal).first();
    auto* yAxis = qobject_cast<QValueAxis*>(chart->axes(Qt::Vertical).first());

    if (!ecg2Pts.isEmpty()) {
        ecg2_ampogram_series->replace(ecg2Pts);
        ecg2_ampogram_series->setPen(QPen(QColor("#FF0000"), 1));
        chart->addSeries(ecg2_ampogram_series);
        ecg2_ampogram_series->attachAxis(xAxis);
        ecg2_ampogram_series->attachAxis(yAxis);
    }
    if (!ecg3Pts.isEmpty()) {
        ecg3_ampogram_series->replace(ecg3Pts);
        ecg3_ampogram_series->setPen(QPen(QColor("#0000FF"), 1));
        chart->addSeries(ecg3_ampogram_series);
        ecg3_ampogram_series->attachAxis(xAxis);
        ecg3_ampogram_series->attachAxis(yAxis);
    }

    // 3. Expand Y range to fit all three
    auto allPts = ecg1Pts + ecg2Pts + ecg3Pts;
    if (yAxis && !allPts.isEmpty()) {
        auto [mi, ma] = std::minmax_element(allPts.begin(), allPts.end(),
            [](const QPointF& a, const QPointF& b) { return a.y() < b.y(); });
        yAxis->setRange(mi->y(), ma->y());
    }
    chart->setTitle("Amp-O-Gram");
    chart->legend()->show();
    const auto markers = chart->legend()->markers(m_ecgCursorBar);
    for (auto* m : markers) m->setVisible(false);
    chart->legend()->setAlignment(Qt::AlignTop);
    chart->legend()->setFont(QFont("Arial", 7));


    create_plot(ui->amp_ppg_axis, ppg_ampogram_series, calculate_amplitude(m_ppg, m_ppgSR), m_ppgCursorBar, COLOR_PPG,
        "PPG Amp-O-Gram");
}

// ============================================================================
// Main signal plot
// ============================================================================

void noise_marking_gui::handle_data_plot() {
    // Clean up existing highlights BEFORE plotSignal clears chart series
    for (auto* area : m_highlights) {
        if (area->chart()) area->chart()->removeSeries(area);
        delete area->upperSeries();
        delete area->lowerSeries();
        delete area;
    }
    m_highlights.clear();

    // Pick the bottommost currently-visible main chart. Its x-axis will show
    // the time labels; every chart above it hides them. This way if ABP, PPG,
    // or any of the ECGs is the last visible row, it gets the time ruler.
    // Each column has its own bottom-most visible chart that gets the time ruler.
    QChartView* xLabelOwnerRight = nullptr;
    QChartView* xLabelOwnerLeft = nullptr;
    {
        const QList<QChartView*> rightCol = {
            ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3,
            ui->ppg_axis, ui->accel_or_abg_axis
        };
        for (auto* cv : rightCol) {
            if (cv && cv->isVisible()) xLabelOwnerRight = cv;
        }

        const QList<QChartView*> leftCol = {
            ui->resp_cvp_axis, ui->resp_cvp_axis
        };
        for (auto* cv : leftCol) {
            if (cv && cv->isVisible()) xLabelOwnerLeft = cv;
        }
    }

    // Belt-and-suspenders: wipe every managed chart's series + axes up front.
    // This prevents stale series from a previous plot mode (line vs. scatter)
    // from being left behind on a chart that the current redraw path doesn't
    // visit (e.g. a channel that just became inactive, or CVP/RESP getting
    // hidden by a toggle).
    auto wipeChart = [](QChartView* cv, const QList<QAbstractSeries*>& keep) {
        if (!cv || !cv->chart()) return;
        QChart* chart = cv->chart();
        const auto serieses = chart->series();
        for (auto* s : serieses) {
            if (keep.contains(s)) continue;
            chart->removeSeries(s);
            delete s;
        }
        const auto axes = chart->axes();
        for (auto* a : axes) { chart->removeAxis(a); delete a; }
        };
    // Preserve each channel's active start marker (we don't want to delete
    // an in-progress marking just because of a redraw).
    QList<QAbstractSeries*> keepECG1 = m_markState_ecg1.startMarkerLine
        ? QList<QAbstractSeries*>{ m_markState_ecg1.startMarkerLine } : QList<QAbstractSeries*>{};
    QList<QAbstractSeries*> keepECG2 = m_markState_ecg2.startMarkerLine
        ? QList<QAbstractSeries*>{ m_markState_ecg2.startMarkerLine } : QList<QAbstractSeries*>{};
    QList<QAbstractSeries*> keepECG3 = m_markState_ecg3.startMarkerLine
        ? QList<QAbstractSeries*>{ m_markState_ecg3.startMarkerLine } : QList<QAbstractSeries*>{};
    QList<QAbstractSeries*> keepPPG = m_markState_ppg.startMarkerLine
        ? QList<QAbstractSeries*>{ m_markState_ppg.startMarkerLine } : QList<QAbstractSeries*>{};
    wipeChart(ui->ecg_axis_1, keepECG1);
    wipeChart(ui->ecg_axis_2, keepECG2);
    wipeChart(ui->ecg_axis_3, keepECG3);
    wipeChart(ui->ppg_axis, keepPPG);
    // resp_cvp_axis is owned here for windowed CVP only when there are
    // no sleep stages. When sleep stages are present, setupHypnogram() owns
    // this chart and we must leave it alone.
    {
        bool sleepPresent = !m_sleepStages.isEmpty()
            && !(m_sleepStages.size() == 1 && m_sleepStages[0] == -1.0);
        if (!sleepPresent) wipeChart(ui->resp_cvp_axis, {});
    }

    // ----------------------------------------------------------------------
    // Plot a single markable channel.
    //
    // Layering (bottom to top):
    //   1. In Line mode: upsampled samples drawn as a colored line.
    //      In Scatter mode: nothing drawn for the upsampled data.
    //   2. Raw samples drawn as a solid black scatter, always whenever
    //      rawData is non-empty. In Scatter mode this is the ONLY visible
    //      data. In Line mode it sits on top of the colored line so the
    //      marker can see where the native-rate samples land vs. what the
    //      upsampler inserted between them.
    //
    // The Y range is computed across BOTH overlays so the raw scatter is
    // never clipped by a tighter upsampled envelope.
    //
    // Returns the (min, max) of the combined sample set for downstream use.
    // ----------------------------------------------------------------------
    auto plotSignal = [&](QChartView* view,
        const QVector<double>& data, double sr,
        const QVector<QPointF>& rawData,
        QLineSeries* marker, double markerPos,
        const QColor& color,
        const QString& signalName = QString())
        -> std::pair<double, double> {
        if (!view || !view->chart()) return { 1e9, -1e9 };
        QChart* chart = view->chart();

        for (auto* s : chart->series()) { if (s != marker) { chart->removeSeries(s); delete s; } }
        for (auto* a : chart->axes()) { chart->removeAxis(a); delete a; }

        // Small colored title in the top-left of the chart, with native
        // (un-upsampled) Hz and the current pixel/second resolution
        // appended to the right of the channel name.
        double nativeHz = sr;
        if (rawData.size() >= 2
            && rawData.last().x() > rawData.first().x()
            && !(rawData.size() == 1 && rawData[0].x() == -1.0)) {
            nativeHz = (rawData.size() - 1)
                / (rawData.last().x() - rawData.first().x());
        }
        const double pxPerSec = (m_windowDuration > 0.0)
            ? chart->plotArea().width() / m_windowDuration : 0.0;
        // Stash on the view so plotAreaChanged can refresh the title with a
        // fresh pxPerSec after the dialog finishes its initial layout.
        view->setProperty("signalName", signalName);
        view->setProperty("nativeHz", nativeHz);

        if (!signalName.isEmpty()) {
            chart->setTitle(formatChartTitle(signalName, nativeHz, pxPerSec));
            QFont titleFont("Arial", 8, QFont::Bold);
            chart->setTitleFont(titleFont);
            chart->setTitleBrush(color);
        }
        else {
            chart->setTitle(QString());
        }

        auto* xAxis = new QCategoryAxis();
        xAxis->setRange(m_currentStartTime, m_currentStartTime + m_windowDuration);
        double offset = m_currentChunkIndex * CHUNK_DURATION_SEC;
        for (int i = 0; i <= 4; ++i) {
            double val = m_currentStartTime + i * m_windowDuration / 4.0;
            double t = offset + val;
            int h = static_cast<int>(t / 3600);
            int m = static_cast<int>(fmod(t, 3600) / 60);
            double s = fmod(t, 60.0);
            xAxis->append(QString("%1:%2:%3")
                .arg(h, 2, 10, QChar('0'))
                .arg(m, 2, 10, QChar('0'))
                .arg(s, 5, 'f', 2, QChar('0')), val);
        }
        xAxis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
        xAxis->setGridLineVisible(false);
        xAxis->setLabelsFont(QFont("Arial", 7));
        if (view != xLabelOwnerRight) xAxis->setLabelsVisible(false);
        chart->addAxis(xAxis, Qt::AlignBottom);
        chart->setMargins(QMargins(0, 0, 20, 0));

        auto* yAxis = new QValueAxis();
        yAxis->setVisible(false);
        chart->addAxis(yAxis, Qt::AlignLeft);

        if (data.size() < 2 || sr <= 0.0) {
            yAxis->setRange(-1.0, 1.0);
            return { 1e9, -1e9 };
        }

        double lMin = 1e9;
        double lMax = -1e9;

        // --- (1) Upsampled foreground: line in Line mode, nothing in Scatter
        //         mode. In Scatter mode the ONLY visible markers are the raw
        //         (black) dots from step (2); we still scan the upsampled data
        //         here to keep the Y range honest (marker line + noise
        //         highlights depend on it).
        QXYSeries* plotSeries = nullptr;
        if (m_plotMode == PlotMode::Line) {
            auto* ln = new QLineSeries();
            ln->setUseOpenGL(true);
            ln->setPen(QPen(color, 1));
            chart->addSeries(ln);
            plotSeries = ln;
        }

        int startIdx = std::clamp(static_cast<int>(m_currentStartTime * sr),
            0, static_cast<int>(data.size() - 1));
        int endIdx = std::clamp(static_cast<int>((m_currentStartTime + m_windowDuration) * sr),
            0, static_cast<int>(data.size()));

        QList<QPointF> pts;
        pts.reserve(endIdx - startIdx);
        for (int i = startIdx; i < endIdx; ++i) {
            pts.append({ static_cast<double>(i) / sr, data[i] });
            if (data[i] < lMin) lMin = data[i];
            if (data[i] > lMax) lMax = data[i];
        }
        if (plotSeries) {
            plotSeries->replace(pts);
            plotSeries->attachAxis(xAxis);
            plotSeries->attachAxis(yAxis);
        }

        // --- (2) Raw scatter overlay: drawn LAST so the black dots sit on
        //         top of the purple line. Each raw sample carries its own
        //         chunk-local timestamp (from the on-disk (t, v) pair block),
        //         so irregular spacing, dropouts, and jitter are preserved
        //         exactly. Skipped when the raw block is missing or just a
        //         sentinel.
        auto isRawUsable = [](const QVector<QPointF>& v) {
            return v.size() >= 2 && !(v.size() == 1 && v[0].x() == -1.0);
            };
        if (isRawUsable(rawData)) {
            auto* rawScatter = new QScatterSeries();
            rawScatter->setColor(COLOR_RAW_SCATTER);
            rawScatter->setBorderColor(Qt::transparent);
            rawScatter->setMarkerSize(3.0);
            rawScatter->setMarkerShape(QScatterSeries::MarkerShapeCircle);
            rawScatter->setUseOpenGL(true);
            chart->addSeries(rawScatter);

            const double viewStart = m_currentStartTime;
            const double viewEnd = m_currentStartTime + m_windowDuration;

            QList<QPointF> rawPts;
            rawPts.reserve(std::min<int>(rawData.size(), 4096));
            // Raw pairs are time-sorted in the on-disk block, so we can
            // break out early once we pass the right edge of the window.
            for (const QPointF& p : rawData) {
                if (p.x() < viewStart) continue;
                if (p.x() > viewEnd)   break;
                rawPts.append(p);
                if (p.y() < lMin) lMin = p.y();
                if (p.y() > lMax) lMax = p.y();
            }
            rawScatter->replace(rawPts);
            rawScatter->attachAxis(xAxis);
            rawScatter->attachAxis(yAxis);
        }

        // Defensive: if no samples were found in range, fall back to a
        // neutral range rather than leaving (1e9, -1e9) in the axis.
        if (lMin > lMax) { lMin = -1.0; lMax = 1.0; }
        // Pad the axis proportionally to the data range (5% each side) so
        // low-amplitude signals (MESA ECG/PPG in mV) still fill the chart.
        // Fall back to a constant pad when the window is flat (e.g. a
        // missing-channel sentinel holding a single value), which would
        // otherwise collapse the axis to zero height.
        double span = lMax - lMin;
        double pad = (span > 1e-9) ? 0.05 * span : 0.5;
        yAxis->setRange(lMin - pad, lMax + pad);

        if (marker && marker->chart() == chart) {
            marker->replace({ {markerPos, yAxis->min()}, {markerPos, yAxis->max()} });
            marker->attachAxis(xAxis);
            marker->attachAxis(yAxis);
        }
        return { lMin, lMax };
        };

    // Plot a single markable channel if it's active. Uses the ChannelRefs
    // lookup to pull both the upsampled and the raw (scatter overlay) data
    // from one source of truth.
    auto maybePlot = [&](const QString& label) {
        if (!isChannelActive(label)) return;
        ChannelRefs r = channelRefs(label);
        if (!r.chartView || !r.data || !r.state) return;

        double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
        double localMarkerPos = r.state->globalStartTime - globalOffset;

        static const QVector<QPointF> kEmptyRaw;
        const QVector<QPointF>& rawData = r.dataRaw ? *r.dataRaw : kEmptyRaw;
        double sr = (r.sampleRate ? *r.sampleRate : m_ecgSR);

        plotSignal(r.chartView, *r.data, sr, rawData,
            r.state->startMarkerLine, localMarkerPos, r.color, label);
        };

    maybePlot("ECG1");
    maybePlot("ECG2");
    maybePlot("ECG3");
    maybePlot("PPG");

    // ------------------------------------------------------------------
    // Display-only charts (no marking): helper to plot 1..N signals on
    // a chart with the shared time axis and auto-ranged Y axis.
    // Each series may optionally carry a `rawData` pointer of (t, v)
    // pairs in chunk-local seconds; those are drawn as a same-color
    // scatter overlay on top of the upsampled line, matching how the
    // markable channels surface their native-rate samples.
    // ------------------------------------------------------------------
    struct DisplaySeries {
        const QVector<double>* data;
        QColor                 color;
        const QVector<QPointF>* rawData = nullptr;   ///< optional raw (t, v) overlay
    };

    auto plotDisplayChart = [&](QChartView* view, const QString& title,
        const QColor& titleColor,
        const QList<DisplaySeries>& serieses,
        bool forceLine = false) {
            if (!view || !view->chart()) return;
            QChart* chart = view->chart();

            for (auto* s : chart->series()) { chart->removeSeries(s); delete s; }
            for (auto* a : chart->axes()) { chart->removeAxis(a); delete a; }

            chart->setTitle(title);
            chart->setTitleFont(QFont("Arial", 8, QFont::Bold));

            auto* xAxis = new QCategoryAxis();
            xAxis->setRange(m_currentStartTime, m_currentStartTime + m_windowDuration);
            double offset = m_currentChunkIndex * CHUNK_DURATION_SEC;
            for (int i = 0; i <= 4; ++i) {
                double val = m_currentStartTime + i * m_windowDuration / 4.0;
                double t = offset + val;
                int h = static_cast<int>(t / 3600);
                int mn = static_cast<int>(fmod(t, 3600) / 60);
                double s = fmod(t, 60.0);
                xAxis->append(QString("%1:%2:%3")
                    .arg(h, 2, 10, QChar('0'))
                    .arg(mn, 2, 10, QChar('0'))
                    .arg(s, 5, 'f', 2, QChar('0')), val);
            }
            xAxis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
            xAxis->setGridLineVisible(false);
            xAxis->setLabelsFont(QFont("Arial", 7));
            xAxis->setLabelsVisible(view == xLabelOwnerRight || view == xLabelOwnerLeft);
            chart->addAxis(xAxis, Qt::AlignBottom);
            chart->setMargins(QMargins(0, 0, 20, 0));

            auto* yAxis = new QValueAxis();
            yAxis->setVisible(false);
            chart->addAxis(yAxis, Qt::AlignLeft);

            const bool useScatter = (m_plotMode == PlotMode::Scatter) && !forceLine;

            // Same skip rule as the markable path: a 1-element raw vector
            // holding (-1, -1) is the missing-channel sentinel.
            auto isRawUsable = [](const QVector<QPointF>& v) {
                return v.size() >= 2 && !(v.size() == 1 && v[0].x() == -1.0);
                };

            double gMin = 1e9, gMax = -1e9;
            for (const auto& d : serieses) {
                if (!d.data || isMissingSignal(*d.data)) continue;

                // (1) Upsampled foreground (line in Line mode, scatter in
                //     Scatter mode -- same as before).
                QXYSeries* plotSeries = nullptr;
                if (useScatter) {
                    auto* sc = new QScatterSeries();
                    sc->setColor(d.color);
                    sc->setBorderColor(Qt::transparent);   // no border: crisp + faster
                    sc->setMarkerSize(2.0);
                    sc->setMarkerShape(QScatterSeries::MarkerShapeCircle);
                    sc->setUseOpenGL(true);
                    chart->addSeries(sc);
                    plotSeries = sc;
                }
                else {
                    auto* ln = new QLineSeries();
                    ln->setUseOpenGL(true);
                    ln->setPen(QPen(d.color, 1));
                    chart->addSeries(ln);
                    plotSeries = ln;
                }

                int startIdx = std::clamp(static_cast<int>(m_currentStartTime * m_ecgSR),
                    0, static_cast<int>(d.data->size() - 1));
                int endIdx = std::clamp(static_cast<int>((m_currentStartTime + m_windowDuration) * m_ecgSR),
                    0, static_cast<int>(d.data->size()));

                QList<QPointF> pts;
                pts.reserve(endIdx - startIdx);
                for (int i = startIdx; i < endIdx; ++i) {
                    double v = (*d.data)[i];
                    pts.append({ static_cast<double>(i) / m_ecgSR, v });
                    if (v < gMin) gMin = v;
                    if (v > gMax) gMax = v;
                }
                plotSeries->replace(pts);
                plotSeries->attachAxis(xAxis);
                plotSeries->attachAxis(yAxis);

                // (2) Optional raw scatter overlay -- drawn AFTER the line so
                //     the dots sit on top. Same color as the upsampled line so
                //     channel identity is preserved when multiple series share
                //     a chart (X/Y/Z accel). Slightly larger than the markable
                //     channels' overlay (3.5 vs 3.0) to read clearly through
                //     the colored line underneath.
                if (d.rawData && isRawUsable(*d.rawData)) {
                    auto* rawScatter = new QScatterSeries();
                    rawScatter->setColor(d.color);
                    rawScatter->setBorderColor(Qt::transparent);
                    rawScatter->setMarkerSize(3.5);
                    rawScatter->setMarkerShape(QScatterSeries::MarkerShapeCircle);
                    rawScatter->setUseOpenGL(true);
                    chart->addSeries(rawScatter);

                    const double viewStart = m_currentStartTime;
                    const double viewEnd = m_currentStartTime + m_windowDuration;

                    QList<QPointF> rawPts;
                    rawPts.reserve(std::min<int>(d.rawData->size(), 4096));
                    for (const QPointF& p : *d.rawData) {
                        if (p.x() < viewStart) continue;
                        if (p.x() > viewEnd)   break;
                        rawPts.append(p);
                        if (p.y() < gMin) gMin = p.y();
                        if (p.y() > gMax) gMax = p.y();
                    }
                    rawScatter->replace(rawPts);
                    rawScatter->attachAxis(xAxis);
                    rawScatter->attachAxis(yAxis);
                }
            }

            // Same proportional padding rule as the markable path: 5% of the
            // data range, with a constant-pad fallback for flat windows.
            if (gMin < gMax) {
                double span = gMax - gMin;
                double pad = (span > 1e-9) ? 0.05 * span : 0.5;
                yAxis->setRange(gMin - pad, gMax + pad);
            }
            else {
                yAxis->setRange(-1.0, 1.0);
            }
        };

    // --- accel_or_abg_axis: prefer accel; if absent, fall through to ABP ---
    if (ui->accel_or_abg_axis) {
        bool anyAccel = !isMissingSignal(m_accelX)
            || !isMissingSignal(m_accelY)
            || !isMissingSignal(m_accelZ);
        bool hasAbp = !isMissingSignal(m_abp);

        if (anyAccel) {
            plotDisplayChart(ui->accel_or_abg_axis, "ACCEL", COLOR_ACCEL_X, {
                { &m_accelX, COLOR_ACCEL_X, &m_accelXRaw },
                { &m_accelY, COLOR_ACCEL_Y, &m_accelYRaw },
                { &m_accelZ, COLOR_ACCEL_Z, &m_accelZRaw },
                });
        }
        else if (hasAbp) {
            // Route ABP through the markable path so it gets the start marker,
            // noise-highlight overlay, and click-to-mark behavior like ECG/PPG.
            // The new maybePlot reads ABP's upsampled + raw overlays and
            // chart view from channelRefs("ABP") -- all uniform with ECG/PPG.
            maybePlot("ABP");
        }
    }

    // --- resp_cvp_axis: windowed RESP and CVP together (CVP only if present) ---
    bool sleepPresent = !m_sleepStages.isEmpty()
        && !(m_sleepStages.size() == 1 && m_sleepStages[0] == -1.0);
    if (!sleepPresent && ui->resp_cvp_axis &&
        (!isMissingSignal(m_resp) || !isMissingSignal(m_cvp))) {
        plotDisplayChart(ui->resp_cvp_axis, "RESP / CVP", COLOR_RESP, {
            { &m_resp, COLOR_RESP },
            { &m_cvp,  COLOR_CVP  },
            });
    }
    
    updateNoiseHighlights();
}

// ============================================================================
// Marking
// ============================================================================

void noise_marking_gui::finalizeMarking(QChartView* cv, double endX,
    const QString& signalLabel) {
    ChannelMarkingState& state = markStateFor(signalLabel);
    double sr = sampleRateForSignal(signalLabel);

    double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
    double globalEnd = endX + globalOffset;
    double globalStart = state.globalStartTime;

    double snappedS = std::round(std::min(globalStart, globalEnd) * sr) / sr;
    double snappedE = std::round(std::max(globalStart, globalEnd) * sr) / sr;

    m_noiseManager->addSegment(
        static_cast<size_t>(snappedS * sr),
        static_cast<size_t>(snappedE * sr),
        signalLabel.toStdString(),
        m_currentMarkingType.toStdString()
    );

    m_genExc.noiseExc.append({ snappedS, snappedE });
    m_genExc.data_type.append(signalLabel);
    m_genExc.marking_type.append(m_currentMarkingType);

    // Reset state
    state.phase = MarkPhase::Idle;
    clearStartMarker(state);

    // Reset the buttons for this specific channel
    if (QPushButton* startBtn = startButtonForSignal(signalLabel))
        startBtn->setStyleSheet("");
    if (QPushButton* stopBtn = stopButtonForSignal(signalLabel)) {
        stopBtn->setStyleSheet("");
        stopBtn->setEnabled(false);
    }

    // If mark-all mode is active, check if all channels are now idle
    if (m_markAllActive) {
        bool allIdle = true;
        for (const QString& lbl : markableChannelLabels()) {
            if (isChannelActive(lbl) && markStateFor(lbl).phase != MarkPhase::Idle) {
                allIdle = false;
                break;
            }
        }
        if (allIdle) {
            m_markAllActive = false;
            ui->start_all_mark->setStyleSheet("");
            ui->stop_all_mark->setStyleSheet("");
            ui->stop_all_mark->setEnabled(false);
        }
    }

    updateNoiseHighlights();
}

void noise_marking_gui::cancelMarking(const QString& signalLabel) {
    ChannelMarkingState& state = markStateFor(signalLabel);
    state.phase = MarkPhase::Idle;
    clearStartMarker(state);

    if (QPushButton* startBtn = startButtonForSignal(signalLabel))
        startBtn->setStyleSheet("");
    if (QPushButton* stopBtn = stopButtonForSignal(signalLabel)) {
        stopBtn->setStyleSheet("");
        stopBtn->setEnabled(false);
    }
}

// ============================================================================
// Event filter
// ============================================================================

bool noise_marking_gui::eventFilter(QObject* watched, QEvent* event) {
    auto* viewport = qobject_cast<QWidget*>(watched);
    if (!viewport) return QDialog::eventFilter(watched, event);

    auto* cv = qobject_cast<QChartView*>(viewport->parent());
    if (!cv || !cv->chart()) return QDialog::eventFilter(watched, event);

    // --- Mouse press ---
    if (event->type() == QEvent::MouseButtonPress) {
        auto* me = static_cast<QMouseEvent*>(event);
        if (me->button() != Qt::LeftButton)
            return QDialog::eventFilter(watched, event);

        double clickedX = cv->chart()->mapToValue(me->pos()).x();
        double chunkDur = totalChunkDuration();
        clickedX = std::clamp(clickedX, 0.0, chunkDur);

        // --- Signal chart click ---
        QString label = signalLabelForChartView(cv);
        if (!label.isEmpty()) {
            if (!isChannelActive(label))
                return QDialog::eventFilter(watched, event);

            // Convert click to local chunk time
            double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;

            // --- Mark-all mode: propagate click to every active channel ---
            if (m_markAllActive) {
                // Categorise the phases of all active channels
                bool anyWaitingForFirstClick = false; // WaitingForStart only
                bool anyWaitingForEnd = false;        // WaitingForEnd (start placed)
                bool anyWaitingStop = false;          // WaitingForStop
                for (const QString& lbl : markableChannelLabels()) {
                    if (!isChannelActive(lbl)) continue;
                    MarkPhase p = markStateFor(lbl).phase;
                    if (p == MarkPhase::WaitingForStart)
                        anyWaitingForFirstClick = true;
                    if (p == MarkPhase::WaitingForEnd)
                        anyWaitingForEnd = true;
                    if (p == MarkPhase::WaitingForStop)
                        anyWaitingStop = true;
                }

                // 1) Finalize any channels that are WaitingForStop first.
                //    This handles the case where some channels were individually
                //    stopped while others are still in WaitingForEnd.
                if (anyWaitingStop) {
                    for (const QString& lbl : markableChannelLabels()) {
                        if (!isChannelActive(lbl)) continue;
                        ChannelMarkingState& st = markStateFor(lbl);
                        if (st.phase != MarkPhase::WaitingForStop) continue;
                        QChartView* targetCv = chartViewForSignalLabel(lbl);
                        finalizeMarking(targetCv, clickedX, lbl);
                    }
                    // Don't touch channels that are still in other phases;
                    // they keep their existing start markers.
                    return true;
                }

                // 2) Place (or move) the start marker on channels that
                //    haven't had their start set yet, or are having it moved.
                if (anyWaitingForFirstClick || anyWaitingForEnd) {
                    for (const QString& lbl : markableChannelLabels()) {
                        if (!isChannelActive(lbl)) continue;
                        ChannelMarkingState& st = markStateFor(lbl);
                        if (st.phase != MarkPhase::WaitingForStart &&
                            st.phase != MarkPhase::WaitingForEnd)
                            continue;
                        st.globalStartTime = clickedX + globalOffset;
                        QChartView* targetCv = chartViewForSignalLabel(lbl);
                        QPushButton* stopBtn = stopButtonForSignal(lbl);
                        showStartMarker(targetCv, clickedX, st,
                            colorForSignal(lbl), stopBtn);
                        st.phase = MarkPhase::WaitingForEnd;
                    }
                    ui->stop_all_mark->setEnabled(true);
                    return true;
                }
            }

            // --- Single-channel mode (original behavior) ---
            ChannelMarkingState& state = markStateFor(label);

            switch (state.phase) {
            case MarkPhase::WaitingForStart:
            case MarkPhase::WaitingForEnd: {
                state.globalStartTime = clickedX + globalOffset;
                QPushButton* stopBtn = stopButtonForSignal(label);
                showStartMarker(cv, clickedX, state, colorForSignal(label), stopBtn);
                state.phase = MarkPhase::WaitingForEnd;
                return true;
            }
            case MarkPhase::WaitingForStop:
                finalizeMarking(cv, clickedX, label);
                return true;

            case MarkPhase::Idle:
                m_isDragging = true;
                m_dragStartPos = me->pos();
                m_dragSignalLabel = label;
                {
                    state.globalStartTime = clickedX + globalOffset;
                }
                if (!m_draggedViewport) {
                    m_draggedViewport = viewport;
                    m_draggedViewport->grabMouse();
                }
                return true;
            }
        }

        // --- Overview click (navigate) ---
        // Ampograms (amp_ecg1, amp_ppg) are always full-8h overview charts.
        // resp_cvp_axis is a full-8h hypnogram when sleep stages are
        // present. RESP and (CVP-mode) resp_cvp_axis are now windowed
        // plots that move with the marking window, so they're NOT navigable
        // via click -- use the arrow keys or the ECG/PPG ampograms instead.
        bool sleepPresent = !m_sleepStages.isEmpty()
            && !(m_sleepStages.size() == 1 && m_sleepStages[0] == -1.0);
        bool isNavChart =
            (cv == ui->ecg_ampogram_axis)
            || (cv == ui->amp_ppg_axis)
            || (cv == ui->resp_cvp_axis && sleepPresent);
        if (isNavChart) {
            double globalClickX = cv->chart()->mapToValue(me->pos()).x();
            double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
            double localTarget = globalClickX - globalOffset;
            double maxStart = std::max(0.0, chunkDur - m_windowDuration);
            m_currentStartTime = std::clamp(localTarget - m_windowDuration / 2.0,
                0.0, maxStart);
            handle_data_plot();
            updateAmpogramCursor();
            return true;
        }
    }

    // --- Mouse move (drag) ---
    if (event->type() == QEvent::MouseMove && m_isDragging)
        return true;

    // --- Mouse release (end drag) ---
    if (event->type() == QEvent::MouseButtonRelease) {
        auto* me = static_cast<QMouseEvent*>(event);
        if (me->button() == Qt::LeftButton && m_isDragging) {
            if (m_draggedViewport) {
                m_draggedViewport->releaseMouse();
                m_draggedViewport = nullptr;
            }
            m_isDragging = false;

            QString label = signalLabelForChartView(cv);
            if (!label.isEmpty() && label == m_dragSignalLabel) {
                double endX = cv->chart()->mapToValue(me->pos()).x();
                endX = std::clamp(endX, m_currentStartTime,
                    m_currentStartTime + m_windowDuration);
                ChannelMarkingState& state = markStateFor(label);
                double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
                double localStart = state.globalStartTime - globalOffset;
                if (std::abs(endX - localStart) > 0.1)
                    finalizeMarking(cv, endX, label);
            }
            m_dragSignalLabel.clear();
            return true;
        }
    }

    return QDialog::eventFilter(watched, event);
}

// ============================================================================
// Cursor & highlight updates
// ============================================================================

void noise_marking_gui::updateAmpogramCursor() {
    auto draw = [this](QChartView* view, QLineSeries* cursor) {
        if (!view || !cursor) return;
        auto axes = view->chart()->axes(Qt::Vertical);
        if (axes.isEmpty()) return;
        auto* yAxis = qobject_cast<QValueAxis*>(axes.first());
        if (!yAxis) return;

        double x = m_currentChunkIndex * CHUNK_DURATION_SEC
            + m_currentStartTime + m_windowDuration / 2.0;
        cursor->replace({ {x, yAxis->min()}, {x, yAxis->max()} });
        };

    draw(ui->ecg_ampogram_axis, m_ecgCursorBar);
    draw(ui->amp_ppg_axis, m_ppgCursorBar);

    // RESP and CVP are now windowed plots (not 8-hour overviews), so they
    // don't need a cursor bar -- the chart's visible extent IS the window.

    // resp_cvp_axis still needs a cursor when showing the 8-hour hypnogram.
    bool sleepPresent = !m_sleepStages.isEmpty()
        && !(m_sleepStages.size() == 1 && m_sleepStages[0] == -1.0);
    if (sleepPresent)
        draw(ui->resp_cvp_axis, m_hypnoCursorBar);
}

void noise_marking_gui::updateNoiseHighlights() {
    for (auto* area : m_highlights) {
        if (area->chart()) area->chart()->removeSeries(area);
        delete area->upperSeries();
        delete area->lowerSeries();
        delete area;
    }
    m_highlights.clear();

    struct ChartAxes {
        QChart* chart = nullptr;
        QAbstractAxis* xAxis = nullptr;
        QValueAxis* yAxis = nullptr;
    };
    QMap<QString, ChartAxes> axesMap;

    for (const QString& lbl : markableChannelLabels()) {
        if (!isChannelActive(lbl)) continue;
        auto* cv = chartViewForSignalLabel(lbl);
        if (!cv) continue;

        ChartAxes ca;
        ca.chart = cv->chart();
        auto hAxes = ca.chart->axes(Qt::Horizontal);
        auto vAxes = ca.chart->axes(Qt::Vertical);
        ca.xAxis = hAxes.isEmpty() ? nullptr : hAxes.first();
        ca.yAxis = vAxes.isEmpty() ? nullptr : qobject_cast<QValueAxis*>(vAxes.first());
        axesMap[lbl] = ca;
    }

    double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
    double viewStart = m_currentStartTime;
    double viewEnd = viewStart + m_windowDuration;

    for (const auto& seg : m_noiseManager->getSegments()) {
        QString segLabel = QString::fromStdString(seg.label);
        if (!axesMap.contains(segLabel)) continue;

        double sr = sampleRateForSignal(segLabel);
        double segStart = seg.startSample / sr - globalOffset;
        double segEnd = seg.endSample / sr - globalOffset;
        if (segEnd < viewStart || segStart > viewEnd) continue;

        double ds = std::max(segStart, viewStart);
        double de = std::min(segEnd, viewEnd);
        QColor color = MARKING_COLORS.value(
            QString::fromStdString(seg.marking_type), QColor(0, 0, 0, 100));

        const ChartAxes& ca = axesMap[segLabel];
        if (!ca.chart || !ca.xAxis || !ca.yAxis) continue;

        auto* upper = new QLineSeries();
        auto* lower = new QLineSeries();
        upper->append({ {ds, ca.yAxis->max()}, {de, ca.yAxis->max()} });
        lower->append({ {ds, ca.yAxis->min()}, {de, ca.yAxis->min()} });

        auto* area = new QAreaSeries(upper, lower);
        area->setBrush(color);
        area->setPen(Qt::NoPen);
        ca.chart->addSeries(area);
        area->attachAxis(ca.xAxis);
        area->attachAxis(ca.yAxis);
        m_highlights.append(area);
    }
}

// ============================================================================
// Marker helpers
// ============================================================================

void noise_marking_gui::clearStartMarker(ChannelMarkingState& state) {
    if (state.startMarkerLine && state.startMarkerLine->chart()) {
        state.startMarkerLine->chart()->removeSeries(state.startMarkerLine);
        delete state.startMarkerLine;
        state.startMarkerLine = nullptr;
    }
}

void noise_marking_gui::showStartMarker(QChartView* cv, double xValue,
    ChannelMarkingState& state, const QColor& color, QPushButton* stopBtn) {
    clearStartMarker(state);

    state.startMarkerLine = new QLineSeries();
    state.startMarkerLine->setPen(QPen(color, 2, Qt::DashLine));

    auto* yAxis = qobject_cast<QValueAxis*>(cv->chart()->axes(Qt::Vertical).first());
    state.startMarkerLine->append(xValue, yAxis->min());
    state.startMarkerLine->append(xValue, yAxis->max());

    cv->chart()->addSeries(state.startMarkerLine);
    state.startMarkerLine->attachAxis(cv->chart()->axes(Qt::Horizontal).first());
    state.startMarkerLine->attachAxis(yAxis);

    if (stopBtn) stopBtn->setEnabled(true);
}

QString noise_marking_gui::formatTimeLabel(double seconds) {
    int hours = static_cast<int>(seconds) / 3600;
    int minutes = (static_cast<int>(seconds) % 3600) / 60;
    double secs = fmod(seconds, 60.0);
    return QString("%1:%2:%3")
        .arg(hours, 2, 10, QChar('0'))
        .arg(minutes, 2, 10, QChar('0'))
        .arg(secs, 5, 'f', 2, QChar('0'));
}

// ============================================================================
// Restore start markers after chunk change
// ============================================================================

void noise_marking_gui::restoreMarkingMarkers() {
    double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
    double chunkEnd = globalOffset + CHUNK_DURATION_SEC;

    auto restore = [&](const QString& label, ChannelMarkingState& state) {
        if (state.phase != MarkPhase::WaitingForEnd &&
            state.phase != MarkPhase::WaitingForStop)
            return;

        state.startMarkerLine = nullptr;

        QChartView* cv = chartViewForSignalLabel(label);
        if (!cv || !isChannelActive(label)) return;

        if (state.globalStartTime >= globalOffset && state.globalStartTime <= chunkEnd) {
            double localX = state.globalStartTime - globalOffset;
            QPushButton* stopBtn = stopButtonForSignal(label);
            showStartMarker(cv, localX, state, colorForSignal(label), stopBtn);
        }
        };

    restore("ECG1", m_markState_ecg1);
    restore("ECG2", m_markState_ecg2);
    restore("ECG3", m_markState_ecg3);
    restore("PPG", m_markState_ppg);
    restore("ABP", m_markState_abp);
}

// ============================================================================
// Public marking API (called by lower_row_buttons)
// ============================================================================

void noise_marking_gui::beginMarking(const QString& signalLabel) {
    ChannelMarkingState& state = markStateFor(signalLabel);

    if (state.phase != MarkPhase::Idle) {
        cancelMarking(signalLabel);
        return;
    }

    m_currentMarkingType = ui->marking_type->currentText();
    state.phase = MarkPhase::WaitingForStart;

    if (QPushButton* startBtn = startButtonForSignal(signalLabel))
        startBtn->setStyleSheet("background-color: #f39c12; color: white;");
}

void noise_marking_gui::beginMarkingAll() {
    m_markAllActive = true;
    m_currentMarkingType = ui->marking_type->currentText();

    for (const QString& label : markableChannelLabels()) {
        if (!isChannelActive(label)) continue;
        ChannelMarkingState& state = markStateFor(label);
        if (state.phase != MarkPhase::Idle) cancelMarking(label);
        state.phase = MarkPhase::WaitingForStart;
        if (QPushButton* startBtn = startButtonForSignal(label))
            startBtn->setStyleSheet("background-color: #f39c12; color: white;");
    }
    ui->start_all_mark->setStyleSheet("background-color: #f39c12; color: white;");
}

void noise_marking_gui::beginStopPhaseAll() {
    for (const QString& label : markableChannelLabels()) {
        if (!isChannelActive(label)) continue;
        ChannelMarkingState& state = markStateFor(label);
        if (state.phase != MarkPhase::WaitingForEnd) continue;

        state.phase = MarkPhase::WaitingForStop;
        if (QPushButton* stopBtn = stopButtonForSignal(label))
            stopBtn->setStyleSheet("background-color: #e74c3c; color: white;");
        if (QPushButton* startBtn = startButtonForSignal(label))
            startBtn->setStyleSheet("");
    }
    ui->start_all_mark->setStyleSheet("");
    ui->stop_all_mark->setStyleSheet("background-color: #e74c3c; color: white;");
}

void noise_marking_gui::beginStopPhase(const QString& signalLabel) {
    ChannelMarkingState& state = markStateFor(signalLabel);
    if (state.phase != MarkPhase::WaitingForEnd) return;

    state.phase = MarkPhase::WaitingForStop;

    if (QPushButton* stopBtn = stopButtonForSignal(signalLabel))
        stopBtn->setStyleSheet("background-color: #e74c3c; color: white;");
    if (QPushButton* startBtn = startButtonForSignal(signalLabel))
        startBtn->setStyleSheet("");
}

// ============================================================================
// Slot handlers
// ============================================================================

void noise_marking_gui::on_skip_interval_box_editingFinished() {
    m_skipInterval = ui->skip_interval_box->text().toDouble();
    ui->skip_interval_box->setText(QString::number(m_skipInterval, 'f', 1));
    ui->skip_interval_box->clearFocus();
}

void noise_marking_gui::on_skip_interval_box_returnPressed() {
    on_skip_interval_box_editingFinished();
}

void noise_marking_gui::on_marking_type_currentTextChanged(const QString& text) {
    m_currentMarkingType = text;
}