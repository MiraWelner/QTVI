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
 *              a colored line underneath the raw black scatter.
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
 *
 * @author  Mira Welner
 * @email   MEW386@pitt.edu
 * @date    2026-03-22
 */

#include "annealing_to_bin//anneal_handler.hpp"
#include "peak_finding//run_find_r_peaks.hpp"
#include "gui_handler.hpp"
#include "post_process.hpp"

#include <QFutureWatcher>
#include <QtConcurrent>
#include <QCheckBox>
#include <QDoubleSpinBox>
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
#include <QSignalBlocker>
#include <QButtonGroup>
#include <QRadioButton>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QtCharts/QLegendMarker>
#include <algorithm>
#include <cstring>
#include <iostream>

 // ============================================================================
 // Constants
 // ============================================================================

static const QColor COLOR_ECG1 = QColor("#FF0000");
static const QColor COLOR_ECG2 = QColor("#0000FF");
static const QColor COLOR_ECG3 = QColor("#00AA00");
static const QColor COLOR_PPG = QColor("#BF00FF");
static const QColor COLOR_ABP = QColor("#BF00FF");
static const QColor COLOR_ACCEL_X = QColor("#F39C12");
static const QColor COLOR_ACCEL_Y = QColor("#27AE60");
static const QColor COLOR_ACCEL_Z = QColor("#8E44AD");
static const QColor COLOR_RESP = QColor("#16A085");
static const QColor COLOR_CVP = QColor("#2980B9");
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

    // Wipe every series and axis on a chart, optionally preserving a list of
    // series the caller wants to keep (e.g. an in-progress start marker).
    void wipeChartContent(QChart* chart,
        const QList<QAbstractSeries*>& keep = {}) {
        if (!chart) return;
        const auto serieses = chart->series();
        for (auto* s : serieses) {
            if (keep.contains(s)) continue;
            chart->removeSeries(s);
            delete s;
        }
        const auto axes = chart->axes();
        for (auto* a : axes) { chart->removeAxis(a); delete a; }
    }

    QString formatChartTitle(const QString& signalName,
        double nativeHz, double pxPerSample) {
        return QString("%1  -- Original Frequency: %2 Hz -- Pixel Resolution: %3 px/sample")
            .arg(signalName)
            .arg(nativeHz, 0, 'f', 1)
            .arg(pxPerSample, 0, 'f', 3);
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

    // Raw (t, v) overlay is usable only when it has actual data -- a 1-element
    // vector holding (-1, -1) is the missing-channel sentinel.
    bool isRawUsable(const QVector<QPointF>& v) {
        return v.size() >= 2 && !(v.size() == 1 && v[0].x() == -1.0);
    }

    QString formatHMS(double seconds) {
        int h = static_cast<int>(seconds) / 3600;
        int m = (static_cast<int>(seconds) % 3600) / 60;
        double s = std::fmod(seconds, 60.0);
        return QString("%1:%2:%3")
            .arg(h, 2, 10, QChar('0'))
            .arg(m, 2, 10, QChar('0'))
            .arg(s, 5, 'f', 2, QChar('0'));
    }

    // Build the windowed HH:MM:SS x-axis used by both markable and display
    // charts. Start/duration are chunk-local; the labels include the chunk's
    // global time offset so they read as wall-clock time within the recording.
    QCategoryAxis* makeWindowedTimeAxis(double startLocal, double duration,
        double globalOffset, bool labelsVisible) {
        auto* xAxis = new QCategoryAxis();
        xAxis->setRange(startLocal, startLocal + duration);
        for (int i = 0; i <= 4; ++i) {
            double localVal = startLocal + i * duration / 4.0;
            double globalSec = globalOffset + localVal;
            xAxis->append(formatHMS(globalSec), localVal);
        }
        xAxis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
        xAxis->setGridLineVisible(false);
        xAxis->setLabelsFont(QFont("Arial", 7));
        xAxis->setLabelsVisible(labelsVisible);
        return xAxis;
    }

    // Apply the standard 5%-of-span Y range with a constant fallback for flat
    // or empty windows. Used by both signal-plot paths.
    void setPaddedYRange(QValueAxis* yAxis, double yMin, double yMax) {
        if (yMin > yMax) { yMin = -1.0; yMax = 1.0; }
        double span = yMax - yMin;
        double pad = (span > 1e-9) ? 0.05 * span : 0.5;
        yAxis->setRange(yMin - pad, yMax + pad);
    }

    // True when the .bin's sleep-stage block holds real data (not the -1
    // sentinel and not empty).
    bool sleepDataPresent(const QVector<double>& sleepStages) {
        return !sleepStages.isEmpty()
            && !(sleepStages.size() == 1 && sleepStages[0] == -1.0);
    }

}  // namespace

// ============================================================================
// Channel lookup table
// ============================================================================

const QStringList& noise_marking_gui::markableChannelLabels() {
    static const QStringList kLabels{ "ECG1", "ECG2", "ECG3", "PPG", "ABP" };
    return kLabels;
}

noise_marking_gui::ChannelRefs noise_marking_gui::channelRefs(const QString& label) const {
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

noise_marking_gui::ChannelMarkingState&
noise_marking_gui::markStateFor(const QString& label) {
    ChannelRefs r = channelRefs(label);
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

double noise_marking_gui::yScaleForSignal(const QString& label) const {
    QCheckBox* check = nullptr;
    QDoubleSpinBox* gain = nullptr;
    if (label == "ECG1") { check = ui->ecg_1_check; gain = ui->ecg_1_gain; }
    else if (label == "ECG2") { check = ui->ecg_2_check; gain = ui->ecg_2_gain; }
    else if (label == "ECG3") { check = ui->ecg_3_check; gain = ui->ecg_3_gain; }
    else if (label == "PPG") { check = ui->ppg_check;   gain = ui->ppg_gain; }
    else if (label == "ABP" || label == "ACCEL") {
        check = ui->abg_check; gain = ui->abg_gain;
    }
    if (!gain) return 1.0;
    double v = gain->value();
    return (v > 0.0) ? v : 1.0;
}

void noise_marking_gui::resetUnpinnedGains() {
    // For each pair, if the checkbox is unchecked, snap the spinbox back to
    // 1.0 without firing handle_data_plot (we're about to redraw anyway).
    auto reset = [](QCheckBox* check, QDoubleSpinBox* gain) {
        if (!check || !gain) return;
        if (check->isChecked()) return;
        QSignalBlocker block(gain);
        gain->setValue(1.0);
        };
    reset(ui->ecg_1_check, ui->ecg_1_gain);
    reset(ui->ecg_2_check, ui->ecg_2_gain);
    reset(ui->ecg_3_check, ui->ecg_3_gain);
    reset(ui->ppg_check, ui->ppg_gain);
    reset(ui->abg_check, ui->abg_gain);
}

double noise_marking_gui::totalChunkDuration() const {
    if (m_ecg1.size() > 1 && m_ecgSR > 0) return m_ecg1.size() / m_ecgSR;
    if (m_ppg.size() > 1 && m_ppgSR > 0) return m_ppg.size() / m_ppgSR;
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

    if (!isChannelActive(label)) {
        startBtn->setEnabled(false);
        stopBtn->setEnabled(false);
        startBtn->setStyleSheet("color: gray;");
        stopBtn->setStyleSheet("color: gray;");
        return;
    }

    const char* startSheet = "";
    const char* stopSheet = "";
    bool stopEnabled = false;

    switch (markStateFor(label).phase) {
    case MarkPhase::WaitingForStart:
        startSheet = "background-color: #f39c12; color: white;";
        break;
    case MarkPhase::WaitingForEnd:
        startSheet = "background-color: #f39c12; color: white;";
        stopEnabled = true;
        break;
    case MarkPhase::WaitingForStop:
        stopSheet = "background-color: #e74c3c; color: white;";
        stopEnabled = true;
        break;
    case MarkPhase::Idle:
    default:
        break;
    }

    startBtn->setEnabled(true);
    startBtn->setStyleSheet(startSheet);
    stopBtn->setEnabled(stopEnabled);
    stopBtn->setStyleSheet(stopSheet);
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

    bool anyWaitingEnd = false, anyWaitingStop = false;
    for (const QString& lbl : markableChannelLabels()) {
        if (!isChannelActive(lbl)) continue;
        MarkPhase p = markStateFor(lbl).phase;
        if (p == MarkPhase::WaitingForEnd)  anyWaitingEnd = true;
        if (p == MarkPhase::WaitingForStop) anyWaitingStop = true;
    }

    if (anyWaitingStop) {
        ui->start_all_mark->setStyleSheet("");
        ui->stop_all_mark->setEnabled(true);
        ui->stop_all_mark->setStyleSheet("background-color: #e74c3c; color: white;");
    }
    else {
        ui->start_all_mark->setStyleSheet("background-color: #f39c12; color: white;");
        ui->stop_all_mark->setEnabled(anyWaitingEnd);
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
        resetUnpinnedGains();
        handle_data_plot();
        updateAmpogramCursor();
        });

    new QShortcut(QKeySequence(Qt::Key_Right), this, [this]() {
        double maxStart = std::max(0.0, totalChunkDuration() - m_windowDuration);
        m_currentStartTime = std::min(m_currentStartTime + m_skipInterval, maxStart);
        resetUnpinnedGains();
        handle_data_plot();
        updateAmpogramCursor();
        });

    // --- All chart views: scrollbars off, antialiasing on, event filter. ---
    const QList<QChartView*> allCharts = {
        ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3,
        ui->ppg_axis, ui->accel_or_abg_axis,
        ui->ecg_ampogram_axis, ui->amp_ppg_axis, ui->resp_cvp_axis
    };
    for (auto* view : allCharts) {
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


    auto wireGain = [this](QCheckBox* check, QDoubleSpinBox* gain) {
        if (!gain) return;
        gain->setDecimals(2);
        gain->setSingleStep(0.5);
        gain->setRange(0.1, 100.0);
        gain->setValue(1.0);
        // ClickFocus: spinbox takes focus only when explicitly clicked, and
        // arrow-key navigation can't accidentally land in it. Combined with
        // editingFinished -> clearFocus, focus reliably returns to the dialog
        // after the user types or steps a new value.
        gain->setFocusPolicy(Qt::ClickFocus);
        connect(gain, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, [this](double) { handle_data_plot(); });
        connect(gain, &QDoubleSpinBox::editingFinished,
            this, [gain]() { gain->clearFocus(); });
        };
    wireGain(ui->ecg_1_check, ui->ecg_1_gain);
    wireGain(ui->ecg_2_check, ui->ecg_2_gain);
    wireGain(ui->ecg_3_check, ui->ecg_3_gain);
    wireGain(ui->ppg_check, ui->ppg_gain);
    wireGain(ui->abg_check, ui->abg_gain);



    // Set up default chart objects + reflow titles when plot area settles.
    // plotAreaChanged fires after QChart updates layout, so plotArea() is
    // fresh -- unlike Resize, which fires first and reads stale coords.
    for (QChartView* v : allCharts) {
        if (!v) continue;
        setupChartDefaults(v);
        connect(v->chart(), &QChart::plotAreaChanged, v,
            [v, this](const QRectF&) {
                const QString sigName = v->property("signalName").toString();
                if (sigName.isEmpty()) return;
                const double nativeHz = v->property("nativeHz").toDouble();
                const double pxPerSec = (m_windowDuration > 0.0)
                    ? v->chart()->plotArea().width() / m_windowDuration : 0.0;
                const double pxPerSample = (nativeHz > 0.0) ? pxPerSec / nativeHz : 0.0;
                v->chart()->setTitle(formatChartTitle(sigName, nativeHz, pxPerSample));
            });
    }

    // Ampogram series (ECG and PPG only).
    ecg1_ampogram_series = new QLineSeries();   ecg1_ampogram_series->setName("ECG1");
    ecg2_ampogram_series = new QLineSeries();   ecg2_ampogram_series->setName("ECG2");
    ecg3_ampogram_series = new QLineSeries();   ecg3_ampogram_series->setName("ECG3");
    ppg_ampogram_series = new QLineSeries();
    ui->ecg_ampogram_axis->chart()->addSeries(ecg1_ampogram_series);
    ui->ecg_ampogram_axis->chart()->addSeries(ecg2_ampogram_series);
    ui->ecg_ampogram_axis->chart()->addSeries(ecg3_ampogram_series);
    ui->amp_ppg_axis->chart()->addSeries(ppg_ampogram_series);

    // Initially disable all stop buttons.
    ui->stop_ecg1_mark->setEnabled(false);
    ui->stop_ecg2_mark->setEnabled(false);
    ui->stop_ecg3_mark->setEnabled(false);
    ui->stopNoisePPG->setEnabled(false);
    ui->stop_all_mark->setEnabled(false);

    // Cursor bars (black vertical line tracking the currently-viewed window).
    auto addCursor = [](QChartView* view, QLineSeries*& series) {
        series = new QLineSeries();
        series->setPen(QPen(Qt::black, 2));
        view->chart()->addSeries(series);
        };
    addCursor(ui->ecg_ampogram_axis, m_ecgCursorBar);
    addCursor(ui->amp_ppg_axis, m_ppgCursorBar);

    // Hypnogram chart (replaces resp_cvp_axis chart object).
    auto* hypnoChart = new QChart();
    hypnoChart->legend()->hide();
    ui->resp_cvp_axis->setChart(hypnoChart);

    m_hypnoCursorBar = new QLineSeries();
    m_hypnoCursorBar->setPen(QPen(Qt::black, 2));
    hypnoChart->addSeries(m_hypnoCursorBar);

    m_currentMarkingType = ui->marking_type->currentText();

    connect(ui->browse_file_button, &QPushButton::clicked,
        this, &noise_marking_gui::handleBrowseFile);

    // Annealing pipeline button.
    connect(ui->process_button, &QPushButton::clicked, this, [this]() {
        if (m_cfg.bin_file_path.empty() || m_cfg.annealed_data_path.empty()) {
            QMessageBox::warning(this, "Process Output",
                "Config not set or annealedDataPath missing in config.csv.");
            return;
        }

        ui->process_button->setEnabled(false);
        QApplication::setOverrideCursor(Qt::WaitCursor);

        auto* watcher = new QFutureWatcher<int>(this);
        connect(watcher, &QFutureWatcher<int>::finished, this, [this, watcher]() {
            QApplication::restoreOverrideCursor();
            ui->process_button->setEnabled(true);
            QMessageBox::information(this, "Process Output",
                QString("Pipeline complete. Processed %1 files.").arg(watcher->result()));
            watcher->deleteLater();
            });

        config_entry cfgCopy = m_cfg;
        watcher->setFuture(QtConcurrent::run([cfgCopy]() {
            return processDataset(cfgCopy);
            }));
        });
}

noise_marking_gui::~noise_marking_gui() = default;

GenExcStruct noise_marking_gui::getMarkings() const {
    GenExcStruct result = m_genExc;
    result.filePath = m_binFilePath;
    return result;
}

QVector<GenExcStruct> noise_marking_gui::getAllMarkings() const {
    QMap<QString, GenExcStruct> all = m_fileMarkings;
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

    constexpr int kNumChannels = NUM_CHANNELS;
    constexpr int kNumFields = 4 + 3 * kNumChannels + 1;
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

    constexpr int kSizesUpBase = 4;
    constexpr int kSizesRawBase = kSizesUpBase + kNumChannels;
    constexpr int kNativeRatesBase = kSizesRawBase + kNumChannels;
    constexpr int kSleepCountIdx = kNativeRatesBase + kNumChannels;

    for (int i = 0; i < kNumChannels; ++i) {
        m_chanSizes[i] = raw32[kSizesUpBase + i];
        m_chanSizesRaw[i] = raw32[kSizesRawBase + i];
        std::memcpy(&m_chanNativeRates[i],
            &raw32[kNativeRatesBase + i], sizeof(float));
    }
    m_totalSleepSamples = raw32[kSleepCountIdx];

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
                m_genExc.marking_type[i].toStdString());
        }
    }
    else {
        m_genExc = GenExcStruct();
        m_genExc.filePath = filePath;
        m_noiseManager = std::make_unique<NoiseManager>(m_ecgSR);
    }

    m_currentStartTime = 0.0;
    m_markAllActive = false;

    for (const QString& lbl : markableChannelLabels()) cancelMarking(lbl);

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

    // ------------------------------------------------------------------
    // Compute per-channel byte offsets. On disk each slot is laid out as
    // (upsampled block, raw block); raw block is (t, v) pairs so its byte
    // length is 2 * m_chanSizesRaw[i] * sizeof(double).
    // ------------------------------------------------------------------
    uint64_t chanUpOffset[NUM_CHANNELS], chanRawOffset[NUM_CHANNELS];
    uint64_t running = 0;
    for (int i = 0; i < NUM_CHANNELS; ++i) {
        chanUpOffset[i] = running;  running += m_chanSizes[i];
        chanRawOffset[i] = running;  running += m_chanSizesRaw[i] * 2;
    }
    const uint64_t sleepByteOffset = running;

    auto rateForChannel = [this](int chIdx) -> double {
        if (chIdx == CH_MARKER || chIdx == CH_TEMP || chIdx == CH_PACEMAKER)
            return m_boolSR;
        if (chIdx >= CH_EKG_OFF && chIdx <= CH_EEG3_OFF)
            return m_boolSR;
        if (chIdx == CH_OXSTATUS || chIdx == CH_SPO2 || chIdx == CH_HR)
            return m_boolSR;
        return m_ecgSR;
        };

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

    // Stream raw (t, v) pair block, keeping only pairs in the current chunk.
    // Pairs are time-sorted, so we bail out as soon as t >= chunkEnd.
    auto loadRaw = [&](QVector<QPointF>& dest, int chIdx) {
        dest.clear();
        uint64_t totalPairs = m_chanSizesRaw[chIdx];

        // Sentinel: file stores a single (-1.0, -1.0) pair for missing.
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

        constexpr uint64_t BLOCK_PAIRS = 1 << 15;
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

    // Upsampled signals.
    loadSignal(m_ecg1, CH_ECG1);
    loadSignal(m_ecg2, CH_ECG2);
    loadSignal(m_ecg3, CH_ECG3);
    loadSignal(m_ppg, CH_PPG);
    loadSignal(m_accelX, CH_ACCEL_X);
    loadSignal(m_accelY, CH_ACCEL_Y);
    loadSignal(m_accelZ, CH_ACCEL_Z);
    loadSignal(m_cvp, CH_PRES);
    loadSignal(m_resp, CH_RESP);
    loadSignal(m_abp, CH_ABP);

    // Raw (t, v) overlays.
    loadRaw(m_ecg1Raw, CH_ECG1);
    loadRaw(m_ecg2Raw, CH_ECG2);
    loadRaw(m_ecg3Raw, CH_ECG3);
    loadRaw(m_ppgRaw, CH_PPG);
    loadRaw(m_abpRaw, CH_ABP);
    loadRaw(m_accelXRaw, CH_ACCEL_X);
    loadRaw(m_accelYRaw, CH_ACCEL_Y);
    loadRaw(m_accelZRaw, CH_ACCEL_Z);
    loadRaw(m_respRaw, CH_RESP);
    loadRaw(m_cvpRaw, CH_PRES);

    // Sleep stages.
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
    m_currentStartTime = 0;

    // ------------------------------------------------------------------
    // Recompute which channels are active & set chart visibility.
    // ------------------------------------------------------------------
    m_activeChannels.clear();
    auto markActive = [this](const QString& label, const QVector<double>& data) {
        bool missing = isMissingSignal(data);
        if (auto* cv = chartViewForSignalLabel(label))
            cv->setVisible(!missing);
        if (!missing) m_activeChannels.insert(label);
        };
    markActive("ECG1", m_ecg1);
    markActive("ECG2", m_ecg2);
    markActive("ECG3", m_ecg3);
    markActive("PPG", m_ppg);

    bool anyAccel = !isMissingSignal(m_accelX)
        || !isMissingSignal(m_accelY)
        || !isMissingSignal(m_accelZ);
    if (!anyAccel) markActive("ABP", m_abp);

    if (ui->accel_or_abg_axis)
        ui->accel_or_abg_axis->setVisible(!isMissingSignal(m_abp));
    if (ui->amp_ppg_axis)
        ui->amp_ppg_axis->setVisible(!isMissingSignal(m_ppg));

    if (ui->resp_cvp_axis) {
        bool sleepPresent = sleepDataPresent(m_sleepStages);
        bool cvpPresent = !isMissingSignal(m_cvp);
        ui->resp_cvp_axis->setVisible(sleepPresent || cvpPresent || anyAccel);
    }

    updateAllChannelButtonStates();

    handle_ampogram_plot();
    handle_data_plot();
    setupHypnogram();        // must run before cursor update
    updateAmpogramCursor();
    restoreMarkingMarkers();

    uint64_t ecgPerChunk = static_cast<uint64_t>(CHUNK_DURATION_SEC * m_ecgSR);
    ui->prev8hours->setEnabled(chunkIndex > 0);
    ui->next8hours->setEnabled((chunkIndex * ecgPerChunk + m_ecg1.size()) < m_chanSizes[CH_ECG1]);

    return true;
}

void noise_marking_gui::on_next8hours_clicked() { resetUnpinnedGains(); loadChunkFromFile(m_currentChunkIndex + 1); }
void noise_marking_gui::on_prev8hours_clicked() { if (m_currentChunkIndex > 0) { resetUnpinnedGains(); loadChunkFromFile(m_currentChunkIndex - 1); } }

void noise_marking_gui::setupHypnogram() {
    /*
        The hypnogram is specifically for the MESA files.
        It is less essential so it can be small, it has a color code for each sleep state
    */
    if (m_sleepSR <= 0.0) return;
    if (!sleepDataPresent(m_sleepStages)) return;

    auto* chart = ui->resp_cvp_axis->chart();


    if (m_cvpCursorBar && m_cvpCursorBar->chart() == chart)
        chart->removeSeries(m_cvpCursorBar);

    for (auto* s : m_hypnoStageSeries) { chart->removeSeries(s); delete s; }
    m_hypnoStageSeries.clear();

    struct Stage { int value; QColor color; };
    const QList<Stage> stages = {
        {0, Qt::black}, {1, Qt::darkGreen}, {2, Qt::blue}, {3, Qt::cyan}, {4, Qt::red}
    };

    const double dt = 1.0 / m_sleepSR;
    const double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;

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

    clearAxes(chart);

    chart->setTitle("Sleep stages");
    chart->setTitleFont(QFont("Arial", 8, QFont::Bold));
    chart->setTitleBrush(Qt::black);
    chart->setMargins(QMargins(0, 0, 0, 0));


    auto* xAxis = new QCategoryAxis();
    xAxis->setRange(globalOffset, globalOffset + CHUNK_DURATION_SEC);
    const int startHour = static_cast<int>(globalOffset / 3600.0);
    for (int h = 0; h <= 8; h += 2)
        xAxis->append(QString::number(startHour + h) + 'h', globalOffset + h * 3600.0);
    xAxis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
    xAxis->setGridLineVisible(false);
    xAxis->setLabelsFont(QFont("Arial", 7));


    auto* yAxis = new QCategoryAxis();
    for (const auto& st : stages) yAxis->append("", st.value + 0.4);
    yAxis->setRange(-0.5, 4.5);
    yAxis->setReverse(true);
    yAxis->setVisible(false);
    yAxis->setGridLineVisible(false);
    xAxis->setLabelsVisible(true);

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
    // Amplitude variability across `range`-second windows -- highlights noise.
    const double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;

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

    auto create_plot = [globalOffset](
        QChartView* view, QLineSeries* series, const QList<QPointF>& pts,
        QLineSeries* cursor, const QColor& color, const QString& title,
        bool showLabels = false) {
            series->replace(pts);
            series->setPen(QPen(color, 1));

            auto* chart = view->chart();
            clearAxes(chart);
            chart->legend()->hide();

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
            const int startHour = static_cast<int>(globalOffset / 3600.0);
            for (int h = 0; h <= 8; ++h) {
                QString lbl = showLabels
                    ? (QString::number(startHour + h) + 'h')
                    : QString::number(h);
                xAxis->append(lbl, globalOffset + h * 3600.0);
            }
            xAxis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
            xAxis->setGridLineVisible(false);
            xAxis->setLabelsFont(QFont("Arial", 7));
            xAxis->setLabelsVisible(showLabels);

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

    const bool sleepPresent = sleepDataPresent(m_sleepStages);
    const bool ppgAmpHasLabels = !sleepPresent
        && ui->amp_ppg_axis && !ui->amp_ppg_axis->isHidden();
    const bool ecgAmpHasLabels = !sleepPresent && !ppgAmpHasLabels;

    auto* chart = ui->ecg_ampogram_axis->chart();
    if (ecg2_ampogram_series->chart() == chart) chart->removeSeries(ecg2_ampogram_series);
    if (ecg3_ampogram_series->chart() == chart) chart->removeSeries(ecg3_ampogram_series);

    create_plot(ui->ecg_ampogram_axis, ecg1_ampogram_series,
        ecg1Pts, m_ecgCursorBar, COLOR_ECG1, "ECG Amp-O-Gram",
        /*showLabels*/ ecgAmpHasLabels);

    auto* xAxis = chart->axes(Qt::Horizontal).first();
    auto* yAxis = qobject_cast<QValueAxis*>(chart->axes(Qt::Vertical).first());

    auto attachExtra = [&](QLineSeries* series, const QList<QPointF>& pts,
        const QColor& color) {
            if (pts.isEmpty()) return;
            series->replace(pts);
            series->setPen(QPen(color, 1));
            chart->addSeries(series);
            series->attachAxis(xAxis);
            series->attachAxis(yAxis);
        };
    attachExtra(ecg2_ampogram_series, ecg2Pts, COLOR_ECG2);
    attachExtra(ecg3_ampogram_series, ecg3Pts, COLOR_ECG3);

    auto allPts = ecg1Pts + ecg2Pts + ecg3Pts;
    if (yAxis && !allPts.isEmpty()) {
        auto [mi, ma] = std::minmax_element(allPts.begin(), allPts.end(),
            [](const QPointF& a, const QPointF& b) { return a.y() < b.y(); });
        yAxis->setRange(mi->y(), ma->y());
    }

    create_plot(ui->amp_ppg_axis, ppg_ampogram_series,
        calculate_amplitude(m_ppg, m_ppgSR),
        m_ppgCursorBar, COLOR_PPG, "PPG Amp-O-Gram",
        /*showLabels*/ ppgAmpHasLabels);
}

// ============================================================================
// Main signal plot
// ============================================================================

// Small POD describing one series to render on a windowed signal chart.
// Used by both the markable path (one channel per chart) and the display
// path (multi-series accel / RESP / CVP charts).
struct WindowedSeries {
    const QVector<double>* data;       // upsampled samples (required)
    QColor                  color;
    const QVector<QPointF>* rawData = nullptr;  // optional native-rate overlay
};

// ----------------------------------------------------------------------------
// renderWindowedChart
//
// Render N signal series on a chart with:
//   - x-axis: HH:MM:SS labels at 1/4 intervals across the current window
//   - y-axis: auto-ranged with 5%-of-span padding (or constant fallback)
//   - per-series upsampled foreground (line OR scatter, per plotMode)
//   - per-series raw (t, v) scatter overlay drawn last (on top)
//
// Returns the (min, max) of the combined sample set. The caller can use
// this to position a marker line that spans the chart's full y range.
//
// `forceLine` overrides plotMode -- used when the upsampled foreground
// MUST be a line regardless of mode (markable channels use this so the
// black raw-scatter overlay is always layered on top of a colored line).
// `singleRawColor` overrides each series' individual color for the raw
// overlay (markable uses COLOR_RAW_SCATTER for solid black).
// ----------------------------------------------------------------------------
namespace {

    std::pair<double, double> renderWindowedChart(
        QChartView* view,
        const QList<WindowedSeries>& serieses,
        double currentStartTime, double windowDuration,
        double globalOffset, double ecgSR,
        bool labelsVisible,
        bool useScatterMode,
        bool forceLineForUpsampled,
        QColor singleRawColor,
        bool useSingleRawColor,
        double yScale = 1.0)
    {
        if (!view || !view->chart()) return { 1e9, -1e9 };
        QChart* chart = view->chart();

        // Wipe everything except series the caller already preserved by
        // attaching them to other charts (we don't pass keep here -- keep
        // is handled at the higher level by wipeChartContent before this).
        for (auto* s : chart->series()) { chart->removeSeries(s); delete s; }
        for (auto* a : chart->axes()) { chart->removeAxis(a);  delete a; }

        auto* xAxis = makeWindowedTimeAxis(currentStartTime, windowDuration,
            globalOffset, labelsVisible);
        chart->addAxis(xAxis, Qt::AlignBottom);
        chart->setMargins(QMargins(0, 0, 20, 0));

        auto* yAxis = new QValueAxis();
        yAxis->setVisible(false);
        chart->addAxis(yAxis, Qt::AlignLeft);

        double gMin = 1e9;
        double gMax = -1e9;

        for (const auto& d : serieses) {
            if (!d.data || isMissingSignal(*d.data)) continue;

            const bool hasRaw = d.rawData && isRawUsable(*d.rawData);

            // Decide what the upsampled foreground looks like:
            //  - forceLine: always line
            //  - else if scatter mode + raw available: skip foreground entirely
            //    (the raw dots are the only markers; we still scan for Y range)
            //  - else if scatter mode + no raw: scatter the upsampled samples
            //  - else (line mode): line
            QXYSeries* plotSeries = nullptr;
            if (forceLineForUpsampled) {
                auto* ln = new QLineSeries();
                ln->setUseOpenGL(true);
                ln->setPen(QPen(d.color, 1));
                chart->addSeries(ln);
                plotSeries = ln;
            }
            else if (useScatterMode && !hasRaw) {
                auto* sc = new QScatterSeries();
                sc->setColor(d.color);
                sc->setBorderColor(Qt::transparent);
                sc->setMarkerSize(2.0);
                sc->setMarkerShape(QScatterSeries::MarkerShapeCircle);
                sc->setUseOpenGL(true);
                chart->addSeries(sc);
                plotSeries = sc;
            }
            else if (!useScatterMode) {
                auto* ln = new QLineSeries();
                ln->setUseOpenGL(true);
                ln->setPen(QPen(d.color, 1));
                chart->addSeries(ln);
                plotSeries = ln;
            }

            // Upsampled samples (and Y-range scan).
            int startIdx = std::clamp(static_cast<int>(currentStartTime * ecgSR),
                0, static_cast<int>(d.data->size() - 1));
            int endIdx = std::clamp(static_cast<int>((currentStartTime + windowDuration) * ecgSR),
                0, static_cast<int>(d.data->size()));

            QList<QPointF> pts;
            pts.reserve(endIdx - startIdx);
            for (int i = startIdx; i < endIdx; ++i) {
                double raw = (*d.data)[i];
                pts.append({ static_cast<double>(i) / ecgSR, raw * yScale });
                if (raw < gMin) gMin = raw;
                if (raw > gMax) gMax = raw;
            }
            if (plotSeries) {
                plotSeries->replace(pts);
                plotSeries->attachAxis(xAxis);
                plotSeries->attachAxis(yAxis);
            }

            // Raw scatter overlay (drawn after so it sits on top).
            if (hasRaw) {
                auto* rawScatter = new QScatterSeries();
                rawScatter->setColor(useSingleRawColor ? singleRawColor : d.color);
                rawScatter->setBorderColor(Qt::transparent);
                rawScatter->setMarkerSize(useSingleRawColor ? 3.0 : 3.5);
                rawScatter->setMarkerShape(QScatterSeries::MarkerShapeCircle);
                rawScatter->setUseOpenGL(true);
                chart->addSeries(rawScatter);

                const double viewStart = currentStartTime;
                const double viewEnd = currentStartTime + windowDuration;

                QList<QPointF> rawPts;
                rawPts.reserve(std::min<int>(d.rawData->size(), 4096));
                for (const QPointF& p : *d.rawData) {
                    if (p.x() < viewStart) continue;
                    if (p.x() > viewEnd)   break;
                    rawPts.append({ p.x(), p.y() * yScale });
                    if (p.y() < gMin) gMin = p.y();
                    if (p.y() > gMax) gMax = p.y();
                }
                rawScatter->replace(rawPts);
                rawScatter->attachAxis(xAxis);
                rawScatter->attachAxis(yAxis);
            }
        }

        setPaddedYRange(yAxis, gMin, gMax);
        return { gMin, gMax };
    }

}  // namespace

void noise_marking_gui::handle_data_plot() {
    // Clear noise highlights before any chart wipes.
    for (auto* area : m_highlights) {
        if (area->chart()) area->chart()->removeSeries(area);
        delete area->upperSeries();
        delete area->lowerSeries();
        delete area;
    }
    m_highlights.clear();

    // Bottom-most visible chart in each column gets the time-ruler labels.
    QChartView* xLabelOwnerRight = nullptr;
    QChartView* xLabelOwnerLeft = nullptr;
    {
        const QList<QChartView*> rightCol = {
            ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3,
            ui->ppg_axis,   ui->accel_or_abg_axis
        };
        for (auto* cv : rightCol)
            if (cv && cv->isVisible()) xLabelOwnerRight = cv;

        if (ui->resp_cvp_axis && ui->resp_cvp_axis->isVisible())
            xLabelOwnerLeft = ui->resp_cvp_axis;
    }

    // Wipe charts up front so stale series from other plot modes / inactive
    // channels don't linger. Preserve any in-progress start markers.
    auto preservedMarker = [](const ChannelMarkingState& st) -> QList<QAbstractSeries*> {
        return st.startMarkerLine ? QList<QAbstractSeries*>{ st.startMarkerLine }
        : QList<QAbstractSeries*>{};
        };
    wipeChartContent(ui->ecg_axis_1->chart(), preservedMarker(m_markState_ecg1));
    wipeChartContent(ui->ecg_axis_2->chart(), preservedMarker(m_markState_ecg2));
    wipeChartContent(ui->ecg_axis_3->chart(), preservedMarker(m_markState_ecg3));
    wipeChartContent(ui->ppg_axis->chart(), preservedMarker(m_markState_ppg));

    // resp_cvp_axis is owned here only when no sleep stages -- otherwise
    // setupHypnogram() owns it and we leave it alone.
    const bool sleepPresent = sleepDataPresent(m_sleepStages);
    if (!sleepPresent && ui->resp_cvp_axis)
        wipeChartContent(ui->resp_cvp_axis->chart());

    // ----------------------------------------------------------------------
    // Render one markable channel: line + raw black scatter overlay,
    // chart title with native Hz + px/s, and a marker line if active.
    // ----------------------------------------------------------------------
    auto plotMarkable = [&](const QString& label) {
        if (!isChannelActive(label)) return;
        ChannelRefs r = channelRefs(label);
        if (!r.chartView || !r.data || !r.state) return;

        const double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
        const double sr = (r.sampleRate ? *r.sampleRate : m_ecgSR);
        const QVector<QPointF> emptyRaw;
        const QVector<QPointF>& rawData = r.dataRaw ? *r.dataRaw : emptyRaw;

        // Compute native Hz from the raw block timestamps when available.
        double nativeHz = sr;
        if (rawData.size() >= 2 && rawData.last().x() > rawData.first().x()
            && !(rawData.size() == 1 && rawData[0].x() == -1.0)) {
            nativeHz = (rawData.size() - 1) / (rawData.last().x() - rawData.first().x());
        }

        // Title (channel name + native Hz + px/sample).
        const double pxPerSec = (m_windowDuration > 0.0)
            ? r.chartView->chart()->plotArea().width() / m_windowDuration : 0.0;
        const double pxPerSample = (nativeHz > 0.0) ? pxPerSec / nativeHz : 0.0;
        r.chartView->setProperty("signalName", label);
        r.chartView->setProperty("nativeHz", nativeHz);
        r.chartView->chart()->setTitle(formatChartTitle(label, nativeHz, pxPerSample));
        r.chartView->chart()->setTitleFont(QFont("Arial", 8, QFont::Bold));
        r.chartView->chart()->setTitleBrush(r.color);

        // Render. Markable channels always force a colored line for the
        // upsampled foreground; raw overlay is solid black (COLOR_RAW_SCATTER).
        QList<WindowedSeries> serieses = { { r.data, r.color, &rawData } };
        auto [yMin, yMax] = renderWindowedChart(
            r.chartView, serieses,
            m_currentStartTime, m_windowDuration, globalOffset, sr,
            /*labelsVisible*/ r.chartView == xLabelOwnerRight,
            /*useScatterMode*/ false,             // ignored when forceLine=true
            /*forceLineForUpsampled*/ m_plotMode == PlotMode::Line,
            COLOR_RAW_SCATTER, /*useSingleRawColor*/ true,
            /*yScale*/ yScaleForSignal(label));

        // Place the start-marker line if this channel has one in progress.
        if (r.state->startMarkerLine && r.state->startMarkerLine->chart() == r.chartView->chart()) {
            const double localMarkerPos = r.state->globalStartTime - globalOffset;
            auto* yAxis = qobject_cast<QValueAxis*>(
                r.chartView->chart()->axes(Qt::Vertical).first());
            r.state->startMarkerLine->replace(
                { {localMarkerPos, yAxis->min()}, {localMarkerPos, yAxis->max()} });
            r.state->startMarkerLine->attachAxis(
                r.chartView->chart()->axes(Qt::Horizontal).first());
            r.state->startMarkerLine->attachAxis(yAxis);
        }
        };

    plotMarkable("ECG1");
    plotMarkable("ECG2");
    plotMarkable("ECG3");
    plotMarkable("PPG");

    // ----------------------------------------------------------------------
    // Display-only chart: 1..N series, no marker / click. RESP/CVP/accel.
    // ----------------------------------------------------------------------
    auto plotDisplay = [&](QChartView* view, const QString& title,
        const QList<WindowedSeries>& serieses) {
            if (!view || !view->chart()) return;
            view->chart()->setTitle(title);
            view->chart()->setTitleFont(QFont("Arial", 8, QFont::Bold));

            const double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
            const bool labelsVisible = (view == xLabelOwnerRight) || (view == xLabelOwnerLeft);

            renderWindowedChart(
                view, serieses,
                m_currentStartTime, m_windowDuration, globalOffset, m_ecgSR,
                labelsVisible,
                /*useScatterMode*/ m_plotMode == PlotMode::Scatter,
                /*forceLineForUpsampled*/ false,
                QColor(), /*useSingleRawColor*/ false);
        };

    // accel_or_abg_axis: ABP only -- routed through the markable path so
    // it gets the start marker, click-to-mark, and noise-highlight overlay.
    if (ui->accel_or_abg_axis && !isMissingSignal(m_abp))
        plotMarkable("ABP");

    // resp_cvp_axis: accel when present; otherwise windowed RESP/CVP.
    const bool anyAccel = !isMissingSignal(m_accelX)
        || !isMissingSignal(m_accelY)
        || !isMissingSignal(m_accelZ);
    if (anyAccel && ui->resp_cvp_axis) {
        plotDisplay(ui->resp_cvp_axis, "ACCEL", {
            { &m_accelX, COLOR_ACCEL_X, &m_accelXRaw },
            { &m_accelY, COLOR_ACCEL_Y, &m_accelYRaw },
            { &m_accelZ, COLOR_ACCEL_Z, &m_accelZRaw },
            });
    }
    else if (!sleepPresent && ui->resp_cvp_axis &&
        (!isMissingSignal(m_resp) || !isMissingSignal(m_cvp))) {
        plotDisplay(ui->resp_cvp_axis, "RESP / CVP", {
            { &m_resp, COLOR_RESP, &m_respRaw },
            { &m_cvp,  COLOR_CVP,  &m_cvpRaw  },
            });
    }

    updateNoiseHighlights();
}

// ============================================================================
// Marking
// ============================================================================

void noise_marking_gui::finalizeMarking(QChartView* /*cv*/, double endX,
    const QString& signalLabel) {
    ChannelMarkingState& state = markStateFor(signalLabel);
    const double sr = sampleRateForSignal(signalLabel);

    const double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
    const double globalEnd = endX + globalOffset;
    const double globalStart = state.globalStartTime;

    const double snappedS = std::round(std::min(globalStart, globalEnd) * sr) / sr;
    const double snappedE = std::round(std::max(globalStart, globalEnd) * sr) / sr;

    m_noiseManager->addSegment(
        static_cast<size_t>(snappedS * sr),
        static_cast<size_t>(snappedE * sr),
        signalLabel.toStdString(),
        m_currentMarkingType.toStdString());

    m_genExc.noiseExc.append({ snappedS, snappedE });
    m_genExc.data_type.append(signalLabel);
    m_genExc.marking_type.append(m_currentMarkingType);

    state.phase = MarkPhase::Idle;
    clearStartMarker(state);

    if (QPushButton* startBtn = startButtonForSignal(signalLabel))
        startBtn->setStyleSheet("");
    if (QPushButton* stopBtn = stopButtonForSignal(signalLabel)) {
        stopBtn->setStyleSheet("");
        stopBtn->setEnabled(false);
    }

    // Mark-all: when every active channel is back to Idle, drop out of mark-all.
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

    if (event->type() == QEvent::MouseButtonPress)
        return handleMousePress(cv, viewport, static_cast<QMouseEvent*>(event))
        ? true : QDialog::eventFilter(watched, event);

    if (event->type() == QEvent::MouseMove && m_isDragging) return true;

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

// Returns true if the press was consumed.
bool noise_marking_gui::handleMousePress(QChartView* cv, QWidget* viewport,
    QMouseEvent* me) {
    if (me->button() != Qt::LeftButton) return false;

    double clickedX = cv->chart()->mapToValue(me->pos()).x();
    const double chunkDur = totalChunkDuration();
    clickedX = std::clamp(clickedX, 0.0, chunkDur);

    const QString label = signalLabelForChartView(cv);
    if (!label.isEmpty()) {
        if (!isChannelActive(label)) return false;
        const double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;

        // --- Mark-all mode ---
        if (m_markAllActive) {
            bool anyWaitingForFirstClick = false, anyWaitingForEnd = false, anyWaitingStop = false;
            for (const QString& lbl : markableChannelLabels()) {
                if (!isChannelActive(lbl)) continue;
                MarkPhase p = markStateFor(lbl).phase;
                if (p == MarkPhase::WaitingForStart) anyWaitingForFirstClick = true;
                if (p == MarkPhase::WaitingForEnd)   anyWaitingForEnd = true;
                if (p == MarkPhase::WaitingForStop)  anyWaitingStop = true;
            }

            // 1) Finalize any channels in WaitingForStop first.
            if (anyWaitingStop) {
                for (const QString& lbl : markableChannelLabels()) {
                    if (!isChannelActive(lbl)) continue;
                    if (markStateFor(lbl).phase != MarkPhase::WaitingForStop) continue;
                    finalizeMarking(chartViewForSignalLabel(lbl), clickedX, lbl);
                }
                return true;
            }

            // 2) Place / move start markers on channels still waiting.
            if (anyWaitingForFirstClick || anyWaitingForEnd) {
                for (const QString& lbl : markableChannelLabels()) {
                    if (!isChannelActive(lbl)) continue;
                    ChannelMarkingState& st = markStateFor(lbl);
                    if (st.phase != MarkPhase::WaitingForStart &&
                        st.phase != MarkPhase::WaitingForEnd)  continue;
                    st.globalStartTime = clickedX + globalOffset;
                    showStartMarker(chartViewForSignalLabel(lbl), clickedX, st,
                        colorForSignal(lbl), stopButtonForSignal(lbl));
                    st.phase = MarkPhase::WaitingForEnd;
                }
                ui->stop_all_mark->setEnabled(true);
                return true;
            }
        }

        // --- Single-channel mode ---
        ChannelMarkingState& state = markStateFor(label);
        switch (state.phase) {
        case MarkPhase::WaitingForStart:
        case MarkPhase::WaitingForEnd:
            state.globalStartTime = clickedX + globalOffset;
            showStartMarker(cv, clickedX, state, colorForSignal(label),
                stopButtonForSignal(label));
            state.phase = MarkPhase::WaitingForEnd;
            return true;

        case MarkPhase::WaitingForStop:
            finalizeMarking(cv, clickedX, label);
            return true;

        case MarkPhase::Idle:
            m_isDragging = true;
            m_dragStartPos = me->pos();
            m_dragSignalLabel = label;
            state.globalStartTime = clickedX + globalOffset;
            if (!m_draggedViewport) {
                m_draggedViewport = viewport;
                m_draggedViewport->grabMouse();
            }
            return true;
        }
    }

    // --- Overview-chart click (navigate the marking window) ---
    const bool sleepPresent = sleepDataPresent(m_sleepStages);
    const bool isNavChart =
        (cv == ui->ecg_ampogram_axis)
        || (cv == ui->amp_ppg_axis)
        || (cv == ui->resp_cvp_axis && sleepPresent);
    if (isNavChart) {
        const double globalClickX = cv->chart()->mapToValue(me->pos()).x();
        const double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
        const double localTarget = globalClickX - globalOffset;
        const double maxStart = std::max(0.0, chunkDur - m_windowDuration);
        m_currentStartTime = std::clamp(localTarget - m_windowDuration / 2.0,
            0.0, maxStart);
        handle_data_plot();
        updateAmpogramCursor();
        return true;
    }

    return false;
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

    if (sleepDataPresent(m_sleepStages))
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

    const double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
    const double viewStart = m_currentStartTime;
    const double viewEnd = viewStart + m_windowDuration;

    for (const auto& seg : m_noiseManager->getSegments()) {
        QString segLabel = QString::fromStdString(seg.label);
        if (!axesMap.contains(segLabel)) continue;

        const double sr = sampleRateForSignal(segLabel);
        const double segStart = seg.startSample / sr - globalOffset;
        const double segEnd = seg.endSample / sr - globalOffset;
        if (segEnd < viewStart || segStart > viewEnd) continue;

        const double ds = std::max(segStart, viewStart);
        const double de = std::min(segEnd, viewEnd);
        const QColor color = MARKING_COLORS.value(
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
    ChannelMarkingState& state,
    const QColor& color, QPushButton* stopBtn) {
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
    return formatHMS(seconds);
}

// ============================================================================
// Restore start markers after chunk change
// ============================================================================

void noise_marking_gui::restoreMarkingMarkers() {
    const double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
    const double chunkEnd = globalOffset + CHUNK_DURATION_SEC;

    auto restore = [&](const QString& label, ChannelMarkingState& state) {
        if (state.phase != MarkPhase::WaitingForEnd &&
            state.phase != MarkPhase::WaitingForStop) return;

        state.startMarkerLine = nullptr;
        QChartView* cv = chartViewForSignalLabel(label);
        if (!cv || !isChannelActive(label)) return;

        if (state.globalStartTime >= globalOffset && state.globalStartTime <= chunkEnd) {
            const double localX = state.globalStartTime - globalOffset;
            showStartMarker(cv, localX, state, colorForSignal(label),
                stopButtonForSignal(label));
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
    if (m_markAllActive) {
        for (const QString& label : markableChannelLabels())
            if (isChannelActive(label)) cancelMarking(label);
        m_markAllActive = false;
        ui->start_all_mark->setStyleSheet("");
        ui->stop_all_mark->setStyleSheet("");
        ui->stop_all_mark->setEnabled(false);
        return;
    }

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