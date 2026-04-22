/**
 * @file   gui_handler.cpp
 * @brief  Entry point of the noise marking program. Makes a GUI where you can
 *         mark the noise you see in the ECG1/ECG2/ECG3/PPG, and outputs a csv
 *         and a .bin file. Also displays acceleration data (read-only).
 *
 * Input .bin format (160-byte header, 40 × uint32):
 *   Offset  0: signal_rate    (uint32) — upsampled signal rate (1000 Hz)
 *   Offset  4: boolean_rate   (uint32) — boolean channel rate (1 Hz)
 *   Offset  8: pacemaker_rate (uint32) — pacemaker epoch rate
 *   Offset 12: sleep_rate     (uint32) — sleep epoch duration in seconds
 *   Offset 16–156: 36 × uint32 channel sizes (ECG1..Resp)
 *   Offset 156: size_sleep    (uint32) — sleep stage epoch count
 *
 * Signal data (contiguous doubles after header):
 *   36 channels in header order, then sleep stages.
 *   Missing channels are stored as a single -1.0 with size = 1.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-22
 */
#include "gui_handler.h"

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
#include <algorithm>

 // ============================================================================
 // Constants
 // ============================================================================

static const QColor COLOR_ECG1 = QColor("#1ABC9C");
static const QColor COLOR_ECG2 = QColor("#3498DB");
static const QColor COLOR_ECG3 = QColor("#9B59B6");
static const QColor COLOR_PPG = QColor("#E74C3C");
static const QColor COLOR_ACCEL_X = QColor("#F39C12");  // orange
static const QColor COLOR_ACCEL_Y = QColor("#27AE60");  // green
static const QColor COLOR_ACCEL_Z = QColor("#8E44AD");  // purple
static const QColor COLOR_RESP = QColor("#16A085");  // teal
static const QColor COLOR_CVP = QColor("#2980B9");  // blue
static const QColor COLOR_ABP = QColor("#C0392B");  // red

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
// Channel state / mapping helpers
// ============================================================================

noise_marking_gui::ChannelMarkingState& noise_marking_gui::markStateFor(const QString& label) {
    if (label == "ECG1") return m_markState_ecg1;
    if (label == "ECG2") return m_markState_ecg2;
    if (label == "ECG3") return m_markState_ecg3;
    if (label == "ABP")  return m_markState_abp;
    return m_markState_ppg;
}

QString noise_marking_gui::signalLabelForChartView(QChartView* cv) const {
    if (cv == ui->ecg_axis_1) return "ECG1";
    if (cv == ui->ecg_axis_2) return "ECG2";
    if (cv == ui->ecg_axis_3) return "ECG3";
    if (cv == ui->ppg_axis)   return "PPG";

    if (cv == ui->accel_or_abg_axis) {
        bool anyAccel = !isMissingSignal(m_accelX)
            || !isMissingSignal(m_accelY)
            || !isMissingSignal(m_accelZ);
        if (!anyAccel && !isMissingSignal(m_abp)) return "ABP";
    }

    return {};
}

QChartView* noise_marking_gui::chartViewForSignalLabel(const QString& label) const {
    if (label == "ECG1") return ui->ecg_axis_1;
    if (label == "ECG2") return ui->ecg_axis_2;
    if (label == "ECG3") return ui->ecg_axis_3;
    if (label == "PPG")  return ui->ppg_axis;
    if (label == "ABP")  return ui->accel_or_abg_axis;
    return nullptr;
}

double noise_marking_gui::sampleRateForSignal(const QString& label) const {
    return (label == "PPG") ? m_ppgSR : m_ecgSR;
}

QColor noise_marking_gui::colorForSignal(const QString& label) const {
    if (label == "ECG1") return COLOR_ECG1;
    if (label == "ECG2") return COLOR_ECG2;
    if (label == "ECG3") return COLOR_ECG3;
    if (label == "ABP")  return COLOR_ABP;
    return COLOR_PPG;
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
    if (label == "ECG1") return ui->start_ecg1_mark;
    if (label == "ECG2") return ui->start_ecg2_mark;
    if (label == "ECG3") return ui->start_ecg3_mark;
    if (label == "PPG")  return ui->startNoisePPG;
    if (label == "ABP")  return ui->startNoiseABP;
    return nullptr;
}

QPushButton* noise_marking_gui::stopButtonForSignal(const QString& label) const {
    if (label == "ECG1") return ui->stop_ecg1_mark;
    if (label == "ECG2") return ui->stop_ecg2_mark;
    if (label == "ECG3") return ui->stop_ecg3_mark;
    if (label == "PPG")  return ui->stopNoisePPG;
    if (label == "ABP")  return ui->stopNoiseABP;
    return nullptr;
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
    for (const QString& label : { "ECG1", "ECG2", "ECG3", "PPG", "ABP" })
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
    for (const QString& lbl : { "ECG1", "ECG2", "ECG3", "PPG", "ABP" }) {
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
        ui->sleep_or_cvp_axis, ui->amp_ecg1_axis, ui->amp_ppg_axis,
        ui->resp_axis,
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

    ui->rb_10s->setChecked(true);
    m_skipInterval = ui->skip_interval_box->text().toDouble();
    if (m_skipInterval <= 0.0) m_skipInterval = 5.0;

    // --- Plot-style toggle (Line / Scatter) ---
    // Uses the radio buttons defined in noise_marking_gui.ui.
    ui->line_plot->setChecked(true);
    ui->line_plot->setFocusPolicy(Qt::NoFocus);
    ui->scatter_plot->setFocusPolicy(Qt::NoFocus);
    {
        auto* plotModeGroup = new QButtonGroup(this);
        plotModeGroup->setExclusive(true);
        plotModeGroup->addButton(ui->line_plot);
        plotModeGroup->addButton(ui->scatter_plot);

        auto onPlotModeToggle = [this](bool checked) {
            if (!checked) return;   // toggled fires on both sides; only react to the checked one
            m_plotMode = ui->scatter_plot->isChecked() ? PlotMode::Scatter : PlotMode::Line;
            handle_data_plot();
            };
        connect(ui->line_plot, &QRadioButton::toggled, this, onPlotModeToggle);
        connect(ui->scatter_plot, &QRadioButton::toggled, this, onPlotModeToggle);
    }

    QList<QChartView*> defaultCharts = { ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3,
                                         ui->ppg_axis, ui->amp_ecg1_axis, ui->amp_ppg_axis,
                                         ui->resp_axis, ui->accel_or_abg_axis };
    for (int i = 0; i < defaultCharts.size(); ++i) {
        if (defaultCharts[i]) setupChartDefaults(defaultCharts[i]);
    }

    // Ampogram / 8h-overview series
    m_ecgAmpSeries = new QLineSeries();
    m_ppgAmpSeries = new QLineSeries();
    m_respAmpSeries = new QLineSeries();
    ui->amp_ecg1_axis->chart()->addSeries(m_ecgAmpSeries);
    ui->amp_ppg_axis->chart()->addSeries(m_ppgAmpSeries);
    ui->resp_axis->chart()->addSeries(m_respAmpSeries);
    // NOTE: m_cvpAmpSeries and m_cvpCursorBar are NOT pre-created in the
    // ctor, because sleep_or_cvp_axis is multi-purpose (hypnogram OR CVP).
    // They are created/destroyed lazily in handle_ampogram_plot().

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
    addCursor(ui->amp_ecg1_axis, m_ecgCursorBar);
    addCursor(ui->amp_ppg_axis, m_ppgCursorBar);
    addCursor(ui->resp_axis, m_respCursorBar);

    // Hypnogram chart
    auto* hypnoChart = new QChart();
    hypnoChart->legend()->hide();
    hypnoChart->setMargins(QMargins(0, 0, 0, 0));
    ui->sleep_or_cvp_axis->setChart(hypnoChart);

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
    // Stash markings for the file we're leaving
    if (!m_binFilePath.isEmpty()) {
        m_genExc.filePath = m_binFilePath;
        m_fileMarkings[m_binFilePath] = m_genExc;
    }

    m_binFilePath = filePath;
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly)) return;

    uint32_t header[45] = {};
    file.read(reinterpret_cast<char*>(header), sizeof(header));
    file.close();

    m_ecgSR = static_cast<double>(header[0]);
    m_boolSR = static_cast<double>(header[1]);
    m_ppgSR = m_ecgSR;
    double sleepEpoch = static_cast<double>(header[3]);
    m_sleepSR = (sleepEpoch > 0) ? (1.0 / sleepEpoch) : 0;

    for (int i = 0; i < 40; ++i)
        m_chanSizes[i] = header[4 + i];
    m_totalSleepSamples = header[4 + 40];

    // Restore markings if we've seen this file before, otherwise start fresh
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

    cancelMarking("ECG1");
    cancelMarking("ECG2");
    cancelMarking("ECG3");
    cancelMarking("PPG");
    cancelMarking("ABP");

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

    // Compute byte offset of each channel's data after the 180-byte header
    uint64_t chanByteOffsets[40];
    uint64_t running = 0;
    for (int i = 0; i < 40; ++i) {
        chanByteOffsets[i] = running;
        running += m_chanSizes[i];
    }
    uint64_t sleepByteOffset = running;

    // Rate for a given channel index
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
        file.seek(FILE_HEADER_SIZE + (chanByteOffsets[chIdx] + start) * sizeof(double));
        file.read(reinterpret_cast<char*>(dest.data()), count * sizeof(double));
        };

    loadSignal(m_ecg1, CH_ECG1);
    loadSignal(m_ecg2, CH_ECG2);
    loadSignal(m_ecg3, CH_ECG3);
    loadSignal(m_ppg, CH_PPG);
    loadSignal(m_accelX, CH_ACCEL_X);
    loadSignal(m_accelY, CH_ACCEL_Y);
    loadSignal(m_accelZ, CH_ACCEL_Z);
    loadSignal(m_cvp, CH_PRES);       // slot 16 — MESA nasal pressure / Bittium CVP
    loadSignal(m_resp, CH_RESP);      // slot 34 — respiration
    loadSignal(m_abp, CH_ABP);        // slot 35 — Bittium arterial blood pressure

    // Sleep stages
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

    // resp_axis: show if respiration has data
    if (ui->resp_axis) {
        ui->resp_axis->setVisible(!isMissingSignal(m_resp));
    }

    // sleep_or_cvp_axis: show if sleep stages OR CVP has data
    if (ui->sleep_or_cvp_axis) {
        bool sleepPresent = !m_sleepStages.isEmpty()
            && !(m_sleepStages.size() == 1 && m_sleepStages[0] == -1.0);
        bool cvpPresent = !isMissingSignal(m_cvp);
        ui->sleep_or_cvp_axis->setVisible(sleepPresent || cvpPresent);
    }

    updateAllChannelButtonStates();

    handle_ampogram_plot();
    handle_data_plot();
    setupHypnogram();        // may replace sleep_or_cvp_axis contents; must run before cursor update
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

    // No sleep stages in this file � the sleep_or_cvp_axis is being
    // used for CVP (handled in handle_data_plot). Do nothing here.
    bool sleepPresent = !m_sleepStages.isEmpty()
        && !(m_sleepStages.size() == 1 && m_sleepStages[0] == -1.0);
    if (!sleepPresent) return;

    auto* chart = ui->sleep_or_cvp_axis->chart();

    // If we're coming from a previous file that had CVP (no sleep stages),
    // the CVP overview series + its cursor are still on this chart.
    // Remove them before we repaint as a hypnogram.
    if (m_cvpAmpSeries && m_cvpAmpSeries->chart() == chart) {
        chart->removeSeries(m_cvpAmpSeries);
        // don't delete; reusable
    }
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

// ============================================================================
// Ampogram
// ============================================================================

void noise_marking_gui::handle_ampogram_plot(double sampling_length) {
    double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;

    auto calcPoints = [sampling_length, globalOffset](
        const QVector<double>& data, double sr) {
            QList<QPointF> pts;
            if (data.isEmpty() || sr <= 0.0) return pts;
            double duration = data.size() / sr;
            for (double t = 0; t <= duration - sampling_length; t += sampling_length) {
                int s = static_cast<int>(t * sr);
                int e = static_cast<int>((t + sampling_length) * sr);
                auto [mi, ma] = std::minmax_element(data.begin() + s, data.begin() + e);
                pts.append({ globalOffset + t, *ma - *mi });
            }
            return pts;
        };

    // ----------------------------------------------------------------------
    // Min/max-per-pixel envelope: for each of `targetPts` time buckets,
    // emit the sample-min and sample-max in time order. Drawn with a
    // QLineSeries, this produces an envelope visually identical to plotting
    // the raw 28.8M-sample signal, but with a few thousand points instead.
    // Used for the 8h raw-signal overviews (RESP, CVP).
    // ----------------------------------------------------------------------
    auto calcEnvelopePoints = [globalOffset](
        const QVector<double>& data, double sr, int targetPts) {
            QList<QPointF> pts;
            if (data.isEmpty() || sr <= 0.0 || targetPts < 1) return pts;

            const int n = data.size();
            if (n <= targetPts * 2) {
                // Few enough samples to plot directly.
                pts.reserve(n);
                for (int i = 0; i < n; ++i)
                    pts.append({ globalOffset + i / sr, data[i] });
                return pts;
            }

            pts.reserve(targetPts * 2);
            for (int b = 0; b < targetPts; ++b) {
                int s = static_cast<int>((int64_t)b * n / targetPts);
                int e = static_cast<int>((int64_t)(b + 1) * n / targetPts);
                if (e > n) e = n;
                if (s >= e) continue;
                auto [miIt, maIt] = std::minmax_element(data.begin() + s, data.begin() + e);
                double tMin = globalOffset + (s + (miIt - (data.begin() + s))) / sr;
                double tMax = globalOffset + (s + (maIt - (data.begin() + s))) / sr;
                // Emit in time order so the line connecting them looks right.
                if (tMin <= tMax) {
                    pts.append({ tMin, *miIt });
                    pts.append({ tMax, *maIt });
                }
                else {
                    pts.append({ tMax, *maIt });
                    pts.append({ tMin, *miIt });
                }
            }
            return pts;
        };

    auto setupPlot = [this, globalOffset](
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

    setupPlot(ui->amp_ecg1_axis, m_ecgAmpSeries,
        calcPoints(m_ecg1, m_ecgSR), m_ecgCursorBar, COLOR_ECG1,
        "ECG1 amplitude");
    setupPlot(ui->amp_ppg_axis, m_ppgAmpSeries,
        calcPoints(m_ppg, m_ppgSR), m_ppgCursorBar, COLOR_PPG,
        "PPG amplitude");

    // --- RESP 8h overview (raw, min/max-decimated) ---
    if (m_respAmpSeries && !isMissingSignal(m_resp)) {
        const int viewportPx = ui->resp_axis->viewport()
            ? ui->resp_axis->viewport()->width() : 1200;
        setupPlot(ui->resp_axis, m_respAmpSeries,
            calcEnvelopePoints(m_resp, m_ecgSR, std::max(1, viewportPx)),
            m_respCursorBar, COLOR_RESP, "RESP");
    }

    // --- CVP 8h overview (raw, min/max-decimated) ---
    // Only owns sleep_or_cvp_axis when there are no sleep stages.
    // setupHypnogram() runs after us (from loadChunkFromFile) and will
    // wipe the chart if sleep stages are present, so we only lay down
    // the CVP overview here when that path won't overwrite us.
    bool sleepPresent = !m_sleepStages.isEmpty()
        && !(m_sleepStages.size() == 1 && m_sleepStages[0] == -1.0);
    if (!sleepPresent && !isMissingSignal(m_cvp) && ui->sleep_or_cvp_axis) {
        auto* chart = ui->sleep_or_cvp_axis->chart();

        // If we're coming from a previous file that had sleep stages,
        // the hypnogram stage series (and cursor) are still on this chart.
        // Strip them out before we take over.
        for (auto* s : m_hypnoStageSeries) { chart->removeSeries(s); delete s; }
        m_hypnoStageSeries.clear();
        if (m_hypnoCursorBar && m_hypnoCursorBar->chart() == chart) {
            chart->removeSeries(m_hypnoCursorBar);
            // don't delete: it's reusable; will re-add if a future file has sleep stages
        }

        // Lazily create the CVP overview series + cursor (they are NOT
        // pre-created in the ctor because this chart is multi-purpose).
        if (!m_cvpAmpSeries) {
            m_cvpAmpSeries = new QLineSeries();
            chart->addSeries(m_cvpAmpSeries);
        }
        else if (m_cvpAmpSeries->chart() != chart) {
            if (m_cvpAmpSeries->chart())
                m_cvpAmpSeries->chart()->removeSeries(m_cvpAmpSeries);
            chart->addSeries(m_cvpAmpSeries);
        }
        if (!m_cvpCursorBar) {
            m_cvpCursorBar = new QLineSeries();
            m_cvpCursorBar->setPen(QPen(Qt::black, 2));
            chart->addSeries(m_cvpCursorBar);
        }
        else if (m_cvpCursorBar->chart() != chart) {
            if (m_cvpCursorBar->chart())
                m_cvpCursorBar->chart()->removeSeries(m_cvpCursorBar);
            chart->addSeries(m_cvpCursorBar);
        }
        const int viewportPx = ui->sleep_or_cvp_axis->viewport()
            ? ui->sleep_or_cvp_axis->viewport()->width() : 1200;
        setupPlot(ui->sleep_or_cvp_axis, m_cvpAmpSeries,
            calcEnvelopePoints(m_cvp, m_ecgSR, std::max(1, viewportPx)),
            m_cvpCursorBar, COLOR_CVP, "CVP");
    }
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
    QChartView* xLabelOwner = nullptr;
    {
        const QList<QChartView*> order = {
            ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3,
            ui->ppg_axis, ui->accel_or_abg_axis
        };
        for (auto* cv : order) {
            if (cv && cv->isVisible()) xLabelOwner = cv;
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
    wipeChart(ui->accel_or_abg_axis, {});
    // Note: resp_axis and sleep_or_cvp_axis are NOT wiped here.
    //   - resp_axis is owned by handle_ampogram_plot (8h RESP overview).
    //   - sleep_or_cvp_axis is owned by either setupHypnogram() or the
    //     CVP branch of handle_ampogram_plot.
    // Wiping here would delete the series they manage.

    auto plotSignal = [&](QChartView* view, const QVector<double>& data, double sr,
        QLineSeries* marker, double markerPos, const QColor& color,
        const QString& signalName = QString())
        -> std::pair<double, double> {
        if (!view || !view->chart()) return { 1e9, -1e9 };
        QChart* chart = view->chart();

        for (auto* s : chart->series()) { if (s != marker) { chart->removeSeries(s); delete s; } }
        for (auto* a : chart->axes()) { chart->removeAxis(a); delete a; }

        // Show signal name as a small left-aligned title
        if (!signalName.isEmpty()) {
            chart->setTitle(signalName);
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
        if (view != xLabelOwner) xAxis->setLabelsVisible(false);
        chart->addAxis(xAxis, Qt::AlignBottom);
        chart->setMargins(QMargins(0, 0, 20, 0));

        auto* yAxis = new QValueAxis();
        yAxis->setVisible(false);
        chart->addAxis(yAxis, Qt::AlignLeft);

        if (data.size() < 2 || sr <= 0.0) {
            yAxis->setRange(-1.0, 1.0);
            return { 1e9, -1e9 };
        }

        QXYSeries* plotSeries = nullptr;
        if (m_plotMode == PlotMode::Scatter) {
            auto* sc = new QScatterSeries();
            sc->setColor(color);
            sc->setBorderColor(Qt::transparent);   // no border: keeps 2px points crisp and fast
            sc->setMarkerSize(2.0);
            sc->setMarkerShape(QScatterSeries::MarkerShapeCircle);
            sc->setUseOpenGL(true);                // scatter supports OpenGL since Qt 5.9
            chart->addSeries(sc);
            plotSeries = sc;
        }
        else {
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
        double lMin = 1e9, lMax = -1e9;
        for (int i = startIdx; i < endIdx; ++i) {
            pts.append({ static_cast<double>(i) / sr, data[i] });
            if (data[i] < lMin) lMin = data[i];
            if (data[i] > lMax) lMax = data[i];
        }
        plotSeries->replace(pts);
        plotSeries->attachAxis(xAxis);
        plotSeries->attachAxis(yAxis);

        yAxis->setRange(lMin - 0.5, lMax + 0.5);
        if (marker && marker->chart() == chart) {
            marker->replace({ {markerPos, yAxis->min()}, {markerPos, yAxis->max()} });
            marker->attachAxis(xAxis);
            marker->attachAxis(yAxis);
        }
        return { lMin, lMax };
        };

    auto maybePlot = [&](const QString& label, QChartView* view,
        const QVector<double>& data, double sr,
        ChannelMarkingState& state, const QColor& color) {
            if (!isChannelActive(label)) return;
            double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
            double localMarkerPos = state.globalStartTime - globalOffset;
            plotSignal(view, data, sr, state.startMarkerLine, localMarkerPos, color, label);
        };

    maybePlot("ECG1", ui->ecg_axis_1, m_ecg1, m_ecgSR, m_markState_ecg1, COLOR_ECG1);
    maybePlot("ECG2", ui->ecg_axis_2, m_ecg2, m_ecgSR, m_markState_ecg2, COLOR_ECG2);
    maybePlot("ECG3", ui->ecg_axis_3, m_ecg3, m_ecgSR, m_markState_ecg3, COLOR_ECG3);
    maybePlot("PPG", ui->ppg_axis, m_ppg, m_ppgSR, m_markState_ppg, COLOR_PPG);

    // ------------------------------------------------------------------
    // Display-only charts (no marking): helper to plot 1..N signals on
    // a chart with the shared time axis and auto-ranged Y axis.
    // ------------------------------------------------------------------
    struct DisplaySeries {
        const QVector<double>* data;
        QColor                 color;
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
            chart->setTitleBrush(titleColor);

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
            xAxis->setLabelsVisible(view == xLabelOwner);
            chart->addAxis(xAxis, Qt::AlignBottom);
            chart->setMargins(QMargins(0, 0, 20, 0));

            auto* yAxis = new QValueAxis();
            yAxis->setVisible(false);
            chart->addAxis(yAxis, Qt::AlignLeft);

            const bool useScatter = (m_plotMode == PlotMode::Scatter) && !forceLine;

            double gMin = 1e9, gMax = -1e9;
            for (const auto& d : serieses) {
                if (!d.data || isMissingSignal(*d.data)) continue;

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
            }

            if (gMin < gMax) yAxis->setRange(gMin - 0.5, gMax + 0.5);
            else             yAxis->setRange(-1.0, 1.0);
        };

    // --- accel_or_abg_axis: prefer accel; if absent, fall through to ABP ---
    if (ui->accel_or_abg_axis) {
        bool anyAccel = !isMissingSignal(m_accelX)
            || !isMissingSignal(m_accelY)
            || !isMissingSignal(m_accelZ);
        bool hasAbp = !isMissingSignal(m_abp);

        if (anyAccel) {
            plotDisplayChart(ui->accel_or_abg_axis, "ACCEL", COLOR_ACCEL_X, {
                { &m_accelX, COLOR_ACCEL_X },
                { &m_accelY, COLOR_ACCEL_Y },
                { &m_accelZ, COLOR_ACCEL_Z },
                });
        }
        else if (hasAbp) {
            // Route ABP through the markable path so it gets the start marker,
            // noise-highlight overlay, and click-to-mark behavior like ECG/PPG.
            maybePlot("ABP", ui->accel_or_abg_axis, m_abp, m_ecgSR,
                m_markState_abp, COLOR_ABP);
        }
    }

    // resp_axis and sleep_or_cvp_axis are 8-hour overview charts (like the
    // ampograms), owned by handle_ampogram_plot() and setupHypnogram(),
    // not by the 10-second-window plotter. Nothing to do here.

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
        for (const QString& lbl : { "ECG1", "ECG2", "ECG3", "PPG", "ABP" }) {
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
                for (const QString& lbl : { "ECG1", "ECG2", "ECG3", "PPG", "ABP" }) {
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
                    for (const QString& lbl : { "ECG1", "ECG2", "ECG3", "PPG", "ABP" }) {
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
                    for (const QString& lbl : { "ECG1", "ECG2", "ECG3", "PPG", "ABP" }) {
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
        // resp_axis is always a full-8h RESP overview when RESP is present.
        // sleep_or_cvp_axis is an 8h overview in BOTH modes now:
        //   - hypnogram when sleep stages are present
        //   - CVP overview when CVP is present and no sleep stages
        // In all these cases, a click maps global-chart-X to a chunk-local
        // time, which is what we want for scrubbing.
        bool sleepPresent = !m_sleepStages.isEmpty()
            && !(m_sleepStages.size() == 1 && m_sleepStages[0] == -1.0);
        bool respPresent = !isMissingSignal(m_resp);
        bool cvpPresent = !isMissingSignal(m_cvp);
        bool isNavChart =
            (cv == ui->amp_ecg1_axis)
            || (cv == ui->amp_ppg_axis)
            || (cv == ui->resp_axis && respPresent)
            || (cv == ui->sleep_or_cvp_axis && (sleepPresent || cvpPresent));
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

    draw(ui->amp_ecg1_axis, m_ecgCursorBar);
    draw(ui->amp_ppg_axis, m_ppgCursorBar);

    // RESP overview (always full 8h when present).
    if (!isMissingSignal(m_resp))
        draw(ui->resp_axis, m_respCursorBar);

    // sleep_or_cvp_axis is multi-purpose. Use the hypno cursor when
    // sleep stages are showing; otherwise use the CVP cursor.
    bool sleepPresent = !m_sleepStages.isEmpty()
        && !(m_sleepStages.size() == 1 && m_sleepStages[0] == -1.0);
    if (sleepPresent)
        draw(ui->sleep_or_cvp_axis, m_hypnoCursorBar);
    else if (!isMissingSignal(m_cvp))
        draw(ui->sleep_or_cvp_axis, m_cvpCursorBar);
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

    for (const QString& lbl : { "ECG1", "ECG2", "ECG3", "PPG", "ABP" }) {
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

    for (const QString& label : { "ECG1", "ECG2", "ECG3", "PPG", "ABP" }) {
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
    for (const QString& label : { "ECG1", "ECG2", "ECG3", "PPG", "ABP" }) {
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