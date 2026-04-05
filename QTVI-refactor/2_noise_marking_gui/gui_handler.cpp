/**
 * @file   gui_handler.cpp
 * @brief  Entry point of the noise marking program. Makes a GUI where you can
 *         mark the noise you see in the ECG1/ECG2/ECG3/PPG, and outputs a csv
 *         and a .bin file.
 *
 * Input .bin format (88-byte header):
 *   Offset  0: ecgRate       (double)  — ECG sampling rate in Hz
 *   Offset  8: ppgRate       (double)  — PPG sampling rate in Hz
 *   Offset 16: epochSize     (double)  — sleep stage epoch duration in seconds
 *   Offset 24: size1–sizeS   (uint64)  — sample counts per channel
 *
 * Signal data (contiguous doubles after header):
 *   [size1] ECG1, [size2] ECG2, [size3] ECG3, [sizeP] PPG, [sizeS] Sleep
 *
 * Missing channels are stored as a single -1.0 with their size field set to 1.
 *
 * Output .bin format:
 *   [uint64] count, then count × 6 doubles:
 *     startSample, endSample, startSec, endSec, labelId, typeId
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
#include <QtCharts/QValueAxis>
#include <QtCharts/QCategoryAxis>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QShortcut>
#include <QFile>
#include <algorithm>

 // ============================================================================
 // Constants
 // ============================================================================

static const QColor COLOR_ECG1 = QColor("#1ABC9C");
static const QColor COLOR_ECG2 = QColor("#3498DB");
static const QColor COLOR_ECG3 = QColor("#9B59B6");
static const QColor COLOR_PPG = QColor("#E74C3C");

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
        for (auto* axis : chart->axes()) chart->removeAxis(axis);
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
    return m_markState_ppg;
}

QString noise_marking_gui::signalLabelForChartView(QChartView* cv) const {
    if (cv == ui->ecg_axis_1) return "ECG1";
    if (cv == ui->ecg_axis_2) return "ECG2";
    if (cv == ui->ecg_axis_3) return "ECG3";
    if (cv == ui->ppg_axis)   return "PPG";
    return {};
}

QChartView* noise_marking_gui::chartViewForSignalLabel(const QString& label) const {
    if (label == "ECG1") return ui->ecg_axis_1;
    if (label == "ECG2") return ui->ecg_axis_2;
    if (label == "ECG3") return ui->ecg_axis_3;
    if (label == "PPG")  return ui->ppg_axis;
    return nullptr;
}

double noise_marking_gui::sampleRateForSignal(const QString& label) const {
    return (label == "PPG") ? m_ppgSR : m_ecgSR;
}

QColor noise_marking_gui::colorForSignal(const QString& label) const {
    if (label == "ECG1") return COLOR_ECG1;
    if (label == "ECG2") return COLOR_ECG2;
    if (label == "ECG3") return COLOR_ECG3;
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
    return nullptr;
}

QPushButton* noise_marking_gui::stopButtonForSignal(const QString& label) const {
    if (label == "ECG1") return ui->stop_ecg1_mark;
    if (label == "ECG2") return ui->stop_ecg2_mark;
    if (label == "ECG3") return ui->stop_ecg3_mark;
    if (label == "PPG")  return ui->stopNoisePPG;
    return nullptr;
}

void noise_marking_gui::updateButtonStatesForChannel(const QString& label) {
    QPushButton* startBtn = startButtonForSignal(label);
    QPushButton* stopBtn = stopButtonForSignal(label);
    if (!startBtn || !stopBtn) return;

    bool active = isChannelActive(label);

    // Gray out both buttons if channel has no data
    startBtn->setEnabled(active);
    stopBtn->setEnabled(false);

    // Clear styling
    startBtn->setStyleSheet("");
    stopBtn->setStyleSheet("");

    if (!active) {
        startBtn->setStyleSheet("color: gray;");
        stopBtn->setStyleSheet("color: gray;");
    }
}

void noise_marking_gui::updateAllChannelButtonStates() {
    for (const QString& label : { "ECG1", "ECG2", "ECG3", "PPG" })
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
    for (const QString& lbl : { "ECG1", "ECG2", "ECG3", "PPG" }) {
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
        ui->sleep_state_axis, ui->amp_ecg_axis, ui->amp_ppg_axis,
        ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3, ui->ppg_axis
    };
    for (auto* view : charts) {
        view->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        view->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        view->setRenderHint(QPainter::Antialiasing);
        view->setFocusPolicy(Qt::NoFocus);
        view->viewport()->installEventFilter(this);
    }

    ui->rb_10s->setChecked(true);
    m_skipInterval = ui->skip_interval_box->text().toDouble();
    if (m_skipInterval <= 0.0) m_skipInterval = 5.0;

    for (auto* cv : { ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3,
                     ui->ppg_axis, ui->amp_ecg_axis, ui->amp_ppg_axis })
        setupChartDefaults(cv);

    // Ampogram series
    m_ecgAmpSeries = new QLineSeries();
    m_ppgAmpSeries = new QLineSeries();
    ui->amp_ecg_axis->chart()->addSeries(m_ecgAmpSeries);
    ui->amp_ppg_axis->chart()->addSeries(m_ppgAmpSeries);

    // Initially disable all stop buttons
    ui->stop_ecg1_mark->setEnabled(false);
    ui->stop_ecg2_mark->setEnabled(false);
    ui->stop_ecg3_mark->setEnabled(false);
    ui->stopNoisePPG->setEnabled(false);
    ui->stop_all_mark->setEnabled(false);

    // Cursor bars
    auto addCursor = [](QChartView* view, QLineSeries*& series) {
        series = new QLineSeries();
        series->setPen(QPen(Qt::black, 2));
        view->chart()->addSeries(series);
        };
    addCursor(ui->amp_ecg_axis, m_ecgCursorBar);
    addCursor(ui->amp_ppg_axis, m_ppgCursorBar);

    // Hypnogram chart
    auto* hypnoChart = new QChart();
    hypnoChart->legend()->hide();
    hypnoChart->setMargins(QMargins(0, 0, 0, 0));
    ui->sleep_state_axis->setChart(hypnoChart);

    m_hypnoCursorBar = new QLineSeries();
    m_hypnoCursorBar->setPen(QPen(Qt::black, 2));
    hypnoChart->addSeries(m_hypnoCursorBar);

    m_currentMarkingType = ui->marking_type->currentText();
}

noise_marking_gui::~noise_marking_gui() = default;

// ============================================================================
// File Loading
// ============================================================================

void noise_marking_gui::setFileSource(const QString& filePath) {
    m_binFilePath = filePath;
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly)) return;

    double ecgRate = 0, ppgRate = 0, epochSize = 0;
    file.read(reinterpret_cast<char*>(&ecgRate), sizeof(double));
    file.read(reinterpret_cast<char*>(&ppgRate), sizeof(double));
    file.read(reinterpret_cast<char*>(&epochSize), sizeof(double));

    m_ecgSR = ecgRate;
    m_ppgSR = ppgRate;
    m_sleepSR = (epochSize > 0) ? (1.0 / epochSize) : 0;

    file.read(reinterpret_cast<char*>(&m_totalEcg1Samples), sizeof(uint64_t));
    file.read(reinterpret_cast<char*>(&m_totalEcg2Samples), sizeof(uint64_t));
    file.read(reinterpret_cast<char*>(&m_totalEcg3Samples), sizeof(uint64_t));
    file.read(reinterpret_cast<char*>(&m_totalPpgSamples), sizeof(uint64_t));
    file.read(reinterpret_cast<char*>(&m_totalSleepSamples), sizeof(uint64_t));

    m_fileHeaderSize = file.pos();
    file.close();

    m_noiseManager = std::make_unique<NoiseManager>(m_ecgSR);
    loadChunkFromFile(0);
}

bool noise_marking_gui::loadChunkFromFile(uint64_t chunkIndex) {
    QFile file(m_binFilePath);
    if (!file.open(QIODevice::ReadOnly)) return false;

    m_currentChunkIndex = chunkIndex;

    auto loadSignal = [&](QVector<double>& dest, uint64_t totalSamples,
        double sr, uint64_t sampleOffset) {
            uint64_t perChunk = static_cast<uint64_t>(CHUNK_DURATION_SEC * sr);
            uint64_t start = chunkIndex * perChunk;
            uint64_t count = (totalSamples > start)
                ? std::min(perChunk, totalSamples - start) : 0;
            dest.resize(static_cast<int>(count));
            file.seek(m_fileHeaderSize + (sampleOffset + start) * sizeof(double));
            file.read(reinterpret_cast<char*>(dest.data()), count * sizeof(double));
        };

    uint64_t offset = 0;
    loadSignal(m_ecg1, m_totalEcg1Samples, m_ecgSR, offset);  offset += m_totalEcg1Samples;
    loadSignal(m_ecg2, m_totalEcg2Samples, m_ecgSR, offset);  offset += m_totalEcg2Samples;
    loadSignal(m_ecg3, m_totalEcg3Samples, m_ecgSR, offset);  offset += m_totalEcg3Samples;
    loadSignal(m_ppg, m_totalPpgSamples, m_ppgSR, offset);  offset += m_totalPpgSamples;
    loadSignal(m_sleepStages, m_totalSleepSamples, m_sleepSR, offset);
    file.close();

    // Header label
    int startHr = chunkIndex * 8;
    ui->topLabel->setText(
        QString("     Data Range: Hour %1 to Hour %2").arg(startHr).arg(startHr + 8));

    m_currentStartTime = 0;

    // Determine which channels have real data
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

    // Gray out / enable buttons based on which channels are present
    updateAllChannelButtonStates();

    handle_ampogram_plot();
    handle_data_plot();
    updateAmpogramCursor();
    setupHypnogram();
    restoreMarkingMarkers();

    uint64_t ecgPerChunk = static_cast<uint64_t>(CHUNK_DURATION_SEC * m_ecgSR);
    ui->prev8hours->setEnabled(chunkIndex > 0);
    ui->next8hours->setEnabled((chunkIndex * ecgPerChunk + m_ecg1.size()) < m_totalEcg1Samples);

    return true;
}

void noise_marking_gui::on_next8hours_clicked() { loadChunkFromFile(m_currentChunkIndex + 1); }
void noise_marking_gui::on_prev8hours_clicked() { if (m_currentChunkIndex > 0) loadChunkFromFile(m_currentChunkIndex - 1); }

// ============================================================================
// Hypnogram
// ============================================================================

void noise_marking_gui::setupHypnogram() {
    if (m_sleepSR <= 0.0) return;
    auto* chart = ui->sleep_state_axis->chart();

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

    auto setupPlot = [this, globalOffset](
        QChartView* view, QLineSeries* series, const QList<QPointF>& pts,
        QLineSeries* cursor, const QColor& color) {
            series->replace(pts);
            series->setPen(QPen(color, 1));

            auto* chart = view->chart();
            clearAxes(chart);
            chart->legend()->hide();

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

    setupPlot(ui->amp_ecg_axis, m_ecgAmpSeries,
        calcPoints(m_ecg1, m_ecgSR), m_ecgCursorBar, COLOR_ECG1);
    setupPlot(ui->amp_ppg_axis, m_ppgAmpSeries,
        calcPoints(m_ppg, m_ppgSR), m_ppgCursorBar, COLOR_PPG);
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

    auto plotSignal = [&](QChartView* view, const QVector<double>& data, double sr,
        QLineSeries* marker, double markerPos, const QColor& color)
        -> std::pair<double, double> {
        if (!view || !view->chart()) return { 1e9, -1e9 };
        QChart* chart = view->chart();

        for (auto* s : chart->series()) { if (s != marker) { chart->removeSeries(s); delete s; } }
        for (auto* a : chart->axes()) { chart->removeAxis(a); delete a; }

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
        if (view != ui->ppg_axis) xAxis->setLabelsVisible(false);
        chart->addAxis(xAxis, Qt::AlignBottom);
        chart->setMargins(QMargins(0, 0, 20, 0));

        auto* yAxis = new QValueAxis();
        yAxis->setVisible(false);
        chart->addAxis(yAxis, Qt::AlignLeft);

        if (data.size() < 2 || sr <= 0.0) {
            yAxis->setRange(-1.0, 1.0);
            return { 1e9, -1e9 };
        }

        auto* series = new QLineSeries();
        series->setUseOpenGL(true);
        series->setPen(QPen(color, 1));
        chart->addSeries(series);

        int startIdx = std::clamp(static_cast<int>(m_currentStartTime * sr),
            0, static_cast<int>(data.size() - 1));
        int endIdx = std::clamp(static_cast<int>((m_currentStartTime + m_windowDuration) * sr),
            0, static_cast<int>(data.size()));

        QList<QPointF> pts;
        double lMin = 1e9, lMax = -1e9;
        for (int i = startIdx; i < endIdx; ++i) {
            pts.append({ static_cast<double>(i) / sr, data[i] });
            lMin = std::min(lMin, data[i]);
            lMax = std::max(lMax, data[i]);
        }
        series->replace(pts);
        series->attachAxis(xAxis);
        series->attachAxis(yAxis);

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
            plotSignal(view, data, sr, state.startMarkerLine, localMarkerPos, color);
        };

    maybePlot("ECG1", ui->ecg_axis_1, m_ecg1, m_ecgSR, m_markState_ecg1, COLOR_ECG1);
    maybePlot("ECG2", ui->ecg_axis_2, m_ecg2, m_ecgSR, m_markState_ecg2, COLOR_ECG2);
    maybePlot("ECG3", ui->ecg_axis_3, m_ecg3, m_ecgSR, m_markState_ecg3, COLOR_ECG3);
    maybePlot("PPG", ui->ppg_axis, m_ppg, m_ppgSR, m_markState_ppg, COLOR_PPG);

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
        for (const QString& lbl : { "ECG1", "ECG2", "ECG3", "PPG" }) {
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
                // Check if any active channel is in a waiting phase
                bool anyWaitingStart = false;
                bool anyWaitingStop = false;
                for (const QString& lbl : { "ECG1", "ECG2", "ECG3", "PPG" }) {
                    if (!isChannelActive(lbl)) continue;
                    MarkPhase p = markStateFor(lbl).phase;
                    if (p == MarkPhase::WaitingForStart || p == MarkPhase::WaitingForEnd)
                        anyWaitingStart = true;
                    if (p == MarkPhase::WaitingForStop)
                        anyWaitingStop = true;
                }

                if (anyWaitingStart) {
                    // Place start marker on all active channels at this time
                    for (const QString& lbl : { "ECG1", "ECG2", "ECG3", "PPG" }) {
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

                if (anyWaitingStop) {
                    // Finalize all active channels at this time
                    for (const QString& lbl : { "ECG1", "ECG2", "ECG3", "PPG" }) {
                        if (!isChannelActive(lbl)) continue;
                        ChannelMarkingState& st = markStateFor(lbl);
                        if (st.phase != MarkPhase::WaitingForStop) continue;
                        QChartView* targetCv = chartViewForSignalLabel(lbl);
                        finalizeMarking(targetCv, clickedX, lbl);
                    }
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

        // --- Ampogram / hypnogram click (navigate) ---
        if (cv == ui->amp_ecg_axis || cv == ui->amp_ppg_axis
            || cv == ui->sleep_state_axis) {
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

    draw(ui->amp_ecg_axis, m_ecgCursorBar);
    draw(ui->amp_ppg_axis, m_ppgCursorBar);
    draw(ui->sleep_state_axis, m_hypnoCursorBar);
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

    for (const QString& lbl : { "ECG1", "ECG2", "ECG3", "PPG" }) {
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

    for (const QString& label : { "ECG1", "ECG2", "ECG3", "PPG" }) {
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
    for (const QString& label : { "ECG1", "ECG2", "ECG3", "PPG" }) {
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