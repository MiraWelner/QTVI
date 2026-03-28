/**
 * @file   gui_handler.cpp
 * @brief  Entry point of the noise marking program. Makes a GUI where you can mark the noise you see in the
 *         ECG1/ECG2/ECG3/PPG, and outputs a csv and a .bin file
 *
 * Input .bin format (88-byte header):
 *   Offset  0: ecgRate       (double)  — ECG sampling rate in Hz (always 2000.0)
 *   Offset  8: ppgRate       (double)  — PPG sampling rate in Hz (always 2000.0)
 *   Offset 16: epochSize     (double)  — sleep stage epoch duration in seconds
 *   Offset 24: size1         (uint64)  — number of samples in ECG channel 1
 *   Offset 32: size2         (uint64)  — number of samples in ECG channel 2
 *   Offset 40: size3         (uint64)  — number of samples in ECG channel 3
 *   Offset 48: sizeP         (uint64)  — number of samples in PPG channel
 *   Offset 56: sizeS         (uint64)  — number of sleep stage values
 *
 * Signal data (contiguous doubles after header):
 *   [size1]  ECG channel 1
 *   [size2]  ECG channel 2
 *   [size3]  ECG channel 3
 *   [sizeP]  PPG
 *   [sizeS]  Sleep stages
 *
 * Missing channels are stored as a single -1.0 with their size field set to 1.
 *
 * Output .bin format:
 *   Header:
 *     [uint64]  count — number of annotation segments
 *
 *   Followed by count rows of 6 doubles each (48 bytes per row):
 *     [0] startSample   — first sample index of the annotation
 *     [1] endSample     — last sample index of the annotation
 *     [2] startSec      — start time in seconds (startSample / sampleRate)
 *     [3] endSec        — end time in seconds (endSample / sampleRate)
 *     [4] labelId       — signal type: 0=unknown, 1=PPG, 2=ECG1, 3=ECG2, 4=ECG3
 *     [5] typeId        — marking type:
 *                            0=unknown, 1=Noise/Artifact, 2=AF, 3=SVT, 4=VT,
 *                            5=PVC, 6=PAC, 7=Benign Arrhythmia, 8=Significant Arrhythmia
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


static const QColor COLOR_ECG1 = QColor("#1ABC9C");
static const QColor COLOR_ECG2 = QColor("#3498DB");
static const QColor COLOR_ECG3 = QColor("#9B59B6");
static const QColor COLOR_PPG = QColor("#E74C3C");

// Marking type -> highlight color for annotation overlays
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
}

// ============================================================================
// Channel state / mapping helpers
// ============================================================================

noise_marking_gui::ChannelMarkingState& noise_marking_gui::markStateFor(const QString& label) {
    if (label == "ECG1") return m_markState_ecg1;
    if (label == "ECG2") return m_markState_ecg2;
    if (label == "ECG3") return m_markState_ecg3;
    return m_markState_ppg;
}

QString noise_marking_gui::selectedEcgLabel() const {
    // ecg_channel_selector items: "ECG 1", "ECG 2", "ECG 3"
    int idx = ui->ecg_channel_selector->currentIndex();
    if (idx == 1) return "ECG2";
    if (idx == 2) return "ECG3";
    return "ECG1";
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
    if (label == "PPG") return m_ppgSR;
    return m_ecgSR;  // ECG1, ECG2, ECG3 all share the same rate
}

QColor noise_marking_gui::colorForSignal(const QString& label) const {
    if (label == "ECG1") return COLOR_ECG1;
    if (label == "ECG2") return COLOR_ECG2;
    if (label == "ECG3") return COLOR_ECG3;
    return COLOR_PPG;
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

    // Size to 75% x 90% of screen
    if (auto* screen = QGuiApplication::primaryScreen()) {
        const QRect avail = screen->availableGeometry();
        resize(avail.width() * 3 / 4, static_cast<int>(avail.height() * 0.9));
        move(avail.center() - rect().center());
    }

    m_buttonHandler->setupConnections();

    // Prevent buttons from stealing keyboard focus (arrow keys navigate)
    for (auto* btn : findChildren<QPushButton*>()) {
        btn->setFocusPolicy(Qt::NoFocus);
    }
    ui->ecg_channel_selector->setFocusPolicy(Qt::NoFocus);

    // Populate the ECG channel selector if not already done in Designer
    if (ui->ecg_channel_selector->count() == 0) {
        ui->ecg_channel_selector->addItems({ "ECG 1", "ECG 2", "ECG 3" });
    }

    // --- Navigation shortcuts ---
    new QShortcut(QKeySequence(Qt::Key_Left), this, [this]() {
        m_currentStartTime = std::max(0.0, m_currentStartTime - m_skipInterval);
        handle_data_plot();
        updateAmpogramCursor();
        });

    new QShortcut(QKeySequence(Qt::Key_Right), this, [this]() {
        double chunkDur = 0.0;
        if (m_ecg1.size() > 1 && m_ecgSR > 0) chunkDur = m_ecg1.size() / m_ecgSR;
        else if (m_ppg.size() > 1 && m_ppgSR > 0) chunkDur = m_ppg.size() / m_ppgSR;

        m_currentStartTime = std::min(m_currentStartTime + m_skipInterval,
            std::max(0.0, chunkDur - m_windowDuration));
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

    setupChartDefaults(ui->ecg_axis_1);
    setupChartDefaults(ui->ecg_axis_2);
    setupChartDefaults(ui->ecg_axis_3);
    setupChartDefaults(ui->ppg_axis);
    setupChartDefaults(ui->amp_ecg_axis);
    setupChartDefaults(ui->amp_ppg_axis);

    // Ampogram series
    m_ecgAmpSeries = new QLineSeries();
    m_ppgAmpSeries = new QLineSeries();
    ui->amp_ecg_axis->chart()->addSeries(m_ecgAmpSeries);
    ui->amp_ppg_axis->chart()->addSeries(m_ppgAmpSeries);

    // Disable stop buttons initially
    ui->stopNoiseECG->setEnabled(false);
    ui->stopNoisePPG->setEnabled(false);

    // Cursor bars for ampogram and hypnogram
    auto setupCursor = [](QChartView* view, QLineSeries*& series) {
        series = new QLineSeries();
        series->setPen(QPen(Qt::black, 2));
        view->chart()->addSeries(series);
        };
    setupCursor(ui->amp_ecg_axis, m_ecgCursorBar);
    setupCursor(ui->amp_ppg_axis, m_ppgCursorBar);

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
// File Loading (64-byte header)
// ============================================================================

void noise_marking_gui::setFileSource(const QString& filePath) {
    m_binFilePath = filePath;
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly)) return;

    double ecgRate = 0, ppgRate = 0, epochSize = 0;
    file.read(reinterpret_cast<char*>(&ecgRate), sizeof(double));  // offset 0
    file.read(reinterpret_cast<char*>(&ppgRate), sizeof(double));  // offset 8
    file.read(reinterpret_cast<char*>(&epochSize), sizeof(double));  // offset 16

    m_ecgSR = ecgRate;
    m_ppgSR = ppgRate;
    m_sleepSR = (epochSize > 0) ? (1.0 / epochSize) : 0;

    file.read(reinterpret_cast<char*>(&m_totalEcg1Samples), sizeof(uint64_t)); // offset 24
    file.read(reinterpret_cast<char*>(&m_totalEcg2Samples), sizeof(uint64_t)); // offset 32
    file.read(reinterpret_cast<char*>(&m_totalEcg3Samples), sizeof(uint64_t)); // offset 40
    file.read(reinterpret_cast<char*>(&m_totalPpgSamples), sizeof(uint64_t)); // offset 48
    file.read(reinterpret_cast<char*>(&m_totalSleepSamples), sizeof(uint64_t)); // offset 56

    m_fileHeaderSize = file.pos(); // should be 64

    file.close();

    m_noiseManager = std::make_unique<NoiseManager>(m_ecgSR);
    loadChunkFromFile(0);
}

bool noise_marking_gui::loadChunkFromFile(uint64_t chunkIndex) {
    QFile file(m_binFilePath);
    if (!file.open(QIODevice::ReadOnly)) return false;

    m_currentChunkIndex = chunkIndex;

    auto loadSignal = [&](QVector<double>& dest, uint64_t totalSamples, double sr, uint64_t sampleOffset) {
        uint64_t samplesPerChunk = static_cast<uint64_t>(CHUNK_DURATION_SEC * sr);
        uint64_t start = chunkIndex * samplesPerChunk;
        uint64_t count = (totalSamples > start) ? std::min(samplesPerChunk, totalSamples - start) : 0;
        dest.resize(static_cast<int>(count));
        file.seek(m_fileHeaderSize + (sampleOffset + start) * sizeof(double));
        file.read(reinterpret_cast<char*>(dest.data()), count * sizeof(double));
        };

    // Signal data layout after 88-byte header:
    //   ECG1, ECG2, ECG3, PPG, Sleep, Abs1, Abs2, Abs3
    // We load the first 5, skip the abs channels.
    uint64_t offset = 0;
    loadSignal(m_ecg1, m_totalEcg1Samples, m_ecgSR, offset);  offset += m_totalEcg1Samples;
    loadSignal(m_ecg2, m_totalEcg2Samples, m_ecgSR, offset);  offset += m_totalEcg2Samples;
    loadSignal(m_ecg3, m_totalEcg3Samples, m_ecgSR, offset);  offset += m_totalEcg3Samples;
    loadSignal(m_ppg, m_totalPpgSamples, m_ppgSR, offset);  offset += m_totalPpgSamples;
    loadSignal(m_sleepStages, m_totalSleepSamples, m_sleepSR, offset);
    // abs channels intentionally not loaded

    file.close();

    int startHr = chunkIndex * 8;
    ui->topLabel->setText(QString("     Data Range: Hour %1 to Hour %2").arg(startHr).arg(startHr + 8));

    m_currentStartTime = 0;
    handle_ampogram_plot();
    handle_data_plot();
    updateAmpogramCursor();
    setupHypnogram();

    uint64_t ecgSamplesPerChunk = static_cast<uint64_t>(CHUNK_DURATION_SEC * m_ecgSR);
    ui->prev8hours->setEnabled(chunkIndex > 0);
    ui->next8hours->setEnabled((chunkIndex * ecgSamplesPerChunk + m_ecg1.size()) < m_totalEcg1Samples);

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
    const QList<Stage> stages = { {0, Qt::black}, {1, Qt::darkGreen}, {2, Qt::blue}, {3, Qt::cyan}, {4, Qt::red} };

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
                if (static_cast<int>(m_sleepStages[i]) == st.value) {
                    s->append(globalOffset + i * dt + dt / 2.0, st.value);
                }
            }
            chart->addSeries(s);
            m_hypnoStageSeries.append(s);
        }
    }

    clearAxes(chart);

    auto* xAxis = new QCategoryAxis();
    xAxis->setRange(globalOffset, globalOffset + CHUNK_DURATION_SEC);
    for (int h = 0; h <= 8; ++h) xAxis->append(QString::number(h), globalOffset + h * 3600.0);
    xAxis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
    xAxis->setGridLineVisible(false);
    xAxis->setLabelsFont(QFont("Arial", 6));

    auto* yAxis = new QCategoryAxis();
    for (const auto& st : stages) yAxis->append("", st.value + 0.4);
    yAxis->setRange(-0.5, 4.5);
    yAxis->setReverse(true);
    yAxis->setVisible(false);
    yAxis->setGridLineVisible(false);

    chart->addAxis(xAxis, Qt::AlignBottom);
    chart->addAxis(yAxis, Qt::AlignLeft);
    for (auto* s : chart->series()) { s->attachAxis(xAxis); s->attachAxis(yAxis); }

    if (m_hypnoCursorBar) {
        chart->removeSeries(m_hypnoCursorBar);
        chart->addSeries(m_hypnoCursorBar);
        m_hypnoCursorBar->attachAxis(xAxis);
        m_hypnoCursorBar->attachAxis(yAxis);
    }
}

// ============================================================================
// Ampogram (uses ECG1 for the ECG ampogram)
// ============================================================================

void noise_marking_gui::handle_ampogram_plot(double sampling_length) {
    double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;

    auto calcPoints = [sampling_length, globalOffset](const QVector<double>& data, double sr) {
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

    auto setupPlot = [this, globalOffset](QChartView* view, QLineSeries* series,
        const QList<QPointF>& pts, QLineSeries* cursor, const QColor& color) {
            series->replace(pts);
            series->setPen(QPen(color, 1));

            auto* chart = view->chart();
            clearAxes(chart);
            chart->legend()->hide();

            auto* xAxis = new QCategoryAxis();
            xAxis->setRange(globalOffset, globalOffset + CHUNK_DURATION_SEC);
            for (int h = 0; h <= 8; ++h) xAxis->append(QString::number(h), globalOffset + h * 3600.0);
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

    setupPlot(ui->amp_ecg_axis, m_ecgAmpSeries, calcPoints(m_ecg1, m_ecgSR), m_ecgCursorBar, COLOR_ECG1);
    setupPlot(ui->amp_ppg_axis, m_ppgAmpSeries, calcPoints(m_ppg, m_ppgSR), m_ppgCursorBar, COLOR_PPG);
}

// ============================================================================
// Main signal plot
// ============================================================================

void noise_marking_gui::handle_data_plot() {
    m_highlights.clear();

    auto plotSignal = [&](QChartView* view, const QVector<double>& data, double sr,
        QLineSeries* marker, double markerPos, const QColor& color, bool resetAxes) -> std::pair<double, double> {
            if (!view || !view->chart()) return { 1e9, -1e9 };
            QChart* chart = view->chart();

            if (resetAxes) {
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
                    xAxis->append(QString("%1:%2:%3").arg(h, 2, 10, QChar('0'))
                        .arg(m, 2, 10, QChar('0')).arg(s, 5, 'f', 2, QChar('0')), val);
                }
                xAxis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
                xAxis->setGridLineVisible(false);
                xAxis->setLabelsFont(QFont("Arial", 7));
                // Only show time labels on the bottom-most signal chart (PPG)
                if (view != ui->ppg_axis) xAxis->setLabelsVisible(false);
                chart->addAxis(xAxis, Qt::AlignBottom);
                chart->setMargins(QMargins(0, 0, 20, 0));

                auto* yAxis = new QValueAxis();
                yAxis->setVisible(false);
                chart->addAxis(yAxis, Qt::AlignLeft);
            }

            auto hAxes = chart->axes(Qt::Horizontal);
            auto vAxes = chart->axes(Qt::Vertical);
            if (hAxes.isEmpty() || vAxes.isEmpty()) return { 1e9, -1e9 };
            if (data.size() < 2 || sr <= 0.0) return { 1e9, -1e9 };

            auto* series = new QLineSeries();
            series->setUseOpenGL(true);
            series->setPen(QPen(color, 1));
            chart->addSeries(series);

            int startIdx = std::clamp(static_cast<int>(m_currentStartTime * sr), 0, static_cast<int>(data.size() - 1));
            int endIdx = std::clamp(static_cast<int>((m_currentStartTime + m_windowDuration) * sr), 0, static_cast<int>(data.size()));

            QList<QPointF> pts;
            double lMin = 1e9, lMax = -1e9;
            for (int i = startIdx; i < endIdx; ++i) {
                pts.append({ static_cast<double>(i) / sr, data[i] });
                if (data[i] < lMin) lMin = data[i];
                if (data[i] > lMax) lMax = data[i];
            }
            series->replace(pts);
            series->attachAxis(hAxes.first());
            series->attachAxis(vAxes.first());

            if (resetAxes) {
                auto* yAxis = qobject_cast<QValueAxis*>(vAxes.first());
                yAxis->setRange(lMin - 0.5, lMax + 0.5);
                if (marker && marker->chart() == chart) {
                    marker->replace({ {markerPos, yAxis->min()}, {markerPos, yAxis->max()} });
                    marker->attachAxis(hAxes.first());
                    marker->attachAxis(yAxis);
                }
            }
            return { lMin, lMax };
        };

    // Each ECG channel on its own chart view
    plotSignal(ui->ecg_axis_1, m_ecg1, m_ecgSR, m_markState_ecg1.startMarkerLine, m_markState_ecg1.startTimeValue, COLOR_ECG1, true);
    plotSignal(ui->ecg_axis_2, m_ecg2, m_ecgSR, m_markState_ecg2.startMarkerLine, m_markState_ecg2.startTimeValue, COLOR_ECG2, true);
    plotSignal(ui->ecg_axis_3, m_ecg3, m_ecgSR, m_markState_ecg3.startMarkerLine, m_markState_ecg3.startTimeValue, COLOR_ECG3, true);
    plotSignal(ui->ppg_axis, m_ppg, m_ppgSR, m_markState_ppg.startMarkerLine, m_markState_ppg.startTimeValue, COLOR_PPG, true);

    updateNoiseHighlights();
}

// ============================================================================
// Marking finalization
// ============================================================================

void noise_marking_gui::finalizeMarking(QChartView* cv, double endX, const QString& signalLabel) {
    ChannelMarkingState& state = markStateFor(signalLabel);
    double sr = sampleRateForSignal(signalLabel);

    double startX = state.startTimeValue;

    // Snap to nearest sample
    double snappedS = std::round(std::min(startX, endX) * sr) / sr;
    double snappedE = std::round(std::max(startX, endX) * sr) / sr;

    // Convert to global time
    double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
    double globalStart = snappedS + globalOffset;
    double globalEnd = snappedE + globalOffset;

    m_noiseManager->addSegment(
        static_cast<size_t>(globalStart * sr),
        static_cast<size_t>(globalEnd * sr),
        signalLabel.toStdString(),
        m_currentMarkingType.toStdString()
    );

    m_genExc.noiseExc.append({ globalStart, globalEnd });
    m_genExc.data_type.append(signalLabel);
    m_genExc.marking_type.append(m_currentMarkingType);

    state.isWaitingForStart = false;
    state.isWaitingForEnd = false;
    clearStartMarker(state);

    // Reset button styles for the appropriate pair
    if (signalLabel.startsWith("ECG")) {
        ui->startNoiseECG->setStyleSheet("");
        ui->stopNoiseECG->setStyleSheet("");
        ui->stopNoiseECG->setEnabled(false);
    }
    else {
        ui->startNoisePPG->setStyleSheet("");
        ui->stopNoisePPG->setStyleSheet("");
        ui->stopNoisePPG->setEnabled(false);
    }

    updateNoiseHighlights();
}

bool noise_marking_gui::eventFilter(QObject* watched, QEvent* event) {
    auto* viewport = qobject_cast<QWidget*>(watched);
    if (!viewport) return QDialog::eventFilter(watched, event);

    auto* cv = qobject_cast<QChartView*>(viewport->parent());
    if (!cv || !cv->chart()) return QDialog::eventFilter(watched, event);

    if (event->type() == QEvent::MouseButtonPress) {
        auto* me = static_cast<QMouseEvent*>(event);
        if (me->button() != Qt::LeftButton) return QDialog::eventFilter(watched, event);

        double clickedX = cv->chart()->mapToValue(me->pos()).x();
        double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
        double chunkDur = 0.0;
        if (m_ecg1.size() > 1 && m_ecgSR > 0)     chunkDur = m_ecg1.size() / m_ecgSR;
        else if (m_ppg.size() > 1 && m_ppgSR > 0)  chunkDur = m_ppg.size() / m_ppgSR;

        clickedX = std::clamp(clickedX, globalOffset, globalOffset + chunkDur);

        // --- Signal chart click (ECG1, ECG2, ECG3, or PPG) ---
        QString label = signalLabelForChartView(cv);
        if (!label.isEmpty()) {
            ChannelMarkingState& state = markStateFor(label);

            // Button-driven marking workflow
            if (state.isWaitingForStart) {
                state.startTimeValue = clickedX;
                QPushButton* stopBtn = label.startsWith("ECG") ? ui->stopNoiseECG : ui->stopNoisePPG;
                showStartMarker(cv, clickedX, state, colorForSignal(label), stopBtn);
                state.isWaitingForStart = false;
                state.isWaitingForEnd = true;
                return true;
            }
            if (state.isWaitingForEnd) {
                finalizeMarking(cv, clickedX, label);
                return true;
            }

            // Drag-based marking
            if (!m_isDragging) {
                m_isDragging = true;
                m_dragStartPos = me->pos();
                m_dragSignalLabel = label;
                if (!m_draggedViewport) { m_draggedViewport = viewport; m_draggedViewport->grabMouse(); }
                state.startTimeValue = clickedX;
            }
            return true;
        }

        // --- Ampogram / hypnogram click (navigate) ---
        if (cv == ui->amp_ecg_axis || cv == ui->amp_ppg_axis || cv == ui->sleep_state_axis) {
            double localTarget = clickedX - globalOffset;
            m_currentStartTime = std::max(0.0, std::min(localTarget - m_windowDuration / 2.0,
                std::max(0.0, chunkDur - m_windowDuration)));
            handle_data_plot();
            updateAmpogramCursor();
            return true;
        }
    }
    else if (event->type() == QEvent::MouseMove && m_isDragging) {
        return true;
    }
    else if (event->type() == QEvent::MouseButtonRelease) {
        auto* me = static_cast<QMouseEvent*>(event);
        if (me->button() == Qt::LeftButton && m_isDragging) {
            if (m_draggedViewport) { m_draggedViewport->releaseMouse(); m_draggedViewport = nullptr; }
            m_isDragging = false;

            QString label = signalLabelForChartView(cv);
            if (!label.isEmpty() && label == m_dragSignalLabel) {
                double endX = cv->chart()->mapToValue(me->pos()).x();
                endX = std::clamp(endX, m_currentStartTime, m_currentStartTime + m_windowDuration);
                ChannelMarkingState& state = markStateFor(label);
                if (std::abs(endX - state.startTimeValue) > 0.1) finalizeMarking(cv, endX, label);
                return true;
            }
            m_dragSignalLabel.clear();
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

        double x = m_currentChunkIndex * CHUNK_DURATION_SEC + m_currentStartTime + m_windowDuration / 2.0;
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

    auto fetchAxis = [](QChart* chart, Qt::Orientation orient) -> QAbstractAxis* {
        if (!chart) return nullptr;
        auto axes = chart->axes(orient);
        return axes.isEmpty() ? nullptr : axes.first();
        };

    struct ChartAxes {
        QChart* chart = nullptr;
        QAbstractAxis* xAxis = nullptr;
        QValueAxis* yAxis = nullptr;
    };
    QMap<QString, ChartAxes> axesMap;
    for (const QString& lbl : { "ECG1", "ECG2", "ECG3", "PPG" }) {
        QChartView* cv = chartViewForSignalLabel(lbl);
        if (!cv) continue;
        ChartAxes ca;
        ca.chart = cv->chart();
        ca.xAxis = fetchAxis(ca.chart, Qt::Horizontal);
        ca.yAxis = qobject_cast<QValueAxis*>(fetchAxis(ca.chart, Qt::Vertical));
        axesMap[lbl] = ca;
    }

    double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
    double viewStart = m_currentStartTime + globalOffset;
    double viewEnd = viewStart + m_windowDuration;

    auto createHighlight = [&](QChart* chart, QAbstractAxis* xAxis, QValueAxis* yAxis,
        double start, double end, const QColor& color) {
            if (!chart || !xAxis || !yAxis) return;
            auto* upper = new QLineSeries();
            auto* lower = new QLineSeries();
            upper->append({ {start, yAxis->max()}, {end, yAxis->max()} });
            lower->append({ {start, yAxis->min()}, {end, yAxis->min()} });
            auto* area = new QAreaSeries(upper, lower);
            area->setBrush(color);
            area->setPen(Qt::NoPen);
            chart->addSeries(area);
            area->attachAxis(xAxis);
            area->attachAxis(yAxis);
            m_highlights.append(area);
        };

    for (const auto& seg : m_noiseManager->getSegments()) {
        QString segLabel = QString::fromStdString(seg.label);
        double sr = sampleRateForSignal(segLabel);
        double segStart = seg.startSample / sr;
        double segEnd = seg.endSample / sr;

        if (segEnd < viewStart || segStart > viewEnd) continue;

        double ds = std::max(segStart, viewStart);
        double de = std::min(segEnd, viewEnd);
        QColor color = MARKING_COLORS.value(QString::fromStdString(seg.marking_type), QColor(0, 0, 0, 100));

        if (axesMap.contains(segLabel)) {
            const ChartAxes& ca = axesMap[segLabel];
            createHighlight(ca.chart, ca.xAxis, ca.yAxis, ds, de, color);
        }
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
    return QString("%1:%2:%3").arg(hours, 2, 10, QChar('0'))
        .arg(minutes, 2, 10, QChar('0'))
        .arg(secs, 5, 'f', 2, QChar('0'));
}

// ============================================================================
// Slot handlers
// ============================================================================

void noise_marking_gui::on_skip_interval_box_editingFinished() {
    m_skipInterval = ui->skip_interval_box->text().toDouble();
    ui->skip_interval_box->setText(QString::number(m_skipInterval, 'f', 1));
    ui->skip_interval_box->clearFocus();
}

void noise_marking_gui::on_skip_interval_box_returnPressed() { on_skip_interval_box_editingFinished(); }

void noise_marking_gui::on_marking_type_currentTextChanged(const QString& text) { m_currentMarkingType = text; }

void noise_marking_gui::start_marking_button_clicked(const QString& signalLabel) {
    ChannelMarkingState& state = markStateFor(signalLabel);
    state.isWaitingForStart = false;
    state.isWaitingForEnd = false;
    clearStartMarker(state);

    if (signalLabel.startsWith("ECG")) {
        ui->startNoiseECG->setStyleSheet("");
        ui->stopNoiseECG->setStyleSheet("");
        ui->stopNoiseECG->setEnabled(false);
    }
    else {
        ui->startNoisePPG->setStyleSheet("");
        ui->stopNoisePPG->setStyleSheet("");
        ui->stopNoisePPG->setEnabled(false);
    }
}