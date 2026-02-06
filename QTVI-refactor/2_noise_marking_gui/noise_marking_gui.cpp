#include "noise_marking_gui.h"
#include <QtCharts/QAreaSeries>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QValueAxis>
#include <QtCharts/QCategoryAxis>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QTime>
#include <QShortcut>
#include <QFile>
#include <algorithm>
#include <utility>



static const QColor COLOR_ECG = QColor("#1ABC9C");
static const QColor COLOR_PPG = QColor("#E74C3C");


namespace {
    void clearAxes(QChart* chart) {
        if (!chart) return;
        const auto axes = chart->axes();
        for (auto* axis : axes) chart->removeAxis(axis);
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
// CONSTRUCTOR
// ============================================================================
noise_marking_gui::noise_marking_gui(QWidget* parent)
    : QDialog(parent)
    , ui(std::make_unique<Ui::noise_marking_gui>())
    , m_noiseManager(std::make_unique<NoiseManager>(256.0))
    , m_buttonHandler(std::make_unique<lower_row_buttons>(this))
{
    ui->setupUi(this);

    //set large window size
    auto* screen = QGuiApplication::primaryScreen();
    if (screen) {
        const QRect available = screen->availableGeometry();
        const QSize targetSize(available.width() * 3 / 4, available.height() * 0.9);
        resize(targetSize);
        move(available.center() - rect().center());
    }



    m_windowDuration = 10.0;

    m_buttonHandler->setupConnections();

    // Ensure NO buttons get focus when clicked (prevents spacebar/arrow issues)
    const QList<QPushButton*> allButtons = this->findChildren<QPushButton*>();
    for (auto* btn : allButtons) btn->setFocusPolicy(Qt::NoFocus);

    // --- NAVIGATION SHORTCUTS ---
    new QShortcut(QKeySequence(Qt::Key_Left), this, [this]() {
        m_currentStartTime = std::max(0.0, m_currentStartTime - m_skipInterval);
        handle_data_plot();
        updateAmpogramCursor();
        });

    new QShortcut(QKeySequence(Qt::Key_Right), this, [this]() {
        // Calculate the duration of the currently loaded chunk
        double chunkDur = 0.0;
        if (!m_ppg.isEmpty() && m_ppgSR > 0) chunkDur = (double)m_ppg.size() / m_ppgSR;
        else if (!m_ecg.isEmpty() && m_ecgSR > 0) chunkDur = (double)m_ecg.size() / m_ecgSR;

        m_currentStartTime = std::min(m_currentStartTime + m_skipInterval,
            std::max(0.0, chunkDur - m_windowDuration));
        handle_data_plot();
        updateAmpogramCursor();
        });

    // --- CHART VIEWPORT SETUP ---
    const QList<QChartView*> charts = {
        ui->sleep_state_axis, ui->amp_ecg_axis, ui->amp_ppg_axis, ui->ecg_axis, ui->ppg_axis
    };

    for (auto* view : charts) {
        view->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        view->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        view->setRenderHint(QPainter::Antialiasing);
        view->setFocusPolicy(Qt::NoFocus);
        // Important for mouse click navigation
        view->viewport()->installEventFilter(this);
    }

    // --- UI DEFAULTS ---
    ui->rb_10s->setChecked(true);
    m_skipInterval = ui->skip_interval_box->text().toDouble();
    if (m_skipInterval <= 0.0) m_skipInterval = 5.0;

    // --- INITIALIZE CHARTS ---
    setupChartDefaults(ui->ecg_axis);
    setupChartDefaults(ui->ppg_axis);
    setupChartDefaults(ui->amp_ecg_axis);
    setupChartDefaults(ui->amp_ppg_axis);

    // ECG/PPG Ampogram Series
    m_ecgAmpSeries = new QLineSeries();
    m_ppgAmpSeries = new QLineSeries();
    ui->amp_ecg_axis->chart()->addSeries(m_ecgAmpSeries);
    ui->amp_ppg_axis->chart()->addSeries(m_ppgAmpSeries);

    ui->stopNoiseECG->setEnabled(false);
    ui->stopNoisePPG->setEnabled(false);

    // Setup Cursor Bars
    auto setupCursor = [](QChartView* view, QLineSeries*& series) {
        series = new QLineSeries();
        series->setPen(QPen(Qt::black, 2));
        view->chart()->addSeries(series);
        };

    setupCursor(ui->amp_ecg_axis, m_ecgCursorBar);
    setupCursor(ui->amp_ppg_axis, m_ppgCursorBar);

    // Hypnogram Chart
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

void noise_marking_gui::setFileSource(const QString& filePath) {
    m_binFilePath = filePath;
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly)) return;

    double scoring_epoch = 0;
    file.read((char*)&m_ecgSR, sizeof(double));
    file.read((char*)&m_ppgSR, sizeof(double));
    file.read((char*)&scoring_epoch, sizeof(double));
    m_sleepSR = (scoring_epoch > 0) ? (1.0 / scoring_epoch) : 0;

    // Read ALL 5 sample counts from the header
    file.read((char*)&m_totalEcgSamples, sizeof(uint64_t));      // Signal 1 (ECG1)
    file.read((char*)&m_totalSignal2Samples, sizeof(uint64_t));  // Signal 2 (ECG2)
    file.read((char*)&m_totalSignal3Samples, sizeof(uint64_t));  // Signal 3 (ECG3)
    file.read((char*)&m_totalPpgSamples, sizeof(uint64_t));      // PPG samples
    file.read((char*)&m_totalSleepSamples, sizeof(uint64_t));    // Sleep samples

    m_fileHeaderSize = file.pos();
    file.close();

    m_noiseManager = std::make_unique<NoiseManager>(m_ecgSR);
    loadChunkFromFile(0); // Load the first 8 hours immediately
}

bool noise_marking_gui::loadChunkFromFile(uint64_t chunkIndex) {
    QFile file(m_binFilePath);
    if (!file.open(QIODevice::ReadOnly)) return false;

    m_currentChunkIndex = chunkIndex;

    // --- 1. ECG Loading (Signals 1, 2, and 3) ---
    uint64_t ecgIn8Hours = static_cast<uint64_t>(CHUNK_DURATION_SEC * m_ecgSR);
    uint64_t ecg1Start = m_currentChunkIndex * ecgIn8Hours;
    uint64_t ecg1Count = (m_totalEcgSamples > ecg1Start) ? std::min(ecgIn8Hours, m_totalEcgSamples - ecg1Start) : 0;
    m_ecg.resize(static_cast<int>(ecg1Count));
    file.seek(m_fileHeaderSize + (ecg1Start * sizeof(double)));
    file.read((char*)m_ecg.data(), ecg1Count * sizeof(double));

    // ECG Signal 2
    uint64_t ecg2Start = m_currentChunkIndex * ecgIn8Hours;
    uint64_t ecg2Count = (m_totalSignal2Samples > ecg2Start) ? std::min(ecgIn8Hours, m_totalSignal2Samples - ecg2Start) : 0;
    m_ecg2.resize(static_cast<int>(ecg2Count));
    uint64_t ecg2Offset = m_totalEcgSamples;
    file.seek(m_fileHeaderSize + (ecg2Offset + ecg2Start) * sizeof(double));
    file.read((char*)m_ecg2.data(), ecg2Count * sizeof(double));

    // ECG Signal 3
    uint64_t ecg3Start = m_currentChunkIndex * ecgIn8Hours;
    uint64_t ecg3Count = (m_totalSignal3Samples > ecg3Start) ? std::min(ecgIn8Hours, m_totalSignal3Samples - ecg3Start) : 0;
    m_ecg3.resize(static_cast<int>(ecg3Count));
    uint64_t ecg3Offset = ecg2Offset + m_totalSignal2Samples;
    file.seek(m_fileHeaderSize + (ecg3Offset + ecg3Start) * sizeof(double));
    file.read((char*)m_ecg3.data(), ecg3Count * sizeof(double));


    // --- 2. PPG Loading ---
    uint64_t ppgIn8Hours = static_cast<uint64_t>(CHUNK_DURATION_SEC * m_ppgSR);
    uint64_t ppgStart = m_currentChunkIndex * ppgIn8Hours;
    uint64_t ppgCount = (m_totalPpgSamples > ppgStart) ? std::min(ppgIn8Hours, m_totalPpgSamples - ppgStart) : 0;
    m_ppg.resize(static_cast<int>(ppgCount));
    // PPG starts after Signal 1 + Signal 2 + Signal 3
    uint64_t ppgOffset = m_totalEcgSamples + m_totalSignal2Samples + m_totalSignal3Samples;
    file.seek(m_fileHeaderSize + (ppgOffset * sizeof(double)) + (ppgStart * sizeof(double)));
    file.read((char*)m_ppg.data(), ppgCount * sizeof(double));

    // --- 3. Sleep State Loading ---
    uint64_t sleepIn8Hours = static_cast<uint64_t>(CHUNK_DURATION_SEC * m_sleepSR);
    uint64_t sleepStart = m_currentChunkIndex * sleepIn8Hours;
    uint64_t sleepCount = (m_totalSleepSamples > sleepStart) ? std::min(sleepIn8Hours, m_totalSleepSamples - sleepStart) : 0;
    m_sleepStages.resize(static_cast<int>(sleepCount));
    // Sleep starts after ALL signals + PPG
    uint64_t sleepOffset = ppgOffset + m_totalPpgSamples;
    file.seek(m_fileHeaderSize + (sleepOffset * sizeof(double)) + (sleepStart * sizeof(double)));
    file.read((char*)m_sleepStages.data(), sleepCount * sizeof(double));

    file.close();

    // --- 4. Update UI Labels ---
    int startHr = m_currentChunkIndex * 8;
    int endHr = startHr + 8;
    ui->topLabel->setText(QString("     Data Range: Hour %1 to Hour %2").arg(startHr).arg(endHr));

    // --- 5. Refresh Plots ---
    m_currentStartTime = 0;
    handle_ampogram_plot();
    handle_data_plot();
    updateAmpogramCursor();
    setupHypnogram();

    // Navigation buttons
    ui->prev8hours->setEnabled(m_currentChunkIndex > 0);
    ui->next8hours->setEnabled((ecg1Start + ecg1Count) < m_totalEcgSamples);

    return true;
}

void noise_marking_gui::on_next8hours_clicked() { loadChunkFromFile(m_currentChunkIndex + 1); }
void noise_marking_gui::on_prev8hours_clicked() { if (m_currentChunkIndex > 0) loadChunkFromFile(m_currentChunkIndex - 1); }



void noise_marking_gui::setupHypnogram() {
    if (m_sleepSR <= 0.0) return;
    auto* chart = ui->sleep_state_axis->chart();

    for (auto* s : m_hypnoStageSeries) { chart->removeSeries(s); delete s; }
    m_hypnoStageSeries.clear();

    struct Stg { int v; QColor c; };
    const QList<Stg> sts = { {0, Qt::black}, {1, Qt::darkGreen}, {2, Qt::blue}, {3, Qt::cyan}, {4, Qt::red} };

    double dt = 1.0 / m_sleepSR;
    double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;

    // Only draw if we have sleep data loaded for this chunk
    if (!m_sleepStages.isEmpty()) {
        for (const auto& st : sts) {
            auto* s = new QScatterSeries();
            s->setColor(st.c); s->setMarkerSize(3.0); s->setPen(Qt::NoPen);
            s->setMarkerShape(QScatterSeries::MarkerShapeRectangle);
            for (int i = 0; i < m_sleepStages.size(); ++i) {
                if (static_cast<int>(m_sleepStages[i]) == st.v) {
                    s->append(globalOffset + (i * dt + dt / 2.0), st.v);
                }
            }
            chart->addSeries(s);
            m_hypnoStageSeries.append(s);
        }
    }

    clearAxes(chart);
    auto* x_axis = new QCategoryAxis();
    x_axis->setRange(globalOffset, globalOffset + CHUNK_DURATION_SEC);

    // Hour labels 0-8
    for (int h = 0; h <= 8; ++h) {
        x_axis->append(QString::number(h), globalOffset + (h * 3600.0));
    }

    x_axis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
    x_axis->setGridLineVisible(false);
    x_axis->setLabelsFont(QFont("Arial", 6));


    auto* y_axis = new QCategoryAxis();
    for (const auto& st : sts) y_axis->append("", st.v + 0.4);
    y_axis->setRange(-0.5, 4.5);
    y_axis->setReverse(true); y_axis->setVisible(false); y_axis->setGridLineVisible(false);

    chart->addAxis(x_axis, Qt::AlignBottom);
    chart->addAxis(y_axis, Qt::AlignLeft);

    for (auto* s : chart->series()) { s->attachAxis(x_axis); s->attachAxis(y_axis); }

    // Re-attach cursor
    if (m_hypnoCursorBar) {
        chart->removeSeries(m_hypnoCursorBar);
        chart->addSeries(m_hypnoCursorBar);
        m_hypnoCursorBar->attachAxis(x_axis);
        m_hypnoCursorBar->attachAxis(y_axis);
    }
}

void noise_marking_gui::handle_ampogram_plot(double sampling_length) {
    double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
    const double FIXED_8_HOURS = CHUNK_DURATION_SEC;

    auto calcPts = [sampling_length, globalOffset](const QVector<double>& d, double sr) {
        QList<QPointF> p;
        if (d.isEmpty() || sr <= 0.0) return p;
        const double localDur = d.size() / sr;
        for (double t = 0; t <= localDur - sampling_length; t += sampling_length) {
            int s = static_cast<int>(t * sr);
            int e = static_cast<int>((t + sampling_length) * sr);
            auto [mi, ma] = std::minmax_element(d.begin() + s, d.begin() + e);
            p.append({ globalOffset + t, *ma - *mi });
        }
        return p;
        };

    auto design_plot = [this, globalOffset, FIXED_8_HOURS](QChartView* v, QLineSeries* s, const QList<QPointF>& p, QLineSeries* cur) {
        s->replace(p);
        s->setPen(QPen((v == ui->amp_ecg_axis) ? COLOR_ECG : COLOR_PPG, 1));
        auto* chart = v->chart();
        clearAxes(chart);
        chart->legend()->hide();

        auto* x_axis = new QCategoryAxis();
        x_axis->setRange(globalOffset, globalOffset + FIXED_8_HOURS);

        for (int h = 0; h <= 8; ++h) {
            x_axis->append(QString::number(h), globalOffset + (h * 3600.0));
        }

        x_axis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
        x_axis->setGridLineVisible(false);
        x_axis->setLabelsVisible(false);

        chart->addAxis(x_axis, Qt::AlignBottom);
        s->attachAxis(x_axis);
        if (cur) cur->attachAxis(x_axis);

        double yMin = 0, yMax = 1.0;
        if (!p.isEmpty()) {
            auto [mi, ma] = std::minmax_element(p.begin(), p.end(), [](const QPointF& a, const QPointF& b) { return a.y() < b.y(); });
            const double pad = std::max(0.5, (ma->y() - mi->y()) * 0.05);
            yMin = mi->y() - pad; yMax = ma->y() + pad;
        }
        auto* y_axis = new QValueAxis();
        y_axis->setRange(yMin, yMax);
        y_axis->setVisible(false); y_axis->setGridLineVisible(false);
        chart->addAxis(y_axis, Qt::AlignLeft);
        s->attachAxis(y_axis);
        if (cur) cur->attachAxis(y_axis);
        };

    design_plot(ui->amp_ecg_axis, m_ecgAmpSeries, calcPts(m_ecg, m_ecgSR), m_ecgCursorBar);
    design_plot(ui->amp_ppg_axis, m_ppgAmpSeries, calcPts(m_ppg, m_ppgSR), m_ppgCursorBar);
}

bool noise_marking_gui::eventFilter(QObject* watched, QEvent* event) {
    auto* viewport = qobject_cast<QWidget*>(watched);
    if (!viewport) return QDialog::eventFilter(watched, event);

    auto* cv = qobject_cast<QChartView*>(viewport->parent());
    if (!cv || !cv->chart()) return QDialog::eventFilter(watched, event);

    if (event->type() == QEvent::MouseButtonPress) {
        auto* mouseEvent = static_cast<QMouseEvent*>(event);
        if (mouseEvent->button() == Qt::LeftButton) {
            double clickedX = cv->chart()->mapToValue(mouseEvent->pos()).x();

            double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;

            double chunkDur = 0.0;
            if (!m_ppg.isEmpty() && m_ppgSR > 0)
                chunkDur = (double)m_ppg.size() / m_ppgSR;
            else if (!m_ecg.isEmpty() && m_ecgSR > 0)
                chunkDur = (double)m_ecg.size() / m_ecgSR;

            // Clamp clickedX to the global chunk range (prevent clicks before globalOffset, i.e., before chunk start)
            clickedX = std::clamp(clickedX, globalOffset, globalOffset + chunkDur);



            if (cv == ui->ecg_axis || cv == ui->ppg_axis) {
                bool isECG = (cv == ui->ecg_axis);

                // Button-driven workflow (click start / click end)
                if (isECG) {
                    if (m_isWaitingForECGStart) {
                        m_ecgStartTimeValue = clickedX;
                        showStartMarker(cv, clickedX, true);
                        m_isWaitingForECGStart = false;
                        m_isWaitingForECGEnd = true;
                        return true;
                    }
                    else if (m_isWaitingForECGEnd) {
                        finalizeMarking(cv, clickedX, true);
                        return true;
                    }
                }
                else {
                    if (m_isWaitingForPPGStart) {
                        m_ppgStartTimeValue = clickedX;
                        showStartMarker(cv, clickedX, false);
                        m_isWaitingForPPGStart = false;
                        m_isWaitingForPPGEnd = true;
                        return true;
                    }
                    else if (m_isWaitingForPPGEnd) {
                        finalizeMarking(cv, clickedX, false);
                        return true;
                    }
                }

                if (!m_isDragging) {
                    m_isDragging = true;
                    m_dragStartPos = mouseEvent->pos();
                    if (!m_draggedViewport) {
                        m_draggedViewport = viewport;
                        m_draggedViewport->grabMouse();
                    }
                    if (isECG) m_ecgStartTimeValue = clickedX;
                    else m_ppgStartTimeValue = clickedX;
                }
                return true;
            }
            else if (cv == ui->amp_ecg_axis || cv == ui->amp_ppg_axis || cv == ui->sleep_state_axis) {
                double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;

                double chunkDur = 0.0;
                if (!m_ppg.isEmpty() && m_ppgSR > 0)
                    chunkDur = (double)m_ppg.size() / m_ppgSR;
                else if (!m_ecg.isEmpty() && m_ecgSR > 0)
                    chunkDur = (double)m_ecg.size() / m_ecgSR;

                double localTarget = clickedX - globalOffset;

                // Centers the 10s window on the click, while staying within bounds
                m_currentStartTime = std::max(0.0, std::min(localTarget - (m_windowDuration / 2.0),
                    std::max(0.0, chunkDur - m_windowDuration)));

                handle_data_plot();
                updateAmpogramCursor();
                return true;
            }

        }
    }
    else if (event->type() == QEvent::MouseMove && m_isDragging) {
        // Optional: provide visual feedback here if desired
        return true;
    }
    else if (event->type() == QEvent::MouseButtonRelease) {
        auto* mouseEvent = static_cast<QMouseEvent*>(event);
        if (mouseEvent->button() == Qt::LeftButton && m_isDragging) {
            if (m_draggedViewport) {
                m_draggedViewport->releaseMouse();
                m_draggedViewport = nullptr;
            }

            m_isDragging = false;
            if (cv == ui->ecg_axis || cv == ui->ppg_axis) {
                double endX = cv->chart()->mapToValue(mouseEvent->pos()).x();
                // Clamp endX as well
                endX = std::clamp(endX, m_currentStartTime, m_currentStartTime + m_windowDuration);
                bool isECG = (cv == ui->ecg_axis);
                double startX = isECG ? m_ecgStartTimeValue : m_ppgStartTimeValue;

                if (std::abs(endX - startX) > 0.1) {
                    finalizeMarking(cv, endX, isECG);
                }
                return true;
            }
        }
    }

    return QDialog::eventFilter(watched, event);
}

void noise_marking_gui::handle_data_plot() {
    // 1. CRITICAL: Clear highlight list first to prevent dangling pointer crashes 
    // when chart->removeSeries() is called below.
    m_highlights.clear();

    auto plotSignal = [&](QChartView* view, const QVector<double>& data, double sr,
        QLineSeries* marker, double markerPos, const QColor& color, bool resetAxes) {

            if (!view || !view->chart()) return std::make_pair(1e9, -1e9);
            QChart* chart = view->chart();

            if (resetAxes) {
                // Clear old series and axes
                for (auto* s : chart->series()) if (s != marker) { chart->removeSeries(s); delete s; }
                for (auto* a : chart->axes()) { chart->removeAxis(a); delete a; }

                auto* xA = new QCategoryAxis();
                xA->setRange(m_currentStartTime, m_currentStartTime + m_windowDuration);


                double offset = m_currentChunkIndex * CHUNK_DURATION_SEC;

                for (int i = 0; i <= 4; ++i) {
                    double val = m_currentStartTime + (i * m_windowDuration / 4.0);
                    double t = offset + val;

                    // Format hh:mm:ss.ss
                    int h = (int)(t / 3600);
                    int m = (int)(fmod(t, 3600) / 60);
                    double s = fmod(t, 60.0);
                    QString label = QString("%1:%2:%3").arg(h, 2, 10, QChar('0')).arg(m, 2, 10, QChar('0')).arg(s, 5, 'f', 2, QChar('0'));

                    xA->append(label, val);
                }

                xA->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
                xA->setGridLineVisible(false);
                xA->setLabelsFont(QFont("Arial", 7));
                if (view == ui->ecg_axis) xA->setLabelsVisible(false);
                chart->addAxis(xA, Qt::AlignBottom);
                chart->setMargins(QMargins(0, 0, 20, 0));

                auto* yA = new QValueAxis();
                yA->setVisible(false);
                chart->addAxis(yA, Qt::AlignLeft);
            }

            auto hAxes = chart->axes(Qt::Horizontal);
            auto vAxes = chart->axes(Qt::Vertical);
            if (hAxes.isEmpty() || vAxes.isEmpty()) return std::make_pair(1e9, -1e9);

            if (data.size() < 2 || sr <= 0.0) return std::make_pair(1e9, -1e9);

            auto* series = new QLineSeries();
            series->setUseOpenGL(true);
            series->setPen(QPen(color, 1));
            chart->addSeries(series);

            int startIdx = std::clamp((int)(m_currentStartTime * sr), 0, (int)data.size() - 1);
            int endIdx = std::clamp((int)((m_currentStartTime + m_windowDuration) * sr), 0, (int)data.size());

            QList<QPointF> pts;
            double lMin = 1e9, lMax = -1e9;
            for (int i = startIdx; i < endIdx; ++i) {
                pts.append({ (double)i / sr, data[i] });
                if (data[i] < lMin) lMin = data[i];
                if (data[i] > lMax) lMax = data[i];
            }
            series->replace(pts);
            series->attachAxis(hAxes.first());
            series->attachAxis(vAxes.first());

            if (resetAxes) {
                auto* yA = qobject_cast<QValueAxis*>(vAxes.first());
                yA->setRange(lMin - 0.5, lMax + 0.5);
                if (marker && marker->chart() == chart) {
                    marker->replace({ {markerPos, yA->min()}, {markerPos, yA->max()} });
                    marker->attachAxis(hAxes.first()); marker->attachAxis(yA);
                }
            }
            return std::make_pair(lMin, lMax);
        };

    // Redraw everything
    plotSignal(ui->ecg_axis, m_ecg, m_ecgSR, m_ecgStartMarkerLine, m_ecgStartTimeValue, COLOR_ECG, true);
    plotSignal(ui->ecg_axis, m_ecg2, m_ecgSR, nullptr, 0.0, QColor("#3498DB"), false);
    plotSignal(ui->ecg_axis, m_ecg3, m_ecgSR, nullptr, 0.0, QColor("#9B59B6"), false);
    plotSignal(ui->ppg_axis, m_ppg, m_ppgSR, m_ppgStartMarkerLine, m_ppgStartTimeValue, COLOR_PPG, true);

    updateNoiseHighlights();
}


void noise_marking_gui::finalizeMarking(QChartView* cv, double endX, bool isECG) {
    double startX = isECG ? m_ecgStartTimeValue : m_ppgStartTimeValue;
    double s = std::min(startX, endX);
    double e = std::max(startX, endX);
    double sr = isECG ? m_ecgSR : m_ppgSR;
    QString label = isECG ? "ECG" : "PPG";

    // 1. Snap to the nearest sample interval to eliminate floating point drift
    double snappedS = std::round(s * sr) / sr;
    double snappedE = std::round(e * sr) / sr;

    // 2. Add the global offset (current 8-hour chunk position)
    double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
    double finalGlobalStart = snappedS + globalOffset;
    double finalGlobalEnd = snappedE + globalOffset;

    // 3. Save the segment to NoiseManager (used for on-screen highlights)
    // Formula: (Time in Seconds * Sample Rate) = Exact Sample Index
    m_noiseManager->addSegment(
        finalGlobalStart * sr,
        finalGlobalEnd * sr,
        label.toStdString(),
        m_currentMarkingType.toStdString()
    );

    // 4. Save to the persistent list (used for CSV export and Undo)
    m_genExc.noiseExc.append({ finalGlobalStart, finalGlobalEnd });
    m_genExc.data_type.append(label);
    m_genExc.marking_type.append(m_currentMarkingType);

    // 5. Reset the UI state
    if (isECG) {
        m_isWaitingForECGStart = m_isWaitingForECGEnd = false;
        clearECGStartMarker();
        ui->stopNoiseECG->setEnabled(false);
    }
    else {
        m_isWaitingForPPGStart = m_isWaitingForPPGEnd = false;
        clearPPGStartMarker();
        ui->stopNoisePPG->setEnabled(false);
    }

    updateNoiseHighlights();
}


void noise_marking_gui::updateAmpogramCursor() {
    auto drw = [this](QChartView* v, QLineSeries* s) {
        if (!v || !s) return;

        // Ensure the chart has axes before trying to calculate range
        auto axes = v->chart()->axes(Qt::Vertical);
        if (axes.isEmpty()) return;
        auto* y = qobject_cast<QValueAxis*>(axes.first());
        if (!y) return;

        // Calculate Global Position: (Chunk Offset) + (Position within Chunk) + (Center of 10s window)
        double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
        double x = globalOffset + m_currentStartTime + (m_windowDuration / 2.0);

        // Update the line to span the full height of the current Y-axis
        s->replace({ {x, y->min()}, {x, y->max()} });
        };

    drw(ui->amp_ecg_axis, m_ecgCursorBar);
    drw(ui->amp_ppg_axis, m_ppgCursorBar);
    drw(ui->sleep_state_axis, m_hypnoCursorBar);
}

void noise_marking_gui::clearECGStartMarker() {
    // If the start marker cursor/line for ECG exists, remove it from the chart and delete it.
    if (m_ecgStartMarkerLine && m_ecgStartMarkerLine->chart()) {
        m_ecgStartMarkerLine->chart()->removeSeries(m_ecgStartMarkerLine);
        delete m_ecgStartMarkerLine;
        m_ecgStartMarkerLine = nullptr;
    }
}

void noise_marking_gui::clearPPGStartMarker() {
    // If the start marker cursor/line for PPG exists, remove it from the chart and delete it.
    if (m_ppgStartMarkerLine and m_ppgStartMarkerLine->chart()) {
        m_ppgStartMarkerLine->chart()->removeSeries(m_ppgStartMarkerLine);
        delete m_ppgStartMarkerLine; m_ppgStartMarkerLine = nullptr;
    }
}

QString noise_marking_gui::formatTimeLabel(double seconds) {
    \
        /*
        * translates the time in the x axis of the PPG plot from
        * raw seconds to HH:MM:SS.SS
        */
        int hours = static_cast<int>(seconds) / 3600;
    int minutes = (static_cast<int>(seconds) % 3600) / 60;
    double secs = fmod(seconds, 60.0);
    return QString("%1:%2:%3")
        .arg(hours, 2, 10, QChar('0'))
        .arg(minutes, 2, 10, QChar('0'))
        .arg(secs, 5, 'f', 2, QChar('0'));
}

void noise_marking_gui::updateNoiseHighlights() {
    // 1. Clear existing highlights
    for (auto* s : m_highlights) {
        if (s->chart()) {
            s->chart()->removeSeries(s);
        }
        delete s->upperSeries();
        delete s->lowerSeries();
        delete s;
    }
    m_highlights.clear();

    // Helper to safely fetch the first axis from a chart
    auto fetchAxis = [](QChart* chart, Qt::Orientation orient) -> QAbstractAxis* {
        if (!chart) return nullptr;
        auto axes = chart->axes(orient);
        return axes.isEmpty() ? nullptr : axes.first();
        };

    auto* ecgChart = ui->ecg_axis->chart();
    auto* ppgChart = ui->ppg_axis->chart();

    auto* ex = fetchAxis(ecgChart, Qt::Horizontal);
    auto* ey = qobject_cast<QValueAxis*>(fetchAxis(ecgChart, Qt::Vertical));
    auto* px = fetchAxis(ppgChart, Qt::Horizontal);
    auto* py = qobject_cast<QValueAxis*>(fetchAxis(ppgChart, Qt::Vertical));

    double globalOffset = m_currentChunkIndex * CHUNK_DURATION_SEC;
    double viewStartGlobal = m_currentStartTime + globalOffset;
    double viewEndGlobal = viewStartGlobal + m_windowDuration;

    static const QMap<QString, QColor> marking_colors = {
        {"Noise/Artifact",          QColor(255, 255, 0,   30)},
        {"Conduction Delay",        QColor(128, 0,   128, 30)},
        {"AF",                      QColor(255, 0,   0,   30)},
        {"SVT",                     QColor(0,   255, 0,   60)},
        {"VT",                      QColor(0,   0,   255, 60)},
        {"PVC",                     QColor(128, 255, 0,   60)},
        {"PAC",                     QColor(255, 128, 0,   60)},
        {"Benign Arrhythmia",       QColor(255, 128, 255, 60)},
        {"Significant Arrhythmia",  QColor(0,   255, 255, 60)}
    };

    // Define createHighlight lambda with all necessary parameters
    auto createHighlight = [&](QChart* ch, QAbstractAxis* xAxis, QValueAxis* yAxis, double ds, double de, const QColor& color) {
        if (!ch || !xAxis || !yAxis) return;
        auto* upperLine = new QLineSeries();
        auto* lowerLine = new QLineSeries();
        upperLine->append({ {ds, yAxis->max()}, {de, yAxis->max()} });
        lowerLine->append({ {ds, yAxis->min()}, {de, yAxis->min()} });
        auto* area = new QAreaSeries(upperLine, lowerLine);
        area->setBrush(color);
        area->setPen(Qt::NoPen);
        ch->addSeries(area);
        area->attachAxis(xAxis);
        area->attachAxis(yAxis);
        m_highlights.append(area);
        };

    for (const auto& seg : m_noiseManager->getSegments()) {
        double sr = (seg.label == "PPG") ? m_ppgSR : m_ecgSR;
        double segStartGlobal = seg.startSample / sr;
        double segEndGlobal = seg.endSample / sr;

        if (segEndGlobal < viewStartGlobal || segStartGlobal > viewEndGlobal) continue;

        double ds = std::max(segStartGlobal, viewStartGlobal);
        double de = std::min(segEndGlobal, viewEndGlobal);

        QColor highlight_color = marking_colors.value(QString::fromStdString(seg.marking_type),
            QColor(0, 0, 0, 100));

        if (seg.label == "ECG") {
            createHighlight(ecgChart, ex, ey, ds, de, highlight_color);
        }
        else if (seg.label == "PPG") {
            createHighlight(ppgChart, px, py, ds, de, highlight_color);
        }
    }
}


void noise_marking_gui::on_skip_interval_box_editingFinished() {
    m_skipInterval = ui->skip_interval_box->text().toDouble();
    ui->skip_interval_box->setText(QString::number(m_skipInterval, 'f', 1));
    ui->skip_interval_box->clearFocus();
}

void noise_marking_gui::on_skip_interval_box_returnPressed() { on_skip_interval_box_editingFinished(); }


//update the current marking type when the combo box selection changes
void noise_marking_gui::on_marking_type_currentTextChanged(const QString& text) {
    m_currentMarkingType = text;
}

void noise_marking_gui::start_marking_button_clicked(bool isECG) {
    if (isECG) {
        m_isWaitingForECGStart = m_isWaitingForECGEnd = false;
        clearECGStartMarker();
        ui->startNoiseECG->setStyleSheet("");
        ui->stopNoiseECG->setStyleSheet("");
        ui->stopNoiseECG->setEnabled(false);
    }
    else {
        m_isWaitingForPPGStart = m_isWaitingForPPGEnd = false;
        clearPPGStartMarker();
        ui->startNoisePPG->setStyleSheet("");
        ui->stopNoisePPG->setStyleSheet("");
        ui->stopNoisePPG->setEnabled(false);
    }
}

void noise_marking_gui::showStartMarker(QChartView* cv, double xValue, bool isECG) {
    if (isECG) {
        clearECGStartMarker(); // Remove old one if it exists
        m_ecgStartMarkerLine = new QLineSeries();
        m_ecgStartMarkerLine->setPen(QPen(COLOR_ECG, 2, Qt::DashLine));
        auto* ay = qobject_cast<QValueAxis*>(cv->chart()->axes(Qt::Vertical).first());
        m_ecgStartMarkerLine->append(xValue, ay->min());
        m_ecgStartMarkerLine->append(xValue, ay->max());
        cv->chart()->addSeries(m_ecgStartMarkerLine);
        m_ecgStartMarkerLine->attachAxis(cv->chart()->axes(Qt::Horizontal).first());
        m_ecgStartMarkerLine->attachAxis(ay);
        ui->stopNoiseECG->setEnabled(true);
    }
    else {
        clearPPGStartMarker();
        m_ppgStartMarkerLine = new QLineSeries();
        m_ppgStartMarkerLine->setPen(QPen(COLOR_PPG, 2, Qt::DashLine));
        auto* ay = qobject_cast<QValueAxis*>(cv->chart()->axes(Qt::Vertical).first());
        m_ppgStartMarkerLine->append(xValue, ay->min());
        m_ppgStartMarkerLine->append(xValue, ay->max());
        cv->chart()->addSeries(m_ppgStartMarkerLine);
        m_ppgStartMarkerLine->attachAxis(cv->chart()->axes(Qt::Horizontal).first());
        m_ppgStartMarkerLine->attachAxis(ay);
        ui->stopNoisePPG->setEnabled(true);
    }
}
