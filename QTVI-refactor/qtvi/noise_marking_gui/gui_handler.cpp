/**
 * @file   gui_handler.cpp
 * @brief  Construction, channel routing, button state management, and
 *         miscellaneous slot handlers for the noise-marking GUI.
 *
 */

#include "annealing/anneal_handler.hpp"
#include "peak_finding/run_find_r_peaks.hpp"
#include "gui_handler.h"
#include "chart_utils.hpp"
#include "grid_overlay.hpp"
#include "post_process.hpp"
#include "config_loader.hpp"
#include "simple_peak_finder.hpp"

#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QMouseEvent>
#include <QShortcut>
#include <QSignalBlocker>
#include <QApplication>

 // ============================================================================
 // Channel lookup
 // ============================================================================

const QStringList& noise_marking_gui::markableChannelLabels() {
    static const QStringList kLabels{ "ECG1", "ECG2", "ECG3", "PPG", "ABP" };
    return kLabels;
}

noise_marking_gui::data_channel_features
noise_marking_gui::channelRefs(const QString& label) const {
    auto* self = const_cast<noise_marking_gui*>(this);
    data_channel_features r;

    if (label == "ECG1") {
        r.chartView = ui->ecg_axis_1; r.startButton = ui->start_ecg1_mark;
        r.stopButton = ui->stop_ecg1_mark; r.state = &self->m_markState_ecg1;
        r.upsampled_data = &m_ecg1; r.dataRaw = &m_ecg1Raw; r.sampleRate = &m_ecgSR;
        r.color = COLOR_ECG1;
        r.threshold_box = ui->ecg_1_threshold; r.blanking_period_box = ui->ecg_1_blanking_period;
    }
    else if (label == "ECG2") {
        r.chartView = ui->ecg_axis_2; r.startButton = ui->start_ecg2_mark;
        r.stopButton = ui->stop_ecg2_mark; r.state = &self->m_markState_ecg2;
        r.upsampled_data = &m_ecg2; r.dataRaw = &m_ecg2Raw; r.sampleRate = &m_ecgSR;
        r.color = COLOR_ECG2;
        r.threshold_box = ui->ecg_2_threshold; r.blanking_period_box = ui->ecg_2_blanking_period;
    }
    else if (label == "ECG3") {
        r.chartView = ui->ecg_axis_3; r.startButton = ui->start_ecg3_mark;
        r.stopButton = ui->stop_ecg3_mark; r.state = &self->m_markState_ecg3;
        r.upsampled_data = &m_ecg3; r.dataRaw = &m_ecg3Raw; r.sampleRate = &m_ecgSR;
        r.color = COLOR_ECG3;
        r.threshold_box = ui->ecg_3_threshold; r.blanking_period_box = ui->ecg_3_blanking_period;
    }
    else if (label == "PPG") {
        r.chartView = ui->ppg_axis; r.startButton = ui->startNoisePPG;
        r.stopButton = ui->stopNoisePPG; r.state = &self->m_markState_ppg;
        r.upsampled_data = &m_ppg; r.dataRaw = &m_ppgRaw; r.sampleRate = &m_ppgSR;
        r.color = COLOR_PPG;
        r.threshold_box = ui->ppg_threshold; r.blanking_period_box = ui->ppg_blanking_period;
    }
    else if (label == "ABP") {
        r.chartView = ui->accel_or_abp_axis; r.startButton = ui->startNoiseABP;
        r.stopButton = ui->stopNoiseABP; r.state = &self->m_markState_abp;
        r.upsampled_data = &m_abp; r.dataRaw = &m_abpRaw; r.sampleRate = &m_ecgSR;
        r.color = COLOR_ABP;
        r.threshold_box = ui->abp_threshold; r.blanking_period_box = ui->abp_blanking_period;
    }
    return r;
}

noise_marking_gui::ChannelMarkingState&
noise_marking_gui::markStateFor(const QString& label) {
    data_channel_features r = channelRefs(label);
    return r.state ? *r.state : m_markState_ppg;
}

QString noise_marking_gui::signalLabelForChartView(QChartView* cv) const {
    if (cv == ui->ecg_axis_1) return "ECG1";
    if (cv == ui->ecg_axis_2) return "ECG2";
    if (cv == ui->ecg_axis_3) return "ECG3";
    if (cv == ui->ppg_axis)   return "PPG";
    if (cv == ui->accel_or_abp_axis) {
        bool anyAccel = !isMissingSignal(m_accelX)
            || !isMissingSignal(m_accelY) || !isMissingSignal(m_accelZ);
        if (!anyAccel && !isMissingSignal(m_abp)) return "ABP";
    }
    return {};
}

QChartView* noise_marking_gui::chartViewForSignalLabel(const QString& label) const {
    return channelRefs(label).chartView;
}

double noise_marking_gui::sampleRateForSignal(const QString& label) const {
    data_channel_features r = channelRefs(label);
    return r.sampleRate ? *r.sampleRate : m_ecgSR;
}

QColor noise_marking_gui::colorForSignal(const QString& label) const {
    data_channel_features r = channelRefs(label);
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
        check = ui->abp_check; gain = ui->abp_gain;
    }
    if (!gain) return 1.0;
    double v = gain->value();
    return (v > 0.0) ? v : 1.0;
}

// ============================================================================
// Button helpers
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
        startBtn->setEnabled(false); stopBtn->setEnabled(false);
        startBtn->setStyleSheet("color: gray;");
        stopBtn->setStyleSheet("color: gray;");
        return;
    }

    const char* startSheet = ""; const char* stopSheet = "";
    bool stopEnabled = false;

    switch (markStateFor(label).phase) {
    case MarkPhase::WaitingForStart:
        startSheet = "background-color: #f39c12; color: white;"; break;
    case MarkPhase::WaitingForEnd:
        startSheet = "background-color: #f39c12; color: white;";
        stopEnabled = true; break;
    case MarkPhase::WaitingForStop:
        stopSheet = "background-color: #e74c3c; color: white;";
        stopEnabled = true; break;
    default: break;
    }

    startBtn->setEnabled(true); startBtn->setStyleSheet(startSheet);
    stopBtn->setEnabled(stopEnabled); stopBtn->setStyleSheet(stopSheet);
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
// Misc helpers
// ============================================================================

void noise_marking_gui::resetUnpinnedGains() {
    auto reset = [](QCheckBox* check, QDoubleSpinBox* gain) {
        if (!check || !gain || check->isChecked()) return;
        QSignalBlocker block(gain);
        gain->setValue(1.0);
        };
    reset(ui->ecg_1_check, ui->ecg_1_gain);
    reset(ui->ecg_2_check, ui->ecg_2_gain);
    reset(ui->ecg_3_check, ui->ecg_3_gain);
    reset(ui->ppg_check, ui->ppg_gain);
    reset(ui->abp_check, ui->abp_gain);
}

void noise_marking_gui::mousePressEvent(QMouseEvent* event) {
    QWidget* focused = QApplication::focusWidget();
    if (auto* sb = qobject_cast<QDoubleSpinBox*>(focused)) {
        const QPoint global = event->globalPosition().toPoint();
        const QRect sbRect(sb->mapToGlobal(QPoint(0, 0)), sb->size());
        if (!sbRect.contains(global)) setFocus(Qt::MouseFocusReason);
    }
    QDialog::mousePressEvent(event);
}

double noise_marking_gui::totalChunkDuration() const {
    if (m_ecg1.size() > 1 && m_ecgSR > 0) return m_ecg1.size() / m_ecgSR;
    if (m_ppg.size() > 1 && m_ppgSR > 0) return m_ppg.size() / m_ppgSR;
    return 0.0;
}

QString noise_marking_gui::formatTimeLabel(double seconds) {
    return formatHMS(seconds);
}

// ============================================================================
// Constructor / Destructor
// ============================================================================

noise_marking_gui::noise_marking_gui(QWidget* parent)
    : QDialog(parent)
    , ui(std::make_unique<Ui::noise_marking_gui>())
    , m_noiseManager(std::make_unique<annotation_handler>(256.0))
    , m_buttonHandler(std::make_unique<user_control_handler>(this))
{
    ui->setupUi(this);
    m_buttonHandler->setupConnections();

    for (auto* btn : findChildren<QPushButton*>())
        btn->setFocusPolicy(Qt::NoFocus);

    new QShortcut(QKeySequence(Qt::Key_Left), this, [this]() {
        m_currentStartTime = std::max(0.0, m_currentStartTime - m_skipInterval);
        resetUnpinnedGains(); handle_data_plot(); updateAmpogramCursor();
        });
    new QShortcut(QKeySequence(Qt::Key_Right), this, [this]() {
        double maxStart = std::max(0.0, totalChunkDuration() - m_windowDuration);
        m_currentStartTime = std::min(m_currentStartTime + m_skipInterval, maxStart);
        resetUnpinnedGains(); handle_data_plot(); updateAmpogramCursor();
        });

    const QList<QChartView*> allCharts = {
        ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3,
        ui->ppg_axis, ui->accel_or_abp_axis,
        ui->ecg_ampogram_axis, ui->ppg_ampogram_axis,
        ui->hyp_accel_resp_axis, ui->cvp_axis

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
    ui->skip_interval_box->setFocusPolicy(Qt::ClickFocus);

    ui->scatter_line->setCurrentIndex(0);
    ui->scatter_line->setFocusPolicy(Qt::NoFocus);
    ui->window_length_selector->setFocusPolicy(Qt::NoFocus);

    connect(ui->scatter_line, QOverload<int>::of(&QComboBox::currentIndexChanged),
        this, [this](int index) {
            m_plotMode = (index == 1) ? PlotMode::Scatter : PlotMode::Line;
            handle_data_plot();
        });

    ui->show_peaks_check->setChecked(false);
    ui->show_peaks_check->setFocusPolicy(Qt::NoFocus);
    connect(ui->show_peaks_check, &QCheckBox::toggled, this, [this](bool on) {
        m_showPeaks = on; handle_data_plot();
        });
    { QSignalBlocker block(ui->checkBox); ui->checkBox->setChecked(false); }  // default: box unchecked
    m_filterBaselineDrift = !ui->checkBox->isChecked();   // unchecked => drift hidden
    ui->checkBox->setFocusPolicy(Qt::NoFocus);
    connect(ui->checkBox, &QCheckBox::toggled, this, [this](bool on) {
        m_filterBaselineDrift = !on;  // checked = show drift (global); unchecked = hide (window autoscale)
        handle_data_plot();
        });

    auto wire_gain = [this](QCheckBox* /*check*/, QDoubleSpinBox* gain) {
        gain->setDecimals(2); 
        gain->setRange(0.1, 100.0); 
        gain->setValue(1.0);
        gain->setFocusPolicy(Qt::ClickFocus);
        connect(gain, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this](double) { handle_data_plot(); });
        connect(gain, &QDoubleSpinBox::editingFinished, this, [gain]() { gain->clearFocus(); });
        };
    wire_gain(ui->ecg_1_check, ui->ecg_1_gain);
    wire_gain(ui->ecg_2_check, ui->ecg_2_gain);
    wire_gain(ui->ecg_3_check, ui->ecg_3_gain);
    wire_gain(ui->ppg_check, ui->ppg_gain);
    wire_gain(ui->abp_check, ui->abp_gain);

    auto wire_threshold = [this](QDoubleSpinBox* box) {
        box->setDecimals(2);
        box->setRange(0.0, 1.0);
        box->setFocusPolicy(Qt::ClickFocus);
        connect(box, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, [this](double) { handle_data_plot(); });
        connect(box, &QDoubleSpinBox::editingFinished,
            this, [box]() { box->clearFocus(); });
        };
    for (const QString& lbl : markableChannelLabels()) {
        const auto r = channelRefs(lbl);
        wire_threshold(r.threshold_box);
        wire_threshold(r.blanking_period_box);
    }


    for (QChartView* v : allCharts) {
        if (!v) continue;
        auto* chart = new QChart();
        chart->legend()->hide();
        chart->layout()->setContentsMargins(0, 0, 0, 0);
        v->setChart(chart);
        connect(v->chart(), &QChart::plotAreaChanged, v,
            [v, this](const QRectF&) {
                const QString sigName = v->property("signalName").toString();
                if (sigName.isEmpty()) return;
                const double nativeHz = v->property("nativeHz").toDouble();
                const double pxPerSec = (m_windowDuration > 0.0)
                    ? v->chart()->plotArea().width() / m_windowDuration : 0.0;
                const double pxPerSample = (nativeHz > 0.0) ? pxPerSec / nativeHz : 0.0;
                const double bpm = v->property("bpm").isValid()
                    ? v->property("bpm").toDouble() : -1.0;
                v->chart()->setTitle(formatChartTitle(sigName, nativeHz, pxPerSample, bpm));
            });
    }

    ecg1_ampogram_series = new QLineSeries(); ecg1_ampogram_series->setName("ECG1");
    ecg2_ampogram_series = new QLineSeries(); ecg2_ampogram_series->setName("ECG2");
    ecg3_ampogram_series = new QLineSeries(); ecg3_ampogram_series->setName("ECG3");
    ppg_ampogram_series = new QLineSeries();
    ui->ecg_ampogram_axis->chart()->addSeries(ecg1_ampogram_series);
    ui->ecg_ampogram_axis->chart()->addSeries(ecg2_ampogram_series);
    ui->ecg_ampogram_axis->chart()->addSeries(ecg3_ampogram_series);
    ui->ppg_ampogram_axis->chart()->addSeries(ppg_ampogram_series);

    ui->stop_ecg1_mark->setEnabled(false);
    ui->stop_ecg2_mark->setEnabled(false);
    ui->stop_ecg3_mark->setEnabled(false);
    ui->stopNoisePPG->setEnabled(false);
    ui->stop_all_mark->setEnabled(false);

    auto addCursor = [](QChartView* view, QLineSeries*& series) {
        series = new QLineSeries();
        series->setPen(QPen(Qt::black, 2));
        view->chart()->addSeries(series);
        };
    addCursor(ui->ecg_ampogram_axis, m_ecgCursorBar);
    addCursor(ui->ppg_ampogram_axis, m_ppgCursorBar);

    auto* hypnoChart = new QChart();
    hypnoChart->legend()->hide();
    ui->hyp_accel_resp_axis->setChart(hypnoChart);
    m_hypnoCursorBar = new QLineSeries();
    m_hypnoCursorBar->setPen(QPen(Qt::black, 2));
    hypnoChart->addSeries(m_hypnoCursorBar);

    m_currentMarkingType = ui->marking_type->currentText();

    connect(ui->browse_file_button, &QPushButton::clicked, this, &noise_marking_gui::handleBrowseFile);

    m_pulseOverlay = std::make_unique<pulse_overlay>(this, markableChannelLabels());
    m_gapIndicator = std::make_unique<gap_indicator>(this);

    { QSignalBlocker block(ui->show_grid_check); ui->show_grid_check->setChecked(false); }
    ui->show_grid_check->setFocusPolicy(Qt::NoFocus);
    m_pulseOverlay->setEnabled(false);
    connect(ui->show_grid_check, &QCheckBox::toggled, this, [this](bool on) {
        if (!m_pulseOverlay) return;
        m_pulseOverlay->setEnabled(on);
        handle_data_plot();
        });
}

noise_marking_gui::~noise_marking_gui() {
    for (auto& list : m_persistentLines)
        for (auto* ln : list) delete ln;
    m_persistentLines.clear();
}

// ============================================================================
// Accessors
// ============================================================================

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
    for (auto it = all.cbegin(); it != all.cend(); ++it)
        if (!it->noiseExc.isEmpty()) result.append(it.value());
    return result;
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