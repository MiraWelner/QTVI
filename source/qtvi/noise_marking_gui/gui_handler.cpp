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
#include "gui_peak_finder.hpp"
#include "beat_log.hpp"
#include "annotation_eraser.h"

#include <QCheckBox>
#include <QScrollBar>
#include <QDoubleSpinBox>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QMouseEvent>
#include <QShortcut>
#include <QSignalBlocker>
#include <QApplication>
#include <QFileInfo>
#include <QTimer>
#include <QDialog>
#include <QFormLayout>
#include <QLabel>
#include <QDialogButtonBox>
#include <QGraphicsLayout>
#include <QStyle>

 // ============================================================================
 // Channel lookup
 // ============================================================================

const QStringList& noise_marking_gui::markableChannelLabels() {
    static const QStringList lables{ "ECG1", "ECG2", "ECG3", "PPG_ACCEL", "ABP" };
    return lables;
}

noise_marking_gui::data_channel_features
noise_marking_gui::channelRefs(const QString& label) const {
    auto* self = const_cast<noise_marking_gui*>(this);
    data_channel_features r;

    if (label == "ECG1") {
        r.chartView = ui->ecg_axis_1; r.startButton = ui->mark_one_chan;
        r.state = &self->m_markState_ecg1;
        r.upsampled_data = &m_ecg1; r.dataRaw = &m_ecg1Raw; r.sampleRate = &m_ecgSR;
        r.color = COLOR_ECG1;
    }
    else if (label == "ECG2") {
        r.chartView = ui->ecg_axis_2; r.startButton = ui->mark_one_chan;
        r.state = &self->m_markState_ecg2;
        r.upsampled_data = &m_ecg2; r.dataRaw = &m_ecg2Raw; r.sampleRate = &m_ecgSR;
        r.color = COLOR_ECG2;
    }
    else if (label == "ECG3") {
        r.chartView = ui->ecg_axis_3; r.startButton = ui->mark_one_chan;
        r.state = &self->m_markState_ecg3;
        r.upsampled_data = &m_ecg3; r.dataRaw = &m_ecg3Raw; r.sampleRate = &m_ecgSR;
        r.color = COLOR_ECG3;
    }
    else if (label == "PPG_ACCEL") {
        r.chartView = ui->ppg_accel_axis; r.startButton = ui->mark_one_chan;
        r.state = &self->m_markState_ppg;
        if (m_cfg.dataset_type == "BITTIUM") {
            r.upsampled_data = &m_accelX; r.dataRaw = &m_accelXRaw;
            r.sampleRate = &m_ecgSR; r.color = COLOR_ACCEL_X;
        }
        else {
            r.upsampled_data = &m_ppg; r.dataRaw = &m_ppgRaw;
            r.sampleRate = &m_ecgSR; r.color = COLOR_PPG;
        }
    }
    else if (label == "ABP") {
        r.chartView = ui->accel_or_abp_axis; r.startButton = ui->mark_one_chan;
        r.state = &self->m_markState_abp;
        r.upsampled_data = &m_abp; r.dataRaw = &m_abpRaw; r.sampleRate = &m_ecgSR;
        r.color = COLOR_ABP;
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
    if (cv == ui->ppg_accel_axis)   return "PPG_ACCEL";
    if (cv == ui->accel_or_abp_axis) {
        bool anyAccel = !is_missing_signal(m_accelX)
            || !is_missing_signal(m_accelY) || !is_missing_signal(m_accelZ);
        if (!anyAccel && !is_missing_signal(m_abp)) return "ABP";
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
    else if (label == "PPG_ACCEL") { check = ui->ppg_check;   gain = ui->ppg_gain; }
    else if (label == "ABP") { check = ui->abp_check; gain = ui->abp_gain; }
    if (!gain) return 1.0;
    double v = gain->value();
    return (v > 0.0) ? v : 1.0;
}

bool noise_marking_gui::invertedForSignal(const QString& label) const {
    QCheckBox* c = nullptr;
    if (label == "ECG1") c = ui->ecg_1_reverse;
    else if (label == "ECG2") c = ui->ecg_2_reverse;
    else if (label == "ECG3") c = ui->ecg_3_reverse;
    return c && c->isChecked();   // PPG/ABP have no reverse box -> never inverted
}

// ============================================================================
// Button & mark-mode helpers
// ============================================================================

void applyMarkStyle(QPushButton* btn, const char* phase) {
    if (!btn) return;
    btn->setProperty("markphase", (phase && *phase) ? QVariant(phase) : QVariant());
    btn->style()->unpolish(btn);
    btn->style()->polish(btn);
    btn->update();
}

QPushButton* noise_marking_gui::startButtonForSignal(const QString& label) const {
    return channelRefs(label).startButton;
}
QPushButton* noise_marking_gui::stopButtonForSignal(const QString& label) const {
    return channelRefs(label).stopButton;
}

void noise_marking_gui::updateAllChannelButtonStates() {
    for (const QString& label : markableChannelLabels())
        updateEcgMarkButtonStyle();

    const bool anyActive = !m_activeChannels.isEmpty();
    const bool anyEcgActive = isChannelActive("ECG1") || isChannelActive("ECG2") || isChannelActive("ECG3");
    ui->mark_all_chan->setEnabled(anyActive);
    ui->mark_all_ecg->setEnabled(anyEcgActive);

    // Each group button reflects ONLY its own mode.
    auto groupStyle = [&](QPushButton* btn, MarkAllMode mode, const QStringList& chans) {
        if (m_markAllMode != mode) { applyMarkStyle(btn, ""); return; }
        bool anyWaitingStop = false;
        for (const QString& lbl : chans)
            if (isChannelActive(lbl) && markStateFor(lbl).phase == MarkPhase::WaitingForStop)
                anyWaitingStop = true;
        applyMarkStyle(btn, anyWaitingStop ? "waiting" : "armed");
        };
    static const QStringList kEcg{ "ECG1", "ECG2", "ECG3" };
    groupStyle(ui->mark_all_chan, MarkAllMode::All, markableChannelLabels());
    groupStyle(ui->mark_all_ecg, MarkAllMode::Ecg, kEcg);
}

void noise_marking_gui::updateEcgMarkButtonStyle() {
    QPushButton* btn = ui->mark_one_chan;
    if (m_activeChannels.isEmpty()) { btn->setEnabled(false); applyMarkStyle(btn, ""); return; }
    btn->setEnabled(true);

    // Single-marker flow only; stay neutral while a group mode owns the channels.
    if (m_markAllMode != MarkAllMode::None) { applyMarkStyle(btn, ""); return; }

    bool waitingStop = false;
    for (const QString& l : markableChannelLabels())
        if (markStateFor(l).phase == MarkPhase::WaitingForStop) { waitingStop = true; break; }

    if (waitingStop)                    applyMarkStyle(btn, "waiting");
    else if (single_ecg_marker_clicked) applyMarkStyle(btn, "armed");
    else                                applyMarkStyle(btn, "");
}

void noise_marking_gui::exitAllMarkModes() {
    for (const QString& label : markableChannelLabels())
        if (markStateFor(label).phase != MarkPhase::Idle) cancelMarking(label);
    m_markAllMode = MarkAllMode::None;
    single_ecg_marker_clicked = false;
}

void noise_marking_gui::toggleEcgMark() {
    if (single_ecg_marker_clicked) {
        single_ecg_marker_clicked = false;     // toggle off
    }
    else {
        exitAllMarkModes();                     // clear other modes first
        m_currentMarkingType = ui->marking_type->currentText();
        single_ecg_marker_clicked = true;
    }
    updateAllChannelButtonStates();
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

double noise_marking_gui::totalChunkDuration() const {
    if (m_ecg1.size() > 1 && m_ecgSR > 0) return m_ecg1.size() / m_ecgSR;
    if (m_ppg.size() > 1 && m_ecgSR > 0) return m_ppg.size() / m_ecgSR;
    return 0.0;
}

double noise_marking_gui::scrollStepSeconds() const {
    //convert the percent of window size set by user to raw seconds
    return percent_of_window_to_shift_on_click / 100.0 * visible_window_size;
}

void noise_marking_gui::syncChunkScrollBar() {
    if (!ui->chunk_scrollbar) return;
    QScrollBar* sb = ui->chunk_scrollbar;

    const double chunkDur = totalChunkDuration();
    const int page = std::max(1, static_cast<int>(std::lround(visible_window_size)));
    const int maxStart = std::max(0,
        static_cast<int>(std::lround(chunkDur)) - page);

    // Programmatic update: block valueChanged so this doesn't re-enter the
    // redraw path that called us.
    QSignalBlocker block(sb);
    sb->setMinimum(0);
    sb->setMaximum(maxStart);
    sb->setPageStep(page);
    sb->setSingleStep(0.1);
    sb->setValue(std::clamp(
        static_cast<int>(std::lround(current_start_time)), 0, maxStart));
    sb->setEnabled(maxStart > 0);
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
        current_start_time = std::max(0.0, current_start_time - scrollStepSeconds());
        resetUnpinnedGains(); handle_data_plot(); updateAmpogramCursor();
        });
    new QShortcut(QKeySequence(Qt::Key_Right), this, [this]() {
        double maxStart = std::max(0.0, totalChunkDuration() - visible_window_size);
        current_start_time = std::min(current_start_time + scrollStepSeconds(), maxStart);
        resetUnpinnedGains(); handle_data_plot(); updateAmpogramCursor();
        });

    const QList<QChartView*> allCharts = {
       ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3,
       ui->ppg_accel_axis, ui->accel_or_abp_axis,
       ui->ecg_ampogram_axis, ui->ppg_ampogram_axis,
       ui->hyp_resp_axis, ui->cvp_axis, ui->pacemaker_axis
    };

    for (auto* view : allCharts) {
        if (!view) continue;
        view->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        view->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        view->setRenderHint(QPainter::Antialiasing);
        view->setFocusPolicy(Qt::NoFocus);
        view->viewport()->installEventFilter(this);

    }

    percent_of_window_to_shift_on_click = ui->skip_interval_box->text().toDouble();
    ui->skip_interval_box->setFocusPolicy(Qt::ClickFocus);

    ui->scatter_line->setCurrentIndex(0);
    ui->scatter_line->setFocusPolicy(Qt::NoFocus);
    ui->window_length_selector->setFocusPolicy(Qt::NoFocus);
    ui->chunk_scrollbar->setFocusPolicy(Qt::NoFocus);

    connect(ui->scatter_line, QOverload<int>::of(&QComboBox::currentIndexChanged),
        this, [this](int index) {
            m_plotMode = (index == 1) ? PlotMode::Scatter : PlotMode::Line;
            handle_data_plot();
        });

    connect(ui->chunk_scrollbar, &QScrollBar::valueChanged, this, [this](int v) {
        const double maxStart = std::max(0.0, totalChunkDuration() - visible_window_size);
        current_start_time = std::clamp(static_cast<double>(v), 0.0, maxStart);
        handle_data_plot();
        updateAmpogramCursor();
        });

    //handle the checkbox that changes max and min of the scaling to the global max and min
    { QSignalBlocker block(ui->global_scaling); ui->global_scaling->setChecked(false); }
    m_filterBaselineDrift = !ui->global_scaling->isChecked();
    ui->global_scaling->setFocusPolicy(Qt::NoFocus);
    connect(ui->global_scaling, &QCheckBox::toggled, this, [this](bool on) {
        m_filterBaselineDrift = !on;
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
                const double pxPerSec = (visible_window_size > 0.0)
                    ? v->chart()->plotArea().width() / visible_window_size : 0.0;
                const double pxPerSample = (nativeHz > 0.0) ? pxPerSec / nativeHz : 0.0;
                const double bpm = v->property("bpm").isValid()
                    ? v->property("bpm").toDouble() : 0.0;
                v->chart()->setTitle(get_chart_title(sigName, nativeHz, pxPerSample, bpm));
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

    auto addCursor = [](QChartView* view, QLineSeries*& series) {
        series = new QLineSeries();
        series->setPen(QPen(QColor(255, 140, 0), 2));   // bright orange
        view->chart()->addSeries(series);
        };
    addCursor(ui->ecg_ampogram_axis, m_ecgCursorBar);
    addCursor(ui->ppg_ampogram_axis, m_ppgCursorBar);

    auto* hypnoChart = new QChart();
    hypnoChart->legend()->hide();
    ui->hyp_resp_axis->setChart(hypnoChart);
    m_hypnoCursorBar = new QLineSeries();
    m_hypnoCursorBar->setPen(QPen(Qt::black, 2));
    hypnoChart->addSeries(m_hypnoCursorBar);

    m_currentMarkingType = ui->marking_type->currentText();

    connect(ui->browse_file_button, &QPushButton::clicked, this, &noise_marking_gui::handleBrowseFile);

    m_pulseOverlay = std::make_unique<pulse_overlay>(this, markableChannelLabels());
    m_gapIndicator = std::make_unique<gap_indicator>(this);
    m_annotationEraser = std::make_unique<annotation_eraser>(this);

    { QSignalBlocker block(ui->show_grid_check); ui->show_grid_check->setChecked(false); }
    ui->show_grid_check->setFocusPolicy(Qt::NoFocus);
    m_pulseOverlay->setEnabled(false);
    connect(ui->show_grid_check, &QCheckBox::toggled, this, [this](bool on) {
        if (!m_pulseOverlay) return;
        m_pulseOverlay->setEnabled(on);
        handle_data_plot();
        });
    for (QCheckBox* c : { ui->ecg_1_reverse, ui->ecg_2_reverse, ui->ecg_3_reverse }) {
        c->setFocusPolicy(Qt::NoFocus);
        connect(c, &QCheckBox::toggled, this, [this](bool) { handle_data_plot(); });
    }

    // Flush the beat log to disk every 30 s: merge the pending buffer into
    // the table (same-time beats overwrite), write the CSV, and empty the
    // buffer. The final partial interval is committed by main on close.
    m_logFlushTimer = new QTimer(this);
    connect(m_logFlushTimer, &QTimer::timeout, this, [this] {
        if (!m_beatLog) return;
        m_beatLog->flushPending();
        const QString stem = QFileInfo(m_binFilePath).completeBaseName();
        m_beatLog->writeCsv(m_cfg.log_path + "/" + stem.toStdString() + "_log.csv");
        });
    m_logFlushTimer->start(30000);   // 30 s
}

noise_marking_gui::~noise_marking_gui() {
    for (auto& list : m_persistentLines)
        for (auto* ln : list) delete ln;
    m_persistentLines.clear();
    for (auto& list : m_persistentRawScatter)
        for (auto* sc : list) delete sc;
    m_persistentRawScatter.clear();
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
    percent_of_window_to_shift_on_click = ui->skip_interval_box->text().toInt();
    if (percent_of_window_to_shift_on_click <= 0) percent_of_window_to_shift_on_click = 50;
    ui->skip_interval_box->setText(
        QString::number(static_cast<int>(percent_of_window_to_shift_on_click)));
    ui->skip_interval_box->clearFocus();
}

void noise_marking_gui::on_marking_type_currentTextChanged(const QString& text) {
    m_currentMarkingType = text;
}

// ============================================================================
// Parameter (threshold / blanking) editing
// ============================================================================

void noise_marking_gui::enterParamEdit() {
    // Toggle the single param-edit mode on/off.
    m_paramEditMode = (m_paramEditMode == ParamEdit::Active)
        ? ParamEdit::None : ParamEdit::Active;
    updateParamButtonStyles();
}

void noise_marking_gui::updateParamButtonStyles() {
    const char* active = "background-color: #2980b9; color: white;";
    ui->param_change->setStyleSheet(
        m_paramEditMode == ParamEdit::Active ? active : "");
}

double noise_marking_gui::thresholdAt(const QString& label, double globalTime) const {
    double v = m_cfg.height_threshold_percent;
    for (const ParamOverride& o : m_thresholdOverrides)            // last match wins
        if (o.channel == label && globalTime >= o.start && globalTime <= o.end) v = o.value;
    return v;
}

double noise_marking_gui::blankingAt(const QString& label, double globalTime) const {
    double v = m_cfg.blanking_period;
    for (const ParamOverride& o : m_blankingOverrides)
        if (o.channel == label && globalTime >= o.start && globalTime <= o.end) v = o.value;
    return v;
}

void noise_marking_gui::finalizeParamEdit(const QString& label,
    double globalStart, double globalEnd) {
    const double lo = std::min(globalStart, globalEnd);
    const double hi = std::max(globalStart, globalEnd);

    // One popup, two fields: threshold and blanking for the dragged span.
    QDialog dlg(this);
    dlg.setWindowTitle("Set threshold & blanking");
    auto* form = new QFormLayout(&dlg);

    auto* lbl = new QLabel(
        QString("%1 over %2\u2013%3 s").arg(label).arg(lo, 0, 'f', 1).arg(hi, 0, 'f', 1),
        &dlg);
    form->addRow(lbl);

    auto* thrSpin = new QDoubleSpinBox(&dlg);
    thrSpin->setRange(0.0, 1.0); thrSpin->setDecimals(2); thrSpin->setSingleStep(0.05);
    thrSpin->setValue(m_cfg.height_threshold_percent);
    form->addRow("Threshold:", thrSpin);

    auto* blkSpin = new QDoubleSpinBox(&dlg);
    blkSpin->setRange(0.0, 1.0); blkSpin->setDecimals(2); blkSpin->setSingleStep(0.05);
    blkSpin->setValue(m_cfg.blanking_period);
    form->addRow("Blanking:", blkSpin);

    auto* buttons = new QDialogButtonBox(
        QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dlg);
    form->addRow(buttons);
    connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);

    if (dlg.exec() == QDialog::Accepted) {
        const double thrVal = thrSpin->value();
        const double blkVal = blkSpin->value();

        // Replace any existing override of each kind overlapping [lo, hi] on
        // this channel, then add the new one. Both vectors get an entry for
        // the same span -- the single gray rectangle in updateNoiseHighlights
        // represents both.
        auto applyTo = [&](QVector<ParamOverride>& vec, double val) {
            vec.erase(std::remove_if(vec.begin(), vec.end(),
                [&](const ParamOverride& o) {
                    return o.channel == label && o.start <= hi && lo <= o.end;
                }), vec.end());
            vec.append(ParamOverride{ label, lo, hi, val });
            };
        applyTo(m_thresholdOverrides, thrVal);
        applyTo(m_blankingOverrides, blkVal);

        // Detection in [lo, hi] changes; drop this channel's logged peaks there
        // so removed beats don't linger. handle_data_plot() re-logs survivors.
        if (m_beatLog) {
            auto beatCh = [](const QString& l) -> beat_log::ChannelIdx {
                if (l == "ECG1") return beat_log::ECG1;
                if (l == "ECG2") return beat_log::ECG2;
                if (l == "ECG3") return beat_log::ECG3;
                if (l == "PPG_ACCEL")  return beat_log::PPG;
                return beat_log::ABP;
                };
            m_beatLog->removeInRange(beatCh(label), lo, hi);
        }
    }
    m_paramEditMode = ParamEdit::None;
    updateParamButtonStyles();
    clearDragPreview();
    handle_data_plot();
}