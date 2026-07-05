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
#include "annotation_types.hpp"


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
#include <QPixmap>
#include <QPainter>
#include <QIcon>
#include <QRadioButton>

 // ============================================================================
 // Channel lookup
 // ============================================================================

const QStringList& noise_marking_gui::markableChannelLabels() {
    static const QStringList lables{ "ECG1","ECG2","ECG3","PPG_ACCEL","ABP","ART","ART_PULM" };
    return lables;
}

noise_marking_gui::data_channel_features
noise_marking_gui::channelRefs(const QString& label) const {
    auto* self = const_cast<noise_marking_gui*>(this);
    data_channel_features r;

    if (label == "ECG1") {
        r.chartView = ui->ecg_axis_1;
        r.state = &self->m_markState_ecg1;
        r.upsampled_data = &m_ecg1; 
        r.dataRaw = &m_ecg1Raw;
        r.sampleRate = channel_upsampled_rates[CH_ECG1];
        r.color = COLOR_ECG1;
    }
    else if (label == "ECG2") {
        r.chartView = ui->ecg_axis_2;
        r.state = &self->m_markState_ecg2;
        r.upsampled_data = &m_ecg2; 
        r.dataRaw = &m_ecg2Raw; 
        r.sampleRate = channel_upsampled_rates[CH_ECG2];
        r.color = COLOR_ECG2;
    }
    else if (label == "ECG3") {
        r.chartView = ui->ecg_axis_3;
        r.state = &self->m_markState_ecg3;
        r.upsampled_data = &m_ecg3;
        r.dataRaw = &m_ecg3Raw;
        r.sampleRate = channel_upsampled_rates[CH_ECG3];
        r.color = COLOR_ECG3;
    }
    else if (label == "PPG_ACCEL") {
        r.chartView = ui->ppg_accel_axis;
        r.state = &self->m_markState_ppg;
        if (m_cfg.dataset_type == "BITTIUM") {
            r.upsampled_data = &m_accelX; 
            r.dataRaw = &m_accelXRaw;
            r.sampleRate = channel_upsampled_rates[CH_ACCEL_X];
            r.color = COLOR_ACCEL_X;
        }
        else {
            r.upsampled_data = &m_ppg;
            r.dataRaw = &m_ppgRaw;
            r.sampleRate = channel_upsampled_rates[CH_PPG];
            r.color = COLOR_PPG;
        }
    }
    else if (label == "ABP") {
        r.chartView = ui->accel_or_abp_axis;
        r.state = &self->m_markState_abp;
        r.upsampled_data = &m_abp; 
        r.dataRaw = &m_abpRaw;
        r.sampleRate = channel_upsampled_rates[CH_ABP];
        r.color = COLOR_ABP;
    }
    else if (label == "ART") {
        r.chartView = ui->art_axis;
        r.state = &self->m_markState_art;
        r.upsampled_data = &m_art;
        r.dataRaw = &m_artRaw;
        r.sampleRate = channel_upsampled_rates[CH_ART];
        r.color = COLOR_ART;
    }
    else if (label == "ART_PULM") {
        r.chartView = ui->art_pulm_axis;
        r.state = &self->m_markState_art_pulm;
        r.upsampled_data = &m_artPulm;
        r.dataRaw = &m_artPulmRaw;
        r.sampleRate = channel_upsampled_rates[CH_ART_PULM];
        r.color = COLOR_ART_PULM;
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
    if (cv == ui->art_axis) return "ART";
    if (cv == ui->art_pulm_axis) return "ART_PULM";
    return {};
}

QChartView* noise_marking_gui::chartViewForSignalLabel(const QString& label) const {
    return channelRefs(label).chartView;
}

double noise_marking_gui::sampleRateForSignal(const QString& label) const {
    return channelRefs(label).sampleRate;
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
    else if (label == "ART") { check = ui->art_check; gain = ui->art_gain; }
    else if (label == "ART_PULM") { check = ui->art_pulm_check; gain = ui->art_pulm_gain; }
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

void noise_marking_gui::setMarkScope(MarkScope scope) {
    m_markScope = scope;
}

void noise_marking_gui::toggleAnnotationArm() {
    m_markArmed = !m_markArmed;
    if (m_markArmed) {
        m_currentMarkingType = ui->marking_type->currentText();
    }
    else {
        if (m_isDragging) {
            if (m_draggedViewport) { m_draggedViewport->releaseMouse(); m_draggedViewport = nullptr; }
            m_isDragging = false;
        }
        clearDragPreview();
        for (const QString& lbl : markableChannelLabels()) cancelMarking(lbl);
    }
    updateMarkingButtons();
}

void noise_marking_gui::updateMarkingButtons() {
    const bool anyActive = !m_activeChannels.isEmpty();
    const bool anyEcg = isChannelActive("ECG1") || isChannelActive("ECG2") || isChannelActive("ECG3");
    if (!anyActive) m_markArmed = false;

    ui->make_annotation->setEnabled(anyActive);
    ui->mark_one_chan->setEnabled(anyActive);
    ui->mark_all_chan->setEnabled(anyActive);
    ui->mark_ecg->setEnabled(anyEcg);

    applyMarkStyle(ui->make_annotation, m_markArmed ? "armed" : "");
}

QStringList noise_marking_gui::scopeChannels(const QString& clickedLabel) const {
    QStringList out;
    switch (m_markScope) {
    case MarkScope::One:
        if (isChannelActive(clickedLabel)) out << clickedLabel;
        break;
    case MarkScope::Ecg:
        for (const char* l : { "ECG1", "ECG2", "ECG3" })
            if (isChannelActive(l)) out << l;
        break;
    case MarkScope::All:
        for (const QString& l : markableChannelLabels())
            if (isChannelActive(l)) out << l;
        break;
    }
    return out;
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
    reset(ui->art_check, ui->art_gain);
    reset(ui->art_pulm_check, ui->art_pulm_gain);
}

double noise_marking_gui::totalChunkDuration() const {
    if (m_ecg1.size() > 1 && channel_upsampled_rates[CH_ECG1] > 0)
        return m_ecg1.size() / channel_upsampled_rates[CH_ECG1];
    if (m_ppg.size() > 1 && channel_upsampled_rates[CH_PPG] > 0)
        return m_ppg.size() / channel_upsampled_rates[CH_PPG];
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
    , m_noiseManager(std::make_unique<annotation_handler>())
    , m_buttonHandler(std::make_unique<user_control_handler>(this))
{
    ui->setupUi(this);
    m_buttonHandler->setupConnections();
    for (QRadioButton* r : { ui->mark_one_chan, ui->mark_ecg, ui->mark_all_chan })
        r->setFocusPolicy(Qt::NoFocus);
    ui->mark_one_chan->setChecked(true);

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
       ui->hyp_resp_axis, ui->cvp_eeg_axis, ui->pacemaker_axis
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
    wire_gain(ui->art_check, ui->art_gain);
    wire_gain(ui->art_pulm_check, ui->art_pulm_gain);

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
                if (v->property("bpm").isValid())
                    v->chart()->setTitle(get_chart_title(sigName, nativeHz, pxPerSample,
                        v->property("bpm").toDouble()));
                else
                    v->chart()->setTitle(get_chart_title(sigName, nativeHz, pxPerSample));
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
        series->setPen(QPen(QColor(255, 140, 0), 2));
        view->chart()->addSeries(series);
        };
    addCursor(ui->ecg_ampogram_axis, m_ecgCursorBar);
    addCursor(ui->ppg_ampogram_axis, m_ppgCursorBar);

    auto* hypnoChart = new QChart();
    hypnoChart->legend()->hide();
    hypnoChart->layout()->setContentsMargins(0, 0, 0, 0);
    ui->hyp_resp_axis->setChart(hypnoChart);
    m_hypnoCursorBar = new QLineSeries();
    m_hypnoCursorBar->setPen(QPen(QColor(255, 140, 0), 2));
    hypnoChart->addSeries(m_hypnoCursorBar);

    for (int i = 0; i < ui->marking_type->count(); ++i)
        if (const auto* t = annotation_types::find(ui->marking_type->itemText(i))) {
            QPixmap pm(12, 12); pm.fill(Qt::transparent);
            QPainter p(&pm);
            p.setPen(QPen(Qt::black, 1)); p.setBrush(QColor(t->r, t->g, t->b));
            p.drawRect(0, 0, 11, 11); p.end();
            ui->marking_type->setItemIcon(i, QIcon(pm));
        }
    ui->marking_type->setIconSize(QSize(12, 12));

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

double noise_marking_gui::thresholdAt(const QString& label, double globalTime) const {
    double v = m_cfg.threshold;
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

bool noise_marking_gui::promptThresholdBlanking(const QString& header,
    double& thr, double& blk) {
    QDialog dlg(this);
    dlg.setWindowTitle("Set threshold & blanking");
    auto* form = new QFormLayout(&dlg);
    form->addRow(new QLabel(header, &dlg));

    auto* thrSpin = new QDoubleSpinBox(&dlg);
    thrSpin->setRange(0.0, 1.0); thrSpin->setDecimals(2); thrSpin->setSingleStep(0.05);
    thrSpin->setValue(thr);
    form->addRow("Threshold:", thrSpin);

    auto* blkSpin = new QDoubleSpinBox(&dlg);
    blkSpin->setRange(0.0, 2000.0); blkSpin->setDecimals(0); blkSpin->setSingleStep(10.0);
    qDebug() << "blk =" << blk << " cfg.blanking_period =" << m_cfg.blanking_period;\
    blkSpin->setValue(blk);
    form->addRow("Blanking (ms):", blkSpin);

    auto* buttons = new QDialogButtonBox(
        QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dlg);
    form->addRow(buttons);
    connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);

    if (dlg.exec() != QDialog::Accepted) return false;
    thr = thrSpin->value();
    blk = blkSpin->value();
    return true;
}

void noise_marking_gui::applyParamOverrides(const QStringList& channels,
    double lo, double hi, double thrVal, double blkVal) {
    auto beatCh = [](const QString& l) -> beat_log::ChannelIdx {
        if (l == "ECG1") return beat_log::ECG1;
        if (l == "ECG2") return beat_log::ECG2;
        if (l == "ECG3") return beat_log::ECG3;
        if (l == "PPG_ACCEL") return beat_log::PPG;
        if (l == "ART")       return beat_log::ART;
        if (l == "ART_PULM")  return beat_log::ART_PULM;
        return beat_log::ABP;
        };
    for (const QString& label : channels) {
        auto applyTo = [&](QVector<ParamOverride>& vec, double val) {
            vec.erase(std::remove_if(vec.begin(), vec.end(),
                [&](const ParamOverride& o) {
                    return o.channel == label && o.start <= hi && lo <= o.end;
                }), vec.end());
            vec.append(ParamOverride{ label, lo, hi, val });
            };
        applyTo(m_thresholdOverrides, thrVal);
        applyTo(m_blankingOverrides, blkVal);
        if (m_beatLog) m_beatLog->removeInRange(beatCh(label), lo, hi);
    }
    handle_data_plot();
}

void noise_marking_gui::finalizeParamEdit(const QStringList& channels,
    double globalStart, double globalEnd) {
    clearDragPreview();
    if (channels.isEmpty()) return;

    const double lo = std::min(globalStart, globalEnd);
    const double hi = std::max(globalStart, globalEnd);

    double thr = m_cfg.threshold;
    double blk = m_cfg.blanking_period;
    const QString hdr = QString("%1 over %2\u2013%3 s")
        .arg(channels.join(", ")).arg(lo, 0, 'f', 1).arg(hi, 0, 'f', 1);

    if (promptThresholdBlanking(hdr, thr, blk))
        applyParamOverrides(channels, lo, hi, thr, blk);
    else
        handle_data_plot();
}

bool noise_marking_gui::editParamOverrideAt(QChartView* cv, const QPoint& pos) {
    if (!cv || !cv->chart()) return false;
    const QString clicked = signalLabelForChartView(cv);
    if (clicked.isEmpty()) return false;

    const double globalOffset = current_chunk_index * seconds_in_memory_at_once;
    const double t = cv->chart()->mapToValue(pos).x() + globalOffset;

    auto findAt = [&](const QVector<ParamOverride>& v) -> int {
        for (int i = v.size() - 1; i >= 0; --i)          // last match wins (matches lookup order)
            if (v[i].channel == clicked && v[i].start <= t && t <= v[i].end) return i;
        return -1;
        };
    const int ti = findAt(m_thresholdOverrides);
    const int bi = findAt(m_blankingOverrides);
    if (ti < 0 && bi < 0) return false;                  // not on a blank+thresh mark

    const double lo = (ti >= 0) ? m_thresholdOverrides[ti].start : m_blankingOverrides[bi].start;
    const double hi = (ti >= 0) ? m_thresholdOverrides[ti].end : m_blankingOverrides[bi].end;
    double thr = (ti >= 0) ? m_thresholdOverrides[ti].value : m_cfg.threshold;
    double blk = (bi >= 0) ? m_blankingOverrides[bi].value : m_cfg.blanking_period;

    // The first click of the double-click may have started a drag / click-click;
    // undo that state before opening the editor.
    if (m_draggedViewport) { m_draggedViewport->releaseMouse(); m_draggedViewport = nullptr; }
    m_isDragging = false;
    clearDragPreview();
    for (const QString& ch : scopeChannels(clicked)) {
        ChannelMarkingState& st = markStateFor(ch);
        st.phase = MarkPhase::Idle;
        clearStartMarker(st);
    }

    // Edit every channel sharing this exact span (mirrors how the mark was created).
    QStringList channels;
    auto gather = [&](const QVector<ParamOverride>& v) {
        for (const ParamOverride& o : v)
            if (o.start == lo && o.end == hi && !channels.contains(o.channel))
                channels << o.channel;
        };
    gather(m_thresholdOverrides);
    gather(m_blankingOverrides);
    if (channels.isEmpty()) channels << clicked;

    const QString hdr = QString("Edit %1 over %2\u2013%3 s")
        .arg(channels.join(", ")).arg(lo, 0, 'f', 1).arg(hi, 0, 'f', 1);
    if (promptThresholdBlanking(hdr, thr, blk))
        applyParamOverrides(channels, lo, hi, thr, blk);
    else
        handle_data_plot();
    return true;
}