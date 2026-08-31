/**
 * @file   gui_handler.cpp
 * @brief  Construction, channel routing, button state management, and
 *         miscellaneous slot handlers for the noise-marking GUI.
 *
 */

#include "gui_handler.h"
#include "chart_utils.hpp"
#include "grid_overlay.hpp"
#include "config_file_handling/config_loader.hpp"
#include "logging/user_mark_log.hpp"
#include "annotation_eraser.h"
#include "annotation_types.hpp"
#include <algorithm>

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
    static const QStringList lables{ "ECG1","ECG2","ECG3","VCG", "PPG","ACCEL","ABP","ART","ART_PULM" };
    return lables;
}

noise_marking_gui::data_channel_features
noise_marking_gui::channelRefs(const QString& label) const {
    auto* self = const_cast<noise_marking_gui*>(this);
    data_channel_features r;

    if (label == "ECG1") {
        r.chartView = ui->ecg_axis_1;
        r.state = &self->mark_state_ecg1;
        r.upsampled_data = &m_ecg1;
        r.dataRaw = &m_ecg1Raw;
        r.sampleRate = channel_upsampled_rates[CH_ECG1];
        r.nativeRate = channel_native_rates[CH_ECG1];
        r.color = COLOR_ECG1;
    }
    else if (label == "ECG2") {
        r.chartView = ui->ecg_axis_2;
        r.state = &self->mark_state_ecg2;
        r.upsampled_data = &m_ecg2;
        r.dataRaw = &m_ecg2Raw;
        r.sampleRate = channel_upsampled_rates[CH_ECG2];
        r.nativeRate = channel_native_rates[CH_ECG2];
        r.color = COLOR_ECG2;
    }
    else if (label == "ECG3") {
        r.chartView = ui->ecg_axis_3;
        r.state = &self->mark_state_ecg3;
        r.upsampled_data = &m_ecg3;
        r.dataRaw = &m_ecg3Raw;
        r.sampleRate = channel_upsampled_rates[CH_ECG3];
        r.nativeRate = channel_native_rates[CH_ECG3];
        r.color = COLOR_ECG3;
    }
    else if (label == "VCG") {
        r.chartView = ui->kors_matrix;
        r.state = &self->mark_state_vcg;
        r.upsampled_data = &m_vcg;
        r.dataRaw = &m_vcgRaw;
        r.sampleRate = channel_upsampled_rates[CH_ECG1];
        r.nativeRate = channel_native_rates[CH_ECG1];
        r.color = COLOR_VCG;
    }
    else if (label == "PPG") {
        r.chartView = ui->ppg_axis;
        r.state = &self->mark_state_ppg;
        r.upsampled_data = &m_ppg;
        r.dataRaw = &m_ppgRaw;
        r.sampleRate = channel_upsampled_rates[CH_PPG];
        r.nativeRate = channel_native_rates[CH_PPG];
        r.color = COLOR_PPG;
    }
    else if (label == "ACCEL") {
        r.chartView = ui->accel_axis;
        r.state = &self->mark_state_accel;
        r.upsampled_data = &m_accelX;
        r.dataRaw = &m_accelXRaw;
        r.sampleRate = channel_upsampled_rates[CH_ACCEL_X];
        r.nativeRate = channel_native_rates[CH_ACCEL_X];
        r.color = COLOR_ACCEL_X;
    }
    else if (label == "ABP") {
        r.chartView = ui->abp_axis;
        r.state = &self->mark_state_abp;
        r.upsampled_data = &m_abp;
        r.dataRaw = &m_abpRaw;
        r.sampleRate = channel_upsampled_rates[CH_ABP];
        r.nativeRate = channel_native_rates[CH_ABP];
        r.color = COLOR_ABP;
    }
    else if (label == "ART") {
        r.chartView = ui->art_axis;
        r.state = &self->mark_state_art;
        r.upsampled_data = &m_art;
        r.dataRaw = &m_artRaw;
        r.sampleRate = channel_upsampled_rates[CH_ART];
        r.nativeRate = channel_native_rates[CH_ART];
        r.color = COLOR_ART;
    }
    else if (label == "ART_PULM") {
        r.chartView = ui->art_pulm_axis;
        r.state = &self->mark_state_art_pulm;
        r.upsampled_data = &m_artPulm;
        r.dataRaw = &m_artPulmRaw;
        r.sampleRate = channel_upsampled_rates[CH_ART_PULM];
        r.nativeRate = channel_native_rates[CH_ART_PULM];
        r.color = COLOR_ART_PULM;
    }
    return r;
}

noise_marking_gui::ChannelMarkingState&
noise_marking_gui::markStateFor(const QString& label) {
    data_channel_features r = channelRefs(label);
    return r.state ? *r.state : mark_state_ppg;
}

QString noise_marking_gui::signalLabelForChartView(QChartView* cv) const {
    if (cv == ui->ecg_axis_1) return "ECG1";
    if (cv == ui->ecg_axis_2) return "ECG2";
    if (cv == ui->ecg_axis_3) return "ECG3";
    if (cv == ui->ppg_axis)   return "PPG";
    if (cv == ui->accel_axis) return "ACCEL";
    if (cv == ui->kors_matrix) return "VCG";
    if (cv == ui->abp_axis) {
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
    else if (label == "PPG") { check = ui->ppg_check;   gain = ui->ppg_gain; }
    else if (label == "ACCEL") { check = ui->accel_check; gain = ui->accel_gain; }
    else if (label == "ABP") { check = ui->abp_check; gain = ui->abp_gain; }
    else if (label == "ART") { check = ui->art_check; gain = ui->art_gain; }
    else if (label == "ART_PULM") { check = ui->art_pulm_check; gain = ui->art_pulm_gain; }
    else if (label == "VCG") { check = ui->kors_check; gain = ui->kors_gain; }
    if (!gain) return 1.0;
    double v = gain->value();
    return (v > 0.0) ? v : 1.0;
}
void noise_marking_gui::onFixScaleToggled(const QString& label, bool on) {
    /*Called when a channel's "Fix Scale" box toggles. On check: capture the
    axis's current range so it can be frozen. On uncheck: drop it so the
    channel autoscales again.
    */
    if (!on) { m_fixedYRange.remove(label); handle_data_plot(); return; }
    QChartView* cv = chartViewForSignalLabel(label);
    if (cv && cv->chart()) {
        auto vAxes = cv->chart()->axes(Qt::Vertical);
        if (!vAxes.isEmpty())
            if (auto* y = qobject_cast<QValueAxis*>(vAxes.first()))
                m_fixedYRange[label] = { y->min(), y->max() };
    }
    handle_data_plot();
}


bool noise_marking_gui::invertedForSignal(const QString& label) const {
    QCheckBox* c = nullptr;
    if (label == "ECG1") c = ui->ecg_1_reverse;
    else if (label == "ECG2") c = ui->ecg_2_reverse;
    else if (label == "ECG3") c = ui->ecg_3_reverse;
    return c && c->isChecked();   // PPG/ABP have no reverse box -> never inverted
}

void noise_marking_gui::autoDetectLeadPolarity() {
    // Measured, not asked. Kirchhoff's voltage law gives II = I + III
    // exactly for genuine limb leads, regardless of any shared DC offset --
    // see vcg::checkLimbLeadPolarity's own header comment in vcg.hpp for the
    // full reasoning (it also covers why a global flip of all three is
    // undetectable here, and why a large residual at every one of the 8 sign
    // combinations means "not limb leads" rather than "still wrong sign").
    const int n = static_cast<int>(std::min({ m_ecg1.size(), m_ecg2.size(), m_ecg3.size() }));
    if (n < 2) return;

    vcg::OrthoAccumulator acc;
    for (int i = 0; i < n; ++i) {
        const double v[3] = { m_ecg1[i], m_ecg2[i], m_ecg3[i] };
        acc.addSample(v);
    }

    const vcg::PolarityCheckResult pc = vcg::checkLimbLeadPolarity(acc);
    if (!pc.consistentWithLimbLeads) {
        // No sign combination brings the residual down -- either these
        // three channels are not limb leads at all, or there is a
        // magnitude/calibration mismatch a sign flip cannot fix either way.
        // Leave the checkboxes exactly as they already were rather than
        // force a "best of 8" answer this measurement cannot actually give.
        return;
    }

    // Pre-set the checkbox to the measured answer. The checkbox stays the
    // visible, overridable source of truth: invertedForSignal() (and every
    // caller of it) still just reads whatever the checkbox says, so the
    // operator can uncheck this immediately if it is ever wrong.
    if (ui->ecg_1_reverse) ui->ecg_1_reverse->setChecked(pc.sign[0] < 0);
    if (ui->ecg_2_reverse) ui->ecg_2_reverse->setChecked(pc.sign[1] < 0);
    if (ui->ecg_3_reverse) ui->ecg_3_reverse->setChecked(pc.sign[2] < 0);
}

void noise_marking_gui::refreshVcgFromLeadFlags() {
    // Fresh every call, per vcg_lead::rebuild()'s own docstring: a box the
    // operator just checked (or unchecked) must withhold (or restore) the
    // VCG on the VERY NEXT rebuild, whether that rebuild was triggered by a
    // chunk load or by this checkbox click directly.
    m_vcgCfg.leadFlaggedInverted[0] = invertedForSignal("ECG1");
    m_vcgCfg.leadFlaggedInverted[1] = invertedForSignal("ECG2");
    m_vcgCfg.leadFlaggedInverted[2] = invertedForSignal("ECG3");
    // rebuild() leaves m_vcg/m_vcgRaw UNTOUCHED on failure (by design, so a
    // failed rebuild can't half-overwrite a good one) -- which would
    // otherwise leave the PREVIOUS (still-checked-out) trace on screen when
    // this call is withheld. Clear first, so "no VCG" actually shows no VCG.
    m_vcg.clear();
    m_vcgRaw.clear();
    vcg_lead::rebuild(m_ecg1Raw, m_ecg2Raw, m_ecg3Raw,
        m_ecg1.size(), m_vcgCfg, m_vcg, m_vcgRaw, m_vcgStatus);

    // markActive's chart-visibility half (loadChunkFromFile's own
    // markActive lambda isn't reachable from here). Without this, the
    // series data gets wiped on the next handle_data_plot(), but the VCG
    // CHART VIEW WIDGET's setVisible() state doesn't move until something
    // else -- a chunk reload, a resize -- happens to touch it, which is
    // exactly the "goes away eventually, needs a resize" symptom. Setting
    // it here, synchronously with the data clear, removes that lag.
    const bool missing = is_missing_signal(m_vcg);
    if (auto* cv = chartViewForSignalLabel("VCG")) cv->setVisible(!missing);
    if (missing) m_activeChannels.remove("VCG");
    else m_activeChannels.insert("VCG");
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
    reset(ui->accel_check, ui->accel_gain);
    reset(ui->abp_check, ui->abp_gain);
    reset(ui->art_check, ui->art_gain);
    reset(ui->art_pulm_check, ui->art_pulm_gain);
    reset(ui->kors_check, ui->kors_gain);
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

    //asdf for DK
    // Letter-key mirrors: a = Left, d = Right (same window scroll);
    new QShortcut(QKeySequence(Qt::Key_D), this, [this]() {
        double maxStart = std::max(0.0, totalChunkDuration() - visible_window_size);
        current_start_time = std::min(current_start_time + scrollStepSeconds(), maxStart);
        resetUnpinnedGains(); handle_data_plot(); updateAmpogramCursor();
        });
    new QShortcut(QKeySequence(Qt::Key_A), this, [this]() {
        current_start_time = std::max(0.0, current_start_time - scrollStepSeconds());
        resetUnpinnedGains(); handle_data_plot(); updateAmpogramCursor();
        });

    const QList<QChartView*> allCharts = {
       ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3,ui->kors_matrix,
       ui->ppg_axis, ui->accel_axis, ui->abp_axis,
       ui->art_axis, ui->art_pulm_axis,
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

    //Toggles the config-driven
    // powerline notch filter (cfg.notch_filter_hz). Disabled entirely if the
    // config has no notch configured (0 Hz), since there'd be nothing to toggle.
    { QSignalBlocker block(ui->notch_filter); ui->notch_filter->setChecked(false); }
    ui->notch_filter->setFocusPolicy(Qt::NoFocus);
    m_notchFilterEnabled = false;
    connect(ui->notch_filter, &QCheckBox::toggled, this, [this](bool on) {
        m_notchFilterEnabled = on;
        loadChunkFromFile(current_chunk_index, /*resetScroll=*/false);
        handle_data_plot();
        updateAmpogramCursor();
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
    wire_gain(ui->accel_check, ui->accel_gain);
    wire_gain(ui->abp_check, ui->abp_gain);
    wire_gain(ui->art_check, ui->art_gain);
    wire_gain(ui->art_pulm_check, ui->art_pulm_gain);
    wire_gain(ui->kors_check, ui->kors_gain);

    auto wire_fix = [this](QCheckBox* check, const QString& label) {
        connect(check, &QCheckBox::toggled, this,
            [this, label](bool on) { onFixScaleToggled(label, on); });
        };
    wire_fix(ui->ecg_1_check, "ECG1");
    wire_fix(ui->ecg_2_check, "ECG2");
    wire_fix(ui->ecg_3_check, "ECG3");
    wire_fix(ui->ppg_check, "PPG");
    wire_fix(ui->accel_check, "ACCEL");
    wire_fix(ui->abp_check, "ABP");
    wire_fix(ui->art_check, "ART");
    wire_fix(ui->art_pulm_check, "ART_PULM");
    wire_fix(ui->kors_check, "VCG");

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
                const double upHz = v->property("upHz").toDouble();
                const double pxPerSec = (visible_window_size > 0.0)
                    ? v->chart()->plotArea().width() / visible_window_size : 0.0;
                const double pxPerSample = (nativeHz > 0.0) ? pxPerSec / nativeHz : 0.0;
                if (v->property("bpm").isValid())
                    v->chart()->setTitle(get_chart_title(sigName, nativeHz, pxPerSample,
                        v->property("bpm").toDouble(), upHz));            // + upHz
                else
                    v->chart()->setTitle(get_chart_title(sigName, nativeHz, pxPerSample, -1.0, upHz));
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

    // THE DROPDOWN IS POPULATED FROM THE TABLE, in table order. The 13 items
    // used to be typed into noise_marking_gui.ui, which made the labels exist in
    // two places that had to be matched by hand -- and matched silently, because
    // find() returns nullptr on a miss and every helper degrades quietly: a
    // retitled item would have exported as code 0, stopped suppressing R-peak
    // detection in its spans, and drawn in the fallback grey, with a missing
    // colour swatch as the only visible symptom.
    //
    // The <item> blocks are gone from the .ui, so this is now the only list. Row
    // order in annotation_types.hpp is therefore the operator's dropdown order,
    // and the table was reordered once to reproduce exactly what the .ui showed
    // (Invert/Noninvert tenth) so nobody's muscle memory moved.
    //
    // setCurrentIndex(0) rather than leaving it: an empty combo has index -1 and
    // currentText() returns an empty string, which m_currentMarkingType would
    // then hold until the operator touched the control.
    {
        QSignalBlocker block(ui->marking_type);
        ui->marking_type->clear();
        for (const auto& t : annotation_types::noise_types)
            ui->marking_type->addItem(QString::fromUtf8(t.label));
        ui->marking_type->setCurrentIndex(0);
    }

    // The swatch loop below is unchanged and still keyed on itemText(i), but it
    // can no longer fail to resolve -- every item came from the table it looks
    // up in. The `if (!t) continue;` is kept as a guard, not as a live path.
    for (int i = 0; i < ui->marking_type->count(); ++i) {
        const QString itemText = ui->marking_type->itemText(i);
        const auto* t = annotation_types::find(itemText);
        if (!t) continue;
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
        connect(c, &QCheckBox::toggled, this, [this](bool) {
            refreshVcgFromLeadFlags();   // recompute -- handle_data_plot() alone only redraws
            handle_data_plot();
            });
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
    //get all markings - both annotations and threshold/invert overrides
    QMap<QString, GenExcStruct> all = m_fileMarkings;
    GenExcStruct current = m_genExc;
    current.filePath = m_binFilePath;
    // marking_type comes from the annotation table, looked up by the FLAG that
    // defines each of these override kinds -- paramEdit for threshold/blanking,
    // invertEdit for inversion. The labels used to be typed here as literals,
    // which made this function a second place the strings lived and a silent one:
    // a retitled table row would leave these spans resolving to nothing and
    // exporting as code 0, with no error anywhere. Resolved at compile time; see
    // the static_asserts in annotation_types.hpp.
    // THE VALUES TRAVEL WITH THE SPAN. applyParamOverrides() writes a threshold
    // and a blanking period for the same extent on the same channel, so the two
    // vectors hold matched entries and one output row can carry both. Pairing
    // them here by (channel, start, end) rather than by index is deliberate:
    // editParamOverrideAt() can replace one vector's entry without the other's,
    // so equal indices are not a safe assumption.
    //
    // m_blankingOverrides was previously never written out at all -- only the
    // threshold spans were, and only as extents. So a reloaded file lost both
    // numbers and every override reverted to the config defaults, moving the R
    // peaks inside it with nothing to indicate why.
    const QString paramLabel = QString::fromUtf8(annotation_types::kParamEditLabel);
    for (const ParamOverride& o : m_thresholdOverrides) {
        double blk = std::numeric_limits<double>::quiet_NaN();
        for (const ParamOverride& b : m_blankingOverrides)
            if (b.channel == o.channel && b.start == o.start && b.end == o.end) {
                blk = b.value; break;
            }
        current.appendMarking(o.start, o.end, o.channel, paramLabel, o.value, blk);
    }
    // A blanking override with no matching threshold override. Should not happen
    // -- applyParamOverrides always writes both -- but a file that ends up with
    // one would otherwise drop the span entirely on save, and losing an operator
    // edit silently is worse than emitting a row whose threshold reads NaN.
    for (const ParamOverride& b : m_blankingOverrides) {
        bool paired = false;
        for (const ParamOverride& o : m_thresholdOverrides)
            if (o.channel == b.channel && o.start == b.start && o.end == b.end) {
                paired = true; break;
            }
        if (!paired)
            current.appendMarking(b.start, b.end, b.channel, paramLabel,
                std::numeric_limits<double>::quiet_NaN(), b.value);
    }
    for (const ParamOverride& o : m_invertOverrides) {
        current.appendMarking(o.start, o.end, o.channel,
            QString::fromUtf8(annotation_types::kInvertEditLabel));
    }

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
        if (m_beatLog) m_beatLog->removeInRange(beat_log::channelForLabel(label), lo, hi);
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

bool noise_marking_gui::invertedAt(const QString& label, double globalTime) const {
    // Only ECG channels have an invert checkbox; PPG/ABP/ART/ART_PULM can't invert.
    if (!(label == "ECG1" || label == "ECG2" || label == "ECG3"))
        return false;
    bool inv = invertedForSignal(label);
    for (const ParamOverride& o : m_invertOverrides)
        if (o.channel == label && globalTime >= o.start && globalTime <= o.end) { inv = !inv; break; }
    return inv;
}

void noise_marking_gui::applyInvertOverride(const QStringList& channels,
    double globalStart, double globalEnd) {
    const double lo = std::min(globalStart, globalEnd);
    const double hi = std::max(globalStart, globalEnd);
    for (const QString& label : channels) {
        m_invertOverrides.append(ParamOverride{ label, lo, hi, 1.0 });
        if (m_beatLog) m_beatLog->removeInRange(beat_log::channelForLabel(label), lo, hi);   // re-detect/re-log in span
    }
    handle_data_plot();
}