/**
 * @file   lower_row_buttons.cpp
 * @author Mira Welner
 * @date   2026-01-20
 *
 * @brief  Handles the lower toolbar of the noise-marking GUI:
 *           1) Undo    - removes the last annotation
 *           2) Clear   - asks for confirmation, then removes all annotations
 *           3) Skip    - closes the dialog without saving (rejected)
 *           4) Save    - closes the dialog and saves annotations (accepted)
 *           5) ECG start/stop - uses ecg_channel_selector to pick ECG1/2/3
 *           6) PPG start/stop - always marks PPG
 *           7) Window-size radio buttons (10s, 30s, 1m, 5m, 10m)
 */

#include "lower_row_buttons.h"
#include "gui_handler.h"
#include "ui_noise_marking_gui.h"
#include <QMessageBox>

lower_row_buttons::lower_row_buttons(noise_marking_gui* parent)
    : QObject(parent), m_gui(parent) {
}

void lower_row_buttons::setupConnections() {
    connect(m_gui->ui->undo_button, &QPushButton::clicked, this, &lower_row_buttons::handle_undo_button);
    connect(m_gui->ui->clearall_button, &QPushButton::clicked, this, &lower_row_buttons::handle_clearall_button);
    connect(m_gui->ui->finalize_button, &QPushButton::clicked, this, &lower_row_buttons::handle_finalize_button);
    connect(m_gui->ui->skip_button, &QPushButton::clicked, this, &lower_row_buttons::handle_skip_button);

    connect(m_gui->ui->startNoiseECG, &QPushButton::clicked, this, &lower_row_buttons::handle_ecgmarkingstart_button);
    connect(m_gui->ui->stopNoiseECG, &QPushButton::clicked, this, &lower_row_buttons::handle_ecgmarkingstop_button);
    connect(m_gui->ui->startNoisePPG, &QPushButton::clicked, this, &lower_row_buttons::handle_ppgmarkingstart_button);
    connect(m_gui->ui->stopNoisePPG, &QPushButton::clicked, this, &lower_row_buttons::handle_ppgmarkingstop_button);

    connect(m_gui->ui->rb_10s, &QRadioButton::toggled, this, [this](bool c) { handle_window_toggle(c, 10); });
    connect(m_gui->ui->rb_30s, &QRadioButton::toggled, this, [this](bool c) { handle_window_toggle(c, 30); });
    connect(m_gui->ui->rb_1m, &QRadioButton::toggled, this, [this](bool c) { handle_window_toggle(c, 60); });
    connect(m_gui->ui->rb_5m, &QRadioButton::toggled, this, [this](bool c) { handle_window_toggle(c, 300); });
    connect(m_gui->ui->rb_10m, &QRadioButton::toggled, this, [this](bool c) { handle_window_toggle(c, 600); });
}

void lower_row_buttons::handle_finalize_button() { m_gui->accept(); }
void lower_row_buttons::handle_skip_button() { m_gui->reject(); }

void lower_row_buttons::handle_window_toggle(bool checked, double duration) {
    if (checked) {
        m_gui->m_windowDuration = duration;
        m_gui->handle_data_plot();
    }
}

void lower_row_buttons::handle_undo_button() {
    if (m_gui->m_genExc.noiseExc.isEmpty()) return;

    m_gui->m_genExc.noiseExc.removeLast();
    m_gui->m_genExc.data_type.removeLast();
    m_gui->m_genExc.marking_type.removeLast();

    // Rebuild noise manager from remaining annotations
    m_gui->m_noiseManager = std::make_unique<NoiseManager>(m_gui->m_ecgSR);
    for (int i = 0; i < m_gui->m_genExc.noiseExc.size(); ++i) {
        double sr = m_gui->sampleRateForSignal(m_gui->m_genExc.data_type[i]);
        m_gui->m_noiseManager->addSegment(
            static_cast<size_t>(m_gui->m_genExc.noiseExc[i].first * sr),
            static_cast<size_t>(m_gui->m_genExc.noiseExc[i].second * sr),
            m_gui->m_genExc.data_type[i].toStdString(),
            m_gui->m_genExc.marking_type[i].toStdString()
        );
    }
    m_gui->updateNoiseHighlights();
}

void lower_row_buttons::handle_clearall_button() {
    if (QMessageBox::question(m_gui, "Clear", "Clear all markings?") != QMessageBox::Yes) return;

    m_gui->m_genExc.noiseExc.clear();
    m_gui->m_genExc.data_type.clear();
    m_gui->m_genExc.marking_type.clear();
    m_gui->m_noiseManager = std::make_unique<NoiseManager>(m_gui->m_ecgSR);
    m_gui->updateNoiseHighlights();
}

// ============================================================================
// ECG start/stop — channel comes from ecg_channel_selector
// ============================================================================

void lower_row_buttons::handle_ecgmarkingstart_button() {
    QString label = m_gui->selectedEcgLabel();  // "ECG1", "ECG2", or "ECG3"
    auto& state = m_gui->markStateFor(label);

    // If already in marking mode for this channel, cancel it
    if (state.isWaitingForStart || state.isWaitingForEnd || state.startMarkerLine) {
        m_gui->start_marking_button_clicked(label);
    }
    else {
        m_gui->m_currentMarkingType = m_gui->ui->marking_type->currentText();
        state.isWaitingForStart = true;
        m_gui->ui->startNoiseECG->setStyleSheet("background-color: #f39c12; color: white;");
    }
}

void lower_row_buttons::handle_ecgmarkingstop_button() {
    QString label = m_gui->selectedEcgLabel();
    auto& state = m_gui->markStateFor(label);
    if (!state.startMarkerLine) return;

    state.isWaitingForEnd = true;
    state.isWaitingForStart = false;
    m_gui->ui->stopNoiseECG->setStyleSheet("background-color: #e74c3c; color: white;");
    m_gui->ui->startNoiseECG->setStyleSheet("");
}

// ============================================================================
// PPG start/stop — always "PPG"
// ============================================================================

void lower_row_buttons::handle_ppgmarkingstart_button() {
    auto& state = m_gui->markStateFor("PPG");

    if (state.isWaitingForStart || state.isWaitingForEnd || state.startMarkerLine) {
        m_gui->start_marking_button_clicked("PPG");
    }
    else {
        m_gui->m_currentMarkingType = m_gui->ui->marking_type->currentText();
        state.isWaitingForStart = true;
        m_gui->ui->startNoisePPG->setStyleSheet("background-color: #f39c12; color: white;");
    }
}

void lower_row_buttons::handle_ppgmarkingstop_button() {
    auto& state = m_gui->markStateFor("PPG");
    if (!state.startMarkerLine) return;

    state.isWaitingForEnd = true;
    state.isWaitingForStart = false;
    m_gui->ui->stopNoisePPG->setStyleSheet("background-color: #e74c3c; color: white;");
    m_gui->ui->startNoisePPG->setStyleSheet("");
}