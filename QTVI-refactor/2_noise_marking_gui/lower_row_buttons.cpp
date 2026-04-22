/**
 * @file   lower_row_buttons.cpp
 * @brief  Handles the lower toolbar of the noise-marking GUI.
 *
 * @author Mira Welner
 * @date   2026-01-20
 */
#include "lower_row_buttons.h"
#include "gui_handler.h"
#include "ui_noise_marking_gui.h"
#include <QMessageBox>

lower_row_buttons::lower_row_buttons(noise_marking_gui* parent)
    : QObject(parent), m_gui(parent) {
}

void lower_row_buttons::setupConnections() {
    auto* ui = m_gui->ui.get();

    connect(ui->undo_button, &QPushButton::clicked, this, &lower_row_buttons::handle_undo_button);
    connect(ui->clearall_button, &QPushButton::clicked, this, &lower_row_buttons::handle_clearall_button);
    connect(ui->finalize_button, &QPushButton::clicked, this, &lower_row_buttons::handle_finalize_button);
    connect(ui->skip_button, &QPushButton::clicked, this, &lower_row_buttons::handle_skip_button);

    // Per-channel ECG buttons
    connect(ui->start_ecg1_mark, &QPushButton::clicked, this, [this]() { m_gui->beginMarking("ECG1"); });
    connect(ui->stop_ecg1_mark, &QPushButton::clicked, this, [this]() { m_gui->beginStopPhase("ECG1"); });
    connect(ui->start_ecg2_mark, &QPushButton::clicked, this, [this]() { m_gui->beginMarking("ECG2"); });
    connect(ui->stop_ecg2_mark, &QPushButton::clicked, this, [this]() { m_gui->beginStopPhase("ECG2"); });
    connect(ui->start_ecg3_mark, &QPushButton::clicked, this, [this]() { m_gui->beginMarking("ECG3"); });
    connect(ui->stop_ecg3_mark, &QPushButton::clicked, this, [this]() { m_gui->beginStopPhase("ECG3"); });

    // PPG buttons
    connect(ui->startNoisePPG, &QPushButton::clicked, this, &lower_row_buttons::handle_ppgmarkingstart_button);
    connect(ui->stopNoisePPG, &QPushButton::clicked, this, &lower_row_buttons::handle_ppgmarkingstop_button);

    // ABP buttons
    connect(ui->startNoiseABP, &QPushButton::clicked, this,
        [this]() { m_gui->beginMarking("ABP"); });
    connect(ui->stopNoiseABP, &QPushButton::clicked, this,
        [this]() { m_gui->beginStopPhase("ABP"); });

    // Mark All Signals buttons
    connect(ui->start_all_mark, &QPushButton::clicked, this, &lower_row_buttons::handle_allmarkingstart_button);
    connect(ui->stop_all_mark, &QPushButton::clicked, this, &lower_row_buttons::handle_allmarkingstop_button);

    // Window radio buttons
    connect(ui->rb_1s, &QRadioButton::toggled, this, [this](bool c) { handle_window_toggle(c, 1);   });
    connect(ui->rb_3s, &QRadioButton::toggled, this, [this](bool c) { handle_window_toggle(c, 3);   });
    connect(ui->rb_10s, &QRadioButton::toggled, this, [this](bool c) { handle_window_toggle(c, 10);  });
    connect(ui->rb_30s, &QRadioButton::toggled, this, [this](bool c) { handle_window_toggle(c, 30);  });
    connect(ui->rb_1m, &QRadioButton::toggled, this, [this](bool c) { handle_window_toggle(c, 60);  });
    connect(ui->rb_10m, &QRadioButton::toggled, this, [this](bool c) { handle_window_toggle(c, 600); });
}

// ============================================================================
// Dialog actions
// ============================================================================

void lower_row_buttons::handle_finalize_button() { m_gui->accept(); }
void lower_row_buttons::handle_skip_button() { m_gui->reject(); }

void lower_row_buttons::handle_window_toggle(bool checked, double duration) {
    if (!checked) return;
    m_gui->m_windowDuration = duration;
    m_gui->handle_data_plot();
}

// ============================================================================
// Undo / Clear
// ============================================================================

void lower_row_buttons::handle_undo_button() {
    if (m_gui->m_genExc.noiseExc.isEmpty()) return;

    m_gui->m_genExc.noiseExc.removeLast();
    m_gui->m_genExc.data_type.removeLast();
    m_gui->m_genExc.marking_type.removeLast();

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
    if (QMessageBox::question(m_gui, "Clear", "Clear all markings?") != QMessageBox::Yes)
        return;

    m_gui->m_genExc.noiseExc.clear();
    m_gui->m_genExc.data_type.clear();
    m_gui->m_genExc.marking_type.clear();
    m_gui->m_noiseManager = std::make_unique<NoiseManager>(m_gui->m_ecgSR);
    m_gui->updateNoiseHighlights();
}

// ============================================================================
// Mark All Signals start / stop
// ============================================================================

void lower_row_buttons::handle_allmarkingstart_button() {
    m_gui->beginMarkingAll();
}

void lower_row_buttons::handle_allmarkingstop_button() {
    m_gui->beginStopPhaseAll();
}

// ============================================================================
// PPG start / stop
// ============================================================================

void lower_row_buttons::handle_ppgmarkingstart_button() {
    m_gui->beginMarking("PPG");
}

void lower_row_buttons::handle_ppgmarkingstop_button() {
    m_gui->beginStopPhase("PPG");
}