/**
 * @file   user_control_handler.cpp
 * @brief  Implementation of the noise-marking GUI's bottom toolbar.
 *
 * @author Mira Welner
 * @date   2026-01-20
 */
#include "user_control_handler.h"
#include "gui_handler.h"
#include "ui_noise_marking_gui.h"

#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QMessageBox>
#include <QTextStream>
#include <QPixmap>
#include <cmath>
#include <QRadioButton>

user_control_handler::user_control_handler(noise_marking_gui* parent)
    : QObject(parent), m_gui(parent) {
}

// ============================================================================
// Wiring
// ============================================================================

void user_control_handler::setupConnections() {
    auto* ui = m_gui->ui.get();
    connect(ui->clearall_button, &QPushButton::clicked, this, &user_control_handler::handle_clearall_button);
    connect(ui->finalize_button, &QPushButton::clicked, this, &user_control_handler::handle_finalize_button);
    connect(ui->skip_button, &QPushButton::clicked, this, &user_control_handler::handle_skip_button);
    connect(ui->save_current_plot, &QPushButton::clicked, this, &user_control_handler::save_current_plot);
    connect(ui->mark_one_chan, &QRadioButton::toggled, this,
        [this](bool on) { if (on) m_gui->setMarkScope(noise_marking_gui::MarkScope::One); });
    connect(ui->mark_ecg, &QRadioButton::toggled, this,
        [this](bool on) { if (on) m_gui->setMarkScope(noise_marking_gui::MarkScope::Ecg); });
    connect(ui->mark_all_chan, &QRadioButton::toggled, this,
        [this](bool on) { if (on) m_gui->setMarkScope(noise_marking_gui::MarkScope::All); });
    connect(ui->make_annotation, &QPushButton::clicked, this,
        [this] { m_gui->toggleAnnotationArm(); });

    // Window-length selector: combo index -> seconds. Default is index 2 (10 s).
    static constexpr double window_length_options[] = { 1, 3, 10, 30, 60, 120,300 };
    connect(ui->window_length_selector,
        QOverload<int>::of(&QComboBox::currentIndexChanged),
        this,
        [this](int idx) { handle_window_toggle(true, window_length_options[idx]); });
    ui->window_length_selector->setCurrentIndex(3);
}

// ============================================================================
// Dialog actions
// ============================================================================

void user_control_handler::handle_finalize_button() { m_gui->accept(); }
void user_control_handler::handle_skip_button() { m_gui->reject(); }
void user_control_handler::save_current_plot(){
    /*
    Save a .png screenshot of the GUI at this moment
    */
    const QString stem = QFileInfo(m_gui->getFilePath()).completeBaseName();
    const double posSec = m_gui->current_chunk_index * noise_marking_gui::seconds_in_memory_at_once + m_gui->current_start_time;
    const int ti = static_cast<int>(posSec + 0.5);
    const QString stamp = QString("%1h%2m%3s")
        .arg(ti / 3600, 2, 10, QChar('0'))
        .arg((ti % 3600) / 60, 2, 10, QChar('0'))
        .arg(ti % 60, 2, 10, QChar('0'));
    const QString base = QString::fromStdString(m_gui->m_cfg.snapshot_path) + stem + "_" + stamp;
    m_gui->grab().save(base + ".png");
}

void user_control_handler::handle_window_toggle(bool checked, double duration) {
    if (!checked) return;
    m_gui->visible_window_size = duration;
    m_gui->handle_data_plot();
}

void user_control_handler::handle_clearall_button() {
    /*
		When you click the 'clear all' button a confirmation dialog pops up to prevent accidental clearing of markings.
        If the user confirms, all markings are cleared from the current file and the GUI is updated to reflect this change.
    */
    if (QMessageBox::question(m_gui, "Clear", "Clear all markings?") != QMessageBox::Yes)
        return;

    auto& exc = m_gui->m_genExc;
    exc.noiseExc.clear();
    exc.data_type.clear();
    exc.marking_type.clear();
    m_gui->m_noiseManager = std::make_unique<annotation_handler>();
    m_gui->handle_data_plot();
}