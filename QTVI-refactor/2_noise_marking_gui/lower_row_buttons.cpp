// lower_row_buttons.cpp - CORRECTED VERSION
/******************************************************************************
* File:    lower_row_buttons.cpp
* Author:  Mira Welner
* Date:    2026-01-20
* Version: 1.0
*
* Description:
*		The Qt GUI which is created by main.cpp has a lower row of buttons which include:
*			1) Undo - undoes last annotation
*			2) Clear All - asks for conformation from user, then clears all annotations
*			3) Skip - if the loaded file is all noise, or has some other problem, the user can skip it in which case
*			   the GUI closes and opens again having loaded the next file, and the csv file and bin file corresponding
*              to the previous signal are NOT recorded
*			4) Save - once the annotations are completed, they are saved in a CSV file which lists the start and
*              end locations in the time series (a discrete value), the start and end seconds based on the reported
*              sampling rate, the type (ECG vs PPG) and the marking_type (Noise/Arifact, PVC, etc)
*           5) 
*
*
*****************************************************************************/

#include "lower_row_buttons.hpp"
#include "noise_marking_gui.h"
#include "ui_noise_marking_gui.h"
#include <QMessageBox>


lower_row_buttons::lower_row_buttons(noise_marking_gui* parent)
    : QObject(parent)
    , m_gui(parent)
{
}

void lower_row_buttons::setupConnections() {
	// the Qt UI is in the noise_marking_gui class, so we connect the buttons there to their handlers here
    connect(m_gui->ui->undo_button, &QPushButton::clicked, this, &lower_row_buttons::handle_undo_button);
    connect(m_gui->ui->clearall_button, &QPushButton::clicked, this, &lower_row_buttons::handle_clearall_button);
    connect(m_gui->ui->finalize_button, &QPushButton::clicked, this, &lower_row_buttons::handle_finalize_button);
    connect(m_gui->ui->skip_button, &QPushButton::clicked, this, &lower_row_buttons::handle_skip_button);

    connect(m_gui->ui->startNoiseECG, &QPushButton::clicked, this, &lower_row_buttons::handle_ecgmarkingstart_button);
    connect(m_gui->ui->stopNoiseECG, &QPushButton::clicked, this, &lower_row_buttons::handle_ecgmarkingstop_button);
    connect(m_gui->ui->startNoisePPG, &QPushButton::clicked, this, &lower_row_buttons::handle_ppgmarkingstart_button);
    connect(m_gui->ui->stopNoisePPG, &QPushButton::clicked, this, &lower_row_buttons::handle_ppgmarkingstop_button);

    connect(m_gui->ui->rb_10s, &QRadioButton::toggled, this, &lower_row_buttons::handle_10s_window_toggle);
    connect(m_gui->ui->rb_30s, &QRadioButton::toggled, this, &lower_row_buttons::handle_30s_window_toggle);
    connect(m_gui->ui->rb_1m, &QRadioButton::toggled, this, &lower_row_buttons::handle_1m_window_toggle);
    connect(m_gui->ui->rb_5m, &QRadioButton::toggled, this, &lower_row_buttons::handle_5m_window_toggle);
    connect(m_gui->ui->rb_10m, &QRadioButton::toggled, this, &lower_row_buttons::handle_10m_window_toggle);
}


void lower_row_buttons::handle_finalize_button() { m_gui->accept(); } //accept is a QT native function that closes the dialog with an "Accepted" result
void lower_row_buttons::handle_skip_button() { m_gui->reject(); }//reject is a QT native function that closes the dialog with an "Rejected" result

//When the window size radio buttons are toggled, update the window duration and re-plot the data
void lower_row_buttons::handle_10s_window_toggle(bool c) { if (c) { m_gui->m_windowDuration = 10; m_gui->handle_data_plot(); } }
void lower_row_buttons::handle_30s_window_toggle(bool c) { if (c) { m_gui->m_windowDuration = 30; m_gui->handle_data_plot(); } }
void lower_row_buttons::handle_1m_window_toggle(bool c) { if (c) { m_gui->m_windowDuration = 60; m_gui->handle_data_plot(); } }
void lower_row_buttons::handle_5m_window_toggle(bool c) { if (c) { m_gui->m_windowDuration = 300; m_gui->handle_data_plot(); } }
void lower_row_buttons::handle_10m_window_toggle(bool c) { if (c) { m_gui->m_windowDuration = 600; m_gui->handle_data_plot(); } }

void lower_row_buttons::handle_undo_button() {
    // Accessing private m_genExc via the friend relationship
    if (!m_gui->m_genExc.noiseExc.isEmpty()) {
        m_gui->m_genExc.noiseExc.removeLast();
        m_gui->m_genExc.data_type.removeLast();
        m_gui->m_genExc.marking_type.removeLast();

        // Reset and rebuild noise manager
        m_gui->m_noiseManager = std::make_unique<NoiseManager>(m_gui->m_ecgSR);
        for (int i = 0; i < m_gui->m_genExc.noiseExc.size(); ++i) {
            m_gui->m_noiseManager->addSegment(
                m_gui->m_genExc.noiseExc[i].first * m_gui->m_ecgSR,
                m_gui->m_genExc.noiseExc[i].second * m_gui->m_ecgSR,
                m_gui->m_genExc.data_type[i].toStdString(),
                m_gui->m_genExc.marking_type[i].toStdString()
            );
        }
        m_gui->updateNoiseHighlights();
    }
}

void lower_row_buttons::handle_clearall_button() {
    if (QMessageBox::question(m_gui, "Clear", "Clear all markings?") == QMessageBox::Yes) {
        m_gui->m_genExc.noiseExc.clear();
        m_gui->m_genExc.data_type.clear();
        m_gui->m_genExc.marking_type.clear();
        m_gui->m_noiseManager = std::make_unique<NoiseManager>(m_gui->m_ecgSR);
        m_gui->updateNoiseHighlights();
    }
}

void lower_row_buttons::handle_ecgmarkingstart_button() {
    // If already active (waiting for start, waiting for end, or marker exists), toggle off
    if (m_gui->m_isWaitingForECGStart || m_gui->m_isWaitingForECGEnd || m_gui->m_ecgStartMarkerLine) {
        m_gui->start_marking_button_clicked(true);
    }
    else {
        m_gui->m_currentMarkingType = m_gui->ui->marking_type->currentText();
        m_gui->m_isWaitingForECGStart = true;
        m_gui->ui->startNoiseECG->setStyleSheet("background-color: #f39c12; color: white;");
    }
}

void lower_row_buttons::handle_ppgmarkingstart_button() {
    if (m_gui->m_isWaitingForPPGStart || m_gui->m_isWaitingForPPGEnd || m_gui->m_ppgStartMarkerLine) {
        m_gui->start_marking_button_clicked(false);
    }
    else {
        m_gui->m_currentMarkingType = m_gui->ui->marking_type->currentText();
        m_gui->m_isWaitingForPPGStart = true;
        m_gui->ui->startNoisePPG->setStyleSheet("background-color: #f39c12; color: white;");
    }
}

void lower_row_buttons::handle_ecgmarkingstop_button() {
    if (m_gui->m_ecgStartMarkerLine) {
        m_gui->m_isWaitingForECGEnd = true;
        m_gui->m_isWaitingForECGStart = false;
        m_gui->ui->stopNoiseECG->setStyleSheet("background-color: #e74c3c; color: white;");
        m_gui->ui->startNoiseECG->setStyleSheet("");
    }
}

void lower_row_buttons::handle_ppgmarkingstop_button() {
    if (m_gui->m_ppgStartMarkerLine) {
        m_gui->m_isWaitingForPPGEnd = true;
        m_gui->m_isWaitingForPPGStart = false;
        m_gui->ui->stopNoisePPG->setStyleSheet("background-color: #e74c3c; color: white;");
        m_gui->ui->startNoisePPG->setStyleSheet("");
    }
}
