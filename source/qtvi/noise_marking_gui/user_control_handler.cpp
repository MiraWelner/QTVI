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
    connect(ui->mark_one_chan, &QPushButton::clicked, this, [this] { m_gui->toggleEcgMark(); });
    connect(ui->mark_all_chan, &QPushButton::clicked, this, [this] { m_gui->toggleMarkAll(); });
    connect(ui->mark_all_ecg, &QPushButton::clicked, this, [this] { m_gui->toggleMarkEcgAll(); });

    // Window-length selector: combo index -> seconds. Default is index 2 (10 s).
    static constexpr double window_length_options[] = { 1, 3, 10, 30, 60, 120,300 };
    connect(ui->window_length_selector,
        QOverload<int>::of(&QComboBox::currentIndexChanged),
        this,
        [this](int idx) { handle_window_toggle(true, window_length_options[idx]); });
    ui->window_length_selector->setCurrentIndex(3);

    //threshold and blanking buttons allow you to select regions
    connect(ui->param_change, &QPushButton::clicked, this,
        [this] { m_gui->enterParamEdit(); });
}

// ============================================================================
// Dialog actions
// ============================================================================

void user_control_handler::handle_finalize_button() { m_gui->accept(); }
void user_control_handler::handle_skip_button() { m_gui->reject(); }

// Export the currently visible signal window to CSV. The file has two
// blocks: an UPSAMPLED block (one row per 1 kHz sample, one column per
// channel) and a RAW block in long format (channel, time, value), since
// raw channels have different native sample rates and don't share a grid.
//
// Path:  <output_path>/saved_plot_snapshots/<bin_stem>_<timestamp>.csv
//        Timestamp is local time, yyyyMMdd_hhmmss, so each click produces
//        a fresh file rather than overwriting.
void user_control_handler::save_current_plot()
{
    const QString stem = QFileInfo(m_gui->getFilePath()).completeBaseName();
    // In-recording position of the visible window's start, filename-safe
    // (no colons). 1 h 23 m 45 s into the recording -> "01h23m45s".
    const double posSec = m_gui->current_chunk_index * noise_marking_gui::seconds_in_memory_at_once
        + m_gui->current_start_time;
    const int ti = static_cast<int>(posSec + 0.5);
    const QString stamp = QString("%1h%2m%3s")
        .arg(ti / 3600, 2, 10, QChar('0'))
        .arg((ti % 3600) / 60, 2, 10, QChar('0'))
        .arg(ti % 60, 2, 10, QChar('0'));
    const QString base = QString::fromStdString(m_gui->m_cfg.snapshot_path) + stem + "_" + stamp;

    // Screenshot of the whole dialog, same stem + timestamp as the CSV.
    m_gui->grab().save(base + ".png");

    const QString path = base + ".csv";
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return;
    QTextStream out(&f);

    const double t0 = m_gui->current_start_time;
    const double t1 = t0 + m_gui->visible_window_size;

    struct Ch {
        const char* label;
        const QVector<double>* up;
        const QVector<QPointF>* raw;
    };
    const Ch channels[] = {
        { "ecg1",   &m_gui->m_ecg1,   &m_gui->m_ecg1Raw   },
        { "ecg2",   &m_gui->m_ecg2,   &m_gui->m_ecg2Raw   },
        { "ecg3",   &m_gui->m_ecg3,   &m_gui->m_ecg3Raw   },
        { "ppg",    &m_gui->m_ppg,    &m_gui->m_ppgRaw    },
        { "abp",    &m_gui->m_abp,    &m_gui->m_abpRaw    },
        { "accelx", &m_gui->m_accelX, &m_gui->m_accelXRaw },
        { "accely", &m_gui->m_accelY, &m_gui->m_accelYRaw },
        { "accelz", &m_gui->m_accelZ, &m_gui->m_accelZRaw },
        { "resp",   &m_gui->m_resp,   &m_gui->m_respRaw   },
        { "cvp",    &m_gui->m_cvp,    &m_gui->m_cvpRaw    },
    };

    // Upsampled block: shared 1 kHz time column + per-channel value columns.
    out << "# UPSAMPLED\ntime_s";
    for (const auto& c : channels) out << "," << c.label;
    out << "\n";

    const double sr = m_gui->m_ecgSR;
    if (sr > 0.0) {
        const double dtUp = 1.0 / sr;
        for (double t = t0; t < t1; t += dtUp) {
            out << QString::number(t, 'f', 6);
            const int idx = static_cast<int>(std::round(t * sr));
            for (const auto& c : channels) {
                const QVector<double>& up = *c.up;
                // Missing channels are empty / a single sentinel, and the last
                // window can run past the loaded samples -- both would index out
                // of bounds. Write a blank cell instead of reading past the end.
                if (idx >= 0 && idx < up.size())
                    out << "," << QString::number(up[idx], 'g', 8);
                else
                    out << ",";
            }
            out << "\n";
        }
    }

    // Raw block: long format (channel, time, value). Raw vectors are
    // monotone in time, so once we pass t1 we can stop scanning.
    out << "\n# RAW\nchannel,time_s,value\n";
    for (const auto& c : channels) {
        for (const QPointF& p : *c.raw) {
            if (p.x() < t0) continue;
            if (p.x() > t1) break;
            out << c.label << ","
                << QString::number(p.x(), 'f', 6) << ","
                << QString::number(p.y(), 'g', 8) << "\n";
        }
    }
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
    m_gui->m_noiseManager = std::make_unique<annotation_handler>(m_gui->m_ecgSR);
    m_gui->handle_data_plot();
}

// ============================================================================
// Mark-all passthroughs
// ============================================================================

void user_control_handler::handle_allmarkingstart_button() { m_gui->beginMarkingAll(); }
void user_control_handler::handle_allmarkingstop_button() { m_gui->beginStopPhaseAll(); }