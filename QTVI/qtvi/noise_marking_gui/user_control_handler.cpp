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

    // Top-row actions
    connect(ui->undo_button, &QPushButton::clicked, this, &user_control_handler::handle_undo_button);
    connect(ui->clearall_button, &QPushButton::clicked, this, &user_control_handler::handle_clearall_button);
    connect(ui->finalize_button, &QPushButton::clicked, this, &user_control_handler::handle_finalize_button);
    connect(ui->skip_button, &QPushButton::clicked, this, &user_control_handler::handle_skip_button);
    connect(ui->save_current_plot, &QPushButton::clicked, this, &user_control_handler::save_current_plot);

    // Per-channel marking start/stop. Every channel follows the same pattern
    // a dozen near-identical connect() lines.
    struct MarkButtons {
        QPushButton* start;
        QPushButton* stop;
        const char* signalName;
    };
    const MarkButtons markButtons[] = {
        { ui->start_ecg1_mark, ui->stop_ecg1_mark, "ECG1" },
        { ui->start_ecg2_mark, ui->stop_ecg2_mark, "ECG2" },
        { ui->start_ecg3_mark, ui->stop_ecg3_mark, "ECG3" },
        { ui->startNoisePPG,   ui->stopNoisePPG,   "PPG"  },
        { ui->startNoiseABP,   ui->stopNoiseABP,   "ABP"  },
    };
    for (const auto& mb : markButtons) {
        const QString name = QString::fromLatin1(mb.signalName);
        connect(mb.start, &QPushButton::clicked, this, [this, name]() { m_gui->beginMarking(name); });
        connect(mb.stop, &QPushButton::clicked, this, [this, name]() { m_gui->beginStopPhase(name); });
    }

    // "Mark all signals"
    connect(ui->start_all_mark, &QPushButton::clicked, this, &user_control_handler::handle_allmarkingstart_button);
    connect(ui->stop_all_mark, &QPushButton::clicked, this, &user_control_handler::handle_allmarkingstop_button);

    // Window-length selector: combo index -> seconds. Default is index 2 (10 s).
    static constexpr double window_length_options[] = { 1, 3, 10, 30, 60, 120,300 };
    connect(ui->window_length_selector,
        QOverload<int>::of(&QComboBox::currentIndexChanged),
        this,
        [this](int idx) { handle_window_toggle(true, window_length_options[idx]); });
    ui->window_length_selector->setCurrentIndex(2);

    //threshold and blanking buttons allow you to select regions
    connect(ui->change_threshold, &QPushButton::clicked, this,
        [this] { m_gui->enterParamEdit(noise_marking_gui::ParamEdit::Threshold); });
    connect(ui->change_blanking, &QPushButton::clicked, this,
        [this] { m_gui->enterParamEdit(noise_marking_gui::ParamEdit::Blanking); });
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
    const double posSec = m_gui->m_currentChunkIndex * noise_marking_gui::seconds_in_memory_at_once
        + m_gui->m_currentStartTime;
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

    const double t0 = m_gui->m_currentStartTime;
    const double t1 = t0 + m_gui->m_windowDuration;

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
    const double dtUp = 1.0 / sr;
    for (double t = t0; t < t1; t += dtUp) {
        out << QString::number(t, 'f', 6);
        const int idx = static_cast<int>(std::round(t * sr));
        for (const auto& c : channels) {
            out << "," << QString::number((*c.up)[idx], 'g', 8);
        }
        out << "\n";
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
    m_gui->m_windowDuration = duration;
    m_gui->handle_data_plot();
}

// ============================================================================
// Undo / Clear
// ============================================================================

// Rebuild the NoiseManager from the surviving exception list. Undo pops the
// last entry and reconstructs, which is simpler than tracking incremental
// inverse operations and matches the "exception list is source of truth"
// invariant the rest of the GUI relies on.
void user_control_handler::handle_undo_button() {
    auto& exc = m_gui->m_genExc;
    if (exc.noiseExc.isEmpty()) return;

    exc.noiseExc.removeLast();
    exc.data_type.removeLast();
    exc.marking_type.removeLast();

    m_gui->m_noiseManager = std::make_unique<annotation_handler>(m_gui->m_ecgSR);
    for (int i = 0; i < exc.noiseExc.size(); ++i) {
        const double sr = m_gui->sampleRateForSignal(exc.data_type[i]);
        m_gui->m_noiseManager->addSegment(
            static_cast<int>(exc.noiseExc[i].first * sr),
            static_cast<int>(exc.noiseExc[i].second * sr),
            exc.data_type[i].toStdString(),
            exc.marking_type[i].toStdString());
    }
    m_gui->updateNoiseHighlights();
}

void user_control_handler::handle_clearall_button() {
    if (QMessageBox::question(m_gui, "Clear", "Clear all markings?") != QMessageBox::Yes)
        return;

    auto& exc = m_gui->m_genExc;
    exc.noiseExc.clear();
    exc.data_type.clear();
    exc.marking_type.clear();
    m_gui->m_noiseManager = std::make_unique<annotation_handler>(m_gui->m_ecgSR);
    m_gui->updateNoiseHighlights();
}

// ============================================================================
// Mark-all passthroughs
// ============================================================================

void user_control_handler::handle_allmarkingstart_button() { m_gui->beginMarkingAll(); }
void user_control_handler::handle_allmarkingstop_button() { m_gui->beginStopPhaseAll(); }