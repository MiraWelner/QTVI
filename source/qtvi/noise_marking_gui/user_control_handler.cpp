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
#include <QVector>
#include <QPointF>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
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
    connect(ui->save_current_csv, &QPushButton::clicked, this, &user_control_handler::save_current_csv);
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
void user_control_handler::save_current_plot() {
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
    QDir().mkpath(QString::fromStdString(m_gui->m_cfg.snapshot_path));
    m_gui->grab().save(base + ".png");
}

void user_control_handler::save_current_csv() {
    /*
    Snapshot of the noise-marking GUI's currently visible window: raw and
    upsampled samples for every markable channel, x in milliseconds and
    y in mV. Wide format -- one column set per channel per view. Each
    channel/view has its own (t_ms, mv) pair since raw and upsampled land
    at different sample counts; the file therefore has jagged columns
    (blank cells past the shorter series). Filename matches the PNG
    saved from the same click.
    */
    auto* ui = m_gui->ui.get();
    const QString stem = QFileInfo(m_gui->getFilePath()).completeBaseName();
    const double posSec = m_gui->current_chunk_index * noise_marking_gui::seconds_in_memory_at_once + m_gui->current_start_time;
    const int ti = static_cast<int>(posSec + 0.5);
    const QString stamp = QString("%1h%2m%3s")
        .arg(ti / 3600, 2, 10, QChar('0'))
        .arg((ti % 3600) / 60, 2, 10, QChar('0'))
        .arg(ti % 60, 2, 10, QChar('0'));
    const QString outDir = QString::fromStdString(m_gui->m_cfg.snapshot_path);
    QDir().mkpath(outDir);
    const QString path = outDir + stem + "_" + stamp + ".csv";

    std::ofstream f(path.toStdString());
    if (!f) {
        QMessageBox::warning(m_gui, "Save CSV", "Could not open " + path);
        return;
    }
    f << std::setprecision(10);

    const double t0 = m_gui->current_start_time;
    const double t1 = t0 + m_gui->visible_window_size;

    // Build per-channel/view (t_ms, y_mv) sequences for the visible window.
    struct Series {
        QString colName;
        std::vector<double> t_ms;
        std::vector<double> y_mv;
    };
    std::vector<Series> series;

    for (const QString& label : noise_marking_gui::markableChannelLabels()) {
        auto ref = m_gui->channelRefs(label);

        // Upsampled: index i -> t = i / sampleRate (seconds, chunk-local).
        if (ref.upsampled_data && ref.sampleRate > 0.0) {
            Series s;
            s.colName = label + "_upsampled";
            const QVector<double>& v = *ref.upsampled_data;
            const int i0 = std::max(0, static_cast<int>(std::floor(t0 * ref.sampleRate)));
            const int i1 = std::min(static_cast<int>(v.size()),
                static_cast<int>(std::ceil(t1 * ref.sampleRate)));
            s.t_ms.reserve(i1 - i0);
            s.y_mv.reserve(i1 - i0);
            for (int i = i0; i < i1; ++i) {
                s.t_ms.push_back((static_cast<double>(i) / ref.sampleRate) * 1000.0);
                s.y_mv.push_back(v[i]);
            }
            series.push_back(std::move(s));
        }

        // Raw: (t_sec, v) pairs at native rate, filter to visible window.
        if (ref.dataRaw) {
            Series s;
            s.colName = label + "_raw";
            const QVector<QPointF>& pts = *ref.dataRaw;
            for (const QPointF& p : pts) {
                const double tc = p.x();
                if (tc < t0 || tc > t1) continue;
                s.t_ms.push_back(tc * 1000.0);
                s.y_mv.push_back(p.y());
            }
            series.push_back(std::move(s));
        }
    }

    // Header row: two columns per series, "<name>_t_ms" and "<name>_mv".
    for (size_t k = 0; k < series.size(); ++k) {
        if (k) f << ',';
        f << series[k].colName.toStdString() << "_t_ms,"
            << series[k].colName.toStdString() << "_mv";
    }
    f << '\n';

    // Longest series drives the row count. Shorter series leave blanks.
    size_t maxRows = 0;
    for (const Series& s : series) maxRows = std::max(maxRows, s.t_ms.size());

    for (size_t row = 0; row < maxRows; ++row) {
        for (size_t k = 0; k < series.size(); ++k) {
            if (k) f << ',';
            const Series& s = series[k];
            if (row < s.t_ms.size()) {
                f << s.t_ms[row] << ',' << s.y_mv[row];
            }
            else {
                f << ',';   // blank t_ms, blank y_mv
            }
        }
        f << '\n';
    }

    f.close();
    std::cout << "Saved: " << path.toStdString() << std::endl;
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