#include "TemplateViewerWindow.hpp"
#include "ui_TemplateViewerWindow.h"
#include "feature_marks.hpp"
#include "NormalizeFeatures.hpp"
#include <QMessageBox>
#include <QColor>
#include <QTimer>
#include <QPixmap>
#include <QDir>
#include <QFileInfo>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <QFile>
#include <iomanip>
#include <iostream>
#include <cstdio>

// ========================================================================
// Construction
// ========================================================================

TemplateViewerWindow::TemplateViewerWindow(QWidget* parent)
    : QMainWindow(parent), ui(new Ui::TemplateViewerWindow)
{
    ui->setupUi(this);

    auto* moveGroup = new QButtonGroup(this);
    moveGroup->addButton(ui->MoveIndividual);
    moveGroup->addButton(ui->MoveSubsequent);
    moveGroup->setExclusive(true);

    connect(ui->MoveSubsequent, &QRadioButton::toggled, this,
        [this](bool on) { m_moveSubsequent = on; });


    m_moveSubsequent = ui->MoveSubsequent->isChecked();

    // Startup defaults, enforced in code (independent of the .ui's checked
    // attributes): all markers OFF, all waveform traces ON. setChecked here
    // updates the visible checkboxes so they reflect the actual state; the
    // toggled() connections below keep the flags in sync thereafter.
    ui->show_ecg_markers->setChecked(false);
    ui->show_ppg_markers->setChecked(false);
    ui->show_abp_markers->setChecked(false);
    ui->show_art_markers->setChecked(false);
    ui->show_art_pulm_markers->setChecked(false);
    ui->show_ecg->setChecked(true);
    ui->show_ppg->setChecked(true);
    ui->show_abp->setChecked(true);
    ui->show_art->setChecked(true);
    ui->show_art_pulm->setChecked(true);

    m_showEcgMarkers = false;
    m_showPpgMarkers = false;
    m_showAbpMarkers = false;
    m_showArtMarkers = false;
    m_showArtPulmMarkers = false;
    m_showEcgTrace = true;
    m_showPpgTrace = true;
    m_showAbpTrace = true;
    m_showArtTrace = true;
    m_showArtPulmTrace = true;

    connect(ui->show_ecg_markers, &QCheckBox::toggled, this, [this](bool on) {
        m_showEcgMarkers = on; applyMarkerVisibility();
        });
    connect(ui->show_ppg_markers, &QCheckBox::toggled, this, [this](bool on) {
        m_showPpgMarkers = on; applyMarkerVisibility();
        });
    connect(ui->show_abp_markers, &QCheckBox::toggled, this, [this](bool on) {
        m_showAbpMarkers = on; applyMarkerVisibility();
        });
    connect(ui->show_art_markers, &QCheckBox::toggled, this, [this](bool on) {
        m_showArtMarkers = on; applyMarkerVisibility();
        });
    connect(ui->show_art_pulm_markers, &QCheckBox::toggled, this, [this](bool on) {
        m_showArtPulmMarkers = on; applyMarkerVisibility();
        });

    connect(ui->show_ecg, &QCheckBox::toggled, this, [this](bool on) {
        m_showEcgTrace = on; applyMarkerVisibility();
        });
    connect(ui->show_ppg, &QCheckBox::toggled, this, [this](bool on) {
        m_showPpgTrace = on; applyMarkerVisibility();
        });
    connect(ui->show_abp, &QCheckBox::toggled, this, [this](bool on) {
        m_showAbpTrace = on; applyMarkerVisibility();
        });
    connect(ui->show_art, &QCheckBox::toggled, this, [this](bool on) {
        m_showArtTrace = on; applyMarkerVisibility();
        });
    connect(ui->show_art_pulm, &QCheckBox::toggled, this, [this](bool on) {
        m_showArtPulmTrace = on; applyMarkerVisibility();
        });
}

TemplateViewerWindow::~TemplateViewerWindow() { delete ui; }

// ========================================================================
// Helpers
// ========================================================================

std::vector<TemplateViewerWindow::Lead>
TemplateViewerWindow::leadsForBin(const TemplateBin& b) const {
    std::vector<Lead> out;
    if (!b.ch1.ecgTemplate_raw.empty())
        out.push_back({ &b.ch1.ecgTemplate_raw, 0, "Ch1" });
    if (!b.ch2.ecgTemplate_raw.empty())
        out.push_back({ &b.ch2.ecgTemplate_raw, 1, "Ch2" });
    if (!b.ch3.ecgTemplate_raw.empty())
        out.push_back({ &b.ch3.ecgTemplate_raw, 2, "Ch3" });
    return out;
}

std::pair<int, int> TemplateViewerWindow::compactGrid(int n) {
    if (n <= 0) return { 1, 1 };
    int cols = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(n))));
    int rows = static_cast<int>(std::ceil(static_cast<double>(n) / cols));
    if (rows > cols) std::swap(rows, cols);
    return { rows, cols };
}

void TemplateViewerWindow::updatePageControls() {
    ui->pageLabel->setText(QString("%1 / %2").arg(m_currentPage + 1).arg(m_totalPages));
    ui->prevButton->setEnabled(m_currentPage > 0);
    ui->nextButton->setEnabled(m_currentPage < m_totalPages - 1);
}

void TemplateViewerWindow::onNextPage() {
    if (m_currentPage < m_totalPages - 1) {
        ++m_currentPage;
        showPage();
    }
}

void TemplateViewerWindow::onPrevPage() {
    if (m_currentPage > 0) {
        --m_currentPage;
        showPage();
    }
}

// ========================================================================
// Load subject
// ========================================================================

void TemplateViewerWindow::loadSubject(const QString& templatePath, const QString& markingPath, const QString& subjectId, double sampleRateHz) {

    m_markingPath = markingPath;
    m_templateDir = QFileInfo(templatePath).absolutePath();
    m_subjectId = subjectId;
    m_sampleRate = sampleRateHz;
    setWindowTitle(QString("Template Marking - %1").arg(subjectId));
    ui->subjectLabel->setText(subjectId);
    // First pass shows "Finish and Next"; after the Q-align reload it reads
    // "Finish" (m_qAlignPass is set in save_bin_and_csv and preserved across
    // the reload).
    ui->finishButton->setText(m_qAlignPass ? "Finish" : "Finish and Next");

    try {
        m_bins = readTemplateInfoBin(templatePath.toStdString());
    }
    catch (const std::exception& e) {
        QMessageBox::critical(this, "Read error",
            QString("Failed to read %1:\n\n%2").arg(templatePath, e.what()));
        m_bins.clear();
        emit finished();
        return;
    }

    if (m_bins.empty()) {
        QMessageBox::warning(this, "Error", "No bins loaded from " + templatePath);
        emit finished();
        return;
    }

    m_maxLeads = 1;
    for (const auto& b : m_bins) {
        int nl = (int)leadsForBin(b).size();
        if (nl > m_maxLeads) m_maxLeads = nl;
    }
    m_binsPerPage = (m_maxLeads <= 1) ? 16 : 4;

    m_currentPage = 0;
    m_totalPages = std::max(1, static_cast<int>(
        std::ceil(static_cast<double>(m_bins.size()) / m_binsPerPage)));

    // Seed markers for every bin so per-subject global refs (which need
    // R/S and foot/peak positions across all bins) can be computed once
    // and stay stable across paging.
    for (auto& b : m_bins) FeatureMarks::seed_all(b, m_sampleRate);
    computeGlobalRefs();

    m_capturedPages.clear();   // new subject: nothing saved yet
    showPage();
    // Defer the analysis CSV write until after the widget page has painted.
    // Writing it inline blocks the UI thread for ~15s on large templates
    // (7 signals x ~1300 samples x ~800 bins = ~7M formatted writes), which
    // looks like the viewer is frozen at load time. QTimer::singleShot(0)
    // schedules it on the next event-loop cycle so the page shows first.
    QTimer::singleShot(0, this, [this]() { writeAlignedTemplateCsv(); });
}

// Compute per-subject Global_Ref for every ECG channel + every pulse
// channel. Called at loadSubject after every bin has been seeded, so
// R/S positions and pulse foot/peak positions are already set. Results
// go into m_ecgGlobalRef and m_pulseGlobalRef and stay stable across
// pages. See normalize_features.hpp for the definitions.
void TemplateViewerWindow::computeGlobalRefs() {
    for (int c = 0; c < 3; ++c)
        m_ecgGlobalRef[c] = normalize_features::compute_ecg_global_ref(
            m_bins, c, m_sampleRate);
    for (int c = 0; c < 4; ++c)
        m_pulseGlobalRef[c] = normalize_features::compute_pulse_global_ref(
            m_bins, c);
    fprintf(stderr,
        "[refs] ecg=[%.4g,%.4g,%.4g] pulse=[%.4g,%.4g,%.4g,%.4g]\n",
        m_ecgGlobalRef[0], m_ecgGlobalRef[1], m_ecgGlobalRef[2],
        m_pulseGlobalRef[0], m_pulseGlobalRef[1],
        m_pulseGlobalRef[2], m_pulseGlobalRef[3]);
}

std::vector<double>
TemplateViewerWindow::normalizeEcgTrace(const std::vector<double>& raw, int ch) const {
    const double ref = (ch >= 0 && ch < 3) ? m_ecgGlobalRef[ch] : std::nan("");
    if (!std::isfinite(ref) || ref == 0.0) return raw;   // untransformed
    std::vector<double> out(raw.size());
    for (size_t i = 0; i < raw.size(); ++i)
        out[i] = std::isnan(raw[i]) ? raw[i] : raw[i] / ref;
    return out;
}

std::vector<double>
TemplateViewerWindow::normalizePulseTrace(const std::vector<double>& raw,
    int footIdx, int pulseChan) const {
    const double ref = (pulseChan >= 0 && pulseChan < 4)
        ? m_pulseGlobalRef[pulseChan] : std::nan("");
    if (!std::isfinite(ref) || ref == 0.0) return raw;
    if (footIdx < 0 || footIdx >= static_cast<int>(raw.size())) return raw;
    const double footY = raw[footIdx];
    if (std::isnan(footY) || std::abs(footY) < 1e-12) return raw;
    std::vector<double> out(raw.size());
    for (size_t i = 0; i < raw.size(); ++i) {
        if (std::isnan(raw[i])) { out[i] = raw[i]; continue; }
        const double localRatio = 100.0 * (raw[i] - footY) / footY;
        out[i] = localRatio / ref;
    }
    return out;
}

// ========================================================================
// Build / clear plots
// ========================================================================

void TemplateViewerWindow::clearPlots() {
    for (auto* pw : m_allPlots) {
        ui->plotGrid->removeWidget(pw);
        delete pw;
    }
    m_allPlots.clear();
    m_binPlots.clear();
    m_pageGlobalIdx.clear();

    // Drop stretch factors left over from a previous (possibly larger) page
    // so unused rows/columns don't reserve empty space on the next page.
    for (int c = 0; c < ui->plotGrid->columnCount(); ++c)
        ui->plotGrid->setColumnStretch(c, 0);
    for (int r = 0; r < ui->plotGrid->rowCount(); ++r)
        ui->plotGrid->setRowStretch(r, 0);
}

void TemplateViewerWindow::showPage() {
    clearPlots();

    int start = m_currentPage * m_binsPerPage;
    int end = std::min(start + m_binsPerPage, static_cast<int>(m_bins.size()));
    int count = end - start;

    bool compact = (m_maxLeads <= 1);

    int gridRows;
    if (compact) {
        auto [r, c] = compactGrid(count);
        (void)r;
        gridRows = c;
    }
    else {
        gridRows = m_maxLeads;
    }

    m_binPlots.resize(count);
    m_pageGlobalIdx.resize(count);

    int usedRows = 0, usedCols = 0;

    for (int i = 0; i < count; ++i) {
        int gi = start + i;
        m_pageGlobalIdx[i] = gi;

        const TemplateBin& b = m_bins[gi];
        auto leads = leadsForBin(b);
        bool hasPPG = !b.ppgTemplate.empty();

        if (leads.empty())
            leads.push_back({ nullptr, 0, "No ECG" });

        std::vector<BinPlotWidget*> group;

        for (int li = 0; li < (int)leads.size(); ++li) {
            auto* pw = new BinPlotWidget(gi, leads[li].channelIndex, leads[li].label);
            pw->setSampleRate(m_sampleRate);

            static const std::vector<double> empty;
            const auto& ecg = leads[li].ecg ? *leads[li].ecg : empty;
            const auto& ppg = hasPPG ? b.ppgTemplate : empty;

            int c = leads[li].channelIndex;
            // R column = argmax of the ECG template itself. Self-locating,
            // NaN-skipping. Works whether the template is R-anchored by
            // pad-based slicing (Patch B/C) or by explicit alignment; no
            // dependency on avg_r_expand.
            double rPeak = 0.0;
            {
                double best = -std::numeric_limits<double>::infinity();
                for (int i = 0; i < static_cast<int>(ecg.size()); ++i) {
                    if (!std::isnan(ecg[i]) && ecg[i] > best) {
                        best = ecg[i]; rPeak = static_cast<double>(i);
                    }
                }
            }

            // std vectors for this channel + PPG. Empty if the templater
            // didn't compute them for this bin -- the widget treats empty
            // std as "no band, just the line".
            const auto& ecgStd = (c == 0) ? b.ch1.ecgTemplate_raw_std
                : (c == 1) ? b.ch2.ecgTemplate_raw_std
                : b.ch3.ecgTemplate_raw_std;
            const auto& ppgStd = hasPPG ? b.ppgTemplate_std : empty;

            // Per-channel beat count for this widget's lead + the bin's
            // PPG count. Both are 0 when the channel is absent, which
            // suppresses that half of the title suffix.
            const uint64_t nEcgBeats = (c == 0) ? b.ch1_n_beats_raw
                : (c == 1) ? b.ch2_n_beats_raw
                : b.ch3_n_beats_raw;

            // Per-subject normalized traces for on-screen display.
            // ECG:   sample / Global_Ref_ecg(ch)
            // Pulse: ((sample - foot_y) / |foot_y| * 100) / Global_Ref_pulse(chan)
            // If the ref or foot is unusable, the helpers pass the raw
            // trace through unchanged.
            const std::vector<double> ecgN = normalizeEcgTrace(ecg, c);
            const std::vector<double> ppgN = hasPPG
                ? normalizePulseTrace(b.ppgTemplate, b.ppg_onset, 0)
                : empty;
            const std::vector<double> abpN = !b.abpTemplate.empty()
                ? normalizePulseTrace(b.abpTemplate, b.abp_onset, 1)
                : b.abpTemplate;
            const std::vector<double> artN = !b.artTemplate.empty()
                ? normalizePulseTrace(b.artTemplate, b.art_onset, 2)
                : b.artTemplate;
            const std::vector<double> artPN = !b.artPulmTemplate.empty()
                ? normalizePulseTrace(b.artPulmTemplate, b.art_pulm_onset, 3)
                : b.artPulmTemplate;

            pw->setData(ppgN, ppgStd, ecgN, ecgStd,
                b.p_peak_ch[c], b.q_begin_ch[c], b.r_peak_ch[c],
                b.s_end_ch[c], b.t_begin_ch[c], b.t_end_ch[c],
                b.ppg_onset, b.ppg_p50, b.ppg_peak,
                b.ppg_dicrotic, b.ppg_peak2, b.ppg_t80, b.ppg_end,
                rPeak,
                static_cast<int>(nEcgBeats),
                static_cast<int>(b.ppg_n_beats));
            pw->setHasPPG(hasPPG);

            // Faint arterial background-context traces (present-only),
            // foot-anchored like the PPG. Colors mirror the noise-marking GUI:
            // ABP teal, ART dark red, ART_PULM dark blue.
            {
                std::vector<std::pair<std::vector<double>, QColor>> bg;
                if (!abpN.empty())
                    bg.push_back({ abpN,  QColor(0, 95, 105) });
                if (!artN.empty())
                    bg.push_back({ artN,  QColor(150, 40, 40) });
                if (!artPN.empty())
                    bg.push_back({ artPN, QColor(40, 60, 150) });
                pw->setBackgroundTraces(bg);
            }

            // Arterial traces for marker geometry/bounds, plus the markers
            // themselves (shared across leads, like PPG).
            pw->setArterialTraces(abpN, artN, artPN,
                b.abpTemplate_std, b.artTemplate_std, b.artPulmTemplate_std);
            pw->setMarker(BinPlotWidget::AbpOnset, b.abp_onset);
            pw->setMarker(BinPlotWidget::AbpPeak, b.abp_peak);
            pw->setMarker(BinPlotWidget::AbpDicrotic, b.abp_dicrotic);
            pw->setMarker(BinPlotWidget::AbpPeak2, b.abp_peak2);
            pw->setMarker(BinPlotWidget::AbpEnd, b.abp_end);
            pw->setMarker(BinPlotWidget::ArtOnset, b.art_onset);
            pw->setMarker(BinPlotWidget::ArtPeak, b.art_peak);
            pw->setMarker(BinPlotWidget::ArtDicrotic, b.art_dicrotic);
            pw->setMarker(BinPlotWidget::ArtPeak2, b.art_peak2);
            pw->setMarker(BinPlotWidget::ArtEnd, b.art_end);
            pw->setMarker(BinPlotWidget::ArtPulmOnset, b.art_pulm_onset);
            pw->setMarker(BinPlotWidget::ArtPulmPeak, b.art_pulm_peak);
            pw->setMarker(BinPlotWidget::ArtPulmDicrotic, b.art_pulm_dicrotic);
            pw->setMarker(BinPlotWidget::ArtPulmPeak2, b.art_pulm_peak2);
            pw->setMarker(BinPlotWidget::ArtPulmEnd, b.art_pulm_end);

            if (b.ppg_issue == 1)
                pw->setState(BinPlotWidget::State::BadPPG);
            else if (b.bad_r_ch[c])
                pw->setState(BinPlotWidget::State::BadR);

            connect(pw, &BinPlotWidget::markerMoved,
                this, &TemplateViewerWindow::onMarkerMoved);
            connect(pw, &BinPlotWidget::markerDragStarted,
                this, &TemplateViewerWindow::onMarkerDragStarted);
            connect(pw, &BinPlotWidget::badRToggled,
                this, &TemplateViewerWindow::onBadRToggled);
            connect(pw, &BinPlotWidget::badPPGToggled,
                this, &TemplateViewerWindow::onBadPPGToggled);

            if (compact) {
                int row = i % gridRows;
                int col = i / gridRows;
                ui->plotGrid->addWidget(pw, row, col);
                usedRows = std::max(usedRows, row + 1);
                usedCols = std::max(usedCols, col + 1);
            }
            else {
                int rowspan = (li == (int)leads.size() - 1)
                    ? (gridRows - li) : 1;
                ui->plotGrid->addWidget(pw, li, i, rowspan, 1);
                usedRows = std::max(usedRows, li + rowspan);
                usedCols = std::max(usedCols, i + 1);
            }

            m_allPlots.push_back(pw);
            group.push_back(pw);
        }

        m_binPlots[i] = std::move(group);
    }

    // Equal stretch on every used row/column => equal-width, equal-height
    // cells that together fill the whole plot area. Combined with each
    // widget scaling its trace to its cell width, every bin window ends up
    // the same on-screen length.
    for (int c = 0; c < usedCols; ++c) ui->plotGrid->setColumnStretch(c, 1);
    for (int r = 0; r < usedRows; ++r) ui->plotGrid->setRowStretch(r, 1);

    applyMarkerVisibility();
    updatePageControls();

    // Save a screenshot of this page (once) now that it's on screen. Fires on
    // initial load (page 1) and whenever the user scrolls to a new page.
    captureCurrentPage();
}

void TemplateViewerWindow::captureCurrentPage() {
    if (m_capturedPages.count(m_currentPage)) return;   // already saved this page
    QDir outDir(m_templateDir);
    const int page = m_currentPage;
    QTimer::singleShot(60, this, [this, page, outDir]() {
        if (m_capturedPages.count(page)) return;

        const QPixmap shot = ui->scrollContents->grab();
        const QString fn = outDir.filePath(
            QString("%1_templates_page%2.png")
            .arg(m_subjectId)
            .arg(page + 1, 2, 10, QChar('0')));
        const bool ok = shot.save(fn, "PNG");
        if (ok) m_capturedPages.insert(page);
        });
}

// Merge the r_align and q_align aligned-template CSVs into one wide CSV.
// Shared keys (file_id,bin_num,x_ms) appear once; every other column is
// emitted as a pair <name>_r,<name>_q. Rows are paired by line index (both
// passes share bin order and per-bin template length); any length mismatch
// is blank-padded. Returns true on success.
namespace {
    bool mergeAlignedPassCsvs(const std::string& rPath, const std::string& qPath,
        const std::string& outPath) {
        std::ifstream rf(rPath), qf(qPath);
        if (!rf || !qf) return false;
        auto readAll = [](std::ifstream& in) {
            std::vector<std::string> lines; std::string ln;
            while (std::getline(in, ln)) {
                if (!ln.empty() && ln.back() == '\r') ln.pop_back();
                lines.push_back(ln);
            }
            return lines;
            };
        const std::vector<std::string> R = readAll(rf), Q = readAll(qf);
        if (R.empty() || Q.empty()) return false;
        auto split = [](const std::string& s) {
            std::vector<std::string> c; std::string cell; std::stringstream ss(s);
            while (std::getline(ss, cell, ',')) c.push_back(cell);
            return c;
            };
        std::ofstream out(outPath);
        if (!out) return false;
        const std::vector<std::string> qh = split(Q[0]);
        const size_t ncol = qh.size();
        out << "file_id,bin_num,x_ms";
        for (size_t i = 3; i < ncol; ++i) out << ',' << qh[i] << "_r_align" << ',' << qh[i] << "_q_align";
        out << '\n';
        const size_t n = std::max(R.size(), Q.size());
        for (size_t li = 1; li < n; ++li) {
            const std::vector<std::string> rc = (li < R.size()) ? split(R[li]) : std::vector<std::string>();
            const std::vector<std::string> qc = (li < Q.size()) ? split(Q[li]) : std::vector<std::string>();
            const std::vector<std::string>& key = !qc.empty() ? qc : rc;
            out << (key.size() > 0 ? key[0] : "") << ','
                << (key.size() > 1 ? key[1] : "") << ','
                << (key.size() > 2 ? key[2] : "");
            for (size_t i = 3; i < ncol; ++i) {
                out << ',' << (i < rc.size() ? rc[i] : "")
                    << ',' << (i < qc.size() ? qc[i] : "");
            }
            out << '\n';
        }
        return true;
    }
}   // namespace

void TemplateViewerWindow::writeAlignedTemplateCsv() {
    if (m_bins.empty()) return;

    QDir outDir(QDir(m_templateDir).absoluteFilePath("../csv_for_analysis"));
    if (!outDir.exists()) outDir.mkpath(".");
    const QString pass = m_qAlignPass ? "q_align" : "r_align";
    const QString path = outDir.filePath(m_subjectId + "_template_" + pass + ".csv");

    std::ofstream f(path.toStdString());
    if (!f) { fprintf(stderr, "[tmplcsv] cannot open %s\n", path.toStdString().c_str()); return; }

    const double toMs = (m_sampleRate > 0.0) ? 1000.0 / m_sampleRate : 1.0;

    // ---- Header ------------------------------------------------------------
    // Per channel: raw_mv, Normalized_mv, raw_std, normalized_std.
    static const char* CHANS[] = {
        "ch1", "ch2", "ch3", "ppg", "abp", "art", "art_pulm"
    };
    f << "file_id,bin_num,x_ms";
    for (const char* n : CHANS) {
        f << ',' << n << "_raw_mv"
            << ',' << n << "_Normalized"
            << ',' << n << "_raw_std"
            << ',' << n << "_normalized_std";
    }

    // Marker-location columns. For each marker: two columns (autodetect,
    // user). A row has "1" in the column iff its row index equals that
    // marker's sample index for the current bin; blank otherwise.
    //
    // ECG markers per channel (8 each x 3 channels x 2 auto/user = 48):
    //   p_peak, q_begin, q_peak (computed), r_peak, s_peak (computed),
    //   s_end, t_peak, t_end.
    static const char* ECG_MARKERS[] = {
        "p_peak", "q_begin", "q_peak", "r_peak", "s_peak",
        "s_end",  "t_begin",  "t_end"
    };
    for (int c = 1; c <= 3; ++c) {
        for (int k = 0; k < 8; ++k) {
            const char* mname = ECG_MARKERS[k];
            const bool userToo = (k != 3);   // r_peak = auto-only
            f << ',' << mname << "_ch" << c << "_location_autodetect";
            if (userToo) f << ',' << mname << "_ch" << c << "_location_user";
        }
    }
    // Pulse markers (PPG has 6 with p50, arterial has 5).
    static const char* PPG_MARKERS[] = {
        "ppg_onset", "ppg_p50", "ppg_peak",
        "ppg_dicrotic", "ppg_peak2", "ppg_t80", "ppg_end"
    };
    static const char* ABP_MARKERS[] = {
        "abp_onset", "abp_peak", "abp_dicrotic", "abp_peak2", "abp_end"
    };
    static const char* ART_MARKERS[] = {
        "art_onset", "art_peak", "art_dicrotic", "art_peak2", "art_end"
    };
    static const char* ARTP_MARKERS[] = {
        "art_pulm_onset", "art_pulm_peak", "art_pulm_dicrotic",
        "art_pulm_peak2", "art_pulm_end"
    };
    auto emitPulseHeaderGroup = [&](auto const& group) {
        for (const char* m : group) {
            f << ',' << m << "_location_autodetect"
                << ',' << m << "_location_user";
        }
        };
    emitPulseHeaderGroup(PPG_MARKERS);
    emitPulseHeaderGroup(ABP_MARKERS);
    emitPulseHeaderGroup(ART_MARKERS);
    emitPulseHeaderGroup(ARTP_MARKERS);
    // Autodetected computed feature locations (no user bar; from AUTODETECT markers).
    auto emitAutoFeatLocHeader = [&](const char* n) { f << ',' << n << "_autodetect_location"; };
    for (int gc = 1; gc <= 3; ++gc) {
        char gb[64];
        for (const char* g : { "p_wave", "q_onset", "r_wave", "t_peak" }) {
            std::snprintf(gb, sizeof gb, "%s_ch%d", g, gc);
            emitAutoFeatLocHeader(gb);
        }
    }
    for (const char* g : { "ppg_foot", "ppg_p1" })
        emitAutoFeatLocHeader(g);
    f << '\n';
    f << std::setprecision(10);

    // Scale a per-sample std trace by the same linear factor the mean
    // transform applies (see the pulse / ECG normalization rules).
    auto scaleEcgStd = [&](const std::vector<double>& raw_std, int chan) {
        std::vector<double> out = raw_std;
        const double ref = (chan >= 0 && chan < 3) ? m_ecgGlobalRef[chan] : 0.0;
        if (!std::isfinite(ref) || ref == 0.0) return out;
        for (double& s : out) if (!std::isnan(s)) s /= ref;
        return out;
        };
    auto scalePulseStd = [&](const std::vector<double>& raw_std,
        int footIdx, const std::vector<double>& mean_trace, int pulseChan)
        {
            std::vector<double> out = raw_std;
            const double ref = (pulseChan >= 0 && pulseChan < 4) ? m_pulseGlobalRef[pulseChan] : 0.0;
            if (!std::isfinite(ref) || ref == 0.0) return out;
            if (footIdx < 0 || footIdx >= (int)mean_trace.size()) return out;
            const double fy = mean_trace[footIdx];
            if (std::isnan(fy) || std::abs(fy) < 1e-12) return out;
            const double k = 100.0 / std::abs(fy) / std::abs(ref);
            for (double& s : out) if (!std::isnan(s)) s *= k;
            return out;
        };

    // ---- Row loop ----------------------------------------------------------
    for (size_t bi = 0; bi < m_bins.size(); ++bi) {
        const TemplateBin& b = m_bins[bi];

        // Raw + std traces per channel.
        const std::vector<double>& ch1R = b.ch1.ecgTemplate_raw;
        const std::vector<double>& ch2R = b.ch2.ecgTemplate_raw;
        const std::vector<double>& ch3R = b.ch3.ecgTemplate_raw;
        const std::vector<double>& ppgR = b.ppgTemplate;
        const std::vector<double>& abpR = b.abpTemplate;
        const std::vector<double>& artR = b.artTemplate;
        const std::vector<double>& artPR = b.artPulmTemplate;

        const std::vector<double>& ch1S = b.ch1.ecgTemplate_raw_std;
        const std::vector<double>& ch2S = b.ch2.ecgTemplate_raw_std;
        const std::vector<double>& ch3S = b.ch3.ecgTemplate_raw_std;
        const std::vector<double>& ppgS = b.ppgTemplate_std;
        const std::vector<double>& abpS = b.abpTemplate_std;
        const std::vector<double>& artS = b.artTemplate_std;
        const std::vector<double>& artPS = b.artPulmTemplate_std;

        // Normalized mean traces (screen-matching).
        std::vector<double> ch1N = normalizeEcgTrace(ch1R, 0);
        std::vector<double> ch2N = normalizeEcgTrace(ch2R, 1);
        std::vector<double> ch3N = normalizeEcgTrace(ch3R, 2);
        std::vector<double> ppgN = normalizePulseTrace(ppgR, b.ppg_onset, 0);
        std::vector<double> abpN = normalizePulseTrace(abpR, b.abp_onset, 1);
        std::vector<double> artN = normalizePulseTrace(artR, b.art_onset, 2);
        std::vector<double> artPN = normalizePulseTrace(artPR, b.art_pulm_onset, 3);

        // Normalized std traces.
        std::vector<double> ch1NS = scaleEcgStd(ch1S, 0);
        std::vector<double> ch2NS = scaleEcgStd(ch2S, 1);
        std::vector<double> ch3NS = scaleEcgStd(ch3S, 2);
        std::vector<double> ppgNS = scalePulseStd(ppgS, b.ppg_onset, ppgR, 0);
        std::vector<double> abpNS = scalePulseStd(abpS, b.abp_onset, abpR, 1);
        std::vector<double> artNS = scalePulseStd(artS, b.art_onset, artR, 2);
        std::vector<double> artPNS = scalePulseStd(artPS, b.art_pulm_onset, artPR, 3);

        // Computed Q/S peaks (both variants) per ECG channel.
        const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
        EcgFeatures ftAuto[3], ftUser[3];
        for (int c = 0; c < 3; ++c) {
            const auto& ecg = chs[c]->ecgTemplate_raw;
            ftAuto[c] = computeEcgFeatures(ecg,
                b.p_peak_auto_ch[c], b.q_begin_auto_ch[c], b.r_peak_auto_ch[c],
                b.s_end_auto_ch[c], b.t_begin_auto_ch[c], b.t_end_auto_ch[c],
                m_sampleRate);
            ftUser[c] = computeEcgFeatures(ecg,
                b.p_peak_ch[c], b.q_begin_ch[c], b.r_peak_ch[c],
                b.s_end_ch[c], b.t_begin_ch[c], b.t_end_ch[c],
                m_sampleRate);
        }

        // ECG marker positions per channel, aligned with ECG_MARKERS order:
        //   p_peak, q_begin, q_peak(computed), r_peak, s_peak(computed),
        //   s_end,  t_peak, t_end.
        int ecgAuto[3][8], ecgUser[3][8];
        for (int c = 0; c < 3; ++c) {
            ecgAuto[c][0] = b.p_peak_auto_ch[c];
            ecgAuto[c][1] = b.q_begin_auto_ch[c];
            ecgAuto[c][2] = ftAuto[c].q_idx;
            ecgAuto[c][3] = b.r_peak_auto_ch[c];
            ecgAuto[c][4] = ftAuto[c].s_idx;
            ecgAuto[c][5] = b.s_end_auto_ch[c];
            ecgAuto[c][6] = b.t_begin_auto_ch[c];
            ecgAuto[c][7] = b.t_end_auto_ch[c];
            ecgUser[c][0] = b.p_peak_ch[c];
            ecgUser[c][1] = b.q_begin_ch[c];
            ecgUser[c][2] = ftUser[c].q_idx;
            ecgUser[c][3] = b.r_peak_ch[c];
            ecgUser[c][4] = ftUser[c].s_idx;
            ecgUser[c][5] = b.s_end_ch[c];
            ecgUser[c][6] = b.t_begin_ch[c];
            ecgUser[c][7] = b.t_end_ch[c];
        }

        // Pulse marker positions, order matching PPG_MARKERS / ABP_MARKERS etc.
        const int ppgAuto[7] = { b.ppg_onset_auto, b.ppg_p50_auto, b.ppg_peak_auto,
                                 b.ppg_dicrotic_auto, b.ppg_peak2_auto, b.ppg_t80_auto, b.ppg_end_auto };
        const int ppgUser[7] = { b.ppg_onset, b.ppg_p50, b.ppg_peak,
                                 b.ppg_dicrotic, b.ppg_peak2, b.ppg_t80, b.ppg_end };
        const int abpAuto[5] = { b.abp_onset_auto, b.abp_peak_auto,
                                 b.abp_dicrotic_auto, b.abp_peak2_auto, b.abp_end_auto };
        const int abpUser[5] = { b.abp_onset, b.abp_peak,
                                 b.abp_dicrotic, b.abp_peak2, b.abp_end };
        const int artAuto[5] = { b.art_onset_auto, b.art_peak_auto,
                                 b.art_dicrotic_auto, b.art_peak2_auto, b.art_end_auto };
        const int artUser[5] = { b.art_onset, b.art_peak,
                                 b.art_dicrotic, b.art_peak2, b.art_end };
        const int artpAuto[5] = { b.art_pulm_onset_auto, b.art_pulm_peak_auto,
                                  b.art_pulm_dicrotic_auto, b.art_pulm_peak2_auto,
                                  b.art_pulm_end_auto };
        const int artpUser[5] = { b.art_pulm_onset, b.art_pulm_peak,
                                  b.art_pulm_dicrotic, b.art_pulm_peak2, b.art_pulm_end };

        // Row span = longest trace (all start at row 0).
        auto Nof = [](const std::vector<double>& v) { return (int)v.size(); };
        const int hiRow = std::max({ Nof(ch1R), Nof(ch2R), Nof(ch3R),
            Nof(ppgR), Nof(abpR), Nof(artR), Nof(artPR) });
        if (hiRow <= 0) continue;

        struct Col {
            const std::vector<double>* raw;
            const std::vector<double>* norm;
            const std::vector<double>* raw_std;
            const std::vector<double>* norm_std;
        };
        const Col cols[] = {
            { &ch1R,  &ch1N,  &ch1S,  &ch1NS  },
            { &ch2R,  &ch2N,  &ch2S,  &ch2NS  },
            { &ch3R,  &ch3N,  &ch3S,  &ch3NS  },
            { &ppgR,  &ppgN,  &ppgS,  &ppgNS  },
            { &abpR,  &abpN,  &abpS,  &abpNS  },
            { &artR,  &artN,  &artS,  &artNS  },
            { &artPR, &artPN, &artPS, &artPNS },
        };

        auto emitVal = [&](const std::vector<double>& v, int j) {
            f << ',';
            if (j >= 0 && j < (int)v.size() && !std::isnan(v[j])) f << v[j];
            };
        // Emit ",1" if row equals marker index, ",<blank>" otherwise.
        // Autodetected computed feature indices (no user bar; from AUTODETECT markers).
        const ChannelTemplateData* gchs[3] = { &b.ch1, &b.ch2, &b.ch3 };
        FeatureMarks::EcgGlyphs egl[3];
        for (int gc = 0; gc < 3; ++gc)
            egl[gc] = FeatureMarks::compute_ecg_glyphs(
                gchs[gc]->ecgTemplate_raw, b.p_peak_auto_ch[gc], b.q_begin_auto_ch[gc],
                b.s_end_auto_ch[gc], b.t_begin_auto_ch[gc], b.t_end_auto_ch[gc], m_sampleRate);
        const FeatureMarks::PpgGlyphs pgl = FeatureMarks::compute_ppg_glyphs(
            b.ppgTemplate, b.ppg_onset_auto, b.ppg_dicrotic_auto, b.ppg_peak2_auto);
        auto emitLoc = [&](int markerIdx, int row) {
            f << ',';
            if (markerIdx >= 0 && markerIdx == row) f << '1';
            };

        for (int row = 0; row < hiRow; ++row) {
            f << m_subjectId.toStdString() << ',' << bi << ','
                << (row * toMs);
            for (const Col& c : cols) {
                emitVal(*c.raw, row);
                emitVal(*c.norm, row);
                emitVal(*c.raw_std, row);
                emitVal(*c.norm_std, row);
            }
            // ECG location columns (per channel, per marker: auto then user).
            for (int c = 0; c < 3; ++c) {
                for (int k = 0; k < 8; ++k) {
                    emitLoc(ecgAuto[c][k], row);
                    if (k != 3) emitLoc(ecgUser[c][k], row);
                }
            }
            // PPG (6), ABP/ART/ART_PULM (5 each).
            for (int k = 0; k < 7; ++k) { emitLoc(ppgAuto[k], row);  emitLoc(ppgUser[k], row); }
            for (int k = 0; k < 5; ++k) { emitLoc(abpAuto[k], row);  emitLoc(abpUser[k], row); }
            for (int k = 0; k < 5; ++k) { emitLoc(artAuto[k], row);  emitLoc(artUser[k], row); }
            for (int k = 0; k < 5; ++k) { emitLoc(artpAuto[k], row); emitLoc(artpUser[k], row); }
            for (int gc = 0; gc < 3; ++gc) {
                emitLoc(egl[gc].p_wave, row); emitLoc(egl[gc].q_onset, row);
                emitLoc(egl[gc].r_wave, row); emitLoc(egl[gc].t_peak, row);
            }
            emitLoc(pgl.foot, row); emitLoc(pgl.p1, row);
            f << '\n';
        }
    }
    fprintf(stderr, "[tmplcsv] wrote %s\n", path.toStdString().c_str());
    f.close();

    // On the final (q_align) pass, merge the r_align + q_align passes into
    // one wide CSV (columns suffixed _r / _q) and drop the per-pass files.
    if (m_qAlignPass) {
        const QString rPath = outDir.filePath(m_subjectId + "_template_r_align.csv");
        const QString qPath = outDir.filePath(m_subjectId + "_template_q_align.csv");
        const QString mPath = outDir.filePath(m_subjectId + "_template.csv");
        if (mergeAlignedPassCsvs(rPath.toStdString(), qPath.toStdString(),
            mPath.toStdString())) {
            QFile::remove(rPath);
            QFile::remove(qPath);
            fprintf(stderr, "[tmplcsv] merged -> %s\n", mPath.toStdString().c_str());
        }
    }
}

void TemplateViewerWindow::applyMarkerVisibility() {
    for (auto* pw : m_allPlots) {
        pw->setShowEcgMarkers(m_showEcgMarkers);
        pw->setShowPpgMarkers(m_showPpgMarkers);
        pw->setShowAbpMarkers(m_showAbpMarkers);
        pw->setShowArtMarkers(m_showArtMarkers);
        pw->setShowArtPulmMarkers(m_showArtPulmMarkers);
        pw->setShowEcgTrace(m_showEcgTrace);
        pw->setShowPpgTrace(m_showPpgTrace);
        pw->setShowAbpTrace(m_showAbpTrace);
        pw->setShowArtTrace(m_showArtTrace);
        pw->setShowArtPulmTrace(m_showArtPulmTrace);
    }
}

// ========================================================================
// Marker movement
// ========================================================================

void TemplateViewerWindow::refreshBinMarkers(int binIdx) {
    for (int li = 0; li < (int)m_pageGlobalIdx.size(); ++li) {
        if (m_pageGlobalIdx[li] != binIdx) continue;
        const TemplateBin& b = m_bins[binIdx];
        for (auto* pw : m_binPlots[li]) {
            int c = pw->leadIndex();
            pw->setMarker(BinPlotWidget::EcgPPeak, b.p_peak_ch[c]);
            pw->setMarker(BinPlotWidget::EcgQBegin, b.q_begin_ch[c]);
            pw->setMarker(BinPlotWidget::EcgRPeak, b.r_peak_ch[c]);
            pw->setMarker(BinPlotWidget::EcgSEnd, b.s_end_ch[c]);
            pw->setMarker(BinPlotWidget::EcgTBegin, b.t_begin_ch[c]);
            pw->setMarker(BinPlotWidget::EcgTEnd, b.t_end_ch[c]);
            pw->setMarker(BinPlotWidget::PpgOnset, b.ppg_onset);
            pw->setMarker(BinPlotWidget::PpgP50, b.ppg_p50);
            pw->setMarker(BinPlotWidget::PpgPeak, b.ppg_peak);
            pw->setMarker(BinPlotWidget::PpgDicrotic, b.ppg_dicrotic);
            pw->setMarker(BinPlotWidget::PpgPeak2, b.ppg_peak2);
            pw->setMarker(BinPlotWidget::PpgT80, b.ppg_t80);
            pw->setMarker(BinPlotWidget::PpgEnd, b.ppg_end);
            pw->setMarker(BinPlotWidget::AbpOnset, b.abp_onset);
            pw->setMarker(BinPlotWidget::AbpPeak, b.abp_peak);
            pw->setMarker(BinPlotWidget::AbpDicrotic, b.abp_dicrotic);
            pw->setMarker(BinPlotWidget::AbpPeak2, b.abp_peak2);
            pw->setMarker(BinPlotWidget::AbpEnd, b.abp_end);
            pw->setMarker(BinPlotWidget::ArtOnset, b.art_onset);
            pw->setMarker(BinPlotWidget::ArtPeak, b.art_peak);
            pw->setMarker(BinPlotWidget::ArtDicrotic, b.art_dicrotic);
            pw->setMarker(BinPlotWidget::ArtPeak2, b.art_peak2);
            pw->setMarker(BinPlotWidget::ArtEnd, b.art_end);
            pw->setMarker(BinPlotWidget::ArtPulmOnset, b.art_pulm_onset);
            pw->setMarker(BinPlotWidget::ArtPulmPeak, b.art_pulm_peak);
            pw->setMarker(BinPlotWidget::ArtPulmDicrotic, b.art_pulm_dicrotic);
            pw->setMarker(BinPlotWidget::ArtPulmPeak2, b.art_pulm_peak2);
            pw->setMarker(BinPlotWidget::ArtPulmEnd, b.art_pulm_end);
            pw->refreshGlyphs();   // one recapture per widget, after all markers set
        }
        break;
    }
}

void TemplateViewerWindow::onMarkerMoved(int binIdx, int leadIdx,
    int marker, int newIdx)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    TemplateBin& b = m_bins[binIdx];

    if (BinPlotWidget::markerIsEcg(marker)) {
        if (leadIdx < 0 || leadIdx > 2) return;
        auto ecgGet = [&](TemplateBin& tb) -> int {
            switch (marker) {
            case BinPlotWidget::EcgPPeak:  return tb.p_peak_ch[leadIdx];
            case BinPlotWidget::EcgQBegin: return tb.q_begin_ch[leadIdx];
            case BinPlotWidget::EcgRPeak:  return tb.r_peak_ch[leadIdx];
            case BinPlotWidget::EcgSEnd:   return tb.s_end_ch[leadIdx];
            case BinPlotWidget::EcgTBegin:  return tb.t_begin_ch[leadIdx];
            case BinPlotWidget::EcgTEnd:   return tb.t_end_ch[leadIdx];
            }
            return -1;
            };
        auto ecgSet = [&](TemplateBin& tb, int v) {
            switch (marker) {
            case BinPlotWidget::EcgPPeak:  tb.p_peak_ch[leadIdx] = v; break;
            case BinPlotWidget::EcgQBegin: tb.q_begin_ch[leadIdx] = v; break;
            case BinPlotWidget::EcgRPeak:  tb.r_peak_ch[leadIdx] = v; break;
            case BinPlotWidget::EcgSEnd:   tb.s_end_ch[leadIdx] = v; break;
            case BinPlotWidget::EcgTBegin:  tb.t_begin_ch[leadIdx] = v; break;
            case BinPlotWidget::EcgTEnd:   tb.t_end_ch[leadIdx] = v; break;
            }
            };

        const int oldIdx = ecgGet(b);
        ecgSet(b, newIdx);
        const int delta = newIdx - oldIdx;

        if (m_moveSubsequent && oldIdx >= 0) {
            // Shift subsequent bins by the SAME AMOUNT (not to the same time).
            for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
                ChannelTemplateData* chs[3] = {
                    &m_bins[i].ch1, &m_bins[i].ch2, &m_bins[i].ch3
                };
                if (chs[leadIdx]->ecgTemplate_raw.empty()) continue;
                if (m_bins[i].bad_r_ch[leadIdx]) continue;
                const int cur = ecgGet(m_bins[i]);
                if (cur < 0) continue;
                const int n = (int)chs[leadIdx]->ecgTemplate_raw.size();
                ecgSet(m_bins[i], std::clamp(cur + delta, 0, n - 1));
            }
            // Push updated markers to any subsequent bins on screen.
            for (int li = 0; li < (int)m_pageGlobalIdx.size(); ++li) {
                int gi = m_pageGlobalIdx[li];
                if (gi > binIdx) refreshBinMarkers(gi);
            }
        }
        return;
    }

    if (BinPlotWidget::markerIsPpg(marker)) {
        auto ppgGet = [&](TemplateBin& tb) -> int {
            switch (marker) {
            case BinPlotWidget::PpgOnset:    return tb.ppg_onset;
            case BinPlotWidget::PpgP50:      return tb.ppg_p50;
            case BinPlotWidget::PpgPeak:     return tb.ppg_peak;
            case BinPlotWidget::PpgDicrotic: return tb.ppg_dicrotic;
            case BinPlotWidget::PpgPeak2:    return tb.ppg_peak2;
            case BinPlotWidget::PpgT80:      return tb.ppg_t80;
            case BinPlotWidget::PpgEnd:      return tb.ppg_end;
            }
            return -1;
            };
        auto ppgSet = [&](TemplateBin& tb, int v) {
            switch (marker) {
            case BinPlotWidget::PpgOnset:    tb.ppg_onset = v; break;
            case BinPlotWidget::PpgP50:      tb.ppg_p50 = v; break;
            case BinPlotWidget::PpgPeak:     tb.ppg_peak = v; break;
            case BinPlotWidget::PpgDicrotic: tb.ppg_dicrotic = v; break;
            case BinPlotWidget::PpgPeak2:    tb.ppg_peak2 = v; break;
            case BinPlotWidget::PpgT80:      tb.ppg_t80 = v; break;
            case BinPlotWidget::PpgEnd:      tb.ppg_end = v; break;
            }
            };

        const int oldIdx = ppgGet(b);
        ppgSet(b, newIdx);
        refreshBinMarkers(binIdx);
        const int delta = newIdx - oldIdx;

        if (m_moveSubsequent && oldIdx >= 0) {
            for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
                if (m_bins[i].ppg_issue != 0) continue;
                const int cur = ppgGet(m_bins[i]);
                if (cur < 0) continue;
                const int n = (int)m_bins[i].ppgTemplate.size();
                ppgSet(m_bins[i], std::clamp(cur + delta, 0, n - 1));
            }
            for (int li = 0; li < (int)m_pageGlobalIdx.size(); ++li) {
                int gi = m_pageGlobalIdx[li];
                if (gi > binIdx && m_bins[gi].ppg_issue == 0) {
                    refreshBinMarkers(gi);
                }
            }
        }
        return;
    }

    // ---- Arterial markers (ABP / ART / ART_PULM) --------------------------
    // Shared across leads like PPG. Route to the right channel's fields and
    // propagate to subsequent bins when Move-Subsequent is on.
    if (BinPlotWidget::markerIsArterial(marker)) {
        // Select the channel's field pointers, issue flag, and trace by group.
        auto assign = [&](TemplateBin& tb, int mk, int val) {
            switch (mk) {
            case BinPlotWidget::AbpOnset:    tb.abp_onset = val; break;
            case BinPlotWidget::AbpPeak:     tb.abp_peak = val; break;
            case BinPlotWidget::AbpDicrotic: tb.abp_dicrotic = val; break;
            case BinPlotWidget::AbpPeak2:       tb.abp_peak2 = val; break;
            case BinPlotWidget::AbpEnd:      tb.abp_end = val; break;
            case BinPlotWidget::ArtOnset:    tb.art_onset = val; break;
            case BinPlotWidget::ArtPeak:     tb.art_peak = val; break;
            case BinPlotWidget::ArtDicrotic: tb.art_dicrotic = val; break;
            case BinPlotWidget::ArtPeak2:       tb.art_peak2 = val; break;
            case BinPlotWidget::ArtEnd:      tb.art_end = val; break;
            case BinPlotWidget::ArtPulmOnset:    tb.art_pulm_onset = val; break;
            case BinPlotWidget::ArtPulmPeak:     tb.art_pulm_peak = val; break;
            case BinPlotWidget::ArtPulmDicrotic: tb.art_pulm_dicrotic = val; break;
            case BinPlotWidget::ArtPulmPeak2:       tb.art_pulm_peak2 = val; break;
            case BinPlotWidget::ArtPulmEnd:      tb.art_pulm_end = val; break;
            }
            };
        auto artGet = [&](TemplateBin& tb, int mk) -> int {
            switch (mk) {
            case BinPlotWidget::AbpOnset:    return tb.abp_onset;
            case BinPlotWidget::AbpPeak:     return tb.abp_peak;
            case BinPlotWidget::AbpDicrotic: return tb.abp_dicrotic;
            case BinPlotWidget::AbpPeak2:    return tb.abp_peak2;
            case BinPlotWidget::AbpEnd:      return tb.abp_end;
            case BinPlotWidget::ArtOnset:    return tb.art_onset;
            case BinPlotWidget::ArtPeak:     return tb.art_peak;
            case BinPlotWidget::ArtDicrotic: return tb.art_dicrotic;
            case BinPlotWidget::ArtPeak2:    return tb.art_peak2;
            case BinPlotWidget::ArtEnd:      return tb.art_end;
            case BinPlotWidget::ArtPulmOnset:    return tb.art_pulm_onset;
            case BinPlotWidget::ArtPulmPeak:     return tb.art_pulm_peak;
            case BinPlotWidget::ArtPulmDicrotic: return tb.art_pulm_dicrotic;
            case BinPlotWidget::ArtPulmPeak2:    return tb.art_pulm_peak2;
            case BinPlotWidget::ArtPulmEnd:      return tb.art_pulm_end;
            }
            return -1;
            };
        auto channelTrace = [&](TemplateBin& tb, int mk,
            const std::vector<double>*& tr, uint8_t*& iss) {
                if (BinPlotWidget::markerIsAbp(mk)) { tr = &tb.abpTemplate; iss = &tb.abp_issue; }
                else if (BinPlotWidget::markerIsArt(mk)) { tr = &tb.artTemplate; iss = &tb.art_issue; }
                else { tr = &tb.artPulmTemplate; iss = &tb.art_pulm_issue; }
            };

        const int oldIdx = artGet(b, marker);
        assign(b, marker, newIdx);
        refreshBinMarkers(binIdx);
        const int delta = newIdx - oldIdx;

        if (m_moveSubsequent && oldIdx >= 0) {
            for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
                const std::vector<double>* tr = nullptr; uint8_t* iss = nullptr;
                channelTrace(m_bins[i], marker, tr, iss);
                if (!iss || *iss != 0) continue;           // absent/bad channel
                if (!tr) continue;
                const int cur = artGet(m_bins[i], marker);
                if (cur < 0) continue;
                assign(m_bins[i], marker, std::clamp(cur + delta, 0, (int)tr->size() - 1));
            }
            for (int li = 0; li < (int)m_pageGlobalIdx.size(); ++li) {
                int gi = m_pageGlobalIdx[li];
                if (gi > binIdx) refreshBinMarkers(gi);
            }
        }
        return;
    }
}

void TemplateViewerWindow::onMarkerDragStarted(int, int, int) {}

// ========================================================================
// BadR / BadPPG
// ========================================================================

void TemplateViewerWindow::onBadRToggled(int binIdx, int leadIdx, bool bad) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    if (leadIdx < 0 || leadIdx > 2) return;
    m_bins[binIdx].bad_r_ch[leadIdx] = bad;
}

void TemplateViewerWindow::onBadPPGToggled(int binIdx, bool bad) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    m_bins[binIdx].ppg_issue = bad ? 1 : 0;

    if (bad) {
        for (int c = 0; c < 3; ++c)
            m_bins[binIdx].bad_r_ch[c] = false;
    }

    for (int li = 0; li < (int)m_pageGlobalIdx.size(); ++li) {
        if (m_pageGlobalIdx[li] == binIdx) {
            auto s = bad ? BinPlotWidget::State::BadPPG
                : BinPlotWidget::State::Good;
            for (auto* pw : m_binPlots[li])
                pw->setState(s);
            break;
        }
    }
}

void TemplateViewerWindow::save_bin_and_csv() {
    // Do NOT re-seed here. Markers were seeded once at loadSubject and
    // updated by user drags via onMarkerMoved. Re-seeding now would wipe
    // every user edit.

    // .bin goes to markingPath (qtvi folder, config's qtvi_marker_path)
    // -- that's where the viewer LOADS state from on subsequent runs.
    // .csv goes to csv_for_analysis alongside <subject>_templates.csv.
    QDir csvDir(QDir(m_templateDir).absoluteFilePath("../csv_for_analysis"));
    if (!csvDir.exists()) csvDir.mkpath(".");
    QDir binDir(m_markingPath);
    if (!binDir.exists()) binDir.mkpath(".");

    try {
        const QString pass = m_qAlignPass ? "q_align" : "r_align";
        QString outPath = m_markingPath + "/" + m_subjectId + "_template_mark_" + pass + ".bin";
        writeTemplateMarkingsBin(outPath.toStdString(), m_bins);
        std::cout << "Saved: " << outPath.toStdString() << "\n";

        QString csvPath = csvDir.absolutePath() + "/" +
            m_subjectId + "_template_mark_" + pass + ".csv";
        writeTemplateMarkingsCsv(csvPath.toStdString(), m_bins,
            m_subjectId.toStdString(), m_sampleRate);
        std::cout << "Saved: " << csvPath.toStdString() << "\n";
        std::cout.flush();
    }
    catch (const std::exception& e) {
        QMessageBox::critical(this, "Save failed",
            QString("Could not write markings for %1:\n\n%2\n\n"
                "If the CSV is open in Excel, close it and try again.")
            .arg(m_subjectId, e.what()));
        return;   // don't emit finished(); let the user retry
    }

    if (!m_qAlignPass) {
        // First (R-aligned) pass just saved. Switch to the Q-aligned pass:
        // ask the controller to regenerate with Q-alignment and reload, and
        // relabel the button. Do NOT emit finished() yet.
        m_qAlignPass = true;
        ui->finishButton->setText("Finish");
        emit requestQAlignReload();
        return;
    }

    // Final (q_align) pass: merge the r_align + q_align markings CSVs into one
    // wide CSV (columns suffixed _r / _q) and drop the per-pass CSVs. The
    // per-pass .bin files are kept (used to restore state on reload).
    {
        const QString rCsv = csvDir.absolutePath() + "/" + m_subjectId + "_template_mark_r_align.csv";
        const QString qCsv = csvDir.absolutePath() + "/" + m_subjectId + "_template_mark_q_align.csv";
        const QString mCsv = csvDir.absolutePath() + "/" + m_subjectId + "_template_mark.csv";
        if (mergeAlignedPassCsvs(rCsv.toStdString(), qCsv.toStdString(), mCsv.toStdString())) {
            QFile::remove(rCsv);
            QFile::remove(qCsv);
            std::cout << "Merged: " << mCsv.toStdString() << "\n";
        }
    }

    emit finished();
}