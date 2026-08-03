#include "TemplateViewerWindow.hpp"
#include "ui_TemplateViewerWindow.h"
#include "feature_marks.hpp"
#include "NormalizeFeatures.hpp"
#include "peak_finding/FilterUtils.hpp"   // notch_filter for display-time toggle
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
#include <QCheckBox>
#include <QDockWidget>
#include <QVBoxLayout>
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

    //when you are moving a marker, do you move just that marker, or all subsequent markers by the same delta, or all subsequent markers to the same raw index?
    auto* moveGroup = new QButtonGroup(this);
    moveGroup->addButton(ui->move_individual);
    moveGroup->addButton(ui->move_subsequent_delta);
    moveGroup->addButton(ui->move_subsequent_raw);
    moveGroup->setExclusive(true);

    connect(ui->move_individual, &QRadioButton::toggled, this, [this](bool on) {
        if (on) m_moveMode = MoveMode::Individual;
        });
    connect(ui->move_subsequent_delta, &QRadioButton::toggled, this, [this](bool on) {
        if (on) m_moveMode = MoveMode::SubsequentDelta;
        });
    connect(ui->move_subsequent_raw, &QRadioButton::toggled, this, [this](bool on) {
        if (on) m_moveMode = MoveMode::SubsequentRaw;
        });

    // Sync m_moveMode to whichever button Designer has checked by default.
    if (ui->move_individual->isChecked())      m_moveMode = MoveMode::Individual;
    else if (ui->move_subsequent_delta->isChecked()) m_moveMode = MoveMode::SubsequentDelta;
    else if (ui->move_subsequent_raw->isChecked())   m_moveMode = MoveMode::SubsequentRaw;

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

    // Optional display-time notch filter toggle. Wired defensively via
    // findChild so this compiles/runs even if the .ui doesn't (yet) declare
    // a checkbox named "notch_filter"; when the widget is present, ticking
    // it re-runs showPage() with each template pushed through notch_filter
    // before drawing. When absent, this block is silently a no-op.
    if (auto* notchBox = findChild<QCheckBox*>("notch_filter")) {
        connect(notchBox, &QCheckBox::toggled, this, [this](bool on) {
            m_notchFilterOn = on;
            showPage();   // full page redraw; templates re-filtered on the way in
            });
    }

    // B2 focus mode: two stacked panels (QRS and JT views) in a right-side
    // dock. The J-point (S-end) is shared between them, so a J-point
    // selection/edit refreshes both (see refreshFocus). Created in code (not
    // the .ui) so the existing Designer layout is untouched.
    {
        auto* dock = new QDockWidget(QStringLiteral("Focus"), this);
        dock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
        auto* holder = new QWidget(dock);
        auto* vlay = new QVBoxLayout(holder);
        vlay->setContentsMargins(4, 4, 4, 4);
        vlay->setSpacing(6);
        m_focusQrs = new FocusPanelWidget(holder);
        m_focusJt = new FocusPanelWidget(holder);
        vlay->addWidget(m_focusQrs);
        vlay->addWidget(m_focusJt);
        holder->setLayout(vlay);
        dock->setWidget(holder);
        addDockWidget(Qt::RightDockWidgetArea, dock);
    }
}

TemplateViewerWindow::~TemplateViewerWindow() { delete ui; }

void TemplateViewerWindow::setAnchorLabel(const QString& s) {
    m_anchorLabel = s;
    // Anchor shown in the window title bar (the UI's top bar).
    if (!m_subjectId.isEmpty())
        setWindowTitle(QString("Template Marking - %1  [%2]").arg(m_subjectId, s));
    else
        setWindowTitle(QString("Template Marking  [%1]").arg(s));
}

// ========================================================================
// Helpers
// ========================================================================


// Prevents PPG from going over the ECG window when dragged
static int ecgClipLenFor(const TemplateBin& tb) {
    const ChannelTemplateData* chs[3] = { &tb.ch1, &tb.ch2, &tb.ch3 };
    int mn = -1;
    for (const auto* ch : chs) {
        const int l = static_cast<int>(ch->ecgTemplate_raw.size());
        if (l > 0) mn = (mn < 0) ? l : std::min(mn, l);
    }
    return mn;
}


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
        captureCurrentPage();   // snapshot the page we are leaving
        ++m_currentPage;
        showPage();
    }
}

void TemplateViewerWindow::onPrevPage() {
    if (m_currentPage > 0) {
        captureCurrentPage();   // snapshot the page we are leaving
        --m_currentPage;
        showPage();
    }
}

// ========================================================================
// Load subject
// ========================================================================

void TemplateViewerWindow::loadSubject(const QString& templatePath, const QString& markingPath,
    const QString& subjectId, double sampleRateHz,
    double ppgRateHz, double abpRateHz, double artRateHz, double artPulmRateHz,
    int notchFilterHz) {

    m_markingPath = markingPath;
    m_templateDir = QFileInfo(templatePath).absolutePath();
    m_subjectId = subjectId;
    m_sampleRate = sampleRateHz;
    m_ppgRateHz = ppgRateHz;
    m_abpRateHz = abpRateHz;
    m_artRateHz = artRateHz;
    m_artPulmRateHz = artPulmRateHz;
    m_notchFilterHz = notchFilterHz;
    setWindowTitle(QString("Template Marking - %1  [%2]").arg(subjectId, m_anchorLabel));
    ui->subjectLabel->setText(subjectId);
    // Button reads "Finish" only on the final pass of the anchor cycle;
    // every earlier pass reads "Finish and Next".
    const int currentPassIndex = m_anchorStep + 1;
    const bool lastPass = (currentPassIndex + 1 >= m_anchorPassCount);
    ui->finishButton->setText(lastPass ? "Finish" : "Finish and Next");

    try {
        m_bins = readTemplateInfoBin(templatePath.toStdString(), m_currentAnchor);
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
    for (auto& b : m_bins) FeatureMarks::seed_all(b, m_sampleRate, currentAnchor());

    // If this subject was already marked in a previous session, restore
    // those marker positions from the single canonical marking file. If
    // it doesn't exist (never marked before), the fresh auto-seed above
    // stands unchanged.
    const QDir markingDir(m_markingPath);
    const QString canonical = markingDir.filePath(m_subjectId + "_template_markings.bin");
    const QString partial = canonical + ".partial";

    // Source + modality by cycle state:
    //  - mid-cycle (.partial present): restore PULSE only; ECG stays at this
    //    anchor's fresh auto-seed (seed_all above).
    //  - finished subject re-opened (canonical present, no partial): restore
    //    everything.
    //  - neither: fresh auto-seed stands.
    bool markersReloaded = false;
    if (QFile::exists(partial)) {
        markersReloaded = restoreMarkersFrom(partial, /*ecg=*/false, /*pulse=*/true);
    }
    else if (QFile::exists(canonical)) {
        markersReloaded = restoreMarkersFrom(canonical, /*ecg=*/true, /*pulse=*/true);
    }
    fprintf(stderr, "[markers] %s for subject %s (%s pass)\n",
        markersReloaded ? "RELOADED prior markers" : "using FRESH auto-seed (no prior markers applied)",
        m_subjectId.toStdString().c_str(), m_anchorLabel.toStdString().c_str());

    computeGlobalRefs();

    m_capturedPages.clear();   // new subject: nothing saved yet
    showPage();
}

bool TemplateViewerWindow::restoreMarkersFrom(const QString& markingsBinPath,
    bool ecg, bool pulse) {
    try {
        std::vector<TemplateBin> saved = readTemplateMarkingsBin(markingsBinPath.toStdString());
        if (saved.empty()) {
            fprintf(stderr, "[markers] NOT reloaded from %s -- file has 0 bins\n",
                markingsBinPath.toStdString().c_str());
            return false;
        }
        if (saved.size() > m_bins.size()) {
            fprintf(stderr, "[markers] ERROR: %s has %zu bins but current subject has only %zu -- "
                "ignoring the extra %zu bin(s) in the file\n",
                markingsBinPath.toStdString().c_str(), saved.size(), m_bins.size(),
                saved.size() - m_bins.size());
        }
        else if (saved.size() < m_bins.size()) {
            fprintf(stderr, "[markers] ERROR: %s has only %zu bins but current subject has %zu -- "
                "leaving the extra %zu bin(s) at their fresh auto-seed (unedited)\n",
                markingsBinPath.toStdString().c_str(), saved.size(), m_bins.size(),
                m_bins.size() - saved.size());
        }
        const size_t n = std::min(saved.size(), m_bins.size());

        // Clamp every restored marker index against the CURRENT bin's own
        // template length before applying it. Marker indices from the saved
        // file were valid for whatever templates existed when it was
        // written -- if this subject's templates regenerated with different
        // per-bin lengths since then (which the bin-count mismatch above
        // already tells us is possible), blindly applying an old index that
        // now exceeds the current array length causes an out-of-bounds read
        // the first time anything indexes into that array with it (a real
        // crash, not just a cosmetic misplacement). Rejected markers keep
        // whatever seed_all() already put there.
        size_t rejectedCount = 0;
        auto safeIdx = [&](int savedVal, int currentVal, size_t len) -> int {
            if (savedVal >= 0 && static_cast<size_t>(savedVal) < len) return savedVal;
            if (savedVal >= 0) ++rejectedCount;   // only count real (non-sentinel) rejections
            return currentVal;
            };

        for (size_t i = 0; i < n; ++i) {
            TemplateBin& d = m_bins[i];
            const TemplateBin& s = saved[i];
            const size_t ecgLen[3] = {
                d.ch1.ecgTemplate_raw.size(), d.ch2.ecgTemplate_raw.size(), d.ch3.ecgTemplate_raw.size()
            };
            const size_t ppgLen = d.ppgTemplate.size();
            const size_t abpLen = d.abpTemplate.size();
            const size_t artLen = d.artTemplate.size();
            const size_t artPLen = d.artPulmTemplate.size();
            // ECG per-channel user markers (copied raw, bounds-checked). R is
            // NOT copied: it's an auto-only anchor and must come from this
            // pass's own template r_col (seeded above), never from a saved file.
            // Copy every anchor's saved marker set (each independent), bounds-
            // checked per channel. R is NOT copied (auto-only, re-derived).
            if (ecg) {
                for (const auto& kv : s.markers_by_anchor) {
                    const TemplateBin::MarkerSet& sm = kv.second;
                    TemplateBin::MarkerSet& dm = d.marks(kv.first);
                    for (int c = 0; c < 3; ++c) {
                        dm.p_begin_ch[c] = safeIdx(sm.p_begin_ch[c], dm.p_begin_ch[c], ecgLen[c]);
                        dm.p_peak_ch[c] = safeIdx(sm.p_peak_ch[c], dm.p_peak_ch[c], ecgLen[c]);
                        dm.q_begin_ch[c] = safeIdx(sm.q_begin_ch[c], dm.q_begin_ch[c], ecgLen[c]);
                        dm.s_end_ch[c] = safeIdx(sm.s_end_ch[c], dm.s_end_ch[c], ecgLen[c]);
                        dm.t_begin_ch[c] = safeIdx(sm.t_begin_ch[c], dm.t_begin_ch[c], ecgLen[c]);
                        dm.t_end_ch[c] = safeIdx(sm.t_end_ch[c], dm.t_end_ch[c], ecgLen[c]);
                    }
                }
                for (int c = 0; c < 3; ++c) d.bad_r_ch[c] = s.bad_r_ch[c];
            } // if (ecg)

            if (pulse) {
                d.bad_ppg = s.bad_ppg;
                d.ppg_onset = safeIdx(s.ppg_onset, d.ppg_onset, ppgLen);
                d.ppg_t50 = safeIdx(s.ppg_t50, d.ppg_t50, ppgLen);
                d.ppg_t80 = safeIdx(s.ppg_t80, d.ppg_t80, ppgLen);
                d.ppg_peak = safeIdx(s.ppg_peak, d.ppg_peak, ppgLen);
                d.ppg_dicrotic = safeIdx(s.ppg_dicrotic, d.ppg_dicrotic, ppgLen);
                d.ppg_peak2 = safeIdx(s.ppg_peak2, d.ppg_peak2, ppgLen);
                d.ppg_end = safeIdx(s.ppg_end, d.ppg_end, ppgLen);
                d.abp_issue = s.abp_issue;
                d.abp_onset = safeIdx(s.abp_onset, d.abp_onset, abpLen);
                d.abp_peak = safeIdx(s.abp_peak, d.abp_peak, abpLen);
                d.abp_dicrotic = safeIdx(s.abp_dicrotic, d.abp_dicrotic, abpLen);
                d.abp_peak2 = safeIdx(s.abp_peak2, d.abp_peak2, abpLen);
                d.abp_end = safeIdx(s.abp_end, d.abp_end, abpLen);
                d.art_issue = s.art_issue;
                d.art_onset = safeIdx(s.art_onset, d.art_onset, artLen);
                d.art_peak = safeIdx(s.art_peak, d.art_peak, artLen);
                d.art_dicrotic = safeIdx(s.art_dicrotic, d.art_dicrotic, artLen);
                d.art_peak2 = safeIdx(s.art_peak2, d.art_peak2, artLen);
                d.art_end = safeIdx(s.art_end, d.art_end, artLen);
                d.art_pulm_issue = s.art_pulm_issue;
                d.art_pulm_onset = safeIdx(s.art_pulm_onset, d.art_pulm_onset, artPLen);
                d.art_pulm_peak = safeIdx(s.art_pulm_peak, d.art_pulm_peak, artPLen);
                d.art_pulm_dicrotic = safeIdx(s.art_pulm_dicrotic, d.art_pulm_dicrotic, artPLen);
                d.art_pulm_peak2 = safeIdx(s.art_pulm_peak2, d.art_pulm_peak2, artPLen);
                d.art_pulm_end = safeIdx(s.art_pulm_end, d.art_pulm_end, artPLen);
            } // if (pulse)
        }
        if (rejectedCount > 0) {
            fprintf(stderr, "[markers] WARNING: %zu marker(s) from %s were out of range for the "
                "current templates and were rejected (kept at auto-seed) instead of applied\n",
                rejectedCount, markingsBinPath.toStdString().c_str());
        }
        fprintf(stderr, "[markers] reloaded from %s (%zu of %zu bin(s) restored)\n",
            markingsBinPath.toStdString().c_str(), n, m_bins.size());
        return true;
    }
    catch (const std::exception& e) {
        fprintf(stderr, "[markers] NOT reloaded from %s -- %s\n",
            markingsBinPath.toStdString().c_str(), e.what());
        return false;   // No / unreadable file: caller keeps the fresh auto-seed.
    }
}

void TemplateViewerWindow::computeGlobalRefs() {
    /*compute the ecg global reference value for QRS complex height(abs(R) + abs(S))) and pulse global ref for
    PPG / ART / ART_PULM(abs(peak) - abs(foot))*/
    for (int c = 0; c < 3; ++c)
        m_ecgGlobalRef[c] = normalize_features::compute_ecg_global_ref(
            m_bins, c, m_sampleRate);
    for (int c = 0; c < 4; ++c) {
        m_pulseGlobalRef[c] = normalize_features::compute_pulse_global_ref(m_bins, c);
    }
}

std::vector<double> TemplateViewerWindow::normalizeEcgTrace(const std::vector<double>& raw, int ch) const {
    const double ref = (ch >= 0 && ch < 3) ? m_ecgGlobalRef[ch] : std::nan("");
    return normalize_features::normalize_ecg_trace(raw, ref);
}

std::vector<double> TemplateViewerWindow::normalize_ppg_or_similar(const std::vector<double>& raw, int footIdx, int pulseChan) const {
    const double ref = (pulseChan >= 0 && pulseChan < 4) ? m_pulseGlobalRef[pulseChan] : std::nan("");
    return normalize_features::normalize_pulse_trace(raw, footIdx, ref);
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
            pw->setChannelRate(BinPlotWidget::Channel::Ecg, m_sampleRate);
            pw->setChannelRate(BinPlotWidget::Channel::Ppg, m_ppgRateHz);
            pw->setChannelRate(BinPlotWidget::Channel::Abp, m_abpRateHz);
            pw->setChannelRate(BinPlotWidget::Channel::Art, m_artRateHz);
            pw->setChannelRate(BinPlotWidget::Channel::ArtPulm, m_artPulmRateHz);

            static const std::vector<double> empty;
            const auto& ecg = leads[li].ecg ? *leads[li].ecg : empty;
            const auto& ppg = hasPPG ? b.ppgTemplate : empty;

            int c = leads[li].channelIndex;
            const double rPeak = static_cast<double>(b.r_peak_ch[c]);

            const std::vector<double>& ecgIqrRaw = (c == 0) ? b.ch1.ecg_template_raw_iqr
                : (c == 1) ? b.ch2.ecg_template_raw_iqr
                : b.ch3.ecg_template_raw_iqr;
            const double ecgRef = (c >= 0 && c < 3) ? m_ecgGlobalRef[c] : std::nan("");
            // Both ecgIqrRaw (Q3-Q1 of raw amplitude) and b.ppg_template_iqr
            // (Q3-Q1 of each beat's own local perfusion-index ratio, computed
            // at build time -- see CreatePPGTemplates.hpp) are pre-ref-division.
            // The only remaining step, for either channel type, is the same
            // scalar /ref used for the mean trace -- normalization lives only
            // in NormalizeFeatures.hpp, this is just calling it.
            const std::vector<double> ecgIqr = normalize_features::scale_array_by_ref(ecgIqrRaw, ecgRef);
            const std::vector<double> ppgIqr = hasPPG
                ? normalize_features::scale_array_by_ref(b.ppg_template_iqr, m_pulseGlobalRef[0])
                : empty;
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
            //
            // Display-time notch filter toggle: when the `notch_filter`
            // checkbox is on AND a valid notch frequency was passed in via
            // loadSubject, push each channel's template through notch_filter
            // at that channel's own rate BEFORE normalization. Purely a
            // viewing convenience: templates on disk are untouched, and
            // toggling off returns to the raw stored templates on the next
            // showPage() (which the checkbox handler triggers).
            //
            // maybeNotch also REBASES: notch_filter (implemented as
            // x - narrow_bandpass(x)) leaves the waveform shape intact but
            // shifts the DC level by whatever tiny residual mean the
            // bandpass emits. That shift moves foot_y (the value at the foot
            // column), which is what normalize_pulse_trace divides by. For
            // pulse channels stored near-baseline-subtracted, foot_y sits
            // near zero, so even a tiny DC shift is a LARGE relative change
            // to foot_y -- the normalized trace then collapses to the
            // bottom of the y-axis (looks like "PPG went away"). Rebasing
            // to the pre-notch foot value keeps the normalization stable
            // and the visible effect really is just the notch, not a
            // baseline-hunt-induced squash.

            const bool notchActive = m_notchFilterOn && m_notchFilterHz > 0;
            auto maybeNotch = [&](const std::vector<double>& sig, double fs, int footIdx) -> std::vector<double> {
                if (!notchActive || sig.empty() || fs <= 0.0) return sig;
                std::vector<double> out = notch_filter(sig, static_cast<double>(m_notchFilterHz), fs);
                if (footIdx >= 0 && footIdx < (int)sig.size() && footIdx < (int)out.size()) {
                    const double shift = sig[footIdx] - out[footIdx];
                    if (std::isfinite(shift) && shift != 0.0)
                        for (auto& v : out) v += shift;
                }
                return out;
                };
            const std::vector<double> ecgSrc = maybeNotch(ecg, m_sampleRate, -1);   // ECG uses /ref, not a foot; no rebase needed
            const std::vector<double> ppgSrc = hasPPG ? maybeNotch(b.ppgTemplate, m_ppgRateHz, b.ppg_onset) : empty;
            const std::vector<double> abpSrc = !b.abpTemplate.empty() ? maybeNotch(b.abpTemplate, m_abpRateHz, b.abp_onset) : b.abpTemplate;
            const std::vector<double> artSrc = !b.artTemplate.empty() ? maybeNotch(b.artTemplate, m_artRateHz, b.art_onset) : b.artTemplate;
            const std::vector<double> artPSrc = !b.artPulmTemplate.empty() ? maybeNotch(b.artPulmTemplate, m_artPulmRateHz, b.art_pulm_onset) : b.artPulmTemplate;

            const std::vector<double> ecgN = normalizeEcgTrace(ecgSrc, c);
            const std::vector<double> ppgN = hasPPG
                ? normalize_ppg_or_similar(ppgSrc, b.ppg_onset, 0)
                : empty;
            const std::vector<double> abpN = !abpSrc.empty()
                ? normalize_ppg_or_similar(abpSrc, b.abp_onset, 1)
                : abpSrc;
            const std::vector<double> artN = !artSrc.empty()
                ? normalize_ppg_or_similar(artSrc, b.art_onset, 2)
                : artSrc;
            const std::vector<double> artPN = !artPSrc.empty()
                ? normalize_ppg_or_similar(artPSrc, b.art_pulm_onset, 3)
                : artPSrc;

            pw->setData(ppgN, ppgIqr, ecgN, ecgIqr,
                b.marks(currentAnchor()).p_peak_ch[c], b.marks(currentAnchor()).q_begin_ch[c], b.r_peak_ch[c],
                b.marks(currentAnchor()).s_end_ch[c], b.marks(currentAnchor()).t_begin_ch[c], b.marks(currentAnchor()).t_end_ch[c],
                b.ppg_onset, b.ppg_t50, b.ppg_peak,
                b.ppg_dicrotic, b.ppg_peak2, b.ppg_t80, b.ppg_end,
                rPeak,
                static_cast<int>(nEcgBeats),
                static_cast<int>(b.ppg_n_beats),
                b.ppg_onset_auto, b.ppg_peak_auto, b.ppg_peak2_auto,
                b.ppg_peak2_found_auto,
                b.ppg_dicrotic_auto, b.ppg_dicrotic_found_auto,
                b.ppg_end_auto, b.ppg_end_found_auto);

            // Frozen ECG autodetect positions for the own-bar glyphs
            // (P-begin, P-peak, Q-begin, S-end, T-end). These come from the
            // *_auto_ch fields and do NOT move when the user drags a bar.
            pw->setEcgAuto({ b.p_begin_auto_ch[c], b.p_peak_auto_ch[c],
                             b.q_begin_auto_ch[c], b.s_end_auto_ch[c],
                             b.t_end_auto_ch[c] });

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
                b.abpTemplate_iqr, b.artTemplate_iqr, b.artPulmTemplate_iqr);
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

            if (b.bad_ppg == 1)
                pw->setState(BinPlotWidget::State::BadPPG);
            else if (b.bad_r_ch[c])
                pw->setState(BinPlotWidget::State::BadR);

            connect(pw, &BinPlotWidget::markerMoved,
                this, &TemplateViewerWindow::onMarkerMoved);
            connect(pw, &BinPlotWidget::markerDragStarted,
                this, &TemplateViewerWindow::onMarkerDragStarted);
            connect(pw, &BinPlotWidget::landmarkSelected,
                this, &TemplateViewerWindow::onLandmarkSelected);
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

}

void TemplateViewerWindow::captureCurrentPage() {
    // Snapshot the page the user is LEAVING (captures its final edited state).
    // Grabbed synchronously so we snap the page still on screen, not the next
    // one; re-captures on each leave so the latest state wins.
    QDir outDir(m_templateDir);
    const int page = m_currentPage;
    const QPixmap shot = ui->scrollContents->grab();
    const QString pass = "_" + m_anchorLabel;   // per-anchor: _R, _Q_ONSET, ...
    const QString fn = outDir.filePath(
        QString("%1_templates_page%2%3.png")
        .arg(m_subjectId)
        .arg(page + 1, 2, 10, QChar('0'))
        .arg(pass));
    shot.save(fn, "PNG");
    m_capturedPages.insert(page);
}

// Suffix every column in `header` (a single header line) with `suffix`,
// EXCEPT the first three (file_id, bin_num, x_ms), which are row keys and
// must stay un-suffixed so the R and Q sides line up when zipped.
static std::string suffixValueColumns(const std::string& header, const std::string& suffix) {
    std::string out;
    out.reserve(header.size() + 32);
    size_t start = 0;
    int colIdx = 0;
    for (size_t i = 0; i <= header.size(); ++i) {
        if (i == header.size() || header[i] == ',') {
            out.append(header, start, i - start);
            if (colIdx >= 3) out.append(suffix);   // skip file_id/bin_num/x_ms
            if (i < header.size()) out.push_back(',');
            start = i + 1;
            ++colIdx;
        }
    }
    return out;
}

// Zip the R canonical file with the just-produced Q content (both as raw CSV
// text). Header line: R header with un-suffixed keys and _r on values, then
// Q's non-key value columns with _q on them. Row lines: paired 1:1; each
// zipped row uses R's file_id/bin_num/x_ms as the keys, then R's value
// columns, then Q's value columns. Row counts must match -- if they don't we
// abort and leave the R file untouched (worse to write a garbled file than
// none). Returns true on success.
static bool zipCanonicalWithQ(const std::string& canonicalPath,
    const std::string& qContent)
{
    std::ifstream rf(canonicalPath);
    if (!rf) {
        fprintf(stderr, "[tmplcsv] zip: cannot open R canonical %s\n", canonicalPath.c_str());
        return false;
    }
    auto readAll = [](std::istream& in) {
        std::vector<std::string> lines; std::string ln;
        while (std::getline(in, ln)) {
            if (!ln.empty() && ln.back() == '\r') ln.pop_back();
            lines.push_back(ln);
        }
        return lines;
        };
    std::vector<std::string> R = readAll(rf);
    rf.close();
    std::istringstream qs(qContent);
    std::vector<std::string> Q = readAll(qs);
    if (R.empty() || Q.empty()) {
        fprintf(stderr, "[tmplcsv] zip: empty input (R=%zu Q=%zu)\n", R.size(), Q.size());
        return false;
    }
    if (R.size() != Q.size()) {
        fprintf(stderr, "[tmplcsv] zip: row count mismatch (R=%zu Q=%zu) -- keeping R file as-is\n",
            R.size(), Q.size());
        return false;
    }

    // Strip Q's first three columns (the shared keys) from every line so
    // we don't duplicate them in the zipped output.
    auto stripFirstThree = [](const std::string& line) -> std::string {
        int commas = 0;
        for (size_t i = 0; i < line.size(); ++i) {
            if (line[i] == ',') {
                if (++commas == 3) return line.substr(i + 1);
            }
        }
        return {};   // fewer than 3 commas -- nothing to append
        };

    std::ofstream out(canonicalPath, std::ios::trunc);
    if (!out) {
        fprintf(stderr, "[tmplcsv] zip: cannot rewrite %s\n", canonicalPath.c_str());
        return false;
    }
    for (size_t i = 0; i < R.size(); ++i) {
        const std::string qTail = stripFirstThree(Q[i]);
        out << R[i];
        if (!qTail.empty()) out << ',' << qTail;
        out << '\n';
    }
    return out.good();
}

void TemplateViewerWindow::writeAlignedTemplateCsv() {
    if (m_bins.empty()) return;

    QDir outDir(QDir(m_templateDir).absoluteFilePath("../csv_for_analysis"));
    if (!outDir.exists()) outDir.mkpath(".");
    // Single canonical file; both passes contribute to it directly (no
    // intermediate _r_align/_q_align files, no post-hoc merge step).
    const QString path = outDir.filePath(m_subjectId + "_template.csv");

    // Build the pass's CSV content in memory. On the R pass this becomes the
    // canonical file directly (with _r suffixes applied to every value column).
    // On the Q pass we hand this buffer to zipCanonicalWithQ, which reads the
    // R file back and appends the Q columns to each row with _q suffixes.
    std::ostringstream f;
    if (!f) { fprintf(stderr, "[tmplcsv] cannot build in-memory buffer for %s\n", path.toStdString().c_str()); return; }

    const double toMs = (m_sampleRate > 0.0) ? 1000.0 / m_sampleRate : 1.0;

    // ---- Header ------------------------------------------------------------
    // Per channel: raw_mv, Normalized_mv, raw_iqr, normalized_iqr.
    static const char* CHANS[] = {
        "ch1", "ch2", "ch3", "ppg", "abp", "art", "art_pulm"
    };
    f << "file_id,bin_num,x_ms";
    for (const char* n : CHANS) {
        f << ',' << n << "_raw_mv"
            << ',' << n << "_Normalized"
            << ',' << n << "_raw_iqr"
            << ',' << n << "_normalized_iqr";
    }

    // Marker-location columns. For each marker: two columns (autodetect,
    // user). A row has "1" in the column iff its row index equals that
    // marker's sample index for the current bin; blank otherwise.
    //
    // ECG markers per channel (8 each x 3 channels x 2 auto/user = 48):
    //   p_peak, q_begin, q_peak (computed), r_peak, s_peak (computed),
    //   s_end, t_peak, t_end.
    static const char* ECG_MARKERS[] = {
        "p_begin", "p_peak", "q_begin", "q_peak", "r_peak", "s_peak",
        "s_end",  "t_begin",  "t_end"
    };
    for (int c = 1; c <= 3; ++c) {
        for (int k = 0; k < 9; ++k) {
            const char* mname = ECG_MARKERS[k];
            const bool userToo = (k != 4);   // r_peak = auto-only
            f << ',' << mname << "_ch" << c << "_location_autodetect";
            if (userToo) f << ',' << mname << "_ch" << c << "_location_user";
        }
    }
    // Pulse markers (PPG has 6 with p50, arterial has 5).
    static const char* PPG_MARKERS[] = {
        "ppg_onset", "ppg_t50", "ppg_peak",
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
    for (const char* g : { "ppg_foot", "ppg_systolic_peak", "ppg_dicrotic", "ppg_end" })
        emitAutoFeatLocHeader(g);
    f << '\n';
    f << std::setprecision(10);

    // Both ecg_template_raw_iqr and every *_template_iqr / *Template_iqr
    // field are pre-ref-division at build time (see CreatePPGTemplates.hpp /
    // build_pulse_template_pair_windowed for the pulse channels), so the
    // only remaining step for any channel is the same scalar /ref used for
    // its mean trace -- normalize_features::scale_array_by_ref is the one
    // place that division happens.

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

        const std::vector<double>& ch1Iqr = b.ch1.ecg_template_raw_iqr;
        const std::vector<double>& ch2Iqr = b.ch2.ecg_template_raw_iqr;
        const std::vector<double>& ch3Iqr = b.ch3.ecg_template_raw_iqr;
        const std::vector<double>& ppgIqrRaw = b.ppg_template_iqr;
        const std::vector<double>& abpIqrRaw = b.abpTemplate_iqr;
        const std::vector<double>& artIqrRaw = b.artTemplate_iqr;
        const std::vector<double>& artPIqrRaw = b.artPulmTemplate_iqr;

        // Normalized mean traces (screen-matching).
        std::vector<double> ch1N = normalizeEcgTrace(ch1R, 0);
        std::vector<double> ch2N = normalizeEcgTrace(ch2R, 1);
        std::vector<double> ch3N = normalizeEcgTrace(ch3R, 2);
        std::vector<double> ppgN = normalize_ppg_or_similar(ppgR, b.ppg_onset, 0);
        std::vector<double> abpN = normalize_ppg_or_similar(abpR, b.abp_onset, 1);
        std::vector<double> artN = normalize_ppg_or_similar(artR, b.art_onset, 2);
        std::vector<double> artPN = normalize_ppg_or_similar(artPR, b.art_pulm_onset, 3);

        // Normalized std traces -- same scalar /ref step for every channel.
        std::vector<double> ch1IqrN = normalize_features::scale_array_by_ref(ch1Iqr, m_ecgGlobalRef[0]);
        std::vector<double> ch2IqrN = normalize_features::scale_array_by_ref(ch2Iqr, m_ecgGlobalRef[1]);
        std::vector<double> ch3IqrN = normalize_features::scale_array_by_ref(ch3Iqr, m_ecgGlobalRef[2]);
        std::vector<double> ppgIqrN = normalize_features::scale_array_by_ref(ppgIqrRaw, m_pulseGlobalRef[0]);
        std::vector<double> abpIqrN = normalize_features::scale_array_by_ref(abpIqrRaw, m_pulseGlobalRef[1]);
        std::vector<double> artIqrN = normalize_features::scale_array_by_ref(artIqrRaw, m_pulseGlobalRef[2]);
        std::vector<double> artPIqrN = normalize_features::scale_array_by_ref(artPIqrRaw, m_pulseGlobalRef[3]);

        // Computed Q/S peaks (both variants) per ECG channel.
        const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
        EcgFeatures ftAuto[3], ftUser[3];
        for (int c = 0; c < 3; ++c) {
            const auto& ecg = chs[c]->ecgTemplate_raw;
            ftAuto[c] = computeEcgFeatures(ecg,
                b.p_peak_auto_ch[c], b.q_begin_auto_ch[c], b.r_peak_auto_ch[c],
                b.s_end_auto_ch[c], b.t_end_auto_ch[c],
                m_sampleRate);
            const TemplateBin::MarkerSet& fmk = b.marks(currentAnchor());
            ftUser[c] = computeEcgFeatures(ecg,
                fmk.p_peak_ch[c], fmk.q_begin_ch[c], b.r_peak_ch[c],
                fmk.s_end_ch[c], fmk.t_end_ch[c],
                m_sampleRate);
        }

        // ECG marker positions per channel, aligned with ECG_MARKERS order:
        //   p_peak, q_begin, q_peak(computed), r_peak, s_peak(computed),
        //   s_end,  t_peak, t_end.
        int ecgAuto[3][9], ecgUser[3][9];
        const TemplateBin::MarkerSet& umk = b.marks(currentAnchor());
        for (int c = 0; c < 3; ++c) {
            ecgAuto[c][0] = b.p_begin_auto_ch[c];
            ecgAuto[c][1] = b.p_peak_auto_ch[c];
            ecgAuto[c][2] = b.q_begin_auto_ch[c];
            ecgAuto[c][3] = ftAuto[c].q_idx;
            ecgAuto[c][4] = b.r_peak_auto_ch[c];
            ecgAuto[c][5] = ftAuto[c].s_idx;
            ecgAuto[c][6] = b.s_end_auto_ch[c];
            ecgAuto[c][7] = b.t_begin_auto_ch[c];
            ecgAuto[c][8] = b.t_end_auto_ch[c];
            ecgUser[c][0] = umk.p_begin_ch[c];
            ecgUser[c][1] = umk.p_peak_ch[c];
            ecgUser[c][2] = umk.q_begin_ch[c];
            ecgUser[c][3] = ftUser[c].q_idx;
            ecgUser[c][4] = b.r_peak_ch[c];
            ecgUser[c][5] = ftUser[c].s_idx;
            ecgUser[c][6] = umk.s_end_ch[c];
            ecgUser[c][7] = umk.t_begin_ch[c];
            ecgUser[c][8] = umk.t_end_ch[c];
        }

        // Pulse marker positions, order matching PPG_MARKERS / ABP_MARKERS etc.
        const int ppgAuto[7] = { b.ppg_onset_auto, b.ppg_t50_auto, b.ppg_peak_auto,
                                 b.ppg_dicrotic_auto, b.ppg_peak2_auto, b.ppg_t80_auto, b.ppg_end_auto };
        const int ppgUser[7] = { b.ppg_onset, b.ppg_t50, b.ppg_peak,
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
            const std::vector<double>* raw_iqr;
            const std::vector<double>* norm_iqr;
        };
        const Col cols[] = {
            { &ch1R,  &ch1N,  &ch1Iqr,  &ch1IqrN  },
            { &ch2R,  &ch2N,  &ch2Iqr,  &ch2IqrN  },
            { &ch3R,  &ch3N,  &ch3Iqr,  &ch3IqrN  },
            { &ppgR,  &ppgN,  &ppgIqrRaw,  &ppgIqrN  },
            { &abpR,  &abpN,  &abpIqrRaw,  &abpIqrN  },
            { &artR,  &artN,  &artIqrRaw,  &artIqrN  },
            { &artPR, &artPN, &artPIqrRaw, &artPIqrN },
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
        // PPG glyph locations are just the bin's own auto fields now --
        // there's no separate glyph recompute anymore (see
        // FeatureMarks::detect_ppg_fiducials).
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
                emitVal(*c.raw_iqr, row);
                emitVal(*c.norm_iqr, row);
            }
            // ECG location columns (per channel, per marker: auto then user).
            for (int c = 0; c < 3; ++c) {
                for (int k = 0; k < 9; ++k) {
                    emitLoc(ecgAuto[c][k], row);
                    if (k != 4) emitLoc(ecgUser[c][k], row);
                }
            }
            // PPG (6), ABP/ART/ART_PULM (5 each).
            for (int k = 0; k < 7; ++k) { emitLoc(ppgAuto[k], row);  emitLoc(ppgUser[k], row); }
            for (int k = 0; k < 5; ++k) { emitLoc(abpAuto[k], row);  emitLoc(abpUser[k], row); }
            for (int k = 0; k < 5; ++k) { emitLoc(artAuto[k], row);  emitLoc(artUser[k], row); }
            for (int k = 0; k < 5; ++k) { emitLoc(artpAuto[k], row); emitLoc(artpUser[k], row); }
            for (int gc = 0; gc < 3; ++gc) {
                // R is never autodetected: emit the passed-in R (r_peak_auto =
                // r_col), not egl[gc].r_peak_glyph (compute_r_wave argmax).
                emitLoc(egl[gc].p_peak_glyph, row); emitLoc(egl[gc].q_begin_glyph, row);
                emitLoc(b.r_peak_auto_ch[gc], row); emitLoc(egl[gc].t_peak_glyph, row);
            }
            emitLoc(b.ppg_onset_auto, row); emitLoc(b.ppg_peak_auto, row);
            emitLoc(b.ppg_dicrotic_auto, row);
            emitLoc(b.ppg_end_auto, row);
            f << '\n';
        }
    }
    // Serialize the in-memory content, suffix its value columns with the
    // pass tag, and either write it as the fresh canonical file (R pass) or
    // read the existing canonical file back and append the Q columns to
    // each line (Q pass). Row-key columns (file_id/bin_num/x_ms) are never
    // suffixed so the two sides line up cleanly on zip.
    std::string content = f.str();
    if (content.empty()) {
        fprintf(stderr, "[tmplcsv] built empty content for %s -- skipping write\n", path.toStdString().c_str());
        return;
    }
    const size_t nl = content.find('\n');
    if (nl == std::string::npos) {
        fprintf(stderr, "[tmplcsv] malformed content (no newline) for %s\n", path.toStdString().c_str());
        return;
    }
    const std::string header = content.substr(0, nl);
    const std::string body = content.substr(nl);   // includes leading '\n'

    const std::string suffix = "_" + m_anchorLabel.toStdString();  // _R, _Q_ONSET, _J_POINT, ...
    if (m_anchorStep < 0) {
        // R pass (first): overwrite the canonical file with _R-suffixed content.
        std::ofstream out(path.toStdString(), std::ios::trunc);
        if (!out) {
            fprintf(stderr, "[tmplcsv] cannot open %s for R-pass write\n", path.toStdString().c_str());
            return;
        }
        out << suffixValueColumns(header, suffix) << body;
        out.close();
        fprintf(stderr, "[tmplcsv] wrote %s (%s pass)\n",
            path.toStdString().c_str(), m_anchorLabel.toStdString().c_str());
    }
    else {
        // Any later anchor pass: zip this anchor's columns (suffix = anchor
        // name) into the accumulating canonical file. zipCanonicalWithQ widens
        // the file by one group per call, so all anchors accumulate.
        std::string aFull = suffixValueColumns(header, suffix);
        aFull.append(body);
        if (zipCanonicalWithQ(path.toStdString(), aFull))
            fprintf(stderr, "[tmplcsv] zipped %s columns into %s\n",
                suffix.c_str(), path.toStdString().c_str());
        else
            fprintf(stderr, "[tmplcsv] zip failed for %s -- file kept as-is\n", path.toStdString().c_str());
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
            pw->setMarker(BinPlotWidget::EcgPBegin, b.marks(currentAnchor()).p_begin_ch[c]);
            pw->setMarker(BinPlotWidget::EcgPPeak, b.marks(currentAnchor()).p_peak_ch[c]);
            pw->setMarker(BinPlotWidget::EcgQBegin, b.marks(currentAnchor()).q_begin_ch[c]);
            pw->setMarker(BinPlotWidget::EcgRPeak, b.r_peak_ch[c]);
            pw->setMarker(BinPlotWidget::EcgSEnd, b.marks(currentAnchor()).s_end_ch[c]);
            pw->setMarker(BinPlotWidget::EcgTBegin, b.marks(currentAnchor()).t_begin_ch[c]);
            pw->setMarker(BinPlotWidget::EcgTEnd, b.marks(currentAnchor()).t_end_ch[c]);
            pw->setMarker(BinPlotWidget::PpgOnset, b.ppg_onset);
            pw->setMarker(BinPlotWidget::PpgPeak, b.ppg_peak);
            pw->setMarker(BinPlotWidget::PpgDicrotic, b.ppg_dicrotic);
            pw->setMarker(BinPlotWidget::PpgPeak2, b.ppg_peak2);
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
            case BinPlotWidget::EcgPBegin: return tb.marks(currentAnchor()).p_begin_ch[leadIdx];
            case BinPlotWidget::EcgPPeak:  return tb.marks(currentAnchor()).p_peak_ch[leadIdx];
            case BinPlotWidget::EcgQBegin: return tb.marks(currentAnchor()).q_begin_ch[leadIdx];
            case BinPlotWidget::EcgRPeak:  return tb.r_peak_ch[leadIdx];
            case BinPlotWidget::EcgSEnd:   return tb.marks(currentAnchor()).s_end_ch[leadIdx];
            case BinPlotWidget::EcgTBegin:  return tb.marks(currentAnchor()).t_begin_ch[leadIdx];
            case BinPlotWidget::EcgTEnd:   return tb.marks(currentAnchor()).t_end_ch[leadIdx];
            }
            return -1;
            };
        auto ecgSet = [&](TemplateBin& tb, int v) {
            switch (marker) {
            case BinPlotWidget::EcgPBegin: tb.marks(currentAnchor()).p_begin_ch[leadIdx] = v; break;
            case BinPlotWidget::EcgPPeak:  tb.marks(currentAnchor()).p_peak_ch[leadIdx] = v; break;
            case BinPlotWidget::EcgQBegin: tb.marks(currentAnchor()).q_begin_ch[leadIdx] = v; break;
            case BinPlotWidget::EcgRPeak:  tb.r_peak_ch[leadIdx] = v; break;
            case BinPlotWidget::EcgSEnd:   tb.marks(currentAnchor()).s_end_ch[leadIdx] = v; break;
            case BinPlotWidget::EcgTBegin:  tb.marks(currentAnchor()).t_begin_ch[leadIdx] = v; break;
            case BinPlotWidget::EcgTEnd:   tb.marks(currentAnchor()).t_end_ch[leadIdx] = v; break;
            }
            };

        const int oldIdx = ecgGet(b);
        ecgSet(b, newIdx);
        const int delta = newIdx - oldIdx;

        if (m_moveMode != MoveMode::Individual && oldIdx >= 0) {
            ChannelTemplateData* bchs[3] = { &b.ch1, &b.ch2, &b.ch3 };
            const int nDragged = (int)bchs[leadIdx]->ecgTemplate_raw.size();
            const double pct = (nDragged > 1) ? double(newIdx) / (nDragged - 1) : 0.0;

            for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
                ChannelTemplateData* chs[3] = {
                    &m_bins[i].ch1, &m_bins[i].ch2, &m_bins[i].ch3
                };
                if (chs[leadIdx]->ecgTemplate_raw.empty()) continue;
                if (m_bins[i].bad_r_ch[leadIdx]) continue;
                const int cur = ecgGet(m_bins[i]);
                if (cur < 0) continue;
                const int n = (int)chs[leadIdx]->ecgTemplate_raw.size();

                const int target = (m_moveMode == MoveMode::SubsequentDelta)
                    ? cur + delta
                    : (n > 1 ? (int)std::lround(pct * (n - 1)) : 0);
                ecgSet(m_bins[i], std::clamp(target, 0, n - 1));
            }
            for (int li = 0; li < (int)m_pageGlobalIdx.size(); ++li) {
                int gi = m_pageGlobalIdx[li];
                if (gi > binIdx) refreshBinMarkers(gi);
            }
        }
        // B2: an edit to this landmark updates its focus view. For the
        // J-point (S-end) this refreshes BOTH the QRS and JT panels; other
        // landmarks refresh their single panel. refreshFocus() handles the
        // routing.
        refreshFocus(binIdx, leadIdx, marker, newIdx);
        return;
    }

    if (BinPlotWidget::markerIsPpg(marker)) {
        auto ppgGet = [&](TemplateBin& tb) -> int {
            switch (marker) {
            case BinPlotWidget::PpgOnset:    return tb.ppg_onset;
            case BinPlotWidget::PpgPeak:     return tb.ppg_peak;
            case BinPlotWidget::PpgDicrotic: return tb.ppg_dicrotic;
            case BinPlotWidget::PpgPeak2:    return tb.ppg_peak2;
            case BinPlotWidget::PpgEnd:      return tb.ppg_end;
            }
            return -1;
            };
        auto ppgSet = [&](TemplateBin& tb, int v) {
            switch (marker) {
            case BinPlotWidget::PpgOnset:    tb.ppg_onset = v; break;
            case BinPlotWidget::PpgPeak:     tb.ppg_peak = v; break;
            case BinPlotWidget::PpgDicrotic: tb.ppg_dicrotic = v; break;
            case BinPlotWidget::PpgPeak2:    tb.ppg_peak2 = v; break;
            case BinPlotWidget::PpgEnd:      tb.ppg_end = v; break;
            }
            };

        const int oldIdx = ppgGet(b);
        ppgSet(b, newIdx);
        refreshBinMarkers(binIdx);
        const int delta = newIdx - oldIdx;

        if (m_moveMode != MoveMode::Individual && oldIdx >= 0) {
            const int nDragged = (int)b.ppgTemplate.size();
            const double pct = (nDragged > 1) ? double(newIdx) / (nDragged - 1) : 0.0;

            for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
                if (m_bins[i].bad_ppg != 0) continue;
                const int cur = ppgGet(m_bins[i]);
                if (cur < 0) continue;

                const int rawLen = (int)m_bins[i].ppgTemplate.size();
                const int ecgClip = ecgClipLenFor(m_bins[i]);
                const int n = (ecgClip > 0) ? std::min(rawLen, ecgClip) : rawLen;
                if (n <= 0) continue;

                const int target = (m_moveMode == MoveMode::SubsequentDelta)
                    ? cur + delta
                    : (n > 1 ? (int)std::lround(pct * (n - 1)) : 0);
                ppgSet(m_bins[i], std::clamp(target, 0, n - 1));
            }
            for (int li = 0; li < (int)m_pageGlobalIdx.size(); ++li) {
                int gi = m_pageGlobalIdx[li];
                if (gi > binIdx && m_bins[gi].bad_ppg == 0) refreshBinMarkers(gi);
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

        if (m_moveMode != MoveMode::Individual && oldIdx >= 0) {
            const std::vector<double>* trDragged = nullptr; uint8_t* issDragged = nullptr;
            channelTrace(b, marker, trDragged, issDragged);
            const int nDragged = trDragged ? (int)trDragged->size() : 0;
            const double pct = (nDragged > 1) ? double(newIdx) / (nDragged - 1) : 0.0;

            for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
                const std::vector<double>* tr = nullptr; uint8_t* iss = nullptr;
                channelTrace(m_bins[i], marker, tr, iss);
                if (!iss || *iss != 0) continue;
                if (!tr) continue;
                const int cur = artGet(m_bins[i], marker);
                if (cur < 0) continue;

                const int rawLen = (int)tr->size();
                const int ecgClip = ecgClipLenFor(m_bins[i]);
                const int n = (ecgClip > 0) ? std::min(rawLen, ecgClip) : rawLen;
                if (n <= 0) continue;

                const int target = (m_moveMode == MoveMode::SubsequentDelta)
                    ? cur + delta
                    : (n > 1 ? (int)std::lround(pct * (n - 1)) : 0);
                assign(m_bins[i], marker, std::clamp(target, 0, n - 1));
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

// B2 focus mode --------------------------------------------------------------

void TemplateViewerWindow::onLandmarkSelected(int binIdx, int leadIdx,
    int marker, int col)
{
    refreshFocus(binIdx, leadIdx, marker, col);
}

// Rebuild the focus panel(s) for one landmark from the current bin/lead's
// anchored-average stats. Reads mean/sd/n straight from the template the
// viewer already holds:
//   mean = ecgTemplate_raw
//   sd   = ecg_template_raw_iqr  (holds STD, ddof=1 -- despite the _iqr name)
//   n    = ch{1,2,3}_n_beats_raw (per-bin, not per-channel-struct)
// The J-point (S-end) is shared by the QRS and JT views, so selecting/editing
// it refreshes BOTH panels; every other landmark refreshes its own single
// panel.
void TemplateViewerWindow::refreshFocus(int binIdx, int leadIdx,
    int marker, int col)
{
    if (!m_focusQrs && !m_focusJt) return;   // panels not created (nothing to do)
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    TemplateBin& b = m_bins[binIdx];

    // ---- PPG and ARTERIAL landmarks (B2 focus extended to all channels) --
    // These ride their own template (ppgTemplate / abpTemplate / artTemplate
    // / artPulmTemplate) with the matching per-sample std (*_iqr, ddof=1) and
    // the shared pulse beat count (ppg_n_beats -- all pulse channels derive
    // from the same foot-anchored beat set). Routed to the QRS panel as the
    // single focus view for pulse channels (they have no QRS/JT split).
    if (!BinPlotWidget::markerIsEcg(marker)) {
        const std::vector<double>* meanRaw = nullptr;
        const std::vector<double>* iqrRaw = nullptr;
        int pulseChan = -1;   // index into m_pulseGlobalRef: PPG=0,ABP=1,ART=2,ART_PULM=3
        int footIdx = -1;     // this channel's foot/onset column (perfusion-index baseline)
        QString chLabel;
        if (BinPlotWidget::markerIsPpg(marker)) {
            meanRaw = &b.ppgTemplate;      iqrRaw = &b.ppg_template_iqr;      pulseChan = 0; footIdx = b.ppg_onset;      chLabel = "PPG";
        }
        else if (BinPlotWidget::markerIsAbp(marker)) {
            meanRaw = &b.abpTemplate;      iqrRaw = &b.abpTemplate_iqr;      pulseChan = 1; footIdx = b.abp_onset;      chLabel = "ABP";
        }
        else if (BinPlotWidget::markerIsArt(marker)) {
            meanRaw = &b.artTemplate;      iqrRaw = &b.artTemplate_iqr;      pulseChan = 2; footIdx = b.art_onset;      chLabel = "ART";
        }
        else if (BinPlotWidget::markerIsArtPulm(marker)) {
            meanRaw = &b.artPulmTemplate;  iqrRaw = &b.artPulmTemplate_iqr;  pulseChan = 3; footIdx = b.art_pulm_onset; chLabel = "ART_PULM";
        }
        if (!meanRaw || meanRaw->empty()) return;

        // Pulse channels are NOT normalized by a plain scalar (that was the
        // bug -- it left the trace flat). The displayed trace uses a per-
        // sample PERFUSION-INDEX transform relative to the pulse's own foot,
        // then /ref (normalize_ppg_or_similar -> normalize_pulse_trace, see
        // the main plot ~line 508). The mean MUST use that same transform.
        const std::vector<double> mean = normalize_ppg_or_similar(*meanRaw, footIdx, pulseChan);
        // The *_iqr is ALREADY in perfusion-index space (local_ratio_iqr at
        // build time), so it only needs the scalar /ref -- NOT the perfusion
        // transform again (main plot ~line 492). It's a true IQR (Q3-Q1), so
        // convert to an SD estimate (IQR/1.349) for the 95% CI.
        const double ref = (pulseChan >= 0 && pulseChan < 4) ? m_pulseGlobalRef[pulseChan] : std::nan("");
        std::vector<double> sd = normalize_features::scale_array_by_ref(*iqrRaw, ref);
        for (double& s : sd) if (!std::isnan(s)) s /= 1.349;

        const int nBeats = static_cast<int>(b.ppg_n_beats);

        auto pulseLabel = [](int m) -> QString {
            switch (m) {
            case BinPlotWidget::PpgOnset:    return QStringLiteral("Foot");
            case BinPlotWidget::PpgT50:      return QStringLiteral("T50");
            case BinPlotWidget::PpgPeak:     return QStringLiteral("Diastolic Peak");
            case BinPlotWidget::PpgDicrotic: return QStringLiteral("Dicrotic Notch");
            case BinPlotWidget::PpgPeak2:    return QStringLiteral("Systolic Peak");
            case BinPlotWidget::PpgT80:      return QStringLiteral("T80");
            case BinPlotWidget::PpgEnd:      return QStringLiteral("End");
            case BinPlotWidget::AbpOnset: case BinPlotWidget::ArtOnset: case BinPlotWidget::ArtPulmOnset:       return QStringLiteral("onset");
            case BinPlotWidget::AbpPeak: case BinPlotWidget::ArtPeak: case BinPlotWidget::ArtPulmPeak:          return QStringLiteral("peak");
            case BinPlotWidget::AbpDicrotic: case BinPlotWidget::ArtDicrotic: case BinPlotWidget::ArtPulmDicrotic: return QStringLiteral("dicrotic");
            case BinPlotWidget::AbpPeak2: case BinPlotWidget::ArtPeak2: case BinPlotWidget::ArtPulmPeak2:        return QStringLiteral("peak2");
            case BinPlotWidget::AbpEnd: case BinPlotWidget::ArtEnd: case BinPlotWidget::ArtPulmEnd:             return QStringLiteral("end");
            }
            return QStringLiteral("landmark");
            };

        if (m_focusQrs)
            m_focusQrs->setFocus(mean, sd, nBeats, col,
                chLabel + " " + pulseLabel(marker));
        // A pulse landmark has nothing to do with the JT (T-wave) view; clear
        // it so a previously-shown J-point view doesn't stay stuck there.
        if (m_focusJt) m_focusJt->clearFocus();
        return;
    }

    // ---- ECG landmarks ---------------------------------------------------
    if (leadIdx < 0 || leadIdx > 2) return;

    ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
    const ChannelTemplateData& ch = *chs[leadIdx];
    // Templates + iqr are stored PRE-reference-division; scale by the same
    // per-lead ref the displayed ECG trace uses (main plot ~line 484). The
    // ECG *_iqr field already holds a STD (ddof=1) -- CreateEcgTemplates step
    // 7 changed it from IQR to std despite the "iqr" name -- so it feeds the
    // CI directly, NO IQR->SD conversion (unlike the pulse channels, whose
    // *_iqr is a true interquartile range).
    const double eref = m_ecgGlobalRef[leadIdx];
    const std::vector<double> mean = normalize_features::scale_array_by_ref(ch.ecgTemplate_raw, eref);
    const std::vector<double> sd = normalize_features::scale_array_by_ref(ch.ecg_template_raw_iqr, eref);
    const uint64_t nb[3] = { b.ch1_n_beats_raw, b.ch2_n_beats_raw, b.ch3_n_beats_raw };
    const int nBeats = static_cast<int>(nb[leadIdx]);
    if (mean.empty()) return;

    m_lastFocusBin = binIdx;
    m_lastFocusLead = leadIdx;

    auto labelFor = [](int m, const char* side) -> QString {
        switch (m) {
        case BinPlotWidget::EcgPBegin: return QStringLiteral("P onset");
        case BinPlotWidget::EcgPPeak:  return QStringLiteral("P peak");
        case BinPlotWidget::EcgQBegin: return QStringLiteral("Q onset");
        case BinPlotWidget::EcgRPeak:  return QStringLiteral("R peak");
        case BinPlotWidget::EcgSEnd:   return QStringLiteral("J-point (%1)").arg(side);
        case BinPlotWidget::EcgTBegin: return QStringLiteral("T onset");
        case BinPlotWidget::EcgTEnd:   return QStringLiteral("T end");
        }
        return QStringLiteral("landmark");
        };

    if (marker == BinPlotWidget::EcgSEnd) {
        // J-point: shared landmark, but framed differently per panel so the
        // two views aren't identical -- QRS shows it at the RIGHT edge (it
        // ENDS the QRS), JT shows it at the LEFT edge (it STARTS the JT).
        if (m_focusQrs) m_focusQrs->setFocus(mean, sd, nBeats, col, labelFor(marker, "QRS"), 30, -1);
        if (m_focusJt)  m_focusJt->setFocus(mean, sd, nBeats, col, labelFor(marker, "JT"), 30, +1);
        return;
    }

    // Non-J-point ECG landmarks: route to the panel for their segment, and
    // CLEAR the other panel so a stale J-point view (which had set both)
    // doesn't stay stuck when the next selection isn't a J-point.
    const bool isJtSide = (marker == BinPlotWidget::EcgTBegin
        || marker == BinPlotWidget::EcgTEnd);
    if (isJtSide) {
        if (m_focusJt) m_focusJt->setFocus(mean, sd, nBeats, col, labelFor(marker, "JT"));
        if (m_focusQrs) m_focusQrs->clearFocus();
    }
    else {
        if (m_focusQrs) m_focusQrs->setFocus(mean, sd, nBeats, col, labelFor(marker, "QRS"));
        if (m_focusJt) m_focusJt->clearFocus();
    }
}

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
    m_bins[binIdx].bad_ppg = bad ? 1 : 0;

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
    captureCurrentPage();   // snapshot the page being left on Finish

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

    const QString canonicalBin = QDir(m_markingPath).filePath(m_subjectId + "_template_markings.bin");
    const QString partialBin = canonicalBin + ".partial";

    // currentPassIndex: 0 for R, anchor step + 1 for the anchor passes.
    const int currentPassIndex = m_anchorStep + 1;
    const bool finalPass = (currentPassIndex + 1 >= m_anchorPassCount);

    try {
        // Markings .bin: written to a PARTIAL every pass so pulse markers (and
        // in-progress state) carry across anchor reloads. Promoted to the
        // canonical name on the final pass only.
        writeTemplateMarkingsBin(partialBin.toStdString(), m_bins);
        std::cout << "Saved: " << partialBin.toStdString() << "\n";

        // ECG markings CSV: per-anchor suffixed columns, zipped into the
        // canonical file (R pass writes fresh, later passes append). Pulse
        // columns are NOT written here -- they go once, on the final pass.
        const QString csvPath = csvDir.absolutePath() + "/"
            + m_subjectId + "_template_markings.csv";
        const QString tmpPath = csvPath + ".tmp";
        writeTemplateMarkingsCsv(tmpPath.toStdString(), m_bins,
            m_subjectId.toStdString(), m_sampleRate, currentAnchor(),
            MarkingsCsvSection::EcgOnly);
        std::ifstream tin(tmpPath.toStdString());
        std::stringstream tbuf; tbuf << tin.rdbuf();
        tin.close();
        QFile::remove(tmpPath);
        std::string tcontent = tbuf.str();
        const size_t tnl = tcontent.find('\n');
        if (tnl == std::string::npos) {
            throw std::runtime_error("markings CSV writer produced malformed content (no newline)");
        }
        const std::string tHeader = tcontent.substr(0, tnl);
        const std::string tBody = tcontent.substr(tnl);
        const std::string tsuffix = "_" + m_anchorLabel.toStdString();
        if (m_anchorStep < 0) {
            std::ofstream out(csvPath.toStdString(), std::ios::trunc);
            if (!out) throw std::runtime_error("cannot open for write: " + csvPath.toStdString());
            out << suffixValueColumns(tHeader, tsuffix) << tBody;
            out.close();
        }
        else {
            std::string aFull = suffixValueColumns(tHeader, tsuffix);
            aFull.append(tBody);
            if (!zipCanonicalWithQ(csvPath.toStdString(), aFull))
                throw std::runtime_error("could not zip " + tsuffix + " markings CSV into " + csvPath.toStdString());
        }
        std::cout << "Saved: " << csvPath.toStdString() << "\n";
        std::cout.flush();

        // Final pass only: pulse markings written ONCE (no anchor suffix),
        // appended to the canonical CSV; then promote the .bin.
        if (finalPass) {
            const QString pulseTmp = csvPath + ".pulse.tmp";
            writeTemplateMarkingsCsv(pulseTmp.toStdString(), m_bins,
                m_subjectId.toStdString(), m_sampleRate, currentAnchor(),
                MarkingsCsvSection::PulseOnly);
            std::ifstream pin(pulseTmp.toStdString());
            std::stringstream pbuf; pbuf << pin.rdbuf();
            pin.close();
            QFile::remove(pulseTmp);
            // No suffixValueColumns() -> pulse columns stay un-suffixed.
            if (!zipCanonicalWithQ(csvPath.toStdString(), pbuf.str()))
                throw std::runtime_error("could not append pulse markings CSV into " + csvPath.toStdString());

            QFile::remove(canonicalBin);
            if (!QFile::rename(partialBin, canonicalBin))
                fprintf(stderr, "[markers] WARNING: could not promote %s -> %s\n",
                    partialBin.toStdString().c_str(), canonicalBin.toStdString().c_str());
        }
    }
    catch (const std::exception& e) {
        QMessageBox::critical(this, "Save failed",
            QString("Could not write markings for %1:\n\n%2\n\n"
                "If the CSV is open in Excel, close it and try again.")
            .arg(m_subjectId, e.what()));
        return;   // don't emit finished(); let the user retry
    }

    // Aligned-template CSV reflecting the FINAL marker edits.
    writeAlignedTemplateCsv();

    // Anchor cycle: emit reload until the last pass, then finish.
    if (!finalPass) {
        emit requestQAlignReload();
        return;
    }

    emit finished();
}