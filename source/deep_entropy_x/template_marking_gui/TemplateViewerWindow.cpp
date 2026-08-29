#include "TemplateViewerWindow.hpp"
#include "ui_TemplateViewerWindow.h"
#include "feature_marks.hpp"
#include "anchor_fit.hpp"
#include "alignment.hpp"
#include "global_intervals.hpp" 
#include "global_interval_lines.hpp"
#include "vcg_signal_average.hpp"
#include "template_generation/NormalizeFeatures.hpp"
#include "peak_finding/FilterUtils.hpp"
#include <QMessageBox>
#include <QColor>
#include <QTimer>
#include <QPixmap>
#include <QDir>
#include <QFileInfo>
#include <cmath>
#include <algorithm>
#include <limits>
#include <fstream>
#include <sstream>
#include <QFile>
#include <QCheckBox>
#include <QDockWidget>
#include <QVBoxLayout>
#include <iomanip>
#include <iostream>
#include <chrono>
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
    m_showArtPulmMarkers = false;
    m_showPpgDerivMarkers = false;

    connect(ui->show_ecg_markers, &QCheckBox::toggled, this, [this](bool on) {
        m_showEcgMarkers = on; applyMarkerVisibility();
        });
    connect(ui->show_ppg_markers, &QCheckBox::toggled, this, [this](bool on) {
        m_showPpgMarkers = on; applyMarkerVisibility();
        });
    connect(ui->show_ppg_deriv, &QCheckBox::toggled, this, [this](bool on) {
        m_showPpgDerivMarkers = on; applyMarkerVisibility();
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
    for (auto& b : m_bins) FeatureMarks::seed_all(b, m_sampleRate, m_ppgRateHz, currentAnchor());

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

    // ---- TEMPORARY INSTRUMENTATION ------------------------------------
    // Everything above this point has already printed by the time the
    // operator sees a delay ([fast-phases], [ectopic], [markers]), so the
    // hang is in one of the two calls below or in Qt's first layout of what
    // showPage builds. If showPage reports a small number and the window is
    // still slow, the cost is the paint, not this function.
    using clk = std::chrono::steady_clock;
    auto t_prev = clk::now();
    auto lap = [&t_prev](const char* what) {
        const auto now = clk::now();
        fprintf(stderr, "[viewer] %-26s %7lld ms\n", what,
            (long long)std::chrono::duration_cast<std::chrono::milliseconds>(
                now - t_prev).count());
        fflush(stderr);
        t_prev = now;
        };

    computeGlobalRefs();
    lap("computeGlobalRefs");

    showPage();
    lap("showPage");
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

    // ---- Sections 5.2-5.4: CV check + ratio/percentile normalization -----
    // Runs HERE because this is the one place TemplateBin data and a global
    // ref coexist (compute_ecg_global_ref itself is called just above). The
    // functions live in NormalizeFeatures.hpp; this drives them on real
    // per-bin data and persists the results.
    //
    // Decisions the spec left open, made explicit (override as needed):
    //  - Reference FEATURE = per-bin |R|+|S| (Option A's own basis), so the
    //    ratio and its reference are the same quantity. Computed per channel.
    //  - GLOBAL REF = all three options (A=|R|+|S|, B=QRS area, C=spatial)
    //    are written side by side rather than picking one; the ratio/pct
    //    columns use Option A to stay consistent with the reference feature.
    //  - p2/p98 for pct_scale come from THIS subject's own distribution of
    //    the per-bin ratio (per channel), matching the acceptance criterion
    //    "2nd and 98th percentiles at the extremes".
    if (!m_normOutputPath.isEmpty())
        writeNormalizationCsvs();
}

// Section 5.2-5.4 persistence. Split out of computeGlobalRefs for clarity.
void TemplateViewerWindow::writeNormalizationCsvs() {
    const std::string subj = m_subjectId.toStdString();
    const std::string dir = m_normOutputPath.toStdString();

    // Per-channel per-bin QRS reference (|R|+|S|), parallel to m_bins, NaN
    // for bins that don't yield one -- same extraction compute_ecg_global_ref
    // uses internally, just retained per bin instead of reduced to a median.
    auto perBinQrsRef = [&](int ch) {
        std::vector<double> v(m_bins.size(), std::numeric_limits<double>::quiet_NaN());
        for (size_t i = 0; i < m_bins.size(); ++i) {
            const auto& b = m_bins[i];
            if (b.bad_segment || b.bad_r_ch[ch]) continue;
            const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
            const auto& ecg = chs[ch]->ecgTemplate_raw;
            if (ecg.empty()) continue;
            const TemplateBin::MarkerSet& rmk = b.marks(AnchorType::R_PEAK);
            EcgFeatures f = computeEcgFeatures(
                ecg, rmk.p_peak_ch[ch], rmk.q_begin_ch[ch], b.r_peak_ch[ch],
                rmk.s_end_ch[ch], rmk.t_end_ch[ch], m_sampleRate);
            const double ry = normalize_features::sample_y(ecg, f.r_idx);
            const double sy = normalize_features::sample_y(ecg, f.s_idx);
            if (std::isnan(ry) || std::isnan(sy)) continue;
            v[i] = std::abs(ry) + std::abs(sy);
        }
        return v;
        };

    // ---- cv_check.csv : one row per (channel, bin) ----------------------
    // global_ref repeated per row (constant within a channel); cv_flag is
    // the per-channel record-level flag, also repeated. This is exactly the
    // long format the acceptance test reads.
    {
        std::ofstream f(dir + "/" + subj + "_cv_check.csv", std::ios::trunc);
        if (f) {
            f << "subject_id,channel,bin_index,qrs_ref_value,"
                "global_ref_A,global_ref_B,global_ref_C,cv_flag\n";
            for (int ch = 0; ch < 3; ++ch) {
                const std::vector<double> qref = perBinQrsRef(ch);
                const double grefA = normalize_features::compute_ecg_global_ref(m_bins, ch, m_sampleRate);
                const double grefB = normalize_features::compute_ecg_global_ref_area(m_bins, ch, m_sampleRate);
                const double grefC = normalize_features::compute_ecg_global_ref_spatial(m_bins);
                const bool flag = normalize_features::cv_flag(qref, grefA);
                for (size_t i = 0; i < qref.size(); ++i) {
                    f << subj << ",CH" << (ch + 1) << ',' << i << ',';
                    if (!std::isnan(qref[i])) f << qref[i];
                    f << ',';
                    if (!std::isnan(grefA)) f << grefA; f << ',';
                    if (!std::isnan(grefB)) f << grefB; f << ',';
                    if (!std::isnan(grefC)) f << grefC; f << ',';
                    f << (flag ? 1 : 0) << '\n';
                }
            }
        }
    }

    // ---- feature_norm.csv : one row per (channel, bin) ------------------
    // ratio = ratio_norm(qrsRef, grefA); feature_norm = pct_scale against
    // this channel's own p2/p98 of the finite ratios. By construction the
    // p2 and p98 rows land at 0 and 100 (clamped), which is the property the
    // acceptance test checks.
    {
        std::ofstream f(dir + "/" + subj + "_feature_norm.csv", std::ios::trunc);
        if (f) {
            f << "subject_id,channel,bin_index,ratio,feature_norm\n";
            for (int ch = 0; ch < 3; ++ch) {
                const std::vector<double> qref = perBinQrsRef(ch);
                const double grefA = normalize_features::compute_ecg_global_ref(m_bins, ch, m_sampleRate);

                std::vector<double> ratios(qref.size(), std::numeric_limits<double>::quiet_NaN());
                std::vector<double> finite;
                for (size_t i = 0; i < qref.size(); ++i) {
                    if (std::isnan(qref[i])) continue;
                    ratios[i] = normalize_features::ratio_norm(qref[i], grefA);
                    if (!std::isnan(ratios[i])) finite.push_back(ratios[i]);
                }
                // p2 / p98 of this channel's own ratio distribution.
                double p2 = std::nan(""), p98 = std::nan("");
                if (finite.size() >= 2) {
                    std::sort(finite.begin(), finite.end());
                    auto pctl = [&](double q) {
                        const double idx = (q / 100.0) * (finite.size() - 1);
                        const size_t lo = static_cast<size_t>(std::floor(idx));
                        const size_t hi = static_cast<size_t>(std::ceil(idx));
                        const double fr = idx - lo;
                        return finite[lo] * (1.0 - fr) + finite[hi] * fr;
                        };
                    p2 = pctl(2.0);
                    p98 = pctl(98.0);
                }
                for (size_t i = 0; i < ratios.size(); ++i) {
                    f << subj << ",CH" << (ch + 1) << ',' << i << ',';
                    if (!std::isnan(ratios[i])) f << ratios[i];
                    f << ',';
                    const double fn = normalize_features::pct_scale(ratios[i], p2, p98);
                    if (!std::isnan(fn)) f << fn;
                    f << '\n';
                }
            }
        }
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
    // ---- TEMPORARY INSTRUMENTATION ------------------------------------
    using clk = std::chrono::steady_clock;
    auto t_prev = clk::now();
    auto lap = [&t_prev](const char* what) {
        const auto now = clk::now();
        fprintf(stderr, "[showPage] %-24s %7lld ms\n", what,
            (long long)std::chrono::duration_cast<std::chrono::milliseconds>(
                now - t_prev).count());
        fflush(stderr);
        t_prev = now;
        };

    clearPlots();
    lap("clearPlots");

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
        // One extra row for the VCG panel, but only if some bin on this page
        // can actually produce one (all three ECG channels present with an R
        // column). Without this the row count equals the lead count and there
        // is no row index left for the VCG to occupy.
        bool anyVcg = false;
        for (int k = start; k < end && !anyVcg; ++k)
            anyVcg = !vcg_avg::derivedTraceOnChannelAxis(
                m_bins[k], 0, vcg::DerivedLead::VectorMagnitude).empty();
        gridRows = m_maxLeads + (anyVcg ? 1 : 0);
        // Note: this probe throws its trace away and the per-bin loop below
        // recomputes the same thing for the same bins. If this lap is large,
        // that duplication is the first thing to remove.
        lap("anyVcg probe");
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

        // VCG needs all three channels; with fewer, the trace comes back empty
        // and no row is reserved, so the leads keep the full height.
        double vcgRCol = -1.0;
        const std::vector<double> vcgTrace = vcg_avg::derivedTraceOnChannelAxis(
            b, 0, vcg::DerivedLead::VectorMagnitude, vcg::kIdentity, &vcgRCol);
        const bool vcgRowWanted = !compact && !vcgTrace.empty()
            && gridRows > m_maxLeads;

        const auto gi_intervals = global_intervals::computeGlobalIntervals(
            b, m_currentAnchor, m_sampleRate,
            global_intervals::MarkerSource::USER);

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

            // Traces only. Every marker and autodetect column arrives via the
            // single applyBinToWidget() call below.
            pw->setData(ppgN, ppgIqr, ecgN, ecgIqr, rPeak,
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

            // Arterial traces for marker geometry/bounds. The arterial markers
            // themselves come from applyBinToWidget() with all the others.
            pw->setArterialTraces(abpN, artN, artPN,
                b.abpTemplate_iqr, b.artTemplate_iqr, b.artPulmTemplate_iqr);

            // Seed every bar + every autodetect column, in one call, after all
            // traces are in place (the glyph capture needs them).
            applyBinToWidget(pw, b);

            pw->setReferenceLines(global_interval_lines::forChannel(b, gi_intervals, c));

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
                // The last lead no longer swallows the remaining rows: the
                // bottom row belongs to the VCG panel added after this loop.
                int rowspan = (li == (int)leads.size() - 1)
                    ? std::max(1, gridRows - li - (vcgRowWanted ? 1 : 0)) : 1;
                ui->plotGrid->addWidget(pw, li, i, rowspan, 1);
                usedRows = std::max(usedRows, li + rowspan);
                usedCols = std::max(usedCols, i + 1);
            }

            m_allPlots.push_back(pw);
            group.push_back(pw);
        }

        // ------------------------------------------------------------------
        // VCG panel: the bottom row of each bin's column, under lead 3.
        // ------------------------------------------------------------------
        // The trace is laid out on ch1's COLUMN axis by
        // derivedTraceOnChannelAxis, so it shares the x axis of the lead
        // panels above and lines up with them sample for sample -- while every
        // channel is still sampled at its own r_col internally, which is what
        // keeps the combination per-instant.
        //
        // Display only: no markers are set and marker signals are not
        // connected, so nothing here is draggable. The VCG is derived from the
        // three leads, so marking it would create a fourth set of fiducials
        // with no channel of its own to store them in.
        if (vcgRowWanted) {
            const int vcgRow = gridRows - 1;
            auto* vp = new BinPlotWidget(gi, vcgRow, "VCG", this);

            static const std::vector<double> emptyVec;
            vp->setChannelRate(BinPlotWidget::Channel::Ecg, m_sampleRate);
            vp->setData(emptyVec, emptyVec, vcgTrace, emptyVec,
                (vcgRCol >= 0.0) ? vcgRCol : 0.0, 0, 0);
            vp->setHasPPG(false);
            vp->setShowPpgTrace(false);
            vp->setShowEcgMarkers(false);
            vp->setShowPpgMarkers(false);
            vp->setShowPpgDerivMarkers(false);
            vp->setShowAbpMarkers(false);
            vp->setShowArtMarkers(false);
            vp->setShowArtPulmMarkers(false);

            // Same global boundaries as the leads above, converted onto ch1's
            // axis since that is the axis this trace is drawn on.
            vp->setReferenceLines(
                global_interval_lines::forChannel(b, gi_intervals, 0));

            ui->plotGrid->addWidget(vp, vcgRow, i, 1, 1);
            usedRows = std::max(usedRows, vcgRow + 1);
            usedCols = std::max(usedCols, i + 1);

            m_allPlots.push_back(vp);
            group.push_back(vp);
        }

        m_binPlots[i] = std::move(group);
        lap("  one bin's widgets");
    }
    lap("all bin widgets");

    // Equal stretch on every used row/column => equal-width, equal-height
    // cells that together fill the whole plot area. Combined with each
    // widget scaling its trace to its cell width, every bin window ends up
    // the same on-screen length.
    for (int c = 0; c < usedCols; ++c) ui->plotGrid->setColumnStretch(c, 1);
    for (int r = 0; r < usedRows; ++r) ui->plotGrid->setRowStretch(r, 1);

    applyMarkerVisibility();
    updatePageControls();
    lap("stretch+visibility+controls");
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
// Merge an ordered list of sidecar CSV files into one canonical file.
// files[0] is the base (kept whole, with its keys); each subsequent file
// contributes only its value columns (first 3 key columns stripped), appended
// to every row -- same column convention as zipCanonicalWithQ, generalized to
// N inputs and done ONCE at the end instead of a growing read-modify-write per
// pass. Row counts must match across all files. Returns true on success.
static bool mergeSidecarCsvs(const std::string& canonicalPath,
    const std::vector<std::string>& sidecarPaths)
{
    if (sidecarPaths.empty()) return false;

    auto readLines = [](const std::string& p) {
        std::vector<std::string> lines; std::string ln;
        std::ifstream in(p);
        if (!in) return lines;
        while (std::getline(in, ln)) {
            if (!ln.empty() && ln.back() == '\r') ln.pop_back();
            lines.push_back(ln);
        }
        return lines;
        };
    auto stripFirstThree = [](const std::string& line) -> std::string {
        int commas = 0;
        for (size_t i = 0; i < line.size(); ++i)
            if (line[i] == ',' && ++commas == 3) return line.substr(i + 1);
        return {};
        };

    std::vector<std::string> merged = readLines(sidecarPaths[0]);
    if (merged.empty()) {
        fprintf(stderr, "[tmplcsv] merge: base sidecar empty/missing: %s\n", sidecarPaths[0].c_str());
        return false;
    }
    for (size_t s = 1; s < sidecarPaths.size(); ++s) {
        std::vector<std::string> next = readLines(sidecarPaths[s]);
        if (next.size() != merged.size()) {
            fprintf(stderr, "[tmplcsv] merge: row mismatch (%zu vs %zu) for %s -- skipping\n",
                next.size(), merged.size(), sidecarPaths[s].c_str());
            continue;   // skip a bad sidecar rather than abort the whole merge
        }
        for (size_t i = 0; i < merged.size(); ++i) {
            const std::string tail = stripFirstThree(next[i]);
            if (!tail.empty()) { merged[i] += ','; merged[i] += tail; }
        }
    }

    std::ofstream out(canonicalPath, std::ios::trunc);
    if (!out) {
        fprintf(stderr, "[tmplcsv] merge: cannot write %s\n", canonicalPath.c_str());
        return false;
    }
    for (const auto& l : merged) out << l << '\n';
    return out.good();
}

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
    // One row per sample column per bin. Every value column is emitted twice,
    // raw and normalized, plus its IQR band in both forms; marker positions
    // ride along as one-hot flags on the row matching their sample index.
    if (m_bins.empty()) return;

    QDir outDir(m_templateDir);
    if (!outDir.exists()) outDir.mkpath(".");
    // Names the canonical merge target. Used only in the log/error messages
    // below -- this function writes a per-anchor sidecar, never `path` itself.
    const QString path = outDir.filePath(m_subjectId + "_template.csv");

    std::ostringstream f;
    const double toMs = (m_sampleRate > 0.0) ? 1000.0 / m_sampleRate : 1.0;

    // ---- Column lists (shared by the header pass and the row loop) ---------
    static const char* CHANS[] = {
        "ch1", "ch2", "ch3", "ppg", "abp", "art", "art_pulm"
    };
    constexpr int num_chans = static_cast<int>(std::size(CHANS));

    // Per ECG channel, in emission order. r_peak is auto-only: there is no
    // user bar for it, so it contributes one column where the others
    // contribute two. kRPeak is derived from the array rather than written as
    // a literal, because the index was previously hardcoded as `4` in two
    // places and inserting a marker ahead of it would have silently moved the
    // auto-only column onto the wrong landmark.
    static const char* ECG_MARKERS[] = {
        "p_begin", "p_peak", "q_begin", "q_peak", "r_peak", "s_peak",
        "s_end",   "t_begin", "t_end"
    };
    constexpr int kNumEcgMarkers = static_cast<int>(std::size(ECG_MARKERS));
    constexpr int kRPeak = 4;
    static_assert(kNumEcgMarkers == 9, "ECG marker count changed; check kRPeak");

    // Pulse marker groups. Counts differ by channel: PPG carries the two
    // interpolated upslope crossings (t50/t80), the arterial channels do not.
    static const char* PPG_MARKERS[] = {
        "ppg_onset", "ppg_t50", "ppg_peak",
        "ppg_dicr", "ppg_peak2", "ppg_t80", "ppg_end"
    };
    static const char* ABP_MARKERS[] = {
        "abp_onset", "abp_peak", "abp_dicr", "abp_peak2", "abp_end"
    };
    static const char* ART_MARKERS[] = {
        "art_onset", "art_peak", "art_dicr", "art_peak2", "art_end"
    };
    static const char* ARTP_MARKERS[] = {
        "art_pulm_onset", "art_pulm_peak", "art_plm_dicr",
        "art_pulm_peak2", "art_pulm_end"
    };
    constexpr int kNumPpgMarkers = static_cast<int>(std::size(PPG_MARKERS));
    constexpr int kNumArterialMarkers = static_cast<int>(std::size(ABP_MARKERS));

    // Per-ECG-channel autodetect glyph group. p_wave/q_onset/r_wave are the
    // stored autodetect columns, emitted directly so no parallel recompute can
    // disagree with them; t_peak is the only derived one.
    static const char* ECG_GLYPHS[] = { "p_wave", "q_onset", "r_wave", "t_peak" };

    // ---- Header ------------------------------------------------------------
    f << "file_id,bin_num,x_ms";
    for (const char* n : CHANS) {
        f << ',' << n << "_raw_mv"
            << ',' << n << "_norm"
            << ',' << n << "_raw_std"
            << ',' << n << "_norm_std";
    }
    for (int c = 1; c <= 3; ++c) {
        for (int k = 0; k < kNumEcgMarkers; ++k) {
            f << ',' << ECG_MARKERS[k] << "_ch" << c << "_auto";
            if (k != kRPeak)
                f << ',' << ECG_MARKERS[k] << "_ch" << c << "_auto";
        }
    }
    auto emitPulseHeaderGroup = [&](auto const& group) {
        for (const char* m : group)
            f << ',' << m << "_auto"
            << ',' << m << "_user";
        };
    emitPulseHeaderGroup(PPG_MARKERS);
    emitPulseHeaderGroup(ABP_MARKERS);
    emitPulseHeaderGroup(ART_MARKERS);
    emitPulseHeaderGroup(ARTP_MARKERS);

    // Autodetected computed feature locations (no user bar).
    for (int gc = 1; gc <= 3; ++gc) {
        char gb[64];
        for (const char* g : ECG_GLYPHS) {
            std::snprintf(gb, sizeof gb, "%s_ch%d", g, gc);
            f << ',' << gb << "_auto";
        }
    }
    // Driven off ppg_and_artpulse_automated_markers, the same table the row
    // loop emits from. This used to be a hand-written list of four names while
    // the table had grown to fifteen entries, so the file carried eleven
    // unnamed columns and every header-keyed reader was reading the wrong
    // field from here to the end of the row.
    for (const auto& gl : ppg_and_artpulse_automated_markers)
        f << ',' << gl.name << "_auto";
    f << '\n';

    f << std::setprecision(10);

    // Both ecg_template_raw_iqr and every *_template_iqr field are
    // pre-ref-division at build time (see CreatePPGTemplates.hpp /
    // build_pulse_template_pair_windowed), so the only remaining step for any
    // channel is the same scalar /ref used for its mean trace --
    // normalize_features::scale_array_by_ref is the one place that happens.

    // ---- Row loop ----------------------------------------------------------
    for (size_t bi = 0; bi < m_bins.size(); ++bi) {
        const TemplateBin& b = m_bins[bi];
        const TemplateBin::MarkerSet& umk = b.marks(currentAnchor());

        // One description per value column, in CHANS order. The seven channels
        // differ only in which normalizer they take and which global reference
        // scales their IQR, so they are described rather than transcribed:
        // the previous form repeated four near-identical statements per channel
        // and the ordering of CHANS was only implicitly matched.
        struct Src {
            const std::vector<double>* raw;
            const std::vector<double>* rawIqr;
            double ref;      ///< global reference for this channel's scaling
            int    idx;      ///< channel index handed to the normalizer
            bool   isEcg;    ///< selects the normalizer
            int    onset;    ///< pulse channels only: alignment onset
        };
        const Src src[num_chans] = {
            { &b.ch1.ecgTemplate_raw, &b.ch1.ecg_template_raw_iqr, m_ecgGlobalRef[0],   0, true,  0 },
            { &b.ch2.ecgTemplate_raw, &b.ch2.ecg_template_raw_iqr, m_ecgGlobalRef[1],   1, true,  0 },
            { &b.ch3.ecgTemplate_raw, &b.ch3.ecg_template_raw_iqr, m_ecgGlobalRef[2],   2, true,  0 },
            { &b.ppgTemplate,         &b.ppg_template_iqr,         m_pulseGlobalRef[0], 0, false, b.ppg_onset },
            { &b.abpTemplate,         &b.abpTemplate_iqr,          m_pulseGlobalRef[1], 1, false, b.abp_onset },
            { &b.artTemplate,         &b.artTemplate_iqr,          m_pulseGlobalRef[2], 2, false, b.art_onset },
            { &b.artPulmTemplate,     &b.artPulmTemplate_iqr,      m_pulseGlobalRef[3], 3, false, b.art_pulm_onset },
        };

        std::vector<double> norm[num_chans], normIqr[num_chans];
        for (int k = 0; k < num_chans; ++k) {
            norm[k] = src[k].isEcg
                ? normalizeEcgTrace(*src[k].raw, src[k].idx)
                : normalize_ppg_or_similar(*src[k].raw, src[k].onset, src[k].idx);
            normIqr[k] = normalize_features::scale_array_by_ref(*src[k].rawIqr, src[k].ref);
        }

        // Row span = longest trace; all traces start at row 0. Marker indices
        // are NOT considered: a landmark past the end of every trace has no row
        // to flag, which is a template-construction problem rather than
        // something this writer can represent.
        int hiRow = 0;
        for (const Src& s : src) hiRow = std::max(hiRow, (int)s.raw->size());
        if (hiRow <= 0) continue;

        // Computed Q/S peaks per ECG channel, from the auto bars and the user
        // bars separately.
        const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
        EcgFeatures ftAuto[3], ftUser[3];
        for (int c = 0; c < 3; ++c) {
            const auto& ecg = chs[c]->ecgTemplate_raw;
            ftAuto[c] = computeEcgFeatures(ecg,
                (int)std::lround(b.p_peak_auto_ch[c]), (int)std::lround(b.q_begin_auto_ch[c]),
                (int)std::lround(b.r_peak_auto_ch[c]), (int)std::lround(b.s_end_auto_ch[c]),
                (int)std::lround(b.t_end_auto_ch[c]), m_sampleRate);
            ftUser[c] = computeEcgFeatures(ecg,
                umk.p_peak_ch[c], umk.q_begin_ch[c], b.r_peak_ch[c],
                umk.s_end_ch[c], umk.t_end_ch[c], m_sampleRate);
        }

        // ECG marker positions, indices aligned with ECG_MARKERS.
        double ecgAuto[3][kNumEcgMarkers], ecgUser[3][kNumEcgMarkers];
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

        // Pulse marker positions, order matching the *_MARKERS arrays.
        // t50/t80 are reactive: bracketed by onset/peak/end, never stored. The
        // autodetect column brackets with the *_auto bars and the user column
        // with the user bars, both through the same FeatureMarks::reactive_ppg
        // the on-screen glyph uses -- so what is plotted is what is exported.
        const FeatureMarks::ReactivePpg rxPpgAuto = FeatureMarks::reactive_ppg(
            b.ppgTemplate, b.ppg_onset_auto, b.ppg_peak_auto, b.ppg_end_auto);
        const FeatureMarks::ReactivePpg rxPpgUser = FeatureMarks::reactive_ppg(
            b.ppgTemplate, b.ppg_onset, b.ppg_peak, b.ppg_end);
        // double, not int: t50/t80 are interpolated crossings and the stored
        // fields promote without loss. Rounding happens once, in emitLoc.
        const double ppgAuto[kNumPpgMarkers] = {
            (double)b.ppg_onset_auto, rxPpgAuto.t50, (double)b.ppg_peak_auto,
            (double)b.ppg_dicrotic_auto, (double)b.ppg_peak2_auto, rxPpgAuto.t80,
            (double)b.ppg_end_auto };
        const double ppgUser[kNumPpgMarkers] = {
            (double)b.ppg_onset, rxPpgUser.t50, (double)b.ppg_peak,
            (double)b.ppg_dicrotic, (double)b.ppg_peak2, rxPpgUser.t80,
            (double)b.ppg_end };
        const int abpAuto[kNumArterialMarkers] = { b.abp_onset_auto, b.abp_peak_auto,
            b.abp_dicrotic_auto, b.abp_peak2_auto, b.abp_end_auto };
        const int abpUser[kNumArterialMarkers] = { b.abp_onset, b.abp_peak,
            b.abp_dicrotic, b.abp_peak2, b.abp_end };
        const int artAuto[kNumArterialMarkers] = { b.art_onset_auto, b.art_peak_auto,
            b.art_dicrotic_auto, b.art_peak2_auto, b.art_end_auto };
        const int artUser[kNumArterialMarkers] = { b.art_onset, b.art_peak,
            b.art_dicrotic, b.art_peak2, b.art_end };
        const int artpAuto[kNumArterialMarkers] = { b.art_pulm_onset_auto, b.art_pulm_peak_auto,
            b.art_pulm_dicrotic_auto, b.art_pulm_peak2_auto, b.art_pulm_end_auto };
        const int artpUser[kNumArterialMarkers] = { b.art_pulm_onset, b.art_pulm_peak,
            b.art_pulm_dicrotic, b.art_pulm_peak2, b.art_pulm_end };

        // Derived T peak for the autodetect glyph group, bracketed by the AUTO
        // T-begin/T-end so it matches that group's name. Brackets stay
        // fractional; compute_t_peak takes doubles.
        double tPeakAutoGlyph[3];
        for (int gc = 0; gc < 3; ++gc)
            tPeakAutoGlyph[gc] = FeatureMarks::compute_t_peak(
                chs[gc]->ecgTemplate_raw,
                b.t_begin_auto_ch[gc], b.t_end_auto_ch[gc]);

        auto emitVal = [&](const std::vector<double>& v, int j) {
            f << ',';
            if (j >= 0 && j < (int)v.size() && !std::isnan(v[j])) f << v[j];
            };
        // Emit ",1" if the row is this marker's row, ",<blank>" otherwise.
        //
        // ROUNDING BOUNDARY, and an unavoidable one: this column is a one-hot
        // flag per SAMPLE ROW, so a landmark at 104.37 has to be attributed to
        // a row and the format has no fractional representation. Rounded to the
        // nearest row explicitly, in one place, rather than truncated by
        // implicit conversions at each call site. The millisecond columns
        // elsewhere in this file carry the fraction.
        auto emitLoc = [&](double markerIdx, int row) {
            f << ',';
            if (markerIdx >= 0.0 && (int)std::lround(markerIdx) == row) f << '1';
            };
        for (int row = 0; row < hiRow; ++row) {
            f << m_subjectId.toStdString() << ',' << bi << ',' << (row * toMs);
            for (int k = 0; k < num_chans; ++k) {
                emitVal(*src[k].raw, row);
                emitVal(norm[k], row);
                emitVal(*src[k].rawIqr, row);
                emitVal(normIqr[k], row);
            }
            // ECG: per channel, per marker, auto then user (r_peak auto only).
            for (int c = 0; c < 3; ++c)
                for (int k = 0; k < kNumEcgMarkers; ++k) {
                    emitLoc(ecgAuto[c][k], row);
                    if (k != kRPeak) emitLoc(ecgUser[c][k], row);
                }
            // Pulse groups, auto and user interleaved per marker.
            auto emitPair = [&](const auto& a, const auto& u, int n) {
                for (int k = 0; k < n; ++k) { emitLoc(a[k], row); emitLoc(u[k], row); }
                };
            emitPair(ppgAuto, ppgUser, kNumPpgMarkers);
            emitPair(abpAuto, abpUser, kNumArterialMarkers);
            emitPair(artAuto, artUser, kNumArterialMarkers);
            emitPair(artpAuto, artpUser, kNumArterialMarkers);
            // Autodetect glyphs: per ECG channel, then the pulse table.
            for (int gc = 0; gc < 3; ++gc) {
                emitLoc(b.p_peak_auto_ch[gc], row);
                emitLoc(b.q_begin_auto_ch[gc], row);
                emitLoc(b.r_peak_auto_ch[gc], row);
                emitLoc(tPeakAutoGlyph[gc], row);
            }
            for (const auto& gl : ppg_and_artpulse_automated_markers)
                emitLoc(b.*gl.idx, row);
            f << '\n';
        }
    }

    // Serialize, suffix the value columns with the pass tag, and write this
    // pass's sidecar. Row-key columns (file_id/bin_num/x_ms) are never
    // suffixed, so the sidecars line up on the merge in save_bin_and_csv.
    std::string content = f.str();
    const size_t nl = content.find('\n');
    if (nl == std::string::npos) {
        fprintf(stderr, "[tmplcsv] malformed content (no newline) for %s\n",
            path.toStdString().c_str());
        return;
    }
    const std::string header = content.substr(0, nl);
    const std::string body = content.substr(nl);   // includes the leading '\n'
    const std::string suffix = "_" + m_anchorLabel.toStdString();

    const QString sidecar = outDir.filePath(
        m_subjectId + "_template" + suffix.c_str() + ".csv");
    std::ofstream out(sidecar.toStdString(), std::ios::trunc);
    if (!out) {
        fprintf(stderr, "[tmplcsv] cannot open sidecar %s\n",
            sidecar.toStdString().c_str());
        return;
    }
    out << suffixValueColumns(header, suffix) << body;
    out.close();
    fprintf(stderr, "[tmplcsv] wrote sidecar %s\n", sidecar.toStdString().c_str());
}

void TemplateViewerWindow::applyMarkerVisibility() {
    for (auto* pw : m_allPlots) {
        pw->setShowEcgMarkers(m_showEcgMarkers);
        pw->setShowPpgMarkers(m_showPpgMarkers);
        pw->setShowPpgDerivMarkers(m_showPpgDerivMarkers);
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

// ============================================================================
// The single seeding path from a TemplateBin into a plot widget.
//
// Every draggable bar and every frozen autodetect column is written here, from
// one bin, in one pass. There is no second list anywhere: showPage() calls this
// on build and refreshBinMarkers() calls it on every later change, so a bar and
// its glyph can never be seeded from different places (which is how the P-onset
// bar ended up unseeded while its glyph was drawn).
//
// setAuto() MUST stay last: it performs the glyph capture, and the frozen
// snapshot reads the R bar.
// ============================================================================
void TemplateViewerWindow::applyBinToWidget(BinPlotWidget* pw, const TemplateBin& b) {
    const int c = pw->leadIndex();
    const TemplateBin::MarkerSet& mk = b.marks(currentAnchor());

    // ---- draggable bars ----------------------------------------------------
    pw->setMarker(BinPlotWidget::EcgPBegin, mk.p_begin_ch[c]);
    pw->setMarker(BinPlotWidget::EcgPPeak, mk.p_peak_ch[c]);
    pw->setMarker(BinPlotWidget::EcgQBegin, mk.q_begin_ch[c]);
    pw->setMarker(BinPlotWidget::EcgRPeak, b.r_peak_ch[c]);   // auto-only, no bar drawn
    pw->setMarker(BinPlotWidget::EcgSEnd, mk.s_end_ch[c]);
    pw->setMarker(BinPlotWidget::EcgTBegin, mk.t_begin_ch[c]);
    pw->setMarker(BinPlotWidget::EcgTEnd, mk.t_end_ch[c]);

    pw->setMarker(BinPlotWidget::PpgOnset, b.ppg_onset);
    pw->setMarker(BinPlotWidget::PpgPeak, b.ppg_peak);
    pw->setMarker(BinPlotWidget::PpgDicrotic, b.ppg_dicrotic);
    pw->setMarker(BinPlotWidget::PpgPeak2, b.ppg_peak2);
    pw->setMarker(BinPlotWidget::PpgEnd, b.ppg_end);
    // T50/T80 are neither drawn nor draggable -- they're reactive glyphs now.
    // Kept in sync anyway so the enum entries never hold a stale position.
    pw->setMarker(BinPlotWidget::PpgT50, b.ppg_t50);
    pw->setMarker(BinPlotWidget::PpgT80, b.ppg_t80);

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
    pw->setAuto(b);   // last: captures the glyph snapshot and repaints
}

void TemplateViewerWindow::refreshBinMarkers(int binIdx) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    for (int li = 0; li < (int)m_pageGlobalIdx.size(); ++li) {
        if (m_pageGlobalIdx[li] != binIdx) continue;
        for (auto* pw : m_binPlots[li]) applyBinToWidget(pw, m_bins[binIdx]);
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
        // Track the drag: the touched store follows the bar to its final
        // position so confirmedIndex reflects where the operator left it.
        if (newIdx >= 0)
            m_touchedMarks[touchKey(binIdx, leadIdx, marker)] = newIdx;
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
    // Focus activation = the operator clicked this bar. Record it as "touched"
    // at position col; logBoundaryTrainingAtSave reads this to fill
    // confirmedIndex for the matching landmark.
    if (binIdx >= 0 && leadIdx >= 0 && col >= 0)
        m_touchedMarks[touchKey(binIdx, leadIdx, marker)] = col;
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

// AnchorType -> short name for the boundary log's `anchor` column.
static const char* anchorName_boundary(AnchorType a) {
    switch (a) {
    case AnchorType::P_ONSET: return "P_ONSET";
    case AnchorType::P_PEAK:  return "P_PEAK";
    case AnchorType::Q_ONSET: return "Q_ONSET";
    case AnchorType::R_PEAK:  return "R_PEAK";
    case AnchorType::J_POINT: return "J_POINT";
    case AnchorType::T_PEAK:  return "T_PEAK";
    }
    return "UNKNOWN";
}

// Log boundary training data at save time. For EVERY template (bin/lead), log
// each landmark (Q-onset, J-point, P-onset, T-end). auto_detect comes from the
// bin's auto-detected glyph fields (*_auto_ch), which hold every landmark on
// every template regardless of anchor pass. expert_mark = the user marker for
// that landmark iff the operator moved it away from auto (else null).
void TemplateViewerWindow::logBoundaryTrainingAtSave() {
    using boundary_training::Landmark;

    struct Target { Landmark lm; };
    // S-end omitted (same feature as J-point). Each landmark's auto position is
    // read from the bin below; its expert mark (if any) from the current pass.
    const Target targets[] = {
        { Landmark::Q_ONSET  },
        { Landmark::J_POINT  },
        { Landmark::P_ONSET  },
        { Landmark::T_OFFSET },   // T-end
    };

    // auto-detected position (glyph field) for a landmark on a lead. The
    // *_auto_ch fields are double (sub-sample); rounded here since this seeds
    // an integer segment window.
    auto autoPosOf = [](const TemplateBin& tb, Landmark lm, int lead) -> int {
        switch (lm) {
        case Landmark::Q_ONSET:  return (int)std::lround(tb.q_begin_auto_ch[lead]);
        case Landmark::J_POINT:  return (int)std::lround(tb.s_end_auto_ch[lead]);   // J-point == S-end field
        case Landmark::P_ONSET:  return (int)std::lround(tb.p_begin_auto_ch[lead]);
        case Landmark::T_OFFSET: return (int)std::lround(tb.t_end_auto_ch[lead]);
        default: return -1;
        }
        };
    // Map each logged landmark to its BinPlotWidget marker id (for the touch key).
    auto markerIdOf = [](Landmark lm) -> int {
        switch (lm) {
        case Landmark::Q_ONSET:  return BinPlotWidget::EcgQBegin;
        case Landmark::J_POINT:  return BinPlotWidget::EcgSEnd;
        case Landmark::P_ONSET:  return BinPlotWidget::EcgPBegin;
        case Landmark::T_OFFSET: return BinPlotWidget::EcgTEnd;
        default: return -1;
        }
        };

    const int half = std::max(1, (int)std::lround(0.100 * m_sampleRate)); // +/-100 ms
    int written = 0, failed = 0;

    // Slicing puts R at r_col = percent_interval_preceeding_rpeak * RR, so
    // RR_samples = r_col / that fraction. Used for heart rate.
    constexpr double kPreRFrac = alignment::percent_interval_preceeding_rpeak;

    for (int i = 0; i < (int)m_bins.size(); ++i) {
        const TemplateBin& b = m_bins[i];
        const std::vector<double>* chs[3] = {
            &b.ch1.ecgTemplate_raw, &b.ch2.ecgTemplate_raw, &b.ch3.ecgTemplate_raw };
        const int rcol[3] = { b.ch1.r_col_raw, b.ch2.r_col_raw, b.ch3.r_col_raw };
        for (int lead = 0; lead < 3; ++lead) {
            const std::vector<double>& sig = *chs[lead];

            // Heart rate (bpm) = 60000 / RR(ms); RR from this lead's r_col.
            double heartRate = 0.0;
            if (rcol[lead] > 0 && m_sampleRate > 0.0) {
                const double rrSamples = rcol[lead] / kPreRFrac;
                const double rrMs = rrSamples / m_sampleRate * 1000.0;
                if (rrMs > 0.0) heartRate = 60000.0 / rrMs;
            }
            // QRS duration (ms) = distance between q_begin and s_end glyphs.
            double qrsMs = 0.0;
            {
                const int q = (int)std::lround(b.q_begin_auto_ch[lead]);
                const int s = (int)std::lround(b.s_end_auto_ch[lead]);
                if (q >= 0 && s >= 0 && m_sampleRate > 0.0)
                    qrsMs = std::abs(s - q) / m_sampleRate * 1000.0;
            }

            for (const Target& tgt : targets) {
                const int autoPos = autoPosOf(b, tgt.lm, lead);
                if (autoPos < 0 || autoPos >= (int)sig.size()) continue;

                const int lo = std::max(0, autoPos - half);
                const int hi = std::min((int)sig.size(), autoPos + half);
                if (hi - lo < 2) continue;

                // confirmedIndex = the operator's clicked position (segment-
                // relative) if this landmark was touched (focus activated on
                // its bar), else -1 => blank. The row is logged either way.
                int confirmed = -1;
                const int mid = markerIdOf(tgt.lm);
                auto it = m_touchedMarks.find(touchKey(i, lead, mid));
                if (it != m_touchedMarks.end()) confirmed = it->second - lo;

                boundary_training::BoundaryTrainingRecord rec;
                rec.segment.assign(sig.begin() + lo, sig.begin() + hi);
                rec.confirmedIndex = confirmed;
                // fit from anchor_fit: fit-and-select on this landmark's window.
                const anchor_fit::FitResult fit = anchor_fit::selectAnchorModel(sig, lo, hi - 1);
                rec.fitType = fit.type;
                rec.fitRSS = fit.rss;
                rec.individualID = m_subjectId.toStdString();
                rec.bbb = false;                 // left blank for now
                rec.heartRate = heartRate;
                rec.qrsDurationMs = qrsMs;
                if (boundary_training::logBoundary(
                    m_boundaryLog,
                    anchorName_boundary(currentAnchor()),// anchor = this template's alignment
                    tgt.lm, rec)) ++written; else ++failed;
            }
        }
    }
    std::fprintf(stderr,
        "[boundary_log] save: dir='%s' wrote=%d failed=%d\n",
        m_boundaryLog.dir.c_str(), written, failed);
}

void TemplateViewerWindow::save_bin_and_csv() {
    captureCurrentPage();   // snapshot the page being left on Finish

    QDir binDir(m_markingPath);
    QDir csvDir(m_markingPath);
    QDir vcgDir(m_vcgOutputPath);
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
        logBoundaryTrainingAtSave();

        // ECG markings CSV: each pass writes its OWN per-anchor sidecar
        // (<id>_template_markings.<ANCHOR>.csv) with suffixed columns. No
        // growing zip per pass; all sidecars are merged into the canonical
        // <id>_template_markings.csv once, at the final pass.
        // VCG loop features: one row per bin, <id>_vcg.csv, beside this CSV.
        // Written every pass -- it is a standalone file, not a per-anchor
        // sidecar that needs merging, and the markers it measures against
        // change on every pass, so the latest write is the one that matters.
        {
            std::vector<vcg_avg::BinFeatures> vcgRows;
            vcgRows.reserve(m_bins.size());
            for (int vi = 0; vi < (int)m_bins.size(); ++vi) {
                // Prefer the operator's markers. On the R pass (and on any bin
                // not yet marked for this anchor) there are none, so fall back
                // to the autodetected ones rather than emitting an empty row --
                // the features are still measurable, and the source is
                // recorded so a row is never ambiguous about where its
                // boundaries came from.
                auto vgi = global_intervals::computeGlobalIntervals(
                    m_bins[vi], currentAnchor(), m_sampleRate,
                    global_intervals::MarkerSource::USER);
                bool fromAuto = false;
                if (!vgi.valid) {
                    vgi = global_intervals::computeGlobalIntervals(
                        m_bins[vi], currentAnchor(), m_sampleRate,
                        global_intervals::MarkerSource::AUTO);
                    fromAuto = vgi.valid;
                }
                auto vrow = vcg_avg::analyzeBinFromTemplates(
                    vi, m_bins[vi], vgi, m_sampleRate);
                if (vrow.valid && fromAuto) vrow.note = "auto markers";
                vcgRows.push_back(vrow);
            }
            const bool vok = vcg_avg::writeVcgCsv(
                vcgDir.absolutePath().toStdString(),
                m_subjectId.toStdString(), vcgRows);
            std::cout << (vok ? "Saved: " : "FAILED: ")
                << vcgDir.absolutePath().toStdString() << "/"
                << m_subjectId.toStdString() << "_vcg.csv\n";
        }

        const QString csvPath = csvDir.absolutePath() + "/" + m_subjectId + "_template_markings.csv";
        const QString tsuffix = "_" + m_anchorLabel;
        const QString ecgSidecar = csvDir.absolutePath() + "/" + m_subjectId + "_template_markings" + tsuffix + ".csv";
        {
            const QString tmpPath = ecgSidecar + ".tmp";
            writeTemplateMarkingsCsv(tmpPath.toStdString(), m_bins,
                m_subjectId.toStdString(), m_sampleRate, currentAnchor(),
                MarkingsCsvSection::EcgOnly);
            std::ifstream tin(tmpPath.toStdString());
            std::stringstream tbuf; tbuf << tin.rdbuf();
            tin.close();
            QFile::remove(tmpPath);
            std::string tcontent = tbuf.str();
            const size_t tnl = tcontent.find('\n');
            if (tnl == std::string::npos)
                throw std::runtime_error("markings CSV writer produced malformed content (no newline)");
            const std::string tHeader = tcontent.substr(0, tnl);
            const std::string tBody = tcontent.substr(tnl);
            std::ofstream out(ecgSidecar.toStdString(), std::ios::trunc);
            if (!out) throw std::runtime_error("cannot open for write: " + ecgSidecar.toStdString());
            out << suffixValueColumns(tHeader, tsuffix.toStdString()) << tBody;
            out.close();
        }
        std::cout << "Saved sidecar: " << ecgSidecar.toStdString() << "\n";
        std::cout.flush();

        // Final pass: write the pulse sidecar (un-suffixed, once), then MERGE
        // every per-anchor sidecar (R first, then the anchor sequence, then
        // pulse) into the canonical markings CSV in one pass; delete sidecars.
        if (finalPass) {
            const QString pulseSidecar = csvDir.absolutePath() + "/"
                + m_subjectId + "_template_markings_PULSE.csv";
            {
                const QString pulseTmp = pulseSidecar + ".tmp";
                writeTemplateMarkingsCsv(pulseTmp.toStdString(), m_bins,
                    m_subjectId.toStdString(), m_sampleRate, currentAnchor(),
                    MarkingsCsvSection::PulseOnly);
                std::ifstream pin(pulseTmp.toStdString());
                std::stringstream pbuf; pbuf << pin.rdbuf();
                pin.close();
                QFile::remove(pulseTmp);
                std::ofstream out(pulseSidecar.toStdString(), std::ios::trunc);
                if (!out) throw std::runtime_error("cannot open for write: " + pulseSidecar.toStdString());
                out << pbuf.str();   // pulse columns stay un-suffixed
                out.close();
            }

            // Discover the per-anchor sidecars by glob (no dependency on the
            // controller's anchor list). Order R first, then the rest
            // alphabetically, then pulse last.
            const QString stem = m_subjectId + "_template_markings_";
            QDir scDir(csvDir.absolutePath());
            QStringList found = scDir.entryList(
                QStringList{ stem + "*.csv" }, QDir::Files, QDir::Name);
            std::vector<std::string> order;
            // R sidecar first if present.
            const QString rName = stem + "R.csv";
            if (found.contains(rName))
                order.push_back(scDir.filePath(rName).toStdString());
            for (const QString& fn : found) {
                if (fn == rName) continue;
                if (fn == QFileInfo(pulseSidecar).fileName()) continue;  // pulse added last
                order.push_back(scDir.filePath(fn).toStdString());
            }
            order.push_back(pulseSidecar.toStdString());

            if (!mergeSidecarCsvs(csvPath.toStdString(), order))
                throw std::runtime_error("could not merge markings sidecars into " + csvPath.toStdString());
            for (const auto& p : order) QFile::remove(QString::fromStdString(p));
            std::cout << "Merged markings CSV: " << csvPath.toStdString() << "\n";

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

    // Aligned-template CSV: writes a per-anchor sidecar every pass; merge the
    // sidecars into the canonical <id>_template.csv once, at the final pass.
    writeAlignedTemplateCsv();
    if (finalPass) {
        QDir alignedDir(QDir(m_templateDir).absoluteFilePath("../csv_for_analysis"));
        const QString canonical = alignedDir.filePath(m_subjectId + "_template.csv");
        const QString stem = m_subjectId + "_template_";
        QStringList found = alignedDir.entryList(
            QStringList{ stem + "*.csv" }, QDir::Files, QDir::Name);
        // Exclude the canonical itself (<id>_template.csv doesn't match the
        // trailing-underscore stem, but guard anyway) and order R first.
        std::vector<std::string> order;
        const QString rName = stem + "R.csv";
        if (found.contains(rName))
            order.push_back(alignedDir.filePath(rName).toStdString());
        for (const QString& fn : found) {
            if (fn == rName) continue;
            if (fn.startsWith(m_subjectId + "_template_markings")) continue;  // different family
            order.push_back(alignedDir.filePath(fn).toStdString());
        }
        if (!order.empty() && mergeSidecarCsvs(canonical.toStdString(), order)) {
            for (const auto& p : order) QFile::remove(QString::fromStdString(p));
            std::cout << "Merged aligned-template CSV: " << canonical.toStdString() << "\n";
        }
    }

    // Anchor cycle: emit reload until the last pass, then finish.
    if (!finalPass) {
        emit requestQAlignReload();
        return;
    }

    emit finished();
}