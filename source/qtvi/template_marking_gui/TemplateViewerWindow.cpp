#include "TemplateViewerWindow.hpp"
#include "ui_TemplateViewerWindow.h"
#include "markers_automatic_position.hpp"
#include "NormalizeFeatures.hpp"
#include <QMessageBox>
#include <QColor>
#include <QTimer>
#include <QPixmap>
#include <QDir>
#include <QFileInfo>
#include <QApplication>
#include <QLayout>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>

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
        try {
            computeMarkingsForPage();
            showPage();
        }
        catch (const std::exception& e) {
            QMessageBox::critical(this, "Viewer error",
                QString("Failed to prepare page for %1:\n\n%2")
                .arg(m_subjectId, e.what()));
            emit finished();
            return;
        }
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

// Forward-declare: seedBinMarkers lives in an anonymous namespace lower in
// this file (right above the ECG per-channel seed section). loadSubject
// calls it before that namespace appears, so we need the declaration here.
namespace { void seedBinMarkers(TemplateBin& b); }

void TemplateViewerWindow::loadSubject(const QString& templatePath, const QString& markingPath, const QString& subjectId, double sampleRateHz) {

    m_markingPath = markingPath;
    m_templateDir = QFileInfo(templatePath).absolutePath();
    m_subjectId = subjectId;
    m_sampleRate = sampleRateHz;
    setWindowTitle(QString("Template Marking - %1").arg(subjectId));
    ui->subjectLabel->setText(subjectId);

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
    for (auto& b : m_bins) seedBinMarkers(b);
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

// ========================================================================
// Per-bin marker seeding.
//
// Markers detected: ECG (per channel) P-begin, Q-begin, T-begin, T-end;
// PPG Onset, Peak, Dicrotic notch, 50% recovery, End. Detectors return a
// sample index or -1 ("no marker"). Detected indices are clamped into
// the visible portion of their trace -- for ECG, that's the per-trace
// cutoff computed by BinPlotWidget::computeEcgVisibleN (cuts off before
// the next beat's P-wave).
// ========================================================================

namespace {
    // Seed a pulse's markers from the window BETWEEN THE TWO R PEAKS
    // (one cardiac cycle) in the pulse's own sample space. The chain:
    //   peak  = argmax over [rFirst, rSecond]   (systolic peak between the R's)
    //   foot  = min before the peak (anywhere in the template)
    //   end   = min after  the peak (anywhere in the template)
    //   dicrotic / 50% = detected on the [foot, end] single-cycle span.
    //
    // rFirst/rSecond are the two R peaks in this template's sample space.
    // Under Patch B every pulse channel is R-anchored: rFirst = padSamples,
    // rSecond = padSamples + medianRR. Callers with no R window available
    // (arterial channels when ch1.raw is absent) pass rFirst=rSecond=0 and
    // this falls back to the whole trace. No cycleLen / lead-in assumption.
    inline void seedPulse(const std::vector<double>& v, int rFirst, int rSecond,
        int& onset, int& peak, int& dicrotic, int& peak2, int& end) {
        const int n = static_cast<int>(v.size());
        if (n < 1) return;
        auto cl = [&](int x) { return std::clamp(x, 0, n - 1); };

        int lo = cl(rFirst), hi = cl(rSecond);
        if (hi - lo < 3) { lo = 0; hi = n - 1; }   // no R window => whole trace

        // Peak = argmax between the two R peaks (the ONLY R-bounded landmark).
        // Skip NaN and start from -inf so the first real sample wins.
        if (peak < 0) {
            int p = -1; double pv = -std::numeric_limits<double>::infinity();
            for (int i = lo; i <= hi; ++i)
                if (!std::isnan(v[i]) && v[i] > pv) { pv = v[i]; p = i; }
            if (p < 0) return;   // all NaN in window; can't seed
            peak = p;
        }

        // Foot = lowest point before the peak, ANYWHERE in the template.
        // The foot doesn't have to sit inside the R window -- only the peak
        // does. Skip NaN; start from +inf so a real sample always beats it.
        if (onset < 0) {
            int f = -1; double fv = std::numeric_limits<double>::infinity();
            for (int i = 0; i <= peak; ++i)
                if (!std::isnan(v[i]) && v[i] < fv) { fv = v[i]; f = i; }
            if (f < 0) f = std::max(0, peak - 1);
            onset = f;
        }

        // End = lowest point after the peak, ANYWHERE in the template.
        // Same NaN-skipping rule. If everything after the peak is NaN
        // (short template, peak near the tail), fall back to peak+1.
        if (end < 0) {
            int e = -1; double ev = std::numeric_limits<double>::infinity();
            for (int i = peak + 1; i < n; ++i)
                if (!std::isnan(v[i]) && v[i] < ev) { ev = v[i]; e = i; }
            if (e < 0) e = std::min(n - 1, peak + 1);
            end = e;
        }

        // Dicrotic auto-detected on the real single-cycle [foot, end] span.
        if (dicrotic < 0 && cl(end) > cl(onset) + 2) {
            const int base = cl(onset);
            std::vector<double> cyc(v.begin() + base, v.begin() + cl(end) + 1);
            dicrotic = cl(base + ecg_markers::detect_ppg_dicrotic(cyc));
        }
        // peak2 is user-placed. Default: 90% of the way from onset to end
        // so the bar starts near the right of the anchor pulse. User drags.
        if (peak2 < 0 && end > onset) {
            peak2 = cl(onset + (9 * (end - onset)) / 10);
        }
    }


    inline int clampToVisible(int idx, int visN) {
        return std::clamp(idx, 0, visN - 1);
    }

    void seedBinMarkers(TemplateBin& b) {

        // ---- PPG --------------------------------------------------------
        if (b.ppgTemplate.empty()) {
            b.ppg_issue = 2;
            b.ppg_onset = -1;
            b.ppg_peak = -1;
            b.ppg_dicrotic = -1;
            b.ppg_peak2 = -1;
            b.ppg_end = -1;
        }
        else if (b.ppg_issue == 1) {
            b.ppg_onset = -1;
            b.ppg_peak = -1;
            b.ppg_dicrotic = -1;
            b.ppg_peak2 = -1;
            b.ppg_end = -1;
        }
        else {
            // Simple, un-clever rules -- no window bounds, no "first half"
            // heuristics. Each landmark is defined by argmin/argmax over
            // the entire post-anchor trace.
            //
            //   foot     = argmin over [rFirst, N)
            //   peak     = argmax over (foot, N)
            //   end      = argmin over (peak, N)
            //   dicrotic = FIRST interior local min over (peak, end)
            //   peak2    = user-placed; default 90% foot->end
            const std::vector<double>& v = b.ppgTemplate;
            const int N = static_cast<int>(v.size());

            // R sits at column padSamples = 2*avg_r_expand in the PPG frame.
            // avg_r_expand is the rate-derived constant (0.15 * ecgRate) --
            // stable across bins, unaffected by per-bin alignment output.
            double rExp = b.ch1.avg_r_expand_raw;
            if (rExp <= 0) rExp = b.ch2.avg_r_expand_raw;
            if (rExp <= 0) rExp = b.ch3.avg_r_expand_raw;
            const int rFirst = (rExp > 0) ? (int)std::llround(2.0 * rExp) : 0;

            // foot
            int foot = -1;
            {
                double best = std::numeric_limits<double>::infinity();
                const int lo = std::clamp(rFirst, 0, N - 1);
                for (int i = lo; i < N; ++i)
                    if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; foot = i; }
            }
            b.ppg_onset = foot;

            // peak
            int peak = -1;
            if (foot >= 0) {
                double best = -std::numeric_limits<double>::infinity();
                for (int i = foot + 1; i < N; ++i)
                    if (!std::isnan(v[i]) && v[i] > best) { best = v[i]; peak = i; }
            }
            b.ppg_peak = peak;

            // end
            int end = -1;
            if (peak >= 0) {
                double best = std::numeric_limits<double>::infinity();
                for (int i = peak + 1; i < N; ++i)
                    if (!std::isnan(v[i]) && v[i] < best) { best = v[i]; end = i; }
            }
            b.ppg_end = end;

            // Dicrotic notch (movable marker): least-negative slope in the
            // MIDDLE 50% of the downslope. Downslope = [firstMinAfterPeak,
            // end]. If the least-negative slope actually goes non-negative
            // in that middle window, a real inflection exists -- place the
            // marker there. Otherwise place it at the center of the down-
            // slope so the user has a movable bar visible for correction.
            int dic = -1;
            if (peak >= 0 && end > peak + 4) {
                int firstMin = -1;
                for (int i = peak + 1; i < end && i < N - 1; ++i) {
                    if (std::isnan(v[i]) || std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                    if (v[i] <= v[i - 1] && v[i] <= v[i + 1]) { firstMin = i; break; }
                }
                const int slopeStart = (firstMin > 0) ? firstMin : peak;
                const int span = end - slopeStart;
                if (span >= 4) {
                    const int mid_lo = slopeStart + span / 4;
                    const int mid_hi = slopeStart + (3 * span) / 4;
                    int bestIdx = -1;
                    double bestSlope = -std::numeric_limits<double>::infinity();
                    for (int i = std::max(1, mid_lo); i <= std::min(N - 2, mid_hi); ++i) {
                        if (std::isnan(v[i - 1]) || std::isnan(v[i + 1])) continue;
                        const double slope = 0.5 * (v[i + 1] - v[i - 1]);
                        if (slope > bestSlope) { bestSlope = slope; bestIdx = i; }
                    }
                    if (bestIdx >= 0 && bestSlope >= 0.0) {
                        // Real notch (slope turned non-negative).
                        dic = bestIdx;
                    }
                    else {
                        // Subtle / absent notch -- default the movable bar
                        // to the middle of the downslope.
                        dic = (slopeStart + end) / 2;
                    }
                }
            }
            b.ppg_dicrotic = dic;

            // peak2: user-placed. Default 90% foot->end.
            if (foot >= 0 && end > foot)
                b.ppg_peak2 = foot + (9 * (end - foot)) / 10;
            else
                b.ppg_peak2 = std::max(0, N - 1);

            fprintf(stderr,
                "[seed] N=%d foot=%d peak=%d end=%d dic=%d peak2=%d\n",
                N, foot, peak, end, dic, b.ppg_peak2);
        }

        // ---- ECG (per channel) -----------------------------------------
        ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
        for (int c = 0; c < 3; ++c) {
            const auto& ecg = chs[c]->ecgTemplate_raw;
            if (ecg.empty()) {
                b.bad_r_ch[c] = true;
                b.p_begin_ch[c] = -1;
                b.q_begin_ch[c] = -1;
                b.s_end_ch[c] = -1;
                b.t_begin_ch[c] = -1;
                b.t_end_ch[c] = -1;
                continue;
            }
            const int rHint = static_cast<int>(std::round(
                chs[c]->avg_r_expand_raw));

            // T-end first — its detector doesn't need visN, and we need T-end
            // to compute the visible cutoff.
            if (b.t_end_ch[c] < 0)
                b.t_end_ch[c] = ecg_markers::detect_t_end(ecg);

            // Now compute the visible cutoff using T-end.
            const int visN = std::max(static_cast<int>(ecg.size()), 2);

            // Clamp T-end into the visible range now that we know visN.
            b.t_end_ch[c] = clampToVisible(b.t_end_ch[c], visN);

            if (b.p_begin_ch[c] < 0)
                b.p_begin_ch[c] = clampToVisible(
                    ecg_markers::detect_p_begin(ecg), visN);
            if (b.q_begin_ch[c] < 0)
                b.q_begin_ch[c] = clampToVisible(
                    ecg_markers::detect_q_begin(ecg), visN);
            if (b.s_end_ch[c] < 0)
                b.s_end_ch[c] = clampToVisible(
                    ecg_markers::detect_s_end(ecg), visN);
            if (b.t_begin_ch[c] < 0)
                b.t_begin_ch[c] = clampToVisible(
                    ecg_markers::detect_t_begin(ecg), visN);
        }

        // ---- Arterial channels (ABP / ART / ART_PULM) -------------------
        // Same 5-marker set as PPG. Under Patch B, arterial channels are
        // sliced the same way as PPG (R-anchored, driven by ch1.raw), so
        // rFirstArt / rSecondArt land at the same columns in the arterial
        // frame as they do in the PPG frame (assumes same rate; see PPG
        // caveat above). issue==2 => absent, issue==1 => flagged bad
        // (markers cleared), else auto-seed if unset.
        double rExpArt = b.ch1.avg_r_expand_raw;
        if (rExpArt <= 0) rExpArt = b.ch2.avg_r_expand_raw;
        if (rExpArt <= 0) rExpArt = b.ch3.avg_r_expand_raw;
        const int rFirstArt = (rExpArt > 0) ? (int)std::llround(2.0 * rExpArt) : 0;
        const int rSecondArt = rFirstArt + ((rExpArt > 0)
            ? (int)std::llround(5.0 * rExpArt) : 0);

        auto seedArterial = [&](const std::vector<double>& trace, uint8_t& issue,
            int& onset, int& peak, int& dicrotic, int& p50, int& end)
            {
                if (trace.empty()) {
                    issue = 2;
                    onset = peak = dicrotic = p50 = end = -1;
                    return;
                }
                if (issue == 1) {
                    onset = peak = dicrotic = p50 = end = -1;
                    return;
                }
                seedPulse(trace, rFirstArt, rSecondArt,
                    onset, peak, dicrotic, p50, end);
            };
        seedArterial(b.abpTemplate, b.abp_issue,
            b.abp_onset, b.abp_peak, b.abp_dicrotic, b.abp_peak2, b.abp_end);
        seedArterial(b.artTemplate, b.art_issue,
            b.art_onset, b.art_peak, b.art_dicrotic, b.art_peak2, b.art_end);
        seedArterial(b.artPulmTemplate, b.art_pulm_issue,
            b.art_pulm_onset, b.art_pulm_peak, b.art_pulm_dicrotic,
            b.art_pulm_peak2, b.art_pulm_end);
    }
}

void TemplateViewerWindow::computeMarkingsForPage() {
    int start = m_currentPage * m_binsPerPage;
    int end = std::min(start + m_binsPerPage, static_cast<int>(m_bins.size()));
    for (int i = start; i < end; ++i) {
        seedBinMarkers(m_bins[i]);
    }
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
    const double absFoot = std::abs(footY);
    std::vector<double> out(raw.size());
    for (size_t i = 0; i < raw.size(); ++i) {
        if (std::isnan(raw[i])) { out[i] = raw[i]; continue; }
        const double localRatio = (raw[i] - footY) / absFoot * 100.0;
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
            // TEMP DIAGNOSTIC: bypass pulse normalization to test whether
            // the raw PPG/arterial templates are correct.
            const std::vector<double>& ppgN = hasPPG ? b.ppgTemplate : empty;
            const std::vector<double>& abpN = b.abpTemplate;
            const std::vector<double>& artN = b.artTemplate;
            const std::vector<double>& artPN = b.artPulmTemplate;
            /*
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
            */

            pw->setData(ppgN, ppgStd, ecgN, ecgStd,
                b.p_begin_ch[c], b.q_begin_ch[c], b.s_end_ch[c],
                b.t_begin_ch[c], b.t_end_ch[c],
                b.ppg_onset, b.ppg_peak,
                b.ppg_dicrotic, b.ppg_peak2, b.ppg_end,
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
    if (m_bins.empty()) { fprintf(stderr, "[capture] skip: no bins\n"); return; }
    if (m_capturedPages.count(m_currentPage)) return;   // already saved this page
    QDir outDir(m_templateDir);
    const int page = m_currentPage;
    QTimer::singleShot(60, this, [this, page, outDir]() {
        if (page >= m_totalPages) { fprintf(stderr, "[capture] skip: page>=total\n"); return; }
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

void TemplateViewerWindow::writeAlignedTemplateCsv() {
    if (m_bins.empty()) return;

    // All analysis CSVs share one folder: <templateDir>/../csv_for_analysis.
    QDir outDir(QDir(m_templateDir).absoluteFilePath("../csv_for_analysis"));
    if (!outDir.exists()) outDir.mkpath(".");
    const QString path = outDir.filePath(m_subjectId + "_templates.csv");

    std::ofstream f(path.toStdString());
    if (!f) { fprintf(stderr, "[tmplcsv] cannot open %s\n", path.toStdString().c_str()); return; }

    const double toMs = (m_sampleRate > 0.0) ? 1000.0 / m_sampleRate : 1.0;
    constexpr double kEcgRAnchorFrac = 2.0;   // == kEcgPrePFrac / display

    // Column order: shared x, then normalized value + own-peak x for each
    // signal. Values are per-subject normalized (ECG by median(|R|+|S|);
    // pulse by median PI via the local-ratio conversion). Matches the
    // on-screen display.
    f << "file_id,bin_num,x_ms,"
        "ch1_norm,ch1_x_peak_ms,ch2_norm,ch2_x_peak_ms,ch3_norm,ch3_x_peak_ms,"
        "ppg_norm,ppg_x_peak_ms,abp_norm,abp_x_peak_ms,"
        "art_norm,art_x_peak_ms,art_pulm_norm,art_pulm_x_peak_ms\n";
    f << std::setprecision(10);

    auto argmax = [](const std::vector<double>& v) -> int {
        int best_i = 0;
        double best_v = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < static_cast<int>(v.size()); ++i) {
            if (!std::isnan(v[i]) && v[i] > best_v) { best_v = v[i]; best_i = i; }
        }
        return best_i;
        };

    for (size_t bi = 0; bi < m_bins.size(); ++bi) {
        const TemplateBin& b = m_bins[bi];

        // Per-bin normalized traces, matching what's drawn on screen.
        // ECG normalized by that channel's Global_Ref_ecg;
        // pulse normalized by ((sample - foot_y) / |foot_y| * 100) / Global_Ref_pulse.
        // If any ref is unusable the helper returns the raw trace unchanged.
        std::vector<double> ch1N = normalizeEcgTrace(b.ch1.ecgTemplate_raw, 0);
        std::vector<double> ch2N = normalizeEcgTrace(b.ch2.ecgTemplate_raw, 1);
        std::vector<double> ch3N = normalizeEcgTrace(b.ch3.ecgTemplate_raw, 2);
        std::vector<double> ppgN = normalizePulseTrace(b.ppgTemplate, b.ppg_onset, 0);
        std::vector<double> abpN = normalizePulseTrace(b.abpTemplate, b.abp_onset, 1);
        std::vector<double> artN = normalizePulseTrace(b.artTemplate, b.art_onset, 2);
        std::vector<double> artPN = normalizePulseTrace(b.artPulmTemplate, b.art_pulm_onset, 3);

        // One signal = its (normalized) template, its shared-axis start offset,
        // and its own peak index.
        struct Sig {
            const std::vector<double>* v;
            double start;
            int    peak;
        };

        const ChannelTemplateData* ch[3] = { &b.ch1, &b.ch2, &b.ch3 };
        const std::vector<double>* chNorm[3] = { &ch1N, &ch2N, &ch3N };

        std::vector<Sig> sigs;
        for (int c = 0; c < 3; ++c) {
            // ECG R column = argmax of that channel's template. Self-locating.
            sigs.push_back({ chNorm[c], 0.0, argmax(*chNorm[c]) });
        }
        sigs.push_back({ &ppgN,  0.0, argmax(ppgN) });
        sigs.push_back({ &abpN,  0.0, argmax(abpN) });
        sigs.push_back({ &artN,  0.0, argmax(artN) });
        sigs.push_back({ &artPN, 0.0, argmax(artPN) });

        // Row span on the shared axis.
        double loD = 1e300, hiD = -1e300;
        for (const Sig& s : sigs) {
            if (s.v->empty()) continue;
            loD = std::min(loD, s.start);
            hiD = std::max(hiD, s.start + static_cast<double>(s.v->size()));
        }
        if (loD > hiD) continue;   // nothing in this bin
        const int loRow = static_cast<int>(std::floor(loD));
        const int hiRow = static_cast<int>(std::ceil(hiD));

        for (int row = loRow; row < hiRow; ++row) {
            // Shared x: raw timeline, 0 at the top of the ECG templates.
            f << m_subjectId.toStdString() << ',' << bi << ','
                << (row * toMs);
            for (const Sig& s : sigs) {
                const int j = row - static_cast<int>(std::llround(s.start));
                f << ',';
                if (j >= 0 && j < static_cast<int>(s.v->size())
                    && !std::isnan((*s.v)[j]))
                    f << (*s.v)[j];
                f << ',';
                if (j >= 0 && j < static_cast<int>(s.v->size()))
                    f << ((j - s.peak) * toMs);   // 0 at this signal's peak
            }
            f << '\n';
        }
    }
    fprintf(stderr, "[tmplcsv] wrote %s\n", path.toStdString().c_str());
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
            pw->setMarker(BinPlotWidget::EcgP, b.p_begin_ch[c]);
            pw->setMarker(BinPlotWidget::EcgQBegin, b.q_begin_ch[c]);
            pw->setMarker(BinPlotWidget::EcgSEnd, b.s_end_ch[c]);
            pw->setMarker(BinPlotWidget::EcgTBegin, b.t_begin_ch[c]);
            pw->setMarker(BinPlotWidget::EcgTEnd, b.t_end_ch[c]);
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
        switch (marker) {
        case BinPlotWidget::EcgP:      b.p_begin_ch[leadIdx] = newIdx; break;
        case BinPlotWidget::EcgQBegin: b.q_begin_ch[leadIdx] = newIdx; break;
        case BinPlotWidget::EcgSEnd:   b.s_end_ch[leadIdx] = newIdx; break;
        case BinPlotWidget::EcgTBegin: b.t_begin_ch[leadIdx] = newIdx; break;
        case BinPlotWidget::EcgTEnd:   b.t_end_ch[leadIdx] = newIdx; break;
        }

        if (m_moveSubsequent) {
            for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
                // Skip bins where this channel has no ECG or is flagged bad.
                ChannelTemplateData* chs[3] = {
                    &m_bins[i].ch1, &m_bins[i].ch2, &m_bins[i].ch3
                };
                if (chs[leadIdx]->ecgTemplate_raw.empty()) continue;
                if (m_bins[i].bad_r_ch[leadIdx]) continue;
                if (newIdx >= (int)chs[leadIdx]->ecgTemplate_raw.size()) continue;

                switch (marker) {
                case BinPlotWidget::EcgP:
                    m_bins[i].p_begin_ch[leadIdx] = newIdx; break;
                case BinPlotWidget::EcgQBegin:
                    m_bins[i].q_begin_ch[leadIdx] = newIdx; break;
                case BinPlotWidget::EcgSEnd:
                    m_bins[i].s_end_ch[leadIdx] = newIdx; break;
                case BinPlotWidget::EcgTBegin:
                    m_bins[i].t_begin_ch[leadIdx] = newIdx; break;
                case BinPlotWidget::EcgTEnd:
                    m_bins[i].t_end_ch[leadIdx] = newIdx; break;
                }
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
        switch (marker) {
        case BinPlotWidget::PpgOnset:    b.ppg_onset = newIdx; break;
        case BinPlotWidget::PpgPeak:     b.ppg_peak = newIdx; break;
        case BinPlotWidget::PpgDicrotic: b.ppg_dicrotic = newIdx; break;
        case BinPlotWidget::PpgPeak2:       b.ppg_peak2 = newIdx; break;
        case BinPlotWidget::PpgEnd:      b.ppg_end = newIdx; break;
        }
        refreshBinMarkers(binIdx);

        if (m_moveSubsequent) {
            for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
                if (m_bins[i].ppg_issue != 0) continue;
                if (newIdx >= (int)m_bins[i].ppgTemplate.size()) continue;

                switch (marker) {
                case BinPlotWidget::PpgOnset:    m_bins[i].ppg_onset = newIdx; break;
                case BinPlotWidget::PpgPeak:     m_bins[i].ppg_peak = newIdx; break;
                case BinPlotWidget::PpgDicrotic: m_bins[i].ppg_dicrotic = newIdx; break;
                case BinPlotWidget::PpgPeak2:       m_bins[i].ppg_peak2 = newIdx; break;
                case BinPlotWidget::PpgEnd:      m_bins[i].ppg_end = newIdx; break;
                }
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
        // Which trace + issue this marker's channel uses (for bounds/skip).
        auto channelTrace = [&](TemplateBin& tb, int mk,
            const std::vector<double>*& tr, uint8_t*& iss) {
                if (BinPlotWidget::markerIsAbp(mk)) { tr = &tb.abpTemplate; iss = &tb.abp_issue; }
                else if (BinPlotWidget::markerIsArt(mk)) { tr = &tb.artTemplate; iss = &tb.art_issue; }
                else { tr = &tb.artPulmTemplate; iss = &tb.art_pulm_issue; }
            };

        assign(b, marker, newIdx);
        refreshBinMarkers(binIdx);

        if (m_moveSubsequent) {
            for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
                const std::vector<double>* tr = nullptr; uint8_t* iss = nullptr;
                channelTrace(m_bins[i], marker, tr, iss);
                if (!iss || *iss != 0) continue;           // absent/bad channel
                if (!tr || newIdx >= (int)tr->size()) continue;
                assign(m_bins[i], marker, newIdx);
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
    for (auto& b : m_bins)
        seedBinMarkers(b);

    try {
        QString outPath = m_markingPath + "/" + m_subjectId + "_template_markings.bin";
        writeTemplateMarkingsBin(outPath.toStdString(), m_bins);
        std::cout << "Saved: " << outPath.toStdString() << "\n";

        QString csvPath = m_markingPath + "/" + m_subjectId + "_template_markings.csv";
        writeTemplateMarkingsCsv(csvPath.toStdString(), m_bins,
            m_subjectId.toStdString(), m_sampleRate);
        std::cout << "Saved: " << csvPath.toStdString() << "\n";

        // Per-subject amplitude normalization: ECG divided by median(|R|+|S|),
        // pulse channels (PPG + arterial) divided by median PI. See
        // normalize_features.hpp for the full formulas.
        QString normPath = m_markingPath + "/" + m_subjectId + "_template_markings_normalized.csv";
        normalize_features::writeNormalizedCsv(normPath.toStdString(),
            m_bins, m_subjectId.toStdString(), m_sampleRate);
        std::cout << "Saved: " << normPath.toStdString() << "\n";
        std::cout.flush();
    }
    catch (const std::exception& e) {
        QMessageBox::critical(this, "Save failed",
            QString("Could not write markings for %1:\n\n%2\n\n"
                "If the CSV is open in Excel, close it and try again.")
            .arg(m_subjectId, e.what()));
        return;   // don't emit finished(); let the user retry
    }

    emit finished();
}