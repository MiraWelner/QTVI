#include "TemplateViewerWindow.hpp"
#include "ui_TemplateViewerWindow.h"
#include "markers_automatic_position.hpp"
#include <QMessageBox>
#include <cmath>
#include <algorithm>
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

    // Marker-visibility toggles. The checkboxes start checked in the
    // .ui, so read their initial state into our flags rather than
    // hard-coding true.
    m_showEcgMarkers = ui->show_ecg->isChecked();
    m_showPpgMarkers = ui->show_ppg->isChecked();
    connect(ui->show_ecg, &QCheckBox::toggled, this, [this](bool on) {
        m_showEcgMarkers = on;
        applyMarkerVisibility();
        });
    connect(ui->show_ppg, &QCheckBox::toggled, this, [this](bool on) {
        m_showPpgMarkers = on;
        applyMarkerVisibility();
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
        computeMarkingsForPage();
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

void TemplateViewerWindow::loadSubject(const QString& templatePath,
    const QString& markingPath, const QString& subjectId) {
    m_markingPath = markingPath;
    m_subjectId = subjectId;
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

    computeMarkingsForPage();
    showPage();
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
            b.ppg_50 = -1;
            b.ppg_end = -1;
        }
        else if (b.ppg_issue == 1) {
            b.ppg_onset = -1;
            b.ppg_peak = -1;
            b.ppg_dicrotic = -1;
            b.ppg_50 = -1;
            b.ppg_end = -1;
        }
        else {
            const int visN = BinPlotWidget::visiblePpgCount(
                static_cast<int>(b.ppgTemplate.size()));
            if (b.ppg_onset < 0) {
                int raw = ecg_markers::detect_ppg_onset(b.ppgTemplate);
                b.ppg_onset = clampToVisible(raw, visN);
            }
            if (b.ppg_peak < 0) {
                int raw = ecg_markers::detect_ppg_peak(b.ppgTemplate);
                b.ppg_peak = clampToVisible(raw, visN);
            }
            if (b.ppg_dicrotic < 0) {
                int raw = ecg_markers::detect_ppg_dicrotic(b.ppgTemplate);
                b.ppg_dicrotic = clampToVisible(raw, visN);
            }
            if (b.ppg_50 < 0) {
                int raw = ecg_markers::detect_ppg_50(b.ppgTemplate);
                b.ppg_50 = clampToVisible(raw, visN);
            }
            if (b.ppg_end < 0) {
                int raw = ecg_markers::detect_ppg_end(b.ppgTemplate);
                b.ppg_end = clampToVisible(raw, visN);
            }
        }

        // ---- ECG (per channel) -----------------------------------------
        ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
        for (int c = 0; c < 3; ++c) {
            const auto& ecg = chs[c]->ecgTemplate_raw;
            if (ecg.empty()) {
                b.bad_r_ch[c] = true;
                b.p_begin_ch[c] = -1;
                b.q_begin_ch[c] = -1;
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
            const int visN = BinPlotWidget::computeEcgVisibleN(ecg, b.t_end_ch[c]);

            // Clamp T-end into the visible range now that we know visN.
            b.t_end_ch[c] = clampToVisible(b.t_end_ch[c], visN);

            if (b.p_begin_ch[c] < 0)
                b.p_begin_ch[c] = clampToVisible(
                    ecg_markers::detect_p_begin(ecg), visN);
            if (b.q_begin_ch[c] < 0)
                b.q_begin_ch[c] = clampToVisible(
                    ecg_markers::detect_q_begin(ecg), visN);
            if (b.t_begin_ch[c] < 0)
                b.t_begin_ch[c] = clampToVisible(
                    ecg_markers::detect_t_begin(ecg), visN);
        }
    }
}

void TemplateViewerWindow::computeMarkingsForPage() {
    int start = m_currentPage * m_binsPerPage;
    int end = std::min(start + m_binsPerPage, static_cast<int>(m_bins.size()));
    for (int i = start; i < end; ++i) {
        seedBinMarkers(m_bins[i]);
    }
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

            static const std::vector<double> empty;
            const auto& ecg = leads[li].ecg ? *leads[li].ecg : empty;
            const auto& ppg = hasPPG ? b.ppgTemplate : empty;

            int c = leads[li].channelIndex;
            // avg_r_expand_raw = median(RR/5) = R-peak position within the
            // ECG template. Used to start PPG at the correct pixel offset.
            const double rPeak = (c == 0) ? b.ch1.avg_r_expand_raw
                : (c == 1) ? b.ch2.avg_r_expand_raw
                : b.ch3.avg_r_expand_raw;

            // std vectors for this channel + PPG. Empty if the templater
            // didn't compute them for this bin -- the widget treats empty
            // std as "no band, just the line".
            const auto& ecgStd = (c == 0) ? b.ch1.ecgTemplate_raw_std
                : (c == 1) ? b.ch2.ecgTemplate_raw_std
                : b.ch3.ecgTemplate_raw_std;
            const auto& ppgStd = hasPPG ? b.ppgTemplate_std : empty;

            pw->setData(ppg, ppgStd, ecg, ecgStd,
                b.p_begin_ch[c], b.q_begin_ch[c],
                b.t_begin_ch[c], b.t_end_ch[c],
                b.ppg_onset, b.ppg_peak,
                b.ppg_dicrotic, b.ppg_50, b.ppg_end,
                rPeak);
            pw->setHasPPG(hasPPG);

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
            }
            else {
                int rowspan = (li == (int)leads.size() - 1)
                    ? (gridRows - li) : 1;
                ui->plotGrid->addWidget(pw, li, i, rowspan, 1);
            }

            m_allPlots.push_back(pw);
            group.push_back(pw);
        }

        m_binPlots[i] = std::move(group);
    }

    applyMarkerVisibility();
    updatePageControls();
}

void TemplateViewerWindow::applyMarkerVisibility() {
    for (auto* pw : m_allPlots) {
        pw->setShowEcgMarkers(m_showEcgMarkers);
        pw->setShowPpgMarkers(m_showPpgMarkers);
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
            pw->setMarker(BinPlotWidget::EcgTBegin, b.t_begin_ch[c]);
            pw->setMarker(BinPlotWidget::EcgTEnd, b.t_end_ch[c]);
            pw->setMarker(BinPlotWidget::PpgOnset, b.ppg_onset);
            pw->setMarker(BinPlotWidget::PpgPeak, b.ppg_peak);
            pw->setMarker(BinPlotWidget::PpgDicrotic, b.ppg_dicrotic);
            pw->setMarker(BinPlotWidget::Ppg50, b.ppg_50);
            pw->setMarker(BinPlotWidget::PpgEnd, b.ppg_end);
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
        case BinPlotWidget::Ppg50:       b.ppg_50 = newIdx; break;
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
                case BinPlotWidget::Ppg50:       m_bins[i].ppg_50 = newIdx; break;
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

// ========================================================================
// Finish
// ========================================================================

void TemplateViewerWindow::onFinish() {
    for (auto& b : m_bins) {
        seedBinMarkers(b);
    }

    QString outPath = m_markingPath + "/" + m_subjectId + "_template_markings.bin";
    writeTemplateMarkingsBin(outPath.toStdString(), m_bins);
    std::cout << "Saved: " << outPath.toStdString() << "\n";
    std::cout.flush();
    emit finished();
}