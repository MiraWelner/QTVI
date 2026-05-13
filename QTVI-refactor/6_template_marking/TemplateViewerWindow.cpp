#include "TemplateViewerWindow.h"
#include "ui_TemplateViewerWindow.h"
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

    connect(ui->actionMoveSubsequent, &QAction::toggled, this, [this](bool on) {
        m_moveSubsequent = on;
        ui->modeButton->setText(on ? "Mode: Move Subsequent"
            : "Mode: Move Individual");
        });

    addAction(ui->actionMoveSubsequent);
    m_moveSubsequent = ui->actionMoveSubsequent->isChecked();
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
// Per-page setup. No automatic fiducial detection — markers are seeded
// at evenly-spaced positions so the user has something to drag rather
// than a bare trace. Structural facts (no PPG / no ECG for a channel)
// are also recorded.
//
// ECG markers (3) get spaced across the first half of the ECG, which is
// the only portion the viewer draws (the back half is the next beat
// creeping in). PPG markers (2) get spaced across the full PPG length,
// since the viewer shows the whole pulse.
// ========================================================================

namespace {
    // Evenly space `count` markers across [0, visN-1] with equal margins:
    // slot k (0-based) lands at (k+1) / (count+1) of the range.
    inline int spaceSlot(int slot, int count, int visN) {
        if (visN < 2 || count < 1) return 0;
        return static_cast<int>(std::round((slot + 1) * (visN - 1)
            / static_cast<double>(count + 1)));
    }

    // ECG: only the first half of the trace is drawn.
    inline int ecgVisibleN(const std::vector<double>& v) {
        return std::max(static_cast<int>(v.size()) / 2, 2);
    }

    // PPG: the full trace is drawn.
    inline int ppgVisibleN(const std::vector<double>& v) {
        return std::max(static_cast<int>(v.size()), 2);
    }
}

void TemplateViewerWindow::computeMarkingsForPage() {
    int start = m_currentPage * m_binsPerPage;
    int end = std::min(start + m_binsPerPage, static_cast<int>(m_bins.size()));

    for (int i = start; i < end; ++i) {
        TemplateBin& b = m_bins[i];

        // ---- PPG (2 markers across full length) ----
        if (b.ppgTemplate.empty()) {
            b.ppg_issue = 2;
            b.ppg_onset = -1;
            b.ppg_peak = -1;
        }
        else if (b.ppg_issue == 0) {
            int visN = ppgVisibleN(b.ppgTemplate);
            if (b.ppg_onset < 0) b.ppg_onset = spaceSlot(0, 2, visN);
            if (b.ppg_peak < 0) b.ppg_peak = spaceSlot(1, 2, visN);
        }

        // ---- ECG (3 markers across first half, per channel) ----
        ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
        for (int c = 0; c < 3; ++c) {
            const auto& ecg = chs[c]->ecgTemplate_raw;
            if (ecg.empty()) {
                b.bad_r_ch[c] = true;
                b.q_begin_ch[c] = -1;
                b.t_begin_ch[c] = -1;
                b.t_end_ch[c] = -1;
                continue;
            }
            int visN = ecgVisibleN(ecg);
            if (b.q_begin_ch[c] < 0) b.q_begin_ch[c] = spaceSlot(0, 3, visN);
            if (b.t_begin_ch[c] < 0) b.t_begin_ch[c] = spaceSlot(1, 3, visN);
            if (b.t_end_ch[c] < 0) b.t_end_ch[c] = spaceSlot(2, 3, visN);
        }
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
            pw->setData(ppg, ecg,
                b.q_begin_ch[c], b.t_begin_ch[c], b.t_end_ch[c],
                b.ppg_onset, b.ppg_peak);
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

    updatePageControls();
}

// ========================================================================
// Marker movement
// ========================================================================

void TemplateViewerWindow::refreshBinMarkers(int binIdx) {
    // Find this bin's row in the current page (if visible) and push the
    // bin's current markings into every plot showing it.
    for (int li = 0; li < (int)m_pageGlobalIdx.size(); ++li) {
        if (m_pageGlobalIdx[li] != binIdx) continue;
        const TemplateBin& b = m_bins[binIdx];
        for (auto* pw : m_binPlots[li]) {
            int c = pw->leadIndex();
            pw->setMarker(BinPlotWidget::EcgQBegin, b.q_begin_ch[c]);
            pw->setMarker(BinPlotWidget::EcgTBegin, b.t_begin_ch[c]);
            pw->setMarker(BinPlotWidget::EcgTEnd, b.t_end_ch[c]);
            pw->setMarker(BinPlotWidget::PpgOnset, b.ppg_onset);
            pw->setMarker(BinPlotWidget::PpgPeak, b.ppg_peak);
        }
        break;
    }
}

void TemplateViewerWindow::onMarkerMoved(int binIdx, int leadIdx,
    int marker, int newIdx)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    TemplateBin& b = m_bins[binIdx];

    const bool isPpg = BinPlotWidget::markerIsPpg(marker);

    if (BinPlotWidget::markerIsEcg(marker)) {
        if (leadIdx < 0 || leadIdx > 2) return;
        switch (marker) {
        case BinPlotWidget::EcgQBegin: b.q_begin_ch[leadIdx] = newIdx; break;
        case BinPlotWidget::EcgTBegin: b.t_begin_ch[leadIdx] = newIdx; break;
        case BinPlotWidget::EcgTEnd:   b.t_end_ch[leadIdx] = newIdx; break;
        }
        // ECG markings are channel-local; the originating widget has already
        // updated itself. Nothing to propagate.
        return;
    }

    if (isPpg) {
        switch (marker) {
        case BinPlotWidget::PpgOnset: b.ppg_onset = newIdx; break;
        case BinPlotWidget::PpgPeak:  b.ppg_peak = newIdx; break;
        }
        // PPG is shared across this bin's channel plots — update siblings.
        refreshBinMarkers(binIdx);

        // "Move subsequent" applies to PPG markers only (matches the old
        // dicrotic-propagation behavior).
        if (m_moveSubsequent) {
            for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
                if (m_bins[i].ppg_issue != 0) continue;
                if (newIdx < (int)m_bins[i].ppgTemplate.size()) {
                    if (marker == BinPlotWidget::PpgOnset)
                        m_bins[i].ppg_onset = newIdx;
                    else
                        m_bins[i].ppg_peak = newIdx;
                }
            }
            // Push to any of those bins currently on screen.
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
    // Sweep every bin so unvisited ones get the same evenly-spaced seed
    // positions a visited bin would have received from computeMarkingsForPage.
    ChannelTemplateData* chs[3];

    for (auto& b : m_bins) {
        chs[0] = &b.ch1; chs[1] = &b.ch2; chs[2] = &b.ch3;

        for (int c = 0; c < 3; ++c) {
            const auto& ecg = chs[c]->ecgTemplate_raw;
            if (ecg.empty()) {
                b.bad_r_ch[c] = true;
                b.q_begin_ch[c] = -1;
                b.t_begin_ch[c] = -1;
                b.t_end_ch[c] = -1;
                continue;
            }
            int visN = ecgVisibleN(ecg);
            if (b.q_begin_ch[c] < 0) b.q_begin_ch[c] = spaceSlot(0, 3, visN);
            if (b.t_begin_ch[c] < 0) b.t_begin_ch[c] = spaceSlot(1, 3, visN);
            if (b.t_end_ch[c] < 0) b.t_end_ch[c] = spaceSlot(2, 3, visN);
        }

        if (b.ppgTemplate.empty()) {
            b.ppg_issue = 2;
            b.ppg_onset = -1;
            b.ppg_peak = -1;
        }
        else if (b.ppg_issue == 0) {
            int visN = ppgVisibleN(b.ppgTemplate);
            if (b.ppg_onset < 0) b.ppg_onset = spaceSlot(0, 2, visN);
            if (b.ppg_peak < 0) b.ppg_peak = spaceSlot(1, 2, visN);
        }
        else if (b.ppg_issue == 1) {
            // PPG flagged bad — clear positional markings.
            b.ppg_onset = -1;
            b.ppg_peak = -1;
        }
    }

    QString outPath = m_markingPath + "/" + m_subjectId + "_template_markings.bin";
    writeTemplateMarkingsBin(outPath.toStdString(), m_bins);
    std::cout << "Saved: " << outPath.toStdString() << "\n";
    std::cout.flush();
    emit finished();
}