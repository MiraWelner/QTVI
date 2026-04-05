#include "TemplateViewerWindow.h"
#include "ui_TemplateViewerWindow.h"
#include "DumbDicrotic.hpp"
#include <QShortcut>
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

    connect(ui->finishButton, &QPushButton::clicked, this, &TemplateViewerWindow::onFinish);
    connect(ui->modeButton, &QPushButton::clicked, this, &TemplateViewerWindow::onToggleMode);
    connect(ui->prevButton, &QPushButton::clicked, this, &TemplateViewerWindow::onPrevPage);
    connect(ui->nextButton, &QPushButton::clicked, this, &TemplateViewerWindow::onNextPage);

    connect(ui->actionMoveSubsequent, &QAction::toggled, this, [this](bool on) {
        m_moveSubsequent = on;
        ui->modeButton->setChecked(on);
        ui->modeButton->setText(on ? "Mode: Move Subsequent" : "Mode: Move Individual");
        });

    m_moveSubsequent = ui->modeButton->isChecked();

    auto bind = [this](auto key, auto slot) {
        connect(new QShortcut(QKeySequence(key), this), &QShortcut::activated, this, slot);
        };
    bind(Qt::Key_M, &TemplateViewerWindow::onToggleMode);
    bind(Qt::Key_Return, &TemplateViewerWindow::onFinish);
    bind(Qt::Key_Right, &TemplateViewerWindow::onNextPage);
    bind(Qt::Key_Left, &TemplateViewerWindow::onPrevPage);
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

void TemplateViewerWindow::onToggleMode() {
    m_moveSubsequent = !m_moveSubsequent;
    ui->modeButton->setChecked(m_moveSubsequent);
    ui->actionMoveSubsequent->setChecked(m_moveSubsequent);
    ui->modeButton->setText(m_moveSubsequent ? "Mode: Move Subsequent"
        : "Mode: Move Individual");
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

    m_bins = readTemplateInfoBin(templatePath.toStdString());
    if (m_bins.empty()) {
        QMessageBox::warning(this, "Error", "No bins loaded from " + templatePath);
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
// Lazy computation for current page
// ========================================================================

void TemplateViewerWindow::computeMarkingsForPage() {
    int start = m_currentPage * m_binsPerPage;
    int end = std::min(start + m_binsPerPage, static_cast<int>(m_bins.size()));

    for (int i = start; i < end; ++i) {
        TemplateBin& b = m_bins[i];

        // PPG
        if (b.ppgTemplate.empty()) {
            b.ppg_issue = 2;       // no ppg
            b.dicrotic = -1;
            b.peak = -1;
            b.end_idx = -1;
            b.onset = -1;
        }
        else if (b.dicrotic < 0 && b.ppg_issue == 0) {
            b.dicrotic = dumbDicrotic(b.ppgTemplate);
            double mx = b.ppgTemplate[0]; int mxI = 0;
            for (int j = 1; j < (int)b.ppgTemplate.size(); ++j)
                if (b.ppgTemplate[j] > mx) { mx = b.ppgTemplate[j]; mxI = j; }
            b.peak = mxI;
            b.end_idx = (int)b.ppgTemplate.size() - 1;
            b.onset = 0;
        }

        // Per-channel: auto-mark bad_r if channel is empty
        const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
        for (int c = 0; c < 3; ++c) {
            if (chs[c]->ecgTemplate_raw.empty())
                b.bad_r_ch[c] = true;
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

    int gridCols;
    if (compact) {
        auto [r, c] = compactGrid(count);
        (void)r;
        gridCols = c;
    }
    else {
        gridCols = m_maxLeads;
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

            pw->setData(ppg, ecg, b.dicrotic);
            pw->setHasPPG(hasPPG);

            if (b.ppg_issue == 1)
                pw->setState(BinPlotWidget::State::BadPPG);
            else if (b.bad_r_ch[leads[li].channelIndex])
                pw->setState(BinPlotWidget::State::BadR);

            connect(pw, &BinPlotWidget::dicroticMoved,
                this, &TemplateViewerWindow::onDicroticMoved);
            connect(pw, &BinPlotWidget::dicroticDragStarted,
                this, &TemplateViewerWindow::onDicroticDragStarted);
            connect(pw, &BinPlotWidget::badRToggled,
                this, &TemplateViewerWindow::onBadRToggled);
            connect(pw, &BinPlotWidget::badPPGToggled,
                this, &TemplateViewerWindow::onBadPPGToggled);

            if (compact) {
                int row = i / gridCols;
                int col = i % gridCols;
                ui->plotGrid->addWidget(pw, row, col);
            }
            else {
                int colspan = (li == (int)leads.size() - 1)
                    ? (gridCols - li) : 1;
                ui->plotGrid->addWidget(pw, i, li, 1, colspan);
            }

            m_allPlots.push_back(pw);
            group.push_back(pw);
        }

        m_binPlots[i] = std::move(group);
    }

    updatePageControls();
}

// ========================================================================
// Dicrotic
// ========================================================================

void TemplateViewerWindow::onDicroticMoved(int binIdx, int newIdx) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    m_bins[binIdx].dicrotic = newIdx;

    for (int li = 0; li < (int)m_binPlots.size(); ++li) {
        if (m_pageGlobalIdx[li] == binIdx) {
            for (auto* pw : m_binPlots[li])
                pw->setDicrotic(newIdx);
            break;
        }
    }

    if (m_moveSubsequent) {
        for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
            if (m_bins[i].ppg_issue == 0 && newIdx < (int)m_bins[i].ppgTemplate.size())
                m_bins[i].dicrotic = newIdx;
        }
        for (int li = 0; li < (int)m_binPlots.size(); ++li) {
            int gi = m_pageGlobalIdx[li];
            if (gi > binIdx && m_bins[gi].ppg_issue == 0) {
                for (auto* pw : m_binPlots[li])
                    pw->setDicrotic(newIdx);
            }
        }
    }
}

void TemplateViewerWindow::onDicroticDragStarted(int) {}

// ========================================================================
// BadR — per channel
// ========================================================================

void TemplateViewerWindow::onBadRToggled(int binIdx, int leadIdx, bool bad) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    if (leadIdx < 0 || leadIdx > 2) return;
    m_bins[binIdx].bad_r_ch[leadIdx] = bad;
}

// ========================================================================
// BadPPG — per bin, propagate to siblings
// ========================================================================

void TemplateViewerWindow::onBadPPGToggled(int binIdx, bool bad) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    m_bins[binIdx].ppg_issue = bad ? 1 : 0;

    if (bad) {
        for (int c = 0; c < 3; ++c)
            m_bins[binIdx].bad_r_ch[c] = false;
    }

    for (int li = 0; li < (int)m_binPlots.size(); ++li) {
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
    const ChannelTemplateData* chs[3];

    for (auto& b : m_bins) {
        chs[0] = &b.ch1; chs[1] = &b.ch2; chs[2] = &b.ch3;

        for (int c = 0; c < 3; ++c) {
            if (chs[c]->ecgTemplate_raw.empty())
                b.bad_r_ch[c] = true;
        }

        if (b.ppgTemplate.empty()) {
            b.ppg_issue = 2;
            b.onset = -1;
            b.peak = -1;
            b.dicrotic = -1;
            b.end_idx = -1;
        }
        else if (b.ppg_issue == 0) {
            if (b.dicrotic < 0)
                b.dicrotic = dumbDicrotic(b.ppgTemplate);
            b.onset = 0;
            double mx = b.ppgTemplate[0]; int mxI = 0;
            for (int j = 1; j < (int)b.ppgTemplate.size(); ++j)
                if (b.ppgTemplate[j] > mx) { mx = b.ppgTemplate[j]; mxI = j; }
            b.peak = mxI;
            b.end_idx = (int)b.ppgTemplate.size() - 1;
        }
        else {
            // ppg_issue == 1 (user marked bad)
            b.onset = -1;
            b.peak = -1;
            b.dicrotic = -1;
            b.end_idx = -1;
        }
    }

    QString outPath = m_markingPath + "/" + m_subjectId + "_template_markings.bin";
    writeTemplateMarkingsBin(outPath.toStdString(), m_bins);
    std::cout << "Saved: " << outPath.toStdString() << "\n";
    std::cout.flush();
    emit finished();
}