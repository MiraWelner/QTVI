#include "TemplateViewerWindow.h"
#include "ui_TemplateViewerWindow.h"
#include "DumbDicrotic.hpp"
#include <QDir>
#include <QMessageBox>

TemplateViewerWindow::TemplateViewerWindow(QWidget* parent)
    : QMainWindow(parent), ui(new Ui::TemplateViewerWindow) {
    ui->setupUi(this);

    connect(ui->finishButton, &QPushButton::clicked,
        this, &TemplateViewerWindow::onFinish);
    connect(ui->actionMoveSubsequent, &QAction::toggled,
        this, &TemplateViewerWindow::onMoveSubsequentToggled);
}

TemplateViewerWindow::~TemplateViewerWindow() {
    delete ui;
}

void TemplateViewerWindow::onMoveSubsequentToggled(bool on) {
    m_moveSubsequent = on;
    ui->statusLabel->setText(on ? "Mode: Move Subsequent" : "Mode: Move Individual");
}

void TemplateViewerWindow::loadSubject(const QString& templatePath,
    const QString& markingPath,
    const QString& subjectId) {
    m_markingPath = markingPath;
    m_subjectId = subjectId;
    setWindowTitle(QString("Template Marking — %1").arg(subjectId));
    ui->subjectLabel->setText(subjectId);

    m_bins = readTemplateInfoBin(templatePath.toStdString());
    if (m_bins.empty()) {
        QMessageBox::warning(this, "Error", "No bins loaded from " + templatePath);
        return;
    }

    computeDefaultMarkings();
    buildPlots();

    ui->statusLabel->setText(QString("Loaded %1 bins — Mode: %2")
        .arg(m_bins.size())
        .arg(m_moveSubsequent ? "Move Subsequent" : "Move Individual"));
}

void TemplateViewerWindow::computeDefaultMarkings() {
    for (size_t i = 0; i < m_bins.size(); ++i) {
        TemplateBin& b = m_bins[i];
        if (b.ppgTemplate.empty()) {
            b.bad_ppg = true;
            b.templateBad = true;
            b.dicrotic = -1;
            b.peak = -1;
            b.end_idx = -1;
            b.onset = 0;
        }
        else {
            b.dicrotic = dumbDicrotic(b.ppgTemplate);
            double maxVal = b.ppgTemplate[0];
            int maxIdx = 0;
            for (int j = 1; j < (int)b.ppgTemplate.size(); ++j) {
                if (b.ppgTemplate[j] > maxVal) {
                    maxVal = b.ppgTemplate[j];
                    maxIdx = j;
                }
            }
            b.peak = maxIdx;
            b.end_idx = (int)b.ppgTemplate.size() - 1;
            b.onset = 0;
            b.bad_ppg = false;
        }

        if (b.ch1.ecgTemplate_raw.empty()) {
            b.bad_r = true;
        }
    }
}

void TemplateViewerWindow::clearPlots() {
    for (BinPlotWidget* p : m_plots) {
        ui->plotLayout->removeWidget(p);
        delete p;
    }
    m_plots.clear();
}

void TemplateViewerWindow::buildPlots() {
    clearPlots();

    for (size_t i = 0; i < m_bins.size(); ++i) {
        BinPlotWidget* pw = new BinPlotWidget(static_cast<int>(i));
        pw->setFixedHeight(100);

        const TemplateBin& b = m_bins[i];
        pw->setData(b.ppgTemplate,
            b.ch1.ecgTemplate_raw,
            b.ch1.alignment_point_raw,
            b.dicrotic);
        pw->setBadR(b.bad_r);
        pw->setBadPPG(b.bad_ppg);

        connect(pw, &BinPlotWidget::dicroticMoved,
            this, &TemplateViewerWindow::onDicroticMoved);
        connect(pw, &BinPlotWidget::dicroticDragStarted,
            this, &TemplateViewerWindow::onDicroticDragStarted);
        connect(pw, &BinPlotWidget::badRToggled,
            this, &TemplateViewerWindow::onBadRToggled);
        connect(pw, &BinPlotWidget::badPPGToggled,
            this, &TemplateViewerWindow::onBadPPGToggled);

        ui->plotLayout->addWidget(pw);
        m_plots.push_back(pw);
    }
    ui->plotLayout->addStretch();
}

void TemplateViewerWindow::onDicroticMoved(int binIdx, int newIdx) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    m_bins[binIdx].dicrotic = newIdx;

    if (m_moveSubsequent) {
        for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
            if (!m_bins[i].bad_ppg && newIdx < (int)m_bins[i].ppgTemplate.size()) {
                m_bins[i].dicrotic = newIdx;
                m_plots[i]->setDicrotic(newIdx);
            }
        }
    }
}

void TemplateViewerWindow::onDicroticDragStarted(int /*binIdx*/) {
}

void TemplateViewerWindow::onBadRToggled(int binIdx, bool bad) {
    if (binIdx >= 0 && binIdx < (int)m_bins.size())
        m_bins[binIdx].bad_r = bad;
}

void TemplateViewerWindow::onBadPPGToggled(int binIdx, bool bad) {
    if (binIdx >= 0 && binIdx < (int)m_bins.size()) {
        m_bins[binIdx].bad_ppg = bad;
        m_bins[binIdx].templateBad = bad;
    }
}

void TemplateViewerWindow::onFinish() {
    for (size_t i = 0; i < m_bins.size(); ++i) {
        TemplateBin& b = m_bins[i];
        if (!b.bad_ppg && !b.ppgTemplate.empty()) {
            b.onset = 0;
            double maxVal = b.ppgTemplate[0];
            int maxIdx = 0;
            for (int j = 1; j < (int)b.ppgTemplate.size(); ++j) {
                if (b.ppgTemplate[j] > maxVal) {
                    maxVal = b.ppgTemplate[j];
                    maxIdx = j;
                }
            }
            b.peak = maxIdx;
            b.end_idx = (int)b.ppgTemplate.size() - 1;
        }
        else {
            b.onset = -1;
            b.peak = -1;
            b.dicrotic = -1;
            b.end_idx = -1;
        }
    }

    QString outPath = m_markingPath + "/" + m_subjectId + "_template_markings.bin";
    writeTemplateMarkingsBin(outPath.toStdString(), m_bins);

    ui->statusLabel->setText(QString("Saved: %1").arg(outPath));
    emit finished();
}