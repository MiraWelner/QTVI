// ============================================================================
// TemplateViewerWindow.h
// ============================================================================
#pragma once

#include <QMainWindow>
#include <QScrollArea>
#include <QVBoxLayout>
#include <QLabel>
#include <QMenuBar>
#include <QAction>
#include <QPushButton>
#include <vector>
#include "TemplateBinIO.hpp"
#include "BinPlotWidget.h"

namespace Ui { class TemplateViewerWindow; }

class TemplateViewerWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit TemplateViewerWindow(QWidget* parent = nullptr);
    ~TemplateViewerWindow();

    void loadSubject(const QString& templatePath, const QString& markingPath,
        const QString& subjectId);

signals:
    void finished();

private slots:
    void onDicroticMoved(int binIdx, int newIdx);
    void onDicroticDragStarted(int binIdx);
    void onBadRToggled(int binIdx, bool bad);
    void onBadPPGToggled(int binIdx, bool bad);
    void onFinish();
    void onMoveSubsequentToggled(bool on);

private:
    void buildPlots();
    void computeDefaultMarkings();
    void clearPlots();

    Ui::TemplateViewerWindow* ui;

    // Data
    std::vector<TemplateBin> m_bins;
    QString m_markingPath;
    QString m_subjectId;
    std::vector<BinPlotWidget*> m_plots;

    // Mode
    bool m_moveSubsequent = false;
};