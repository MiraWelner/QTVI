#pragma once

#include <QMainWindow>
#include <vector>
#include <QString>
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
    void onBadRToggled(int binIdx, int leadIdx, bool bad);
    void onBadPPGToggled(int binIdx, bool bad);
    void onFinish();
    void onToggleMode();
    void onNextPage();
    void onPrevPage();

private:
    // A lead that has ECG data
    struct Lead {
        const std::vector<double>* ecg;
        int channelIndex;          // 0, 1, or 2
        QString label;
    };

    std::vector<Lead> leadsForBin(const TemplateBin& b) const;

    void showPage();
    void clearPlots();
    void computeMarkingsForPage();
    void updatePageControls();
    static std::pair<int, int> compactGrid(int n);

    Ui::TemplateViewerWindow* ui;

    // Data
    std::vector<TemplateBin> m_bins;
    QString m_markingPath;
    QString m_subjectId;

    // All plot widgets on current page (flat, for cleanup)
    std::vector<BinPlotWidget*> m_allPlots;
    // Grouped by page-local bin index
    std::vector<std::vector<BinPlotWidget*>> m_binPlots;
    // Global bin index for each page-local slot
    std::vector<int> m_pageGlobalIdx;

    // Layout sizing (computed once on load)
    int m_maxLeads = 1;
    int m_binsPerPage = 16;

    // Pagination
    int m_currentPage = 0;
    int m_totalPages = 1;

    // Mode
    bool m_moveSubsequent = true;
};