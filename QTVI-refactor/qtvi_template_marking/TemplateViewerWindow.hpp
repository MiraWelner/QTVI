#pragma once

#include <QMainWindow>
#include <vector>
#include <QString>
#include "TemplateBinIO.hpp"
#include "BinPlotWidget.hpp"

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

public slots:
    // Wired in Designer via <connections>
    void onFinish();
    void onNextPage();
    void onPrevPage();

private slots:
    void onMarkerMoved(int binIdx, int leadIdx, int marker, int newIdx);
    void onMarkerDragStarted(int binIdx, int leadIdx, int marker);
    void onBadRToggled(int binIdx, int leadIdx, bool bad);
    void onBadPPGToggled(int binIdx, bool bad);

private:
    struct Lead {
        const std::vector<double>* ecg;
        int channelIndex;
        QString label;
    };

    std::vector<Lead> leadsForBin(const TemplateBin& b) const;

    void showPage();
    void clearPlots();
    void computeMarkingsForPage();
    void updatePageControls();
    static std::pair<int, int> compactGrid(int n);

    // Pushes the current bin's markings into every plot showing it.
    // Used when a PPG marker drags (which propagates across channels).
    void refreshBinMarkers(int binIdx);

    Ui::TemplateViewerWindow* ui;

    std::vector<TemplateBin> m_bins;
    QString m_markingPath;
    QString m_subjectId;

    std::vector<BinPlotWidget*> m_allPlots;
    std::vector<std::vector<BinPlotWidget*>> m_binPlots;
    std::vector<int> m_pageGlobalIdx;

    int m_maxLeads = 1;
    int m_binsPerPage = 16;

    int m_currentPage = 0;
    int m_totalPages = 1;

    bool m_moveSubsequent = true;
};