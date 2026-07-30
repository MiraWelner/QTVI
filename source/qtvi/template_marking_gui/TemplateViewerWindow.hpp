#pragma once

#include <QMainWindow>
#include <vector>
#include <set>
#include <cmath>
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
        const QString& subjectId, double sampleRateHz);

    // Anchor cycle state. m_anchorStep: -1 = R pass; 0..N-1 index into the
    // controller's anchor sequence. m_qAlignPass is kept as a derived flag
    // (true for any non-R pass) so existing pass-dependent code still works.
    void setAnchorStep(int s) { m_anchorStep = s; m_qAlignPass = (s >= 0); }
    void setAnchorLabel(const QString& s);
    // Total passes in the cycle (R + all anchors); controller sets this so the
    // finish handler knows when to stop emitting reloads. Default 2 = old R/Q.
    void setAnchorPassCount(int n) { m_anchorPassCount = n; }
    // Which AnchorType this pass is marking. Controller sets it before
    // loadSubject; all per-anchor marker reads/writes use it. Defaults to R.
    void setCurrentAnchor(AnchorType a) { m_currentAnchor = a; }
    AnchorType currentAnchor() const { return m_currentAnchor; }

signals:
    void finished();
    // Emitted after the FIRST finish (R-aligned pass). The controller should
    // re-run template generation with Q-alignment enabled and reload the
    // viewer; the button then reads "Finish" and the next click emits
    // finished().
    void requestQAlignReload();

public slots:
    // Wired in Designer via <connections>
    void save_bin_and_csv();
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
    // Save a PNG of the CURRENT page (markers are hidden by default) into the
    // folder holding templates.bin/.csv. Called from showPage(), so each page
    // is captured once, when the user first scrolls to it. No flicker: the
    // page is already on screen and painted.
    void captureCurrentPage();
    // Write templates.csv from the viewer at load time, using the SAME
    // per-bin anchoring the widgets use to draw (rPeak, ppgDelay, ppg_onset),
    // so PPG/arterial rows are shifted relative to ECG exactly as displayed.
    // One shared x_ms column (0 at ch1 R) plus a per-signal *_x_peak_ms column
    // (0 at that signal's own peak). Written once per subject.
    void writeAlignedTemplateCsv();
    void updatePageControls();
    static std::pair<int, int> compactGrid(int n);

    // Pushes the current bin's markings into every plot showing it.
    // Used when a PPG marker drags (which propagates across channels).
    void refreshBinMarkers(int binIdx);

    Ui::TemplateViewerWindow* ui;

    std::vector<TemplateBin> m_bins;
    QString m_markingPath;
    QString m_templateDir;   // folder containing templates.bin/.csv (screenshot target)
    std::set<int> m_capturedPages;   // pages already screenshotted this subject
    QString m_subjectId;
    double  m_sampleRate = 0.0;

    std::vector<BinPlotWidget*> m_allPlots;
    std::vector<std::vector<BinPlotWidget*>> m_binPlots;
    std::vector<int> m_pageGlobalIdx;

    int m_maxLeads = 1;
    int m_binsPerPage = 16;

    int m_currentPage = 0;
    int m_totalPages = 1;

    enum class MoveMode { Individual, SubsequentDelta, SubsequentRaw };
    MoveMode m_moveMode = MoveMode::SubsequentDelta;
    // false = first (R-aligned) pass, button reads "Finish and Next";
    // true  = any anchor pass, button reads per remaining-anchor logic.
    // Derived from m_anchorStep via setAnchorStep; kept for pass-dependent code.
    bool m_qAlignPass = false;
    int  m_anchorStep = -1;        // -1 = R pass; 0..N-1 = anchor index
    int  m_anchorPassCount = 2;    // R + anchors; controller overrides
    QString m_anchorLabel = "R";   // shown in the top bar next to the subject id
    AnchorType m_currentAnchor = AnchorType::R_PEAK;   // which anchor's markers this pass edits
    bool m_showEcgMarkers = false;
    bool m_showPpgMarkers = false;
    bool m_showAbpMarkers = false;
    bool m_showArtMarkers = false;
    bool m_showArtPulmMarkers = false;
    bool m_showEcgTrace = true;
    bool m_showPpgTrace = true;
    bool m_showAbpTrace = true;
    bool m_showArtTrace = true;
    bool m_showArtPulmTrace = true;

    // Per-subject global references, computed once at loadSubject after
    // marker seeding. Used to normalize traces both on-screen and in
    // _templates.csv output. NaN => channel not usable / not normalizable
    // (falls back to raw trace).
    double m_ecgGlobalRef[3] = { std::nan(""), std::nan(""), std::nan("") };
    double m_pulseGlobalRef[4] = { std::nan(""), std::nan(""), std::nan(""), std::nan("") };
    void computeGlobalRefs();
    // Copies the persisted per-bin marker fields (everything a user can
    // drag/toggle -- NOT r_peak_ch, which is always auto-derived from that
    // pass's own template r_col) from a previously-saved
    // _template_markings_*.bin into m_bins, when the bin counts match.
    // Returns false (leaving m_bins untouched) if the file is missing,
    // unreadable, or sized differently than the current template set.
    bool restoreMarkersFrom(const QString& markingsBinPath);
    // Normalize a copy of `raw` according to the rules in normalize_features.hpp.
    // ECG: sample / globalRef.
    // Pulse: 100*(sample - footY) / footY / globalRef.
    // If globalRef or footY is not usable, returns raw unchanged.
    std::vector<double> normalizeEcgTrace(const std::vector<double>& raw, int ch) const;
    std::vector<double> normalize_ppg_or_similar(const std::vector<double>& raw,
        int footIdx, int pulseChan) const;

    // Push m_showEcgMarkers / m_showPpgMarkers into every visible plot.
    void applyMarkerVisibility();
};