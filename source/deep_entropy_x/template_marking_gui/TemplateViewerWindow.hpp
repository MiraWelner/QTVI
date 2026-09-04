#pragma once

#include <QMainWindow>
#include <vector>
#include <utility>
#include <map>
#include <cmath>
#include <QString>
#include "template_marking_bin_io.hpp"
#include "BinPlotWidget.hpp"
#include "FocusPanelWidget.hpp"
#include "template_anchoring\anchor_view.hpp"
#include "logging/boundary_training_log.hpp"

namespace Ui { class TemplateViewerWindow; }

class TemplateViewerWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit TemplateViewerWindow(QWidget* parent = nullptr);
    ~TemplateViewerWindow();

    // Destination for operator-confirmed boundary training data (from
    // cfg.training_log). Set once before/at loadSubject; the log is
    // (re)constructed here so the directory is created up front.
    void setBoundaryTrainingDir(const QString& dir) {
        m_boundaryLog = boundary_training::BoundaryTrainingLog(dir.toStdString());
    }
    void setVcgOutputDir(const QString& dir) { m_vcgOutputPath = dir; }
    // Where <id>_feature_norm.csv and <id>_cv_check.csv are written when
    // computeGlobalRefs() runs. Set from cfg (bin_archive_path or a
    // dedicated norm dir) at viewer setup, same as setVcgOutputDir. Empty =
    // don't write.
    void setNormOutputDir(const QString& dir) { m_normOutputPath = dir; }

    void loadSubject(const QString& templatePath, const QString& markingPath,
        const QString& subjectId, double sampleRateHz,
        // Per-pulse-channel upsample rates (Hz). Default 0.0 = unknown, in
        // which case BinPlotWidget::rateRatio() falls back to 1.0 -- the
        // historical behavior from when every channel shared one rate.
        double ppgRateHz = 0.0, double abpRateHz = 0.0,
        double artRateHz = 0.0, double artPulmRateHz = 0.0,
        // Notch frequency for the viewer's display-time notch toggle
        // (see the `notch_filter` checkbox). 0 disables the toggle entirely
        // regardless of the checkbox state; typically comes from
        // cfg.notch_filter_hz so the display filter matches whatever was
        // (or would have been) applied at build time.
        int notchFilterHz = 0);

    // ---- THE ANCHOR CYCLE IS GONE ---------------------------------------
    //
    // setAnchorStep / setAnchorPassCount / setCurrentAnchor / currentAnchor and
    // the requestQAlignReload signal are all deleted. They existed so a
    // controller could drive this window four times, once per alignment,
    // regenerating _templates.bin between openings. The window now loads every
    // alignment at once (TemplateBin::anchored) and each bar reads and writes
    // its own (anchor_view::anchorFor), so there is no pass, no step and no
    // reload.
    //
    // CALLERS MUST BE UPDATED. Whatever owns regenerateWithAnchor /
    // prepareViewerJob must call alignTemplatesFromCache once per
    // anchor_view::kAllAnchors entry BEFORE writing _templates.bin, then open
    // this window once. If it keeps the old sequence, the reload connection
    // simply never fires -- but raw_anchors will hold only the one alignment
    // that run generated, chFor falls back to R for the rest, and every bar's
    // close-up shows the R template. That is the symptom to look for.

signals:
    void finished();

public slots:
    // Wired in Designer via <connections>
    void save_bin_and_csv();
    void onNextPage();
    void onPrevPage();

private slots:
    void onMarkerMoved(int binIdx, int leadIdx, int marker, int newIdx);
    void onMarkerDragStarted(int binIdx, int leadIdx, int marker);
    void onBadRToggled(int binIdx, int leadIdx, int templateIdx, bool bad);
    // Helpers for the two above; declared here so both can find them.
    tbank::BankTemplate* slotFor(int binIdx, int leadIdx, int templateIdx);
    void repaintPanel(int binIdx, int leadIdx, int templateIdx,
        BinPlotWidget::State st);
    void onBadPPGToggled(int binIdx, int templateIdx, bool bad);
    // B2 focus mode: operator selected a landmark in some bin/lead.
    void onLandmarkSelected(int binIdx, int leadIdx, int templateIdx,
        int marker, int col);

private:
    struct Lead {
        const std::vector<double>* ecg;
        int channelIndex;
        QString label;
        // Beats in THIS bank member, not in the bin. The panel title reports
        // it, so it has to travel with the lead rather than being looked up
        // from the bin later -- a bin-level count is a different number and
        // showing one where the other is meant is not visibly wrong.
        // 0 means "unknown": the pre-bank fallback path, where only the bin
        // total exists.
        int nMembers = 0;
    };

    std::vector<Lead> leadsForBin(const TemplateBin& b) const;

    // Leads for one BANK MEMBER of a bin. templateIdx 0 is the sinus seed and
    // returns the chN_raw templates exactly as leadsForBin() does, so a bin
    // with no bank behaves identically to before. Higher indices return that
    // channel's bank template, and a channel whose bank is shorter contributes
    // no lead -- channels are allowed to disagree on template count, so the
    // grid is ragged by design.
    std::vector<Lead> leadsForBinTemplate(const TemplateBin& b,
        int templateIdx) const;

    // Number of grid columns a bin needs: 1 + the highest bank index that earns
    // a column on any channel. Both pagination and layout must call this same
    // function or the two disagree about where a bin's columns end.
    // Bank slots in this bin that want landmark marking, sparse and ascending,
    // always starting with slot 0. Replaces columnsForBin(), which returned a
    // count and could therefore only describe a contiguous prefix.
    std::vector<int> markingSlotsForBin(const TemplateBin& b) const;

    // Section 4.6 class confirmation, from BinPlotWidget::classConfirmRequested.
    // Turns one operator click into tbank::propagateLabel() across all three
    // channels' banks, then rebuilds the page so the label, the subtype the
    // bank issued, and the changed marking eligibility all become visible at
    // once. This is the call site Section 4.6 bullets 3 and 4 were written for
    // and which did not previously exist -- propagateLabel() was reachable from
    // nowhere, so no template in any record had ever been confirmed.
    void onClassConfirmRequested(int binIndex, int leadIndex, int templateIdx,
        int annotationCode);

    // Bars for a bank-template column, from that template's own
    // BankMarkerSet (seeded lazily from its own median). Sub-templates had no
    // bars at all before this: applyBinToWidget draws the BIN's marker set,
    // which describes sinus, so it was correctly applied to slot 0 only.
    void applyBankTemplateToWidget(BinPlotWidget* pw, TemplateBin& b,
        int channel, int templateIdx);

    // Slot-aware marker write-back. Routes to the bin's marker set for slot 0
    // and to the bank template's for every other column.
    void onMarkerMovedOnTemplate(int binIdx, int leadIdx, int templateIdx,
        int marker, int newIdx);

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
    // ONE SIDECAR PER ALIGNMENT, called once per anchor_view::kAllAnchors
    // entry inside a single save. Writes that alignment's own averages and
    // names the file from anchor_view::label(anchor).
    void writeAlignedTemplateCsv(AnchorType anchor);
    void updatePageControls();
    static std::pair<int, int> compactGrid(int n);

    // Pushes the current bin's markings into every plot showing it.
    // Used when a PPG marker drags (which propagates across channels).
    // The ONE place a bin's marker positions and autodetect columns are
    // pushed into a plot widget. Both the initial page build (showPage) and
    // every later refresh go through it, so the bars and the glyphs are always
    // written from the same TemplateBin in the same call.
    void applyBinToWidget(BinPlotWidget* pw, const TemplateBin& b);
    void refreshBinMarkers(int binIdx);

    // Bank-column counterpart of refreshBinMarkers. Repaints only the columns
    // showing (binIdx, templateIdx), from that template's own BankMarkerSet.
    // templateIdx 0 forwards to refreshBinMarkers, since slot 0 is the one
    // column that does carry the bin's marker set.
    void refreshBankMarkers(int binIdx, int templateIdx);

    // Focus mode: ONE panel, showing the clicked landmark's OWN alignment
    // magnified around it -- P-aligned under the P-onset bar, Q under Q-onset,
    // R under the J point, T under T-end (anchor_view::anchorFor). Stats come
    // from that alignment's average via TemplateBin::chFor: mean =
    // ecgTemplate_raw, sd = ecg_template_raw_iqr [ddof=1 std], n = n_beats_raw.
    //
    // There were two panels (QRS above, JT below) because the J point had to
    // appear in both, framed to opposite edges: with only one alignment on
    // screen, neither panel alone could show both of its neighbourhoods. That
    // is what the alignment switch replaces.
    void refreshFocus(int binIdx, int leadIdx, int templateIdx,
        int marker, int col);
    FocusPanelWidget* m_focus = nullptr;
    // Remember the last J-point column per (bin,lead) so a QRS-side or
    // JT-side edit can refresh the other view against the same landmark.
    int m_lastFocusBin = -1, m_lastFocusLead = -1;

    // Operator-confirmed boundary training data (Section 9.10). Destination
    // set via setBoundaryTrainingDir (from cfg.training_log). logBoundary is
    // called at the B2 focus-mode confirmation point.
    boundary_training::BoundaryTrainingLog m_boundaryLog;

    // Operator-touched landmark positions, keyed by (binIdx, leadIdx, marker)
    // packed into a single int, value = the bar position at the click. Filled
    // in onLandmarkSelected (focus activation = bar click). logBoundaryTrainingAtSave
    // reads this to fill confirmedIndex; landmarks never clicked stay blank.
    std::map<long long, int> m_touchedMarks;
    static long long touchKey(int binIdx, int leadIdx, int marker) {
        return ((long long)binIdx * 100 + leadIdx) * 100 + marker;
    }

    // Log boundary training data for all landmarks at save (auto_detect from
    // the bin's *_auto_ch glyph fields; expert_mark from the user marks).
    void logBoundaryTrainingAtSave();

    Ui::TemplateViewerWindow* ui;

    std::vector<TemplateBin> m_bins;
    QString m_markingPath;
    QString m_templateDir;   // folder containing templates.bin/.csv (screenshot target)
    QString m_vcgOutputPath;   // cfg.vcg_output; <id>_vcg.csv lands here
    QString m_normOutputPath;  // feature_norm / cv_check CSVs land here
    QString m_subjectId;
    double  m_sampleRate = 0.0;    // ECG rate; also feeds ECG-only feature/ms code below
    double  m_ppgRateHz = 0.0;
    double  m_abpRateHz = 0.0;
    double  m_artRateHz = 0.0;
    double  m_artPulmRateHz = 0.0;

    std::vector<BinPlotWidget*> m_allPlots;
    std::vector<std::vector<BinPlotWidget*>> m_binPlots;
    std::vector<int> m_pageGlobalIdx;
    // Parallel to m_pageGlobalIdx: which bank member each column shows. A
    // column is now a (bin, template) pair, so every lookup that used to key on
    // the bin index alone has to consult both.
    std::vector<int> m_pageTemplateIdx;

    int m_maxLeads = 1;
    // Bins per page is now a CEILING, not the page size. Pages are packed by
    // COLUMN COUNT instead -- see buildPages(). A bin contributes one column per
    // markable template, so a record with several morphologies per bin used to
    // put 20+ panels on a page and squeeze each to a few pixels wide. The panel
    // width is what makes a template readable, so the page count gives way
    // instead: more pages, each with panels the same size as on a clean record.
    int m_binsPerPage = 16;

    // Columns a page may hold. Chosen so a panel keeps a usable width at the
    // window sizes this tool is used at; a bin whose own column count exceeds
    // it gets a page to itself and is the only case that still compresses,
    // which is also the case the columnsForBin diagnostic is about (three or
    // more markable templates in one bin means the bank over-segmented).
    int m_maxColsPerPage = 8;

    // (first bin, bin count) per page, packed by column budget. Rebuilt whenever
    // marking eligibility changes, because confirming a template's class can add
    // or remove a column and therefore move every later page boundary.
    std::vector<std::pair<int, int>> m_pages;
    void buildPages();

    int m_currentPage = 0;
    int m_totalPages = 1;

    enum class MoveMode { Individual, SubsequentDelta, SubsequentRaw };
    MoveMode m_moveMode = MoveMode::SubsequentDelta;
    // (m_qAlignPass / m_anchorStep / m_anchorPassCount / m_anchorLabel /
    //  m_currentAnchor removed with the cycle. No member holds "the current
    //  alignment" any more, deliberately: whichever alignment a read or write
    //  concerns is a property of the MARKER, answered by
    //  anchor_view::anchorFor, and a member shadowing that is exactly how a
    //  drag on one bar used to land in another alignment's set.)
    void setTitleForSubject();
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
    bool m_showPpgDerivMarkers = false;

    // Display-time notch filter toggle. When on, every template array is
    // pushed through notch_filter() at draw time (in showPage) before it
    // reaches normalize_* and the plot widget. Does NOT rebuild templates
    // or touch saved files -- purely a viewing convenience. m_notchHz is
    // seeded from cfg.notch_filter_hz at loadSubject; a value <= 0 also
    // means "disabled" (the checkbox becomes a no-op).
    bool m_notchFilterOn = false;
    int  m_notchFilterHz = 0;

    // Per-subject global references, computed once at loadSubject after
    // marker seeding. Used to normalize traces both on-screen and in
    // _templates.csv output. NaN => channel not usable / not normalizable
    // (falls back to raw trace).
    double m_ecgGlobalRef[3] = { std::nan(""), std::nan(""), std::nan("") };
    double m_pulseGlobalRef[4] = { std::nan(""), std::nan(""), std::nan(""), std::nan("") };
    void computeGlobalRefs();
    // Sections 5.2-5.4: writes <id>_cv_check.csv and <id>_feature_norm.csv
    // into m_normOutputPath. Called by computeGlobalRefs when that path is
    // set. No-op target otherwise.
    void writeNormalizationCsvs();
    // Copies the persisted per-bin marker fields (everything a user can
    // drag/toggle -- NOT r_peak_ch, which is always auto-derived from that
    // pass's own template r_col) from a previously-saved
    // _template_markings_*.bin into m_bins, when the bin counts match.
    // Returns false (leaving m_bins untouched) if the file is missing,
    // unreadable, or sized differently than the current template set.
    bool restoreMarkersFrom(const QString& markingsBinPath, bool ecg, bool pulse);
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