#pragma once

/*
* template_viewer.hpp
* @brief the TemplateViewer class is defined here. It contains:
*
*
*
*/

#include <QMainWindow>
#include <QEvent>
#include <QPointer>
#include <QElapsedTimer> 
#include <vector>
#include <utility>
#include <map>
#include <set>
#include <array>
#include <unordered_map>
#include <cmath>
#include <QString>
#include <QMessageBox>
#include <QRadioButton>
#include <QPushButton>
#include <QButtonGroup>
#include <QColor>
#include <QPixmap>
#include <QDir>
#include <QFileInfo>
#include <QFile>
#include <QCheckBox>
#include <QSpinBox> 
#include <QDockWidget>
#include <QVBoxLayout>
#include <QShortcut>
#include <QKeyEvent>
#include <QApplication>
#include <QGuiApplication>
#include <QStatusBar>
#include <QStringList>
#include <QtConcurrent/QtConcurrentMap>
#include <algorithm>
#include <limits>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cstdio>
#include <cassert>

#include "fiducial_marker_finding/template_marking_bin_io.hpp"
#include "fiducial_marker_finding/bin_plot_widget.hpp"
#include "fiducial_marker_finding/focus_panel_widget.hpp"
#include "fiducial_marker_finding/anchor_view.hpp"
#include "fiducial_marker_finding/feature_marks.hpp"
#include "fiducial_marker_finding/alignment.hpp"
#include "fiducial_marker_finding/global_intervals.hpp"
#include "fiducial_marker_finding/global_interval_lines.hpp"
#include "fiducial_marker_finding/vcg_signal_average.hpp"
#include "fiducial_marker_finding/ppg_derivative.hpp"
#include "fiducial_marker_finding/subsample_refine.hpp"
#include "template_generation/ppg_realign.hpp"
#include "fiducial_marker_finding/curve_fit.hpp"
#include "ui_template_viewer.h"

#include "template_generation/normalize_template_amplitude.hpp"
#include "peak_finding/FilterUtils.hpp"
#include "logging/boundary_training_log.hpp"

//forward declares the earliest onset and latest offset of qrs complex
namespace global_intervals { struct GlobalIntervals; }

class QVBoxLayout;
class QRadioButton;


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
    // ---- THE BANKS AS THE OPERATOR LEFT THEM ----------------------------
    //
    // Read by main.cpp after the window closes, to copy the banks back into
    // the TemplateFile before templates.bin is rewritten.
    //
    // WHY IT HAS TO BE COPIED BACK. showPage sets confirmed_by_operator on
    // THESE bins as each panel is built, and marked_invalid_template lands
    // here on a right-click -- but prepareViewerJob wrote templates.bin before
    // the window opened, from a separate TemplateFile that nothing touches
    // afterwards. So the file's `confirmed` column could only ever read
    // "presumed", for every template in every record, however much marking
    // had been done.
    const std::vector<time_bin>& bins() const { return m_bins; }

    void set_vcg_output_dir(const QString& dir) { m_vcgOutputPath = dir; }
    void setNormOutputDir(const QString& dir) { m_normOutputPath = dir; } //write <id>_feature_norm.csv and <id>_cv_check.csv 

    // THE OPERATOR'S PER-CHANNEL "Lead Reversed" ANSWER, from the noise-
    // marking stage via AnalysisJob::ecg{1,2,3}_inverted. MUST BE CALLED
    // BEFORE loadSubject: loadSubject stamps it onto every bin and the
    // seeding pass reads it from there, so a setter called afterwards would
    // leave the whole record detected as upright -- and unlike the empty-path
    // setters above, the default here is silently WRONG rather than visibly
    // absent. Same trap as set_vcg_output_dir, worse consequence.
    void setLeadPolarity(const LeadPolarity& pol) { m_polarity = pol; }

    void setBeats(const template_structs::BeatsFile* beats) { m_beatsInMemory = beats; }

    void loadSubject(const template_structs::TemplateFile& tf,
        const QString& templateDir, const QString& markingPath,
        const QString& subjectId, double sampleRateHz,
        double ppgRateHz = 0.0, double abpRateHz = 0.0,
        double artRateHz = 0.0, double artPulmRateHz = 0.0,
        double notchFilterHz = 0.0);

signals:
    void finished();

public slots:
    // Wired in Designer via <connections>
    void save_bin_and_csv();
    // Builds every pulse variant that the data allows, so the CSV's _F / _P
    // blocks do not depend on which pages were paged through.
    void buildAllPulseVariants();
    void onNextPage();
    void onPrevPage();

private slots:
    void onMarkerMovedOnTemplate(int binIdx, int leadIdx, int templateIdx, int marker, int newIdx);
    void onMarkerReleasedOnTemplate(int binIdx, int leadIdx, int templateIdx, int marker, int newIdx);
    void movePpgMarker(int binIdx, int leadIdx, int templateIdx, int marker, int newIdx);
    void moveEcgMarker(int binIdx, int leadIdx, int templateIdx, int marker, int newIdx);
    void onMarkerMoved(int binIdx, int leadIdx, int marker, int newIdx);
    void resetMarks();
    void onMarkerDragStarted(int binIdx, int leadIdx, int marker);
    void onBadRToggled(int binIdx, int leadIdx, int templateIdx, bool bad);
    // Helpers for the two above; declared here so both can find them.
    tbank::template_of_all_signals* slotFor(int binIdx, int leadIdx, int templateIdx);
    // One panel's combined bad-ECG / bad-PPG state, from both flag sources.
    BinPlotWidget::State panelState(int binIdx, int leadIdx,
        int templateIdx) const;
    void repaintPanel(int binIdx, int leadIdx, int templateIdx, BinPlotWidget::State st);
    void onBadPPGToggled(int binIdx, int templateIdx, bool bad);
    // `col` IS A DOUBLE, matching BinPlotWidget::landmarkSelected. Qt connects
    // a signal to a slot by parameter type; a mismatch here connects at runtime
    // and then silently never fires, so the two must change together.
    void user_clicked_on_bar(int binIdx, int leadIdx, int templateIdx, int marker, double col); //focus mode - the focus is open in sidebar

private:
    struct Lead {
        const std::vector<double>* ecg;
        const std::vector<double>* ecg_std = nullptr;
        int channelIndex;
        QString label;
        int nMembers = 0;
    };

    std::vector<Lead> leadsForBin(const time_bin& b) const;
    std::vector<Lead> leadsForBinTemplate(const time_bin& b, int templateIdx) const;
    std::vector<int> markingSlotsForBin(const time_bin& b) const;
    std::vector<std::pair<int, int>> pageColumns(int start, int count) const;
    int pageGridRows(bool compact, const std::vector<std::pair<int, int>>& cols) const;
    void addVcgPanel(int gi, int column, int gridRows, const time_bin& b, const std::vector<double>& vcgTrace, double vcgRCol, const global_intervals::GlobalIntervals& intervals, std::vector<BinPlotWidget*>& group, int& usedRows, int& usedCols);
    bool unionEcgFrameSeconds(const time_bin& b, int lead, int templateIdx, double& tMinSec, double& tMaxSec) const;
    void onClassConfirmRequested(int binIndex, int leadIndex, int templateIdx, int annotationCode);
    void applyTemplateToWidget(BinPlotWidget* pw, time_bin& b, int channel, int templateIdx);
    void initAfterBinsLoaded();

    std::vector<double> maybeNotchTrace(const std::vector<double>& sig, double fs, double footIdx) const;
    bool pulseTraceForSlot(tbank::template_of_all_signals& slot, std::vector<double>& outTrace, std::vector<double>& out_std, double& outFootIdx);
    void showPage();
    QString beatsBinPath() const;
    std::set<int> m_ppgRealigned;
    std::map<int, std::pair<std::vector<double>, std::vector<double>>> m_ppgBuilt;
    tbank::PulseAnchor m_ppgViewVariant = tbank::PulseAnchor::Foot;

    // Which variant a pulse landmark's value lives in. ONE DEFINITION, so the
    // drag path, the compose step and the CSV cannot disagree:
    //
    //   PpgOnset                    -> Foot   (it defines that alignment)
    //   PpgDicrotic, PpgEnd         -> Peak   (both sit on the decay)
    //   everything else             -> Peak   (peak, peak2, t50, t80: two of
    //                                          the three brackets are Peak's)
    static tbank::PulseAnchor pulseVariantForMarker(int marker);

    // Build ONE variant of one slot's pulse from the AS-BUILT member rows, and
    // detect that variant's own auto marks on the result. `pct` is ignored for
    // Peak, which is always 100.
    //
    // Idempotent for fixed inputs, which is the entire point: Foot's inputs are
    // the anchor column and the percentage, Peak's are neither. Nothing here
    // reads slot.tmpl, so no build is ever fed a previous build's output.
    bool buildPulseVariant(int binIdx, int templateIdx,
        tbank::PulseAnchor v, double pct, bool announce = false);

    // Install a built variant as the waveform the panels draw, compose
    // pulse_marks from both variants, and push. Builds on demand if the
    // variant is not there yet.
    bool showPulseVariant(int binIdx, int templateIdx, tbank::PulseAnchor v,
        bool announce = false);

    // pulse_marks <- Foot's foot + Peak's notch/end/glyphs. The composed set is
    // what every existing reader sees; pulse_by_variant is the source.
    // ppgRateHz IS A PARAMETER because this is static: the geometric
    // fallbacks at the end of it need a rate (dicrotic sits 0.12 s past the
    // peak) and a static member has no m_ppgRateHz to read. Both call sites
    // are member functions and pass it.
    static void composePulseMarks(tbank::template_of_all_signals& slot,
        double ppgRateHz);

    // The variant the page draws: the forced positions say it outright, Auto
    // defers to m_ppgViewVariant.
    tbank::PulseAnchor pulseVariantToShow() const;

    // Auto's view follows the last pulse BAR clicked. Returns whether the view
    // changed (and therefore whether the page was re-shown).
    bool followPulseBarSelection(int marker);

    // The average the build produced, before any re-level: the stash when one
    // has been taken, the live template when none has. Every anchor and every
    // variant is measured against THIS, so no computation is ever fed its own
    // output.
    const std::vector<double>& ppgAsBuiltTmpl(int binIdx, int templateIdx) const;

    void stashBuiltPulse(int binIdx, int templateIdx,
        const tbank::template_of_all_signals& slot);
    bool restorePulseAsBuilt(int binIdx, int templateIdx);

    // Normalize one slot's pulse through pulseTraceForSlot and push it into
    // every panel of its column, in place. Shared by the re-stack and the
    // restore, so the two cannot come to normalize against different feet.
    void pushPulseToPanels(int binIdx, int templateIdx,
        bool alsoFocus = true);

    // ---- THE FOOT BAR'S OWN GESTURE, ON THE OTHER AXIS ------------------
    //
    // Level every member row to a common baseline at the column the operator
    // dragged the foot bar to, then re-median. NO sample moves sideways: the
    // foot is the pulse's vertical reference (normalize_ppg_or_similar
    // subtracts and divides by it), so "re-do the foot alignment" is a
    // levelling -- the per-beat-crossing reference belongs to the alignment
    // group instead (whose .ui label still reads "Align PPG Horizontal").
    //
    // Returns whether the waveform was replaced; `announce` as
    // relevelPulseAtPct.
    bool relevelPulseAtFoot(int binIdx, int templateIdx, double footCol, bool announce = true);

    const template_structs::BeatsFile* m_beatsInMemory = nullptr;
    // KEYED ON (bin, channel), not bin alone: the ECG re-stack reads "CH1".."CH3"
    // out of the same BeatsFile the pulse path reads "PPG" from, and a cache
    // ignoring the channel would serve lead 0's rows for lead 2.
    int m_beatsCacheBin = -1;
    std::string m_beatsCacheChan;
    ppg_realign::BinBeats m_beatsCache;
    const ppg_realign::BinBeats& beatsForBin(int binIdx,
        const char* channel = "PPG");
    void clearBeatsCache();

    void clearPlots();
    void captureCurrentPage();
    std::string buildAlignedTemplateCsv(AnchorType anchor);
    // Companion file <id>_landmark_fits.csv, beside <id>_bins.csv: per bin/
    // channel/landmark, the fitted curve TYPE and PARAMS (peaks: the weighted
    // quadratic from symmetricExtremumFit; onsets/offsets: the curve_fit model
    // selected by BIC). Recomputed at save, same window convention as the
    // boundary log.
    void writeLandmarkFitsCsv(const std::string& dir);
    std::array<std::vector<std::array<FeatureMarks::TemplateLandmarks, 3>>, 4>  m_exportLm;

    static int anchorSlot(AnchorType a);
    void primeExportLandmarks();
    const FeatureMarks::TemplateLandmarks&
        exportLandmarks(std::size_t bi, int lead, AnchorType a) const;
    void updatePageControls();
    static std::pair<int, int> compactGrid(int n);
    void applyBinCommonToWidget(BinPlotWidget* pw, const time_bin& b);
    void refreshBinMarkers(int binIdx);

    // Bank-column counterpart of refreshBinMarkers. Repaints only the columns
    // showing (binIdx, templateIdx), from that template's own BankMarkerSet.
    // templateIdx 0 forwards to refreshBinMarkers, since slot 0 is the one
    // column that does carry the bin's marker set.
    void refreshBankMarkers(int binIdx, int templateIdx);
    void refreshFocus(BinPlotWidget* pw, int binIdx, int leadIdx, int templateIdx, int marker, double col);
    void pageIn();   // showPage, re-seeding this page with active fit modes first

    FocusPanelWidget* zoomed_in_section_top = nullptr; //for most close ups, they only use focus top
    FocusPanelWidget* zoomed_in_section_bottom = nullptr;   // J point only - the bottom panel is used to show the JT segment (top is QRS)
    bool m_forceAlign = false;
    AnchorType m_forcedAlign = AnchorType::R_PEAK;
    void wireAlignButtons();

    // ---- ONE PLACE THAT CHANGES THE ALIGNMENT ---------------------------
    //
    // The radio buttons, the letter hotkeys and the Tab cycle all land here.
    // Each of them used to carry its own copy of "set the two members, call
    // showPage, then refreshFocus if a landmark is selected", which is three
    // chances for one of them to forget the focus refresh.
    void applyAlignmentSelection(bool force, AnchorType a);

    enum class PpgAlign { Auto, Foot, Percent, Peak };
    PpgAlign m_ppgAlignMode = PpgAlign::Auto;
    int      m_ppgAlignPercent = 0;


    //sometimes the foot is noisy, and so rather than vertical align around foot, it vertical aligns a bit ahead.
    static constexpr double max_foot_std_to_trigger_looking_ahead = 0.01;
    static constexpr int    percent_up_upstroke_to_go_if_noisy = 10;
    static constexpr double region_around_foot_to_measure_std = 0.030;   // +-30 ms about the foot

    //right now either p or q
    int which_alignment_fiducial_marker_should_auto_use(int binIdx, int templateIdx, const tbank::template_of_all_signals& slot, double footCol) const;

    //without this, the percentage box is never read
    double percentage_for_aligning() const {
        switch (m_ppgAlignMode) {
        case PpgAlign::Foot:    return 0.0;
        case PpgAlign::Percent: return static_cast<double>(m_ppgAlignPercent);
        case PpgAlign::Peak:    return 100.0;
        case PpgAlign::Auto:
        default:                return 0.0;
        }
    }

    // ONE PLACE THAT CHANGES THE PULSE ALIGNMENT, as applyAlignmentSelection
    // is for the ECG one: sets the members, syncs the controls, re-stacks the
    // page. The radios and the spin box all land here.
    void applyPpgAlignSelection(PpgAlign mode, int pct);

    // Mode only, NO re-stack. The path a foot DRAG uses to make the group say
    // "Foot": that gesture's re-stack is its own, on the dragged column alone,
    // and a whole-page one from here would re-stack every other column on the
    // page as a side effect of touching one bar.
    void setPpgAlignMode(PpgAlign mode);

    // Re-stack every pulse column on the CURRENT PAGE at the selected
    // alignment. The page and not the record: a record is thousands of columns
    // and re-stacking all of them on a radio click would read every bin's beat
    // matrix to produce waveforms nobody may ever look at.
    void realignAllVisiblePulses();

    void wirePpgAlignButtons();
    // Push the members back into the radio group and the spin box with their
    // signals blocked, so the checked button always names the alignment the
    // panels are actually drawn on.
    void syncPpgAlignControls();

    // The anchor the grid draws in when NOT in a forced alignment (i.e. in
    // Automatic). Set to the last ECG bar the operator clicked, so Automatic
    // follows that bar the same way the focus panel does. Its own member,
    // deliberately: the grid's alignment and the focus zoom's are two
    // separate things and reading one from the other's state is how they end
    // up quietly coupled.
    AnchorType m_autoGridAnchor = AnchorType::R_PEAK;
    QPointer<BinPlotWidget> m_focusWidget;
    curve_fit::FitMode     m_onOffsetFitMode = curve_fit::FitMode::Auto;
    curve_fit::PeakFitMode m_peakFitMode = curve_fit::PeakFitMode::Auto;
    // Last focused landmark, so a fit-mode change re-runs refreshFocus in place.
    int    m_focusBin = -1, m_focusLead = -1, m_focusSlot = -1, m_focusMarker = -1;
    double m_focusCol = -1.0;

    // WHAT THE GRID (and the glyphs, and the R glyph column) IS CURRENTLY
    // DRAWN IN. Forced -> m_forcedAlign; Automatic -> m_autoGridAnchor. ONE
    // place, so leadsForBinTemplate, applyBankTemplateToWidget and
    // applyBinToWidget cannot disagree about what "automatic" means.
    AnchorType currentGridAnchor() const {
        return m_forceAlign ? m_forcedAlign : m_autoGridAnchor;
    }
    void reskinGridForAnchor(int onlyBin = -1, int onlySlot = -1);
    void cycleAlignment(int step); //move to next alignement if user presses tab
    bool eventFilter(QObject* obj, QEvent* ev) override;

    // Last focus selection, so a button press redraws the same landmark
    // instead of waiting for the next click.
    int m_lastFocusBinIdx = -1;
    int m_lastFocusLeadIdx = -1;
    int m_lastFocusTemplateIdx = 0;
    int m_lastFocusMarker = -1;
    double m_lastFocusCol = -1.0;

    QVBoxLayout* m_focusLay = nullptr;
    void setFocusSplit(bool split);
    boundary_training::BoundaryTrainingLog m_boundaryLog;

    std::map<long long, double> m_touchedMarks;
    static long long touchKey(int binIdx, int leadIdx, int marker, AnchorType a) {
        return (((long long)binIdx * 100 + leadIdx) * 100 + marker) * 8 + static_cast<int>(a);
    }

    // Log boundary training data for all landmarks at save (auto_detect from
    // the bin's *_auto_ch glyph fields; expert_mark from the user marks).
    void logBoundaryTrainingAtSave();

    Ui::TemplateViewerWindow* ui;

    std::vector<time_bin> m_bins;
    QString m_markingPath;
    QString m_templateDir;   // folder containing templates.bin/.csv (screenshot target)
    QString m_vcgOutputPath;   // cfg.vcg_output; <id>_vcg.csv lands here
    QString m_normOutputPath;  // feature_norm / cv_check CSVs land here
    // Stamped onto every time_bin in loadSubject; the bins are what every
    // detector call reads, so this member is only the inbound copy.
    LeadPolarity m_polarity;
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

    // (bin, slot) -> page column, rebuilt at the end of showPage. A pair
    // occupies exactly one column, so this is exact, not a first-match cache.
    // Replaces the linear m_pageGlobalIdx scan that the propagation loops ran
    // once per propagated column.
    std::unordered_map<int, int> m_pageColOf;

    // The one key function for a (bin, slot) pair; the stride was a bare 64 in
    // four places. Slots run to max_templates_per_bin * 4, hence the assert.
    static constexpr int kSlotKeyStride = 64;
    static constexpr int slotKey(int bin, int slot) {
        return bin * kSlotKeyStride + slot;
    }
    static_assert(tbank::max_templates_per_bin * 4 <= kSlotKeyStride,
        "slotKey would collide: a bin's slot range no longer fits the stride");

    // This column's panels, or nullptr when the (bin, slot) is not on the page.
    const std::vector<BinPlotWidget*>* panelsForColumn(int binIdx,
        int templateIdx) const;

    // A bin's visible x-axis span in seconds, the denominator every
    // equal-screen-distance propagation divides by. Was an identical lambda in
    // each of the propagation paths. -1 when the bin has no drawable extent.
    double binSpanSeconds(int binIdx) const;

    void clearFocusPanels();
    tbank::BankMarkerSet barsForPanel(const BinPlotWidget* pw, const time_bin& b, int lead, int slot) const;

    struct SdMsModel {
        std::vector<double>  sdMs;        // per column, NaN where floored
        std::vector<uint8_t> floorMask;   // 1 where the floor engaged (shaded)
        std::vector<double>  absSlope;    // per column |dV/dt|, amp/sample
        double               floor = 0.0;
    };
    struct SdKey {
        int bin = -1, lead = -1, slot = -1;
        AnchorType anchor = AnchorType::R_PEAK;
        bool operator==(const SdKey& o) const {
            return bin == o.bin && lead == o.lead && slot == o.slot
                && anchor == o.anchor;
        }
    };
    static SdMsModel sd_in_msec(const std::vector<double>& mean,
        const std::vector<double>& sd, double fs);
    SdMsModel m_sdCache;
    SdKey     m_sdCacheKey;
    bool      m_sdCacheValid = false;
    struct FocusKey {
        int bin = -1, lead = -1, slot = -1, marker = -1;
        AnchorType anchor = AnchorType::R_PEAK;
        const void* panel = nullptr;
        curve_fit::FitMode onOffset = curve_fit::FitMode::Auto;
        curve_fit::PeakFitMode peak = curve_fit::PeakFitMode::Auto;
        bool operator==(const FocusKey& o) const {
            return bin == o.bin && lead == o.lead && slot == o.slot
                && marker == o.marker && anchor == o.anchor
                && panel == o.panel && onOffset == o.onOffset
                && peak == o.peak;
        }
    };
    FocusKey m_focusFastKey;
    bool     m_focusFastValid = false;
    int      m_focusFastLen = 0;      // mean.size() the last full pass saw
    bool     m_focusFastSplit = false;// whether the bottom panel is in use
    // Anything that changes what the panels hold without changing the key
    // above calls this. clearFocusPanels does, so every path that blanks the
    // panels is covered.
    void invalidateFocusFastPath() { m_focusFastValid = false; }

    // THE TWO CHANNEL PATHS. refreshFocus is the dispatcher: it records the
    // focus for the re-fire paths and hands off. pw IS THE PANEL, in both --
    // the detector's position for the focused landmark comes from
    // pw->detectedLandmarks() / pw->detectedPulse(), which is the answer the
    // X glyph is drawn at. Both accept nullptr (a re-fire that could not name
    // its panel) and leave the fiducial absent.
    void focusEcg(BinPlotWidget* pw, time_bin& b, int binIdx, int leadIdx,
        int templateIdx, int marker, double col);
    void focusPulse(BinPlotWidget* pw, time_bin& b, int templateIdx,
        int marker, double col);

    int max_leads = 1;

    int m_binsPerPage = 16;

    // Columns a page may hold. Chosen so a panel keeps a usable width at the
    // window sizes this tool is used at; a bin whose own column count exceeds
    // it gets a page to itself and is the only case that still compresses,
    // which is also the case the columnsForBin diagnostic is about (three or
    // more markable templates in one bin means the bank over-segmented).
    int m_maxColsPerPage = 8;

    // EVERY COLUMN IN THE RECORD, in draw order: (global bin index, template
    // slot). Built by buildPages, which is the only thing that writes it.
    // Pages index into this.
    std::vector<std::pair<int, int>> m_columnTable;

    // (first column, column count) per page -- NOT bins. Every page holds a
    // full grid's worth of columns except the last, which is what packing by
    // bin could not do. Rebuilt whenever marking eligibility changes, because
    // confirming a template's class can add or remove a column and therefore
    // move every later page boundary.
    std::vector<std::pair<int, int>> m_pages;
    void buildPages();

    /// Columns a page may hold: the full compact grid when panels wrap, one
    /// column per bin-width otherwise (the leads occupy the rows there). One
    /// function so buildPages and the grid layout cannot disagree.
    int pageColumnBudget() const;

    int m_currentPage = 0;
    int m_totalPages = 1;

    enum class MoveMode { Individual, SubsequentDelta, SubsequentRaw };
    MoveMode m_moveMode = MoveMode::SubsequentDelta;
    std::map<int, int> original_location_of_bar;//helps ensure that the subsequent bars are moved by the same delta as the first bar
    int originFor(int col, int cur) const;
    // The dragged bar's own position when the drag began, set on the first move
    // event and reset by onMarkerDragStarted. The propagated shift is measured
    // from this, not from the previous event, so the total percentage is
    // computed and rounded once per drag instead of once per mouse-move.
    double m_dragStartIdx = -1.0;

    // ---- PROPAGATED COLUMNS REPAINT AT A RATE, NOT AT EVENT RATE --------
    //
    // Move-Subsequent pushes the dragged bar's shift into every column right
    // of it, and a repaint of one panel is four traces, four bands, ~eight
    // bars and a re-bracketing of P and T peak (the bar set changed, so the
    // reactive cache misses). At twelve columns of three leads that was ~36 of
    // those per drag pixel, which is what made the bar lag the cursor while
    // the arithmetic behind it was never the problem.
    //
    // So the propagated bars are STORED on every move (setMarkerQuiet, so a
    // hit test or an export in between reads the new value) and PAINTED on a
    // clock: the dragged panel stays live at event rate, the rest follow at
    // kDragRepaintMs, and the release flushes whatever is outstanding so the
    // page is never left showing a stale column.
    static constexpr int kDragRepaintMs = 40;   // ~25 fps for the followers
    std::set<int> m_dragDirtyCols;              // page column indices

    // ---- WHICH COLUMNS THIS GESTURE ACTUALLY WROTE A BAR INTO -----------
    //
    // Page column indices, and deliberately NOT m_dragDirtyCols. That one is a
    // repaint queue: flushDragRepaints empties it several times during a drag,
    // so by mouse-up it holds only the last interval's columns. This one
    // accumulates for the whole gesture and is cleared by onMarkerDragStarted,
    // because the RELEASE has to re-stack every column the propagation moved --
    // a bar that moved without its waveform following it is exactly the defect
    // the re-stack exists to fix, and it is no less wrong on a propagated
    // column than on the one under the cursor.
    std::set<int> m_dragPropCols;
    QElapsedTimer m_dragPaintClock;
    // force = the end of a gesture: paint now regardless of the clock.
    void flushDragRepaints(bool force);

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

    // From loadSubject: > 0 means on. No companion bool -- there is no
    // checkbox here, so nothing to disagree with it.
    double m_notchFilterHz = 0.0;

    //Global references - earliest QRS onset, latest QRS offset, etc
    double m_pulseGlobalRef[4] = { std::nan(""), std::nan(""), std::nan(""), std::nan("") };
    void compute_global_refs();

    // One bin's load-time landmark seeding. const because it is called
    // concurrently from initAfterBinsLoaded and must not touch this window's
    // state -- it reads the rates and fit modes and writes only into `b`.
    void seedOneBin(time_bin& b) const;


    void writeNormalizationCsvs(); // Writes <id>_cv_check.csv and <id>_feature_norm.csv

    bool restoreMarkersFrom(const QString& markingsBinPath, bool ecg, bool pulse); //attemps to reload markers from previous markings.bin file, returns true if successful, false otherwise
    // Normalize a copy of `raw` according to the rules in normalize_features.hpp.
    // ECG: sample / globalRef.
    // Pulse: 100*(sample - footY) / footY / globalRef.
    // If globalRef or footY is not usable, returns raw unchanged.
    // This time bin's own |R|+|S|, for this slot. Not a record-wide median:
    // slot indices name the same morphology only within a bin.
    double ecgRefFor(const time_bin& b, int ch, int slot) const;
    std::vector<double> normalizeEcgTrace(const std::vector<double>& raw, const time_bin& b, int ch, int slot) const;

    std::vector<double> normalize_ppg_or_similar(const std::vector<double>& raw, double footIdx, int pulseChan) const;

    // Push m_showEcgMarkers / m_showPpgMarkers into every visible plot.
    void applyMarkerVisibility();
};