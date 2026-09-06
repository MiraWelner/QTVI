#include <QMessageBox>
#include <QRadioButton>
#include <QColor>
#include <QPixmap>
#include <QDir>
#include <QFileInfo>
#include <cmath>
#include <algorithm>
#include <limits>
#include <fstream>
#include <sstream>
#include <QFile>
#include <QCheckBox>
#include <QDockWidget>
#include <QVBoxLayout>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <cstdio>
#include <QRadioButton>
#include <QStatusBar>
#include <cassert>
#include <QStringList>

#include "TemplateViewerWindow.hpp"
#include "ui_TemplateViewerWindow.h"
#include "feature_marks.hpp"
#include "template_anchoring\anchor_view.hpp"
#include "template_anchoring\anchor_fit.hpp"
#include "alignment.hpp"
#include "global_intervals.hpp" 
#include "global_interval_lines.hpp"
#include "vcg_signal_average.hpp"
#include "template_generation/NormalizeFeatures.hpp"
#include "peak_finding/FilterUtils.hpp"

namespace {
    constexpr int n_template_cols = 3;
    constexpr int n_template_rows = 4;
    constexpr int max_templates_per_page = n_template_cols * n_template_rows;
}

// ========================================================================
// Construction
// ========================================================================

TemplateViewerWindow::TemplateViewerWindow(QWidget* parent)
    : QMainWindow(parent), ui(new Ui::TemplateViewerWindow)
{
    ui->setupUi(this);

    //when you are moving a marker, do you move just that marker, or all subsequent markers by the same delta, or all subsequent markers to the same raw index?
    auto* moveGroup = new QButtonGroup(this);
    moveGroup->addButton(ui->move_individual);
    moveGroup->addButton(ui->move_subsequent_delta);
    moveGroup->addButton(ui->move_subsequent_raw);
    moveGroup->setExclusive(true);

    connect(ui->move_individual, &QRadioButton::toggled, this, [this](bool on) {
        if (on) m_moveMode = MoveMode::Individual;
        });
    connect(ui->move_subsequent_delta, &QRadioButton::toggled, this, [this](bool on) {
        if (on) m_moveMode = MoveMode::SubsequentDelta;
        });
    connect(ui->move_subsequent_raw, &QRadioButton::toggled, this, [this](bool on) {
        if (on) m_moveMode = MoveMode::SubsequentRaw;
        });

    // Sync m_moveMode to whichever button Designer has checked by default.
    if (ui->move_individual->isChecked())      m_moveMode = MoveMode::Individual;
    else if (ui->move_subsequent_delta->isChecked()) m_moveMode = MoveMode::SubsequentDelta;
    else if (ui->move_subsequent_raw->isChecked())   m_moveMode = MoveMode::SubsequentRaw;

    // Startup defaults, enforced in code (independent of the .ui's checked
    // attributes): all markers OFF, all waveform traces ON. setChecked here
    // updates the visible checkboxes so they reflect the actual state; the
    // toggled() connections below keep the flags in sync thereafter.
    ui->show_ecg_markers->setChecked(false);
    ui->show_ppg_markers->setChecked(false);
    ui->show_abp_markers->setChecked(false);
    ui->show_art_markers->setChecked(false);
    ui->show_art_pulm_markers->setChecked(false);
    ui->show_ecg->setChecked(true);
    ui->show_ppg->setChecked(true);
    ui->show_abp->setChecked(true);
    ui->show_art->setChecked(true);
    ui->show_art_pulm->setChecked(true);

    m_showEcgMarkers = false;
    m_showPpgMarkers = false;
    m_showAbpMarkers = false;
    m_showArtMarkers = false;
    m_showArtPulmMarkers = false;
    m_showEcgTrace = true;
    m_showPpgTrace = true;
    m_showAbpTrace = true;
    m_showArtTrace = true;
    m_showArtPulmTrace = true;
    m_showArtPulmMarkers = false;
    m_showPpgDerivMarkers = false;

    connect(ui->show_ecg_markers, &QCheckBox::toggled, this, [this](bool on) {
        m_showEcgMarkers = on; applyMarkerVisibility();
        });
    connect(ui->show_ppg_markers, &QCheckBox::toggled, this, [this](bool on) {
        m_showPpgMarkers = on; applyMarkerVisibility();
        });
    connect(ui->show_ppg_deriv, &QCheckBox::toggled, this, [this](bool on) {
        m_showPpgDerivMarkers = on; applyMarkerVisibility();
        });
    connect(ui->show_abp_markers, &QCheckBox::toggled, this, [this](bool on) {
        m_showAbpMarkers = on; applyMarkerVisibility();
        });
    connect(ui->show_art_markers, &QCheckBox::toggled, this, [this](bool on) {
        m_showArtMarkers = on; applyMarkerVisibility();
        });
    connect(ui->show_art_pulm_markers, &QCheckBox::toggled, this, [this](bool on) {
        m_showArtPulmMarkers = on; applyMarkerVisibility();
        });

    connect(ui->show_ecg, &QCheckBox::toggled, this, [this](bool on) {
        m_showEcgTrace = on; applyMarkerVisibility();
        });
    connect(ui->show_ppg, &QCheckBox::toggled, this, [this](bool on) {
        m_showPpgTrace = on; applyMarkerVisibility();
        });
    connect(ui->show_abp, &QCheckBox::toggled, this, [this](bool on) {
        m_showAbpTrace = on; applyMarkerVisibility();
        });
    connect(ui->show_art, &QCheckBox::toggled, this, [this](bool on) {
        m_showArtTrace = on; applyMarkerVisibility();
        });
    connect(ui->show_art_pulm, &QCheckBox::toggled, this, [this](bool on) {
        m_showArtPulmTrace = on; applyMarkerVisibility();
        });

    // Optional display-time notch filter toggle. Wired defensively via
    // findChild so this compiles/runs even if the .ui doesn't (yet) declare
    // a checkbox named "notch_filter"; when the widget is present, ticking
    // it re-runs showPage() with each template pushed through notch_filter
    // before drawing. When absent, this block is silently a no-op.
    if (auto* notchBox = findChild<QCheckBox*>("notch_filter")) {
        connect(notchBox, &QCheckBox::toggled, this, [this](bool on) {
            m_notchFilterOn = on;
            showPage();   // full page redraw; templates re-filtered on the way in
            });
    }

    // Focus mode: one panel in a right-side dock, showing the clicked
    // landmark's own alignment magnified around it (see refreshFocus).
    // Created in code (not the .ui) so the existing Designer layout is
    // untouched.
    {
        auto* dock = new QDockWidget(QStringLiteral("Focus"), this);
        dock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
        auto* holder = new QWidget(dock);
        auto* vlay = new QVBoxLayout(holder);
        vlay->setContentsMargins(4, 4, 4, 4);
        vlay->setSpacing(6);
        // TWO PANELS AND A SPACER, LAID OUT IN THIRDS. The J point needs both
        // (see refreshFocus); every other landmark uses the top one and the
        // spacer holds the rest, which is what makes a single view occupy the
        // TOP THIRD instead of stretching over the whole dock.
        zoomed_in_section_top = new FocusPanelWidget(holder);
        zoomed_in_section_bottom = new FocusPanelWidget(holder);
        vlay->addWidget(zoomed_in_section_top, 1);
        vlay->addWidget(zoomed_in_section_bottom, 1);
        vlay->addStretch(1);
        m_focusLay = vlay;
        zoomed_in_section_bottom->hide();
        setFocusSplit(false);
        holder->setLayout(vlay);
        dock->setWidget(holder);
        addDockWidget(Qt::RightDockWidgetArea, dock);
    }
    wireAlignButtons();
}

TemplateViewerWindow::~TemplateViewerWindow() { delete ui; }

// The title used to carry the pass's alignment, because the window WAS one
// alignment. It now holds all four at once, so there is nothing to name here:
// the alignment on show is whichever bar was last clicked, and the focus panel
// header says which that is.
void TemplateViewerWindow::setTitleForSubject() {
    if (!m_subjectId.isEmpty())
        setWindowTitle(QString("Template Marking - %1").arg(m_subjectId));
    else
        setWindowTitle(QStringLiteral("Template Marking"));
}

// ========================================================================
// Helpers
// ========================================================================


// Prevents PPG from going over the ECG window when dragged
static int ecgClipLenFor(const TemplateBin& tb) {
    const ChannelTemplateData* chs[3] = { &tb.ch1, &tb.ch2, &tb.ch3 };
    int mn = -1;
    for (const auto* ch : chs) {
        const int l = static_cast<int>(ch->ecgTemplate_raw.size());
        if (l > 0) mn = (mn < 0) ? l : std::min(mn, l);
    }
    return mn;
}


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

std::vector<TemplateViewerWindow::Lead>
TemplateViewerWindow::leadsForBinTemplate(const TemplateBin& b,
    int templateIdx) const {
    std::vector<Lead> out;
    static const char* kNames[3] = { "Ch1", "Ch2", "Ch3" };

    for (int c = 0; c < 3; ++c) {
        const tbank::TemplateBank& bank = b.ecg_bank[c];

        // Slot 0 falls back to the chN_raw template when no bank reached this
        // bin, so a pre-bank file renders exactly as it always did.
        const std::vector<double>* trace = nullptr;
        int nMembers = 0;
        uint8_t labelCode = tbank::kUnlabeled;

        // tooFewBeats agrees with markingSlotsForBin, which already refused
        // this template a column. Repeated here because the two are reached
        // independently -- refreshBankMarkers and the focus path call this
        // function directly -- and a lead supplied for a slot with no column
        // would draw a trace nothing else on the page accounts for.
        // The pulse side of the same gate. markingSlotsForBin has already
        // refused this slot a column; repeated because this function is also
        // reached directly by refreshBankMarkers and the focus path, and a lead
        // supplied for a slot with no column would draw a trace nothing else on
        // the page accounts for.
        // NO PULSE, NO LEAD. UNCONDITIONALLY.
        //
        // This used to be skipped whenever the bin had no pulse cohort ANYWHERE
        // -- the reasoning being that an absent channel should not suppress
        // markable ECG. That exemption is what put ECG-only panels on screen: a
        // bin whose pulse filter produced nothing took the exemption and every
        // one of its ECG columns was drawn with no pulse to mark.
        //
        // A panel is a (bin, template) pair for marking BOTH faces of a group.
        // A group with no pulse is not markable, whether the pulse is missing
        // because the cohort was thin or because there was none to begin with.
        // The threshold is min_beats_template_ppg from config.csv, via
        // tooFewBeats().
        const bool pulseThin =
            (templateIdx >= b.ppg_bank.size())
            || b.ppg_bank.templates[templateIdx].tmpl.empty()
            || b.ppg_bank.templates[templateIdx].memberCount() <= 0
            || b.ppg_bank.templates[templateIdx].tooFewBeats(/*is_ppg=*/true);
        if (!pulseThin
            && templateIdx < bank.size()
            && !bank.templates[templateIdx].tmpl.empty()
            && !bank.templates[templateIdx].tooFewBeats(/*is_ppg=*/false)
            && (templateIdx == 0
                || bank.templates[templateIdx].wantsLandmarkMarking())) {
            const tbank::BankTemplate& t = bank.templates[templateIdx];
            trace = &t.tmpl;
            nMembers = t.memberCount();
            labelCode = t.label_code;
            // subtype is no longer read here: tbank::letterRanks applies the
            // confirmed-subtype rule itself, from the same BankTemplate.
        }
        else if (templateIdx == 0 && !pulseThin) {
            // PRE-BANK FALLBACK, AND IT MUST STILL RESPECT THE PULSE MINIMUM.
            //
            // A file with no bank has no BankTemplate to measure, so slot 0
            // renders from the bin's own chN_raw -- that is what this branch is
            // for. But it was unconditional, so slot 0 also caught every bin
            // whose bank DID exist and whose slot 0 the test above had just
            // refused. pulseThin was computed and then bypassed here, which is
            // why panels with no PPG kept appearing however tightly
            // markingSlotsForBin and showPage were gated: the panel never
            // reached those gates, it took this exit.
            //
            // pulseThin carries min_beats_template_ppg from config.csv via
            // tooFewBeats(), and is false when the bin has no pulse cohort at
            // all -- an ABSENT channel, which must not suppress markable ECG.
            const std::vector<double>* raw =
                (c == 0) ? &b.ch1.ecgTemplate_raw
                : (c == 1) ? &b.ch2.ecgTemplate_raw : &b.ch3.ecgTemplate_raw;
            if (raw->empty()) continue;
            trace = raw;
        }
        else {
            continue;   // ragged: this channel's bank is shorter
        }

        // NAMING. Bin, then the class, then an underscore and the letter the
        // ALGORITHM assigned when it separated the morphologies: A, B, C in
        // bank order. PQRST is the name for a template no operator has
        // confirmed yet -- it says "an ordinary complex", not "this is normal
        // and not ectopic", which is still the operator's call. Once a beat in
        // the template is confirmed, the class replaces PQRST and the letter
        // tracks the subtype index the bank issued (PVC_A, PVC_B).
        QString cls = "PQRST";
        if (labelCode != tbank::kUnlabeled) {
            switch (labelCode) {
            case tbank::kCodePvc:       cls = "PVC";   break;
            case tbank::kCodePac:       cls = "PAC";   break;
            case tbank::kCodeVt:        cls = "VT";    break;
            case tbank::kCodeMinorNoise:cls = "NOISE"; break;
            default: cls = QString("CODE%1").arg(labelCode); break;
            }
        }

        // THE LETTER COMES FROM tbank::letterRanks, the same function the
        // morphology CSV/bin writers use. It used to be the raw bank slot
        // (`letterIdx = templateIdx`), which skipped letters whenever a lower
        // slot was empty -- a bin showing slots 0, 1, 5 read A, B, F while the
        // CSV called that third template C. One template, two names. The
        // shared function also handles the confirmed-subtype case, so the
        // local subtype override is gone with it.
        int letterIdx = 0;
        if (templateIdx < bank.size()) {
            const std::vector<uint8_t> letters = tbank::letterRanks(bank);
            if (static_cast<size_t>(templateIdx) < letters.size())
                letterIdx = letters[static_cast<size_t>(templateIdx)];
        }
        const QChar letter = QChar('A' + (letterIdx % 26));

        // The count is NOT appended to the label any more. It used to read
        // "Ch1 PQRST_A n=365", which put a per-TEMPLATE number inside the
        // name and left the panel's separate beat figure to report the whole
        // bin -- two counts of different things, one line, neither labelled
        // as to which. The name is now just the name; the count travels on
        // Lead::nMembers and the widget prints it as the ECG beat count.
        const QString lbl = QString("%1 %2_%3").arg(kNames[c]).arg(cls).arg(letter);

        out.push_back({ trace, c, lbl, nMembers });
    }
    return out;
}

// Pack bins into pages so no page exceeds m_maxColsPerPage COLUMNS.
//
// WHY BY COLUMNS AND NOT BY BINS. A page used to hold a fixed m_binsPerPage
// bins, and each bin contributes one column per markable template, so the panel
// count was whatever the banks happened to produce -- four bins with five
// markable templates each is twenty panels sharing the window, a few pixels
// apiece and unreadable. Panel width is what makes a template markable at all,
// so it is the page COUNT that gives way: bins per page falls until the columns
// fit, and a clean record still shows m_binsPerPage bins because its bins
// contribute one column each.
//
// A bin whose own column count exceeds the budget gets a page to itself and is
// the only case that still compresses. Splitting one bin across two pages would
// avoid even that, but a bin's panels are edited and snapshotted as a unit
// (captureCurrentPage), so half a bin per page would mean half a snapshot.
//
// Called again whenever eligibility changes -- confirming a class can add or
// remove a column and therefore move every later boundary.
void TemplateViewerWindow::buildPages() {
    m_pages.clear();
    const int nBins = static_cast<int>(m_bins.size());
    if (nBins == 0) { m_pages.push_back({ 0, 0 }); m_totalPages = 1; return; }

    const int budget = std::max(1, m_maxColsPerPage);
    int i = 0;
    while (i < nBins) {
        int cols = 0, taken = 0;
        while (i + taken < nBins && taken < m_binsPerPage) {
            const int c = static_cast<int>(
                markingSlotsForBin(m_bins[i + taken]).size());
            // Always take at least one bin, even if it alone blows the budget:
            // a page must make progress or paging never terminates. A bin
            // contributing ZERO columns (every template below the configured
            // minimum) is taken for the same reason -- skipping it would leave
            // `taken` at 0 and the outer loop would never advance.
            if (taken > 0 && c > 0 && cols + c > budget) break;
            cols += c;
            ++taken;
        }
        m_pages.push_back({ i, taken });
        i += taken;
    }
    m_totalPages = static_cast<int>(m_pages.size());

    int widest = 0;
    for (const auto& pg : m_pages) {
        int cols = 0;
        for (int k = 0; k < pg.second; ++k)
            cols += static_cast<int>(markingSlotsForBin(m_bins[pg.first + k]).size());
        widest = std::max(widest, cols);
    }
    std::fprintf(stderr,
        "  [pages] %d bins -> %d pages (budget %d cols/page, widest %d)\n",
        nBins, m_totalPages, budget, widest);
    std::fflush(stderr);
}

std::vector<int> TemplateViewerWindow::markingSlotsForBin(const TemplateBin& b) const {
    // ONLY CATEGORY 1 TEMPLATES GET A COLUMN. Landmark marking exists to feed
    // feature extraction, and only category 1 beats do that -- a P-onset on a
    // PVC template has nothing downstream to consume it, and a PVC's QT is not
    // comparable to a sinus QT, so putting them in one feature column would be
    // worse than leaving it empty. Ectopic and noise templates take a class
    // label from the operator and are not landmark-marked at all.
    //
    // Category 1 now means what Section 4.6 says it means: slot 0, which is the
    // Phase 1 sinus template by construction, plus any template an operator has
    // CONFIRMED as normal. Nothing is presumed. So a bin presents one markable
    // column until an operator says otherwise, and the six-column pages came
    // from the presumption, not from the bank.
    //
    // RETURNS THE ELIGIBLE SLOTS, NOT A COUNT. This used to return an int, and a
    // count can only describe a contiguous prefix 0..n-1: `cols = max(cols, t+1)`
    // meant one eligible slot 5 created columns for slots 1-4 as well, whatever
    // the per-slot test had said about them. leadsForBinTemplate() then correctly
    // refused to supply a lead for those slots, so they rendered as empty panels
    // -- the gate was being applied one layer too late to affect layout. A sparse
    // list cannot express the wrong thing.
    //
    // TREAT THREE OR MORE AS A DIAGNOSTIC. Being asked to mark three or more
    // normal templates in one bin means the bank has over-segmented a single
    // morphology. Do not raise the match threshold to hide it -- that merges the
    // fragments and removes the symptom while leaving the cause. Check
    // normalization, baseline drift, and false R detections first.
    // NOT NAMED `slots`. Qt's qobjectdefs.h defines `slots` as an empty macro
    // unless QT_NO_KEYWORDS is set, so a local of that name vanishes at
    // preprocess time: the declaration becomes `std::vector<int> { 0 };` and
    // every use becomes a bare `.push_back(...)`. MSVC reports it as
    // "C2059: syntax error: '.'" several lines from the actual cause, plus a
    // spurious "function must return a value" because `return slots;` collapses
    // to `return ;`. `signals`, `emit` and `foreach` are the same trap.
    // ---- AND NOT IF IT HAS TOO FEW BEATS ------------------------------
    //
    // tbank::minBeatsEcg(), from config.csv. A template below it is flagged
    // too_few_beats in _templates.csv and _templates.bin and gets no column
    // here -- the archive keeps it, the operator is not asked to mark it.
    //
    // THIS APPLIES TO SLOT 0 TOO. A seed built from three beats is no more
    // markable than a spawned template built from three, and exempting it
    // would be the same "it is slot 0 so it must be fine" assumption that
    // put 2-beat waveforms in front of the operator in the first place. So a
    // bin can now contribute ZERO columns, and the callers have to tolerate
    // that: buildPages takes a bin with no columns without stalling, and
    // showPage reports which bins vanished rather than letting them disappear
    // quietly.
    // EITHER CHANNEL BEING THIN INVALIDATES THE TEMPLATE. A group is one set
    // of beats seen on four channels, so a slot whose pulse cohort is below
    // min_beats_template_ppg is not a half-good template -- it is a template
    // whose PPG face rests on too few beats to be a reference, and marking its
    // ECG face would attach landmarks to a group that cannot be measured
    // jointly. The two minimums are one gate.
    //
    // A CHANNEL WITH NO DATA IN THIS BIN DOES NOT VOTE. If the pulse filter
    // produced no template for the bin at all, its pulse bank holds no members
    // anywhere -- that is the channel being ABSENT, not thin, and it must not
    // suppress ECG columns that are perfectly markable. The test is therefore
    // conditional on the bin having some pulse cohort to speak of; the same
    // rule an absent CH2 or CH3 already gets by having no bank.
    auto shown = [](const TemplateBin& bb, int t) {
        {
            // UNCONDITIONAL, for the reason given at leadsForBinTemplate's
            // pulseThin: the "no pulse cohort anywhere is an ABSENT channel"
            // exemption that used to wrap this is exactly what let ECG-only
            // columns through.
            if (t >= bb.ppg_bank.size()) return false;
            const tbank::BankTemplate& pt = bb.ppg_bank.templates[t];
            // STRUCTURAL GATES FIRST, THEN THE CONFIGURED ONE. tooFewBeats
            // compares cleanCount() against minBeatsPpg(), and returns FALSE
            // when that minimum is 0 -- which is the default. So on a config
            // without min_beats_template_ppg this branch tested nothing, and a
            // slot with an EMPTY pulse template or ZERO members passed straight
            // through to become an ECG-only panel.
            //
            // These two are not thresholds and are not configurable: a slot with
            // no waveform and a slot with no members are not thin, they are
            // absent, and the same three tests are what showPage's hasPPG
            // applies at the draw site. The two must agree or a column is
            // created that cannot draw its pulse.
            if (pt.tmpl.empty()) return false;
            if (pt.memberCount() <= 0) return false;
            if (pt.tooFewBeats(/*is_ppg=*/true)) return false;
        }
        for (int c = 0; c < 3; ++c) {
            const tbank::TemplateBank& bank = bb.ecg_bank[c];
            if (t >= bank.size()) continue;
            const tbank::BankTemplate& tp = bank.templates[t];
            if (tp.tmpl.empty()) continue;
            if (tp.tooFewBeats(/*is_ppg=*/false)) continue;
            if (t == 0 || tp.wantsLandmarkMarking()) return true;
        }
        return false;
        };

    std::vector<int> eligible;
    // Slot 0 on a bank-less bin has no BankTemplate to measure, so there is
    // nothing to suppress and it keeps its column: the chN_raw template IS
    // the whole bin, and its beat count is the bin's.
    bool anyBank = false;
    for (int c = 0; c < 3; ++c) if (b.ecg_bank[c].size() > 0) anyBank = true;
    if (!anyBank || shown(b, 0)) eligible.push_back(0);

    for (int t = 1; t < tbank::kDefaultMaxTemplatesPerBin * 4; ++t)
        if (shown(b, t)) eligible.push_back(t);
    return eligible;
}

// Rows x columns for n compact panels, inside the kMaxGridRows x kMaxGridCols
// box. Fill widthwise first (a row of 5 before a second row starts), so a page
// with few panels stays one wide row rather than a tall thin stack.
//
// RETURNS ROWS FIRST, AND THE CALLER USES THE ROW COUNT. This used to be a
// square-ish grid that put the LARGER dimension second and showPage() then
// read that second value as the row count -- so a 16-panel page laid out 4x4
// either way and the mislabelling never showed. It shows the moment the box
// stops being square, hence the explicit note.
std::pair<int, int> TemplateViewerWindow::compactGrid(int n) {
    if (n <= 0) return { 1, 1 };
    const int cols = std::min(n_template_cols, n);
    const int rows = std::min(n_template_rows, (n + cols - 1) / cols);
    return { rows, cols };
}

void TemplateViewerWindow::updatePageControls() {
    ui->pageLabel->setText(QString("%1 / %2").arg(m_currentPage + 1).arg(m_totalPages));
    ui->prevButton->setEnabled(m_currentPage > 0);
    ui->nextButton->setEnabled(m_currentPage < m_totalPages - 1);
}

void TemplateViewerWindow::onNextPage() {
    if (m_currentPage < m_totalPages - 1) {
        captureCurrentPage();   // snapshot the page we are leaving
        ++m_currentPage;
        showPage();
    }
}

void TemplateViewerWindow::onPrevPage() {
    if (m_currentPage > 0) {
        captureCurrentPage();   // snapshot the page we are leaving
        --m_currentPage;
        showPage();
    }
}

// ========================================================================
// Load subject
// ========================================================================

void TemplateViewerWindow::loadSubject(const QString& templatePath, const QString& markingPath,
    const QString& subjectId, double sampleRateHz,
    double ppgRateHz, double abpRateHz, double artRateHz, double artPulmRateHz,
    int notchFilterHz) {

    m_markingPath = markingPath;
    m_templateDir = QFileInfo(templatePath).absolutePath();
    m_subjectId = subjectId;
    m_sampleRate = sampleRateHz;
    m_ppgRateHz = ppgRateHz;
    m_abpRateHz = abpRateHz;
    m_artRateHz = artRateHz;
    m_artPulmRateHz = artPulmRateHz;
    m_notchFilterHz = notchFilterHz;
    setTitleForSubject();
    ui->subjectLabel->setText(subjectId);
    // ONE PASS. The button used to read "Finish and Next" for three of four
    // openings, each one regenerating templates and reloading the window on a
    // different alignment.
    ui->finishButton->setText("Finish");

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
    // Compact mode wraps its panels, so a page may hold the whole box
    // (kMaxGridRows x kMaxGridCols panels). Lead-per-row mode spends a whole
    // column on each bin and stacks its leads down the rows, so its budget is
    // the column cap itself.
    m_binsPerPage = (m_maxLeads <= 1) ? max_templates_per_page : n_template_cols;
    m_maxColsPerPage = (m_maxLeads <= 1) ? max_templates_per_page : n_template_cols;

    m_currentPage = 0;
    buildPages();

    // Seed markers for every bin so per-subject global refs (which need
    // R/S and foot/peak positions across all bins) can be computed once
    // and stay stable across paging.
    // ---- FOUR SEEDING PASSES, ONE LOAD --------------------------------
    //
    // seed_all detects on ch1..3, so each alignment's templates take that slot
    // for the length of its own call, and R goes last -- leaving the flat
    // ch1..3 and the flat *_auto_ch fields holding R, which is what the grid
    // draws and what every existing consumer of those fields expects.
    //
    // Every alignment's glyphs are kept (auto_by_anchor) because a glyph is
    // reported on all four; every alignment's bar is seeded into its own
    // markers_by_anchor entry, and the exclusive mask in maskFor means each
    // pass seeds only the one bar it owns.
    //
    // This loop is what four separate windows used to be.
    for (auto& b : m_bins) {
        const std::array<ChannelTemplateData, 3> savedR = { b.ch1, b.ch2, b.ch3 };
        for (AnchorType a : anchor_view::kAllAnchors) {
            if (a == AnchorType::R_PEAK) {
                b.ch1 = savedR[0]; b.ch2 = savedR[1]; b.ch3 = savedR[2];
            }
            else {
                auto it = b.anchored.find(static_cast<int>(a));
                if (it == b.anchored.end()) continue;   // no block -> nothing to seed
                b.ch1 = it->second[0]; b.ch2 = it->second[1]; b.ch3 = it->second[2];
            }
            FeatureMarks::seed_all(b, m_sampleRate, m_ppgRateHz, a);
            // Captures the flat fields seed_all just wrote: autoFor falls back
            // to them when this alignment has no entry yet, which is exactly
            // the capture wanted.
            b.auto_by_anchor[static_cast<int>(a)] = b.autoFor(a);
        }
        // R last so the flat state the grid reads is R's.
        b.ch1 = savedR[0]; b.ch2 = savedR[1]; b.ch3 = savedR[2];
        FeatureMarks::seed_all(b, m_sampleRate, m_ppgRateHz, AnchorType::R_PEAK);
        for (int c = 0; c < 3; ++c) b.syncReactiveGlyphs(c, 0, m_sampleRate);
    }

    // If this subject was already marked in a previous session, restore
    // those marker positions from the single canonical marking file. If
    // it doesn't exist (never marked before), the fresh auto-seed above
    // stands unchanged.
    const QDir markingDir(m_markingPath);
    const QString canonical = markingDir.filePath(m_subjectId + "_template_markings.bin");

    // ONE SOURCE. The .partial branch is gone with the cycle: mid-cycle used to
    // be a real state -- three of four openings of this window found a partial
    // file and had to restore PULSE only, leaving ECG at that alignment's fresh
    // auto-seed, because the partial's ECG marks belonged to a DIFFERENT
    // alignment than the one about to be shown. Every alignment is present at
    // once now, so a marking file is either a finished subject's (restore
    // everything) or absent (the fresh auto-seed stands).
    //
    // A leftover .partial from a pre-change session is deliberately ignored
    // rather than migrated: its ECG marks are keyed by anchor tag and would
    // restore correctly, but it may hold a half-finished cycle whose later
    // alignments were never marked, and silently presenting that as a restored
    // subject hides which bars are actually the operator's.
    bool markersReloaded = false;
    if (QFile::exists(canonical)) {
        markersReloaded = restoreMarkersFrom(canonical, /*ecg=*/true, /*pulse=*/true);
    }
    fprintf(stderr, "[markers] %s for subject %s (all %zu alignments)\n",
        markersReloaded ? "RELOADED prior markers" : "using FRESH auto-seed (no prior markers applied)",
        m_subjectId.toStdString().c_str(), anchor_view::kAllAnchors.size());

    // ---- TEMPORARY INSTRUMENTATION ------------------------------------
    // Everything above this point has already printed by the time the
    // operator sees a delay ([fast-phases], [ectopic], [markers]), so the
    // hang is in one of the two calls below or in Qt's first layout of what
    // showPage builds. If showPage reports a small number and the window is
    // still slow, the cost is the paint, not this function.
    using clk = std::chrono::steady_clock;
    auto t_prev = clk::now();
    auto lap = [&t_prev](const char* what) {
        const auto now = clk::now();
        fprintf(stderr, "[viewer] %-26s %7lld ms\n", what,
            (long long)std::chrono::duration_cast<std::chrono::milliseconds>(
                now - t_prev).count());
        fflush(stderr);
        t_prev = now;
        };

    computeGlobalRefs();
    lap("computeGlobalRefs");

    showPage();
    lap("showPage");
}

bool TemplateViewerWindow::restoreMarkersFrom(const QString& markingsBinPath,
    bool ecg, bool pulse) {
    try {
        std::vector<TemplateBin> saved = readTemplateMarkingsBin(markingsBinPath.toStdString());
        if (saved.empty()) {
            fprintf(stderr, "[markers] NOT reloaded from %s -- file has 0 bins\n",
                markingsBinPath.toStdString().c_str());
            return false;
        }
        if (saved.size() > m_bins.size()) {
            fprintf(stderr, "[markers] ERROR: %s has %zu bins but current subject has only %zu -- "
                "ignoring the extra %zu bin(s) in the file\n",
                markingsBinPath.toStdString().c_str(), saved.size(), m_bins.size(),
                saved.size() - m_bins.size());
        }
        else if (saved.size() < m_bins.size()) {
            fprintf(stderr, "[markers] ERROR: %s has only %zu bins but current subject has %zu -- "
                "leaving the extra %zu bin(s) at their fresh auto-seed (unedited)\n",
                markingsBinPath.toStdString().c_str(), saved.size(), m_bins.size(),
                m_bins.size() - saved.size());
        }
        const size_t n = std::min(saved.size(), m_bins.size());

        // Clamp every restored marker index against the CURRENT bin's own
        // template length before applying it. Marker indices from the saved
        // file were valid for whatever templates existed when it was
        // written -- if this subject's templates regenerated with different
        // per-bin lengths since then (which the bin-count mismatch above
        // already tells us is possible), blindly applying an old index that
        // now exceeds the current array length causes an out-of-bounds read
        // the first time anything indexes into that array with it (a real
        // crash, not just a cosmetic misplacement). Rejected markers keep
        // whatever seed_all() already put there.
        size_t rejectedCount = 0;
        auto safeIdx = [&](int savedVal, int currentVal, size_t len) -> int {
            if (savedVal >= 0 && static_cast<size_t>(savedVal) < len) return savedVal;
            if (savedVal >= 0) ++rejectedCount;   // only count real (non-sentinel) rejections
            return currentVal;
            };

        for (size_t i = 0; i < n; ++i) {
            TemplateBin& d = m_bins[i];
            const TemplateBin& s = saved[i];
            const size_t ecgLen[3] = {
                d.ch1.ecgTemplate_raw.size(), d.ch2.ecgTemplate_raw.size(), d.ch3.ecgTemplate_raw.size()
            };
            const size_t ppgLen = d.ppgTemplate.size();
            const size_t abpLen = d.abpTemplate.size();
            const size_t artLen = d.artTemplate.size();
            const size_t artPLen = d.artPulmTemplate.size();
            // ECG per-channel user markers (copied raw, bounds-checked). R is
            // NOT copied: it's an auto-only anchor and must come from this
            // pass's own template r_col (seeded above), never from a saved file.
            // Copy every anchor's saved marker set (each independent), bounds-
            // checked per channel. R is NOT copied (auto-only, re-derived).
            if (ecg) {
                // LEAD, SLOT, ANCHOR -- the shape the marking file now has.
                //
                // The saved bins carry their marks on placeholder bank slots
                // (see readTemplateMarkingsBin), so the walk is over the SAVED
                // slot count and each position is bounds-checked against the
                // DESTINATION template's own length. A saved slot with no
                // counterpart in this pass's bank is dropped: the bank is
                // rebuilt every run and its slot count can shrink, and a
                // landmark has no meaning without the waveform it sat on.
                for (int c = 0; c < 3; ++c) {
                    const int nSaved = s.markedSlotCount(c);
                    const int nNow =
                        static_cast<int>(d.ecg_bank[c].templates.size());
                    for (int slot = 0; slot < nSaved && slot < nNow; ++slot) {
                        // Slot 0's length is the bin's own template; a deeper
                        // slot's is its bank template's, which can be shorter.
                        // Using the bin's for both would admit a position past a
                        // sub-template's end.
                        const size_t len = (slot == 0)
                            ? ecgLen[c]
                            : d.ecg_bank[c].templates[slot].tmpl.size();
                        for (const auto& kv :
                            s.ecg_bank[c].templates[slot].markers_by_anchor) {
                            const tbank::BankMarkerSet& sm = kv.second;
                            tbank::BankMarkerSet& dm =
                                d.ecg_bank[c].templates[slot].marks(kv.first);
                            dm.p_begin = safeIdx(sm.p_begin, dm.p_begin, len);
                            dm.p_peak = safeIdx(sm.p_peak, dm.p_peak, len);
                            dm.q_begin = safeIdx(sm.q_begin, dm.q_begin, len);
                            dm.s_end = safeIdx(sm.s_end, dm.s_end, len);
                            dm.t_end = safeIdx(sm.t_end, dm.t_end, len);
                        }
                    }
                }
                for (int c = 0; c < 3; ++c) d.bad_r_ch[c] = s.bad_r_ch[c];
            } // if (ecg)

            if (pulse) {
                d.bad_ppg = s.bad_ppg;
                d.ppg_onset = safeIdx(s.ppg_onset, d.ppg_onset, ppgLen);
                d.ppg_t50 = safeIdx(s.ppg_t50, d.ppg_t50, ppgLen);
                d.ppg_t80 = safeIdx(s.ppg_t80, d.ppg_t80, ppgLen);
                d.ppg_peak = safeIdx(s.ppg_peak, d.ppg_peak, ppgLen);
                d.ppg_dicrotic = safeIdx(s.ppg_dicrotic, d.ppg_dicrotic, ppgLen);
                d.ppg_peak2 = safeIdx(s.ppg_peak2, d.ppg_peak2, ppgLen);
                d.ppg_end = safeIdx(s.ppg_end, d.ppg_end, ppgLen);
                d.abp_issue = s.abp_issue;
                d.abp_onset = safeIdx(s.abp_onset, d.abp_onset, abpLen);
                d.abp_peak = safeIdx(s.abp_peak, d.abp_peak, abpLen);
                d.abp_dicrotic = safeIdx(s.abp_dicrotic, d.abp_dicrotic, abpLen);
                d.abp_peak2 = safeIdx(s.abp_peak2, d.abp_peak2, abpLen);
                d.abp_end = safeIdx(s.abp_end, d.abp_end, abpLen);
                d.art_issue = s.art_issue;
                d.art_onset = safeIdx(s.art_onset, d.art_onset, artLen);
                d.art_peak = safeIdx(s.art_peak, d.art_peak, artLen);
                d.art_dicrotic = safeIdx(s.art_dicrotic, d.art_dicrotic, artLen);
                d.art_peak2 = safeIdx(s.art_peak2, d.art_peak2, artLen);
                d.art_end = safeIdx(s.art_end, d.art_end, artLen);
                d.art_pulm_issue = s.art_pulm_issue;
                d.art_pulm_onset = safeIdx(s.art_pulm_onset, d.art_pulm_onset, artPLen);
                d.art_pulm_peak = safeIdx(s.art_pulm_peak, d.art_pulm_peak, artPLen);
                d.art_pulm_dicrotic = safeIdx(s.art_pulm_dicrotic, d.art_pulm_dicrotic, artPLen);
                d.art_pulm_peak2 = safeIdx(s.art_pulm_peak2, d.art_pulm_peak2, artPLen);
                d.art_pulm_end = safeIdx(s.art_pulm_end, d.art_pulm_end, artPLen);
            } // if (pulse)
        }
        if (rejectedCount > 0) {
            fprintf(stderr, "[markers] WARNING: %zu marker(s) from %s were out of range for the "
                "current templates and were rejected (kept at auto-seed) instead of applied\n",
                rejectedCount, markingsBinPath.toStdString().c_str());
        }
        fprintf(stderr, "[markers] reloaded from %s (%zu of %zu bin(s) restored)\n",
            markingsBinPath.toStdString().c_str(), n, m_bins.size());
        return true;
    }
    catch (const std::exception& e) {
        fprintf(stderr, "[markers] NOT reloaded from %s -- %s\n",
            markingsBinPath.toStdString().c_str(), e.what());
        return false;   // No / unreadable file: caller keeps the fresh auto-seed.
    }
}

void TemplateViewerWindow::computeGlobalRefs() {
    /*compute the ecg global reference value for QRS complex height(abs(R) + abs(S))) and pulse global ref for
    PPG / ART / ART_PULM(abs(peak) - abs(foot))*/
    for (int c = 0; c < 3; ++c)
        m_ecgGlobalRef[c] = normalize_features::compute_ecg_global_ref(
            m_bins, c, m_sampleRate);
    for (int c = 0; c < 4; ++c) {
        m_pulseGlobalRef[c] = normalize_features::compute_pulse_global_ref(m_bins, c);
    }

    // ---- Sections 5.2-5.4: CV check + ratio/percentile normalization -----
    // Runs HERE because this is the one place TemplateBin data and a global
    // ref coexist (compute_ecg_global_ref itself is called just above). The
    // functions live in NormalizeFeatures.hpp; this drives them on real
    // per-bin data and persists the results.
    //
    // Decisions the spec left open, made explicit (override as needed):
    //  - Reference FEATURE = per-bin |R|+|S| (Option A's own basis), so the
    //    ratio and its reference are the same quantity. Computed per channel.
    //  - GLOBAL REF = all three options (A=|R|+|S|, B=QRS area, C=spatial)
    //    are written side by side rather than picking one; the ratio/pct
    //    columns use Option A to stay consistent with the reference feature.
    //  - p2/p98 for pct_scale come from THIS subject's own distribution of
    //    the per-bin ratio (per channel), matching the acceptance criterion
    //    "2nd and 98th percentiles at the extremes".
    if (!m_normOutputPath.isEmpty())
        writeNormalizationCsvs();
}

// Section 5.2-5.4 persistence. Split out of computeGlobalRefs for clarity.
void TemplateViewerWindow::writeNormalizationCsvs() {
    const std::string subj = m_subjectId.toStdString();
    const std::string dir = m_normOutputPath.toStdString();

    // Per-channel per-bin QRS reference (|R|+|S|), parallel to m_bins, NaN
    // for bins that don't yield one -- same extraction compute_ecg_global_ref
    // uses internally, just retained per bin instead of reduced to a median.
    auto perBinQrsRef = [&](int ch) {
        std::vector<double> v(m_bins.size(), std::numeric_limits<double>::quiet_NaN());
        for (size_t i = 0; i < m_bins.size(); ++i) {
            const auto& b = m_bins[i];
            if (b.bad_segment || b.bad_r_ch[ch]) continue;
            // THE R BASE, and this function has no alignment parameter to
            // choose otherwise. The per-bin QRS reference divides normalized
            // amplitudes subject-wide, so it must be one alignment for the
            // whole file -- the same reason writeTemplateMarkingsCsv pins its
            // ecgRef to R even inside a non-R sidecar.
            const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
            const auto& ecg = chs[ch]->ecgTemplate_raw;
            if (ecg.empty()) continue;
            // slotMarks selects the lead, so the per-lead subscripts below
            // are gone. Slot 0: the reference is measured on the sinus seed.
            const tbank::BankMarkerSet& rmk =
                b.slotMarks(ch, 0, AnchorType::R_PEAK);
            EcgFeatures f = computeEcgFeatures(
                ecg, rmk.p_peak, rmk.q_begin, b.r_peak_ch[ch],
                rmk.s_end, rmk.t_end, m_sampleRate);
            const double ry = normalize_features::sample_y(ecg, f.r_idx);
            const double sy = normalize_features::sample_y(ecg, f.s_idx);
            if (std::isnan(ry) || std::isnan(sy)) continue;
            v[i] = std::abs(ry) + std::abs(sy);
        }
        return v;
        };

    // ---- cv_check.csv : one row per (channel, bin) ----------------------
    // global_ref repeated per row (constant within a channel); cv_flag is
    // the per-channel record-level flag, also repeated. This is exactly the
    // long format the acceptance test reads.
    {
        std::ofstream f(dir + "/" + subj + "_cv_check.csv", std::ios::trunc);
        if (f) {
            f << "subject_id,channel,bin_index,qrs_ref_value,"
                "global_ref_A,global_ref_B,global_ref_C,cv_flag\n";
            for (int ch = 0; ch < 3; ++ch) {
                const std::vector<double> qref = perBinQrsRef(ch);
                const double grefA = normalize_features::compute_ecg_global_ref(m_bins, ch, m_sampleRate);
                const double grefB = normalize_features::compute_ecg_global_ref_area(m_bins, ch, m_sampleRate);
                const double grefC = normalize_features::compute_ecg_global_ref_spatial(m_bins);
                const bool flag = normalize_features::cv_flag(qref, grefA);
                for (size_t i = 0; i < qref.size(); ++i) {
                    f << subj << ",CH" << (ch + 1) << ',' << i << ',';
                    if (!std::isnan(qref[i])) f << qref[i];
                    f << ',';
                    if (!std::isnan(grefA)) f << grefA; f << ',';
                    if (!std::isnan(grefB)) f << grefB; f << ',';
                    if (!std::isnan(grefC)) f << grefC; f << ',';
                    f << (flag ? 1 : 0) << '\n';
                }
            }
        }
    }

    // ---- feature_norm.csv : one row per (channel, bin) ------------------
    // ratio = ratio_norm(qrsRef, grefA); feature_norm = pct_scale against
    // this channel's own p2/p98 of the finite ratios. By construction the
    // p2 and p98 rows land at 0 and 100 (clamped), which is the property the
    // acceptance test checks.
    {
        std::ofstream f(dir + "/" + subj + "_feature_norm.csv", std::ios::trunc);
        if (f) {
            f << "subject_id,channel,bin_index,ratio,feature_norm\n";
            for (int ch = 0; ch < 3; ++ch) {
                const std::vector<double> qref = perBinQrsRef(ch);
                const double grefA = normalize_features::compute_ecg_global_ref(m_bins, ch, m_sampleRate);

                std::vector<double> ratios(qref.size(), std::numeric_limits<double>::quiet_NaN());
                std::vector<double> finite;
                for (size_t i = 0; i < qref.size(); ++i) {
                    if (std::isnan(qref[i])) continue;
                    ratios[i] = normalize_features::ratio_norm(qref[i], grefA);
                    if (!std::isnan(ratios[i])) finite.push_back(ratios[i]);
                }
                // p2 / p98 of this channel's own ratio distribution.
                double p2 = std::nan(""), p98 = std::nan("");
                if (finite.size() >= 2) {
                    std::sort(finite.begin(), finite.end());
                    auto pctl = [&](double q) {
                        const double idx = (q / 100.0) * (finite.size() - 1);
                        const size_t lo = static_cast<size_t>(std::floor(idx));
                        const size_t hi = static_cast<size_t>(std::ceil(idx));
                        const double fr = idx - lo;
                        return finite[lo] * (1.0 - fr) + finite[hi] * fr;
                        };
                    p2 = pctl(2.0);
                    p98 = pctl(98.0);
                }
                for (size_t i = 0; i < ratios.size(); ++i) {
                    f << subj << ",CH" << (ch + 1) << ',' << i << ',';
                    if (!std::isnan(ratios[i])) f << ratios[i];
                    f << ',';
                    const double fn = normalize_features::pct_scale(ratios[i], p2, p98);
                    if (!std::isnan(fn)) f << fn;
                    f << '\n';
                }
            }
        }
    }
}

std::vector<double> TemplateViewerWindow::normalizeEcgTrace(const std::vector<double>& raw, int ch) const {
    const double ref = (ch >= 0 && ch < 3) ? m_ecgGlobalRef[ch] : std::nan("");
    return normalize_features::normalize_ecg_trace(raw, ref);
}

std::vector<double> TemplateViewerWindow::normalize_ppg_or_similar(const std::vector<double>& raw, int footIdx, int pulseChan) const {
    const double ref = (pulseChan >= 0 && pulseChan < 4) ? m_pulseGlobalRef[pulseChan] : std::nan("");
    return normalize_features::normalize_pulse_trace(raw, footIdx, ref);
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
    m_pageTemplateIdx.clear();
    m_pageGlobalIdx.clear();

    // Drop stretch factors left over from a previous (possibly larger) page
    // so unused rows/columns don't reserve empty space on the next page.
    for (int c = 0; c < ui->plotGrid->columnCount(); ++c)
        ui->plotGrid->setColumnStretch(c, 0);
    for (int r = 0; r < ui->plotGrid->rowCount(); ++r)
        ui->plotGrid->setRowStretch(r, 0);
}

void TemplateViewerWindow::showPage() {
    // ---- TEMPORARY INSTRUMENTATION ------------------------------------
    using clk = std::chrono::steady_clock;
    auto t_prev = clk::now();
    auto lap = [&t_prev](const char* what) {
        const auto now = clk::now();
        fprintf(stderr, "[showPage] %-24s %7lld ms\n", what,
            (long long)std::chrono::duration_cast<std::chrono::milliseconds>(
                now - t_prev).count());
        fflush(stderr);
        t_prev = now;
        };

    clearPlots();
    lap("clearPlots");

    // Page bounds come from the packed table, not from multiplication: pages
    // hold a variable number of bins so the COLUMN count stays bounded.
    if (m_pages.empty()) buildPages();
    m_currentPage = std::clamp(m_currentPage, 0, (int)m_pages.size() - 1);
    int start = m_pages[m_currentPage].first;
    int count = m_pages[m_currentPage].second;
    int end = start + count;

    bool compact = (m_maxLeads <= 1);

    // ---- expand this page's bins into (bin, template) COLUMNS -------------
    // A column used to be a bin. It is now a bin plus a bank member, so a bin
    // holding three morphologies occupies three adjacent columns and the sinus
    // seed stays leftmost. Columns per page therefore floats with the record's
    // ectopy: pages stay bin-aligned rather than column-count-aligned, because
    // splitting a bin across a page boundary would put its sinus template on
    // one page and its PVC template on the next.
    std::vector<std::pair<int, int>> cols;   // (global bin index, template index)
    for (int i = 0; i < count; ++i) {
        const int gi = start + i;
        for (int t : markingSlotsForBin(m_bins[gi])) cols.push_back({ gi, t });
    }
    const int nCols = static_cast<int>(cols.size());

    // WHICH BINS HAVE NO COLUMNS AT ALL. With a minimum-beats threshold set,
    // every template in a bin can fall below it, and that bin then shows
    // nothing. Printed rather than left to be noticed, because a bin that
    // silently has no panel is indistinguishable from a paging bug.
    if (tbank::minBeatsEcg() > 0 || tbank::minBeatsPpg() > 0) {
        std::string gone;
        for (int i = 0; i < count; ++i) {
            const int gi = start + i;
            if (markingSlotsForBin(m_bins[gi]).empty())
                gone += (gone.empty() ? "" : ",") + std::to_string(gi);
        }
        if (!gone.empty())
            std::fprintf(stderr,
                "  [min-beats] bins with NO displayable template "
                "(every slot below ECG %d or PPG %d clean beats): %s\n",
                tbank::minBeatsEcg(), tbank::minBeatsPpg(), gone.c_str());
    }

    // Says out loud whether the banks survived the trip from generation to the
    // viewer. nCols == count means every bin resolved to a single column, which
    // is what an EMPTY bank looks like -- and an empty bank is indistinguishable
    // on screen from a record with no ectopy, so the distinction has to be
    // printed rather than inferred from the plots.
    {
        int withBanks = 0;
        for (int i = 0; i < count; ++i)
            for (int c = 0; c < 3; ++c)
                if (m_bins[start + i].ecg_bank[c].size() > 1) { ++withBanks; break; }
        std::fprintf(stderr, "[bank-view] page bins=%d columns=%d "
            "bins_with_multi_template_bank=%d\n", count, nCols, withBanks);
    }

    // ---- row count -------------------------------------------------------
    // Computed HERE, not before the expansion above, because compact mode
    // wraps PANELS and a bin can contribute several. Sizing the wrap from the
    // bin count let a page with banks overflow past kMaxGridCols columns while
    // still reporting a legal grid.
    int gridRows;
    if (compact) {
        // Not `auto [rows, cols]`: `cols` is already the (bin, template)
        // vector above, and a structured binding would shadow it.
        const std::pair<int, int> rc = compactGrid(nCols);
        gridRows = std::max(1, rc.first);   // .second is implied by the wrap
    }
    else {
        // One extra row for the VCG panel, but only if some bin on this page
        // can actually produce one (all three ECG channels present with an R
        // column). Without this the row count equals the lead count and there
        // is no row index left for the VCG to occupy.
        bool anyVcg = false;
        for (int k = start; k < end && !anyVcg; ++k)
            anyVcg = !vcg_avg::derivedTraceOnChannelAxis(
                m_bins[k], 0, vcg::DerivedLead::VectorMagnitude).empty();
        // Capped like the compact path. Three leads plus VCG already sits at
        // the cap; a fourth lead would drop the VCG row rather than grow the
        // stack, and vcgRowWanted below (gridRows > m_maxLeads) agrees with
        // that on its own, so the two cannot disagree about which row exists.
        gridRows = std::min(n_template_rows, m_maxLeads + (anyVcg ? 1 : 0));
        // Note: this probe throws its trace away and the per-bin loop below
        // recomputes the same thing for the same bins. If this lap is large,
        // that duplication is the first thing to remove.
        lap("anyVcg probe");
    }

    m_binPlots.resize(nCols);
    m_pageGlobalIdx.resize(nCols);
    m_pageTemplateIdx.resize(nCols);

    int usedRows = 0, usedCols = 0;

    for (int i = 0; i < nCols; ++i) {
        int gi = cols[i].first;
        const int ti = cols[i].second;
        m_pageGlobalIdx[i] = gi;
        m_pageTemplateIdx[i] = ti;

        const TemplateBin& b = m_bins[gi];
        auto leads = leadsForBinTemplate(b, ti);

        // Only shout when no bank arrived at all: an empty bank and a genuinely
        // single-morphology bin render identically, so the difference has to be
        // written somewhere visible rather than inferred from the plots.
        if (ti == 0)
            for (auto& L : leads)
                if (L.channelIndex >= 0 && L.channelIndex < 3
                    && b.ecg_bank[L.channelIndex].size() == 0)
                    L.label += " bank=0";

        // ---- THE PULSE WAVEFORM, ITS BAND AND ITS COUNT ARE ALL PER GROUP -
        //
        // ppg_bank slot i IS group i: projectToChannel walks the joint groups
        // in order for every channel, so the pulse face is parallel to the
        // three ECG faces and sits on the same axis as the bin's own pulse
        // template (both are built from the same foot-anchored ppg_kept
        // beats). The bin's pulse marker columns stay valid against it.
        //
        // ALL THREE USED TO COME FROM THE BIN -- b.ppgTemplate, its
        // b.ppg_template_iqr, and b.ppg_n_beats -- drawn identically under
        // every column. So a panel could report one group's pulse count beside
        // a band computed from the whole bin, both morphologies and the full
        // respiratory amplitude swing, which is why a column with no pulses of
        // its own still showed a wide band.
        //
        // THERE IS NO FALLBACK. Not to the bin's template, not to its count,
        // not on a file whose ppg_bank is absent. A group with no pulse cohort
        // draws no pulse: no trace, no band, no count. Substituting the bin's
        // numbers is what made an unmeasured value look measured.
        //
        // ONE SOURCE FOR THE COUNT AND THE TRACE: the number printed is the
        // member count of the very template the drawn waveform is the median
        // of, so the two cannot disagree.
        //
        // tooFewBeats carries the PPG minimum from config.csv, so a cohort
        // below it is not drawn either -- a band fitted to two pulses is as
        // meaningless as a waveform built from two beats.
        tbank::BankTemplate* ppgSlot =
            (ti >= 0 && ti < b.ppg_bank.size())
            ? &m_bins[gi].ppg_bank.templates[ti] : nullptr;
        const bool hasPPG = ppgSlot && !ppgSlot->tmpl.empty()
            && ppgSlot->memberCount() > 0
            && !ppgSlot->tooFewBeats(/*is_ppg=*/true);
        const int nPpgForColumn = hasPPG ? ppgSlot->memberCount() : -1;

        // NO PULSE, NO PANEL -- AND THE TEST HAS TO BE HERE, not only in
        // markingSlotsForBin.
        //
        // shown() decides which slots are ELIGIBLE; this function builds and
        // populates the widgets, and it did so regardless of what shown()
        // returned. So tightening shown() alone left the same ECG-only panels on
        // screen: the column had already been created by the time hasPPG came
        // out false, and the only consequence was that no pulse trace was drawn.
        //
        // An ECG-only panel is worse than an absent one. It is
        // indistinguishable from a bin whose record genuinely has no pulse
        // channel, and the operator cannot place the pulse landmarks the panel
        // exists to collect.
        if (!hasPPG) continue;   // unconditional -- see leadsForBinTemplate

        if (leads.empty())
            leads.push_back({ nullptr, 0, "No ECG" });

        // VCG needs all three channels; with fewer, the trace comes back empty
        // and no row is reserved, so the leads keep the full height.
        double vcgRCol = -1.0;
        const std::vector<double> vcgTrace = vcg_avg::derivedTraceOnChannelAxis(
            b, 0, vcg::DerivedLead::VectorMagnitude, vcg::kIdentity, &vcgRCol);
        const bool vcgRowWanted = !compact && !vcgTrace.empty()
            && gridRows > m_maxLeads;

        // R frame: the reference lines these produce are drawn on the
        // R-aligned panels, and leadMarkersFor's USER path assembles the bars
        // from all four alignments into that frame.
        const auto gi_intervals = global_intervals::computeGlobalIntervals(
            b, AnchorType::R_PEAK, m_sampleRate,
            global_intervals::MarkerSource::USER);

        std::vector<BinPlotWidget*> group;

        for (int li = 0; li < (int)leads.size(); ++li) {
            auto* pw = new BinPlotWidget(gi, leads[li].channelIndex, leads[li].label);
            pw->setChannelRate(BinPlotWidget::Channel::Ecg, m_sampleRate);
            pw->setChannelRate(BinPlotWidget::Channel::Ppg, m_ppgRateHz);
            pw->setChannelRate(BinPlotWidget::Channel::Abp, m_abpRateHz);
            pw->setChannelRate(BinPlotWidget::Channel::Art, m_artRateHz);
            pw->setChannelRate(BinPlotWidget::Channel::ArtPulm, m_artPulmRateHz);

            static const std::vector<double> empty;
            const auto& ecg = leads[li].ecg ? *leads[li].ecg : empty;
            const auto& ppg = hasPPG ? ppgSlot->tmpl : empty;

            int c = leads[li].channelIndex;
            const double rPeak = static_cast<double>(b.r_peak_ch[c]);

            const std::vector<double>& ecgIqrRaw = (c == 0) ? b.ch1.ecg_template_raw_iqr
                : (c == 1) ? b.ch2.ecg_template_raw_iqr
                : b.ch3.ecg_template_raw_iqr;
            const double ecgRef = (c >= 0 && c < 3) ? m_ecgGlobalRef[c] : std::nan("");
            // Both ecgIqrRaw (Q3-Q1 of raw amplitude) and b.ppg_template_iqr
            // (Q3-Q1 of each beat's own local perfusion-index ratio, computed
            // at build time -- see CreatePPGTemplates.hpp) are pre-ref-division.
            // The only remaining step, for either channel type, is the same
            // scalar /ref used for the mean trace -- normalization lives only
            // in NormalizeFeatures.hpp, this is just calling it.
            const std::vector<double> ecgIqr = normalize_features::scale_array_by_ref(ecgIqrRaw, ecgRef);
            // The GROUP's own spread. b.ppg_template_iqr is the interquartile
            // range across every pulse in the bin, so under a group holding one
            // morphology it drew a band sized by the difference BETWEEN
            // morphologies.
            std::vector<double> ppgIqr = empty;
            int ppgFootIdx = -1;
            if (hasPPG) {
                if (!ppgSlot->hasDetectedPulseMarks())
                    FeatureMarks::seed_pulse_bank_template(ppgSlot->tmpl,
                        m_ppgRateHz, ppgSlot->pulse_marks);
                ppgFootIdx = ppgSlot->pulse_marks.onset;
                ppgIqr = normalize_features::scale_pulse_spread_by_ref(
                    ppgSlot->tmpl_iqr,
                    normalize_features::sample_y(ppgSlot->tmpl, ppgFootIdx),
                    m_pulseGlobalRef[0]);
            };
            // ECG beats for THIS TEMPLATE, not for the bin. Lead::nMembers is
            // the bank member's own beat count; it is 0 only on the pre-bank
            // fallback path (slot 0 of a bin whose bank never arrived), where
            // the bin total is the only count that exists and the template IS
            // the whole bin, so the two agree.
            //
            // The pulse count is per group too, resolved above from the pulse
            // bank slot. Nothing on this line derives it.
            const uint64_t nEcgBinTotal = (c == 0) ? b.ch1_n_beats_raw
                : (c == 1) ? b.ch2_n_beats_raw
                : b.ch3_n_beats_raw;
            const uint64_t nEcgBeats = (leads[li].nMembers > 0)
                ? static_cast<uint64_t>(leads[li].nMembers)
                : nEcgBinTotal;

            // Per-subject normalized traces for on-screen display.
            // ECG:   sample / Global_Ref_ecg(ch)
            // Pulse: ((sample - foot_y) / |foot_y| * 100) / Global_Ref_pulse(chan)
            // If the ref or foot is unusable, the helpers pass the raw
            // trace through unchanged.
            //
            // Display-time notch filter toggle: when the `notch_filter`
            // checkbox is on AND a valid notch frequency was passed in via
            // loadSubject, push each channel's template through notch_filter
            // at that channel's own rate BEFORE normalization. Purely a
            // viewing convenience: templates on disk are untouched, and
            // toggling off returns to the raw stored templates on the next
            // showPage() (which the checkbox handler triggers).
            //
            // maybeNotch also REBASES: notch_filter (implemented as
            // x - narrow_bandpass(x)) leaves the waveform shape intact but
            // shifts the DC level by whatever tiny residual mean the
            // bandpass emits. That shift moves foot_y (the value at the foot
            // column), which is what normalize_pulse_trace divides by. For
            // pulse channels stored near-baseline-subtracted, foot_y sits
            // near zero, so even a tiny DC shift is a LARGE relative change
            // to foot_y -- the normalized trace then collapses to the
            // bottom of the y-axis (looks like "PPG went away"). Rebasing
            // to the pre-notch foot value keeps the normalization stable
            // and the visible effect really is just the notch, not a
            // baseline-hunt-induced squash.

            const bool notchActive = m_notchFilterOn && m_notchFilterHz > 0;
            auto maybeNotch = [&](const std::vector<double>& sig, double fs, int footIdx) -> std::vector<double> {
                if (!notchActive || sig.empty() || fs <= 0.0) return sig;
                std::vector<double> out = notch_filter(sig, static_cast<double>(m_notchFilterHz), fs);
                if (footIdx >= 0 && footIdx < (int)sig.size() && footIdx < (int)out.size()) {
                    const double shift = sig[footIdx] - out[footIdx];
                    if (std::isfinite(shift) && shift != 0.0)
                        for (auto& v : out) v += shift;
                }
                return out;
                };
            const std::vector<double> ecgSrc = maybeNotch(ecg, m_sampleRate, -1);   // ECG uses /ref, not a foot; no rebase needed
            const std::vector<double> ppgSrc = hasPPG
                ? maybeNotch(ppgSlot->tmpl, m_ppgRateHz, b.ppg_onset) : empty;
            const std::vector<double> abpSrc = !b.abpTemplate.empty() ? maybeNotch(b.abpTemplate, m_abpRateHz, b.abp_onset) : b.abpTemplate;
            const std::vector<double> artSrc = !b.artTemplate.empty() ? maybeNotch(b.artTemplate, m_artRateHz, b.art_onset) : b.artTemplate;
            const std::vector<double> artPSrc = !b.artPulmTemplate.empty() ? maybeNotch(b.artPulmTemplate, m_artPulmRateHz, b.art_pulm_onset) : b.artPulmTemplate;

            const std::vector<double> ecgN = normalizeEcgTrace(ecgSrc, c);
            // ppgFootIdx, not b.ppg_onset: the trace being normalized is
            // ppgSlot->tmpl, and b.ppg_onset was measured on b.ppgTemplate. The
            // band above uses this same foot, so the two agree by construction.
            const std::vector<double> ppgN = hasPPG
                ? normalize_ppg_or_similar(ppgSrc, ppgFootIdx, 0)
                : empty;
            const std::vector<double> abpN = !abpSrc.empty()
                ? normalize_ppg_or_similar(abpSrc, b.abp_onset, 1)
                : abpSrc;
            const std::vector<double> artN = !artSrc.empty()
                ? normalize_ppg_or_similar(artSrc, b.art_onset, 2)
                : artSrc;
            const std::vector<double> artPN = !artPSrc.empty()
                ? normalize_ppg_or_similar(artPSrc, b.art_pulm_onset, 3)
                : artPSrc;

            // (No plot-width call here any more. The frame is derived from the
            // traces themselves in BinPlotWidget::recomputeFrame -- the union of
            // every channel's own drawn extent, in seconds relative to R -- so
            // there is no width to hand in ahead of setData. median_rr_samples
            // was never populated anyway: copyMethod does not set it and
            // readMethod does not deserialize it, which is why the widget had
            // been falling back to the array length.)

            // Traces only. Every marker and autodetect column arrives via the
            // single applyBinToWidget() call below. nPpgForColumn was resolved
            // once per column above, from the pulse bank slot this panel draws.

            pw->setData(ppgN, ppgIqr, ecgN, ecgIqr, rPeak,
                static_cast<int>(nEcgBeats), nPpgForColumn);

            pw->setHasPPG(hasPPG);

            // Faint arterial background-context traces (present-only),
            // foot-anchored like the PPG. Colors mirror the noise-marking GUI:
            // ABP teal, ART dark red, ART_PULM dark blue.
            {
                std::vector<std::pair<std::vector<double>, QColor>> bg;
                if (!abpN.empty())
                    bg.push_back({ abpN,  QColor(0, 95, 105) });
                if (!artN.empty())
                    bg.push_back({ artN,  QColor(150, 40, 40) });
                if (!artPN.empty())
                    bg.push_back({ artPN, QColor(40, 60, 150) });
                pw->setBackgroundTraces(bg);
            }

            // Arterial traces for marker geometry/bounds. The arterial markers
            // themselves come from applyBinToWidget() with all the others.
            pw->setArterialTraces(abpN, artN, artPN,
                b.abpTemplate_iqr, b.artTemplate_iqr, b.artPulmTemplate_iqr);

            // Seed every bar + every autodetect column, in one call, after all
            // traces are in place (the glyph capture needs them).
            //
            // SLOT 0 ONLY. TemplateBin::marks() holds one MarkerSet per anchor
            // for the bin, which describes the sinus template. A PVC's Q-onset
            // sits at a different column than sinus's, so applying those bars to
            // a bank column would draw landmarks that are simply wrong -- and a
            // drag would then write them back. Bank members carry their own
            // BankMarkerSet (tbank::BankTemplate::markers_by_anchor); until the
            // marking path is threaded through to them, an ectopic column shows
            // its waveform with no bars, which is honest.
            // Slot 0 gets the bin's marker set, which describes the sinus
            // template. Every other slot gets its OWN BankMarkerSet, seeded
            // from its own median waveform -- not a copy of the bin's, because a
            // PVC's Q-onset sits at a different column than sinus's and drawing
            // sinus's bars there would be wrong in a way a drag would then
            // persist. That was the reason sub-templates had no bars at all.
            if (ti == 0) applyBinToWidget(pw, b);
            // m_bins[gi] rather than `b`: the loop binds `b` as const, and the
            // lazy seed writes the marker set it just computed back into the
            // template so the next repaint and any drag see the same positions.
            else         applyBankTemplateToWidget(pw, m_bins[gi], c, ti);

            pw->setReferenceLines(global_interval_lines::forChannel(b, gi_intervals, c));

            // ---- RESTORE THIS PANEL'S OWN MARK -----------------------
            // operator_state is per template. It used to read b.bad_ppg and
            // b.bad_r_ch[c], which are per BIN, so a rebuild painted every
            // panel of a marked bin -- the marks spread on a page turn.
            //
            // Slot 0 falls back to the bin flags when it has no state of its
            // own, so a record marked before operator_state existed, or one
            // whose flags feature_marks set automatically, still shows them.
            {
                uint8_t st = 0;
                const tbank::TemplateBank& bkq =
                    (c >= 0 && c <= 2) ? b.ecg_bank[c] : b.ppg_bank;
                if (ti >= 0 && ti < bkq.size())
                    st = bkq.templates[ti].operator_state;
                if (st == 0 && ti == 0)
                    st = (b.bad_ppg == 1) ? 2u
                    : ((c >= 0 && c <= 2 && b.bad_r_ch[c]) ? 1u : 0u);
                if (st == 2)      pw->setState(BinPlotWidget::State::BadPPG);
                else if (st == 1) pw->setState(BinPlotWidget::State::BadR);
            }

            connect(pw, &BinPlotWidget::markerMovedOnTemplate,
                this, &TemplateViewerWindow::onMarkerMovedOnTemplate);
            connect(pw, &BinPlotWidget::markerDragStarted,
                this, &TemplateViewerWindow::onMarkerDragStarted);
            connect(pw, &BinPlotWidget::landmarkSelected,
                this, &TemplateViewerWindow::onLandmarkSelected);
            connect(pw, &BinPlotWidget::badRToggled,
                this, &TemplateViewerWindow::onBadRToggled);
            pw->setTemplateIndex(ti);
            connect(pw, &BinPlotWidget::classConfirmRequested,
                this, &TemplateViewerWindow::onClassConfirmRequested);
            connect(pw, &BinPlotWidget::badPPGToggled,
                this, &TemplateViewerWindow::onBadPPGToggled);

            if (compact) {
                int row = i % gridRows;
                int col = i / gridRows;
                ui->plotGrid->addWidget(pw, row, col);
                usedRows = std::max(usedRows, row + 1);
                usedCols = std::max(usedCols, col + 1);
            }
            else {
                // The last lead no longer swallows the remaining rows: the
                // bottom row belongs to the VCG panel added after this loop.
                int rowspan = (li == (int)leads.size() - 1)
                    ? std::max(1, gridRows - li - (vcgRowWanted ? 1 : 0)) : 1;
                ui->plotGrid->addWidget(pw, li, i, rowspan, 1);
                usedRows = std::max(usedRows, li + rowspan);
                usedCols = std::max(usedCols, i + 1);
            }

            m_allPlots.push_back(pw);
            group.push_back(pw);
        }

        // ------------------------------------------------------------------
        // VCG panel: the bottom row of each bin's column, under lead 3.
        // ------------------------------------------------------------------
        // The trace is laid out on ch1's COLUMN axis by
        // derivedTraceOnChannelAxis, so it shares the x axis of the lead
        // panels above and lines up with them sample for sample -- while every
        // channel is still sampled at its own r_col internally, which is what
        // keeps the combination per-instant.
        //
        // Display only: no markers are set and marker signals are not
        // connected, so nothing here is draggable. The VCG is derived from the
        // three leads, so marking it would create a fourth set of fiducials
        // with no channel of its own to store them in.
        if (vcgRowWanted) {
            const int vcgRow = gridRows - 1;
            auto* vp = new BinPlotWidget(gi, vcgRow, "VCG", this);

            static const std::vector<double> emptyVec;
            vp->setChannelRate(BinPlotWidget::Channel::Ecg, m_sampleRate);
            vp->setData(emptyVec, emptyVec, vcgTrace, emptyVec,
                (vcgRCol >= 0.0) ? vcgRCol : 0.0, 0, 0);
            vp->setHasPPG(false);
            vp->setShowPpgTrace(false);
            vp->setShowEcgMarkers(false);
            vp->setShowPpgMarkers(false);
            vp->setShowPpgDerivMarkers(false);
            vp->setShowAbpMarkers(false);
            vp->setShowArtMarkers(false);
            vp->setShowArtPulmMarkers(false);

            // Same global boundaries as the leads above, converted onto ch1's
            // axis since that is the axis this trace is drawn on.
            vp->setReferenceLines(
                global_interval_lines::forChannel(b, gi_intervals, 0));

            ui->plotGrid->addWidget(vp, vcgRow, i, 1, 1);
            usedRows = std::max(usedRows, vcgRow + 1);
            usedCols = std::max(usedCols, i + 1);

            m_allPlots.push_back(vp);
            group.push_back(vp);
        }

        m_binPlots[i] = std::move(group);
        lap("  one bin's widgets");
    }
    lap("all bin widgets");

    // Equal stretch on every used row/column => equal-width, equal-height
    // cells that together fill the whole plot area. Combined with each
    // widget scaling its trace to its cell width, every bin window ends up
    // the same on-screen length.
    for (int c = 0; c < usedCols; ++c) ui->plotGrid->setColumnStretch(c, 1);
    for (int r = 0; r < usedRows; ++r) ui->plotGrid->setRowStretch(r, 1);

    applyMarkerVisibility();
    updatePageControls();
    lap("stretch+visibility+controls");
}

void TemplateViewerWindow::captureCurrentPage() {
    // Snapshot the page the user is LEAVING (captures its final edited state).
    // Grabbed synchronously so we snap the page still on screen, not the next
    // one; re-captures on each leave so the latest state wins.
    QDir outDir(m_templateDir);
    const int page = m_currentPage;
    const QPixmap shot = ui->scrollContents->grab();
    // No alignment suffix: there is one screenshot per page now, not one per
    // page per pass. The grid it captures is R-aligned on every panel.
    const QString fn = outDir.filePath(
        QString("%1_templates_page%2.png")
        .arg(m_subjectId)
        .arg(page + 1, 2, 10, QChar('0')));
    shot.save(fn, "PNG");
}

// Suffix every column in `header` (a single header line) with `suffix`,
// EXCEPT the first three (file_id, bin_num, x_ms), which are row keys and
// must stay un-suffixed so the R and Q sides line up when zipped.
static std::string suffixValueColumns(const std::string& header, const std::string& suffix) {
    std::string out;
    out.reserve(header.size() + 32);
    size_t start = 0;
    int colIdx = 0;
    for (size_t i = 0; i <= header.size(); ++i) {
        if (i == header.size() || header[i] == ',') {
            out.append(header, start, i - start);
            if (colIdx >= 3) out.append(suffix);   // skip file_id/bin_num/x_ms
            if (i < header.size()) out.push_back(',');
            start = i + 1;
            ++colIdx;
        }
    }
    return out;
}

// Zip the R canonical file with the just-produced Q content (both as raw CSV
// text). Header line: R header with un-suffixed keys and _r on values, then
// Q's non-key value columns with _q on them. Row lines: paired 1:1; each
// zipped row uses R's file_id/bin_num/x_ms as the keys, then R's value
// columns, then Q's value columns. Row counts must match -- if they don't we
// abort and leave the R file untouched (worse to write a garbled file than
// none). Returns true on success.
// Merge an ordered list of sidecar CSV files into one canonical file.
// files[0] is the base (kept whole, with its keys); each subsequent file
// contributes only its value columns (first 3 key columns stripped), appended
// to every row -- same column convention as zipCanonicalWithQ, generalized to
// N inputs and done ONCE at the end instead of a growing read-modify-write per
// pass. Row counts must match across all files. Returns true on success.
static bool mergeSidecarCsvs(const std::string& canonicalPath,
    const std::vector<std::string>& sidecarPaths)
{
    if (sidecarPaths.empty()) return false;

    auto readLines = [](const std::string& p) {
        std::vector<std::string> lines; std::string ln;
        std::ifstream in(p);
        if (!in) return lines;
        while (std::getline(in, ln)) {
            if (!ln.empty() && ln.back() == '\r') ln.pop_back();
            lines.push_back(ln);
        }
        return lines;
        };
    auto stripFirstThree = [](const std::string& line) -> std::string {
        int commas = 0;
        for (size_t i = 0; i < line.size(); ++i)
            if (line[i] == ',' && ++commas == 3) return line.substr(i + 1);
        return {};
        };

    std::vector<std::string> merged = readLines(sidecarPaths[0]);
    if (merged.empty()) {
        fprintf(stderr, "[tmplcsv] merge: base sidecar empty/missing: %s\n", sidecarPaths[0].c_str());
        return false;
    }
    for (size_t s = 1; s < sidecarPaths.size(); ++s) {
        std::vector<std::string> next = readLines(sidecarPaths[s]);
        if (next.size() != merged.size()) {
            fprintf(stderr, "[tmplcsv] merge: row mismatch (%zu vs %zu) for %s -- skipping\n",
                next.size(), merged.size(), sidecarPaths[s].c_str());
            continue;   // skip a bad sidecar rather than abort the whole merge
        }
        for (size_t i = 0; i < merged.size(); ++i) {
            const std::string tail = stripFirstThree(next[i]);
            if (!tail.empty()) { merged[i] += ','; merged[i] += tail; }
        }
    }

    std::ofstream out(canonicalPath, std::ios::trunc);
    if (!out) {
        fprintf(stderr, "[tmplcsv] merge: cannot write %s\n", canonicalPath.c_str());
        return false;
    }
    for (const auto& l : merged) out << l << '\n';
    return out.good();
}

static bool zipCanonicalWithQ(const std::string& canonicalPath,
    const std::string& qContent)
{
    std::ifstream rf(canonicalPath);
    if (!rf) {
        fprintf(stderr, "[tmplcsv] zip: cannot open R canonical %s\n", canonicalPath.c_str());
        return false;
    }
    auto readAll = [](std::istream& in) {
        std::vector<std::string> lines; std::string ln;
        while (std::getline(in, ln)) {
            if (!ln.empty() && ln.back() == '\r') ln.pop_back();
            lines.push_back(ln);
        }
        return lines;
        };
    std::vector<std::string> R = readAll(rf);
    rf.close();
    std::istringstream qs(qContent);
    std::vector<std::string> Q = readAll(qs);
    if (R.empty() || Q.empty()) {
        fprintf(stderr, "[tmplcsv] zip: empty input (R=%zu Q=%zu)\n", R.size(), Q.size());
        return false;
    }
    if (R.size() != Q.size()) {
        fprintf(stderr, "[tmplcsv] zip: row count mismatch (R=%zu Q=%zu) -- keeping R file as-is\n",
            R.size(), Q.size());
        return false;
    }

    // Strip Q's first three columns (the shared keys) from every line so
    // we don't duplicate them in the zipped output.
    auto stripFirstThree = [](const std::string& line) -> std::string {
        int commas = 0;
        for (size_t i = 0; i < line.size(); ++i) {
            if (line[i] == ',') {
                if (++commas == 3) return line.substr(i + 1);
            }
        }
        return {};   // fewer than 3 commas -- nothing to append
        };

    std::ofstream out(canonicalPath, std::ios::trunc);
    if (!out) {
        fprintf(stderr, "[tmplcsv] zip: cannot rewrite %s\n", canonicalPath.c_str());
        return false;
    }
    for (size_t i = 0; i < R.size(); ++i) {
        const std::string qTail = stripFirstThree(Q[i]);
        out << R[i];
        if (!qTail.empty()) out << ',' << qTail;
        out << '\n';
    }
    return out.good();
}

void TemplateViewerWindow::writeAlignedTemplateCsv(AnchorType anchor) {
    // One row per sample column per bin. Every value column is emitted twice,
    // raw and normalized, plus its IQR band in both forms; marker positions
    // ride along as one-hot flags on the row matching their sample index.
    if (m_bins.empty()) return;

    QDir outDir(m_templateDir);
    if (!outDir.exists()) outDir.mkpath(".");
    // Names the canonical merge target. Used only in the log/error messages
    // below -- this function writes a per-anchor sidecar, never `path` itself.
    const QString path = outDir.filePath(m_subjectId + "_template.csv");

    std::ostringstream f;
    const double toMs = (m_sampleRate > 0.0) ? 1000.0 / m_sampleRate : 1.0;

    // ---- Column lists (shared by the header pass and the row loop) ---------
    static const char* CHANS[] = {
        "ch1", "ch2", "ch3", "ppg", "abp", "art", "art_pulm"
    };
    constexpr int num_chans = static_cast<int>(std::size(CHANS));

    // Per ECG channel, in emission order. r_peak is auto-only: there is no
    // user bar for it, so it contributes one column where the others
    // contribute two. kRPeak is derived from the array rather than written as
    // a literal, because the index was previously hardcoded as `4` in two
    // places and inserting a marker ahead of it would have silently moved the
    // auto-only column onto the wrong landmark.
    // t_begin removed with the marker (see BankMarkerSet): its two columns
    // were structurally blank, because nothing ever set it.
    static const char* ECG_MARKERS[] = {
        "p_begin", "p_peak", "q_begin", "q_peak", "r_peak", "s_peak",
        "s_end",   "t_end"
    };
    constexpr int kNumEcgMarkers = static_cast<int>(std::size(ECG_MARKERS));
    constexpr int kRPeak = 4;
    static_assert(kNumEcgMarkers == 8, "ECG marker count changed; check kRPeak");
    // kRPeak cannot be pinned by static_assert: anchor_view::hasUserColumn
    // compares strings with std::strcmp and ECG_MARKERS is an array of
    // pointers, so neither is a constant expression. The count assert above is
    // the guard -- inserting a marker changes it, which forces a look at
    // kRPeak. (writeTemplateMarkingsCsv has no such gap: it keys off the
    // column NAME through hasUserColumn at run time, so it needs no index.)
    assert(!anchor_view::hasUserColumn(ECG_MARKERS[kRPeak])
        && "kRPeak must index the one marker with no user column");

    // Pulse marker groups. Counts differ by channel: PPG carries the two
    // interpolated upslope crossings (t50/t80), the arterial channels do not.
    static const char* PPG_MARKERS[] = {
        "ppg_onset", "ppg_t50", "ppg_peak",
        "ppg_dicr", "ppg_peak2", "ppg_t80", "ppg_end"
    };
    static const char* ABP_MARKERS[] = {
        "abp_onset", "abp_peak", "abp_dicr", "abp_peak2", "abp_end"
    };
    static const char* ART_MARKERS[] = {
        "art_onset", "art_peak", "art_dicr", "art_peak2", "art_end"
    };
    static const char* ARTP_MARKERS[] = {
        "art_pulm_onset", "art_pulm_peak", "art_plm_dicr",
        "art_pulm_peak2", "art_pulm_end"
    };
    constexpr int kNumPpgMarkers = static_cast<int>(std::size(PPG_MARKERS));
    constexpr int kNumArterialMarkers = static_cast<int>(std::size(ABP_MARKERS));

    // Per-ECG-channel autodetect glyph group. p_wave/q_onset/r_wave are the
    // stored autodetect columns, emitted directly so no parallel recompute can
    // disagree with them; t_peak is the only derived one.
    static const char* ECG_GLYPHS[] = { "p_wave", "q_onset", "r_wave", "t_peak" };

    // ---- Header ------------------------------------------------------------
    f << "file_id,bin_num,x_ms";
    for (const char* n : CHANS) {
        f << ',' << n << "_raw_mv"
            << ',' << n << "_norm"
            << ',' << n << "_raw_std"
            << ',' << n << "_norm_std";
    }
    for (int c = 1; c <= 3; ++c) {
        for (int k = 0; k < kNumEcgMarkers; ++k) {
            f << ',' << ECG_MARKERS[k] << "_ch" << c << "_auto";
            // WAS "_auto" A SECOND TIME -- a copy-paste that gave every
            // landmark two identically-named columns and no user column at
            // all, so a reader keying on name got whichever pandas kept.
            if (k != kRPeak)
                f << ',' << ECG_MARKERS[k] << "_ch" << c << "_user";
        }
    }
    auto emitPulseHeaderGroup = [&](auto const& group) {
        for (const char* m : group)
            f << ',' << m << "_auto"
            << ',' << m << "_user";
        };
    emitPulseHeaderGroup(PPG_MARKERS);
    emitPulseHeaderGroup(ABP_MARKERS);
    emitPulseHeaderGroup(ART_MARKERS);
    emitPulseHeaderGroup(ARTP_MARKERS);

    // Autodetected computed feature locations (no user bar).
    for (int gc = 1; gc <= 3; ++gc) {
        char gb[64];
        for (const char* g : ECG_GLYPHS) {
            std::snprintf(gb, sizeof gb, "%s_ch%d", g, gc);
            f << ',' << gb << "_auto";
        }
    }
    // Driven off ppg_and_artpulse_automated_markers, the same table the row
    // loop emits from. This used to be a hand-written list of four names while
    // the table had grown to fifteen entries, so the file carried eleven
    // unnamed columns and every header-keyed reader was reading the wrong
    // field from here to the end of the row.
    for (const auto& gl : ppg_and_artpulse_automated_markers)
        f << ',' << gl.name << "_auto";
    f << '\n';

    f << std::setprecision(10);

    // Both ecg_template_raw_iqr and every *_template_iqr field are
    // pre-ref-division at build time (see CreatePPGTemplates.hpp /
    // build_pulse_template_pair_windowed), so the only remaining step for any
    // channel is the same scalar /ref used for its mean trace --
    // normalize_features::scale_array_by_ref is the one place that happens.

    // ---- Row loop ----------------------------------------------------------
    for (size_t bi = 0; bi < m_bins.size(); ++bi) {
        const TemplateBin& b = m_bins[bi];

        // One description per value column, in CHANS order. The seven channels
        // differ only in which normalizer they take and which global reference
        // scales their IQR, so they are described rather than transcribed:
        // the previous form repeated four near-identical statements per channel
        // and the ordering of CHANS was only implicitly matched.
        struct Src {
            const std::vector<double>* raw;
            const std::vector<double>* rawIqr;
            double ref;      ///< global reference for this channel's scaling
            int    idx;      ///< channel index handed to the normalizer
            bool   isEcg;    ///< selects the normalizer
            int    onset;    ///< pulse channels only: alignment onset
        };
        const Src src[num_chans] = {
            // THIS SIDECAR'S ALIGNMENT. One file per alignment, each holding
            // that alignment's own averages -- which is what makes the merged
            // <id>_template.csv four aligned views of the same subject rather
            // than four copies of one.
            { &b.chFor(0, anchor).ecgTemplate_raw, &b.chFor(0, anchor).ecg_template_raw_iqr, m_ecgGlobalRef[0],   0, true,  0 },
            { &b.chFor(1, anchor).ecgTemplate_raw, &b.chFor(1, anchor).ecg_template_raw_iqr, m_ecgGlobalRef[1],   1, true,  0 },
            { &b.chFor(2, anchor).ecgTemplate_raw, &b.chFor(2, anchor).ecg_template_raw_iqr, m_ecgGlobalRef[2],   2, true,  0 },
            { &b.ppgTemplate,         &b.ppg_template_iqr,         m_pulseGlobalRef[0], 0, false, b.ppg_onset },
            { &b.abpTemplate,         &b.abpTemplate_iqr,          m_pulseGlobalRef[1], 1, false, b.abp_onset },
            { &b.artTemplate,         &b.artTemplate_iqr,          m_pulseGlobalRef[2], 2, false, b.art_onset },
            { &b.artPulmTemplate,     &b.artPulmTemplate_iqr,      m_pulseGlobalRef[3], 3, false, b.art_pulm_onset },
        };

        std::vector<double> norm[num_chans], normIqr[num_chans];
        for (int k = 0; k < num_chans; ++k) {
            norm[k] = src[k].isEcg
                ? normalizeEcgTrace(*src[k].raw, src[k].idx)
                : normalize_ppg_or_similar(*src[k].raw, src[k].onset, src[k].idx);
            normIqr[k] = normalize_features::scale_array_by_ref(*src[k].rawIqr, src[k].ref);
        }

        // Row span = longest trace; all traces start at row 0. Marker indices
        // are NOT considered: a landmark past the end of every trace has no row
        // to flag, which is a template-construction problem rather than
        // something this writer can represent.
        int hiRow = 0;
        for (const Src& s : src) hiRow = std::max(hiRow, (int)s.raw->size());
        if (hiRow <= 0) continue;

        // Computed Q/S peaks per ECG channel, from the auto bars and the user
        // bars separately.
        // THIS SIDECAR'S ALIGNMENT, matching the traces above. Reading
        // ch1..3 here would measure this alignment's landmarks against the
        // R-aligned average.
        const ChannelTemplateData* chs[3] = {
            &b.chFor(0, anchor), &b.chFor(1, anchor), &b.chFor(2, anchor) };
        EcgFeatures ftAuto[3], ftUser[3];
        for (int c = 0; c < 3; ++c) {
            const auto& ecg = chs[c]->ecgTemplate_raw;
            const auto aaF = b.autoFor(anchor);
            ftAuto[c] = computeEcgFeatures(ecg,
                (int)std::lround(aaF.p_peak[c]), (int)std::lround(aaF.q_begin[c]),
                (int)std::lround(aaF.r_peak[c]), (int)std::lround(aaF.s_end[c]),
                (int)std::lround(aaF.t_end[c]), m_sampleRate);
            // Per lead, because slotMarks selects the lead -- the old bin-wide
            // MarkerSet held all three leads in one object and was fetched once
            // per bin.
            const tbank::BankMarkerSet umk = b.userMarks(c, 0, anchor);
            // userMarks returns BARS ONLY, so umk.p_peak is -1 by
            // construction: P peak is a glyph. Recompute it from the two bars
            // that bracket it rather than passing the -1 through, which would
            // have handed computeEcgFeatures an absent P on every lead.
            const FeatureMarks::ReactiveEcg rxF = FeatureMarks::reactive_ecg(
                ecg, umk.p_begin, umk.q_begin, umk.s_end, umk.t_end, m_sampleRate);
            ftUser[c] = computeEcgFeatures(ecg,
                (int)std::lround(rxF.p_peak), umk.q_begin, b.r_peak_ch[c],
                umk.s_end, umk.t_end, m_sampleRate);
        }

        // ECG marker positions, indices aligned with ECG_MARKERS.
        double ecgAuto[3][kNumEcgMarkers], ecgUser[3][kNumEcgMarkers];
        for (int c = 0; c < 3; ++c) {
            // Glyphs from THIS alignment's detections; bars assembled across
            // all four alignments and expressed in this one's frame, because
            // no single alignment's marker set holds a whole beat any more.
            const auto aa = b.autoFor(anchor);
            const tbank::BankMarkerSet umk = b.userMarks(c, 0, anchor);
            const std::vector<double>& ecgA = b.chFor(c, anchor).ecgTemplate_raw;
            const FeatureMarks::ReactiveEcg rxA = FeatureMarks::reactive_ecg(
                ecgA, (int)std::lround(aa.p_begin[c]), (int)std::lround(aa.q_begin[c]),
                (int)std::lround(aa.s_end[c]), (int)std::lround(aa.t_end[c]), m_sampleRate);
            const FeatureMarks::ReactiveEcg rxU = FeatureMarks::reactive_ecg( ecgA, umk.p_begin, umk.q_begin, umk.s_end, umk.t_end, m_sampleRate);

            ecgAuto[c][0] = aa.p_begin[c];
            ecgAuto[c][1] = rxA.p_peak;          // reactive glyph, detector brackets
            ecgAuto[c][2] = aa.q_begin[c];
            ecgAuto[c][3] = ftAuto[c].q_idx;
            ecgAuto[c][4] = aa.r_peak[c];
            ecgAuto[c][5] = ftAuto[c].s_idx;
            ecgAuto[c][6] = aa.s_end[c];
            ecgAuto[c][7] = aa.t_end[c];
            ecgUser[c][0] = umk.p_begin;
            ecgUser[c][1] = rxU.p_peak;          // reactive glyph, operator brackets
            ecgUser[c][2] = umk.q_begin;
            ecgUser[c][3] = ftUser[c].q_idx;
            ecgUser[c][4] = b.r_peak_ch[c];
            ecgUser[c][5] = ftUser[c].s_idx;
            ecgUser[c][6] = umk.s_end;
            ecgUser[c][7] = umk.t_end;
        }

        // Pulse marker positions, order matching the *_MARKERS arrays.
        // t50/t80 are reactive: bracketed by onset/peak/end, never stored. The
        // autodetect column brackets with the *_auto bars and the user column
        // with the user bars, both through the same FeatureMarks::reactive_ppg
        // the on-screen glyph uses -- so what is plotted is what is exported.
        const FeatureMarks::ReactivePpg rxPpgAuto = FeatureMarks::reactive_ppg(
            b.ppgTemplate, b.ppg_onset_auto, b.ppg_peak_auto, b.ppg_dicrotic_auto, b.ppg_end_auto);
        const FeatureMarks::ReactivePpg rxPpgUser = FeatureMarks::reactive_ppg(
            b.ppgTemplate, b.ppg_onset, b.ppg_peak, b.ppg_dicrotic, b.ppg_end);
        // double, not int: t50/t80 are interpolated crossings and the stored
        // fields promote without loss. Rounding happens once, in emitLoc.
        const double ppgAuto[kNumPpgMarkers] = {
            (double)b.ppg_onset_auto, rxPpgAuto.t50, (double)b.ppg_peak_auto,
            (double)b.ppg_dicrotic_auto, (double)b.ppg_peak2_auto, rxPpgAuto.t80,
            (double)b.ppg_end_auto };
        const double ppgUser[kNumPpgMarkers] = {
            (double)b.ppg_onset, rxPpgUser.t50, (double)b.ppg_peak,
            (double)b.ppg_dicrotic, (double)b.ppg_peak2, rxPpgUser.t80,
            (double)b.ppg_end };
        const int abpAuto[kNumArterialMarkers] = { b.abp_onset_auto, b.abp_peak_auto,
            b.abp_dicrotic_auto, b.abp_peak2_auto, b.abp_end_auto };
        const int abpUser[kNumArterialMarkers] = { b.abp_onset, b.abp_peak,
            b.abp_dicrotic, b.abp_peak2, b.abp_end };
        const int artAuto[kNumArterialMarkers] = { b.art_onset_auto, b.art_peak_auto,
            b.art_dicrotic_auto, b.art_peak2_auto, b.art_end_auto };
        const int artUser[kNumArterialMarkers] = { b.art_onset, b.art_peak,
            b.art_dicrotic, b.art_peak2, b.art_end };
        const int artpAuto[kNumArterialMarkers] = { b.art_pulm_onset_auto, b.art_pulm_peak_auto,
            b.art_pulm_dicrotic_auto, b.art_pulm_peak2_auto, b.art_pulm_end_auto };
        const int artpUser[kNumArterialMarkers] = { b.art_pulm_onset, b.art_pulm_peak,
            b.art_pulm_dicrotic, b.art_pulm_peak2, b.art_pulm_end };

        // Derived T peak for the autodetect glyph group, bracketed by the AUTO J-point and T-end. The J-point is the left bracket because no T-onset
        double tPeakAutoGlyph[3];
        for (int gc = 0; gc < 3; ++gc)
            tPeakAutoGlyph[gc] = FeatureMarks::compute_t_peak(
                chs[gc]->ecgTemplate_raw,
                b.s_end_auto_ch[gc], b.t_end_auto_ch[gc]);

        auto emitVal = [&](const std::vector<double>& v, int j) {
            f << ',';
            if (j >= 0 && j < (int)v.size() && !std::isnan(v[j])) f << v[j];
            };
        // Emit ",1" if the row is this marker's row, ",<blank>" otherwise.
        //
        // ROUNDING BOUNDARY, and an unavoidable one: this column is a one-hot
        // flag per SAMPLE ROW, so a landmark at 104.37 has to be attributed to
        // a row and the format has no fractional representation. Rounded to the
        // nearest row explicitly, in one place, rather than truncated by
        // implicit conversions at each call site. The millisecond columns
        // elsewhere in this file carry the fraction.
        auto emitLoc = [&](double markerIdx, int row) {
            f << ',';
            if (markerIdx >= 0.0 && (int)std::lround(markerIdx) == row) f << '1';
            };
        for (int row = 0; row < hiRow; ++row) {
            f << m_subjectId.toStdString() << ',' << bi << ',' << (row * toMs);
            for (int k = 0; k < num_chans; ++k) {
                emitVal(*src[k].raw, row);
                emitVal(norm[k], row);
                emitVal(*src[k].rawIqr, row);
                emitVal(normIqr[k], row);
            }
            // ECG: per channel, per marker, auto then user (r_peak auto only).
            for (int c = 0; c < 3; ++c)
                for (int k = 0; k < kNumEcgMarkers; ++k) {
                    emitLoc(ecgAuto[c][k], row);
                    if (k != kRPeak) emitLoc(ecgUser[c][k], row);
                }
            // Pulse groups, auto and user interleaved per marker.
            auto emitPair = [&](const auto& a, const auto& u, int n) {
                for (int k = 0; k < n; ++k) { emitLoc(a[k], row); emitLoc(u[k], row); }
                };
            emitPair(ppgAuto, ppgUser, kNumPpgMarkers);
            emitPair(abpAuto, abpUser, kNumArterialMarkers);
            emitPair(artAuto, artUser, kNumArterialMarkers);
            emitPair(artpAuto, artpUser, kNumArterialMarkers);
            // Autodetect glyphs: per ECG channel, then the pulse table.
            for (int gc = 0; gc < 3; ++gc) {
                emitLoc(b.p_peak_auto_ch[gc], row);
                emitLoc(b.q_begin_auto_ch[gc], row);
                emitLoc(b.r_peak_auto_ch[gc], row);
                emitLoc(tPeakAutoGlyph[gc], row);
            }
            for (const auto& gl : ppg_and_artpulse_automated_markers)
                emitLoc(b.*gl.idx, row);
            f << '\n';
        }
    }

    // Serialize, suffix the value columns with the pass tag, and write this
    // pass's sidecar. Row-key columns (file_id/bin_num/x_ms) are never
    // suffixed, so the sidecars line up on the merge in save_bin_and_csv.
    std::string content = f.str();
    const size_t nl = content.find('\n');
    if (nl == std::string::npos) {
        fprintf(stderr, "[tmplcsv] malformed content (no newline) for %s\n",
            path.toStdString().c_str());
        return;
    }
    const std::string header = content.substr(0, nl);
    const std::string body = content.substr(nl);   // includes the leading '\n'
    // From the ARGUMENT, not from window state: this function is called once
    // per alignment inside a single save.
    const std::string suffix = std::string("_") + anchor_view::label(anchor);

    const QString sidecar = outDir.filePath(
        m_subjectId + "_template" + suffix.c_str() + ".csv");
    std::ofstream out(sidecar.toStdString(), std::ios::trunc);
    if (!out) {
        fprintf(stderr, "[tmplcsv] cannot open sidecar %s\n",
            sidecar.toStdString().c_str());
        return;
    }
    out << suffixValueColumns(header, suffix) << body;
    out.close();
    fprintf(stderr, "[tmplcsv] wrote sidecar %s\n", sidecar.toStdString().c_str());
}

void TemplateViewerWindow::applyMarkerVisibility() {
    for (auto* pw : m_allPlots) {
        pw->setShowEcgMarkers(m_showEcgMarkers);
        pw->setShowPpgMarkers(m_showPpgMarkers);
        pw->setShowPpgDerivMarkers(m_showPpgDerivMarkers);
        pw->setShowAbpMarkers(m_showAbpMarkers);
        pw->setShowArtMarkers(m_showArtMarkers);
        pw->setShowArtPulmMarkers(m_showArtPulmMarkers);
        pw->setShowEcgTrace(m_showEcgTrace);
        pw->setShowPpgTrace(m_showPpgTrace);
        pw->setShowAbpTrace(m_showAbpTrace);
        pw->setShowArtTrace(m_showArtTrace);
        pw->setShowArtPulmTrace(m_showArtPulmTrace);
    }
}

// ========================================================================
// Marker movement
// ========================================================================

// ============================================================================
// The single seeding path from a TemplateBin into a plot widget.
//
// Every draggable bar and every frozen autodetect column is written here, from
// one bin, in one pass. There is no second list anywhere: showPage() calls this
// on build and refreshBinMarkers() calls it on every later change, so a bar and
// its glyph can never be seeded from different places (which is how the P-onset
// bar ended up unseeded while its glyph was drawn).
//
// setAuto() MUST stay last: it performs the glyph capture, and the frozen
// snapshot reads the R bar.
// ============================================================================
void TemplateViewerWindow::onClassConfirmRequested(int binIndex, int leadIndex,
    int templateIdx, int annotationCode)
{
    if (binIndex < 0 || binIndex >= (int)m_bins.size()) return;
    TemplateBin& b = m_bins[binIndex];
    if (leadIndex < 0 || leadIndex >= 3) return;

    // KEYED ON THE SLOT, because there is one partition now. These banks come
    // from jbank::projectToChannel, which walks the joint groups in order, so
    // template i is group i on all three leads AND on PPG -- the slot is the
    // shared handle.
    //
    // IT USED TO BE KEYED ON A BEAT, because three independently built banks
    // disagreed about which slot held which morphology. That is now not just
    // unnecessary but wrong: the beat index it passed was a member of the
    // clicked lead's bank, and the projection fills `members` with CHANNEL
    // LOCAL ROW indices. Row 12 of CH1 and row 12 of CH2 are different
    // heartbeats, so the lookup labeled whichever unrelated morphology happened
    // to occupy that row number on the other leads. It also never touched the
    // pulse bank, so a confirmed PVC left its PPG cohort unlabeled.
    const tbank::TemplateBank& bank = b.ecg_bank[leadIndex];
    if (templateIdx < 0 || templateIdx >= bank.size()) return;
    const auto& members = bank.templates[templateIdx].members;
    if (members.empty()) {
        // A seeded slot 0 with no assigned beats has no beat to key on. That is
        // a real state (an empty bin), not an error, and silently doing nothing
        // would look like a dead menu item.
        statusBar()->showMessage(
            tr("No beats assigned to this template yet - nothing to confirm."),
            4000);
        return;
    }

    const uint8_t code = static_cast<uint8_t>(annotationCode);

    const tbank::PropagationResult pr =
        tbank::propagateLabelBySlot(b.ecg_bank, b.ppg_bank, templateIdx, code);

    // Say what happened, per channel. A confirmation that reached one channel
    // and not the others means that beat was unscorable in those channels --
    // worth knowing at the moment of the click, not later from a CSV.
    QStringList reached;
    for (int c = 0; c < 3; ++c)
        if (pr.labeled_template[c] >= 0)
            reached << tr("CH%1 slot %2 (subtype %3)")
            .arg(c + 1).arg(pr.labeled_template[c]).arg(pr.subtype[c]);

    if (pr.ppg_labeled_template >= 0)
        reached << tr("PPG slot %1 (subtype %2)")
        .arg(pr.ppg_labeled_template).arg(pr.ppg_subtype);

    if (reached.isEmpty()) {
        statusBar()->showMessage(
            tr("Confirmation did not reach any bank - slot %1 is empty "
                "on every channel.").arg(templateIdx), 5000);
        return;
    }
    statusBar()->showMessage(
        tr("Confirmed: %1 beats relabeled across %2")
        .arg(pr.beats_relabeled).arg(reached.join(", ")), 6000);

    // The label changes marking eligibility -- a template confirmed as PVC
    // stops wanting landmarks, one confirmed as normal starts -- and
    // markingSlotsForBin() reads exactly that. So the page is rebuilt rather than
    // repainted, and the column count can legitimately change under the
    // operator's hands. Snapshot first, as page navigation does, or the marker
    // edits made on this page are lost to the rebuild.
    captureCurrentPage();
    // Eligibility changed, so column counts changed, so page boundaries moved.
    buildPages();
    showPage();
}

// Bars for a bank template column, from the template's own BankMarkerSet.
//
// Seeded LAZILY, on first display, rather than during the build: the anchor is
// an interactive choice, so a template needs a marker set per anchor and the
// build has no idea which ones the operator will visit. Seeding on demand also
// means an existing archive gains bars without being regenerated.
//
// ECG ONLY. A bank template is a per-channel ECG median; it has no PPG or
// arterial waveform of its own, so those markers stay at -1 and their bars are
// simply absent -- which is honest, and better than showing the bin's PPG bars
// against an ECG-only column.
void TemplateViewerWindow::applyBankTemplateToWidget(BinPlotWidget* pw,
    TemplateBin& b, int channel, int templateIdx)
{
    if (channel < 0 || channel > 2) return;
    tbank::TemplateBank& bank = b.ecg_bank[channel];
    if (templateIdx < 0 || templateIdx >= bank.size()) return;
    tbank::BankTemplate& tp = bank.templates[templateIdx];
    if (tp.tmpl.empty()) return;

    // Retained only to key marks() below for the *display* fetch; every
    // alignment is seeded above and the bars are assembled from all four.
    // SEED EVERY ALIGNMENT, not just the one this window was opened for.
    //
    // SEEDED ON CONTENT, NOT ON KEY PRESENCE. hasDetectedMarks tests for a
    // VALUE; testing markers_by_anchor.find(tag) == end() does not work,
    // because marks() is operator[] and so any read of this template's marker
    // set -- from the subsequent-move propagation loop, or anywhere else that
    // touches a bin it is not displaying -- inserts an empty all -1 entry. That
    // test then reported "already seeded", auto-detection never ran, the
    // template drew with no bars at all, and serialization persisted the empty
    // set so a reload did not recover. Slot 0 was unaffected because it draws
    // from the bin's marker set, which seed_all always populates -- which is
    // exactly why it only ever showed up on the non-_A columns. With four
    // alignments per slot there are now four chances to trip it, so the guard
    // matters more, not less.
    //
    // All four are seeded from this slot's ONE R-aligned waveform: raw_anchors
    // is per (bin, channel), not per (bin, channel, slot), so a sub-template
    // has no aligned averages of its own to detect on, and its close-up falls
    // back to the bin's alignment template via chFor. Per-slot alignments would
    // mean teaching alignTemplatesFromCache to walk the banks -- a separate job.
    //
    // r_col comes from the template when it has one, otherwise from the bin's
    // channel: every beat in the bank was aligned on the same R column by
    // construction, so the bin's value is correct rather than a guess when the
    // template's own field was never filled.
    for (AnchorType a4 : anchor_view::kAllAnchors) {
        const int tag4 = static_cast<int>(a4);
        if (tp.hasDetectedMarks(tag4)) continue;
        const int rc4 = (tp.r_col >= 0) ? tp.r_col : b.r_peak_ch[channel];
        FeatureMarks::seed_bank_template(tp.tmpl, rc4, m_sampleRate, a4, tp.marks(tag4));
    }
    b.syncReactiveGlyphs(channel, templateIdx, m_sampleRate);
    // ASSEMBLED, NOT FETCHED. Each bar lives in its owning alignment's set;
    // userMarks pulls all four and translates them into the R frame the panel
    // draws in. Fetching one alignment's set here is what limited a session to
    // editing one alignment's bars.
    const tbank::BankMarkerSet mk = b.userMarks(channel, templateIdx, AnchorType::R_PEAK);
    // Bin first, for the arterial markers (a bank slot has no ABP/ART waveform
    // of its own). The seven ECG bars are overridden below; the PULSE bars and
    // glyphs are overridden here, from this slot's own waveform.
    //
    // The old comment claimed "there is ONE PPG template per bin, so every
    // column must show the SAME pulse markers." That stopped being true when
    // ppg_bank arrived: showPage draws ppg_bank.templates[templateIdx], so the
    // bin's marks were indices into a different pulse. That is why the foot sat
    // nowhere near a minimum.
    applyBinToWidget(pw, b);

    if (templateIdx >= 0 && templateIdx < b.ppg_bank.size()) {
        tbank::BankTemplate& ps = b.ppg_bank.templates[templateIdx];
        if (!ps.tmpl.empty()) {
            if (!ps.hasDetectedPulseMarks())
                FeatureMarks::seed_pulse_bank_template(ps.tmpl, m_ppgRateHz,
                    ps.pulse_marks);
            const tbank::BankPulseMarkerSet& pm = ps.pulse_marks;
            pw->setMarker(BinPlotWidget::PpgOnset, pm.onset);
            pw->setMarker(BinPlotWidget::PpgPeak, pm.peak);
            pw->setMarker(BinPlotWidget::PpgDicrotic, pm.dicrotic);
            pw->setMarker(BinPlotWidget::PpgPeak2, pm.peak2);
            pw->setMarker(BinPlotWidget::PpgEnd, pm.end);
            pw->setMarker(BinPlotWidget::PpgT50, pm.t50);
            pw->setMarker(BinPlotWidget::PpgT80, pm.t80);
            pw->overridePulseGlyphs(pm);
        }
    }

    pw->setMarker(BinPlotWidget::EcgPBegin, mk.p_begin);
    pw->setMarker(BinPlotWidget::EcgPPeak, mk.p_peak);
    pw->setMarker(BinPlotWidget::EcgQBegin, mk.q_begin);
    pw->setMarker(BinPlotWidget::EcgRPeak,
        (tp.r_col >= 0) ? tp.r_col : b.r_peak_ch[channel]);
    pw->setMarker(BinPlotWidget::EcgSEnd, mk.s_end);
    pw->setMarker(BinPlotWidget::EcgTEnd, mk.t_end);
}

void TemplateViewerWindow::applyBinToWidget(BinPlotWidget* pw, const TemplateBin& b) {
    const int c = pw->leadIndex();
    // ALL FOUR BARS, EACH FROM ITS OWN ALIGNMENT, translated into the R frame
    // this widget draws in. The grid is R-aligned on every panel; only the
    // close-up switches waveform.
    const tbank::BankMarkerSet mk = b.userMarks(c, 0, AnchorType::R_PEAK);

    // Bank members are their own COLUMNS now (see the (bin, template) expansion
    // in the layout loop), so nothing is overlaid here -- drawing them again as
    // dashed traces under slot 0 would duplicate what the neighbouring columns
    // already show.

    // ---- draggable bars ----------------------------------------------------
    pw->setMarker(BinPlotWidget::EcgPBegin, mk.p_begin);
    pw->setMarker(BinPlotWidget::EcgPPeak, mk.p_peak);
    pw->setMarker(BinPlotWidget::EcgQBegin, mk.q_begin);
    pw->setMarker(BinPlotWidget::EcgRPeak, b.r_peak_ch[c]);   // auto-only, no bar drawn
    pw->setMarker(BinPlotWidget::EcgSEnd, mk.s_end);
    pw->setMarker(BinPlotWidget::EcgTEnd, mk.t_end);

    pw->setMarker(BinPlotWidget::PpgOnset, b.ppg_onset);
    pw->setMarker(BinPlotWidget::PpgPeak, b.ppg_peak);
    pw->setMarker(BinPlotWidget::PpgDicrotic, b.ppg_dicrotic);
    pw->setMarker(BinPlotWidget::PpgPeak2, b.ppg_peak2);
    pw->setMarker(BinPlotWidget::PpgEnd, b.ppg_end);
    // T50/T80 are neither drawn nor draggable -- they're reactive glyphs now.
    // Kept in sync anyway so the enum entries never hold a stale position.
    pw->setMarker(BinPlotWidget::PpgT50, b.ppg_t50);
    pw->setMarker(BinPlotWidget::PpgT80, b.ppg_t80);

    pw->setMarker(BinPlotWidget::AbpOnset, b.abp_onset);
    pw->setMarker(BinPlotWidget::AbpPeak, b.abp_peak);
    pw->setMarker(BinPlotWidget::AbpDicrotic, b.abp_dicrotic);
    pw->setMarker(BinPlotWidget::AbpPeak2, b.abp_peak2);
    pw->setMarker(BinPlotWidget::AbpEnd, b.abp_end);
    pw->setMarker(BinPlotWidget::ArtOnset, b.art_onset);
    pw->setMarker(BinPlotWidget::ArtPeak, b.art_peak);
    pw->setMarker(BinPlotWidget::ArtDicrotic, b.art_dicrotic);
    pw->setMarker(BinPlotWidget::ArtPeak2, b.art_peak2);
    pw->setMarker(BinPlotWidget::ArtEnd, b.art_end);
    pw->setMarker(BinPlotWidget::ArtPulmOnset, b.art_pulm_onset);
    pw->setMarker(BinPlotWidget::ArtPulmPeak, b.art_pulm_peak);
    pw->setMarker(BinPlotWidget::ArtPulmDicrotic, b.art_pulm_dicrotic);
    pw->setMarker(BinPlotWidget::ArtPulmPeak2, b.art_pulm_peak2);
    pw->setMarker(BinPlotWidget::ArtPulmEnd, b.art_pulm_end);
    pw->setAuto(b);   // last: captures the glyph snapshot and repaints
}

void TemplateViewerWindow::refreshBinMarkers(int binIdx) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    for (int li = 0; li < (int)m_pageGlobalIdx.size(); ++li) {
        if (m_pageGlobalIdx[li] != binIdx) continue;

        // Only slot-0 columns carry the bin's MarkerSet, so bank columns are
        // skipped here for the same reason showPage() skips them.
        //
        // And no `break`: a bin now owns one column per bank member, and the
        // original loop stopped after the first match -- which would silently
        // refresh only the sinus column of a polymorphic bin.
        if (li < (int)m_pageTemplateIdx.size() && m_pageTemplateIdx[li] != 0)
            continue;
        for (auto* pw : m_binPlots[li]) applyBinToWidget(pw, m_bins[binIdx]);
    }
}

// The bank-column counterpart of refreshBinMarkers. Repaints only the columns
// showing (binIdx, templateIdx), because a bank slot's bars live in that
// template's own BankMarkerSet and applyBinToWidget would draw the BIN's set
// over them -- the sinus landmarks, on an ectopic waveform.
void TemplateViewerWindow::refreshBankMarkers(int binIdx, int templateIdx) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    if (templateIdx <= 0) { refreshBinMarkers(binIdx); return; }

    for (int li = 0; li < (int)m_pageGlobalIdx.size(); ++li) {
        if (m_pageGlobalIdx[li] != binIdx) continue;
        if (li >= (int)m_pageTemplateIdx.size()) continue;
        if (m_pageTemplateIdx[li] != templateIdx) continue;
        // Each widget in the column is one lead, and a bank is per lead, so the
        // widget's own leadIndex() selects the bank to draw from -- not the
        // dragged lead, which would paint lead 1's bars onto lead 2's panel.
        for (auto* pw : m_binPlots[li])
            applyBankTemplateToWidget(pw, m_bins[binIdx],
                pw->leadIndex(), templateIdx);
    }
}

// Slot-aware entry point. A drag on a sub-template column must land in that
// template's BankMarkerSet, not in the bin's -- the bin's set describes sinus,
// and writing a PVC's Q-onset into it would corrupt the sinus landmarks for
// every other panel showing the same bin.
void TemplateViewerWindow::onMarkerMovedOnTemplate(int binIdx, int leadIdx,
    int templateIdx, int marker, int newIdx)
{
    if (templateIdx == 0) { onMarkerMoved(binIdx, leadIdx, marker, newIdx); return; }
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    if (leadIdx < 0 || leadIdx > 2) return;
    if (!BinPlotWidget::markerIsEcg(marker)) return;   // bank templates are ECG only

    tbank::TemplateBank& bank = m_bins[binIdx].ecg_bank[leadIdx];
    if (templateIdx < 0 || templateIdx >= bank.size()) return;

    // THE BAR CHOOSES THE ALIGNMENT, not the pass. This was
    // static_cast<int>(currentAnchor()) -- whichever alignment the window had
    // been opened for -- so a drag could only ever land in one alignment's set
    // and the other three were reachable only by restarting the session.
    const AnchorType ownerAnchor = anchor_view::anchorFor(marker);
    const int anchor = static_cast<int>(ownerAnchor);

    // Read/write one landmark of one bank slot. Taken as lambdas over the bank
    // rather than written inline, because the subsequent-move loop below has to
    // perform the same access on a DIFFERENT bin's bank and the two must not
    // drift -- the same shape onMarkerMoved uses for the bin's marker set.
    // SLOT IS A PARAMETER, not templateIdx captured. The propagation loop below
    // targets the later bin's BAND-MATCHING slot, which is generally a different
    // number from the dragged one -- so these two hardcoding templateIdx was
    // half of why propagation only ever touched same-numbered columns.
    auto bankGet = [&](tbank::TemplateBank& bk, int slot) -> int {
        tbank::BankMarkerSet& m = bk.templates[slot].marks(anchor);
        switch (marker) {
        case BinPlotWidget::EcgPBegin: return m.p_begin;
        case BinPlotWidget::EcgPPeak:  return m.p_peak;
        case BinPlotWidget::EcgQBegin: return m.q_begin;
        case BinPlotWidget::EcgSEnd:   return m.s_end;
        case BinPlotWidget::EcgTEnd:   return m.t_end;
        }
        return -1;
        };
    auto bankSet = [&](tbank::TemplateBank& bk, int slot, int v) {
        tbank::BankMarkerSet& m = bk.templates[slot].marks(anchor);
        switch (marker) {
        case BinPlotWidget::EcgPBegin: m.p_begin = v; break;
        case BinPlotWidget::EcgPPeak:  m.p_peak = v; break;
        case BinPlotWidget::EcgQBegin: m.q_begin = v; break;
        case BinPlotWidget::EcgSEnd:   m.s_end = v; break;
        case BinPlotWidget::EcgTEnd:   m.t_end = v; break;
            // EcgRPeak has no bar and is not draggable: it is the alignment
            // anchor, and moving it would invalidate every other landmark's
            // coordinate.
        }
        };

    // Anything not in the six cases above is not a bank-draggable landmark, so
    // return before writing rather than silently no-op'ing inside bankSet and
    // then running the propagation loop for a marker nobody moved.
    // Bars only. P peak used to be in this list, which let a drag write a
    // stored P peak that the reactive recompute then contradicted; it is a
    // glyph and markerAtX does not hand it out in the first place.
    if (!anchor_view::isBar(marker)) return;

    const int oldIdx = bankGet(bank, templateIdx);
    if (m_dragStartIdx < 0) m_dragStartIdx = oldIdx;   // first move of this drag
    bankSet(bank, templateIdx, newIdx);

    // ---- subsequent-move, BY SLOT INDEX ---------------------------------
    // Bank columns used to ignore m_moveMode entirely: this handler wrote one
    // field and returned, so only slot-0 drags ever propagated and a drag on a
    // B or C column moved exactly one bar. It now follows the same two modes
    // the bin's markers do.
    //
    // SUBSEQUENT MEANS EVERY PANEL AFTER THIS ONE ON THE PAGE, ALL LETTERS.
    //
    // The set is exactly the columns on screen, in the order they are drawn:
    // m_pageGlobalIdx and m_pageTemplateIdx are parallel and ARE the page. So a
    // drag on _A moves _E, and it moves nothing the operator cannot see.
    //
    // TWO THINGS THIS FIXES.
    //
    // It used to propagate by slot INDEX -- bk.templates[templateIdx] in each
    // later bin -- so a drag on _A touched nothing but other _A columns, and a
    // drag on a deep slot touched nothing at all in any bin whose bank was
    // shorter than that index. That is the "only moves ones of the same letter"
    // symptom.
    //
    // And it used to walk EVERY BIN IN THE RECORD, not the page. Bins beyond the
    // page have no panel, the operator never saw them, and markers_by_anchor IS
    // serialized -- so a single drag persisted landmarks across the whole record
    // and exported them as though someone had placed each one. Walking the page
    // is what makes "subsequent" mean what the button says.
    //
    // Iterating the page also removes the need to re-test what a column is:
    // markingSlotsForBin already applied every eligibility rule to build it, so
    // there is no second copy of those tests here to fall out of step.
    if (m_moveMode != MoveMode::Individual && oldIdx >= 0) {
        const int nDragged =
            static_cast<int>(bank.templates[templateIdx].tmpl.size());

        // BOTH MODES ARE FRACTIONAL. A raw sample delta is the wrong unit across
        // bins: each bin's frame is its own width (alignment frames on that bin's
        // longest RR), so the same +12 samples is a different physical shift in
        // every one, and a much smaller fraction of the beat on a bin holding a
        // pause. A drag of "10% later" should stay 10% later.
        //
        //   Delta : cur + deltaFrac * (n - 1)   -- a relative nudge
        //   Raw   : posFrac  * (n - 1)          -- an absolute position
        // FROM THE DRAG START, not the previous mouse-move. Per-event, each
        // shift was rounded on its own, so a one-sample step scaled to a
        // fraction and lround'd to ZERO -- dragging smoothly moved the other
        // bars not at all, a fast sweep jumped them, and out-and-back never
        // returned. One total shift, rounded once.
        const double deltaFrac = (nDragged > 1)
            ? static_cast<double>(newIdx - m_dragStartIdx) / (nDragged - 1) : 0.0;
        const double posFrac = (nDragged > 1)
            ? static_cast<double>(newIdx) / (nDragged - 1) : 0.0;

        // One body, applied to a (bin, slot) pair drawn from the page.
        auto applyTo = [&](int binI, int slot) {
            if (binI < 0 || binI >= (int)m_bins.size()) return;
            if (m_bins[binI].bad_r_ch[leadIdx]) return;
            tbank::TemplateBank& bk = m_bins[binI].ecg_bank[leadIdx];
            if (slot < 0 || slot >= bk.size()) return;
            if (bk.templates[slot].tmpl.empty()) return;
            const int n = static_cast<int>(bk.templates[slot].tmpl.size());
            if (n <= 0) return;

            // SEED BEFORE READING. bankGet goes through marks(), which
            // inserts, so reading an unseeded template here is what used to
            // suppress its auto-detection forever. Seeding first makes the read
            // harmless and gives the propagated delta a real bar to move.
            //
            // hasDetectedMarks is on the TEMPLATE, not the bank: marker sets
            // are per template, so there is no bank-wide answer to "are the
            // marks detected".
            tbank::BankTemplate& tgt = bk.templates[slot];
            if (!tgt.hasDetectedMarks(anchor)) {
                const int rc = (tgt.r_col >= 0)
                    ? tgt.r_col : m_bins[binI].r_peak_ch[leadIdx];
                FeatureMarks::seed_bank_template(tgt.tmpl, rc, m_sampleRate,
                    ownerAnchor, tgt.marks(anchor));
            }
            const int cur = bankGet(bk, slot);
            if (cur < 0) return;   // this landmark was not found on this one

            // This panel's own start. Keyed on (bin, slot): applyTo is a lambda
            // over those two, and the page column index is not in scope here.
            // After the cur < 0 guard, so an absent landmark never stores -1.
            const int key = binI * 64 + slot;
            if (!original_location_of_bar.count(key))
                original_location_of_bar[key] = cur;

            const int target = (m_moveMode == MoveMode::SubsequentDelta)
                ? originFor(key, cur) + (int)std::lround(deltaFrac * (n - 1))
                : (int)std::lround(posFrac * (n - 1));
            if (target < 0 || target > n - 1) return;
            bankSet(bk, slot, target);
            };

        // Walk the PAGE from the column after the dragged one to the end. The
        // two arrays are parallel, so index li identifies a panel; the dragged
        // panel is found by matching both halves, because a bin contributes
        // several columns and a slot number recurs across bins.
        int dragCol = -1;
        for (int li = 0; li < (int)m_pageGlobalIdx.size()
            && li < (int)m_pageTemplateIdx.size(); ++li) {
            if (m_pageGlobalIdx[li] == binIdx
                && m_pageTemplateIdx[li] == templateIdx) {
                dragCol = li; break;
            }
        }
        for (int li = dragCol + 1; li < (int)m_pageGlobalIdx.size()
            && li < (int)m_pageTemplateIdx.size(); ++li)
            applyTo(m_pageGlobalIdx[li], m_pageTemplateIdx[li]);

        // Refresh exactly the panels that were written, by the same page walk.
        // This used to refresh (gi > binIdx, templateIdx) -- the wrong slot on
        // the wrong set of bins, now that propagation is per column.
        for (int li = dragCol + 1; li < (int)m_pageGlobalIdx.size()
            && li < (int)m_pageTemplateIdx.size(); ++li)
            refreshBankMarkers(m_pageGlobalIdx[li], m_pageTemplateIdx[li]);
    }
}

void TemplateViewerWindow::onMarkerMoved(int binIdx, int leadIdx,
    int marker, int newIdx)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    TemplateBin& b = m_bins[binIdx];

    if (BinPlotWidget::markerIsEcg(marker)) {
        if (leadIdx < 0 || leadIdx > 2) return;
        // SLOT IS A PARAMETER NOW, AND THAT IS THE WHOLE POINT.
        //
        // This function used to have no slot dimension at all: it read and wrote
        // the bin's own MarkerSet, which held slot 0 and nothing else. It is also
        // the function a slot-0 drag runs (onMarkerMovedOnTemplate routes
        // templateIdx == 0 here). So a drag on _A could not move _E however the
        // propagation loop was written -- the storage it reached had nowhere to
        // put a sub-template's landmark. One home fixed that; this is where it
        // pays off.
        // THE BAR CHOOSES THE ALIGNMENT. This was currentAnchor() -- the
        // pass -- so a drag could only ever land in one alignment's set, and
        // the other three were reachable only by finishing the window and
        // reopening it on the next alignment.
        //
        // Reads and writes are in the OWNER'S frame; the widget hands us a
        // column in the R frame, so the caller converts. P peak is gone from
        // these switches: it is a glyph, recomputed by syncReactiveGlyphs from
        // the two bars that bracket it, so a drag must not write it.
        if (!anchor_view::isBar(marker) && marker != BinPlotWidget::EcgRPeak) return;
        const AnchorType owner = anchor_view::anchorFor(marker);
        auto ecgGet = [&](TemplateBin& tb, int slot) -> int {
            if (marker == BinPlotWidget::EcgRPeak) return tb.r_peak_ch[leadIdx];
            const tbank::BankMarkerSet& m =
                tb.slotMarks(leadIdx, slot, owner);
            switch (marker) {
            case BinPlotWidget::EcgPBegin: return m.p_begin;
            case BinPlotWidget::EcgQBegin: return m.q_begin;
            case BinPlotWidget::EcgSEnd:   return m.s_end;
            case BinPlotWidget::EcgTEnd:   return m.t_end;
            }
            return -1;
            };
        auto ecgSet = [&](TemplateBin& tb, int slot, int v) {
            // R is the alignment anchor: auto-only, and per BIN rather than per
            // slot, because every template in the bank is aligned on it.
            if (marker == BinPlotWidget::EcgRPeak) { tb.r_peak_ch[leadIdx] = v; return; }
            tbank::BankMarkerSet& m = tb.slotMarks(leadIdx, slot, owner);
            switch (marker) {
            case BinPlotWidget::EcgPBegin: m.p_begin = v; break;
            case BinPlotWidget::EcgQBegin: m.q_begin = v; break;
            case BinPlotWidget::EcgSEnd:   m.s_end = v; break;
            case BinPlotWidget::EcgTEnd:   m.t_end = v; break;
            }
            };

        // newIdx is a column the WIDGET measured, and the grid is R-aligned on
        // every panel. The bar is stored in its own alignment's frame, so the
        // column converts here -- a no-op for the J point (owner is R) and a
        // sample or two for the other three.
        const int stored = newIdx + b.frameShift(leadIdx, AnchorType::R_PEAK, owner);
        const int oldIdx = ecgGet(b, 0);
        if (m_dragStartIdx < 0) m_dragStartIdx = oldIdx;
        ecgSet(b, 0, stored);
        // P peak follows its brackets. Recomputed for every alignment, so the
        // value written to the bin and the CSV is the one drawn on screen.
        b.syncReactiveGlyphs(leadIdx, 0, m_sampleRate);
        // Track the drag: the touched store follows the bar to its final
        // position so confirmedIndex reflects where the operator left it.
        if (newIdx >= 0)
            m_touchedMarks[touchKey(binIdx, leadIdx, marker)] = newIdx;
        if (m_moveMode != MoveMode::Individual && oldIdx >= 0) {
            // ONE PROPAGATION LOOP, OVER THE PAGE'S OWN COLUMNS.
            //
            // Subsequent means every panel after this one in drawing order, all
            // letters. m_pageGlobalIdx and m_pageTemplateIdx are parallel and ARE
            // the page, so _A moves _E -- and nothing off-page is touched, which
            // matters because landmarks are serialized and this loop used to run
            // to the end of the record, persisting marks into bins the operator
            // never saw and exporting them as though someone had placed each one.
            //
            // BOTH MODES ARE FRACTIONAL. A raw sample delta is the wrong unit
            // across bins: each bin's frame is its own width, so the same +12
            // samples is a different physical shift in every one, and a much
            // smaller fraction of the beat on a bin holding a pause. A drag of
            // "10% later" stays 10% later.
            ChannelTemplateData* bchs[3] = { &b.ch1, &b.ch2, &b.ch3 };
            const int nDragged = (int)bchs[leadIdx]->ecgTemplate_raw.size();
            // From the drag start -- see the note in the bank handler.
            const double deltaFrac = (nDragged > 1)
                ? double(newIdx - m_dragStartIdx) / (nDragged - 1) : 0.0;
            const double posFrac = (nDragged > 1)
                ? double(newIdx) / (nDragged - 1) : 0.0;

            // Slot 0's length is the bin's own template; a deeper slot's is its
            // bank template's, which can be shorter.
            auto lenOf = [&](int binI, int slot) -> int {
                if (slot == 0) {
                    ChannelTemplateData* chs[3] = {
                        &m_bins[binI].ch1, &m_bins[binI].ch2, &m_bins[binI].ch3
                    };
                    return (int)chs[leadIdx]->ecgTemplate_raw.size();
                }
                const tbank::TemplateBank& bk = m_bins[binI].ecg_bank[leadIdx];
                if (slot >= (int)bk.templates.size()) return 0;
                return (int)bk.templates[slot].tmpl.size();
                };

            // The dragged panel, matched on BOTH halves: a slot number recurs
            // across bins and a bin contributes several columns, so matching
            // either alone lands on the wrong one.
            int dragCol = -1;
            for (int li = 0; li < (int)m_pageGlobalIdx.size()
                && li < (int)m_pageTemplateIdx.size(); ++li) {
                if (m_pageGlobalIdx[li] == binIdx
                    && m_pageTemplateIdx[li] == 0) {
                    dragCol = li; break;
                }
            }

            for (int li = dragCol + 1; li < (int)m_pageGlobalIdx.size()
                && li < (int)m_pageTemplateIdx.size(); ++li) {
                const int gi = m_pageGlobalIdx[li];
                const int slot = m_pageTemplateIdx[li];
                if (gi < 0 || gi >= (int)m_bins.size()) continue;
                if (m_bins[gi].bad_r_ch[leadIdx]) continue;
                const int n = lenOf(gi, slot);
                if (n <= 0) continue;
                const int cur = ecgGet(m_bins[gi], slot);
                if (cur < 0) continue;

                // This panel's own start. `li` IS the page column here.
                if (!original_location_of_bar.count(li))
                    original_location_of_bar[li] = cur;

                const int target = (m_moveMode == MoveMode::SubsequentDelta)
                    ? originFor(li, cur) + (int)std::lround(deltaFrac * (n - 1))
                    : (int)std::lround(posFrac * (n - 1));
                if (target < 0 || target > n - 1) continue;
                ecgSet(m_bins[gi], slot, target);
            }
            // Refresh exactly the columns that were written.
            for (int li = dragCol + 1; li < (int)m_pageGlobalIdx.size()
                && li < (int)m_pageTemplateIdx.size(); ++li)
                refreshBankMarkers(m_pageGlobalIdx[li], m_pageTemplateIdx[li]);
        }
        // B2: an edit to this landmark updates its focus view. For the
        // J-point (S-end) this refreshes BOTH the QRS and JT panels; other
        // landmarks refresh their single panel. refreshFocus() handles the
        // routing.
        // Slot 0: onMarkerMovedOnTemplate routes every other slot to its own
        // handler, so a drag reaching here is by definition on the seed column.
        refreshFocus(binIdx, leadIdx, 0, marker, newIdx);
        return;
    }

    if (BinPlotWidget::markerIsPpg(marker)) {
        auto ppgGet = [&](TemplateBin& tb) -> int {
            switch (marker) {
            case BinPlotWidget::PpgOnset:    return tb.ppg_onset;
            case BinPlotWidget::PpgPeak:     return tb.ppg_peak;
            case BinPlotWidget::PpgDicrotic: return tb.ppg_dicrotic;
            case BinPlotWidget::PpgPeak2:    return tb.ppg_peak2;
            case BinPlotWidget::PpgEnd:      return tb.ppg_end;
            }
            return -1;
            };
        auto ppgSet = [&](TemplateBin& tb, int v) {
            switch (marker) {
            case BinPlotWidget::PpgOnset:    tb.ppg_onset = v; break;
            case BinPlotWidget::PpgPeak:     tb.ppg_peak = v; break;
            case BinPlotWidget::PpgDicrotic: tb.ppg_dicrotic = v; break;
            case BinPlotWidget::PpgPeak2:    tb.ppg_peak2 = v; break;
            case BinPlotWidget::PpgEnd:      tb.ppg_end = v; break;
            }
            };

        const int oldIdx = ppgGet(b);
        ppgSet(b, newIdx);
        refreshBinMarkers(binIdx);
        const int delta = newIdx - oldIdx;

        if (m_moveMode != MoveMode::Individual && oldIdx >= 0) {
            const int nDragged = (int)b.ppgTemplate.size();
            const double pct = (nDragged > 1) ? double(newIdx) / (nDragged - 1) : 0.0;

            for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
                if (m_bins[i].bad_ppg != 0) continue;
                const int cur = ppgGet(m_bins[i]);
                if (cur < 0) continue;

                const int rawLen = (int)m_bins[i].ppgTemplate.size();
                const int ecgClip = ecgClipLenFor(m_bins[i]);
                const int n = (ecgClip > 0) ? std::min(rawLen, ecgClip) : rawLen;
                if (n <= 0) continue;

                const int target = (m_moveMode == MoveMode::SubsequentDelta)
                    ? cur + delta
                    : (n > 1 ? (int)std::lround(pct * (n - 1)) : 0);
                if (target < 0 || target > n - 1) continue;
                ppgSet(m_bins[i], target);
            }
            for (int li = 0; li < (int)m_pageGlobalIdx.size(); ++li) {
                int gi = m_pageGlobalIdx[li];
                if (gi > binIdx && m_bins[gi].bad_ppg == 0) refreshBinMarkers(gi);
            }
        }
        return;
    }

    // ---- Arterial markers (ABP / ART / ART_PULM) --------------------------
    // Shared across leads like PPG. Route to the right channel's fields and
    // propagate to subsequent bins when Move-Subsequent is on.
    if (BinPlotWidget::markerIsArterial(marker)) {
        // Select the channel's field pointers, issue flag, and trace by group.
        auto assign = [&](TemplateBin& tb, int mk, int val) {
            switch (mk) {
            case BinPlotWidget::AbpOnset:    tb.abp_onset = val; break;
            case BinPlotWidget::AbpPeak:     tb.abp_peak = val; break;
            case BinPlotWidget::AbpDicrotic: tb.abp_dicrotic = val; break;
            case BinPlotWidget::AbpPeak2:       tb.abp_peak2 = val; break;
            case BinPlotWidget::AbpEnd:      tb.abp_end = val; break;
            case BinPlotWidget::ArtOnset:    tb.art_onset = val; break;
            case BinPlotWidget::ArtPeak:     tb.art_peak = val; break;
            case BinPlotWidget::ArtDicrotic: tb.art_dicrotic = val; break;
            case BinPlotWidget::ArtPeak2:       tb.art_peak2 = val; break;
            case BinPlotWidget::ArtEnd:      tb.art_end = val; break;
            case BinPlotWidget::ArtPulmOnset:    tb.art_pulm_onset = val; break;
            case BinPlotWidget::ArtPulmPeak:     tb.art_pulm_peak = val; break;
            case BinPlotWidget::ArtPulmDicrotic: tb.art_pulm_dicrotic = val; break;
            case BinPlotWidget::ArtPulmPeak2:       tb.art_pulm_peak2 = val; break;
            case BinPlotWidget::ArtPulmEnd:      tb.art_pulm_end = val; break;
            }
            };
        auto artGet = [&](TemplateBin& tb, int mk) -> int {
            switch (mk) {
            case BinPlotWidget::AbpOnset:    return tb.abp_onset;
            case BinPlotWidget::AbpPeak:     return tb.abp_peak;
            case BinPlotWidget::AbpDicrotic: return tb.abp_dicrotic;
            case BinPlotWidget::AbpPeak2:    return tb.abp_peak2;
            case BinPlotWidget::AbpEnd:      return tb.abp_end;
            case BinPlotWidget::ArtOnset:    return tb.art_onset;
            case BinPlotWidget::ArtPeak:     return tb.art_peak;
            case BinPlotWidget::ArtDicrotic: return tb.art_dicrotic;
            case BinPlotWidget::ArtPeak2:    return tb.art_peak2;
            case BinPlotWidget::ArtEnd:      return tb.art_end;
            case BinPlotWidget::ArtPulmOnset:    return tb.art_pulm_onset;
            case BinPlotWidget::ArtPulmPeak:     return tb.art_pulm_peak;
            case BinPlotWidget::ArtPulmDicrotic: return tb.art_pulm_dicrotic;
            case BinPlotWidget::ArtPulmPeak2:    return tb.art_pulm_peak2;
            case BinPlotWidget::ArtPulmEnd:      return tb.art_pulm_end;
            }
            return -1;
            };
        auto channelTrace = [&](TemplateBin& tb, int mk,
            const std::vector<double>*& tr, uint8_t*& iss) {
                if (BinPlotWidget::markerIsAbp(mk)) { tr = &tb.abpTemplate; iss = &tb.abp_issue; }
                else if (BinPlotWidget::markerIsArt(mk)) { tr = &tb.artTemplate; iss = &tb.art_issue; }
                else { tr = &tb.artPulmTemplate; iss = &tb.art_pulm_issue; }
            };

        const int oldIdx = artGet(b, marker);
        assign(b, marker, newIdx);
        refreshBinMarkers(binIdx);
        const int delta = newIdx - oldIdx;

        if (m_moveMode != MoveMode::Individual && oldIdx >= 0) {
            const std::vector<double>* trDragged = nullptr; uint8_t* issDragged = nullptr;
            channelTrace(b, marker, trDragged, issDragged);
            const int nDragged = trDragged ? (int)trDragged->size() : 0;
            const double pct = (nDragged > 1) ? double(newIdx) / (nDragged - 1) : 0.0;

            for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
                const std::vector<double>* tr = nullptr; uint8_t* iss = nullptr;
                channelTrace(m_bins[i], marker, tr, iss);
                if (!iss || *iss != 0) continue;
                if (!tr) continue;
                const int cur = artGet(m_bins[i], marker);
                if (cur < 0) continue;

                const int rawLen = (int)tr->size();
                const int ecgClip = ecgClipLenFor(m_bins[i]);
                const int n = (ecgClip > 0) ? std::min(rawLen, ecgClip) : rawLen;
                if (n <= 0) continue;

                const int target = (m_moveMode == MoveMode::SubsequentDelta)
                    ? cur + delta
                    : (n > 1 ? (int)std::lround(pct * (n - 1)) : 0);
                if (target < 0 || target > n - 1) continue;
                assign(m_bins[i], marker, target);
            }
            for (int li = 0; li < (int)m_pageGlobalIdx.size(); ++li) {
                int gi = m_pageGlobalIdx[li];
                if (gi > binIdx) refreshBinMarkers(gi);
            }
        }
        return;
    }
}
int TemplateViewerWindow::originFor(int col, int cur) const {
    auto it = original_location_of_bar.find(col);
    return (it != original_location_of_bar.end()) ? it->second : cur;
}

// A drag begins. BinPlotWidget emits this once from mousePressEvent, so it is
// the one place that can define "where the drag started".
//
// BOTH HALVES of the shift are anchored to that instant. The dragged bar has
// moved a PERCENTAGE of its plot; applying that total percentage to a bar
// already moved by earlier events in the same drag counts it twice. So the
// numerator (how far the dragged bar has come) and the base (where each
// propagated bar began) both date from here and hold for the whole drag --
// which is what makes dragging back to the start put every bar back exactly.
void TemplateViewerWindow::onMarkerDragStarted(int, int, int) {
    m_dragStartIdx = -1;              // dragged bar's start; set on the first move
    original_location_of_bar.clear(); // per-panel starts; filled lazily below
}
// B2 focus mode --------------------------------------------------------------

void TemplateViewerWindow::onLandmarkSelected(int binIdx, int leadIdx,
    int templateIdx, int marker, int col)
{
    // Focus activation = the operator clicked this bar. Record it as "touched"
    // at position col; logBoundaryTrainingAtSave reads this to fill
    // confirmedIndex for the matching landmark.
    if (binIdx >= 0 && leadIdx >= 0 && col >= 0)
        m_touchedMarks[touchKey(binIdx, leadIdx, marker)] = col;
    refreshFocus(binIdx, leadIdx, templateIdx, marker, col);
}

// Rebuild the focus panel(s) for one landmark from the current bin/lead's
// anchored-average stats. Reads mean/sd/n straight from the template the
// viewer already holds:
//   mean = ecgTemplate_raw
//   sd   = ecg_template_raw_iqr  (holds STD, ddof=1 -- despite the _iqr name)
//   n    = ch{1,2,3}_n_beats_raw (per-bin, not per-channel-struct)
// The J-point (S-end) is shared by the QRS and JT views, so selecting/editing
// it refreshes BOTH panels; every other landmark refreshes its own single
// panel.
// Move the stretch between the second panel and the trailing spacer so the
// panels always sit on the SAME third-height grid.
//
//   split=false -> panel 1/3, panel(hidden) 0, spacer 2/3
//   split=true  -> panel 1/3, panel        1/3, spacer 1/3
//
// A hidden widget contributes no stretch, so without moving it into the spacer
// a lone visible panel would expand to fill half the dock -- and the same
// landmark would then be drawn at one scale on its own and another right after
// the J point had been selected. The spacer holds the leftover.
void TemplateViewerWindow::setFocusSplit(bool split) {
    if (!m_focusLay || !zoomed_in_section_bottom) return;
    zoomed_in_section_bottom->setVisible(split);
    m_focusLay->setStretch(1, split ? 1 : 0);   // second panel
    m_focusLay->setStretch(2, split ? 1 : 2);   // trailing spacer absorbs the rest
}
// Radio group that pins the focus panel to one alignment.
//
// findChild rather than ui->r_align_button on purpose: a button that is not in
// the .ui yet simply leaves that option unavailable instead of failing the
// build, so they can be added one at a time.
void TemplateViewerWindow::wireAlignButtons() {
    struct Btn { const char* name; bool force; AnchorType a; };
    static const Btn kBtns[] = {
        { "r_align_button",         true,  AnchorType::R_PEAK  },
        { "p_align_button",         true,  AnchorType::P_ONSET },
        { "q_align_button",         true,  AnchorType::Q_ONSET },
        { "j_point_align_button",         true,  AnchorType::J_POINT },
        { "automatic_align_button", false, AnchorType::R_PEAK  },
    };
    for (const Btn& b : kBtns) {
        QRadioButton* rb = findChild<QRadioButton*>(QString::fromLatin1(b.name));
        if (!rb) continue;
        const bool force = b.force;
        const AnchorType a = b.a;
        connect(rb, &QRadioButton::toggled, this, [this, force, a](bool on) {
            if (!on) return;                  // only the newly-checked one acts
            m_forceAlign = force;
            m_forcedAlign = a;
            if (m_lastFocusMarker >= 0)
                refreshFocus(m_lastFocusBinIdx, m_lastFocusLeadIdx,
                    m_lastFocusTemplateIdx, m_lastFocusMarker, m_lastFocusCol);
            });
        if (rb->isChecked()) { m_forceAlign = force; m_forcedAlign = a; }
    }
}

void TemplateViewerWindow::refreshFocus(int binIdx, int leadIdx,
    int templateIdx, int marker, int col)
{
    if (!zoomed_in_section_top) return;   // panels not created (nothing to do)
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    TemplateBin& b = m_bins[binIdx];

    // ---- PPG and ARTERIAL landmarks (B2 focus extended to all channels) --
    // These ride their own template (ppgTemplate / abpTemplate / artTemplate
    // / artPulmTemplate) with the matching per-sample std (*_iqr, ddof=1) and
    // the shared pulse beat count (ppg_n_beats -- all pulse channels derive
    // from the same foot-anchored beat set). Routed to the QRS panel as the
    // single focus view for pulse channels (they have no QRS/JT split).
    if (!BinPlotWidget::markerIsEcg(marker)) {
        const std::vector<double>* meanRaw = nullptr;
        const std::vector<double>* iqrRaw = nullptr;
        int pulseChan = -1;   // index into m_pulseGlobalRef: PPG=0,ABP=1,ART=2,ART_PULM=3
        // -1 until a branch sets it. The arterial channels have no bank, so
        // they keep the bin-wide pulse beat count, which for them IS the whole
        // population -- there is no per-group arterial cohort to get wrong.
        int nPulseBeats = -1;
        int footIdx = -1;     // this channel's foot/onset column (perfusion-index baseline)
        QString chLabel;
        if (BinPlotWidget::markerIsPpg(marker)) {
            // THE GROUP'S PULSE, NOT THE BIN'S. ppg_bank slot i is group i, on
            // the same axis as the bin's pulse template. This path read
            // b.ppgTemplate / b.ppg_template_iqr / b.ppg_n_beats
            // unconditionally, so clicking a pulse landmark on ANY column
            // showed the bin's mean, its bin-wide spread and its bin-wide beat
            // count -- the same defect the main panel had, one layer over.
            //
            // NO FALLBACK. A group with no pulse cohort has no focus view:
            // both panels are cleared and the function returns. Showing the
            // bin's waveform there would be a measurement attributed to beats
            // that are not in this template.
            const tbank::BankTemplate* ps =
                (templateIdx >= 0 && templateIdx < b.ppg_bank.size())
                ? &b.ppg_bank.templates[templateIdx] : nullptr;
            if (!ps || ps->tmpl.empty() || ps->memberCount() <= 0) {
                setFocusSplit(false);
                if (zoomed_in_section_top) zoomed_in_section_top->clearFocus();
                if (zoomed_in_section_bottom) zoomed_in_section_bottom->clearFocus();
                return;
            }
            meanRaw = &ps->tmpl;           iqrRaw = &ps->tmpl_iqr;
            pulseChan = 0; footIdx = b.ppg_onset;  chLabel = "PPG";
            nPulseBeats = ps->memberCount();
        }
        else if (BinPlotWidget::markerIsAbp(marker)) {
            meanRaw = &b.abpTemplate;      iqrRaw = &b.abpTemplate_iqr;      pulseChan = 1; footIdx = b.abp_onset;      chLabel = "ABP";
        }
        else if (BinPlotWidget::markerIsArt(marker)) {
            meanRaw = &b.artTemplate;      iqrRaw = &b.artTemplate_iqr;      pulseChan = 2; footIdx = b.art_onset;      chLabel = "ART";
        }
        else if (BinPlotWidget::markerIsArtPulm(marker)) {
            meanRaw = &b.artPulmTemplate;  iqrRaw = &b.artPulmTemplate_iqr;  pulseChan = 3; footIdx = b.art_pulm_onset; chLabel = "ART_PULM";
        }
        if (!meanRaw || meanRaw->empty()) return;

        // Pulse channels are NOT normalized by a plain scalar (that was the
        // bug -- it left the trace flat). The displayed trace uses a per-
        // sample PERFUSION-INDEX transform relative to the pulse's own foot,
        // then /ref (normalize_ppg_or_similar -> normalize_pulse_trace, see
        // the main plot ~line 508). The mean MUST use that same transform.
        const std::vector<double> mean = normalize_ppg_or_similar(*meanRaw, footIdx, pulseChan);
        // The *_iqr is ALREADY in perfusion-index space (local_ratio_iqr at
        // build time), so it only needs the scalar /ref -- NOT the perfusion
        // transform again (main plot ~line 492). It's a true IQR (Q3-Q1), so
        // convert to an SD estimate (IQR/1.349) for the 95% CI.
        const double ref = (pulseChan >= 0 && pulseChan < 4) ? m_pulseGlobalRef[pulseChan] : std::nan("");
        std::vector<double> sd = normalize_features::scale_array_by_ref(*iqrRaw, ref);
        for (double& s : sd) if (!std::isnan(s)) s /= 1.349;

        const int nBeats = (nPulseBeats >= 0)
            ? nPulseBeats : static_cast<int>(b.ppg_n_beats);

        auto pulseLabel = [](int m) -> QString {
            switch (m) {
            case BinPlotWidget::PpgOnset:    return QStringLiteral("Foot");
            case BinPlotWidget::PpgT50:      return QStringLiteral("T50");
                // SWAPPED. PpgPeak is the SYSTOLIC peak -- the forward-wave maximum
                // the whole pulse is anchored on -- and PpgPeak2 is the DIASTOLIC
                // peak, the reflected wave arriving after the dicrotic notch. Every
                // other reference in the tree agrees: ppg_peak2_color is commented
                // "(2nd/diastolic peak)", and feature_marks builds peak2 as "first
                // local max after the notch". Only these two labels disagreed, and
                // they disagreed with each other in a way that made the diastolic
                // bar look like a missing systolic one.
            case BinPlotWidget::PpgPeak:     return QStringLiteral("Systolic Peak");
            case BinPlotWidget::PpgDicrotic: return QStringLiteral("Dicrotic Notch");
            case BinPlotWidget::PpgPeak2:    return QStringLiteral("Diastolic Peak");
            case BinPlotWidget::PpgT80:      return QStringLiteral("T80");
            case BinPlotWidget::PpgEnd:      return QStringLiteral("End");
            case BinPlotWidget::AbpOnset: case BinPlotWidget::ArtOnset: case BinPlotWidget::ArtPulmOnset:       return QStringLiteral("onset");
            case BinPlotWidget::AbpPeak: case BinPlotWidget::ArtPeak: case BinPlotWidget::ArtPulmPeak:          return QStringLiteral("peak");
            case BinPlotWidget::AbpDicrotic: case BinPlotWidget::ArtDicrotic: case BinPlotWidget::ArtPulmDicrotic: return QStringLiteral("dicrotic");
            case BinPlotWidget::AbpPeak2: case BinPlotWidget::ArtPeak2: case BinPlotWidget::ArtPulmPeak2:        return QStringLiteral("peak2");
            case BinPlotWidget::AbpEnd: case BinPlotWidget::ArtEnd: case BinPlotWidget::ArtPulmEnd:             return QStringLiteral("end");
            }
            return QStringLiteral("landmark");
            };

        // Pulse channels have no alignment dimension: they are foot-anchored,
        // once, and raw_anchors is ECG-only. One panel, no suffix.
        // Pulse landmarks bound one part of the wave: top third only.
        setFocusSplit(false);
        if (zoomed_in_section_bottom) zoomed_in_section_bottom->clearFocus();
        if (zoomed_in_section_top)
            zoomed_in_section_top->setFocus(mean, sd, nBeats, col,
                chLabel + " " + pulseLabel(marker));
        return;
    }

    // ---- ECG landmarks ---------------------------------------------------
    if (leadIdx < 0 || leadIdx > 2) return;

    // ================= THE ALIGNMENT SWITCH =============================
    //
    // This is the whole feature. The grid draws the R-aligned average on every
    // panel; this panel draws the alignment the clicked landmark is measured
    // on -- P-aligned under the P-onset bar, Q-aligned under Q-onset,
    // R-aligned under the J point, T-aligned under T-end. Same alignment the
    // bar is stored in and reported under, so what the operator is looking at
    // while placing a bar is the waveform its column is a column OF.
    //
    // A glyph has no alignment of its own (it is measured on all four), so a
    // click that resolves to one shows the R-aligned average -- the same
    // waveform the grid does, magnified. markerAtX does not hand out glyphs
    // anyway; this is just what anchorFor returns for them.
    m_lastFocusBinIdx = binIdx;
    m_lastFocusLeadIdx = leadIdx;
    m_lastFocusTemplateIdx = templateIdx;
    m_lastFocusMarker = marker;
    m_lastFocusCol = col;

    const AnchorType focusAnchor = m_forceAlign
        ? m_forcedAlign
        : anchor_view::anchorFor(marker);    // NO FALLBACK: chForStrict returns nullptr when this alignment is absent
    // from the file, and the panel is cleared rather than showing the R-aligned
    // average under this bar's header. chFor's fallback did the latter, which
    // is what made every bar look identical.
    const ChannelTemplateData* chP = b.chForStrict(leadIdx, focusAnchor);
    if (!chP) {
        setFocusSplit(false);
        if (zoomed_in_section_top) zoomed_in_section_top->clearFocus();
        if (zoomed_in_section_bottom) zoomed_in_section_bottom->clearFocus();
        fprintf(stderr, "[focus] bin=%d lead=%d marker=%d anchor=%s NOT IN FILE"
            " -- regenerate templates\n",
            binIdx, leadIdx, marker, anchor_view::label(focusAnchor));
        return;
    }
    const ChannelTemplateData& ch = *chP;


    // Templates + iqr are stored PRE-reference-division; scale by the same
    // per-lead ref the displayed ECG trace uses (main plot ~line 484). The
    // ECG *_iqr field already holds a STD (ddof=1) -- CreateEcgTemplates step
    // 7 changed it from IQR to std despite the "iqr" name -- so it feeds the
    // CI directly, NO IQR->SD conversion (unlike the pulse channels, whose
    // *_iqr is a true interquartile range).
    const double eref = m_ecgGlobalRef[leadIdx];

    // ---- WHICH WAVEFORM THIS PANEL IS SHOWING --------------------------
    // Slot 0 is the bin's chN template. Every other slot is a bank member with
    // its own median, its own spread and its own beat count, and this path used
    // ch.ecgTemplate_raw and the bin's chN_n_beats_raw regardless -- so the
    // focus view for a PVC column plotted the sinus average and reported the
    // whole bin's beat count under it.
    //
    // NO FALLBACK for a slot that exists but is empty: the panels are cleared
    // and the function returns, because a focus view of the wrong morphology is
    // worse than none.
    const std::vector<double>* meanRawEcg = &ch.ecgTemplate_raw;
    const std::vector<double>* sdRawEcg = &ch.ecg_template_raw_iqr;
    // Set when a sub-template panel had to fall back to the slot's own
    // unaligned average; changes the header so it cannot claim an alignment
    // the data is not in.
    bool slotUnaligned = false;
    const uint64_t nb[3] = { b.ch1_n_beats_raw, b.ch2_n_beats_raw, b.ch3_n_beats_raw };
    int nBeats = static_cast<int>(nb[leadIdx]);
    if (templateIdx > 0) {
        const tbank::TemplateBank& bank = b.ecg_bank[leadIdx];
        if (templateIdx >= bank.size()
            || bank.templates[templateIdx].tmpl.empty()) {
            setFocusSplit(false);
            if (zoomed_in_section_top) zoomed_in_section_top->clearFocus();
            if (zoomed_in_section_bottom) zoomed_in_section_bottom->clearFocus();
            return;
        }
        const tbank::BankTemplate& tp = bank.templates[templateIdx];
        meanRawEcg = &tp.tmpl;
        sdRawEcg = &tp.tmpl_iqr;
        nBeats = tp.cleanCount();

        // ---- THIS SLOT'S OWN ALIGNED AVERAGE ---------------------------
        //
        // THIS IS THE BUG THIS BRANCH HAD. focusAnchor was computed and
        // chFor fetched the right aligned waveform two dozen lines above --
        // and then the two assignments right here threw both away, because
        // tp.tmpl is the slot's single R-aligned average with no anchor
        // dimension at all. So on an _A column the panel showed the same trace
        // and the same spread under every bar, while the header named whichever
        // alignment the bar belonged to.
        //
        // Per-slot aligned averages now exist (v6 section, computed in
        // alignTemplatesFromCache as a reduction over the beats it already
        // aligns). nullptr means the templates file predates that section, in
        // which case the header says so rather than naming an alignment the
        // data is not in.
        if (const AnchoredBankSlot* asl =
            b.bankSlotFor(leadIdx, templateIdx, focusAnchor)) {
            meanRawEcg = &asl->tmpl;
            sdRawEcg = &asl->tmpl_iqr;
        }
        else {
            slotUnaligned = true;
        }
    }
    const std::vector<double> mean = normalize_features::scale_array_by_ref(*meanRawEcg, eref);
    const std::vector<double> sd = normalize_features::scale_array_by_ref(*sdRawEcg, eref);
    if (mean.empty()) return;


    // No `side` argument any more: "J-point (QRS)" vs "(JT)" was the
    // two-panel distinction, and the panel now names the ALIGNMENT instead,
    // which is the thing that actually differs between views.
    auto labelFor = [](int m) -> QString {
        switch (m) {
        case BinPlotWidget::EcgPBegin: return QStringLiteral("P onset");
        case BinPlotWidget::EcgPPeak:  return QStringLiteral("P peak");
        case BinPlotWidget::EcgQBegin: return QStringLiteral("Q onset");
        case BinPlotWidget::EcgRPeak:  return QStringLiteral("R peak");
        case BinPlotWidget::EcgSEnd:   return QStringLiteral("J-point");
        case BinPlotWidget::EcgTEnd:   return QStringLiteral("T end");
        }
        return QStringLiteral("landmark");
        };

    // `col` arrived in the R frame -- the widget's, because the grid is
    // R-aligned. This panel plots THIS alignment's average, so the bar has to
    // be placed in its columns. Same conversion and same direction as the drag
    // path in onMarkerMoved, so a bar dragged on the grid lands under the
    // crosshair here.
    // CLAMPED. frameShift returned 0 for the life of this code until
    // sub-sample alignment landed -- every anchor shared one r_col, so the
    // translation was a no-op. Now the anchors' r_cols genuinely differ, so
    // this is the first build where colHere can fall outside the trace, and
    // an out-of-range landmark column is not something the panel should be
    // asked to render.
    int colHere = col + b.frameShift(leadIdx, AnchorType::R_PEAK, focusAnchor);
    if (mean.empty()) return;
    if (colHere < 0) colHere = 0;
    if (colHere >= static_cast<int>(mean.size()))
        colHere = static_cast<int>(mean.size()) - 1;

    // ONE PANEL. There were two -- QRS above, JT below -- because the J point
    // had to appear in both, framed right-edge in one and left-edge in the
    // other. That split existed to compensate for a single-alignment view: the
    // J point is where the QRS view ends and the T view begins, so neither
    // panel alone could show both of its neighbourhoods. With the alignment
    // switching per bar there is one waveform per landmark and nothing to
    // reconcile.
    //
    // Onsets frame toward the LEFT edge (they start their segment), offsets
    // toward the RIGHT (they end it).
    const int bias = BinPlotWidget::markerIsBegin(marker) ? +1 : -1;

    if (zoomed_in_section_top)
    {
        // Mean SD over the same window the panel draws. On screen because the
        // panel autoscales y, so the BAND looks about the same width whatever
        // the spread actually is -- the number is the only way to compare two
        // alignments. In a P-aligned average the T region is the most smeared
        // part of the trace, so P-aligned-at-T should read larger than
        // T-aligned-at-T.
        double acc = 0.0; int cnt = 0;
        for (int k = colHere - 30; k <= colHere + 30; ++k)
            if (k >= 0 && k < (int)sd.size() && !std::isnan(sd[k])) { acc += sd[k]; ++cnt; }
        const double sdMean = cnt
            ? acc / cnt : std::numeric_limits<double>::quiet_NaN();

        // Two states only: this alignment's own template, or -- on a
        // sub-template column whose file predates the per-slot section -- the
        // slot's unaligned average, said plainly. chForStrict has already
        // returned early if the alignment itself is absent, so the label can
        // never name an alignment the data is not in.
        const QString tag = slotUnaligned
            ? QStringLiteral(" [slot avg, NOT aligned]")
            : QStringLiteral(" [%1-aligned]")
            .arg(QString::fromLatin1(anchor_view::label(focusAnchor)));

        // THE J POINT GETS BOTH PANELS. It is the one landmark bounding two
        // segments, so a single framing always hides one neighbourhood: above
        // it is framed to the RIGHT edge, where it ends the QRS; below to the
        // LEFT, where it starts the JT. Same waveform in both -- differing
        // only in which side of the landmark is shown.
        const QString head = QString("%1%2  sd=%3")
            .arg(labelFor(marker), tag, QString::number(sdMean, 'f', 5));
        if (marker == BinPlotWidget::EcgSEnd) {
            setFocusSplit(true);
            if (zoomed_in_section_top)
                zoomed_in_section_top->setFocus(mean, sd, nBeats, colHere,
                    head + QStringLiteral("  (QRS)"), 30, -1);
            if (zoomed_in_section_bottom)
                zoomed_in_section_bottom->setFocus(mean, sd, nBeats, colHere,
                    head + QStringLiteral("  (JT)"), 30, +1);
        }
        else {
            // One segment -> top third only.
            setFocusSplit(false);
            if (zoomed_in_section_bottom) zoomed_in_section_bottom->clearFocus();
            if (zoomed_in_section_top)
                zoomed_in_section_top->setFocus(mean, sd, nBeats, colHere, head, 30, bias);
        }
    }
}

// ========================================================================
// BadR / BadPPG
// ========================================================================

// ============================================================================
// QUALITY MARKS: PER PANEL, WITH SLOT 0 OWNING THE BIN
//
// A right-click marks the panel it was made on and nothing else. It used to
// record against the BIN, so with a bin occupying one panel per morphology, one
// click crossed out up to six columns -- and the extra ones appeared on the next
// page rebuild rather than at the moment of the click, which made it look like
// the click had landed somewhere else entirely.
//
// WHY THE BIN-LEVEL FLAGS SURVIVE. TemplateBin::bad_r_ch and bad_ppg are
// consumed: NormalizeFeatures skips a bin's channel on bad_r_ch when building
// the global reference, template_marking_bin_io serializes them and exports them
// per bin, and feature_marks sets them automatically. They have to keep meaning
// "this bin's lead is untrustworthy", so the SEED panel owns them and a
// sub-template panel does not touch them. Marking a 2-beat junk column must not
// exclude a whole bin from the feature reference.
// ============================================================================

// The bank slot behind a panel, for a given bin and lead. Null when the slot is
// out of range, which is a real state -- a bin whose bank did not survive the
// trip from generation has one panel and no templates.
tbank::BankTemplate* TemplateViewerWindow::slotFor(int binIdx, int leadIdx,
    int templateIdx)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return nullptr;
    if (templateIdx < 0) return nullptr;
    TemplateBin& b = m_bins[binIdx];
    tbank::TemplateBank& bk = (leadIdx >= 0 && leadIdx <= 2)
        ? b.ecg_bank[leadIdx] : b.ppg_bank;
    if (templateIdx >= bk.size()) return nullptr;
    return &bk.templates[templateIdx];
}

// Repaints exactly the panel that was clicked.
void TemplateViewerWindow::repaintPanel(int binIdx, int leadIdx, int templateIdx,
    BinPlotWidget::State st)
{
    for (size_t ci = 0; ci < m_pageGlobalIdx.size(); ++ci) {
        if (m_pageGlobalIdx[ci] != binIdx) continue;
        if (ci >= m_pageTemplateIdx.size() || m_pageTemplateIdx[ci] != templateIdx)
            continue;
        if (ci >= m_binPlots.size()) continue;
        for (auto* pw : m_binPlots[ci])
            if (pw && (leadIdx < 0 || pw->leadIndex() == leadIdx))
                pw->setState(st);
        return;
    }
}

void TemplateViewerWindow::onBadRToggled(int binIdx, int leadIdx,
    int templateIdx, bool bad) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    if (leadIdx < 0 || leadIdx > 2) return;

    if (tbank::BankTemplate* t = slotFor(binIdx, leadIdx, templateIdx))
        t->operator_state = bad ? 1u : 0u;

    // SLOT 0 ONLY writes the bin-level flag. See the header note above.
    if (templateIdx == 0) m_bins[binIdx].bad_r_ch[leadIdx] = bad;

    repaintPanel(binIdx, leadIdx, templateIdx,
        bad ? BinPlotWidget::State::BadR : BinPlotWidget::State::Good);
}

void TemplateViewerWindow::onBadPPGToggled(int binIdx, int templateIdx,
    bool bad) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;

    // The pulse verdict is recorded on the PULSE bank's slot, not on the ECG
    // lead's -- it is a statement about the pulse waveform in this panel.
    if (tbank::BankTemplate* t = slotFor(binIdx, -1, templateIdx))
        t->operator_state = bad ? 2u : 0u;

    if (templateIdx == 0) {
        m_bins[binIdx].bad_ppg = bad ? 1 : 0;
        // bad_ppg supersedes bad_r on the same bin, as before: the two are
        // alternatives in the right-click cycle, not independent flags.
        if (bad)
            for (int c = 0; c < 3; ++c) m_bins[binIdx].bad_r_ch[c] = false;
    }

    // Every lead of THIS panel: a bad pulse is not a per-lead judgement, and the
    // panel shows the same pulse trace under each lead.
    repaintPanel(binIdx, -1, templateIdx,
        bad ? BinPlotWidget::State::BadPPG : BinPlotWidget::State::Good);
}

// AnchorType -> short name for the boundary log's `anchor` column.
static const char* anchorName_boundary(AnchorType a) {
    switch (a) {
    case AnchorType::P_ONSET: return "P_ONSET";
    case AnchorType::Q_ONSET: return "Q_ONSET";
    case AnchorType::R_PEAK:  return "R_PEAK";
    case AnchorType::J_POINT: return "J_POINT";
    }
    return "UNKNOWN";
}

// Log boundary training data at save time. For EVERY template (bin/lead), log
// each landmark (Q-onset, J-point, P-onset, T-end). auto_detect comes from the
// bin's auto-detected glyph fields (*_auto_ch), which hold every landmark on
// every template regardless of anchor pass. expert_mark = the user marker for
// that landmark iff the operator moved it away from auto (else null).
void TemplateViewerWindow::logBoundaryTrainingAtSave() {
    using boundary_training::Landmark;

    struct Target { Landmark lm; };
    // S-end omitted (same feature as J-point). Each landmark's auto position is
    // read from the bin below; its expert mark (if any) from the current pass.
    const Target targets[] = {
        { Landmark::Q_ONSET  },
        { Landmark::J_POINT  },
        { Landmark::P_ONSET  },
        { Landmark::T_OFFSET },   // T-end
    };

    // auto-detected position (glyph field) for a landmark on a lead. The
    // *_auto_ch fields are double (sub-sample); rounded here since this seeds
    // an integer segment window.
    auto autoPosOf = [](const TemplateBin& tb, Landmark lm, int lead) -> int {
        switch (lm) {
        case Landmark::Q_ONSET:  return (int)std::lround(tb.q_begin_auto_ch[lead]);
        case Landmark::J_POINT:  return (int)std::lround(tb.s_end_auto_ch[lead]);   // J-point == S-end field
        case Landmark::P_ONSET:  return (int)std::lround(tb.p_begin_auto_ch[lead]);
        case Landmark::T_OFFSET: return (int)std::lround(tb.t_end_auto_ch[lead]);
        default: return -1;
        }
        };
    // Map each logged landmark to its BinPlotWidget marker id (for the touch key).
    auto markerIdOf = [](Landmark lm) -> int {
        switch (lm) {
        case Landmark::Q_ONSET:  return BinPlotWidget::EcgQBegin;
        case Landmark::J_POINT:  return BinPlotWidget::EcgSEnd;
        case Landmark::P_ONSET:  return BinPlotWidget::EcgPBegin;
        case Landmark::T_OFFSET: return BinPlotWidget::EcgTEnd;
        default: return -1;
        }
        };

    const int half = std::max(1, (int)std::lround(0.100 * m_sampleRate)); // +/-100 ms
    int written = 0, failed = 0;

    // Slicing puts R at r_col = percent_interval_preceeding_rpeak * RR, so
    // RR_samples = r_col / that fraction. Used for heart rate.
    constexpr double kPreRFrac = alignment::percent_interval_preceeding_rpeak;

    for (int i = 0; i < (int)m_bins.size(); ++i) {
        const TemplateBin& b = m_bins[i];
        const std::vector<double>* chs[3] = {
            &b.ch1.ecgTemplate_raw, &b.ch2.ecgTemplate_raw, &b.ch3.ecgTemplate_raw };
        const int rcol[3] = { b.ch1.r_col_raw, b.ch2.r_col_raw, b.ch3.r_col_raw };
        for (int lead = 0; lead < 3; ++lead) {
            const std::vector<double>& sig = *chs[lead];

            // Heart rate (bpm) = 60000 / RR(ms); RR from this lead's r_col.
            double heartRate = 0.0;
            if (rcol[lead] > 0 && m_sampleRate > 0.0) {
                const double rrSamples = rcol[lead] / kPreRFrac;
                const double rrMs = rrSamples / m_sampleRate * 1000.0;
                if (rrMs > 0.0) heartRate = 60000.0 / rrMs;
            }
            // QRS duration (ms) = distance between q_begin and s_end glyphs.
            double qrsMs = 0.0;
            {
                const int q = (int)std::lround(b.q_begin_auto_ch[lead]);
                const int s = (int)std::lround(b.s_end_auto_ch[lead]);
                if (q >= 0 && s >= 0 && m_sampleRate > 0.0)
                    qrsMs = std::abs(s - q) / m_sampleRate * 1000.0;
            }

            for (const Target& tgt : targets) {
                const int autoPos = autoPosOf(b, tgt.lm, lead);
                if (autoPos < 0 || autoPos >= (int)sig.size()) continue;

                const int lo = std::max(0, autoPos - half);
                const int hi = std::min((int)sig.size(), autoPos + half);
                if (hi - lo < 2) continue;

                // confirmedIndex = the operator's clicked position (segment-
                // relative) if this landmark was touched (focus activated on
                // its bar), else -1 => blank. The row is logged either way.
                int confirmed = -1;
                const int mid = markerIdOf(tgt.lm);
                auto it = m_touchedMarks.find(touchKey(i, lead, mid));
                if (it != m_touchedMarks.end()) confirmed = it->second - lo;

                boundary_training::BoundaryTrainingRecord rec;
                rec.segment.assign(sig.begin() + lo, sig.begin() + hi);
                rec.confirmedIndex = confirmed;
                // fit from anchor_fit: fit-and-select on this landmark's window.
                const anchor_fit::FitResult fit = anchor_fit::selectAnchorModel(sig, lo, hi - 1);
                rec.fitType = fit.type;
                rec.fitRSS = fit.rss;
                rec.individualID = m_subjectId.toStdString();
                rec.bbb = false;                 // left blank for now
                rec.heartRate = heartRate;
                rec.qrsDurationMs = qrsMs;
                if (boundary_training::logBoundary(
                    m_boundaryLog,
                    // PER LANDMARK, not per session. This column records which
                    // alignment the expert mark was placed on, and that is now
                    // a property of the landmark: the P onset was placed on the
                    // P-aligned average, the J point on the R-aligned one. It
                    // was constant for a whole file when a file was one pass.
                    anchorName_boundary(anchor_view::anchorFor(markerIdOf(tgt.lm))),
                    tgt.lm, rec)) ++written; else ++failed;
            }
        }
    }
    std::fprintf(stderr,
        "[boundary_log] save: dir='%s' wrote=%d failed=%d\n",
        m_boundaryLog.dir.c_str(), written, failed);
}

void TemplateViewerWindow::save_bin_and_csv() {
    captureCurrentPage();   // snapshot the page being left on Finish

    QDir binDir(m_markingPath);
    QDir csvDir(m_markingPath);
    QDir vcgDir(m_vcgOutputPath);
    const QString canonicalBin = QDir(m_markingPath).filePath(m_subjectId + "_template_markings.bin");

    try {
        // STRAIGHT TO THE CANONICAL NAME. There used to be a .partial written
        // on every pass and promoted on the last one, because the window was
        // torn down and rebuilt between alignments and the pulse markers (which
        // have no alignment dimension) had to survive the gap. One session, one
        // write: the file already carries every alignment's marker set, keyed
        // by anchor tag, so there is nothing in flight to stage.
        writeTemplateMarkingsBin(canonicalBin.toStdString(), m_bins);
        std::cout << "Saved: " << canonicalBin.toStdString() << "\n";
        logBoundaryTrainingAtSave();

        // ECG markings CSV: each pass writes its OWN per-anchor sidecar
        // (<id>_template_markings.<ANCHOR>.csv) with suffixed columns. No
        // growing zip per pass; all sidecars are merged into the canonical
        // <id>_template_markings.csv once, at the final pass.
        // VCG loop features: one row per bin, <id>_vcg.csv, beside this CSV.
        // Written every pass -- it is a standalone file, not a per-anchor
        // sidecar that needs merging, and the markers it measures against
        // change on every pass, so the latest write is the one that matters.
        {
            std::vector<vcg_avg::BinFeatures> vcgRows;
            vcgRows.reserve(m_bins.size());
            for (int vi = 0; vi < (int)m_bins.size(); ++vi) {
                // Prefer the operator's markers. On the R pass (and on any bin
                // not yet marked for this anchor) there are none, so fall back
                // to the autodetected ones rather than emitting an empty row --
                // the features are still measurable, and the source is
                // recorded so a row is never ambiguous about where its
                // boundaries came from.
                // R frame: global intervals compare landmarks ACROSS leads on
                // one axis, and the reference lines drawn from them are drawn
                // on the R-aligned panels. leadMarkersFor's USER path assembles
                // the bars from all four alignments into this frame.
                auto vgi = global_intervals::computeGlobalIntervals(
                    m_bins[vi], AnchorType::R_PEAK, m_sampleRate,
                    global_intervals::MarkerSource::USER);
                bool fromAuto = false;
                if (!vgi.valid) {
                    vgi = global_intervals::computeGlobalIntervals(
                        m_bins[vi], AnchorType::R_PEAK, m_sampleRate,
                        global_intervals::MarkerSource::AUTO);
                    fromAuto = vgi.valid;
                }
                auto vrow = vcg_avg::analyzeBinFromTemplates(
                    vi, m_bins[vi], vgi, m_sampleRate);
                if (vrow.valid && fromAuto) vrow.note = "auto markers";
                vcgRows.push_back(vrow);
            }
            const bool vok = vcg_avg::writeVcgCsv(
                vcgDir.absolutePath().toStdString(),
                m_subjectId.toStdString(), vcgRows);
            std::cout << (vok ? "Saved: " : "FAILED: ")
                << vcgDir.absolutePath().toStdString() << "/"
                << m_subjectId.toStdString() << "_vcg.csv\n";
        }

        // ---- FOUR SIDECARS, ONE SAVE ----------------------------------
        //
        // Each sidecar holds one alignment's glyph columns and the one bar it
        // owns, with its value columns suffixed _R / _P / _Q / _T; the merge
        // joins them on (file_id, bin_index) exactly as it did when the four
        // were written one window apart. Pulse last, un-suffixed, once -- it
        // has no alignment dimension.
        //
        // The glob-discovery this replaces existed only to find sidecars left
        // by PREVIOUS sessions of the same subject. There is one session now,
        // so the order is kAllAnchors, and a stale sidecar from an older build
        // can no longer sneak into the merge.
        const QString csvPath = csvDir.absolutePath() + "/" + m_subjectId + "_template_markings.csv";
        const QString pulseSidecar = csvDir.absolutePath() + "/"
            + m_subjectId + "_template_markings_PULSE.csv";

        // Write one sidecar: CSV to a temp file, then re-emit with every value
        // column suffixed. An empty suffix skips the renaming (pulse).
        auto writeSidecar = [&](const QString& path, const QString& suffix,
            AnchorType anchor, MarkingsCsvSection section) {
                const QString tmpPath = path + ".tmp";
                writeTemplateMarkingsCsv(tmpPath.toStdString(), m_bins,
                    m_subjectId.toStdString(), m_sampleRate, anchor, section);
                std::ifstream tin(tmpPath.toStdString());
                std::stringstream tbuf; tbuf << tin.rdbuf();
                tin.close();
                QFile::remove(tmpPath);
                const std::string tcontent = tbuf.str();
                const size_t tnl = tcontent.find('\n');
                if (tnl == std::string::npos)
                    throw std::runtime_error("markings CSV writer produced malformed content (no newline)");
                const std::string tHeader = tcontent.substr(0, tnl);
                const std::string tBody = tcontent.substr(tnl);
                std::ofstream out(path.toStdString(), std::ios::trunc);
                if (!out) throw std::runtime_error("cannot open for write: " + path.toStdString());
                out << (suffix.isEmpty()
                    ? tHeader
                    : suffixValueColumns(tHeader, suffix.toStdString()))
                    << tBody;
                out.close();
                std::cout << "Saved sidecar: " << path.toStdString() << "\n";
            };

        std::vector<std::string> order;
        for (AnchorType a : anchor_view::kAllAnchors) {
            const QString suffix = QString("_") + anchor_view::label(a);
            const QString sidecar = csvDir.absolutePath() + "/"
                + m_subjectId + "_template_markings" + suffix + ".csv";
            writeSidecar(sidecar, suffix, a, MarkingsCsvSection::EcgOnly);
            order.push_back(sidecar.toStdString());
        }
        writeSidecar(pulseSidecar, QString(), AnchorType::R_PEAK,
            MarkingsCsvSection::PulseOnly);
        order.push_back(pulseSidecar.toStdString());
        std::cout.flush();

        {
            if (!mergeSidecarCsvs(csvPath.toStdString(), order))
                throw std::runtime_error("could not merge markings sidecars into " + csvPath.toStdString());
            for (const auto& p : order) QFile::remove(QString::fromStdString(p));
            std::cout << "Merged markings CSV: " << csvPath.toStdString() << "\n";

        }
    }
    catch (const std::exception& e) {
        QMessageBox::critical(this, "Save failed",
            QString("Could not write markings for %1:\n\n%2\n\n"
                "If the CSV is open in Excel, close it and try again.")
            .arg(m_subjectId, e.what()));
        return;   // don't emit finished(); let the user retry
    }

    // Aligned-template CSV: one sidecar per alignment, holding that
    // alignment's own averages, merged into the canonical <id>_template.csv.
    // Same restructuring as the markings sidecars above, same reason -- the
    // glob was reading files left by earlier passes of the same subject, and
    // there are no earlier passes now.
    {
        QDir alignedDir(m_templateDir);
        const QString canonical = alignedDir.filePath(m_subjectId + "_template.csv");
        std::vector<std::string> order;
        for (AnchorType a : anchor_view::kAllAnchors) {
            writeAlignedTemplateCsv(a);
            order.push_back(alignedDir.filePath(
                m_subjectId + "_template_" + anchor_view::label(a) + ".csv")
                .toStdString());
        }
        if (!order.empty() && mergeSidecarCsvs(canonical.toStdString(), order)) {
            for (const auto& p : order) QFile::remove(QString::fromStdString(p));
            std::cout << "Merged aligned-template CSV: " << canonical.toStdString() << "\n";
        }
    }

    // ONE SAVE, ONE FINISH. There is no cycle to advance: requestQAlignReload
    // and the pass counter are gone, and the four alignments were all written
    // above.
    emit finished();
}