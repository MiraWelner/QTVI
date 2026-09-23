// ========================================================================
// template_viewer_page.cpp
// The bin grid: paging, column layout, panel construction and refresh.
//
// One of five units; TemplateViewerWindow is declared in template_viewer.hpp.
// ========================================================================

#include "template_viewer.hpp"


// ---- SHARED PULSE FIELD TABLE -------------------------------------------
// Two units walk it: applyBinCommonToWidget pushes the arterial entries into a
// panel (_page) and movePpgMarker resolves a dragged marker back to its field
// (_marks). `inline` in a named namespace, with identical definitions in both,
// so the two TUs name one object.
namespace tv_detail {

    struct PulseField { int marker; double TemplateBin::* field; };

    inline constexpr PulseField kPulseFields[] = {
        { BinPlotWidget::PpgOnset,           &TemplateBin::ppg_onset },
        { BinPlotWidget::PpgPeak,            &TemplateBin::ppg_peak },
        { BinPlotWidget::PpgDicrotic,        &TemplateBin::ppg_dicrotic },
        { BinPlotWidget::PpgPeak2,           &TemplateBin::ppg_peak2 },
        { BinPlotWidget::PpgEnd,             &TemplateBin::ppg_end },
        // T50/T80 are reactive glyphs: neither drawn from here nor draggable.
        // Pushed anyway so the enum entries never hold a stale position.
        { BinPlotWidget::PpgT50,             &TemplateBin::ppg_t50 },
        { BinPlotWidget::PpgT80,             &TemplateBin::ppg_t80 },
        { BinPlotWidget::AbpOnset,           &TemplateBin::abp_onset },
        { BinPlotWidget::AbpPeak,            &TemplateBin::abp_peak },
        { BinPlotWidget::AbpDicrotic,        &TemplateBin::abp_dicrotic },
        { BinPlotWidget::AbpPeak2,           &TemplateBin::abp_peak2 },
        { BinPlotWidget::AbpEnd,             &TemplateBin::abp_end },
        { BinPlotWidget::ArtOnset,           &TemplateBin::art_onset },
        { BinPlotWidget::ArtPeak,            &TemplateBin::art_peak },
        { BinPlotWidget::ArtDicrotic,        &TemplateBin::art_dicrotic },
        { BinPlotWidget::ArtPeak2,           &TemplateBin::art_peak2 },
        { BinPlotWidget::ArtEnd,             &TemplateBin::art_end },
        { BinPlotWidget::ArtPulmOnset,       &TemplateBin::art_pulm_onset },
        { BinPlotWidget::ArtPulmPeak,        &TemplateBin::art_pulm_peak },
        { BinPlotWidget::ArtPulmDicrotic,    &TemplateBin::art_pulm_dicrotic },
        { BinPlotWidget::ArtPulmPeak2,       &TemplateBin::art_pulm_peak2 },
        { BinPlotWidget::ArtPulmEnd,         &TemplateBin::art_pulm_end },
    };

    inline double* pulseField(TemplateBin& tb, int marker) {
        for (const PulseField& f : kPulseFields)
            if (f.marker == marker) return &(tb.*f.field);
        return nullptr;
    }

}   // namespace tv_detail

using tv_detail::PulseField;
using tv_detail::kPulseFields;
using tv_detail::pulseField;

// Grid geometry. Only the page build reads these.
namespace {

}

namespace {
    constexpr int n_template_cols = 3;
    constexpr int n_template_rows = 4;
    constexpr int max_templates_per_page = n_template_cols * n_template_rows;
}

std::vector<TemplateViewerWindow::Lead>
TemplateViewerWindow::leadsForBin(const TemplateBin& b) const {
    std::vector<Lead> out;
    // Shared with the markings CSV, so it cannot report a channel this
    // function says does not exist.
    static const char* kNames[3] = { "Ch1", "Ch2", "Ch3" };
    const ChannelTemplateData* ch[3] = { &b.ch1, &b.ch2, &b.ch3 };
    for (int c = 0; c < 3; ++c)
        if (ecgChannelPresent(b, c))
            out.push_back({ &ch[c]->ecgTemplate_raw,
                            &ch[c]->ecg_template_raw_iqr, c, kNames[c] });
    return out;
}

std::vector<TemplateViewerWindow::Lead>
TemplateViewerWindow::leadsForBinTemplate(const TemplateBin& b,
    int templateIdx) const {
    std::vector<Lead> out;
    static const char* kNames[3] = { "Ch1", "Ch2", "Ch3" };

    // WHICH ALIGNMENT THE GRID DRAWS. The panels used to be R-aligned always,
    // so pressing P/Q/J switched the focus panel's waveform while the grid
    // behind it kept showing the R-aligned average. Same selection now drives
    // both -- including in Automatic, where the grid follows whichever bar
    // was last clicked (currentGridAnchor), same as the focus panel, rather
    // than sitting on R regardless of what the operator clicked.
    const AnchorType gridAnchor = currentGridAnchor();

    for (int c = 0; c < 3; ++c) {
        const tbank::TemplateBank& bank = b.ecg_bank[c];

        // Slot 0 falls back to the chN_raw template when no bank reached this
        // bin, so a pre-bank file renders exactly as it always did.
        const std::vector<double>* trace = nullptr;
        // The spread for THAT trace. Paired with it from the same source, so
        // the band and the tail trim describe the waveform they are drawn
        // around.
        const std::vector<double>* traceIqr = nullptr;
        int nMembers = 0;
        uint8_t labelCode = tbank::kUnlabeled;
        uint8_t splitSource = tbank::kSplitUnknown;

        const bool pulseThin = templateIdx < b.ppg_bank.size() && b.ppg_bank.templates[templateIdx].tooFewBeats(/*is_ppg=*/true); //is there fewer ppgs than the given limit
        if (!pulseThin
            && templateIdx < bank.size()
            && !bank.templates[templateIdx].tmpl.empty()
            && !bank.templates[templateIdx].tooFewBeats(/*is_ppg=*/false)
            // wantsLandmarkMarking FOR EVERY SLOT. Slot 0 used to bypass it
            // (`templateIdx == 0 ||`), so the _A column appeared whatever the
            // template said about wanting landmarks while B and C obeyed the
            // rule. Either the predicate means something or it does not; it
            // cannot mean something for two columns out of three.
            && bank.templates[templateIdx].wantsLandmarkMarking()) {
            const tbank::BankTemplate& t = bank.templates[templateIdx];

            // ---- ONE SOURCE: THIS SLOT, THIS ALIGNMENT -------------------
            //
            // There used to be three waveforms reachable here, chosen by a
            // fallback chain, with nothing on screen saying which one you got:
            //
            //   asl->tmpl                        this slot, this alignment,
            //                                    averaged over members_clean
            //   t.tmpl                           this slot, R-ALIGNED only
            //   chFor(c, a).ecgTemplate_raw      the WHOLE BIN, no exclusions,
            //                                    built before the partition
            //                                    exists
            //
            // Three populations under one label. The third is the worst of
            // them: create_ecg_templates builds it from every baseline-valid
            // beat, so a dropped R detection -- whose "RR" is two cardiac
            // cycles and whose slice therefore spans four -- is in it, is the
            // only slice with samples past 1.8x the real RR, and the column
            // median out there IS its later complexes. That is the four-QRS
            // panel: bin 5 slot A, 986 members, a 2.84 s axis, beside slots B
            // and C at 1.40 and 1.30 s drawing t.tmpl correctly.
            //
            // NO FALLBACK NOW. If this slot has no average for this alignment
            // the panel is ABSENT, which is a visible gap rather than a silent
            // substitution of a different population. prepareViewerJob now
            // aligns all four anchors including R, so every (slot, anchor) has
            // one -- R was missing, which is why Automatic, the alignment the
            // operator starts on, was the one view with nothing to draw.
            //
            // EMPTY IS NOT THE ONLY WAY TO BE ABSENT. bankSlotFor only tests
            // tmpl.empty(), and alignTemplatesFromCache sizes each slot's
            // average to W NaNs and then fills the columns it has members for
            // -- so a slot with no members, or whose member indices fall
            // outside the aligned matrix, comes back correctly sized and
            // entirely NaN. Checked for a finite sample instead.
            // THROUGH slotView, the shared accessor -- so this trace and the
            // array ecgDetect measures on are the same object by construction.
            // This block used to call bankSlotFor itself and check for a finite
            // sample; bankSlotFor already rejects empty and all-NaN, and
            // slotView pairs the waveform with its r_col so the two cannot be
            // taken from different places.
            const SlotView svT = slotView(b, c, templateIdx, gridAnchor);
            if (!svT.valid || !svT.tmpl) {
                const char* why = "no-per-slot-average";
                fprintf(stderr, "[bank-trace] bin %llu lead %d slot %d anchor %s:"
                    " NO PER-SLOT AVERAGE (%s) -- no panel\n",
                    (unsigned long long)b.index, c, templateIdx,
                    anchor_view::label(gridAnchor), why);
                continue;
            }
            trace = svT.tmpl;
            traceIqr = svT.iqr;
            nMembers = t.memberCount();
            labelCode = t.label_code;
            splitSource = t.split_source;
            // subtype is no longer read here: tbank::letterRanks applies the
            // confirmed-subtype rule itself, from the same BankTemplate.
        }
        else {
            // No chN_raw fallback and no unconditional slot 0: both drew a
            // panel whose waveform was a different population from the
            // template it claimed to be. Covers every old exit -- ragged bank,
            // too few beats, thin pulse cohort, no average for this alignment.
            continue;
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
        QString lbl = QString("%1 %2_%3").arg(kNames[c]).arg(cls).arg(letter);

        // WHICH SIGNAL SPLIT THIS TEMPLATE OFF. A joint partition spawns a
        // group when one channel rejects the best existing one, and whether
        // that channel was an ECG lead or the pulse is the difference between
        // "this bin holds two morphologies" and "the pulse was noisy here".
        // The bin used to report that only as a per-bin tally in the log.
        //
        // ABSENT ON THE SEED, and on anything from an archive written before
        // the field existed: splitSourceLabel returns nullptr and nothing is
        // appended, rather than a word standing in for a value nobody
        // measured.
        if (const char* src = tbank::splitSourceLabel(splitSource))
            lbl += QString("  split %1").arg(src);

        out.push_back({ trace, traceIqr, c, lbl, nMembers });
    }
    return out;
}

// PACKED BY COLUMN, NOT BY BIN. A column is a (bin, template slot) pair, so
// a bin holding three morphologies contributes three of them and no bin's
// column count is known until its bank is. Taking whole bins therefore left
// SHORT PAGES: with a 3-column bin next in line a 12-cell grid stopped at 10
// and the operator saw two empty cells with no way to tell a paging rule from
// missing data. The grid is filled to the budget every page now and a bin's
// columns SPLIT across the boundary when they have to -- its _A template can
// end one page and its _B template open the next.
//
// THE COST IS REAL AND IS ACCEPTED: siblings from one bin are no longer
// guaranteed to be side by side, so comparing them can mean paging. Every
// panel still names its own bin (drawn in black, ahead of the gray template
// name), which is what makes a split column legible rather than confusing.
void TemplateViewerWindow::buildPages() {
    m_pages.clear();
    m_columnTable.clear();
    const int nBins = static_cast<int>(m_bins.size());

    // THE WHOLE RECORD'S COLUMNS, ONCE, IN DRAW ORDER. Pages are windows onto
    // this; markingSlotsForBin is asked exactly once per bin here instead of
    // once per packing trial, once per widest-page scan and again per page
    // build.
    for (int gi = 0; gi < nBins; ++gi)
        for (int t : markingSlotsForBin(m_bins[gi]))
            m_columnTable.push_back({ gi, t });

    // WHICH BINS HAVE NO COLUMNS AT ALL. With a minimum-beats threshold set,
    // every template in a bin can fall below it, and that bin then shows
    // nothing. Printed rather than left to be noticed, because a bin that
    // silently has no panel is indistinguishable from a paging bug. Reported
    // here, over the whole record: per page it printed the same bins again on
    // every visit and could not say how many there were in total.
    if (tbank::minBeatsEcg() > 0 || tbank::minBeatsPpg() > 0) {
        std::string gone;
        int nGone = 0;
        for (int gi = 0; gi < nBins; ++gi)
            if (markingSlotsForBin(m_bins[gi]).empty()) {
                ++nGone;
                if (nGone <= 40)
                    gone += (gone.empty() ? "" : ",") + std::to_string(gi);
            }
        if (nGone > 0)
            std::fprintf(stderr,
                "  [min-beats] %d bin(s) with NO displayable template "
                "(every slot below ECG %d or PPG %d clean beats): %s%s\n",
                nGone, tbank::minBeatsEcg(), tbank::minBeatsPpg(),
                gone.c_str(), (nGone > 40) ? ",..." : "");
    }

    // Cells a page holds: the whole grid when one lead is drawn per panel
    // (compact wraps), one column per lead-stack otherwise.
    const int budget = std::max(1, (max_leads <= 1)
        ? n_template_cols * n_template_rows : n_template_cols);
    const int nCols = static_cast<int>(m_columnTable.size());

    // No columns is a real state -- an empty record, or every template below
    // the minimum -- and it gets one empty page rather than none, so
    // showPage's clamp has something to land on.
    if (nCols == 0) { m_pages.push_back({ 0, 0 }); m_totalPages = 1; return; }

    for (int i = 0; i < nCols; i += budget)
        m_pages.push_back({ i, std::min(budget, nCols - i) });
    m_totalPages = static_cast<int>(m_pages.size());

    // How many pages hold a bin that continues onto the next one. Not a
    // warning: it is the cost of full pages, and it is worth knowing how
    // often it happens on a given record.
    int nStraddling = 0;
    for (size_t pg = 0; pg + 1 < m_pages.size(); ++pg) {
        const int lastOfPage = m_pages[pg].first + m_pages[pg].second - 1;
        if (m_columnTable[lastOfPage].first
            == m_columnTable[lastOfPage + 1].first) ++nStraddling;
    }
    std::fprintf(stderr,
        "  [pages] %d bins -> %d columns -> %d pages"
        " (%d cols/page, %d page(s) split a bin)\n",
        nBins, nCols, m_totalPages, budget, nStraddling);
    std::fflush(stderr);
}

bool TemplateViewerWindow::unionEcgFrameSeconds(const TemplateBin& b, int lead,
    int templateIdx, double& tMinSec, double& tMaxSec) const
{
    if (!(m_sampleRate > 0.0)) return false;
    if (lead < 0 || lead > 2) return false;

    // x = 0 is R for every alignment (the panel's origin is b.r_peak_ch, the
    // same column regardless of anchor), so an extent in columns maps to
    // seconds by (col - rPeak) / fs with one shared rPeak and rate.
    const double rPeak = static_cast<double>(b.r_peak_ch[lead]);

    // The same four alignments the ring walks. Kept local so this does not
    // depend on the anonymous-namespace kAlignRing defined further down.
    static const AnchorType kFour[4] = {
        AnchorType::P_ONSET, AnchorType::Q_ONSET,
        AnchorType::R_PEAK,  AnchorType::J_POINT,
    };

    int loCol = std::numeric_limits<int>::max();
    int hiCol = -1;
    for (AnchorType a : kFour) {
        // THROUGH slotView, like every other reader. The chN_raw fallback for
        // slot 0 is gone with it: sizing the frame from the bin's average while
        // the panel draws the per-slot one gave slot 0 an x-window belonging to
        // a different waveform.
        const SlotView svU = slotView(b, lead, templateIdx, a);
        if (!svU.valid || !svU.tmpl) continue;
        const std::vector<double>* trace = svU.tmpl;

        int f = -1, l = -1;
        for (int i = 0; i < (int)trace->size(); ++i)
            if (!std::isnan((*trace)[i])) { f = i; break; }
        for (int i = (int)trace->size() - 1; i >= 0; --i)
            if (!std::isnan((*trace)[i])) { l = i; break; }
        if (f < 0 || l <= f) continue;

        loCol = std::min(loCol, f);
        hiCol = std::max(hiCol, l);
    }

    if (hiCol < 0 || loCol == std::numeric_limits<int>::max() || hiCol <= loCol)
        return false;

    tMinSec = (static_cast<double>(loCol) - rPeak) / m_sampleRate;
    tMaxSec = (static_cast<double>(hiCol) - rPeak) / m_sampleRate;
    return true;
}

std::vector<int> TemplateViewerWindow::markingSlotsForBin(const TemplateBin& b) const {
    // ONE SOURCE, shared with the markings CSV writer, so the file's row set
    // IS this column set. Also per ALIGNMENT, which the old local lambda was
    // not: a slot can have an average for R and none for P.
    return visibleSlots(b, currentGridAnchor());
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
        pageIn();
    }
}

void TemplateViewerWindow::onPrevPage() {
    if (m_currentPage > 0) {
        captureCurrentPage();   // snapshot the page we are leaving
        --m_currentPage;
        pageIn();
    }
}

// Show the current page. No re-seed step: the fit modes are inputs to the
// detection barsForPanel reads, so a mode change is already reflected in what
// the next paint produces.
void TemplateViewerWindow::pageIn() {
    showPage();
}

// ========================================================================
// Build / clear plots
// ========================================================================

void TemplateViewerWindow::clearPlots() {
    // deleteLater, NOT delete. showPage() reaches this from a panel's own
    // signal -- mousePressEvent opens a context menu and emits
    // classConfirmRequested from inside its exec() loop -- so a plain delete
    // freed a panel whose handler was still on the stack, and the QMenu lost
    // its parent mid-loop. That presented as a permanent hang, not a crash.
    // deleteLater cannot do it: the event loop is only idle when no handler is
    // mid-call. disconnect is narrow (sender -> this) because the wildcard form
    // also severs destroyed(), which Qt wires internally.
    for (auto* pw : m_allPlots) {
        // NARROW, NOT pw->disconnect(). The no-argument form is a wildcard and
        // severs the panel's destroyed() connections too, which Qt itself wires
        // up for parent and layout bookkeeping -- hence
        // "wildcard call disconnects from destroyed signal" on every page
        // change. All that is wanted here is that a signal already emitted
        // cannot reach THIS window's slots while the panel is on its way out.
        disconnect(pw, nullptr, this, nullptr);
        ui->plotGrid->removeWidget(pw);
        pw->hide();
        pw->setParent(nullptr);
        pw->deleteLater();
    }
    m_allPlots.clear();
    m_binPlots.clear();
    m_pageTemplateIdx.clear();
    m_pageGlobalIdx.clear();
    m_pageColOf.clear();   // indexes panels that have just been deleted

    // Drop stretch factors left over from a previous (possibly larger) page
    // so unused rows/columns don't reserve empty space on the next page.
    for (int c = 0; c < ui->plotGrid->columnCount(); ++c)
        ui->plotGrid->setColumnStretch(c, 0);
    for (int r = 0; r < ui->plotGrid->rowCount(); ++r)
        ui->plotGrid->setRowStretch(r, 0);
}

// ---- expand this page's bins into (bin, template) COLUMNS ----------------
// A column is a bin plus a bank member, so a bin holding three morphologies
// occupies three adjacent columns and its seed stays leftmost. The page is a
// RANGE OF COLUMNS in m_columnTable, so this is a slice and nothing here
// recomputes eligibility -- buildPages settled which columns exist, and a
// second opinion formed here is how the page table and the drawn grid used to
// be able to disagree about how many panels a page holds.
std::vector<std::pair<int, int>>
TemplateViewerWindow::pageColumns(int start, int count) const
{
    const int n = static_cast<int>(m_columnTable.size());
    const int lo = std::clamp(start, 0, n);
    const int hi = std::clamp(start + std::max(0, count), lo, n);
    return std::vector<std::pair<int, int>>(m_columnTable.begin() + lo,
        m_columnTable.begin() + hi);
}

// ---- row count -----------------------------------------------------------
// Computed from the COLUMN count, not the bin count, because compact mode
// wraps PANELS and a bin can contribute several. Sizing the wrap from the bin
// count let a page with banks overflow past kMaxGridCols columns while still
// reporting a legal grid.
int TemplateViewerWindow::pageGridRows(bool compact,
    const std::vector<std::pair<int, int>>& cols) const
{
    const int nCols = static_cast<int>(cols.size());
    if (compact) {
        // Not `auto [rows, cols]`: the caller's `cols` is the (bin, template)
        // vector, and a structured binding there would shadow it.
        const std::pair<int, int> rc = compactGrid(nCols);
        return std::max(1, rc.first);   // .second is implied by the wrap
    }

    // One extra row for the VCG panel, but only if some bin on this page can
    // actually produce one (all three ECG channels present with an R column).
    // Without this the row count equals the lead count and there is no row
    // index left for the VCG to occupy.
    // THE BINS ON THIS PAGE COME FROM ITS COLUMNS. A page is a column range,
    // so there is no [start, end) over bins to walk any more -- and the first
    // and last bin of a page may each be showing only some of their columns.
    // Consecutive columns of one bin are adjacent, so tracking the previous
    // bin index is enough to probe each bin once.
    bool anyVcg = false;
    int prevBin = -1;
    for (const std::pair<int, int>& col : cols) {
        if (anyVcg) break;
        if (col.first == prevBin) continue;
        prevBin = col.first;
        anyVcg = !vcg_avg::derivedTraceOnChannelAxis(
            m_bins[col.first], 0, vcg::DerivedLead::VectorMagnitude).empty();
    }

    // Capped like the compact path. Three leads plus VCG already sits at the
    // cap; a fourth lead would drop the VCG row rather than grow the stack,
    // and vcgRowWanted in showPage (gridRows > max_leads) agrees with that on
    // its own, so the two cannot disagree about which row exists.
    //
    // STILL A THROWAWAY PROBE: the per-column loop recomputes this trace for
    // the same bins. Cheap (one interpolation pass per channel plus a 3x3
    // reconstruction, microseconds) and it needs the answer before the loop
    // that would produce it, so it stays -- but it is a real duplication and
    // this is where it lives if it ever matters.
    return std::min(n_template_rows, max_leads + (anyVcg ? 1 : 0));
}

// ---- VCG panel: the bottom row of one column, under lead 3 ---------------
// The trace is laid out on ch1's COLUMN axis by derivedTraceOnChannelAxis, so
// it shares the x axis of the lead panels above and lines up with them sample
// for sample -- while every channel is still sampled at its own r_col
// internally, which is what keeps the combination per-instant.
//
// Display only: no markers are set and no marker signals are connected, so
// nothing here is draggable. The VCG is derived from the three leads, so
// marking it would create a fourth set of fiducials with no channel of its own
// to store them in. That is also why it takes no template index.
void TemplateViewerWindow::addVcgPanel(int gi, int column, int gridRows,
    const TemplateBin& b,
    const std::vector<double>& vcgTrace, double vcgRCol,
    const global_intervals::GlobalIntervals& intervals,
    std::vector<BinPlotWidget*>& group, int& usedRows, int& usedCols)
{
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

    // Same global boundaries as the leads above, converted onto ch1's axis
    // since that is the axis this trace is drawn on.
    vp->setReferenceLines(
        global_interval_lines::forChannel(b, intervals, 0));

    ui->plotGrid->addWidget(vp, vcgRow, column, 1, 1);
    usedRows = std::max(usedRows, vcgRow + 1);
    usedCols = std::max(usedCols, column + 1);

    m_allPlots.push_back(vp);
    group.push_back(vp);
}

void TemplateViewerWindow::showPage() {
    clearPlots();

    // Page bounds are COLUMN indices from the packed table: every page holds
    // a full grid's worth except the last, and a bin's columns may continue
    // onto the next page. There is no bin range for the page.
    if (m_pages.empty()) buildPages();
    m_currentPage = std::clamp(m_currentPage, 0, (int)m_pages.size() - 1);
    const int start = m_pages[m_currentPage].first;
    const int count = m_pages[m_currentPage].second;

    bool compact = (max_leads <= 1);

    const std::vector<std::pair<int, int>> cols = pageColumns(start, count);
    const int nCols = static_cast<int>(cols.size());

    const int gridRows = pageGridRows(compact, cols);

    m_binPlots.resize(nCols);
    m_pageGlobalIdx.resize(nCols);
    m_pageTemplateIdx.resize(nCols);

    int usedRows = 0, usedCols = 0;

    for (int i = 0; i < nCols; ++i) {
        int gi = cols[i].first;
        const int template_index = cols[i].second;
        m_pageGlobalIdx[i] = gi;
        m_pageTemplateIdx[i] = template_index;

        const TemplateBin& b = m_bins[gi];
        auto leads = leadsForBinTemplate(b, template_index);

        // Only shout when no bank arrived at all: an empty bank and a genuinely
        // single-morphology bin render identically, so the difference has to be
        // written somewhere visible rather than inferred from the plots.
        if (template_index == 0)
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
            (template_index >= 0 && template_index < b.ppg_bank.size())
            ? &m_bins[gi].ppg_bank.templates[template_index] : nullptr;
        const bool hasPPG = ppgSlot && !ppgSlot->tmpl.empty()
            && ppgSlot->memberCount() > 0
            && !ppgSlot->tooFewBeats(/*is_ppg=*/true);
        const int nPpgForColumn = hasPPG ? ppgSlot->memberCount() : -1;

        if (leads.empty())
            leads.push_back({ nullptr, nullptr, 0, "No ECG" });

        // VCG needs all three channels; with fewer, the trace comes back empty
        // and no row is reserved, so the leads keep the full height.
        double vcgRCol = -1.0;
        const std::vector<double> vcgTrace = vcg_avg::derivedTraceOnChannelAxis(
            b, 0, vcg::DerivedLead::VectorMagnitude, vcg::kIdentity, &vcgRCol);
        const bool vcgRowWanted = !compact && !vcgTrace.empty()
            && gridRows > max_leads;

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

            int lead_index = leads[li].channelIndex;
            // The panel's x origin: the bin's R column, the SAME for every
            // alignment. Deriving it from the selected anchor instead moved the
            // axis rather than the waveform -- x=0 landed on that anchor's R
            // column, so the frame slid under a trace that had not changed
            // position, and two alignments could not be compared against a
            // common axis. The traces differ because the per-anchor averages
            // differ; the frame they are drawn in should not.
            const double rPeak = static_cast<double>(b.r_peak_ch[lead_index]);

            // FROM THE LEAD, not by channel. This picked
            // b.chN.ecg_template_raw_iqr -- the whole bin's spread -- while the
            // trace came from a bank slot, so the band around a slot's waveform
            // described a different population, and recomputeFrame's tail trim
            // read that band to decide where the trace stopped.
            static const std::vector<double> emptyIqr;
            const std::vector<double>& ecgIqrRaw =
                leads[li].ecgIqr ? *leads[li].ecgIqr : emptyIqr;
            const double ecgRef = (lead_index >= 0 && lead_index < 3) ? m_ecgGlobalRef[lead_index] : std::nan("");
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
            const uint64_t nEcgBinTotal = (lead_index == 0) ? b.ch1_n_beats_raw
                : (lead_index == 1) ? b.ch2_n_beats_raw
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

            const std::vector<double> ecgN = normalizeEcgTrace(ecgSrc, lead_index);
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

            // Pin the ECG x-window to the union of the four alignments' extents
            // so switching anchor (Automatic bar click, or a forced P/Q/R/J)
            // does not rescale the axis. Anchor-independent, so it survives the
            // in-place re-skin without being re-set. Only the ECG leads (0..2);
            // the VCG panel has no alignment dimension and keeps sizing itself.
            {
                double fLo = 0.0, fHi = 0.0;
                if (lead_index >= 0 && lead_index <= 2
                    && unionEcgFrameSeconds(b, lead_index, template_index, fLo, fHi))
                    pw->setEcgFrame(fLo, fHi);
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
            // ONE CALL FOR EVERY COLUMN. This was a templateIdx == 0 fork into
            // applyBinToWidget, which is what left slot 0 unseeded. m_bins[gi]
            // rather than `b`: the loop binds `b` as const, and the lazy seed
            // writes the marker set it just computed back into the template so
            // the next repaint and any drag see the same positions.
            applyTemplateToWidget(pw, m_bins[gi], lead_index, template_index);

            pw->setReferenceLines(global_interval_lines::forChannel(b, gi_intervals, lead_index));

            // ---- RESTORE THIS PANEL'S OWN MARK -----------------------
            // operator_state is per template. It used to read b.bad_ppg and
            // b.bad_r_ch[c], which are per BIN, so a rebuild painted every
            // panel of a marked bin -- the marks spread on a page turn.
            //
            // Slot 0 falls back to the bin flags when it has no state of its
            // own, so a record marked before operator_state existed, or one
            // whose flags feature_marks set automatically, still shows them.
            // The SAME function the right-click handlers repaint with, so a page
            // rebuild cannot disagree with a click. (This was inline, and it
            // tested marked_invalid_template == 2u for the pulse -- the field is
            // a bool, so that was never true and BadPPG never restored.)
            pw->setState(panelState(gi, lead_index, template_index));

            connect(pw, &BinPlotWidget::markerMovedOnTemplate, this, &TemplateViewerWindow::onMarkerMovedOnTemplate);
            connect(pw, &BinPlotWidget::markerDragStarted, this, &TemplateViewerWindow::onMarkerDragStarted);
            connect(pw, &BinPlotWidget::landmarkSelected, this, &TemplateViewerWindow::user_clicked_on_bar);
            // Glyph click: focus only -- refresh the panel, do NOT record a
            // touch or re-align (the glyph is the detector's answer, not a bar).
            connect(pw, &BinPlotWidget::landmarkFocusOnly, this,
                [this, pw](int binIdx, int leadIdx, int templateIdx, int marker, double col) {
                    refreshFocus(pw, binIdx, leadIdx, templateIdx, marker, col);
                });
            // (The two R-overlay connects lived here: rMarkerDragStarted,
            //  which forced the view to R so the drag happened on the R
            //  trace, and rMarkerMoved, which wrote r_bars_ch and propagated
            //  it across the page. Forced-R now marks the R cell through the
            //  ordinary bar path, so both are gone with the overlay -- and
            //  with them the last site that clamped an R-frame column against
            //  a drawn-frame wall.)
            connect(pw, &BinPlotWidget::badRToggled, this, &TemplateViewerWindow::onBadRToggled);
            pw->setTemplateIndex(template_index);

            if (lead_index >= 0 && lead_index <= 2 && template_index >= 0 && template_index < m_bins[gi].ecg_bank[lead_index].size())
            {
                m_bins[gi].ecg_bank[lead_index].templates[template_index].confirmed_by_operator = true;
            }
            if (template_index >= 0 && template_index < m_bins[gi].ppg_bank.size())
            {
                m_bins[gi].ppg_bank.templates[template_index].confirmed_by_operator = true;
            }
            // QUEUED, DELIBERATELY. This is emitted from inside
            // BinPlotWidget::mousePressEvent's QMenu::exec() loop, and
            // onClassConfirmRequested calls showPage() -> clearPlots(),
            // which destroys the emitting panel. Queued defers the slot
            // until mousePressEvent has returned and the menu has closed,
            // so the page is rebuilt with nothing of the old one on the
            // stack. Direct connection here is what hung the GUI.
            connect(pw, &BinPlotWidget::classConfirmRequested, this,
                &TemplateViewerWindow::onClassConfirmRequested,
                Qt::QueuedConnection);
            connect(pw, &BinPlotWidget::badPPGToggled, this, &TemplateViewerWindow::onBadPPGToggled);

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

        if (vcgRowWanted)
            addVcgPanel(gi, i, gridRows, b, vcgTrace, vcgRCol, gi_intervals,
                group, usedRows, usedCols);

        m_binPlots[i] = std::move(group);
    }

    // Equal stretch on every used row/column => equal-width, equal-height
    // cells that together fill the whole plot area. Combined with each
    // widget scaling its trace to its cell width, every bin window ends up
    // the same on-screen length.
    for (int c = 0; c < usedCols; ++c) ui->plotGrid->setColumnStretch(c, 1);
    for (int r = 0; r < usedRows; ++r) ui->plotGrid->setRowStretch(r, 1);

    // The column index, once, rather than a linear scan per lookup.
    m_pageColOf.clear();
    for (int li = 0; li < (int)m_pageGlobalIdx.size()
        && li < (int)m_pageTemplateIdx.size(); ++li)
        m_pageColOf[slotKey(m_pageGlobalIdx[li], m_pageTemplateIdx[li])] = li;

    applyMarkerVisibility();
    updatePageControls();
}

// ===========================================================================
// THE BARS THIS PANEL DRAWS -- AND THE ONLY PLACE THEY COME FROM.
//
// ONE DETECTION IN THE SYSTEM. pw->detectedLandmarks() is the panel's own
// m_det, which captureGlyphSnapshot also reads for every ECG glyph -- so a bar
// IS its glyph's column, not a copy that can drift from it.
//
// A CELL HOLDS OPERATOR EDITS AND NOTHING ELSE: -1 means untouched, and an
// untouched bar falls back to the detection on every paint. Clearing a cell is
// therefore the whole of "reset", and a fit-mode change needs no re-seed.
//
// showsBar decides which bars an alignment carries, so P shows one, Q three,
// R all four, J one -- unchanged.
tbank::BankMarkerSet TemplateViewerWindow::barsForPanel(const BinPlotWidget* pw,
    const TemplateBin& b, int lead, int slot) const
{
    tbank::BankMarkerSet out;   // all -1
    if (!pw) return out;

    const AnchorType a = currentGridAnchor();
    const FeatureMarks::TemplateLandmarks& lm = pw->detectedLandmarks();

    // showsBar GATES A FORCED VIEW ONLY. Automatic draws all four bars whatever
    // alignment it is anchored on -- and it re-anchors on a bar CLICK, so
    // gating it here made clicking the P bar hide the three Q bars.
    // EACH BAR'S EDIT COMES FROM THE CELL THE DRAG WROTE IT TO -- its OWNER
    // anchor (anchor_view::anchorFor), not the alignment on screen. Reading
    // them all from currentGridAnchor() worked only in Automatic, where a bar
    // click re-anchors the grid to that bar's own anchor so the two coincide;
    // force any alignment afterwards and every bar read an empty cell and fell
    // back to the detection, which looked like the edit being thrown away.
    //
    // The stored value is in the owner's columns, so it is converted into the
    // drawn frame here -- the inverse of the view -> owner conversion
    // moveEcgMarker applies on the way in.
    // WHICH CELL AN EDIT LIVES IN depends on the alignment, and that is the
    // point -- see anchor_view::ownsCanonicalBar.
    //
    //   forced R or J : this alignment's OWN cell, in its own columns, so no
    //                   conversion -- it was measured on the waveform drawn.
    //   otherwise     : the bar's canonical owner cell (P for p_begin, Q for
    //                   the rest), converted into the drawn frame. Automatic
    //                   and forced P/Q therefore show and edit one shared bar.
    const bool ownBars = m_forceAlign && anchor_view::hasOwnBars(a);

    const auto field = [](const tbank::BankMarkerSet& s, int marker) {
        if (marker == anchor_view::p_begin) return s.p_begin;
        if (marker == anchor_view::q_begin) return s.q_onset;
        if (marker == anchor_view::j_point)   return s.s_end;
        if (marker == anchor_view::t_end)   return s.t_end;
        return -1.0;
        };

    const auto pick = [&](int marker, double detected) -> double {
        if (m_forceAlign && !anchor_view::showsBar(a, marker)) return -1.0;
        if (ownBars) {
            const double edit = field(b.slotMarks(lead, slot, a), marker);
            if (edit >= 0.0) return edit;
        }
        else {
            const AnchorType owner = anchor_view::anchorFor(marker);
            const double edit = field(b.slotMarks(lead, slot, owner), marker);
            if (edit >= 0.0) return edit + b.frameShift(lead, owner, a);
        }
        return lm.valid ? detected : -1.0;     // otherwise the detection
        };
    out.p_begin = pick(anchor_view::p_begin, lm.p_begin);
    out.q_onset = pick(anchor_view::q_begin, lm.q_onset);
    out.s_end = pick(anchor_view::j_point, lm.s_end);
    out.t_end = pick(anchor_view::t_end, lm.t_end);
    return out;
}

// This column's panels, or nullptr when the (bin, slot) is not on this page.
const std::vector<BinPlotWidget*>* TemplateViewerWindow::panelsForColumn(
    int binIdx, int templateIdx) const
{
    const auto it = m_pageColOf.find(slotKey(binIdx, templateIdx));
    if (it == m_pageColOf.end()) return nullptr;
    if (it->second < 0 || it->second >= (int)m_binPlots.size()) return nullptr;
    return &m_binPlots[it->second];
}

// A bin's visible x-axis span in seconds. The sample rate cancels in the
// ratio, so this is the only term an equal-screen-distance shift needs.
double TemplateViewerWindow::binSpanSeconds(int binIdx) const {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return -1.0;
    double lo = 0.0, hi = 0.0;
    return unionEcgFrameSeconds(m_bins[binIdx], 0, 0, lo, hi) ? (hi - lo) : -1.0;
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
void TemplateViewerWindow::applyTemplateToWidget(BinPlotWidget* pw,
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
    // ---- THIS SLOT'S AVERAGE FOR ONE ALIGNMENT ---------------------------
    //
    // R is tp.tmpl: that is what an R-aligned slot average IS, and
    // alignTemplatesFromCache does not accumulate R into bank_anchors -- only
    // the re-aligned anchors go there. Every other alignment comes from
    // bankSlotFor, and a null is a WRITER GAP, not a state to render.
    //
    // NO FALLBACK TO tp.tmpl FOR A NON-R ANCHOR. There is one format version
    // and every section is written unconditionally (see template_io.cpp), so
    // "the file predates this section" is not a case. Substituting the
    // R-aligned average would put a bar in the R frame under a non-R tag, and
    // userMarks would then translate it OUT of a frame it was never in --
    // displacing it by r_col(R) - r_col(anchor), which is the defect this
    // whole function is being fixed for.
    //
    // r_col comes from the template when it has one, otherwise from the bin's
    // channel: every beat in the bank was aligned on the same R column by
    // construction, so the bin's value is correct rather than a guess when the
    // template's own field was never filled. detect_template_landmarks refines
    // it only locally -- symmetricExtremum is clamped to +/-7 samples and drops
    // to a five-point parabola when the residual guard trips on a sharp R -- so
    // this is not a rough hint, it is very nearly the answer, and it is the
    // search origin every other finder brackets on.
    const int rColR = (tp.r_col >= 0)
        ? tp.r_col
        : static_cast<int>(std::lround(b.r_peak_ch[channel]));

    // slotView is shared, so the bar-seeding sites elsewhere use the same
    // answer -- they used to seed from the R-aligned average under another
    // alignment's tag because this logic was a lambda they could not call.
    auto slotWaveform = [&](AnchorType a, const std::vector<double>*& w,
        int& rc) -> bool {
            const SlotView sv = slotView(b, channel, templateIdx, a);
            if (!sv.valid) return false;
            w = sv.tmpl; rc = sv.r_col;
            return true;
        };

    for (AnchorType a4 : anchor_view::anchor_array) {
        const int tag4 = static_cast<int>(a4);
        if (tp.hasDetectedMarks(tag4)) continue;

        // THIS ANCHOR'S WAVEFORM AND THIS ANCHOR'S R COLUMN. Both used to be
        // the R-frame ones on all four passes -- every bar measured on tp.tmpl
        // and stored under its owner's tag -- so a sub-template's Q-onset bar
        // came back from userMarks displaced by r_col(R) - r_col(Q).
        //
        // The per-slot anchored averages are row subsets of the bin's aligned
        // matrix, so they share the bin's frame and the bin's per-anchor r_col
        // is the right seed for them.
        const std::vector<double>* w4 = nullptr;
        int rc4 = -1;
        if (!slotWaveform(a4, w4, rc4)) {
            // Reported, not papered over: this is build_templates failing to
            // write a per-slot average for an anchor it aligned. The bar stays
            // absent, which is the same answer refreshFocus gives for the same
            // gap.
            fprintf(stderr, "[bank-seed] bin=%llu lead=%d slot=%d anchor=%s"
                " NO PER-SLOT ALIGNED AVERAGE -- bar not seeded\n",
                (unsigned long long)b.index, channel, templateIdx,
                anchor_view::label(a4));
            continue;
        }
        FeatureMarks::seed_bank_template(*w4, rc4, m_sampleRate, 1.0, a4,
            tp.marks(tag4));
    }
    // (no glyph sync: p_peak is derived at every read now -- see the
    //  reactive_ecg call in applyBankTemplateToWidget below.)
    // ASSEMBLED, NOT FETCHED. Each bar lives in its owning alignment's set;
    // userMarks pulls all four and translates them into the R frame the panel
    // draws in. Fetching one alignment's set here is what limited a session to
    // editing one alignment's bars.
    // setAuto FIRST, BEFORE THE BARS ARE READ. It is what sets m_bin and
    // m_frame on the panel, and barsForPanel asks the panel for its detection
    // -- so with setAuto later, the bars came from a detection made against
    // the panel's PREVIOUS alignment. On a reskin (forced Q, then back)
    // setEcgData has already dropped the cache, so that stale-frame detection
    // was recomputed on the wrong average and P came back out of range, which
    // keep() folds to -1: the P bar simply vanished.
    pw->setAuto(b, currentGridAnchor());

    const tbank::BankMarkerSet mk =
        barsForPanel(pw, b, channel, templateIdx);
    // Bin first, for the arterial markers (a bank slot has no ABP/ART waveform
    // of its own). The seven ECG bars are overridden below; the PULSE bars are
    // overridden here, from this slot's own waveform.
    //
    // BARS ONLY. A glyph push used to follow them (overridePulseGlyphs), and
    // the panel paints the pulse glyphs from its own detectedPulse() now -- on
    // this same ps.tmpl, so the bars and the X marks are one measurement.
    applyBinCommonToWidget(pw, b);

    if (templateIdx >= 0 && templateIdx < b.ppg_bank.size()
        && !b.ppg_bank.templates[templateIdx].tmpl.empty()) {
        tbank::BankTemplate& ps = b.ppg_bank.templates[templateIdx];
        {
            if (!ps.hasDetectedPulseMarks())
                FeatureMarks::seed_pulse_bank_template(ps.tmpl, m_ppgRateHz,
                    ps.pulse_marks);
            const tbank::BankPulseMarkerSet& pm = ps.pulse_marks;
            // THREE BARS, THE REST DERIVED. peak / peak2 / t50 / t80 left
            // BankPulseMarkerSet -- markerAtX hands none of them out, and all
            // of them come back from reactive_ppg bracketed by the bars plus
            // the detector's own systolic peak.
            const FeatureMarks::ReactivePpg rp = FeatureMarks::reactive_ppg(
                ps.tmpl, pm.onset, pm.peak_auto, pm.dicrotic, pm.end);
            pw->setMarker(BinPlotWidget::PpgOnset, pm.onset);
            pw->setMarker(BinPlotWidget::PpgPeak, pm.peak_auto);
            pw->setMarker(BinPlotWidget::PpgDicrotic, pm.dicrotic);
            pw->setMarker(BinPlotWidget::PpgPeak2, rp.peak2);
            pw->setMarker(BinPlotWidget::PpgEnd, pm.end);
            pw->setMarker(BinPlotWidget::PpgT50, rp.t50);
            pw->setMarker(BinPlotWidget::PpgT80, rp.t80);
        }
    }
    else {
        // NO PULSE FOR THIS SLOT, so no pulse bars. Left unset they would keep
        // whatever the previous occupant of this widget had -- the bin's values
        // used to be pushed here unconditionally, which hid it.
        pw->setMarker(BinPlotWidget::PpgOnset, -1.0);
        pw->setMarker(BinPlotWidget::PpgPeak, -1.0);
        pw->setMarker(BinPlotWidget::PpgDicrotic, -1.0);
        pw->setMarker(BinPlotWidget::PpgPeak2, -1.0);
        pw->setMarker(BinPlotWidget::PpgEnd, -1.0);
        pw->setMarker(BinPlotWidget::PpgT50, -1.0);
        pw->setMarker(BinPlotWidget::PpgT80, -1.0);
    }

    // P PEAK IS DERIVED, not stored: reactive_ecg on this slot's own bars,
    // against this slot's own waveform. Same function BinPlotWidget::
    // reactiveGlyphs calls, so the bar set, the X on screen and the CSV column
    // cannot disagree.
    // ON THE WAVEFORM THIS PANEL DRAWS. tp.tmpl is the slot's R-aligned
    // average, so bracketing it with bars in the P/Q/J frame measured the P
    // peak on one waveform with another's columns. slotView pairs the trace
    // with its own R column for the alignment on screen.
    const SlotView svDraw = slotView(b, channel, templateIdx, currentGridAnchor());
    const std::vector<double>& drawnEcg = svDraw.valid ? *svDraw.tmpl : tp.tmpl;
    const FeatureMarks::ReactiveEcg reBank = FeatureMarks::reactive_ecg(
        drawnEcg, mk.p_begin, mk.q_onset, mk.s_end, mk.t_end, m_sampleRate);
    pw->setMarker(BinPlotWidget::EcgPBegin, mk.p_begin);
    pw->setMarker(BinPlotWidget::EcgPPeak, reBank.p_peak);
    pw->setMarker(BinPlotWidget::EcgQBegin, mk.q_onset);
    // R FROM THE DETECTION, NOT THE BOOKKEPT COLUMN.
    //
    // svDraw.r_col is alignment.hpp's R_anchor + median(applied shifts) with a
    // +-5 ms snap: a bookkeeping reference column, never refined and never put
    // through a peak fit. The GLYPH draws the refined position
    // (m_glyphs.ecgRPeak, from detect_template_landmarks), so setting the bar
    // from r_col left the two a few samples apart -- and on a P/Q/J average,
    // where the snap window cannot reach the smeared apex, several samples
    // apart. That is the bar sitting off the peak while the focus panel, which
    // re-fits, sits on it.
    //
    // bin_plot_widget.cpp's glyph hit-test already documented the mismatch:
    // "NOT m_markers[EcgRPeak], which holds r_col_raw and sits a few samples
    // off the drawn cross -- that mismatch is why clicking the R glyph used to
    // do nothing." Taking both from one detection closes that too: the bar and
    // its glyph are now the same column, so the X is clickable.
    //
    // setAuto above has already run, so detectedLandmarks() is populated for
    // this bin and this alignment. r_col remains the fallback for a template
    // whose detection did not produce an R at all.
    {
        const FeatureMarks::TemplateLandmarks& lmR = pw->detectedLandmarks();
        const double rFallback =
            svDraw.valid ? static_cast<double>(svDraw.r_col) : rColR;
        pw->setMarker(BinPlotWidget::EcgRPeak,
            (lmR.valid && lmR.r_peak >= 0.0) ? lmR.r_peak : rFallback);
    }
    pw->setMarker(BinPlotWidget::EcgSEnd, mk.s_end);
    pw->setMarker(BinPlotWidget::EcgTEnd, mk.t_end);

    // ---- THIS SLOT'S OWN GLYPHS, ON THE TRACE THIS PANEL DRAWS -----------
    //
    // (NOTHING PUSHES GLYPHS ANY MORE, on either channel. setAuto names the
    //  (bin, slot, alignment) and the panel detects that waveform itself, once,
    //  when something asks -- the paint that draws the X marks, the hit test
    //  that makes them clickable, the R bar above, and the focus panel, all off
    //  one answer. Two override passes lived here, each re-detecting on the raw
    //  slot array in Auto fit modes and overwriting the detection setAuto had
    //  just made in the operator's modes: a full detector run per panel per
    //  apply, discarded, and the surviving answer was the one that ignored the
    //  radios.)
}

// ---------------------------------------------------------------------------
// THE BIN-WIDE HALF, SHARED BY EVERY COLUMN OF THE BIN.
//
// This was applyBinToWidget, "the slot 0 path", and it set the ECG bars, the
// pulse bars and the arterial bars. That made slot 0 the one column nobody
// seeded, the one column whose pulse bars came from the bin instead of from the
// pulse it draws, and the one column whose P peak was measured on
// chFor(c, R_PEAK) -- hardcoded R -- with bars in the current frame. Those were
// three separate bugs with one cause: slot 0 was special.
//
// It now does only what is genuinely PER BIN: the arterial channels (one ABP /
// ART / ART_PULM trace per bin, no bank to hang them off), the alignment badge
// and the glyph snapshot. Everything per-template is in
// applyTemplateToWidget, which every column goes through, slot 0 included.
// ---------------------------------------------------------------------------
void TemplateViewerWindow::applyBinCommonToWidget(BinPlotWidget* pw,
    const TemplateBin& b) {
    // ARTERIAL ONLY out of the field table. The pulse entries are per slot now
    // and applyTemplateToWidget sets them from this column's own pulse_marks;
    // pushing the bin's values here first would put a bar from another
    // waveform on screen for however long it took to be overwritten.
    for (const PulseField& f : kPulseFields) {
        if (BinPlotWidget::markerIsPpg(f.marker)) continue;
        pw->setMarker(static_cast<BinPlotWidget::Marker>(f.marker), b.*f.field);
    }
    // BADGE FOR R AND J ONLY: dotted, one colour, labelled (the `overlay`
    // branch in paintEvent). Those two alignments carry their own bars, so the
    // different style says something true -- the operator is editing this
    // alignment's measurement, not the shared one. P and Q show the canonical
    // bars and keep the per-landmark palette.
    const bool ownBadge = m_forceAlign && anchor_view::hasOwnBars(m_forcedAlign);
    pw->setAlignmentBadge(ownBadge ? anchor_view::label(m_forcedAlign)
        : nullptr);

    // Glyphs in the frame of the alignment the grid is drawing.

}

// In-place ECG re-skin of every panel on the current page to
// currentGridAnchor(). Mirrors the ECG half of showPage()'s per-column body
// -- per-anchor average from leadsForBinTemplate, the same optional notch,
// the same /ref normalization -- but pushes it through setEcgData (ECG trace,
// band, R column and count only) and re-applies the glyphs, leaving the pulse
// channels, the layout and the widgets themselves alone. That is what lets it
// run mid-click without breaking a drag.
void TemplateViewerWindow::reskinGridForAnchor() {
    const bool notchActive = m_notchFilterOn && m_notchFilterHz > 0;

    for (int i = 0; i < (int)m_binPlots.size()
        && i < (int)m_pageGlobalIdx.size()
        && i < (int)m_pageTemplateIdx.size(); ++i)
    {
        const int gi = m_pageGlobalIdx[i];
        const int templateIndex = m_pageTemplateIdx[i];
        if (gi < 0 || gi >= (int)m_bins.size()) continue;
        const TemplateBin& b = m_bins[gi];

        const auto leads = leadsForBinTemplate(b, templateIndex);

        for (auto* pw : m_binPlots[i]) {
            if (!pw) continue;
            const int lead_index = pw->leadIndex();
            // 0..2 only: the VCG panel (a fourth "lead" row) is derived from
            // all three R-aligned leads and has no alignment of its own, so it
            // neither needs nor wants a re-anchor.
            if (lead_index < 0 || lead_index > 2) continue;

            // The per-anchor average for THIS lead, picked exactly as showPage
            // does. If this lead has no entry (e.g. absent channel), skip it.
            const Lead* L = nullptr;
            for (const auto& cand : leads)
                if (cand.channelIndex == lead_index) { L = &cand; break; }
            if (!L || !L->ecg) continue;

            static const std::vector<double> emptyIqr;
            const std::vector<double>& ecgRaw = *L->ecg;
            const std::vector<double>& ecgIqrRaw = L->ecgIqr ? *L->ecgIqr : emptyIqr;

            // ECG notch, footIdx = -1 (no rebase; ECG normalizes by /ref, not a
            // foot) -- the same call showPage makes.
            std::vector<double> ecgSrc = ecgRaw;
            if (notchActive && !ecgSrc.empty())
                ecgSrc = notch_filter(ecgSrc,
                    static_cast<double>(m_notchFilterHz), m_sampleRate);

            const std::vector<double> ecgN = normalizeEcgTrace(ecgSrc, lead_index);
            const double ecgRef = m_ecgGlobalRef[lead_index];
            const std::vector<double> ecgIqr =
                normalize_features::scale_array_by_ref(ecgIqrRaw, ecgRef);

            const double rPeak = static_cast<double>(b.r_peak_ch[lead_index]);
            const uint64_t nEcgBinTotal = (lead_index == 0) ? b.ch1_n_beats_raw
                : (lead_index == 1) ? b.ch2_n_beats_raw : b.ch3_n_beats_raw;
            const uint64_t nEcgBeats = (L->nMembers > 0)
                ? static_cast<uint64_t>(L->nMembers) : nEcgBinTotal;

            pw->setEcgData(ecgN, ecgIqr, rPeak, static_cast<int>(nEcgBeats));

            // Glyphs + the R glyph column follow the anchor too; these helpers
            // already read currentGridAnchor(). The draggable bars they also
            // re-apply are read from their stored R-framed positions, which at
            // click time (before any move) equal what is on screen, so nothing
            // the operator is holding jumps.
            applyTemplateToWidget(pw, m_bins[gi], lead_index, templateIndex);
        }
    }
}

// Repaints exactly the panel that was clicked.
// A PANEL'S STATE IS A FUNCTION OF BOTH VERDICTS, not of whichever one just
// changed. Each handler used to repaint with its own -- BadR or BadPPG -- so the
// third right-click set BadBoth in the widget and then had it overwritten by
// whichever signal was handled last, landing back on BadR.
//
// marked_invalid_template IS A BOOL (template_bank.hpp), so the ECG/PPG
// distinction is NOT in the value -- it is in WHICH SLOT carries the flag: the
// ECG lead's slot for ECG, the pulse slot for pulse. Testing the value against
// 1 or 2 cannot work; a bool never equals 2.
BinPlotWidget::State TemplateViewerWindow::panelState(int binIdx, int leadIdx,
    int templateIdx) const
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size())
        return BinPlotWidget::State::Good;
    const TemplateBin& b = m_bins[binIdx];

    bool ecgBad = false, ppgBad = false;
    if (leadIdx >= 0 && leadIdx <= 2 && templateIdx >= 0 && templateIdx < b.ecg_bank[leadIdx].size())
        ecgBad = b.ecg_bank[leadIdx].templates[templateIdx].marked_invalid_template;
    if (templateIdx >= 0 && templateIdx < b.ppg_bank.size()) {
        ppgBad = b.ppg_bank.templates[templateIdx].marked_invalid_template;
    }
    if (templateIdx == 0 && leadIdx >= 0 && leadIdx <= 2 && b.bad_r_ch[leadIdx]) {
        ecgBad = true;
    }
    return (ecgBad && ppgBad) ? BinPlotWidget::State::BadBoth
        : ppgBad ? BinPlotWidget::State::BadPPG
        : ecgBad ? BinPlotWidget::State::BadR
        : BinPlotWidget::State::Good;
}

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