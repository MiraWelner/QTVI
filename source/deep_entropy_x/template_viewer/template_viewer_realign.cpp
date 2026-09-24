// ========================================================================
// template_viewer_realign.cpp
// Operator-guided pulse re-stacking, plus the two display helpers the page
// and the re-stack share.
//
// THE GESTURE. Drag the pulse foot bar to where the foot really is, release.
// On release this re-finds each member beat's own trough near that column,
// re-medians the stack about it, and pushes the new trace into the panels of
// that column in place. Nothing else on the page moves.
//
// AND THE SELECTION, which is the same operation with a different anchor. The
// "Align PPG Horizontal" group picks WHAT the stack lines up on: the build's
// own choice (Auto), the foot, or a percentage up the upstroke -- in
// amplitude, so 0 is the foot and 10 is a tenth of the way up the systolic
// rise. Choosing one re-stacks every pulse column on the page; choosing Auto
// puts them back. A foot drag moves the group to Foot, because after that drag
// the foot IS the alignment and a control saying otherwise is a false claim
// about the trace beside it.
//
// WHAT IT DOES NOT DO, AND WHY THAT IS WRITTEN HERE AS WELL AS IN
// ppg_realign.hpp: the cohort is not re-selected. The beats averaged are
// exactly the ones the build chose, under the build's own alignment and its
// correlation filter. So the waveform that results is anchored one way and was
// populated under another. That is a deliberate, bounded compromise -- re-
// filtering would change group membership, which changes the morphology split,
// which is the pipeline's decision -- and it is why the status line says the
// bin was re-stacked and m_ppgRealigned remembers which ones.
//
// One of six units; TemplateViewerWindow is declared in template_viewer.hpp.
// ========================================================================

#include "template_viewer.hpp"
#include "template_generation/ecg_realign.hpp"

// ========================================================================
// Shared display preparation
// ========================================================================

// The display-time notch, with the rebase. Lifted verbatim out of showPage's
// lambda so the re-stack path cannot notch a trace differently from the page.
//
// WHY THE REBASE. notch_filter is x - narrow_bandpass(x): the shape survives
// but the DC level moves by whatever residual mean the bandpass emits. That
// shift lands on foot_y, which is what normalize_pulse_trace DIVIDES BY, and
// for a pulse stored near baseline foot_y sits near zero -- so a tiny DC shift
// is a large relative change and the normalized trace collapses to the bottom
// of the axis, looking like the channel vanished. Rebasing to the pre-notch
// foot value keeps the divisor stable.
//
// footIdx < 0 means "no foot" (the ECG's case, and the detector's honest answer
// on a pulse it could not resolve): filter, do not rebase.
std::vector<double> TemplateViewerWindow::maybeNotchTrace(
    const std::vector<double>& sig, double fs, double footIdx) const
{
    if (!(m_notchFilterOn && m_notchFilterHz > 0)) return sig;
    if (sig.empty() || fs <= 0.0) return sig;

    std::vector<double> out =
        notch_filter(sig, static_cast<double>(m_notchFilterHz), fs);

    const int fi = static_cast<int>(std::lround(footIdx));
    if (footIdx >= 0.0 && fi >= 0
        && fi < (int)sig.size() && fi < (int)out.size()) {
        const double shift = sig[fi] - out[fi];
        if (std::isfinite(shift) && shift != 0.0)
            for (auto& v : out) v += shift;
    }
    return out;
}

// Trace + band + the foot both are measured against, for one pulse-bank slot.
//
// THE FOOT IS THE SLOT'S OWN pulse_marks.onset, seeded on first use. Not the
// bin's ppg_onset: the waveform here is this slot's median, and the bin field
// was measured on b.ppgTemplate -- a different pulse. The band is scaled about
// the same foot, so trace and band cannot disagree about where the baseline is.
bool TemplateViewerWindow::pulseTraceForSlot(tbank::BankTemplate& slot,
    std::vector<double>& outTrace,
    std::vector<double>& outIqr,
    double& outFootIdx)
{
    outTrace.clear();
    outIqr.clear();
    outFootIdx = -1.0;
    if (slot.tmpl.empty()) return false;

    if (!slot.hasDetectedPulseMarks())
        FeatureMarks::seed_pulse_bank_template(slot.tmpl, m_ppgRateHz,
            slot.pulse_marks);
    outFootIdx = slot.pulse_marks.onset;

    const std::vector<double> src =
        maybeNotchTrace(slot.tmpl, m_ppgRateHz, outFootIdx);

    outTrace = normalize_ppg_or_similar(src, outFootIdx, 0);
    outIqr = normalize_features::scale_pulse_spread_by_ref(
        slot.tmpl_iqr,
        normalize_features::sample_y(slot.tmpl, outFootIdx),
        m_pulseGlobalRef[0]);
    return true;
}

// ========================================================================
// The re-stack
// ========================================================================

// ---- ADOPT A RE-AVERAGED PAIR WITHOUT CHANGING THE PANEL'S EXTENT --------
//
// The beat matrix is padded to the bin's LONGEST slice; the bank slot the panel
// draws is trimmed to roughly where its own membership ends. Those are
// different widths -- bin 1 drew an 880-column slot off a 1483-column matrix,
// and two slots of bin 5 drew 747 and 819 off 1466 -- so a re-average that
// emits the full matrix width hands the panel a longer waveform than the one it
// replaced. The x-frame stretches, and the columns it gained are the far tail
// where only the few longest beats contribute.
//
// ppg_realign's per-column occupancy gate already NaNs most of that out. This
// clips what is left to the incoming length, so a re-stack or a re-level can
// only ever change a template's VALUES, never how much axis it occupies. A
// result SHORTER than the old template is left alone: that is the gate saying
// this stacking supports fewer columns, which is true and worth seeing.
//
// ---- THE SPREAD IS RECOMPUTED, AND THAT IS THE POINT ----------------------
//
// The corridor is the EVIDENCE the alignment took. If every beat really does
// share a landmark at the operator's column then the cross-beat spread must be
// small there and grow away from it -- so a band that tightens at the bar is
// how the operator sees the correction worked, and one that does not is how
// they see it did not. Keeping the build's spread would pair a measurement
// with a waveform it no longer describes and remove the only readout this
// gesture has.
//
// SO THE PINCH AT THE BAR IS THE SIGNAL, NOT A DEFECT. local_ratio_iqr maps
// each beat to (sample - footY_i)/footY_i, which is identically 0 at the column
// footY_i was read from: the spread there is zero by construction and SHOULD
// be. The build has the same property at each beat's own foot; a shared column
// simply puts it in one place, where it can be read.
static void adoptPulsePair(tbank::BankTemplate& slot,
    std::vector<double> tmpl, std::vector<double> iqr)
{
    const std::size_t keep = slot.tmpl.size();
    if (keep > 0 && tmpl.size() > keep) {
        tmpl.resize(keep);
        if (iqr.size() > keep) iqr.resize(keep);
    }
    slot.tmpl = std::move(tmpl);
    slot.tmpl_iqr = std::move(iqr);
    slot.band_lo.clear();
    slot.band_hi.clear();
}


QString TemplateViewerWindow::beatsBinPath() const {
    // SAME DIRECTORY AND STEM morphology_csv::set was handed in
    // analysis_job::prepare, which is what wrote the file. Derived on each call
    // rather than cached at load, so it cannot survive a subject change and
    // point at the previous record's beats.
    if (m_templateDir.isEmpty() || m_subjectId.isEmpty()) return {};
    return m_templateDir + "/" + m_subjectId + "_beats.bin";
}

// Release of a dragged bar. TWO KINDS OF BAR MATTER HERE: the four ECG landmark
// bars, whose anchor's per-slot average is re-stacked on the operator's column,
// and the pulse foot, whose cohort is re-levelled on it. Everything else is
// reactive and was already applied during the drag, so this returns immediately
// for those rather than doing work per gesture.
//
// ON RELEASE, NOT ON MOVE, AND THAT IS THE WHOLE REASON THIS SLOT EXISTS. A
// re-stack reads the bin's beat matrix off disk and re-medians a few hundred
// rows; hanging that off markerMoved would run it once per mouse-move event
// for the entire drag, so the bar would lurch and every intermediate column
// would produce a template nobody asked for. The gesture is "the landmark is
// HERE", and a gesture is only complete on mouse-up.
//
// AND IT IS A GESTURE, NOT A PANEL. In Move-Subsequent the operator's claim
// covers every column right of the dragged one -- m_dragPropCols is the list
// the propagation wrote -- so each of those is re-stacked too, on its own
// propagated column. Doing only the dragged panel left the other eleven with a
// bar in one place and a stack built around another, which is the exact defect
// this slot exists to remove, hidden on the panels nobody was watching.
void TemplateViewerWindow::onMarkerReleasedOnTemplate(int binIdx, int leadIdx,
    int templateIdx, int marker, int newIdx)
{
    // FIRST, BEFORE ANY BRANCH RETURNS. The propagated columns were painted on
    // a clock during the drag (see flushDragRepaints) and the last move is
    // usually inside the interval, so the gesture is not finished on screen
    // until this runs -- and it has to run for every marker, including the
    // ones this function then ignores.
    flushDragRepaints(/*force=*/true);

    // ---- THE ECG BARS: RE-STACK THAT ANCHOR'S AVERAGE ------------------
    //
    // A landmark bar's own alignment is the one built by aligning every beat on
    // that landmark, and until now the operator could move the BAR but not the
    // STACK -- so the two disagreed about where the landmark was, on a waveform
    // whose target column came from a detector that (per the [anchors] summary)
    // frequently could not find it at all. R is excluded because it is the base
    // frame every other anchor aligns FROM, not an alignment of its own.
    if (BinPlotWidget::markerIsEcg(marker) && anchor_view::isBar(marker)
        && marker != BinPlotWidget::EcgRPeak) {
        const AnchorType a = anchor_view::anchorFor(marker);
        if (a == AnchorType::R_PEAK) return;

        // ---- EVERY COLUMN THE GESTURE MOVED, NOT JUST THE ONE UNDER THE
        //      CURSOR ---------------------------------------------------
        //
        // Move-Subsequent means "this landmark is HERE on this column, and by
        // the same screen distance on every column right of it". That claim is
        // about the waveforms, not about the bars: a propagated bar that lands
        // on a new column while its average stays stacked on the old one
        // re-creates, on eleven columns, precisely the disagreement this
        // gesture exists to remove -- and it is worse there than on the dragged
        // column, because the operator is not looking at those panels and has
        // no cue that the trace under the bar did not follow.
        //
        // READ EVERY BAR FIRST, THEN RE-STACK. realignEcgFromBar reskins the
        // column it writes, which re-pushes bars and re-runs a detection, so
        // collecting the targets up front keeps the values being acted on the
        // ones the drag actually produced.
        //
        // ASCENDING COLUMN ORDER IS ALSO CACHE ORDER: a bin owns one column per
        // bank slot and those columns are adjacent on the page, so sibling
        // slots of one bin are re-stacked back to back and beatsForBin's
        // one-deep cache serves all of them from a single read.
        struct Pending { int bin; int slot; double col; };
        std::vector<Pending> also;
        also.reserve(m_dragPropCols.size());
        for (const int li : m_dragPropCols) {
            if (li < 0 || li >= (int)m_pageGlobalIdx.size()
                || li >= (int)m_pageTemplateIdx.size()) continue;
            const int gi = m_pageGlobalIdx[li];
            const int sl = m_pageTemplateIdx[li];
            if (gi < 0 || gi >= (int)m_bins.size() || sl < 0) continue;
            if (gi == binIdx && sl == templateIdx) continue;   // the dragged one
            // THE BAR AS DRAWN, on the panel for the DRAGGED LEAD -- an ECG bar
            // is per lead and the propagation wrote only that lead's cells, so
            // a column whose panel for this lead is absent has nothing to
            // re-stack.
            double col = -1.0;
            if (li < (int)m_binPlots.size())
                for (auto* pw : m_binPlots[li])
                    if (pw && pw->leadIndex() == leadIdx) {
                        col = pw->marker(
                            static_cast<BinPlotWidget::Marker>(marker));
                        break;
                    }
            if (col >= 0.0) also.push_back(Pending{ gi, sl, col });
        }
        m_dragPropCols.clear();

        // ONE MESSAGE FOR THE WHOLE GESTURE when there is more than one column
        // in it: the per-column message would be overwritten eleven times and
        // the operator would read whichever column happened to be last.
        const bool many = !also.empty();
        const bool draggedOk = realignEcgFromBar(binIdx, leadIdx, templateIdx,
            a, static_cast<double>(newIdx), /*announce=*/!many);
        if (!many) return;

        // A CURSOR, BECAUSE THIS IS NOW O(COLUMNS) BEAT MATRICES. One column's
        // re-stack is a locator run per member beat (parallelised inside
        // realignAt); a page of them is that many times over, plus a file read
        // per distinct bin that is not already in memory.
        QGuiApplication::setOverrideCursor(Qt::WaitCursor);
        int ok = 0;
        for (const Pending& pcol : also)
            if (realignEcgFromBar(pcol.bin, leadIdx, pcol.slot, a, pcol.col,
                /*announce=*/false)) ++ok;
        QGuiApplication::restoreOverrideCursor();

        statusBar()->showMessage(tr("%1 re-stacked on %2 of %3 columns "
            "(the dragged column %4). A refused column keeps the average it "
            "had and its bar still moved. Cohorts unchanged - not re-filtered.")
            .arg(anchor_view::label(a))
            .arg(ok + (draggedOk ? 1 : 0))
            .arg(static_cast<int>(also.size()) + 1)
            .arg(draggedOk ? tr("included") : tr("refused")), 8000);
        return;
    }

    (void)leadIdx;   // the pulse is per (bin, slot); leads share it
    // [ppg-dbg]
    fprintf(stderr, "[ppg-dbg] released: bin=%d slot=%d marker=%d newIdx=%d "
        "(want PpgOnset=%d)\n", binIdx, templateIdx, marker, newIdx,
        (int)BinPlotWidget::PpgOnset);
    if (marker != BinPlotWidget::PpgOnset) return;

    // THE FOOT BAR IS THE VERTICAL REFERENCE, so this re-LEVELS the cohort at
    // the operator's column and re-medians it. It does NOT re-time anything.
    //
    // It used to call realignPulseFromFoot, which shifts beats sideways -- the
    // wrong axis for this bar, and on rows that arrive already up50
    // time-aligned that computed a shift near zero for every one of them, so
    // the gesture appeared to do nothing at all. The radio group owns the
    // horizontal axis; the foot bar owns the vertical one. Two gestures, two
    // axes, and this one no longer touches the radio group.
    // ---- AND EVERY PROPAGATED FOOT, FOR THE SAME REASON ----------------
    //
    // movePpgMarker re-normalizes the DISPLAY of each propagated column during
    // the drag (pushPulseToPanels), which is the cheap half; the cohort
    // re-level off _beats.bin is the half that only a mouse-up can afford, and
    // it was being done for the dragged column alone. So every other column's
    // trace was drawn against a foot its member beats had never been levelled
    // on. The foot is read from each column's own pulse_marks -- pulse marks
    // are per (bin, slot) and every lead of a column shares them.
    std::vector<std::pair<int, int> > alsoPulse;   // (bin, slot)
    alsoPulse.reserve(m_dragPropCols.size());
    for (const int li : m_dragPropCols) {
        if (li < 0 || li >= (int)m_pageGlobalIdx.size()
            || li >= (int)m_pageTemplateIdx.size()) continue;
        const int gi = m_pageGlobalIdx[li];
        const int sl = m_pageTemplateIdx[li];
        if (gi < 0 || gi >= (int)m_bins.size() || sl < 0) continue;
        if (gi == binIdx && sl == templateIdx) continue;   // the dragged one
        alsoPulse.push_back(std::make_pair(gi, sl));
    }
    m_dragPropCols.clear();

    const bool manyPulse = !alsoPulse.empty();
    const bool footOk = relevelPulseAtFoot(binIdx, templateIdx,
        static_cast<double>(newIdx), /*announce=*/!manyPulse);
    if (!manyPulse) return;

    QGuiApplication::setOverrideCursor(Qt::WaitCursor);
    int okPulse = 0;
    for (const std::pair<int, int>& pc : alsoPulse) {
        const TemplateBin& tb = m_bins[pc.first];
        if (pc.second >= (int)tb.ppg_bank.size()) continue;
        const double foot =
            tb.ppg_bank.templates[pc.second].pulse_marks.onset;
        if (!(foot >= 0.0)) continue;
        if (relevelPulseAtFoot(pc.first, pc.second, foot, /*announce=*/false))
            ++okPulse;
    }
    QGuiApplication::restoreOverrideCursor();

    statusBar()->showMessage(tr("Pulse re-levelled on %1 of %2 columns "
        "(the dragged column %3). A refused column keeps the average it had "
        "and its bar still moved. Cohorts unchanged - not re-filtered.")
        .arg(okPulse + (footOk ? 1 : 0))
        .arg(static_cast<int>(alsoPulse.size()) + 1)
        .arg(footOk ? tr("included") : tr("refused")), 8000);
}

// ---- ONE BIN'S BEAT MATRIX, CACHED ONE DEEP --------------------------------
//
// The original note here said: read per gesture, and if it ever shows up as a
// delay, cache the LAST bin read -- not all of them. It has: a percent change
// re-stacks every column on the page, and a page is up to eight columns of
// which several are usually sibling slots of the SAME bin, so an uncached read
// re-walks the file once per column to return the same matrix.
//
// ONE ENTRY, and it is a cache of the file rather than of the alignment -- the
// rows are what the build sliced and no re-stack ever writes to them, so it
// cannot go stale within a subject. It is keyed only by bin, so it MUST be
// dropped when the subject changes (initAfterBinsLoaded does that); bin 3 of
// the next record would otherwise be served bin 3 of this one.
const ppg_realign::BinBeats& TemplateViewerWindow::beatsForBin(int binIdx,
    const char* channel)
{
    static const ppg_realign::BinBeats kNone;
    if (binIdx < 0 || !channel) return kNone;
    if (m_beatsCacheBin == binIdx && m_beatsCacheChan == channel)
        return m_beatsCache;

    // ---- MEMORY FIRST, AND NORMALLY THERE IS NO SECOND ------------------
    //
    // per_channel_beats["PPG"][bin] IS this bin's kept pulses, [beat][sample],
    // in slice order -- which is the LOCAL ROW SPACE members_clean indexes, so
    // there is no join to do. loadBin has to reconstruct that space by
    // counting became_beat columns in file order; here it is simply the vector.
    //
    // This is also why the disk path was the wrong design and not just badly
    // ordered: the file is written by a deferred task on the thread running
    // finalize(), started at the same instant this window opens, so it does
    // not exist when the operator first drags a bar. See setBeats.
    if (m_beatsInMemory) {
        const auto it = m_beatsInMemory->per_channel_beats.find(channel);
        if (it != m_beatsInMemory->per_channel_beats.end()
            && binIdx < (int)it->second.size()) {
            const std::vector<std::vector<double>>& rows = it->second[binIdx];
            ppg_realign::BinBeats out;
            out.rows = rows;   // one bin, cached one deep below
            for (const auto& row : out.rows)
                out.width = std::max(out.width, (int)row.size());
            m_beatsCache = std::move(out);
            m_beatsCacheBin = binIdx;
            m_beatsCacheChan = channel;
            return m_beatsCache;
        }
    }

    // FALLBACK: the file. Reachable from the path overload of loadSubject,
    // which has no BeatsFile to be handed, and from a re-run where the archive
    // from a previous pass is on disk.
    const QString path = beatsBinPath();
    if (path.isEmpty()) return kNone;

    m_beatsCache = ppg_realign::loadBin(path.toStdString(),
        static_cast<uint32_t>(binIdx), channel);
    m_beatsCacheBin = binIdx;
    m_beatsCacheChan = channel;
    return m_beatsCache;
}

void TemplateViewerWindow::clearBeatsCache()
{
    m_beatsCacheBin = -1;
    m_beatsCacheChan.clear();
    m_beatsCache = ppg_realign::BinBeats{};
}

// ---- THE WAVEFORM THE BUILD PRODUCED, KEPT ---------------------------------
//
// A re-stack overwrites slot.tmpl in place, which is what makes the panel
// update; it also means the build's own alignment is gone the moment the first
// gesture lands. That was acceptable while the only control was a drag -- there
// is no "un-drag" -- but "Auto" is a POSITION OF A RADIO GROUP, and a control
// the operator can return to has to return to something. So the built pair is
// stashed on first overwrite and restored by Auto.
//
// VIEWER-ONLY, like m_ppgRealigned, and for the same reason: it is a copy of
// what the pipeline already wrote, held so this window can put it back.
void TemplateViewerWindow::stashBuiltPulse(int binIdx, int templateIdx,
    const tbank::BankTemplate& slot)
{
    const int key = slotKey(binIdx, templateIdx);
    if (m_ppgBuilt.count(key)) return;   // FIRST overwrite only
    m_ppgBuilt[key] = { slot.tmpl, slot.tmpl_iqr };
}

// NO CALLER AS OF THE Auto REWRITE, AND KEPT ON PURPOSE. "Auto" used to mean
// "the stacking the build produced" and this was how it got back there; Auto
// now decides an alignment per column like every other position of the group,
// so nothing asks for the built pair any more. The pair is still STASHED on
// every first overwrite (stashBuiltPulse) -- autoPctForSlot depends on it to
// keep its verdict stable -- so the state this restores is still there, and a
// fourth "As built" radio position would need only to call this.
bool TemplateViewerWindow::restorePulseAsBuilt(int binIdx, int templateIdx)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return false;
    TemplateBin& b = m_bins[binIdx];
    if (templateIdx < 0 || templateIdx >= (int)b.ppg_bank.size()) return false;

    const auto it = m_ppgBuilt.find(slotKey(binIdx, templateIdx));
    if (it == m_ppgBuilt.end()) return false;   // never re-stacked: nothing to undo

    tbank::BankTemplate& slot = b.ppg_bank.templates[templateIdx];
    slot.tmpl = it->second.first;
    slot.tmpl_iqr = it->second.second;
    // The corridor is dropped on the way out of a re-stack and is not part of
    // the stash, so it stays dropped: it described neither stacking and is
    // rebuilt by the next recomputeTemplate.
    slot.band_lo.clear();
    slot.band_hi.clear();

    m_ppgRealigned.erase(slotKey(binIdx, templateIdx));
    m_ppgBuilt.erase(it);
    pushPulseToPanels(binIdx, templateIdx);
    return true;
}

// The display half of a re-stack: normalize this slot through the shared
// pulse-display path and push it into every panel of its column, in place.
// Factored out because the restore path needs exactly the same push and
// copying it is how the two would come to normalize against different feet.
void TemplateViewerWindow::pushPulseToPanels(int binIdx, int templateIdx,
    bool alsoFocus)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    TemplateBin& b = m_bins[binIdx];
    if (templateIdx < 0 || templateIdx >= (int)b.ppg_bank.size()) return;
    tbank::BankTemplate& slot = b.ppg_bank.templates[templateIdx];

    std::vector<double> trace, iqr;
    double footIdx = -1.0;
    if (!pulseTraceForSlot(slot, trace, iqr, footIdx)) return;

    const std::vector<BinPlotWidget*>* col =
        panelsForColumn(binIdx, templateIdx);
    // [ppg-dbg] A null column is the silent failure that looks exactly like
    // "nothing happened": the re-level can succeed, write slot.tmpl, and have
    // nowhere on screen to push it. footIdx < 0 is the other one -- no foot
    // means normalize_ppg_or_similar has no divisor.
    fprintf(stderr, "[ppg-dbg] push: bin=%d slot=%d footIdx=%.1f trace=%zu "
        "iqr=%zu panels=%d\n", binIdx, templateIdx, footIdx, trace.size(),
        iqr.size(), col ? (int)col->size() : -1);
    if (col) {
        for (auto* pw : *col)
            if (pw) pw->setPpgData(trace, iqr, slot.memberCount());
    }

    // The focus panel, if it is showing this pulse, is re-run from the new
    // trace rather than left displaying the old stacking magnified.
    //
    // alsoFocus false for the LIVE path: movePpgMarker runs its own
    // refreshFocus at the end of the same event, and two per mouse-move is one
    // wasted magnified re-detect per pixel of drag.
    if (alsoFocus && m_focusMarker >= 0
        && m_focusBin == binIdx && m_focusSlot == templateIdx)
        refreshFocus(m_focusWidget, m_focusBin, m_focusLead, m_focusSlot,
            m_focusMarker, m_focusCol);
}

bool TemplateViewerWindow::realignPulseFromFoot(int binIdx, int templateIdx,
    double footCol, bool announce, double pctOverride)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return false;
    TemplateBin& b = m_bins[binIdx];
    if (templateIdx < 0 || templateIdx >= (int)b.ppg_bank.size()) return false;
    tbank::BankTemplate& slot = b.ppg_bank.templates[templateIdx];
    if (slot.tmpl.empty()) return false;

    // THE COHORT: members_clean when it exists, because that is the set the
    // waveform on screen was actually averaged over -- `members` includes the
    // premature and Tukey-excluded beats that were deliberately kept in the
    // record but out of the average. Re-averaging `members` would quietly
    // re-admit them, which would look like the re-alignment had changed the
    // morphology.
    const std::vector<uint32_t>& cohort = !slot.members_clean.empty()
        ? slot.members_clean : slot.members;
    if (cohort.empty()) {
        if (auto* sb = announce ? statusBar() : nullptr)
            sb->showMessage(tr("No pulse members recorded for this template; "
                "nothing to re-stack."), 5000);
        return false;
    }

    const QString path = beatsBinPath();
    if (path.isEmpty()) return false;

    // One seeking pass per bin, held one deep. See beatsForBin.
    const ppg_realign::BinBeats& beats = beatsForBin(binIdx);
    if (beats.empty()) {
        // The usual cause is a run with morphology_csv::setWriteBeatsBin(false),
        // which truncates the file on purpose. Say which file, because the
        // alternative is a control that silently does nothing.
        if (auto* sb = announce ? statusBar() : nullptr)
            sb->showMessage(tr("No per-beat pulse data in %1 - cannot re-stack "
                "(the bar was still moved).").arg(path), 6000);
        return false;
    }

    // ---- WHAT THE RADIO GROUP SELECTED --------------------------------
    //
    // Foot is Percent(0) and is not a second code path: upstrokePctCol
    // returns the foot itself for pct <= 0, so the only thing the Foot
    // position does is pin the percentage to zero.
    //
    // THE PERCENT POINT IS MEASURED ON THE AVERAGE, not on a beat, and that is
    // the whole reason it is measured here and passed in rather than derived
    // inside realignAt. slot.tmpl is this cohort's median: if the upstroke is
    // resolvable anywhere it is resolvable there, whereas a single beat is the
    // noisy measurement the percent alignment exists to get away from. The
    // beats are then each located individually and shifted ONTO that column.
    //
    // AN OVERRIDE WINS OVER THE CONTROLS, and only the Auto path passes one:
    // Auto's percentage is a per-column verdict (autoPctForSlot) and there is
    // no control holding it, so it arrives as an argument rather than through
    // m_ppgAlignPercent. Negative means "ask the radio group", which is what
    // the foot drag and both explicit positions do.
    const double pct = (pctOverride >= 0.0)
        ? std::clamp(pctOverride, 0.0, 100.0)
        : ((m_ppgAlignMode == PpgAlign::Percent)
            ? std::clamp(static_cast<double>(m_ppgAlignPercent), 0.0, 100.0)
            : 0.0);

    double anchorCol = footCol;
    if (pct > 0.0) {
        anchorCol = ppg_realign::upstrokePctCol(slot.tmpl,
            static_cast<int>(std::lround(footCol)), pct);
        if (!(anchorCol >= 0.0)) {
            // REFUSED, AND THE OLD TEMPLATE STANDS -- the same rule as a
            // refused re-stack below. No upstroke means the percentage has no
            // referent on this pulse, and stacking on the foot instead would
            // silently give this column a different alignment from every
            // other one on the page while the radio said otherwise.
            if (auto* sb = announce ? statusBar() : nullptr)
                sb->showMessage(tr("No resolvable upstroke above that foot - "
                    "cannot align this pulse at %1%% (the bar was still moved).")
                    .arg(pct), 6000);
            return false;
        }
    }

    const ppg_realign::Result res = ppg_realign::realignAt(
        beats, cohort, footCol, ppg_realign::searchHalfWinFor(m_ppgRateHz),
        pct, anchorCol);
    if (!res.ok) {
        // REFUSED, AND THE OLD TEMPLATE STANDS. A re-median over one or two
        // beats is noise, and replacing a real template with it while reporting
        // success is the worst available outcome.
        if (auto* sb = announce ? statusBar() : nullptr)
            sb->showMessage(tr("Pulse re-stack refused: %1.")
                .arg(QString::fromStdString(res.why)), 6000);
        return false;
    }

    // ---- ADOPT THE NEW WAVEFORM ---------------------------------------
    //
    // The bar is NOT moved: the operator put it where the foot is and the
    // re-stack anchored every beat to that column, so the trace has come to
    // the bar rather than the reverse. pulse_marks.onset was already written by
    // the drag handler (movePpgMarker), which is why nothing here touches it.
    //
    // THE BUILD'S OWN PAIR IS KEPT FIRST, before it is overwritten -- this is
    // the only moment it still exists. See stashBuiltPulse: Auto is a position
    // the operator can return to.
    stashBuiltPulse(binIdx, templateIdx, slot);
    adoptPulsePair(slot, res.tmpl, res.iqr);

    m_ppgRealigned.insert(slotKey(binIdx, templateIdx));

    // ---- PUSH IT INTO THE LIVE PANELS, IN PLACE -----------------------
    //
    // setPpgData, not showPage(): a page rebuild would destroy the panel this
    // gesture came from, drop the focus view, and need a captureCurrentPage
    // first to avoid losing marker edits. The pulse average is the only thing
    // that changed, and every panel in this column draws the same one.
    pushPulseToPanels(binIdx, templateIdx);

    // SAID OUT LOUD, including the part that is a compromise. n_used vs
    // n_members is how much of the cohort had an anchor near the new column: a
    // large gap means the correction disagrees with most of the beats, which is
    // worth seeing before trusting the result.
    //
    // WHICH POINT IT STACKED ON IS PART OF THE SENTENCE. "re-stacked on col
    // 214" is ambiguous the moment the alignment can be something other than
    // the foot, and the difference between a foot alignment and a 10% one is
    // invisible on a thumbnail-sized panel.
    if (auto* sb = announce ? statusBar() : nullptr) {
        const QString where = (res.pct > 0.0)
            ? tr("%1% up the upstroke (col %2, foot col %3)")
            .arg(res.pct).arg(res.anchor_col).arg(res.foot_col)
            : tr("the foot (col %1)").arg(res.foot_col);
        sb->showMessage(tr("Pulse re-stacked on %1: %2 of %3 beats "
            "(median shift %4, max %5 samples). Cohort unchanged - not re-filtered.")
            .arg(where)
            .arg(res.n_used).arg(res.n_members)
            .arg(res.median_shift).arg(res.max_shift), 8000);
    }
    return true;
}

// ========================================================================
// The "Align PPG Horizontal" selection
// ========================================================================

// ---- WHAT "Auto" DECIDES, PER COLUMN --------------------------------------
//
// The foot is the worst-conditioned landmark on a pulse -- zero slope through
// it by definition -- so the percent alignment exists for the columns where
// stacking on it cannot be trusted. Auto asks one question of each column: how
// tight is the stack at the foot? Under kAutoFootIqrMax it takes
// kAutoFallbackPct up the upstroke instead; above it, the foot.
//
// THE BUILD'S SPREAD, NOT THE SLOT'S CURRENT ONE. A re-stack overwrites
// slot.tmpl_iqr with the spread of its own output, so reading that field would
// have the second Auto pass over a column judge the result of the first -- and
// a column could flip between the foot and 10% on nothing but a page turn.
// m_ppgBuilt holds the build's pair from the first overwrite onward; before
// that there has been no overwrite and the slot's own field IS the build's.
// Either way the verdict for a given column is the same every time it is
// asked, which is what makes "Auto" a position rather than a history.
int TemplateViewerWindow::autoPctForSlot(int binIdx, int templateIdx,
    const tbank::BankTemplate& slot, double footCol) const
{
    const auto it = m_ppgBuilt.find(slotKey(binIdx, templateIdx));
    const std::vector<double>& iqr = (it != m_ppgBuilt.end())
        ? it->second.second
        : slot.tmpl_iqr;

    const int halfwin = std::max(1, static_cast<int>(
        std::lround(kAutoIqrWindowSec * m_ppgRateHz)));
    const double spread = ppg_realign::iqrAbout(iqr, footCol, halfwin);

    // NO SPREAD MEASURABLE -> THE FOOT. An absent or all-zero band is not
    // evidence for the fallback: it means this slot's spread was never
    // written, or fewer than two beats reached the columns around the foot.
    // The foot is where the operator's bar is and what they would expect.
    if (!(spread >= 0.0)) return 0;

    return (spread < kAutoFootIqrMax) ? kAutoFallbackPct : 0;
}

// ONE PLACE THAT CHANGES THE PULSE ALIGNMENT, for the same reason
// applyAlignmentSelection is the one place that changes the ECG one: the radio
// buttons, the spin box and the foot-drag override all land here, and each of
// them carrying its own copy of "set the members, re-stack the page, sync the
// controls" is three chances for one of them to leave the group checked on an
// alignment that is not on screen.
void TemplateViewerWindow::applyPpgAlignSelection(PpgAlign mode, int pct)
{
    // [ppg-dbg]
    fprintf(stderr, "[ppg-dbg] align selection: mode=%d pct=%d\n",
        (int)mode, pct);
    m_ppgAlignMode = mode;
    // CLAMPED ON THE WAY IN, not trusted from the widget. The spin box's range
    // is set in wirePpgAlignButtons, but m_ppgAlignPercent is also read by the
    // re-stack path, and a percentage outside [0, 100] there means an anchor
    // outside the foot-to-peak span -- which upstrokePctCol would answer with
    // the peak column, silently.
    m_ppgAlignPercent = std::clamp(pct, 0, 100);
    syncPpgAlignControls();
    realignAllVisiblePulses();
}

void TemplateViewerWindow::setPpgAlignMode(PpgAlign mode)
{
    // NO RE-STACK, and that is the difference from applyPpgAlignSelection.
    // This is the path a foot DRAG uses to make the group say "Foot": the
    // re-stack for that gesture is the drag's own, on the dragged column
    // alone, and running a whole-page one from here would re-stack every other
    // column on the page as a side effect of touching one bar.
    m_ppgAlignMode = mode;
    syncPpgAlignControls();
}

// ---- RE-STACK WHAT IS ON SCREEN -------------------------------------------
//
// THE PAGE, NOT THE RECORD. A record is thousands of columns; re-stacking all
// of them on a radio click would read every bin's beat matrix and freeze the
// window for minutes, to produce waveforms for columns the operator may never
// look at. The page is what a selection has to be true of to be believed, and
// showPage calls this so paging forward brings the new columns into line
// (see the tail of showPage).
//
// EACH COLUMN FROM ITS OWN FOOT. The percent point is measured up from the
// foot, so the hint is this slot's own pulse_marks.onset -- the operator's bar
// where they have moved one, the detector's seed where they have not. There is
// deliberately no single column shared across the page: the bins have
// different rates and different RR, so one column is a different instant in
// each of them.
void TemplateViewerWindow::realignAllVisiblePulses()
{
    int nDone = 0, nSkipped = 0;
    // Auto's split, which is the only thing worth reporting about it: how many
    // columns it judged stackable on their own foot and how many it moved up
    // the upstroke. (nRestored went with restorePulseAsBuilt -- Auto no longer
    // puts the build's stacking back, it decides an alignment like every other
    // position of the group.)
    int nAutoFoot = 0, nAutoPct = 0;
    // [ppg-dbg] An empty page table means no column was considered at all,
    // which is a different problem from every column refusing.
    fprintf(stderr, "[ppg-dbg] realign page: mode=%d pct=%d cols=%zu/%zu\n",
        (int)m_ppgAlignMode, m_ppgAlignPercent,
        m_pageGlobalIdx.size(), m_pageTemplateIdx.size());

    for (int li = 0; li < (int)m_pageGlobalIdx.size()
        && li < (int)m_pageTemplateIdx.size(); ++li) {
        const int gi = m_pageGlobalIdx[li];
        const int slotIdx = m_pageTemplateIdx[li];
        if (gi < 0 || gi >= (int)m_bins.size() || slotIdx < 0) continue;

        TemplateBin& b = m_bins[gi];
        if (slotIdx >= (int)b.ppg_bank.size()) continue;
        tbank::BankTemplate& slot = b.ppg_bank.templates[slotIdx];
        if (slot.tmpl.empty()) continue;
        // A column the operator has already called bad is left exactly as it
        // is: re-stacking it would spend a file read to improve a waveform
        // that has been excluded.
        if (b.bad_ppg == 1 || slot.marked_invalid_template) continue;

        // Seed before reading, as every other pulse-mark reader does: a column
        // that has never been displayed has no foot yet, and -1 is not a hint.
        //
        // FOR EVERY MODE INCLUDING AUTO. Auto used to return above this, put
        // the build's waveform back and never look at the pulse marks at all;
        // it now stacks on the foot like the other positions, so it needs the
        // same hint they need.
        if (!slot.hasDetectedPulseMarks())
            FeatureMarks::seed_pulse_bank_template(slot.tmpl, m_ppgRateHz,
                slot.pulse_marks);
        const double foot = slot.pulse_marks.onset;
        if (!(foot >= 0.0)) { ++nSkipped; continue; }

        // AUTO DECIDES HERE, PER COLUMN, and passes its answer down as an
        // override -- there is no control holding it. Every other mode leaves
        // this negative and realignPulseFromFoot reads the radio group.
        double pctOverride = -1.0;
        if (m_ppgAlignMode == PpgAlign::Auto) {
            pctOverride = static_cast<double>(
                autoPctForSlot(gi, slotIdx, slot, foot));
            if (pctOverride > 0.0) ++nAutoPct; else ++nAutoFoot;
        }

        // THE RETURN VALUE, not a membership test on m_ppgRealigned. A slot
        // already in that set stays in it when a later re-stack is REFUSED, so
        // asking the set would report a refusal as a success.
        if (realignPulseFromFoot(gi, slotIdx, foot, /*announce=*/false,
            pctOverride)) ++nDone;
        else ++nSkipped;
    }

    // ONE LINE FOR THE WHOLE PAGE. Per-column messages would each overwrite
    // the last and the operator would see only whichever column happened to be
    // drawn last -- and the number that matters here is how many columns the
    // selection actually reached, not what happened to any one of them.
    if (auto* sb = statusBar()) {
        if (m_ppgAlignMode == PpgAlign::Auto) {
            // THE SPLIT, NOT A TOTAL. "Auto - 6 columns re-stacked" would hide
            // the only thing the operator cannot see on a thumbnail: which
            // columns were judged too tight at the foot to stack on it.
            sb->showMessage(tr("Pulse alignment: auto - %1 on the foot, "
                "%2 at %3% (foot IQR below %4), %5 unchanged. "
                "Cohorts unchanged - not re-filtered.")
                .arg(nAutoFoot).arg(nAutoPct).arg(kAutoFallbackPct)
                .arg(kAutoFootIqrMax).arg(nSkipped), 8000);
        }
        else {
            const QString where = (m_ppgAlignMode == PpgAlign::Percent
                && m_ppgAlignPercent > 0)
                ? tr("%1% up the upstroke").arg(m_ppgAlignPercent)
                : tr("the foot");
            sb->showMessage(tr("Pulse alignment: %1 - %2 column(s) re-stacked, "
                "%3 unchanged. Cohorts unchanged - not re-filtered.")
                .arg(where).arg(nDone).arg(nSkipped), 8000);
        }
    }
}

// ========================================================================
// The foot bar: re-level the cohort where the operator put it
// ========================================================================
//
// THE VERTICAL COUNTERPART OF realignPulseFromFoot, and deliberately its own
// function rather than a mode of it: the two share the cohort rule and the
// status line and nothing else. One shifts rows sideways onto an instant; this
// one shifts them up and down onto a baseline. A single function with an axis
// flag would be two disjoint bodies under one name.
//
// WHAT IT DOES NOT DO, as ever: the cohort is not re-selected and the pulse QC
// filter is not re-run. See ppg_realign.hpp.
bool TemplateViewerWindow::relevelPulseAtFoot(int binIdx, int templateIdx,
    double footCol, bool announce)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return false;
    TemplateBin& b = m_bins[binIdx];
    if (templateIdx < 0 || templateIdx >= (int)b.ppg_bank.size()) return false;
    tbank::BankTemplate& slot = b.ppg_bank.templates[templateIdx];
    if (slot.tmpl.empty()) return false;

    // members_clean when it exists -- the set the waveform on screen was
    // averaged over. Re-averaging `members` would re-admit the premature and
    // Tukey-excluded beats and look like the correction had changed the
    // morphology.
    const std::vector<uint32_t>& cohort = !slot.members_clean.empty()
        ? slot.members_clean : slot.members;
    if (cohort.empty()) {
        if (auto* sb = announce ? statusBar() : nullptr)
            sb->showMessage(tr("No pulse members recorded for this template; "
                "nothing to re-level."), 5000);
        return false;
    }

    const QString path = beatsBinPath();
    if (path.isEmpty()) return false;
    const ppg_realign::BinBeats& beats = beatsForBin(binIdx);
    // [ppg-dbg] width/rows say whether the file was found, whether the PPG
    // block was found inside it, and whether this bin has any rows in it --
    // three different failures that all end up as beats.empty().
    fprintf(stderr, "[ppg-dbg] relevel: bin=%d slot=%d cohort=%zu foot=%.1f "
        "beats.width=%d beats.rows=%zu path=%s\n",
        binIdx, templateIdx, cohort.size(), footCol, beats.width,
        beats.rows.size(), path.toStdString().c_str());
    if (beats.empty()) {
        if (auto* sb = announce ? statusBar() : nullptr)
            sb->showMessage(tr("No per-beat pulse data in %1 - cannot re-level "
                "(the bar was still moved).").arg(path), 6000);
        return false;
    }

    const ppg_realign::RelevelResult res = ppg_realign::relevelAt(
        beats, cohort, footCol, ppg_realign::levelHalfWinFor(m_ppgRateHz));
    if (!res.ok) {
        // REFUSED, AND THE OLD TEMPLATE STANDS.
        if (auto* sb = announce ? statusBar() : nullptr)
            sb->showMessage(tr("Pulse re-level refused: %1.")
                .arg(QString::fromStdString(res.why)), 6000);
        return false;
    }

    // [ppg-dbg] THE NUMBERS THAT SEPARATE "BAND TOO BIG" FROM "TRACE TOO
    // FLAT". tmpl[foot] is what normalize_ppg_or_similar DIVIDES BY, so a large
    // one collapses the trace toward zero and makes any band look enormous
    // beside it. old vs new iqr at the bar and 40 samples on says whether the
    // corridor actually tightened where it should.
    {
        const int f = std::clamp(res.foot_col, 0,
            (int)std::max<std::size_t>(1, res.tmpl.size()) - 1);
        const int m = std::min<int>((int)res.tmpl.size() - 1, f + 40);
        auto at = [](const std::vector<double>& v, int i) {
            return (i >= 0 && i < (int)v.size()) ? v[i]
                : std::numeric_limits<double>::quiet_NaN();
            };
        fprintf(stderr, "[ppg-dbg] relevel adopt: baseline=%.4f "
            "tmplOld[foot]=%.4f tmplNew[foot]=%.4f | iqrOld foot=%.4f +40=%.4f"
            " | iqrNew foot=%.4f +40=%.4f\n",
            res.baseline, at(slot.tmpl, f), at(res.tmpl, f),
            at(slot.tmpl_iqr, f), at(slot.tmpl_iqr, m),
            at(res.iqr, f), at(res.iqr, m));
    }

    // The bar is NOT moved: pulse_marks.onset was written by movePpgMarker
    // during the drag and the rows have come to it. Stash the build's pair
    // first -- this is the last moment it exists. See stashBuiltPulse.
    stashBuiltPulse(binIdx, templateIdx, slot);
    adoptPulsePair(slot, res.tmpl, res.iqr);

    m_ppgRealigned.insert(slotKey(binIdx, templateIdx));
    pushPulseToPanels(binIdx, templateIdx);

    // n_used vs n_members is how much of the cohort had samples at the
    // operator's column at all; median and max offset say how far the rows had
    // to move to agree there, which is the honest measure of how much the
    // correction disagreed with the build's own levelling.
    if (auto* sb = announce ? statusBar() : nullptr)
        sb->showMessage(tr("Pulse re-levelled on foot col %1 (baseline %2): "
            "%3 of %4 beats, median offset %5, max %6. "
            "Cohort unchanged - not re-filtered.")
            .arg(res.foot_col).arg(res.baseline)
            .arg(res.n_used).arg(res.n_members)
            .arg(res.median_level).arg(res.max_level), 8000);
    return true;
}

// ========================================================================
// The ECG bars: re-stack one anchor's per-slot average
// ========================================================================

// THE WRITE TARGET, MUTABLY. TemplateBin::bankSlotFor is const and also
// REJECTS an all-NaN slot as absent, which is right for a reader and wrong for
// a writer: the whole point of this gesture is to fill a slot whose average the
// build could not produce. So the store is reached directly -- it is a public
// member and bankSlotKey is a public static -- and grown to fit the slot.
static AnchoredBankSlot* anchoredSlotMut(TemplateBin& b, int lead, int slot,
    AnchorType a)
{
    if (lead < 0 || lead > 2 || slot < 0) return nullptr;
    std::vector<AnchoredBankSlot>& v =
        b.anchored_bank[TemplateBin::bankSlotKey(lead, a)];
    if ((int)v.size() <= slot) v.resize(slot + 1);
    return &v[slot];
}

// ---- ONE FRAME'S COLUMNS INTO ANOTHER'S ------------------------------------
//
// A whole-array shift, NaN-filled at the end it moves away from, which is the
// same thing align_beat_matrix does to a row and the same thing frameShift
// describes for a single column. Used to convert a re-stacked average out of
// the beat matrix's R frame and into the columns of the anchored slot it is
// about to replace.
static void shiftColumnsInPlace(std::vector<double>& v, int delta)
{
    if (delta == 0 || v.empty()) return;
    const int n = static_cast<int>(v.size());
    std::vector<double> out(v.size(),
        std::numeric_limits<double>::quiet_NaN());
    for (int i = 0; i < n; ++i) {
        const int j = i + delta;
        if (j >= 0 && j < n) out[static_cast<std::size_t>(j)] = v[static_cast<std::size_t>(i)];
    }
    v.swap(out);
}

bool TemplateViewerWindow::realignEcgFromBar(int binIdx, int leadIdx,
    int templateIdx, AnchorType anchor, double barCol, bool announce)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return false;
    if (leadIdx < 0 || leadIdx > 2) return false;
    if (anchor == AnchorType::R_PEAK) return false;   // the base frame, not an alignment
    TemplateBin& b = m_bins[binIdx];

    tbank::TemplateBank& bank = b.ecg_bank[leadIdx];
    if (templateIdx < 0 || templateIdx >= bank.size()) return false;
    tbank::BankTemplate& slot = bank.templates[templateIdx];

    // members_clean when it exists -- the set the average on screen was taken
    // over. `members` also holds the premature and Tukey-excluded beats the
    // record keeps but the average excludes, and re-admitting them would look
    // like the re-alignment had changed the morphology.
    const std::vector<uint32_t>& cohort = !slot.members_clean.empty()
        ? slot.members_clean : slot.members;
    if (cohort.empty()) {
        if (auto* sb = announce ? statusBar() : nullptr)
            sb->showMessage(tr("No beats recorded for this template; "
                "nothing to re-stack."), 5000);
        return false;
    }

    const char* chan = ecg_realign::channelKey(leadIdx);
    const ppg_realign::BinBeats& beats = beatsForBin(binIdx, chan);
    if (beats.empty()) {
        if (auto* sb = announce ? statusBar() : nullptr)
            sb->showMessage(tr("No per-beat %1 data for this bin - cannot "
                "re-stack (the bar was still moved).")
                .arg(QString::fromLatin1(chan)), 6000);
        return false;
    }

    // ---- INTO THE FRAME THE ROWS ARE IN --------------------------------
    //
    // barCol arrives as a DRAWN-frame column: moveEcgMarker clamps `placed`
    // against the drawn wall and stores it through setView, which adds
    // frameShift(grid -> owner). The rows here are R-framed slices -- the
    // matrix alignTemplatesFromCache hands to align_beat_matrix, framed on
    // blk.r_col -- so the bar converts by frameShift(grid -> R).
    //
    // NOT frameShift(owner -> R), which is only the same number when the grid
    // happens to be drawing this bar's own alignment.
    const int rCol = b.chFor(leadIdx, AnchorType::R_PEAK).r_col_raw;
    const double targetR = barCol
        + b.frameShift(leadIdx, currentGridAnchor(), AnchorType::R_PEAK);

    const ecg_realign::Result res = ecg_realign::realignAt(
        beats, cohort, targetR, rCol, m_sampleRate, anchor,
        ecg_realign::fenceHalfWinFor(m_sampleRate));
    if (!res.ok) {
        // REFUSED, AND THE OLD AVERAGE STANDS.
        if (auto* sb = announce ? statusBar() : nullptr)
            sb->showMessage(tr("%1 re-stack refused: %2.")
                .arg(anchor_view::label(anchor))
                .arg(QString::fromStdString(res.why)), 6000);
        return false;
    }

    AnchoredBankSlot* dst = anchoredSlotMut(b, leadIdx, templateIdx, anchor);
    if (!dst) return false;

    // CLIPPED TO THE EXTENT IT REPLACES, exactly as adoptPulsePair does and for
    // the same reason: the beat matrix is padded to the bin's longest slice
    // while the anchored slot is trimmed to where its own membership ends, so
    // emitting the full width stretches the panel's x-frame into the far tail.
    std::vector<double> tmpl = res.tmpl;
    std::vector<double> iqr = res.iqr;

    // ---- BACK OUT OF THE ROW FRAME, BEFORE ANYTHING ELSE ----------------
    //
    // realignAt worked on _beats.bin's slices, which are R-framed, and it put
    // every member's landmark on targetR -- an R-frame column. The slot it is
    // about to replace is read through slotView, which pairs it with
    // chFor(lead, anchor).r_col_raw: its columns are the ANCHOR'S, and the
    // panel plots them straight onto its axis while the BAR is converted into
    // the drawn frame from its own cell (barsForPanel). Storing an R-framed
    // array under the anchor's tag therefore displaced the whole waveform by
    // r_col(anchor) - r_col(R) from the bar that produced it: the trace slid
    // sideways on mouse-up and the landmark came to rest somewhere other than
    // under the bar the operator had just placed.
    //
    // ZERO ON EXACTLY THE BINS THAT ALREADY LOOKED RIGHT. frameShift is 0 when
    // either r_col is missing, and 0 when this anchor's alignment came out
    // byte-identical to R -- the NOT ALIGNED bins the [anchors] summary counts,
    // which ecg_realign.hpp's own header says is roughly half of them. That is
    // why the gesture appeared to work on some bins and to throw the waveform
    // off on others.
    const int frameDelta = b.frameShift(leadIdx, AnchorType::R_PEAK, anchor);
    shiftColumnsInPlace(tmpl, frameDelta);
    shiftColumnsInPlace(iqr, frameDelta);

    const std::size_t keep = dst->tmpl.size();
    if (keep > 0 && tmpl.size() > keep) {
        tmpl.resize(keep);
        if (iqr.size() > keep) iqr.resize(keep);
    }
    dst->tmpl = std::move(tmpl);
    // RAW-AMPLITUDE SD, ddof=1 -- what build_bins.hpp writes into this field
    // and what focusEcg draws through scale_array_by_ref (a plain divide by the
    // global reference, no foot involved). Unlike the pulse spread, nothing
    // divides this a second time, so there is no unit trap here.
    dst->tmpl_iqr = std::move(iqr);
    dst->n_members = static_cast<uint32_t>(res.n_used);

    // ---- REDRAW, WITHOUT REBUILDING THE PAGE ---------------------------
    //
    // The panels read this average through slotView every time they are
    // re-skinned, so re-anchoring the grid to the alignment on screen is enough
    // and leaves the layout, the focus view and any drag in progress alone.
    // A showPage() here would destroy the panel the gesture came from.
    //
    // SCOPED TO THIS COLUMN: one slot's average changed, and the unscoped pass
    // re-detects glyphs on every panel of every column on the page.
    reskinGridForAnchor(binIdx, templateIdx);
    // THE WAVEFORM CHANGED UNDER AN UNCHANGED FOCUS KEY -- the one case the
    // column-only fast path cannot detect for itself, since (bin, lead, slot,
    // anchor, marker, panel) all still name this view. Say so, or the refresh
    // below would move the crosshair on the OLD average.
    invalidateFocusFastPath();
    m_sdCacheValid = false;   // same reason: the msec model was built on it
    if (m_focusMarker >= 0 && m_focusBin == binIdx && m_focusSlot == templateIdx)
        refreshFocus(m_focusWidget, m_focusBin, m_focusLead, m_focusSlot,
            m_focusMarker, m_focusCol);

    // n_used vs n_members is how much of the cohort had this landmark findable
    // near the operator's column; a large gap means the correction disagrees
    // with most of the beats, which is worth seeing before trusting it.
    if (auto* sb = announce ? statusBar() : nullptr)
        sb->showMessage(tr("%1 re-stacked on col %2: %3 of %4 beats "
            "(median shift %5, max %6 samples). Cohort unchanged - not re-filtered.")
            .arg(anchor_view::label(anchor))
            .arg(res.anchor_col + frameDelta)   // the anchored slot's own columns
            .arg(res.n_used).arg(res.n_members)
            .arg(res.median_shift).arg(res.max_shift), 8000);
    return true;
}