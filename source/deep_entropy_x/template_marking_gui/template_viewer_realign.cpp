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

#include <algorithm>   // stable_sort, for grouping the release columns by bin

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
    if (!(m_notchFilterHz > 0.0)) return sig;   // 0 = the operator left Notch off
    if (sig.empty() || fs <= 0.0) return sig;

    std::vector<double> out =
        notch_filter(sig, m_notchFilterHz, fs);

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
bool TemplateViewerWindow::pulseTraceForSlot(tbank::template_of_all_signals& slot,
    std::vector<double>& outTrace,
    std::vector<double>& outIqr,
    double& outFootIdx)
{
    outTrace.clear();
    outIqr.clear();
    outFootIdx = -1.0;
    if (slot.tmpl.empty()) return false;

    // Seeded for every slot in seedOneBin, so there is nothing to detect here.
    outFootIdx = slot.pulse_marks.onset;

    const std::vector<double> src =
        maybeNotchTrace(slot.tmpl, m_ppgRateHz, outFootIdx);

    outTrace = normalize_ppg_or_similar(src, outFootIdx, 0);
    outIqr = normalize_features::scale_pulse_spread_by_ref(
        slot.tmpl_std,
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
static void adoptPulsePair(tbank::template_of_all_signals& slot,
    std::vector<double> tmpl, std::vector<double> iqr)
{
    const std::size_t keep = slot.tmpl.size();
    if (keep > 0 && tmpl.size() > keep) {
        tmpl.resize(keep);
        if (iqr.size() > keep) iqr.resize(keep);
    }
    slot.tmpl = std::move(tmpl);
    slot.tmpl_std = std::move(iqr);
    slot.band_lo.clear();
    slot.band_hi.clear();
    // THE MARKS ARE NOT TOUCHED HERE ANY MORE.
    //
    // This used to re-run seed_pulse_bank_template on the waveform it had just
    // installed, on the reasoning that the marks are derived from tmpl. They
    // are -- but from the variant they were measured on, not from whichever
    // variant happens to be installed for display. Re-detecting here put the
    // foot back on a peak-aligned average, where in a messy bin it landed half
    // a second away; that foot was then both the display's normalisation
    // reference and the next alignment's search hint, so one gesture moved the
    // whole trace and the next one moved it somewhere else again.
    //
    // Each variant's marks are detected once, on that variant, by
    // buildPulseVariant. composePulseMarks assembles the set the panel reads.
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

    (void)leadIdx;   // the pulse is per (bin, slot); leads share it

    // ---- ONLY THE FOOT RE-STACKS ---------------------------------------
    //
    // The foot bar IS the alignment landmark: saying where the foot is says
    // how the cohort should line up, so a re-level follows from it.
    //
    // THE DICROTIC NOTCH AND THE END ARE MEASUREMENTS ON the waveform, not
    // statements about how to build it. They briefly triggered a peak-aligned
    // re-level here, which was wrong twice over: the column they were dropped
    // on was never an input to the arithmetic, so the gesture was a button
    // wearing a bar's clothes; and the re-level was not repeatable, so the
    // average moved on every drag. Peak alignment belongs to the Peak
    // position of the align group, where it is chosen deliberately and once.
    if (marker != BinPlotWidget::PpgOnset) return;

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

    // ---- BY BIN, SO THE ONE-DEEP MATRIX CACHE HITS ---------------------
    //
    // m_dragPropCols is in PANEL order, which walks the page left to right and
    // therefore alternates bins; beatsForBin caches exactly one bin, so that
    // order missed on nearly every column. Sorting by bin makes every slot of
    // a bin consecutive, and the second and later slots of each bin cost
    // nothing. It matters most on the disk fallback, where a miss re-walks the
    // whole _beats.bin.
    //
    // STABLE, so slots within a bin keep their left-to-right order and the
    // status line still reports them in the order the operator sees.
    std::stable_sort(alsoPulse.begin(), alsoPulse.end(),
        [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            return a.first < b.first;
        });

    const bool manyPulse = !alsoPulse.empty();

    // relevelAt, not relevelAtOwnCrossing: the operator's column IS the answer
    // for this bar, and newIdx is where they dropped it.
    relevelPulseAtFoot(binIdx, templateIdx,
        static_cast<double>(newIdx), /*announce=*/!manyPulse);
    if (!manyPulse) return;

    QGuiApplication::setOverrideCursor(Qt::WaitCursor);
    for (const std::pair<int, int>& pc : alsoPulse) {
        const time_bin& tb = m_bins[pc.first];
        if (pc.second >= (int)tb.ppg_bank.size()) continue;
        // THE PROPAGATED COLUMN'S OWN BAR, which movePpgMarker wrote during
        // the drag -- so each column re-levels on the column Move-Subsequent
        // gave it, not on the dragged column's.
        const double foot =
            tb.ppg_bank.templates[pc.second].pulse_marks.onset;
        if (!(foot >= 0.0)) continue;
        relevelPulseAtFoot(pc.first, pc.second, foot, /*announce=*/false);
    }
    QGuiApplication::restoreOverrideCursor();


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
            // BORROWED, NOT COPIED. `out.rows = rows` here was a deep copy of
            // the whole bin -- megabytes and hundreds of allocations per
            // gesture, and per COLUMN under Move-Subsequent, because the cache
            // below is one deep and keyed by bin. per_channel_beats belongs to
            // the BeatsFile this window was handed and outlives the cache, so
            // pointing at it is safe; see BinBeats.
            out.external = &rows;
            for (const auto& row : rows)
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
// ---- THE AS-BUILT AVERAGE, AND THE ANCHOR MEASURED ON IT -----------------
//
// Every re-level reads its hint from slot.pulse_marks.onset, which is stored;
// this is kept for the spread band autoPctForSlot measures.
const std::vector<double>& TemplateViewerWindow::ppgAsBuiltTmpl(
    int binIdx, int templateIdx) const
{
    static const std::vector<double> kEmpty;
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return kEmpty;
    const time_bin& b = m_bins[binIdx];
    if (templateIdx < 0 || templateIdx >= (int)b.ppg_bank.size()) return kEmpty;

    // THE STASH WHEN THERE IS ONE. stashBuiltPulse keeps the FIRST overwrite
    // only, so its entry is the build's own pair however many re-levels have
    // run since; before the first one the live template still is it.
    const auto it = m_ppgBuilt.find(slotKey(binIdx, templateIdx));
    if (it != m_ppgBuilt.end()) return it->second.first;
    return b.ppg_bank.templates[templateIdx].tmpl;
}

// ---- LANDMARK -> VARIANT, ONE DEFINITION ---------------------------------
tbank::PulseAnchor TemplateViewerWindow::pulseVariantForMarker(int marker)
{
    return (marker == BinPlotWidget::PpgOnset)
        ? tbank::PulseAnchor::Foot
        : tbank::PulseAnchor::Peak;
}

// ---- pulse_marks <- THE TWO VARIANTS -------------------------------------
//
// The foot from Foot, everything else from Peak, per pulseVariantForMarker.
// A variant that has not been built contributes nothing rather than -1, so a
// half-built slot keeps whatever it already had on screen instead of losing
// its bars.
void TemplateViewerWindow::composePulseMarks(tbank::template_of_all_signals& slot)
{
    const tbank::PulseVariant& F =
        slot.pulseVariant(tbank::PulseAnchor::Foot);
    const tbank::PulseVariant& P =
        slot.pulseVariant(tbank::PulseAnchor::Peak);

    if (F.ok()) {
        slot.pulse_marks.onset = F.marks.onset;
        slot.pulse_marks.onset_auto = F.marks.onset_auto;
    }
    if (P.ok()) {
        slot.pulse_marks.dicrotic = P.marks.dicrotic;
        slot.pulse_marks.end = P.marks.end;
        slot.pulse_marks.dicrotic_auto = P.marks.dicrotic_auto;
        slot.pulse_marks.end_auto = P.marks.end_auto;
        // THE PEAK AND ITS DERIVATIVES COME FROM Peak, which is the alignment
        // that sharpens them: the apex is what its rows were levelled on.
        slot.pulse_marks.peak_auto = P.marks.peak_auto;
        slot.pulse_marks.peak2_auto = P.marks.peak2_auto;
        slot.pulse_marks.notch_found = P.marks.notch_found;
    }
}

bool TemplateViewerWindow::buildPulseVariant(int binIdx, int templateIdx,
    tbank::PulseAnchor v, double pct, bool announce)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return false;
    time_bin& b = m_bins[binIdx];
    if (templateIdx < 0 || templateIdx >= (int)b.ppg_bank.size()) return false;
    tbank::template_of_all_signals& slot = b.ppg_bank.templates[templateIdx];
    if (slot.tmpl.empty()) return false;

    // members_clean when it exists -- the set the build averaged. Re-averaging
    // `members` would re-admit the premature and Tukey-excluded beats and look
    // like the alignment had changed the morphology.
    const std::vector<uint32_t>& cohort = !slot.members_clean.empty()
        ? slot.members_clean : slot.members;
    if (cohort.empty()) {
        if (auto* sb = announce ? statusBar() : nullptr)
            sb->showMessage(tr("No pulse members recorded for this template; "
                "nothing to align."), 5000);
        return false;
    }

    // ---- THE HINT IS FIXED PER VARIANT ---------------------------------
    //
    // Foot: the operator's anchor column, which only a foot drag changes.
    // Peak: the as-built detected foot, which nothing changes.
    //
    // Neither is slot.pulse_marks.onset. That field is re-detected per variant
    // and reading it here is what made the old path iterate on its own output.
    // FOOT TAKES THE STORED BAR. Re-detecting it each session landed it a
    // sample out on the waveform bank_reload had restored, which moved
    // relevelAtOwnCrossing's search window, which moved out.baseline, which
    // offset every column of the result.
    //
    // PEAK DOES NOT, and must not: pulse_marks.onset is written by
    // composePulseMarks out of the Foot variant, so reading it here makes this
    // variant depend on whether Foot was built first. ppgAsBuiltTmpl is a pure
    // function of an array nothing mutates.
    const double hint = (v == tbank::PulseAnchor::Foot)
        ? slot.pulse_marks.onset
        : static_cast<double>(FeatureMarks::detect_ppg_onset(
            ppgAsBuiltTmpl(binIdx, templateIdx)));
    if (!(hint >= 0.0)) return false;

    // Peak is 100 by definition. Foot takes the height already resolved for
    // this slot when there is one, so a rebuild repeats it rather than asking
    // the threshold again; otherwise the caller's, kept for next time.
    double usePct = 100.0;
    if (v != tbank::PulseAnchor::Peak) {
        usePct = (slot.pulse_marks.foot_pct >= 0.0)
            ? slot.pulse_marks.foot_pct : std::clamp(pct, 0.0, 100.0);
        slot.pulse_marks.foot_pct = usePct;
    }

    if (beatsBinPath().isEmpty()) return false;
    const ppg_realign::BinBeats& beats = beatsForBin(binIdx);
    if (beats.empty()) return false;

    const ppg_realign::RelevelResult res = ppg_realign::relevelAtOwnCrossing(
        beats, cohort, hint,
        ppg_realign::searchHalfWinFor(m_ppgRateHz),
        ppg_realign::levelHalfWinFor(m_ppgRateHz), usePct);
    if (!res.ok) {
        if (auto* sb = announce ? statusBar() : nullptr)
            sb->showMessage(tr("Pulse alignment refused: %1.")
                .arg(QString::fromStdString(res.why)), 6000);
        // BUILT, AND EMPTY. The refusal is recorded so the next paint does not
        // retry it; PulseVariant::ok() is false either way.
        tbank::PulseVariant& out = slot.pulseVariant(v);
        out.built = true;
        out.tmpl.clear();
        out.tmpl_iqr.clear();
        return false;
    }

    tbank::PulseVariant& out = slot.pulseVariant(v);
    // CLIPPED TO THE SLOT'S OWN LENGTH here rather than in adoptPulsePair, so
    // both variants are stored at the same width and a column index means the
    // same instant in each -- which is the property that lets a foot read off
    // Foot sit on the same axis as a notch read off Peak.
    out.tmpl = res.tmpl;
    out.tmpl_iqr = res.iqr;
    const std::size_t keep = slot.tmpl.size();
    if (keep > 0 && out.tmpl.size() > keep) {
        out.tmpl.resize(keep);
        if (out.tmpl_iqr.size() > keep) out.tmpl_iqr.resize(keep);
    }

    // ---- THIS VARIANT'S OWN AUTO MARKS, ON THIS VARIANT ----------------
    //
    // Detected once, here, and never again -- so the foot detector never runs
    // on a peak-aligned average and the notch detector never runs on a
    // foot-aligned one. The operator's bars are NOT touched: a variant that
    // already carries a placed bar keeps it.
    tbank::BankPulseMarkerSet fresh;
    FeatureMarks::seed_pulse_bank_template(out.tmpl, m_ppgRateHz, fresh);
    const double keepOnset = out.marks.onset;
    const double keepDicr = out.marks.dicrotic;
    const double keepEnd = out.marks.end;
    out.marks = fresh;
    if (keepOnset >= 0.0) out.marks.onset = keepOnset;
    if (keepDicr >= 0.0) out.marks.dicrotic = keepDicr;
    if (keepEnd >= 0.0) out.marks.end = keepEnd;
    // FOOT'S BAR IS THE ANCHOR, not the detection: the operator's column is
    // what this variant was built to, so it is what its foot bar reads.
    if (v == tbank::PulseAnchor::Foot) out.marks.onset = hint;

    out.built = true;
    return true;
}

bool TemplateViewerWindow::showPulseVariant(int binIdx, int templateIdx,
    tbank::PulseAnchor v, bool announce)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return false;
    time_bin& b = m_bins[binIdx];
    if (templateIdx < 0 || templateIdx >= (int)b.ppg_bank.size()) return false;
    tbank::template_of_all_signals& slot = b.ppg_bank.templates[templateIdx];

    if (!slot.pulseVariant(v).built)
        buildPulseVariant(binIdx, templateIdx, v, percentage_for_aligning(), announce);
    const tbank::PulseVariant& pv = slot.pulseVariant(v);
    if (!pv.ok()) return false;

    // The build's own pair, kept once, so "as built" remains reachable and
    // ppgAsBuiltTmpl has something stable to measure hints against.
    stashBuiltPulse(binIdx, templateIdx, slot);

    adoptPulsePair(slot, pv.tmpl, pv.tmpl_iqr);
    composePulseMarks(slot);
    m_ppgRealigned.insert(slotKey(binIdx, templateIdx));
    pushPulseToPanels(binIdx, templateIdx);
    return true;
}

void TemplateViewerWindow::stashBuiltPulse(int binIdx, int templateIdx,
    const tbank::template_of_all_signals& slot)
{
    const int key = slotKey(binIdx, templateIdx);
    if (m_ppgBuilt.count(key)) return;   // FIRST overwrite only
    m_ppgBuilt[key] = { slot.tmpl, slot.tmpl_std };
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
    time_bin& b = m_bins[binIdx];
    if (templateIdx < 0 || templateIdx >= (int)b.ppg_bank.size()) return false;

    const auto it = m_ppgBuilt.find(slotKey(binIdx, templateIdx));
    if (it == m_ppgBuilt.end()) return false;   // never re-stacked: nothing to undo

    tbank::template_of_all_signals& slot = b.ppg_bank.templates[templateIdx];
    slot.tmpl = it->second.first;
    slot.tmpl_std = it->second.second;
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
    time_bin& b = m_bins[binIdx];
    if (templateIdx < 0 || templateIdx >= (int)b.ppg_bank.size()) return;
    tbank::template_of_all_signals& slot = b.ppg_bank.templates[templateIdx];

    std::vector<double> trace, iqr;
    double footIdx = -1.0;
    if (!pulseTraceForSlot(slot, trace, iqr, footIdx)) return;

    const std::vector<BinPlotWidget*>* col = panelsForColumn(binIdx, templateIdx);

    if (col) {
        for (auto* pw : *col) {
            if (!pw) continue;
            pw->setPpgData(trace, iqr, slot.memberCount());
            // AND THE THREE BARS, NOT JUST THE TRACE. adoptPulsePair has
            // re-seeded pulse_marks on the new waveform, but the panels still
            // hold the columns applyBankTemplateToWidget pushed BEFORE the
            // re-level -- showPage builds the panels and then calls
            // realignAllVisiblePulses at its tail, so on the initial page the
            // bars were a page older than the trace under them. It took a bar
            // click to reskin the column and pick the new cells up, which is
            // the whole "not until I click another bar" symptom.
            const tbank::BankPulseMarkerSet& pm = slot.pulse_marks;
            const FeatureMarks::ReactivePpg rp = FeatureMarks::reactive_ppg(
                slot.tmpl, pm.onset, pm.peak_auto, pm.dicrotic, pm.end);
            pw->setMarker(BinPlotWidget::PpgOnset, pm.onset);
            pw->setMarker(BinPlotWidget::PpgPeak, pm.peak_auto);
            pw->setMarker(BinPlotWidget::PpgDicrotic, pm.dicrotic);
            pw->setMarker(BinPlotWidget::PpgPeak2, rp.peak2);
            pw->setMarker(BinPlotWidget::PpgEnd, pm.end);
            pw->setMarker(BinPlotWidget::PpgT50, rp.t50);
            pw->setMarker(BinPlotWidget::PpgT80, rp.t80);
        }
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

int TemplateViewerWindow::which_alignment_fiducial_marker_should_auto_use(int binIdx, int templateIdx,
    const tbank::template_of_all_signals& slot, double footCol) const
{
    const auto it = m_ppgBuilt.find(slotKey(binIdx, templateIdx));
    const std::vector<double>& iqr = (it != m_ppgBuilt.end())
        ? it->second.second
        : slot.tmpl_std;

    const int halfwin = std::max(1, static_cast<int>(
        std::lround(region_around_foot_to_measure_std * m_ppgRateHz)));
    const double spread = ppg_realign::std_about(iqr, footCol, halfwin);

    // NO SPREAD MEASURABLE -> THE FOOT. An absent or all-zero band is not
    // evidence for the fallback: it means this slot's spread was never
    // written, or fewer than two beats reached the columns around the foot.
    // The foot is where the operator's bar is and what they would expect.
    if (!(spread >= 0.0)) return 0;

    return (spread < max_foot_std_to_trigger_looking_ahead) ? percent_up_upstroke_to_go_if_noisy : 0;
}

// ONE PLACE THAT CHANGES THE PULSE ALIGNMENT, for the same reason
// applyAlignmentSelection is the one place that changes the ECG one: the radio
// buttons, the spin box and the foot-drag override all land here, and each of
// them carrying its own copy of "set the members, re-stack the page, sync the
// controls" is three chances for one of them to leave the group checked on an
// alignment that is not on screen.
void TemplateViewerWindow::applyPpgAlignSelection(PpgAlign mode, int pct)
{
    m_ppgAlignMode = mode;
    m_ppgAlignPercent = std::clamp(pct, 0, 100);
    syncPpgAlignControls();
    realignAllVisiblePulses();
}

void TemplateViewerWindow::setPpgAlignMode(PpgAlign mode)
{
    m_ppgAlignMode = mode;
    syncPpgAlignControls();
}

// ---- BUILD BOTH VARIANTS, SHOW THE SELECTED ONE --------------------------
//
// This used to re-level slot.tmpl in place, once per column, per selection.
// Now each column's two alignments are computed from fixed inputs and stored,
// and the selection picks one. Re-selecting the same position therefore shows
// the identical waveform instead of computing a new one, which is what the
// operator means by "it is already peak aligned".
void TemplateViewerWindow::realignAllVisiblePulses()
{
    int nDone = 0, nSkipped = 0;
    int nAutoFoot = 0, nAutoPct = 0;

    for (int li = 0; li < (int)m_pageGlobalIdx.size()
        && li < (int)m_pageTemplateIdx.size(); ++li) {
        const int gi = m_pageGlobalIdx[li];
        const int slotIdx = m_pageTemplateIdx[li];
        if (gi < 0 || gi >= (int)m_bins.size() || slotIdx < 0) continue;

        time_bin& b = m_bins[gi];
        if (slotIdx >= (int)b.ppg_bank.size()) continue;
        tbank::template_of_all_signals& slot = b.ppg_bank.templates[slotIdx];
        if (slot.tmpl.empty()) continue;
        // This column's own verdict: a column already called bad is left as it
        // is. The bin-level flag reaches these bits through the load seed.
        if (slot.badPulseMarked() || slot.marked_invalid_template) continue;

        // ---- THE PERCENTAGE _F IS BUILT AT ----------------------------
        //
        // Auto decides per column (autoPctForSlot: the foot, or
        // kAutoFallbackPct up the upstroke where the foot's own spread is too
        // tight to level on). Foot / Percent / Peak read the group through
        // pctForAlignMode. Either way it is _F's definition -- _P is 100 and
        // reads nothing.
        double pctF = percentage_for_aligning();
        if (slot.pulse_marks.foot_pct >= 0.0) {
            pctF = slot.pulse_marks.foot_pct;   // resolved once, then kept
        }
        else if (m_ppgAlignMode == PpgAlign::Auto) {
            pctF = static_cast<double>(which_alignment_fiducial_marker_should_auto_use(
                gi, slotIdx, slot, slot.pulse_marks.onset));
            if (pctF > 0.0) ++nAutoPct; else ++nAutoFoot;
        }

        // BOTH, ALWAYS. The composed marker set needs Peak's notch and end
        // even when Foot is the trace on screen, and Peak costs nothing after
        // the first build because nothing invalidates it.
        buildPulseVariant(gi, slotIdx, tbank::PulseAnchor::Foot, pctF);
        buildPulseVariant(gi, slotIdx, tbank::PulseAnchor::Peak, 100.0);

        if (showPulseVariant(gi, slotIdx, pulseVariantToShow())) ++nDone;
        else ++nSkipped;
    }
}

// The variant the page draws. Forced positions say it outright; Auto follows
// the last pulse bar the operator clicked (m_ppgViewVariant), which starts at
// Foot.
tbank::PulseAnchor TemplateViewerWindow::pulseVariantToShow() const
{
    switch (m_ppgAlignMode) {
    case PpgAlign::Peak: return tbank::PulseAnchor::Peak;
    case PpgAlign::Foot:
    case PpgAlign::Percent: return tbank::PulseAnchor::Foot;
    case PpgAlign::Auto:
    default: return m_ppgViewVariant;
    }
}

// The operator clicked a pulse bar. In Auto, show the average that bar was
// measured on -- the ECG side does the same with m_autoGridAnchor, and for the
// same reason. Returns whether the view changed.
bool TemplateViewerWindow::followPulseBarSelection(int marker)
{
    if (m_ppgAlignMode != PpgAlign::Auto) return false;
    if (!BinPlotWidget::markerIsPpg(marker)) return false;
    // GLYPHS LEAVE IT ALONE. t50, t80, peak and peak2 are crossings and
    // detections, not a statement about which alignment is being read.
    if (marker != BinPlotWidget::PpgOnset
        && marker != BinPlotWidget::PpgDicrotic
        && marker != BinPlotWidget::PpgEnd) return false;

    const tbank::PulseAnchor want = pulseVariantForMarker(marker);
    if (want == m_ppgViewVariant) return false;
    m_ppgViewVariant = want;
    realignAllVisiblePulses();
    return true;
}

// ---- THE FOOT DRAG: REDEFINE _F, AND SHOW IT -----------------------------
//
// The operator's column becomes the anchor _F is built from, for this slot,
// permanently -- so dragging the foot away and back returns the original
// waveform exactly. _P is not touched: it has no bar among its inputs.
//
// THE VIEW FOLLOWS. Having just placed the foot, the operator is looking at
// the foot, so the foot-aligned average is what they are shown.
bool TemplateViewerWindow::relevelPulseAtFoot(int binIdx, int templateIdx,
    double footCol, bool announce)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return false;
    time_bin& b = m_bins[binIdx];
    if (templateIdx < 0 || templateIdx >= (int)b.ppg_bank.size()) return false;
    if (b.ppg_bank.templates[templateIdx].tmpl.empty()) return false;
    if (!(footCol >= 0.0)) return false;

    // No anchor memo: buildPulseVariant reads slot.pulse_marks.onset, which
    // this drag has already written, so the bar IS the anchor.
    // A NEW HEIGHT TOO -- the operator moved the foot, so the alignment that
    // was resolved against the old one is stale.
    if (binIdx >= 0 && binIdx < static_cast<int>(m_bins.size())
        && templateIdx >= 0
        && templateIdx < static_cast<int>(m_bins[binIdx].ppg_bank.size()))
        m_bins[binIdx].ppg_bank.templates[templateIdx].pulse_marks.foot_pct = -1.0;

    // REBUILT, NOT PATCHED. The stored _F is a function of the anchor, so a
    // new anchor means a new _F -- computed from the as-built rows, never from
    // the _F it replaces.
    tbank::PulseVariant& fv =
        b.ppg_bank.templates[templateIdx].pulseVariant(tbank::PulseAnchor::Foot);
    fv.built = false;
    if (!buildPulseVariant(binIdx, templateIdx, tbank::PulseAnchor::Foot,
        percentage_for_aligning(), announce)) return false;

    if (m_ppgAlignMode == PpgAlign::Auto)
        m_ppgViewVariant = tbank::PulseAnchor::Foot;
    return showPulseVariant(binIdx, templateIdx, pulseVariantToShow(), announce);
}