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

QString TemplateViewerWindow::beatsBinPath() const {
    // SAME DIRECTORY AND STEM morphology_csv::set was handed in
    // analysis_job::prepare, which is what wrote the file. Derived on each call
    // rather than cached at load, so it cannot survive a subject change and
    // point at the previous record's beats.
    if (m_templateDir.isEmpty() || m_subjectId.isEmpty()) return {};
    return m_templateDir + "/" + m_subjectId + "_beats.bin";
}

// Release of a dragged bar. ONE BAR MATTERS HERE -- the pulse foot. Every other
// marker's consequences are reactive and were already applied during the drag,
// so this returns immediately for them rather than doing work per gesture.
void TemplateViewerWindow::onMarkerReleasedOnTemplate(int binIdx, int leadIdx,
    int templateIdx, int marker, int newIdx)
{
    (void)leadIdx;   // the pulse is per (bin, slot); leads share it
    if (marker != BinPlotWidget::PpgOnset) return;
    realignPulseFromFoot(binIdx, templateIdx, static_cast<double>(newIdx));
}

void TemplateViewerWindow::realignPulseFromFoot(int binIdx, int templateIdx,
    double footCol)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    TemplateBin& b = m_bins[binIdx];
    if (templateIdx < 0 || templateIdx >= (int)b.ppg_bank.size()) return;
    tbank::BankTemplate& slot = b.ppg_bank.templates[templateIdx];
    if (slot.tmpl.empty()) return;

    // THE COHORT: members_clean when it exists, because that is the set the
    // waveform on screen was actually averaged over -- `members` includes the
    // premature and Tukey-excluded beats that were deliberately kept in the
    // record but out of the average. Re-averaging `members` would quietly
    // re-admit them, which would look like the re-alignment had changed the
    // morphology.
    const std::vector<uint32_t>& cohort = !slot.members_clean.empty()
        ? slot.members_clean : slot.members;
    if (cohort.empty()) {
        if (auto* sb = statusBar())
            sb->showMessage(tr("No pulse members recorded for this template; "
                "nothing to re-stack."), 5000);
        return;
    }

    const QString path = beatsBinPath();
    if (path.isEmpty()) return;

    // READ PER GESTURE, NOT CACHED. One bin's pulse rows is a few MB at most
    // and the read is a single seeking pass (ppg_realign::loadBin); holding
    // every bin's matrix for the life of the window would cost hundreds of MB
    // to serve a gesture that happens a handful of times per record. If this
    // ever shows up as a delay, cache the LAST bin read -- not all of them.
    const ppg_realign::BinBeats beats =
        ppg_realign::loadBin(path.toStdString(), static_cast<uint32_t>(binIdx));
    if (beats.empty()) {
        // The usual cause is a run with morphology_csv::setWriteBeatsBin(false),
        // which truncates the file on purpose. Say which file, because the
        // alternative is a control that silently does nothing.
        if (auto* sb = statusBar())
            sb->showMessage(tr("No per-beat pulse data in %1 - cannot re-stack "
                "(the bar was still moved).").arg(path), 6000);
        return;
    }

    const ppg_realign::Result res = ppg_realign::realign(
        beats, cohort, footCol, ppg_realign::searchHalfWinFor(m_ppgRateHz));
    if (!res.ok) {
        // REFUSED, AND THE OLD TEMPLATE STANDS. A re-median over one or two
        // beats is noise, and replacing a real template with it while reporting
        // success is the worst available outcome.
        if (auto* sb = statusBar())
            sb->showMessage(tr("Pulse re-stack refused: %1.")
                .arg(QString::fromStdString(res.why)), 6000);
        return;
    }

    // ---- ADOPT THE NEW WAVEFORM ---------------------------------------
    //
    // The bar is NOT moved: the operator put it where the foot is and the
    // re-stack anchored every beat to that column, so the trace has come to
    // the bar rather than the reverse. pulse_marks.onset was already written by
    // the drag handler (movePpgMarker), which is why nothing here touches it.
    slot.tmpl = res.tmpl;
    slot.tmpl_iqr = res.iqr;
    // The corridor described the previous stacking. Cleared rather than kept,
    // so nothing draws a band or scores a match against a shape that is no
    // longer there; it is rebuilt by the next recomputeTemplate that runs.
    slot.band_lo.clear();
    slot.band_hi.clear();

    m_ppgRealigned.insert(slotKey(binIdx, templateIdx));

    // ---- PUSH IT INTO THE LIVE PANELS, IN PLACE -----------------------
    //
    // setPpgData, not showPage(): a page rebuild would destroy the panel this
    // gesture came from, drop the focus view, and need a captureCurrentPage
    // first to avoid losing marker edits. The pulse average is the only thing
    // that changed, and every panel in this column draws the same one.
    std::vector<double> trace, iqr;
    double footIdx = -1.0;
    if (!pulseTraceForSlot(slot, trace, iqr, footIdx)) return;

    if (const std::vector<BinPlotWidget*>* col =
        panelsForColumn(binIdx, templateIdx)) {
        for (auto* pw : *col)
            if (pw) pw->setPpgData(trace, iqr, slot.memberCount());
    }

    // The focus panel, if it is showing this pulse, is re-run from the new
    // trace rather than left displaying the old stacking magnified.
    if (m_focusMarker >= 0 && m_focusBin == binIdx && m_focusSlot == templateIdx)
        refreshFocus(m_focusWidget, m_focusBin, m_focusLead, m_focusSlot,
            m_focusMarker, m_focusCol);

    // SAID OUT LOUD, including the part that is a compromise. n_used vs
    // n_members is how much of the cohort had a foot near the new column: a
    // large gap means the correction disagrees with most of the beats, which is
    // worth seeing before trusting the result.
    if (auto* sb = statusBar())
        sb->showMessage(tr("Pulse re-stacked on foot col %1: %2 of %3 beats "
            "(median shift %4, max %5 samples). Cohort unchanged - not re-filtered.")
            .arg(res.anchor_col)
            .arg(res.n_used).arg(res.n_members)
            .arg(res.median_shift).arg(res.max_shift), 8000);
}
