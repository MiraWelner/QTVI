// ========================================================================
// template_viewer_focus.cpp
// The focus panels: landmark close-ups for ECG and pulse.
//
// One of five units; TemplateViewerWindow is declared in template_viewer.hpp.
// ========================================================================

#include "template_viewer.hpp"


// Both panels to "nothing selected".
void TemplateViewerWindow::clearFocusPanels() {
    setFocusSplit(false);
    if (zoomed_in_section_top) zoomed_in_section_top->clearFocus();
    if (zoomed_in_section_bottom) zoomed_in_section_bottom->clearFocus();
}

static std::vector<double> get_savitzky_golay_derivative_for_every_sample_in_vector(const std::vector<double>& t, double fs) {
    /*to get the standard deviation for the focus panel, the variance is divided by the savitsky golay derivative.
    This function takes in a vector and outputs the vector's savitsky golay derivative for each sample*/
    const int N = (int)t.size();
    std::vector<double> s(N, std::numeric_limits<double>::quiet_NaN());
    const int h = (fs > 0.0) ? std::max(1, (int)std::lround(4.0 * fs / 1000.0)) : 1;   // +/- 4 ms
    const int order = std::min(4, 2 * h);          // need 2h+1 >= order+1
    const std::vector<double> k = ppg_deriv::sgCoeffs(h, order, /*deriv=*/1);

    for (int i = h; i < N - h; ++i) {
        double d = 0.0;
        bool ok = true;
        for (int j = -h; j <= h; ++j) {
            const double v = t[i + j];
            if (std::isnan(v)) { ok = false; break; }
            d += k[j + h] * v;
        }
        if (ok) s[i] = std::fabs(d);   // per-sample derivative (deriv! = 1)
    }
    return s;
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

// ---- PULSE AND ARTERIAL FOCUS ------------------------------------------
//
// These ride their own template (ppgTemplate / abpTemplate / artTemplate /
// artPulmTemplate) with the matching per-sample std (*_iqr, ddof=1) and the
// shared pulse beat count (ppg_n_beats -- all pulse channels derive from the
// same foot-anchored beat set). One panel, not two: pulse channels are
// foot-anchored once and have no QRS/JT split and no alignment dimension.
//
// Split out of refreshFocus, which was 500 lines doing three unrelated jobs.
void TemplateViewerWindow::focusPulse(TemplateBin& b, int templateIdx,
    int marker, double col)
{
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
            clearFocusPanels();
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
        zoomed_in_section_top->setFocus(mean, sd, nBeats, col, chLabel + " " + pulseLabel(marker));
    zoomed_in_section_top->setFitKind(FocusPanelWidget::FitKind::Transition);
    return;
}

void TemplateViewerWindow::refreshFocus(int binIdx, int leadIdx,
    int templateIdx, int marker, double col)
{
    if (!zoomed_in_section_top) return;   // panels not created (nothing to do)
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    // Remember this focus so a fit-mode radio change can replay it in place.
    m_focusBin = binIdx; m_focusLead = leadIdx; m_focusSlot = templateIdx;
    m_focusMarker = marker; m_focusCol = col;
    TemplateBin& b = m_bins[binIdx];

    // Pulse channels have no alignment dimension, so they take a separate,
    // self-contained path (focusPulse) rather than threading through the ECG
    // slot and anchor selection below.
    if (!BinPlotWidget::markerIsEcg(marker)) {
        focusPulse(b, templateIdx, marker, col);
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

    // A BAR TAKES ITS OWN ALIGNMENT; A GLYPH TAKES THE ONE ON SCREEN.
    //
    // anchorFor returns R_PEAK for every glyph -- deliberately, since a glyph
    // is measured on all four averages and has no alignment of its own. Using
    // that as the FOCUS alignment meant clicking the P peak showed the
    // R-aligned average while clicking the P-onset bar a few pixels away showed
    // the P-aligned one. The glyph was drawn on the waveform currently
    // displayed, so that is the waveform its close-up must magnify.
    //
    // Bars keep anchorFor rather than currentGridAnchor() because of ordering:
    // user_clicked_on_bar calls refreshFocus BEFORE it moves m_autoGridAnchor
    // and re-skins, so at this point currentGridAnchor() is still the alignment
    // being left. A glyph click moves no anchor, so for glyphs it is current.
    const AnchorType focusAnchor = m_forceAlign
        ? m_forcedAlign
        : (anchor_view::isBar(marker) ? anchor_view::anchorFor(marker)
            : currentGridAnchor());    // NO FALLBACK: chForStrict returns nullptr when this alignment is absent
    // from the file, and the panel is cleared rather than showing the R-aligned
    // average under this bar's header. chFor's fallback did the latter, which
    // is what made every bar look identical.
    const ChannelTemplateData* chP = b.chForStrict(leadIdx, focusAnchor);
    if (!chP) {
        clearFocusPanels();
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
    const uint64_t nb[3] = { b.ch1_n_beats_raw, b.ch2_n_beats_raw, b.ch3_n_beats_raw };
    int nBeats = static_cast<int>(nb[leadIdx]);
    if (templateIdx > 0) {
        const tbank::TemplateBank& bank = b.ecg_bank[leadIdx];
        if (templateIdx >= bank.size()
            || bank.templates[templateIdx].tmpl.empty()) {
            clearFocusPanels();
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
        // Per-slot aligned averages are computed in alignTemplatesFromCache,
        // as a reduction over the beats it already aligns.
        //
        // NO FALLBACK. Showing the slot's R-aligned average under a header
        // naming another alignment is the defect this branch existed to
        // label, and a label is not a fix. There is one format version and
        // every section is written unconditionally (template_io.cpp), so a
        // null here is not an old file -- it is build_templates failing to
        // write the per-slot average for this anchor, which is a writer bug to
        // go and fix, not a state to render.
        const AnchoredBankSlot* asl =
            b.bankSlotFor(leadIdx, templateIdx, focusAnchor);
        if (!asl) {
            clearFocusPanels();
            fprintf(stderr, "[focus] bin=%d lead=%d slot=%d anchor=%s"
                " NO PER-SLOT ALIGNED AVERAGE -- build_templates did not write"
                " it\n", binIdx, leadIdx, templateIdx,
                anchor_view::label(focusAnchor));
            return;
        }
        meanRawEcg = &asl->tmpl;
        sdRawEcg = &asl->tmpl_iqr;
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
        case BinPlotWidget::EcgQPeak:  return QStringLiteral("Q peak");
        case BinPlotWidget::EcgTPeak:  return QStringLiteral("T peak");
        case BinPlotWidget::EcgSEnd:   return QStringLiteral("J-point");
        case BinPlotWidget::EcgTEnd:   return QStringLiteral("T end");
        }
        return QStringLiteral("landmark");
        };

    // `col` arrived in the R frame -- the widget's, because the grid is
    // R-aligned. This panel plots THIS alignment's average, so the bar has to
    // be placed in its columns. Same conversion and same direction as the drag
    // path in onMarkerMoved, so a bar dragged on the grid lands under the
    // crosshair here -- and, crucially, the panel FOLLOWS the bar as the
    // operator drags it, because colHere is the live bar position (col), not a
    // frozen auto-detected column.
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

    if (zoomed_in_section_top) //if the zoomed in top section is activated (it will always be with any focus)                                                                                                                                                   
    {

        const std::vector<double> absSlope = get_savitzky_golay_derivative_for_every_sample_in_vector(mean, m_sampleRate);
        // Floor at 5% of the template's OWN max slope: the divide can't blow up
        // in flat regions, and the threshold scales with the waveform instead
        // of being a fixed amplitude constant.
        double maxSlope = 0.0;
        for (double s : absSlope) if (std::isfinite(s) && s > maxSlope) maxSlope = s;
        const double floor = 0.05 * maxSlope;
        const double msPerSample = (m_sampleRate > 0.0) ? 1000.0 / m_sampleRate : 0.0;
        const double NaNv = std::numeric_limits<double>::quiet_NaN();
        std::vector<double>  sdMs(sd.size(), NaNv);
        std::vector<uint8_t> floorMask(sd.size(), 0u);
        for (size_t k = 0; k < sd.size() && k < absSlope.size(); ++k) {
            if (std::isnan(sd[k]) || !std::isfinite(absSlope[k])) continue;
            double slope = absSlope[k];
            if (slope < floor) { floorMask[k] = 1u; continue; }   // leaves NaN
            sdMs[k] = sd[k] / slope * msPerSample;
        }

        // (The sd at the bar is no longer read here: FocusPanelWidget already
        //  holds m_sdMs and prints it on its own two lines below the plot, so a
        //  copy in the header would be the same number twice.)

        // Two states only: this alignment's own template, or -- on a
        // sub-template column whose file predates the per-slot section -- the
        // slot's unaligned average, said plainly. chForStrict has already
        // returned early if the alignment itself is absent, so the label can
        // never name an alignment the data is not in.
        const QString tag = QStringLiteral(" [%1-aligned]")
            .arg(QString::fromLatin1(anchor_view::label(focusAnchor)));

        // THE J POINT GETS BOTH PANELS. It is the one landmark bounding two
        // segments, so a single framing always hides one neighbourhood: above
        // it is framed to the RIGHT edge, where it ends the QRS; below to the
        // LEFT, where it starts the JT. Same waveform in both -- differing
        // only in which side of the landmark is shown.
        // THE LANDMARK AND ITS ALIGNMENT, nothing else. The sd moved to the
        // panel's own footer lines; the model name moved to a second header
        // line that the panel draws itself.
        const QString head = QString("%1%2").arg(labelFor(marker), tag);
        // EXACT transition candidates: recompute the detector's fits on the
        // SAME displayed average. detect_template_landmarks is scale- and
        // position-invariant (the BIC argmin and the landmark positions do not
        // move under the eref amplitude scale), so this reproduces the winner
        // and positions that placed the mark, and the curves overlay `mean`.
        // Only the transition bars carry candidates; peaks fit locally in the
        // panel and ignore an invalid set.
        subsample_refine::TransitionCandidates transCand;
        // The detector's own position for the focused landmark, in `mean`'s
        // columns. -1 until the block below supplies it.
        double detFid = -1.0;
        // Hoisted out of the block below, where lm lives: the fit-kind choice
        // needs it and lm does not reach that far. True by default so a path
        // that never detects keeps the old behaviour.
        bool qOnsetFound = true;
        {
            // Peaks need this block as well: it is where the detector runs on
            // `mean`, and its peak fields are the only way the panel can know
            // where the mark actually is. The transition CANDIDATES are
            // selected by the `marker` switch below; the `isTrans` flag that
            // used to be computed here was never read -- it stopped gating the
            // detection and nothing replaced its use.
            if (!mean.empty()) {
                // The detector's fit depends on the TEMPLATE, not on where the
                // operator drags the bar -- so it is identical on every mouse-
                // move of the same landmark. refreshFocus fires on each drag
                // move; recomputing detect_template_landmarks every time was
                // the drag lag. Cache by (bin,slot,lead,marker,anchor) and only
                // re-detect when the focused landmark actually changes.
                const long long tkey =
                    (((((static_cast<long long>(binIdx) * 64 + templateIdx) * 4
                        + leadIdx) * 32 + marker) * 8)
                        + static_cast<int>(focusAnchor));
                if (tkey == m_lastTransKey) {
                    transCand = m_lastTransCand;
                }
                else {
                    // Seed with the DISPLAYED alignment's own R column -- the
                    // same one the main-plot glyphs use (chFor(frame).r_col_raw)
                    // -- NOT r_peak_ch + frameShift, a computed value that
                    // drifted and put the P/T detection window on the wrong part
                    // of the wave. detect runs on the displayed (e.g. P-aligned)
                    // mean, which is correct; only the R seed was off.
                    const int rColInMean = std::clamp(
                        b.chFor(leadIdx, focusAnchor).r_col_raw,
                        0, static_cast<int>(mean.size()) - 1);
                    const FeatureMarks::TemplateLandmarks lm =
                        FeatureMarks::detect_template_landmarks(mean, rColInMean, m_sampleRate,
                            curve_fit::FitMode::Auto, curve_fit::PeakFitMode::Auto);
                    qOnsetFound = lm.q_onset_found;
                    switch (marker) {
                    case BinPlotWidget::EcgPBegin: transCand = lm.p_begin_cand; break;
                    case BinPlotWidget::EcgQBegin: transCand = lm.q_onset_cand; break;
                    case BinPlotWidget::EcgSEnd:   transCand = lm.s_end_cand;   break;
                    case BinPlotWidget::EcgTEnd:   transCand = lm.t_end_cand;   break;
                    default: break;
                    }
                    // The expensive half, on the same key. ecgDetect pairs the
                    // trace with its own R column, so this panel cannot pair
                    // them differently from the seeding or the grid glyphs.
                    m_lastDet = ecgDetect(b, leadIdx, templateIdx, focusAnchor,
                        m_sampleRate, m_onOffsetFitMode, m_peakFitMode);
                    m_lastTransKey = tkey;
                    m_lastTransCand = transCand;
                }

                // ---- THE FIDUCIAL, RE-BRACKETED EVERY CALL ---------------
                //
                // Outside the cache: the peaks follow the bars, which move
                // without changing the key.
                //
                // userMarks, NOT slotMarks. The P peak is bracketed by the
                // P-onset and Q-onset bars, which live in DIFFERENT anchors'
                // sets -- one anchor's set holds only the bars it admits, so
                // the bracket arrived half-empty and compute_p_peak clamped it
                // to the trace edge, putting the fiducial outside the window.
                // userMarks assembles all four and translates them into
                // focusAnchor's columns: the same set the grid draws its X
                // from. It is also the const accessor -- slotMarks' non-const
                // overload inserts through operator[].
                const EcgFiducials fid = ecgFiducialsFrom(m_lastDet,
                    m_sampleRate, m_peakFitMode,
                    b.userMarks(leadIdx, templateIdx, focusAnchor));
                switch (marker) {
                case BinPlotWidget::EcgPPeak:  detFid = fid.p_peak;  break;
                    // T PEAK WAS MISSING FROM THIS SWITCH, so the T-peak focus
                    // never received a fiducial at all and the panel fell back to
                    // the bar column for its dotted line.
                case BinPlotWidget::EcgTPeak:  detFid = fid.t_peak;  break;
                case BinPlotWidget::EcgQPeak:  detFid = fid.q_peak;  break;
                case BinPlotWidget::EcgRPeak:  detFid = fid.r_peak;  break;
                case BinPlotWidget::EcgPBegin: detFid = fid.p_begin; break;
                case BinPlotWidget::EcgQBegin: detFid = fid.q_onset; break;
                case BinPlotWidget::EcgSEnd:   detFid = fid.s_end;   break;
                case BinPlotWidget::EcgTEnd:   detFid = fid.t_end;   break;
                default: break;
                }
                // The candidates are computed in Auto (a stable window). The
                // radio only recolors WHICH is green -- override the winner
                // index here, without re-fitting, so switching the model never
                // reshapes the curves (that was the "sigmoid goes flat" bug).
                if (transCand.valid) {
                    switch (m_onOffsetFitMode) {
                    case curve_fit::FitMode::Linear:      transCand.winner = 0; break;
                    case curve_fit::FitMode::Sigmoid:     transCand.winner = 1; break;
                    case curve_fit::FitMode::FracPoly:    transCand.winner = 2; break;
                    case curve_fit::FitMode::CubicSpline: transCand.winner = 3; break;
                    case curve_fit::FitMode::Cubic:       transCand.winner = 4; break;
                    case curve_fit::FitMode::Auto:
                    default:                              break;   // keep BIC winner
                    }
                }
            }
        }

        if (marker == BinPlotWidget::EcgSEnd) {
            setFocusSplit(true);
            if (zoomed_in_section_top) {
                zoomed_in_section_top->setFocus(mean, sd, nBeats, colHere, head + QStringLiteral("  (QRS)"), 100, -1);
                zoomed_in_section_top->setFitKind(FocusPanelWidget::FitKind::Transition);
                zoomed_in_section_top->setTransitionCandidates(transCand);
                zoomed_in_section_top->setDetectorFiducial(detFid);
                zoomed_in_section_top->setSdMs(sdMs, floorMask, absSlope, floor);
            }
            if (zoomed_in_section_bottom) {
                zoomed_in_section_bottom->setFocus(mean, sd, nBeats, colHere, head + QStringLiteral("  (JT)"), 100, +1);
                zoomed_in_section_bottom->setFitKind(FocusPanelWidget::FitKind::Transition);
                zoomed_in_section_bottom->setTransitionCandidates(transCand);
                zoomed_in_section_bottom->setDetectorFiducial(detFid);
                zoomed_in_section_bottom->setSdMs(sdMs, floorMask, absSlope, floor);
            }
        }
        else {
            // Peaks (R/P/Q/T) draw their tested candidates; onsets/offsets draw
            // theirs. peakSigma MUST match the detector's per-landmark sigma so
            // the drawn quadratic/cubic are the same fits that placed the mark.
            const bool isPeak =
                (marker == BinPlotWidget::EcgRPeak
                    || marker == BinPlotWidget::EcgPPeak
                    || marker == BinPlotWidget::EcgQPeak
                    || marker == BinPlotWidget::EcgTPeak);
            // NO Q TROUGH => NO FIT. compute_q_onset falls back to its
            // R-upstroke branch on a monophasic-R beat, and lm.q_onset_found is
            // false -- the same flag the grid uses to draw a CIRCLE there
            // instead of an X, and the reason lm.q_peak comes back -1. Neither
            // the Q onset nor the Q peak was placed by a curve, so the panel
            // must not draw candidates for them.
            const bool qUnfound =
                (marker == BinPlotWidget::EcgQBegin
                    || marker == BinPlotWidget::EcgQPeak)
                && !qOnsetFound;
            const FocusPanelWidget::FitKind fk = qUnfound
                ? FocusPanelWidget::FitKind::None
                : (isPeak ? FocusPanelWidget::FitKind::PeakQuadratic
                    : FocusPanelWidget::FitKind::Transition);
            // Only read on the peak branch below, so the initialisers are
            // unreachable placeholders rather than a fallback for a real
            // landmark -- every peak marker sets both explicitly.
            double peakSigma = 0.0;
            int peakHalfWidth = 0;
            switch (marker) {
            case BinPlotWidget::EcgRPeak:
                peakSigma = subsample_refine::peak_sigma::R;
                peakHalfWidth = subsample_refine::peak_halfwidth::R; break;
            case BinPlotWidget::EcgPPeak:
                peakSigma = subsample_refine::peak_sigma::P;
                peakHalfWidth = subsample_refine::peak_halfwidth::P; break;
            case BinPlotWidget::EcgTPeak:
                peakSigma = subsample_refine::peak_sigma::T;
                peakHalfWidth = subsample_refine::peak_halfwidth::T; break;
            case BinPlotWidget::EcgQPeak:
                peakSigma = subsample_refine::peak_sigma::Q;
                peakHalfWidth = subsample_refine::peak_halfwidth::Q; break;
            default: break;
            }
            setFocusSplit(false);
            if (zoomed_in_section_bottom) zoomed_in_section_bottom->clearFocus();
            if (zoomed_in_section_top) {
                zoomed_in_section_top->setFocus(mean, sd, nBeats, colHere, head, 100, 0);
                // Peaks carry a sigma and a window; transitions have neither,
                // and the one-argument overload is the one that says so.
                if (fk == FocusPanelWidget::FitKind::PeakQuadratic
                    || fk == FocusPanelWidget::FitKind::PeakCubic)
                    zoomed_in_section_top->setPeakFitKind(fk, peakSigma,
                        peakHalfWidth);
                else
                    zoomed_in_section_top->setFitKind(fk);
                zoomed_in_section_top->setPeakFitMode(m_peakFitMode);
                zoomed_in_section_top->setTransitionCandidates(transCand);   // invalid for peaks -> ignored
                zoomed_in_section_top->setDetectorFiducial(detFid);
                zoomed_in_section_top->setSdMs(sdMs, floorMask, absSlope, floor);
            }
        }
    }
}