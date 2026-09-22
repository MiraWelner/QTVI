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

// The sd-in-msec model: each column's amplitude SD divided by the template's
// own |dV/dt| there, floored at 5% of that template's max slope so the divide
// cannot blow up in a flat region and the threshold scales with the waveform
// instead of being an absolute constant.
//
// ONE COPY, BOTH CHANNELS. This was inline in the ECG branch, which is the
// whole reason the pulse panel printed "sd = --": there was nothing wrong with
// the model on a pulse, it just lived somewhere a pulse could not reach.
//
// SdMsModel is declared in template_viewer.hpp, because focusEcg caches one --
// which is also why this is a static member rather than a file-scope function.
TemplateViewerWindow::SdMsModel TemplateViewerWindow::sd_in_msec(
    const std::vector<double>& mean, const std::vector<double>& sd, double fs)
{
    SdMsModel m;
    m.absSlope =
        get_savitzky_golay_derivative_for_every_sample_in_vector(mean, fs);
    double maxSlope = 0.0;
    for (double v : m.absSlope)
        if (std::isfinite(v) && v > maxSlope) maxSlope = v;
    m.floor = 0.05 * maxSlope;

    const double msPerSample = (fs > 0.0) ? 1000.0 / fs : 0.0;
    const double NaNv = std::numeric_limits<double>::quiet_NaN();
    m.sdMs.assign(sd.size(), NaNv);
    m.floorMask.assign(sd.size(), 0u);
    for (size_t k = 0; k < sd.size() && k < m.absSlope.size(); ++k) {
        if (std::isnan(sd[k]) || !std::isfinite(m.absSlope[k])) continue;
        const double slope = m.absSlope[k];
        if (slope < m.floor) { m.floorMask[k] = 1u; continue; }   // leaves NaN
        m.sdMs[k] = sd[k] / slope * msPerSample;
    }
    return m;
}

// Both panels sit on the SAME third-height grid whether one or two are shown:
//
//   split=false -> panel 1/3, panel(hidden) 0, spacer 2/3
//   split=true  -> panel 1/3, panel        1/3, spacer 1/3
//
// A hidden widget contributes no stretch, so the spacer has to absorb the
// leftover. Without that, a lone visible panel expands to fill half the dock
// and the same landmark is drawn at one scale alone and another once the J
// point has been selected.
void TemplateViewerWindow::setFocusSplit(bool split) {
    if (!m_focusLay || !zoomed_in_section_bottom) return;
    zoomed_in_section_bottom->setVisible(split);
    m_focusLay->setStretch(1, split ? 1 : 0);   // second panel
    m_focusLay->setStretch(2, split ? 1 : 2);   // trailing spacer absorbs the rest
}

// ---- PULSE AND ARTERIAL FOCUS ------------------------------------------
//
// PPG rides its bank slot's own pulse average; the arterial channels have no
// bank and ride the bin's. One panel, not two: pulse channels are foot-anchored
// once, so there is no QRS/JT split and no alignment dimension.
void TemplateViewerWindow::focusPulse(BinPlotWidget* pw, TemplateBin& b,
    int templateIdx, int marker, double col)
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
        // THE GROUP'S PULSE, NOT THE BIN'S: ppg_bank slot i is group i, on
        // the same axis as the bin's pulse template.
        //
        // NO FALLBACK. A group with no pulse cohort has no focus view -- both
        // panels are cleared and the function returns -- because the bin's
        // waveform there would be a measurement attributed to beats that are
        // not in this template.
        const tbank::BankTemplate* ps =
            (templateIdx >= 0 && templateIdx < b.ppg_bank.size())
            ? &b.ppg_bank.templates[templateIdx] : nullptr;
        if (!ps || ps->tmpl.empty() || ps->memberCount() <= 0) {
            clearFocusPanels();
            return;
        }
        meanRaw = &ps->tmpl;           iqrRaw = &ps->tmpl_iqr;
        // THE SLOT'S OWN FOOT, not the bin's. The perfusion transform divides
        // by the value AT this column, and b.ppg_onset was measured on
        // b.ppgTemplate -- a different pulse. showPage normalizes the same
        // trace by ps->pulse_marks.onset, so this is the foot the trace on
        // screen was divided by.
        pulseChan = 0; footIdx = ps->pulse_marks.onset;  chLabel = "PPG";
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
    if (!meanRaw || meanRaw->empty()) { clearFocusPanels(); return; }

    // Pulse channels are NOT normalized by a plain scalar (that was the
    // bug -- it left the trace flat). The displayed trace uses a per-
    // sample PERFUSION-INDEX transform relative to the pulse's own foot,
    // then /ref (normalize_ppg_or_similar -> normalize_pulse_trace, see
    // the main plot ~line 508). The mean MUST use that same transform.
    const std::vector<double> mean = normalize_ppg_or_similar(*meanRaw, footIdx, pulseChan);
    // THE SPREAD GETS THE SPREAD TRANSFORM -- scale_pulse_spread_by_ref,
    // 100/|foot_y|/|ref|, the scale factor of the affine perfusion transform
    // applied to the mean above. It must be the same call showPage makes for
    // the main plot's band, or this panel's band is a different size from the
    // one the same beats draw on the panel it is zooming into.
    const double ref = (pulseChan >= 0 && pulseChan < 4) ? m_pulseGlobalRef[pulseChan] : std::nan("");
    const std::vector<double> sd = normalize_features::scale_pulse_spread_by_ref(
        *iqrRaw, normalize_features::sample_y(*meanRaw, footIdx), ref);

    const int nBeats = (nPulseBeats >= 0)
        ? nPulseBeats : static_cast<int>(b.ppg_n_beats);

    auto pulseLabel = [](int m) -> QString {
        switch (m) {
        case BinPlotWidget::PpgOnset:    return QStringLiteral("Foot");
        case BinPlotWidget::PpgT50:      return QStringLiteral("T50");
            // PpgPeak is the SYSTOLIC peak (the forward-wave maximum the whole
            // pulse is anchored on); PpgPeak2 is the DIASTOLIC peak, the
            // reflected wave after the dicrotic notch.
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

    // ---- THE DETECTOR'S POSITION FOR THIS LANDMARK --------------------
    //
    // From the panel, exactly as focusEcg takes detFid from
    // pw->detectedLandmarks(): detectedPulse() is the detection the X glyphs
    // are drawn at, so the dotted fiducial and the X are one number rather than
    // two positions for one landmark.
    //
    // pw == nullptr on a re-fire path that could not name its panel: leave the
    // fiducial absent rather than substitute a second measurement.
    double detFid = -1.0;
    // The fits behind the three landmarks that HAVE one, for the panel to draw.
    // Invalid for every other pulse marker, which the panel then renders as the
    // detector's position and no curves.
    subsample_refine::PeakCandidates peakCand;
    FocusPanelWidget::FitKind fk = FocusPanelWidget::FitKind::None;
    double peakSigma = 0.0;
    int    peakHalfWidth = 0;
    if (pw) {
        const FeatureMarks::PpgFiducials& pf = pw->detectedPulse();
        // T50 and T80 come from the REACTIVE set, because that is where their
        // X marks come from: they are bar-bracketed crossings, the pulse twin
        // of the ECG P and T peak. Reading pf.t50 here would be the detector's
        // own bracket, which is not the glyph on screen the moment a bar moves.
        const BinPlotWidget::Reactive rx = pw->reactiveGlyphs();
        switch (marker) {
        case BinPlotWidget::PpgOnset:    detFid = pf.onset;    break;
        case BinPlotWidget::PpgPeak:     detFid = pf.peak;     break;
        case BinPlotWidget::PpgDicrotic: detFid = pf.dicrotic; break;
        case BinPlotWidget::PpgPeak2:    detFid = pf.peak2;    break;
        case BinPlotWidget::PpgEnd:      detFid = pf.end;      break;
        case BinPlotWidget::PpgT50:      detFid = rx.ppgT50;   break;
        case BinPlotWidget::PpgT80:      detFid = rx.ppgT80;   break;
        default: break;   // arterial: no per-channel detection to read
        }

        // ---- THE CURVE THAT PLACED IT -----------------------------------
        //
        // The peak is a weighted QUADRATIC's vertex; the foot and end are a
        // CUBIC's. Everything else on a pulse -- the notch (a fixed offset
        // placeholder), the diastolic peak (bracketed search), T50/T80
        // (interpolated crossings) -- is not placed by a polynomial and stays
        // FitKind::None.
        switch (marker) {
        case BinPlotWidget::PpgOnset:
            peakCand = pf.onset_cand;
            fk = FocusPanelWidget::FitKind::PeakCubic;
            peakSigma = subsample_refine::pulse_sigma::Foot;
            peakHalfWidth = subsample_refine::pulse_halfwidth::Foot;
            break;
        case BinPlotWidget::PpgEnd:
            peakCand = pf.end_cand;
            fk = FocusPanelWidget::FitKind::PeakCubic;
            peakSigma = subsample_refine::pulse_sigma::Foot;
            peakHalfWidth = subsample_refine::pulse_halfwidth::Foot;
            break;
        case BinPlotWidget::PpgPeak:
            peakCand = pf.peak_cand;
            fk = FocusPanelWidget::FitKind::PeakQuadratic;
            peakSigma = subsample_refine::pulse_sigma::Peak;
            peakHalfWidth = subsample_refine::pulse_halfwidth::Peak;
            break;
        default: break;
        }

        // ---- INTO THE UNITS THE PANEL DRAWS -----------------------------
        //
        // The fits were made on the RAW pulse average; this panel plots the
        // perfusion-index copy, which is AFFINE in the raw amplitude --
        //   norm(y) = A*y + B,   A = 100 / (|foot_y| * |ref|),
        //                        B = -100 * foot_y / (|foot_y| * |ref|)
        // -- so the constant term takes the offset and the higher powers take
        // the scale only. focusEcg divides all four coefficients by eref
        // because its transform is a bare scale; that would give the right
        // shape at the wrong height here. Positions are COLUMNS: no conversion.
        if (peakCand.valid) {
            const double footY =
                normalize_features::sample_y(*meanRaw, footIdx);
            if (std::isfinite(ref) && ref != 0.0
                && std::isfinite(footY) && std::abs(footY) >= 1e-12) {
                const double A = 100.0 / (std::abs(footY) * std::abs(ref));
                const double B = -100.0 * footY / (std::abs(footY) * std::abs(ref));
                for (auto& f : peakCand.draw) {
                    for (double& cf : f.coeff) cf *= A;
                    f.coeff[0] += B;
                }
            }
            else {
                peakCand = subsample_refine::PeakCandidates{};   // cannot place a curve
            }
        }
    }

    // Pulse channels have no alignment dimension: they are foot-anchored,
    // once, and raw_anchors is ECG-only. One panel, no suffix, and no frame
    // conversion -- `col` is already in this trace's columns.
    setFocusSplit(false);
    if (zoomed_in_section_bottom) zoomed_in_section_bottom->clearFocus();
    if (!zoomed_in_section_top) return;

    const SdMsModel sm = sd_in_msec(mean, sd, m_ppgRateHz);
    zoomed_in_section_top->setFocus(mean, sd, nBeats, static_cast<int>(col),
        chLabel + " " + pulseLabel(marker));
    // Peaks carry a sigma and a window; the rest carry neither, and the
    // one-argument overload is the one that says so -- same split as focusEcg.
    if (fk == FocusPanelWidget::FitKind::PeakQuadratic
        || fk == FocusPanelWidget::FitKind::PeakCubic)
        zoomed_in_section_top->setPeakFitKind(fk, peakSigma, peakHalfWidth);
    else
        zoomed_in_section_top->setFitKind(fk);
    zoomed_in_section_top->setPeakCandidates(peakCand);
    zoomed_in_section_top->setDetectorFiducial(detFid);
    zoomed_in_section_top->setSdMs(sm.sdMs, sm.floorMask, sm.absSlope, sm.floor);
}

// pw may be null on a re-fire path; m_focusWidget keeps the last real one.
void TemplateViewerWindow::refreshFocus(BinPlotWidget* pw, int binIdx, int leadIdx,
    int templateIdx, int marker, double col)
{
    if (!zoomed_in_section_top) return;   // panels not created
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;

    // Remember this focus so a fit-mode radio or an alignment change can replay
    // it in place. THE PANEL TOO, but only when a real one was named: the
    // re-fire paths pass m_focusWidget back in, so overwriting it with their own
    // nullptr would lose the panel on the first replay. m_focusWidget is a
    // QPointer, so a panel destroyed by a page rebuild reads back as null rather
    // than as garbage.
    if (pw) m_focusWidget = pw;
    m_focusBin = binIdx; m_focusLead = leadIdx; m_focusSlot = templateIdx;
    m_focusMarker = marker; m_focusCol = col;

    TemplateBin& b = m_bins[binIdx];
    if (BinPlotWidget::markerIsEcg(marker))
        focusEcg(pw, b, binIdx, leadIdx, templateIdx, marker, col);
    else
        focusPulse(pw, b, templateIdx, marker, col);
}

// ---- ECG FOCUS ---------------------------------------------------------
//
// Rides one alignment's anchored average for one bank slot, with that
// alignment's own per-sample spread and the slot's own beat count. The J point
// is the one landmark bounding two segments, so it gets both panels; every
// other landmark gets the top one.
//
// Split out of refreshFocus for the same reason focusPulse was: the two
// channels share nothing but the panels they write to, and one function doing
// both was 500 lines in which the pulse path could quietly inherit half the ECG
// path's state.
void TemplateViewerWindow::focusEcg(BinPlotWidget* pw, TemplateBin& b,
    int binIdx, int leadIdx, int templateIdx, int marker, double col)
{
    if (leadIdx < 0 || leadIdx > 2) return;

    // THE ALIGNMENT SWITCH, which is the whole feature. The grid draws one
    // alignment on every panel; this panel draws the alignment the clicked
    // landmark is MEASURED on -- P-aligned under the P-onset bar, Q-aligned
    // under Q-onset, and so on -- which is the alignment the bar is stored in
    // and reported under. So the waveform the operator places a bar against is
    // the waveform its column is a column OF.
    m_lastFocusBinIdx = binIdx;
    m_lastFocusLeadIdx = leadIdx;
    m_lastFocusTemplateIdx = templateIdx;
    m_lastFocusMarker = marker;
    m_lastFocusCol = col;

    // A BAR TAKES ITS OWN ALIGNMENT; A GLYPH TAKES THE ONE ON SCREEN.
    //
    // anchorFor returns R_PEAK for every glyph, deliberately -- a glyph is
    // measured on all four averages and has no alignment of its own -- so a
    // glyph must magnify the waveform it was DRAWN on, which is the one
    // currently displayed.
    //
    // Bars use anchorFor rather than currentGridAnchor() because of ORDERING:
    // user_clicked_on_bar calls refreshFocus BEFORE it moves m_autoGridAnchor
    // and re-skins, so currentGridAnchor() here is still the alignment being
    // left. A glyph click moves no anchor, so for glyphs it is current.
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
    // Templates and spreads are stored PRE-reference-division; scale by the
    // same per-lead ref the displayed ECG trace uses. The ECG spread is
    // already an SD, so it feeds the band directly -- no IQR->SD conversion,
    // unlike the pulse channels, whose bank spread is a true interquartile
    // range.
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
    // ---- ONE ACCESSOR, EVERY SLOT ---------------------------------------
    //
    // slotView is what leadsForBinTemplate draws and what ecgDetect measures
    // on, so taking the focus waveform from it means the panel shows the same
    // array the grid shows and the same array the landmarks were found on.
    //
    // NO FALLBACK. A null is build_templates failing to write the per-slot
    // average for this anchor: a writer bug to go and fix, not a state to
    // render under a header naming an alignment it is not.
    const SlotView svF = slotView(b, leadIdx, templateIdx, focusAnchor);
    if (!svF.valid || !svF.tmpl) {
        clearFocusPanels();
        fprintf(stderr, "[focus] bin=%d lead=%d slot=%d anchor=%s"
            " NO PER-SLOT ALIGNED AVERAGE -- build_templates did not write it\n",
            binIdx, leadIdx, templateIdx, anchor_view::label(focusAnchor));
        return;
    }
    const std::vector<double>* meanRawEcg = svF.tmpl;
    const std::vector<double>* sdRawEcg = svF.iqr;

    // The count belongs to the waveform above: the bank slot's clean member
    // count when the slot exists, the bin total only on a file with no bank.
    const uint64_t nb[3] = { b.ch1_n_beats_raw, b.ch2_n_beats_raw, b.ch3_n_beats_raw };
    int nBeats = static_cast<int>(nb[leadIdx]);
    {
        const tbank::TemplateBank& bank = b.ecg_bank[leadIdx];
        if (templateIdx >= 0 && templateIdx < bank.size())
            nBeats = bank.templates[templateIdx].cleanCount();
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

    // `col` arrived in the R frame (the widget's, because the grid is
    // R-aligned); this panel plots THIS alignment's average, so the bar has to
    // be translated into its columns. Same conversion and direction as the drag
    // path in onMarkerMoved, so a bar dragged on the grid lands under the
    // crosshair here -- and the panel FOLLOWS the bar during a drag, because
    // colHere is the live bar position, not a frozen detected column.
    //
    // CLAMPED: the anchors' r_cols genuinely differ, so colHere can fall
    // outside the trace.
    int colHere = col + b.frameShift(leadIdx, AnchorType::R_PEAK, focusAnchor);
    if (mean.empty()) return;
    if (colHere < 0) colHere = 0;
    if (colHere >= static_cast<int>(mean.size()))
        colHere = static_cast<int>(mean.size()) - 1;


    if (zoomed_in_section_top) //if the zoomed in top section is activated (it will always be with any focus)                                                                                                                                                   
    {

        // CACHED ON THE TEMPLATE, NOT RECOMPUTED PER MOUSE-MOVE. A drag
        // re-fires this function on every move with a new `col` and the SAME
        // mean and sd, and sd_in_msec is a Savitzky-Golay pass plus three
        // N-length allocations. (bin, lead, slot, anchor) determines mean and
        // sd completely, so it determines this.
        const SdKey sdKey{ binIdx, leadIdx, templateIdx, focusAnchor };
        if (!(m_sdCacheValid && m_sdCacheKey == sdKey)) {
            m_sdCache = sd_in_msec(mean, sd, m_sampleRate);
            m_sdCacheKey = sdKey;
            m_sdCacheValid = true;
        }
        const SdMsModel& sm = m_sdCache;
        const std::vector<double>& absSlope = sm.absSlope;
        const std::vector<double>& sdMs = sm.sdMs;
        const std::vector<uint8_t>& floorMask = sm.floorMask;
        const double                floor = sm.floor;

        // (The sd at the bar is no longer read here: FocusPanelWidget already
        //  holds m_sdMs and prints it on its own two lines below the plot, so a
        //  copy in the header would be the same number twice.)

        // chForStrict has already returned early if the alignment is absent,
        // so this label can never name an alignment the data is not in.
        const QString tag = QStringLiteral(" [%1-aligned]")
            .arg(QString::fromLatin1(anchor_view::label(focusAnchor)));

        // Header: THE LANDMARK AND ITS ALIGNMENT, nothing else. The sd is on
        // the panel's own footer lines and the model name on its second header
        // line.
        const QString head = QString("%1%2").arg(labelFor(marker), tag);
        // Only the transition bars carry candidates; the peak branch below
        // supplies its own and ignores an invalid set.
        subsample_refine::TransitionCandidates transCand;
        // The detector's own position for the focused landmark, in `mean`'s
        // columns. -1 until the block below supplies it.
        double detFid = -1.0;
        // Hoisted out of the block below, where lm lives: the fit-kind choice
        // needs it and lm does not reach that far. True by default so a path
        // that never detects keeps the old behaviour.
        bool qOnsetFound = true;
        {
            // ---- THE COLUMNS THE GLYPHS ARE DRAWN AT -------------------
            //
            // Taken from the panel, not re-detected. BinPlotWidget runs exactly
            // one detection per (bin, slot, alignment) and draws every X glyph
            // from it; detectedLandmarks() hands that same answer out and
            // reactiveGlyphs() gives the two bar-bracketed peaks. Reading them
            // is what makes the focus mark and the glyph one number rather than
            // two computations to keep in step.
            //
            // pw == nullptr on a path that could not name its panel: leave the
            // fiducial absent rather than substitute a second measurement; the
            // panel then falls back to the bar column.
            if (!mean.empty() && pw) {
                const FeatureMarks::TemplateLandmarks& lm = pw->detectedLandmarks();
                const BinPlotWidget::Reactive rx = pw->reactiveGlyphs();
                qOnsetFound = lm.q_onset_found;

                switch (marker) {
                    // Transitions carry candidate curves as well as a position.
                case BinPlotWidget::EcgPBegin:
                    transCand = lm.p_begin_cand; detFid = lm.p_begin; break;
                case BinPlotWidget::EcgQBegin:
                    transCand = lm.q_onset_cand; detFid = lm.q_onset; break;
                case BinPlotWidget::EcgSEnd:
                    transCand = lm.s_end_cand;   detFid = lm.s_end;   break;
                case BinPlotWidget::EcgTEnd:
                    transCand = lm.t_end_cand;   detFid = lm.t_end;   break;
                    // Detected peaks.
                case BinPlotWidget::EcgRPeak: detFid = lm.r_peak; break;
                case BinPlotWidget::EcgQPeak: detFid = lm.q_peak; break;
                    // Reactive peaks: bracketed by the bars, so they come from
                    // reactiveGlyphs rather than the detection.
                case BinPlotWidget::EcgPPeak: detFid = rx.ecgPPeak; break;
                case BinPlotWidget::EcgTPeak: detFid = rx.ecgTPeak; break;
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

                // ---- INTO THE UNITS THE PANEL DRAWS ---------------------
                //
                // The transition curves come back as CLOSURES OVER THE ARRAY
                // THE DETECTOR FITTED, which is not the array this panel
                // plots, in two independent ways. Both are pure y-axis
                // scalings, so the crossings in cross[] are unaffected and
                // only the drawn curve moves.
                //
                //  1. /eref. ecgDetect measures on the RAW stored slot
                //     average; `mean` above is that array divided by the
                //     per-lead reference. The peak branch below already
                //     divides its coefficients by eref for exactly this
                //     reason -- the transition branch never did, so its
                //     curves were drawn at raw amplitudes on a normalized
                //     axis and sat off the trace by a factor of eref.
                //
                //  2. SIGN. compute_p_begin, compute_q_onset and
                //     compute_j_point all fit `u`, which is -v when the QRS
                //     is negative in this lead, so on those leads the curve
                //     is the mirror of the trace it is drawn over.
                //     compute_t_end fits v directly and needs no flip, but it
                //     rides the same closures, so the sign is resolved per
                //     LANDMARK, not per lead.
                //
                // Negating the closure is exact and is not a re-fit: it is the
                // same fitted model, read in v's units instead of u's.
                if (transCand.valid && std::isfinite(eref) && eref != 0.0) {
                    const bool fitsInverted =
                        (marker == BinPlotWidget::EcgPBegin
                            || marker == BinPlotWidget::EcgQBegin
                            || marker == BinPlotWidget::EcgSEnd)
                        && !FeatureMarks::qrs_positive_at(*meanRawEcg, svF.r_col);
                    const double amp = (fitsInverted ? -1.0 : 1.0) / eref;
                    if (amp != 1.0)
                        for (auto& fn : transCand.curve)
                            if (fn)
                                fn = [inner = fn, amp](double x) {
                                return inner(x) * amp;
                                };
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
            // Fit window for this landmark, from the shared table. Only read
            // on the peak branch below, so the 0 defaults are unreachable
            // placeholders rather than a fallback for a real landmark.
            double peakSigma = 0.0;
            int peakHalfWidth = 0;
            bool isEcgPeakMarker = true;
            EcgPeak whichPeak = EcgPeak::R;
            subsample_refine::PeakCandidates peakCand;   // for the panel to DRAW
            switch (marker) {
            case BinPlotWidget::EcgRPeak: whichPeak = EcgPeak::R; break;
            case BinPlotWidget::EcgPPeak: whichPeak = EcgPeak::P; break;
            case BinPlotWidget::EcgTPeak: whichPeak = EcgPeak::T; break;
            case BinPlotWidget::EcgQPeak: whichPeak = EcgPeak::Q; break;
            default: isEcgPeakMarker = false; break;
            }
            if (isEcgPeakMarker) {
                peakSigma = peakSigmaFor(whichPeak);
                peakHalfWidth = peakHalfWidthFor(whichPeak);

                // THE CURVES ONLY. detFid is already set from the glyph
                // detection above and nothing here may change it, or the focus
                // mark and the glyph part company again. The panel still needs
                // the fitted curves, on `mean`, so their amplitudes are in the
                // units it draws.
                //
                // T-peak excluded: it is a reactive glyph, defined by the
                // S-end / T-end brackets rather than by a fit of its own.
                //
                // pw IS TESTED AGAIN HERE, not inherited from the block above:
                // without it a destroyed panel reaches a `->`, which is the
                // access violation on switching alignment -- the switch tears
                // the panels down and the re-fire hands this function the
                // stale pointer.
                if (pw && whichPeak != EcgPeak::T) {
                    peakCand = pw->peakCandidatesFor(whichPeak);
                    // Coefficients are in the RAW array's amplitude units while
                    // this panel draws the /ref-normalized copy; positions are
                    // columns and need no conversion.
                    if (peakCand.valid && eref > 0.0 && eref != 1.0)
                        for (auto& f : peakCand.draw)
                            for (double& cf : f.coeff) cf /= eref;
                }
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
                // The fits themselves, not the mode: the panel draws them.
                zoomed_in_section_top->setPeakCandidates(peakCand);
                zoomed_in_section_top->setTransitionCandidates(transCand);   // invalid for peaks -> ignored
                zoomed_in_section_top->setDetectorFiducial(detFid);
                zoomed_in_section_top->setSdMs(sdMs, floorMask, absSlope, floor);
            }
        }
    }
}