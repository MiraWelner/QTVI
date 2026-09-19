// ========================================================================
// template_viewer_marks.cpp
// Marker editing: drags, propagation, reseeding and per-bin flags.
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
        for (auto* pw : m_binPlots[li])
            applyTemplateToWidget(pw, m_bins[binIdx], pw->leadIndex(), 0);
    }
}

// The bank-column counterpart of refreshBinMarkers. Repaints only the columns
// showing (binIdx, templateIdx), because a bank slot's bars live in that
// template's own BankMarkerSet and applyBinToWidget would draw the BIN's set
// over them -- the sinus landmarks, on an ectopic waveform.
void TemplateViewerWindow::refreshBankMarkers(int binIdx, int templateIdx) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    if (templateIdx <= 0) { refreshBinMarkers(binIdx); return; }

    // One column, looked up rather than searched for.
    const std::vector<BinPlotWidget*>* col = panelsForColumn(binIdx, templateIdx);
    if (!col) return;
    // Each widget in the column is one lead, and a bank is per lead, so the
    // widget's own leadIndex() selects the bank to draw from -- not the
    // dragged lead, which would paint lead 1's bars onto lead 2's panel.
    for (auto* pw : *col)
        applyTemplateToWidget(pw, m_bins[binIdx],
            pw->leadIndex(), templateIdx);
}

void TemplateViewerWindow::onMarkerMovedOnTemplate(int binIdx, int leadIdx,
    int templateIdx, int marker, int newIdx)
{
    // DISPATCH BY CHANNEL. One signal for every bar; the three channel families
    // want different storage and propagation, so each has its own handler.
    // There is deliberately no templateIdx==0 redirect -- routing slot 0
    // elsewhere is what let the slot-0 and bank paths drift apart historically.
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;

    if (BinPlotWidget::markerIsPpg(marker)) {
        movePpgMarker(binIdx, leadIdx, templateIdx, marker, newIdx);
        return;
    }
    // Arterial (ABP/ART/ART_PULM) still routes to the bin-level handler.
    if (!BinPlotWidget::markerIsEcg(marker)) {
        onMarkerMoved(binIdx, leadIdx, marker, newIdx);
        return;
    }
    moveEcgMarker(binIdx, leadIdx, templateIdx, marker, newIdx);
}

// ---------------------------------------------------------------------------
// PPG drag: SLOT-AWARE, PER-COLUMN. Each morphology column has its OWN pulse
// marks -- slot 0 in the bin-level b.ppg_* fields, slot N in
// ppg_bank.templates[N].pulse_marks -- and that is what each column DRAWS. The
// old path forwarded PPG to onMarkerMoved, which dropped the slot and edited
// only the bin-level fields, so a drag on a bank column (1_b) wrote a field
// only slot 0 (1_a) read: 1_a jumped, 1_c did not. This edits the DRAGGED
// column's own store and propagates to SUBSEQUENT page columns, equal screen
// distance.
// ---------------------------------------------------------------------------
void TemplateViewerWindow::movePpgMarker(int binIdx, int leadIdx, int templateIdx, int marker, int newIdx)
{
    if (templateIdx < 0) return;

    // Per-slot pulse-mark storage for the three draggable bars.
    // SLOT 0 IS A SLOT LIKE ANY OTHER. These four lambdas each had a
    // `slot == 0` branch reading and writing the BIN-level ppg_onset /
    // ppg_dicrotic / ppg_end, while every column DRAWS its own pulse
    // (ppg_bank.templates[slot].tmpl) -- so slot 0's bars were indices into a
    // different waveform from the one under them, and the markings file now
    // serializes pulse_marks per slot, which the bin fields would never reach.
    auto ppgGet = [&](int gi, int slot) -> double {
        TemplateBin& tb = m_bins[gi];
        if (slot >= 0 && slot < (int)tb.ppg_bank.size()) {
            const tbank::BankPulseMarkerSet& pm =
                tb.ppg_bank.templates[slot].pulse_marks;
            switch (marker) {
            case BinPlotWidget::PpgOnset:    return pm.onset;
            case BinPlotWidget::PpgDicrotic: return pm.dicrotic;
            case BinPlotWidget::PpgEnd:      return pm.end;
            }
        }
        return -1.0;
        };
    auto ppgSet = [&](int gi, int slot, double v) {
        TemplateBin& tb = m_bins[gi];
        if (slot >= 0 && slot < (int)tb.ppg_bank.size()) {
            tbank::BankPulseMarkerSet& pm =
                tb.ppg_bank.templates[slot].pulse_marks;
            switch (marker) {
            case BinPlotWidget::PpgOnset:    pm.onset = v; break;
            case BinPlotWidget::PpgDicrotic: pm.dicrotic = v; break;
            case BinPlotWidget::PpgEnd:      pm.end = v; break;
            }
        }
        };
    // Every slot seeds its pulse marks lazily; seed before reading so a
    // never-displayed column still has a real bar to move from. Slot 0 was
    // exempt because its bars came from the bin; it no longer is.
    auto ppgSeed = [&](int gi, int slot) {
        if (slot < 0) return;
        TemplateBin& tb = m_bins[gi];
        if (slot >= (int)tb.ppg_bank.size()) return;
        tbank::BankTemplate& ps = tb.ppg_bank.templates[slot];
        if (ps.tmpl.empty() || ps.hasDetectedPulseMarks()) return;
        FeatureMarks::seed_pulse_bank_template(ps.tmpl, m_ppgRateHz, ps.pulse_marks);
        };
    // Drawn length of this column's pulse (clipped to the ECG window).
    auto ppgLen = [&](int gi, int slot) -> int {
        TemplateBin& tb = m_bins[gi];
        // This slot's own pulse, not the bin's -- it is the trace the column
        // draws and therefore the one its bars are columns of.
        const int rawLen = (slot >= 0 && slot < (int)tb.ppg_bank.size())
            ? (int)tb.ppg_bank.templates[slot].tmpl.size() : 0;
        const int ecgClip = ecgClipLenFor(tb);
        return (ecgClip > 0) ? std::min(rawLen, ecgClip) : rawLen;
        };
    // O(1), from the index showPage builds.
    int dragCol = -1;
    {
        const auto it = m_pageColOf.find(slotKey(binIdx, templateIdx));
        if (it != m_pageColOf.end()) dragCol = it->second;
    }

    // The dragged bar itself.
    ppgSeed(binIdx, templateIdx);
    const double oldIdx = ppgGet(binIdx, templateIdx);
    int placed = newIdx;
    const int dragLen = ppgLen(binIdx, templateIdx);
    if (dragLen > 0) placed = std::clamp(placed, 0, dragLen - 1);
    ppgSet(binIdx, templateIdx, placed);
    // The dragged column's sibling leads share these marks, so they need the
    // position -- but not a re-seed and re-detect.
    if (dragCol >= 0 && dragCol < (int)m_binPlots.size())
        for (auto* pw : m_binPlots[dragCol])
            if (pw)
                pw->setMarker(static_cast<BinPlotWidget::Marker>(marker),
                    static_cast<double>(placed));

    if (m_moveMode == MoveMode::Individual || oldIdx < 0) {
        refreshFocus(binIdx, leadIdx, templateIdx, marker, placed);
        return;
    }

    // Propagate to every LATER column on the page, equal screen distance.
    const double delta = placed - oldIdx;
    const double dragSpan = binSpanSeconds(binIdx);
    for (int li = dragCol + 1; li < (int)m_pageGlobalIdx.size()
        && li < (int)m_pageTemplateIdx.size(); ++li) {
        const int gi = m_pageGlobalIdx[li];
        const int slot = m_pageTemplateIdx[li];
        if (gi < 0 || slot < 0) continue;
        if (m_bins[gi].bad_ppg != 0) continue;
        ppgSeed(gi, slot);
        const double cur = ppgGet(gi, slot);
        if (cur < 0.0) continue;
        const int n = ppgLen(gi, slot);
        if (n <= 0) continue;
        const double tgtSpan = binSpanSeconds(gi);
        if (!(dragSpan > 0.0) || !(tgtSpan > 0.0)) continue;
        const double target = cur + delta * (tgtSpan / dragSpan);
        if (target < 0.0 || target > n - 1) continue;
        ppgSet(gi, slot, target);

        // PUSH ONE BAR -- as the ECG path, and likewise the whole update.
        // Pulse marks are per (bin, slot), shared across a column's leads, so
        // every panel in the column gets the push.
        if (li >= 0 && li < (int)m_binPlots.size())
            for (auto* pw : m_binPlots[li])
                if (pw)
                    pw->setMarker(static_cast<BinPlotWidget::Marker>(marker),
                        target);
    }
    refreshFocus(binIdx, leadIdx, templateIdx, marker, placed);
}

// ---------------------------------------------------------------------------
// ONE ECG DRAG HANDLER, FOR EVERY SLOT. slotMarks() makes the storage uniform
// (slot 0 is a bank slot like any other), so there is one path, not two that
// drift. Reads/writes are in the OWNER alignment's frame; the widget measures
// in the R frame, so every write converts (view -> owner) and every read
// converts back. Bounds come from the DRAWN wall (lastDrawnSample), not the
// array end (which runs past it after the tail trim). Out-of-range targets are
// CLAMPED, not skipped, so one panel can't abandon the rest of the page. The
// per-drag origin map is keyed binI*64+slot (never the page column, which
// collides). m_touchedMarks is recorded for every slot.
// ---------------------------------------------------------------------------
void TemplateViewerWindow::moveEcgMarker(int binIdx, int leadIdx,
    int templateIdx, int marker, int newIdx)
{
    if (leadIdx < 0 || leadIdx > 2) return;
    if (templateIdx < 0) return;

    TemplateBin& b = m_bins[binIdx];

    // R PEAK IS PER BIN, NOT PER SLOT: every template in the bank is aligned
    // on it. Handled and returned here, because none of the bar machinery
    // below applies -- no BankMarkerSet field, no frame conversion (R IS the
    // frame), no propagation. markerAtX does not hand R out, so this is
    // defensive rather than reachable.
    if (marker == BinPlotWidget::EcgRPeak) {
        if (templateIdx == 0) b.r_peak_ch[leadIdx] = newIdx;
        return;
    }
    if (!anchor_view::isBar(marker)) return;   // glyphs are not draggable

    // WHICH CELL THIS DRAG EDITS: the bar's OWNER, always, forced or not.
    // A bar has one home -- p_begin in P_ONSET, q_onset / s_end / t_end in
    // Q_ONSET -- and every view reads that home through barsForPanel, so a
    // drag has to write it or the edit lands somewhere nothing reads.
    // WHICH CELL THIS DRAG EDITS, and it matches barsForPanel exactly.
    //
    // Forced P or Q: m_forcedAlign IS the bar's owner, so this writes the
    // canonical cell -- the same one Automatic reads. Forced R or J: that
    // alignment's own cell, which only it displays. Automatic: the owner.
    // THE SAME CELL barsForPanel READS. Forced R or J edits that alignment's
    // own bar; everything else edits the canonical one (anchorFor). See
    // anchor_view::ownsCanonicalBar for why those two are different things.
    const AnchorType owner =
        (m_forceAlign && anchor_view::hasOwnBars(m_forcedAlign))
        ? m_forcedAlign
        : anchor_view::anchorFor(marker);
    if (!anchor_view::showsBar(owner, marker)) return;

    // ---- storage: one accessor pair, any (bin, slot) ----------------------
    auto get = [&](TemplateBin& tb, int slot) -> double {
        const tbank::BankMarkerSet& m = tb.slotMarks(leadIdx, slot, owner);
        switch (marker) {
        case BinPlotWidget::EcgPBegin: return m.p_begin;
        case BinPlotWidget::EcgQBegin: return m.q_onset;
        case BinPlotWidget::EcgSEnd:   return m.s_end;
        case BinPlotWidget::EcgTEnd:   return m.t_end;
        }
        return -1.0;
        };
    auto set = [&](TemplateBin& tb, int slot, double v) {
        tbank::BankMarkerSet& m = tb.slotMarks(leadIdx, slot, owner);
        switch (marker) {
        case BinPlotWidget::EcgPBegin: m.p_begin = v; break;
        case BinPlotWidget::EcgQBegin: m.q_onset = v; break;
        case BinPlotWidget::EcgSEnd:   m.s_end = v; break;
        case BinPlotWidget::EcgTEnd:   m.t_end = v; break;
        }
        };

    // ---- the panel's waveform and its DRAWN wall --------------------------
    auto wallAt = [&](int li) -> int {
        if (li < 0 || li >= (int)m_binPlots.size()) return -1;
        for (auto* pw : m_binPlots[li])
            if (pw && pw->leadIndex() == leadIdx)
                return pw->lastDrawnSample(BinPlotWidget::Channel::Ecg);
        return -1;
        };
    // The VISIBLE x-axis span, in seconds (equal-screen-distance propagation).
    auto spanAt = [&](int li) -> double {
        if (li < 0 || li >= (int)m_binPlots.size()) return -1.0;
        for (auto* pw : m_binPlots[li])
            if (pw && pw->leadIndex() == leadIdx)
                return pw->frameTMax() - pw->frameTMin();   // seconds
        return -1.0;
        };

    // O(1), from the index showPage builds.
    int dragCol = -1;
    {
        const auto it = m_pageColOf.find(slotKey(binIdx, templateIdx));
        if (it != m_pageColOf.end()) dragCol = it->second;
    }

    // SEED BEFORE READING, for every slot (marks() is operator[]; a read
    // inserts, which used to suppress a slot's auto-detection permanently).
    auto seedIfNeeded = [&](int binI, int slot) {
        // No slot-0 exclusion: one seeding call for every slot.
        tbank::TemplateBank& bk = m_bins[binI].ecg_bank[leadIdx];
        if (slot >= (int)bk.templates.size()) return;
        tbank::BankTemplate& tgt = bk.templates[slot];
        if (tgt.hasDetectedMarks(static_cast<int>(owner))) return;
        // This alignment's average, not the R-aligned one: tgt.tmpl/tgt.r_col
        // are the R-frame pair, and filing that under `owner` is what put an
        // R-frame number in the P-onset bar.
        const SlotView sv = slotView(m_bins[binI], leadIdx, slot, owner);
        if (!sv.valid) return;
        FeatureMarks::seed_bank_template(*sv.tmpl, sv.r_col, m_sampleRate,
            owner, tgt.marks(static_cast<int>(owner)));
        };

    // ---- ONE FRAME FOR ALL THE ARITHMETIC --------------------------------
    //
    // THE INVARIANT: bars are STORED in the owner alignment's columns and
    // MEASURED in the frame the grid is drawing, so the conversion belongs on
    // BOTH sides of the accessor. It used to be on the write only, which left
    // every propagated bar reading owner-framed, shifted a second time on the
    // way back, and clamped against a drawn-frame wall. frameShift is 0 when
    // the grid draws the bar's own alignment, so it surfaced as forced-P plus
    // Move-Subsequent.
    //
    // Below this line, every column is a DRAWN-frame column.
    auto getView = [&](TemplateBin& tb, int slot) -> double {
        const double v = get(tb, slot);
        return (v < 0.0) ? -1.0
            : v - tb.frameShift(leadIdx, currentGridAnchor(), owner);
        };
    auto setView = [&](TemplateBin& tb, int slot, double v) {
        set(tb, slot, v + tb.frameShift(leadIdx, currentGridAnchor(), owner));
        };

    // ---- the dragged bar --------------------------------------------------
    seedIfNeeded(binIdx, templateIdx);
    const double oldIdx = getView(b, templateIdx);
    if (m_dragStartIdx < 0) m_dragStartIdx = oldIdx;   // first move of this drag

    const int dragWall = wallAt(dragCol);

    int placed = newIdx;
    if (dragWall > 0) placed = std::clamp(placed, 0, dragWall);
    // Stored through setView (view -> owner), so userMarks reads it back at
    // `placed` in whatever frame the grid is drawing.
    setView(b, templateIdx, placed);

    if (placed >= 0)
        m_touchedMarks[touchKey(binIdx, leadIdx, marker, owner)] = placed;

    if (m_moveMode == MoveMode::Individual || oldIdx < 0) {
        refreshFocus(binIdx, leadIdx, templateIdx, marker, placed);
        return;
    }

    // ---- propagation, over the PAGE'S OWN COLUMNS (equal screen distance) --
    const double dragSpan = spanAt(dragCol);

    for (int li = dragCol + 1; li < (int)m_pageGlobalIdx.size()
        && li < (int)m_pageTemplateIdx.size(); ++li) {
        const int gi = m_pageGlobalIdx[li];
        const int slot = m_pageTemplateIdx[li];
        if (gi < 0 || gi >= (int)m_bins.size()) continue;
        if (slot < 0) continue;
        if (m_bins[gi].bad_r_ch[leadIdx]) continue;

        const int wall = wallAt(li);
        if (wall <= 0) continue;

        seedIfNeeded(gi, slot);
        const double cur = getView(m_bins[gi], slot);
        if (cur < 0) continue;   // this landmark was not found on this one

        const int key = slotKey(gi, slot);
        if (!original_location_of_bar.count(key))
            original_location_of_bar[key] = (int)std::lround(cur);

        const double tgtSpan = spanAt(li);
        if (!(dragSpan > 0.0) || !(tgtSpan > 0.0)) continue;
        const int sampleShift = (int)std::lround(
            double(placed - m_dragStartIdx) * (tgtSpan / dragSpan));
        int target = originFor(key, (int)std::lround(cur)) + sampleShift;

        // CLAMP, DON'T SKIP, AND KEEP A BUFFER (edges are ungrabbable).
        const int kEdgeBuffer = 5;
        if (wall <= 2 * kEdgeBuffer) continue;
        target = std::clamp(target, kEdgeBuffer, wall - kEdgeBuffer);

        setView(m_bins[gi], slot, target);

        // PUSH ONE BAR, NOT THE WHOLE COLUMN, AND NOTHING ON RELEASE. A bar
        // move changes nothing that refreshBankMarkers recomputes (seeding,
        // pulse marks, glyph detection), and the glyphs that DO follow a bar
        // are reactive at paint time -- so this push is the whole update, not
        // a cheap stand-in for one. `target` is already a drawn-frame column,
        // which is what setMarker wants.
        for (auto* pw : m_binPlots[li])
            if (pw && pw->leadIndex() == leadIdx)
                pw->setMarker(static_cast<BinPlotWidget::Marker>(marker),
                    static_cast<double>(target));
    }

    // The dragged panel's own focus view (J-point refreshes both QRS and JT).
    refreshFocus(binIdx, leadIdx, templateIdx, marker, placed);
}

void TemplateViewerWindow::onMarkerMoved(int binIdx, int leadIdx,
    int marker, int newIdx)
{
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    TemplateBin& b = m_bins[binIdx];

    // ECG AND PULSE FORWARD; ONLY THE ARTERIAL BODY LIVES HERE.
    //
    // This function is reachable from exactly one place:
    // onMarkerMovedOnTemplate's arterial dispatch. BinPlotWidget::markerMoved,
    // the signal it was written for, is declared but never emitted anywhere in
    // the tree -- markerMovedOnTemplate replaced it -- so the ~150 lines of ECG
    // and PPG handling that used to sit here were unreachable, and were a
    // second, slot-less implementation of what moveEcgMarker and movePpgMarker
    // already do. The PPG copy in particular wrote the bin-level ppg_* fields
    // with no slot dimension, which is the defect movePpgMarker exists to fix,
    // so keeping it around was keeping the old bug one connect() away.
    //
    // Both kinds now forward to the slot-aware handlers with slot 0 (a drag
    // arriving without a slot is by definition on the bin's first column), so
    // if the legacy signal is ever wired up it behaves correctly instead of
    // diverging. Arterial must NOT forward: onMarkerMovedOnTemplate routes it
    // back here, which would recurse.
    if (BinPlotWidget::markerIsEcg(marker) || BinPlotWidget::markerIsPpg(marker)) {
        onMarkerMovedOnTemplate(binIdx, leadIdx, /*templateIdx=*/0, marker, newIdx);
        return;
    }

    // ---- Arterial markers (ABP / ART / ART_PULM) --------------------------
    // Shared across leads like PPG. Route to the right channel's fields and
    // propagate to subsequent bins when Move-Subsequent is on.
    if (BinPlotWidget::markerIsArterial(marker)) {
        // Field access is the shared table; only the trace and issue flag
        // still need a per-channel branch.
        auto channelTrace = [&](TemplateBin& tb, int mk,
            const std::vector<double>*& tr, uint8_t*& iss) {
                if (BinPlotWidget::markerIsAbp(mk)) { tr = &tb.abpTemplate; iss = &tb.abp_issue; }
                else if (BinPlotWidget::markerIsArt(mk)) { tr = &tb.artTemplate; iss = &tb.art_issue; }
                else { tr = &tb.artPulmTemplate; iss = &tb.art_pulm_issue; }
            };

        double* dragged = pulseField(b, marker);
        if (!dragged) return;
        const int oldIdx = (int)*dragged;
        *dragged = newIdx;
        refreshBinMarkers(binIdx);
        const int delta = newIdx - oldIdx;

        if (m_moveMode != MoveMode::Individual && oldIdx >= 0) {
            const std::vector<double>* trDragged = nullptr; uint8_t* issDragged = nullptr;
            channelTrace(b, marker, trDragged, issDragged);
            // PERCENT OF THE X-AXIS MOVED, same rule as the ECG/PPG paths.
            const int rawLenD = trDragged ? (int)trDragged->size() : 0;
            const int ecgClipD = ecgClipLenFor(b);
            const int draggedN = (ecgClipD > 0 && rawLenD > 0)
                ? std::min(rawLenD, ecgClipD) : rawLenD;
            const double axisFrac = (draggedN > 1) ? double(delta) / (draggedN - 1) : 0.0;

            for (int i = binIdx + 1; i < (int)m_bins.size(); ++i) {
                const std::vector<double>* tr = nullptr; uint8_t* iss = nullptr;
                channelTrace(m_bins[i], marker, tr, iss);
                if (!iss || *iss != 0) continue;
                if (!tr) continue;
                double* fld = pulseField(m_bins[i], marker);
                if (!fld) continue;
                const int cur = (int)*fld;
                if (cur < 0) continue;

                const int rawLen = (int)tr->size();
                const int ecgClip = ecgClipLenFor(m_bins[i]);
                const int n = (ecgClip > 0) ? std::min(rawLen, ecgClip) : rawLen;
                if (n <= 0) continue;

                const int target = cur + (int)std::lround(axisFrac * (n - 1));
                if (target < 0 || target > n - 1) continue;
                *fld = target;
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

void TemplateViewerWindow::user_clicked_on_bar(int binIdx, int leadIdx, int templateIdx, int marker, double col)
{
    //the user clicked a bar. Record that it was confirmed, align automatically accordingly, and load/refresh the focus panel on the side
    // Against the alignment on screen: that is the cell whose bar was clicked.
    if (binIdx >= 0 && leadIdx >= 0 && col >= 0)
        m_touchedMarks[touchKey(binIdx, leadIdx, marker, currentGridAnchor())] = col;
    refreshFocus(binIdx, leadIdx, templateIdx, marker, col);
    if (!m_forceAlign && BinPlotWidget::markerIsEcg(marker) && anchor_view::isBar(marker))
    {
        const AnchorType a = anchor_view::anchorFor(marker);
        if (a != m_autoGridAnchor) {
            m_autoGridAnchor = a;
            reskinGridForAnchor();
        }
    }
}
// Clear every operator-placed bar on the CURRENT PAGE, so each one falls back
// to the detection again.
//
// THAT IS THE WHOLE OF A RESET NOW. A cell holds operator edits and nothing
// else, and barsForPanel returns the panel's detection for any bar whose cell
// is -1 -- so emptying the cells restores the detected positions without
// re-running anything. (This function used to re-detect and write the results
// back in; it was briefly an empty stub, which reset nothing at all.)
void TemplateViewerWindow::resetMarks() {
    for (size_t i = 0; i < m_pageGlobalIdx.size(); ++i) {
        const int gi = m_pageGlobalIdx[i];
        const int slot = (i < m_pageTemplateIdx.size())
            ? m_pageTemplateIdx[i] : 0;
        if (gi < 0 || gi >= static_cast<int>(m_bins.size())) continue;
        TemplateBin& b = m_bins[gi];

        for (int lead = 0; lead < 3; ++lead) {
            for (AnchorType a : anchor_view::kAllAnchors) {
                tbank::BankMarkerSet& m = b.slotMarks(lead, slot, a);
                m.p_begin = -1.0;
                m.q_onset = -1.0;
                m.s_end = -1.0;
                m.t_end = -1.0;
                // The bar is no longer operator-confirmed either, so the
                // boundary log must not report it as ground truth.
                m_touchedMarks.erase(
                    touchKey(gi, lead, anchor_view::kPBegin, a));
                m_touchedMarks.erase(
                    touchKey(gi, lead, anchor_view::kQBegin, a));
                m_touchedMarks.erase(
                    touchKey(gi, lead, anchor_view::kSEnd, a));
                m_touchedMarks.erase(
                    touchKey(gi, lead, anchor_view::kTEnd, a));
            }
        }
    }

    // Re-push through the one path that reads the cells.
    for (size_t i = 0; i < m_binPlots.size(); ++i) {
        const int gi = (i < m_pageGlobalIdx.size()) ? m_pageGlobalIdx[i] : -1;
        const int slot = (i < m_pageTemplateIdx.size())
            ? m_pageTemplateIdx[i] : 0;
        if (gi < 0 || gi >= static_cast<int>(m_bins.size())) continue;
        for (BinPlotWidget* pw : m_binPlots[i]) {
            if (!pw) continue;
            const int lead = pw->leadIndex();
            if (lead < 0 || lead > 2) continue;
            applyTemplateToWidget(pw, m_bins[gi], lead, slot);
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

void TemplateViewerWindow::onBadRToggled(int binIdx, int leadIdx,
    int templateIdx, bool bad) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;
    if (leadIdx < 0 || leadIdx > 2) return;

    if (tbank::BankTemplate* t = slotFor(binIdx, leadIdx, templateIdx))
        // BOOL. The 1u/2u looked like a tag for which verdict this is, but the
        // field is a bool -- which slot holds it IS the distinction: this one is
        // the ECG lead's slot.
        t->marked_invalid_template = bad;

    // SLOT 0 ONLY writes the bin-level flag. See the header note above.
    if (templateIdx == 0) m_bins[binIdx].bad_r_ch[leadIdx] = bad;

    repaintPanel(binIdx, leadIdx, templateIdx,
        panelState(binIdx, leadIdx, templateIdx));
}

void TemplateViewerWindow::onBadPPGToggled(int binIdx, int templateIdx,
    bool bad) {
    if (binIdx < 0 || binIdx >= (int)m_bins.size()) return;

    // The pulse verdict is recorded on the PULSE bank's slot, not on the ECG
    // lead's -- it is a statement about the pulse waveform in this panel.
    if (tbank::BankTemplate* t = slotFor(binIdx, -1, templateIdx))
        // BOOL -- see onBadRToggled. This is the PULSE slot, which is what
        // makes it the pulse verdict.
        t->marked_invalid_template = bad;

    if (templateIdx == 0) {
        m_bins[binIdx].bad_ppg = bad ? 1 : 0;
        // NO LONGER CLEARS bad_r. The two were alternatives in the old
        // three-step right-click cycle; the cycle now has a both-bad step, so
        // marking the pulse must leave the ECG verdict alone.
    }

    // Every lead of THIS panel: a bad pulse is not a per-lead judgement, and the
    // panel shows the same pulse trace under each lead.
    // Per lead: the pulse verdict is shared across a panel's leads, the ECG one
    // is not, so each lead's combined state can differ.
    for (int lead = 0; lead < 3; ++lead)
        repaintPanel(binIdx, lead, templateIdx,
            panelState(binIdx, lead, templateIdx));
}