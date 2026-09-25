// ========================================================================
// template_viewer_setup.cpp
// Construction, UI wiring, subject load and alignment-radio plumbing.
//
// One of five units; TemplateViewerWindow is declared in template_viewer.hpp.
// ========================================================================

#include "template_viewer.hpp"



// ========================================================================
// Construction
// ========================================================================

TemplateViewerWindow::TemplateViewerWindow(QWidget* parent)
    : QMainWindow(parent), ui(new Ui::TemplateViewerWindow)
{
    ui->setupUi(this);

    //when you are moving a marker, do you move just that marker, or all subsequent markers by the same screen distance?
    auto* moveGroup = new QButtonGroup(this);
    moveGroup->addButton(ui->move_individual);
    moveGroup->addButton(ui->move_subsequent_delta);
    moveGroup->setExclusive(true);

    connect(ui->move_individual, &QRadioButton::toggled, this, [this](bool on) {
        if (on) m_moveMode = MoveMode::Individual;
        });
    connect(ui->move_subsequent_delta, &QRadioButton::toggled, this, [this](bool on) {
        if (on) m_moveMode = MoveMode::SubsequentDelta;
        });

    // Sync m_moveMode to whichever button Designer has checked by default.
    if (ui->move_individual->isChecked())      m_moveMode = MoveMode::Individual;
    else if (ui->move_subsequent_delta->isChecked()) m_moveMode = MoveMode::SubsequentDelta;

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
    // Fit-model radios. Picking a model forces it (Auto = the BIC contest);
    // the focus fit follows immediately by re-running the current focus, and
    // the invalidated cache makes the next landmark use it too.
    auto applyFitMode = [this]() {
        // Switching the model: re-detect the X-mark glyphs at the selected
        // model's crossing. setFitModes stores the mode + invalidates the glyph
        // snapshot; setAuto re-captures it (that is the ONLY path that re-runs
        // detect_template_landmarks). Both are needed -- without setAuto only
        // the reactive P-onset glyph moved, since the detected Q/S/T glyphs are
        // frozen in the snapshot. In place, per panel: no page rebuild, no
        // shift, dragged BARS untouched. Then the focus panel recolors.
        // (m_lastTransKey = -1 was here: the focus path no longer caches a
        //  detection of its own, so there is nothing to invalidate --
        //  setFitModes below clears the panel caches that do exist.)
        const AnchorType frame = currentGridAnchor();
        for (int li = 0; li < (int)m_binPlots.size()
            && li < (int)m_pageGlobalIdx.size(); ++li) {
            const int gi = m_pageGlobalIdx[li];
            if (gi < 0 || gi >= (int)m_bins.size()) continue;
            for (auto* pw : m_binPlots[li]) {
                if (!pw) continue;
                pw->setFitModes(m_onOffsetFitMode, m_peakFitMode);
                pw->setAuto(m_bins[gi], frame);   // re-capture glyphs with the mode
            }
        }
        if (m_focusMarker >= 0)
            refreshFocus(m_focusWidget, m_focusBin, m_focusLead, m_focusSlot, m_focusMarker, m_focusCol);
        };
    auto wireOnOffset = [this, applyFitMode](const char* name, curve_fit::FitMode mode) {
        if (auto* rb = findChild<QRadioButton*>(name))
            connect(rb, &QRadioButton::toggled, this,
                [this, applyFitMode, mode](bool on) {
                    if (!on) return; m_onOffsetFitMode = mode; applyFitMode();
                });
        else   // SAY SO. findChild returning null here is a radio that silently
               // does nothing when clicked, with no diagnostic anywhere.
            fprintf(stderr, "[fit-onoffset] NOT WIRED: %s\n", name);
        };
    wireOnOffset("auto_fit_onoffset", curve_fit::FitMode::Auto);
    wireOnOffset("linear_fit_onoffset", curve_fit::FitMode::Linear);
    wireOnOffset("cubic_spline_fit_onoffset", curve_fit::FitMode::CubicSpline);
    wireOnOffset("cubic_fit_onoffset", curve_fit::FitMode::Cubic);
    wireOnOffset("sigmoid_fit_onoffset", curve_fit::FitMode::Sigmoid);
    wireOnOffset("fracpoly_fit_onoffset", curve_fit::FitMode::FracPoly);

    // The initial selection is set in the .ui (fit_peaks_auto is checked), so it
    // happens inside setupUi -- before these connects -- and therefore does not
    // fire applyFitMode during construction, when no page or panel exists yet.
    auto wirePeak = [this, applyFitMode](const char* name, curve_fit::PeakFitMode mode) {
        if (auto* rb = findChild<QRadioButton*>(name))
            connect(rb, &QRadioButton::toggled, this,
                [this, applyFitMode, mode](bool on) {
                    if (!on) return; m_peakFitMode = mode; applyFitMode();
                });
        else
            fprintf(stderr, "[fit-peaks] NOT WIRED: %s\n", name);
        };
    wirePeak("fit_peaks_auto", curve_fit::PeakFitMode::Auto);
    wirePeak("fit_peaks_cubic", curve_fit::PeakFitMode::Cubic);
    wirePeak("fit_peaks_parabola", curve_fit::PeakFitMode::Parabola);
    wirePeak("fit_peaks_5pt", curve_fit::PeakFitMode::FivePoint);

    if (auto* resetBtn = findChild<QPushButton*>("reset_marks"))
    {
        connect(resetBtn, &QPushButton::clicked, this, &TemplateViewerWindow::resetMarks);
    }
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
    wirePpgAlignButtons();
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
// Load subject
// ========================================================================

void TemplateViewerWindow::loadSubject(const QString& templatePath, const QString& markingPath,
    const QString& subjectId, double sampleRateHz,
    double ppgRateHz, double abpRateHz, double artRateHz, double artPulmRateHz,
    double notchFilterHz) {

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

    // ONE CALL, TWO ENTRY POINTS. Everything below used to be the back half
    // of loadSubject, and it moved out when the in-memory overload arrived --
    // rather than being copied into it, because the four-pass seeding loop is
    // the ONE place all four alignments are detected and two copies of it
    // would drift.
    initAfterBinsLoaded();
}

// ---------------------------------------------------------------------------
// THE IN-MEMORY OVERLOAD
// ---------------------------------------------------------------------------
//
// Takes the TemplateFile post_process already holds, rather than a path. The
// disk round trip it replaces produced nothing but a filename -- and it was
// the reason a half-populated templates.bin had to exist, because
// prepareViewerJob wrote one before the squared/absval blocks were built. That
// file looked complete and was not, which is what the _templates.partial.bin
// and its remove+rename promote were papering over. templates.bin is now
// written once, at the end of the anchor cycle, so its existence means
// complete.
void TemplateViewerWindow::loadSubject(const template_io::TemplateFile& tf,
    const QString& templateDir, const QString& markingPath,
    const QString& subjectId, double sampleRateHz,
    double ppgRateHz, double abpRateHz, double artRateHz, double artPulmRateHz,
    double notchFilterHz) {

    m_markingPath = markingPath;
    // EXPLICIT, not derived. The path overload takes it from the filename;
    // there is no filename here, and captureCurrentPage and the bins CSV both
    // write into it.
    m_templateDir = templateDir;
    m_subjectId = subjectId;
    m_sampleRate = sampleRateHz;
    m_ppgRateHz = ppgRateHz;
    m_abpRateHz = abpRateHz;
    m_artRateHz = artRateHz;
    m_artPulmRateHz = artPulmRateHz;
    m_notchFilterHz = notchFilterHz;
    setTitleForSubject();
    ui->subjectLabel->setText(subjectId);
    ui->finishButton->setText("Finish");

    // NO try/catch. The path overload guards a FILE READ; there is no read
    // here, so the only failure left is an empty set, checked below.
    m_bins = binsFromTemplateFile(tf);

    if (m_bins.empty()) {
        QMessageBox::warning(this, "Error",
            "No bins in the template data for " + subjectId);
        emit finished();
        return;
    }

    initAfterBinsLoaded();
}

void TemplateViewerWindow::initAfterBinsLoaded() {
    //various bookeeping after the bins are loaded

    // PER-SUBJECT STATE, DROPPED HERE AND NOWHERE ELSE. All three are keyed on
    // (bin, slot) or on bin alone, so carrying them across a subject change
    // serves bin 3 of the next record with bin 3 of this one: a beat matrix
    // from another patient, and a "built" waveform that would be restored over
    // a template it never described.
    clearBeatsCache();
    m_ppgBuilt.clear();
    m_ppgRealigned.clear();

    for (TemplateBin& b : m_bins) b.polarity = m_polarity;

    max_leads = 1;
    for (const auto& b : m_bins) {
        int nl = (int)leadsForBin(b).size();
        if (nl > max_leads) max_leads = nl;
    }
    m_currentPage = 0;
    buildPages();
    QGuiApplication::setOverrideCursor(Qt::WaitCursor);
    QtConcurrent::blockingMap(m_bins,
        [this](TemplateBin& b) { seedOneBin(b); });
    QGuiApplication::restoreOverrideCursor();
    const QDir markingDir(m_markingPath);
    const QString canonical = markingDir.filePath(m_subjectId + "_template_markings.bin");

    bool markersReloaded = false;
    if (QFile::exists(canonical)) {
        markersReloaded = restoreMarkersFrom(canonical, /*ecg=*/true, /*pulse=*/true);
        if (markersReloaded) {
            for (auto& b : m_bins) b.syncReactivePpg();
        }
    }

    for (TemplateBin& b : m_bins) {
        if (b.bad_ppg != 1) continue;
        for (tbank::BankTemplate& t : b.ppg_bank.templates)
            t.marked_invalid_template = true;
    }
    compute_global_refs();
    showPage();
}

// Everything that has to happen to ONE bin at load, in the order it has to
// happen in. Called concurrently, one bin per call, from initAfterBinsLoaded.
// Reads this window's rates and fit modes and writes nothing but `b`, which is
// what makes the concurrency safe -- see the note at the call site.
void TemplateViewerWindow::seedOneBin(TemplateBin& b) const
{
    const std::array<ChannelTemplateData, 3> savedR = { b.ch1, b.ch2, b.ch3 };

    for (AnchorType a : anchor_view::anchor_array) {
        // R IS SEEDED ONCE, AFTER THIS LOOP, NOT IN IT. kAllAnchors puts
        // R_PEAK first (anchor_view.hpp), so the old code seeded R here and
        // then seeded it AGAIN at the end -- five passes per bin for four
        // alignments' worth of result. The trailing pass is the one that has
        // to stay, because it is what leaves the flat ch1..3 and *_auto_ch
        // fields holding R, which is what the grid draws; this one was pure
        // duplicate work.
        if (a == AnchorType::R_PEAK) continue;

        auto it = b.anchored.find(static_cast<int>(a));
        if (it == b.anchored.end()) continue;   // no block -> nothing to seed
        b.ch1 = it->second[0]; b.ch2 = it->second[1]; b.ch3 = it->second[2];
        FeatureMarks::seed_all(b, m_sampleRate, m_ppgRateHz, a, b.polarity);
    }

    // R last so the flat state the grid reads is R's.
    b.ch1 = savedR[0]; b.ch2 = savedR[1]; b.ch3 = savedR[2];
    FeatureMarks::seed_all(b, m_sampleRate, m_ppgRateHz, AnchorType::R_PEAK, b.polarity);

    // NO ECG GLYPH SYNC: p_peak is not stored any more, so there is nothing to
    // cache. The PPG reactive values ARE cached (t50 / t80 / t80_rise / pw80 /
    // peak2), so they are rederived here from the bars the seed just wrote
    // plus the auto-detected systolic peak.
    b.syncReactivePpg();

    // NO BAR SEEDING ANY MORE. barsForPanel falls back to the panel's own
    // detection for any bar the operator has not placed, so a cell holds
    // operator edits and nothing else and there is nothing to pre-fill. The
    // pass that was here ran a SECOND detection over every (lead, slot,
    // alignment) -- the one that could, and did, disagree with the detection
    // the glyphs are drawn from.
    //
    // seed_all above is untouched: it fills the flat *_auto_ch columns, which
    // are what the CSVs report and are not bars.
}

bool TemplateViewerWindow::restoreMarkersFrom(const QString& markingsBinPath, bool ecg, bool pulse) {
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
        // DOUBLE, and the bound compared as one: casting a fractional saved
        // position to size_t truncates, which would let len-0.5 through.
        auto safeIdx = [&](double savedVal, double currentVal, size_t len) -> double {
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
                    const int nSaved = s.slotCount(c);
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
                            // BARS ONLY. p_peak is not in the record and not on
                            // the struct: readers call reactive_ecg on these.
                            dm.p_begin = safeIdx(sm.p_begin, dm.p_begin, len);
                            dm.q_onset = safeIdx(sm.q_onset, dm.q_onset, len);
                            dm.s_end = safeIdx(sm.s_end, dm.s_end, len);
                            dm.t_end = safeIdx(sm.t_end, dm.t_end, len);
                        }
                    }
                }
                for (int c = 0; c < 3; ++c) d.bad_r_ch[c] = s.bad_r_ch[c];
            } // if (ecg)

            if (pulse) {
                d.bad_ppg = s.bad_ppg;
                // BARS ONLY. t50 / peak / peak2 / t80 are auto-only glyphs and
                // are not in the record; the caller calls syncReactivePpg()
                // once the merge is done. (ppg_t80_rise / ppg_pw80 were never
                // merged here even when the file carried them.)
                d.ppg_onset = safeIdx(s.ppg_onset, d.ppg_onset, ppgLen);
                d.ppg_dicrotic = safeIdx(s.ppg_dicrotic, d.ppg_dicrotic, ppgLen);
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

// ---- ONE PLACE THAT CHANGES THE ALIGNMENT ---------------------------------
//
// showPage() because leadsForBinTemplate reads the same selection, so redrawing
// the page swaps every panel's trace to this alignment's average; refreshFocus
// because the close-up is keyed on the same choice and would otherwise sit on
// the previous alignment until the next click.
void TemplateViewerWindow::applyAlignmentSelection(bool force, AnchorType a) {
    m_forceAlign = force;
    m_forcedAlign = a;
    showPage();
    if (m_lastFocusMarker >= 0)
        refreshFocus(m_focusWidget, m_lastFocusBinIdx, m_lastFocusLeadIdx,
            m_lastFocusTemplateIdx, m_lastFocusMarker, m_lastFocusCol);
}

// P -> Q -> R -> J -> P. The ring is spelled out rather than taken from
// anchor_view::kAllAnchors, which is in CSV-merge order (R, P, Q, J): that is
// the order columns are written in and is not the order an operator wants to
// walk a beat in. Changing one must not change the other.
namespace {
    constexpr std::array<AnchorType, 4> kAlignRing = {
        AnchorType::P_ONSET,
        AnchorType::Q_ONSET,
        AnchorType::R_PEAK,
        AnchorType::J_POINT,
    };
    const char* alignRingButton(AnchorType a) {
        switch (a) {
        case AnchorType::P_ONSET: return "p_align_button";
        case AnchorType::Q_ONSET: return "q_align_button";
        case AnchorType::R_PEAK:  return "r_align_button";
        case AnchorType::J_POINT: return "j_point_align_button";
        }
        return "r_align_button";
    }
}

void TemplateViewerWindow::cycleAlignment(int step) {
    const int n = static_cast<int>(kAlignRing.size());

    // Where the ring currently stands. Automatic is not ON the ring, so Tab out
    // of it enters at P going forwards and at J going backwards, rather than
    // silently treating automatic as R and skipping P.
    int idx = -1;
    if (m_forceAlign)
        for (int k = 0; k < n; ++k)
            if (kAlignRing[k] == m_forcedAlign) { idx = k; break; }

    const int next = (idx < 0)
        ? (step > 0 ? 0 : n - 1)
        : ((idx + step) % n + n) % n;
    const AnchorType a = kAlignRing[next];

    // THROUGH THE BUTTON when there is one. Setting the members directly would
    // leave the checked radio naming the previous alignment, and the operator
    // would be reading a label that disagrees with the trace. The toggled
    // handler calls applyAlignmentSelection for us.
    if (QRadioButton* rb = findChild<QRadioButton*>(
        QString::fromLatin1(alignRingButton(a)))) {
        rb->setChecked(true);
        return;
    }
    applyAlignmentSelection(true, a);
}

// Tab / Shift+Tab. Filtered on the application object rather than bound as a
// QShortcut: Tab is consumed by focus navigation inside whichever child holds
// focus, so a window-context shortcut fires only when focus happens to sit on
// something that does not want Tab. Accepting the event here also stops focus
// from moving as a side effect of changing alignment.
bool TemplateViewerWindow::eventFilter(QObject* obj, QEvent* ev) {
    if (ev->type() == QEvent::KeyPress && isActiveWindow()) {
        auto* ke = static_cast<QKeyEvent*>(ev);
        // Backtab is what Qt delivers for Shift+Tab; Key_Tab with the Shift
        // modifier arrives on some platforms, so both are tested.
        const bool back = ke->key() == Qt::Key_Backtab
            || (ke->key() == Qt::Key_Tab
                && (ke->modifiers() & Qt::ShiftModifier));
        if (ke->key() == Qt::Key_Tab || ke->key() == Qt::Key_Backtab) {
            cycleAlignment(back ? -1 : +1);
            return true;   // eaten: no focus change
        }
    }
    return QMainWindow::eventFilter(obj, ev);
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
            applyAlignmentSelection(force, a);
            });
        if (rb->isChecked()) { m_forceAlign = force; m_forcedAlign = a; }
    }
    // Hotkeys: same mapping as the radio buttons above. These go through the
    // BUTTON, not through applyAlignmentSelection, so the checked radio always
    // names the alignment on screen -- pressing Q used to change the trace and
    // leave the R radio checked.
    struct Key { const char* seq; const char* btn; bool force; AnchorType a; };
    static const Key kKeys[] = {
        { "A", "automatic_align_button", false, AnchorType::R_PEAK  },
        { "P", "p_align_button",         true,  AnchorType::P_ONSET },
        { "Q", "q_align_button",         true,  AnchorType::Q_ONSET },
        { "R", "r_align_button",         true,  AnchorType::R_PEAK  },
        { "J", "j_point_align_button",   true,  AnchorType::J_POINT },
    };
    for (const Key& k : kKeys) {
        const char* btn = k.btn;
        const bool force = k.force;
        const AnchorType a = k.a;
        auto* sc = new QShortcut(QKeySequence(QString::fromLatin1(k.seq)),
            this, nullptr, nullptr, Qt::WindowShortcut);
        connect(sc, &QShortcut::activated, this, [this, btn, force, a]() {
            if (QRadioButton* rb = findChild<QRadioButton*>(
                QString::fromLatin1(btn))) {
                rb->setChecked(true);
                return;
            }
            applyAlignmentSelection(force, a);
            });
    }
    // TAB IS NOT A QShortcut. See eventFilter.
    if (qApp) qApp->installEventFilter(this);
}

// ---- THE PULSE ALIGNMENT RADIO GROUP --------------------------------------
//
// These four widgets existed in the .ui with NOTHING behind them: three radios
// and a spin box that changed no state and ran no code, while the pulse was
// always stacked on whatever the build chose. A control that does nothing is
// worse than an absent one -- it is a claim about the waveform beside it.
//
// findChild rather than ui->ppg_foot_align, following wireAlignButtons: a
// widget missing from the .ui leaves that option unavailable instead of
// failing the build, and says so on stderr rather than silently.
void TemplateViewerWindow::wirePpgAlignButtons() {
    // KEYBOARD TRACKING OFF, which is what makes the spin box affordable.
    // valueChanged fires per KEYSTROKE by default, so typing "10" would run a
    // whole-page re-stack at 1 and again at 10 -- the first one on a
    // percentage the operator never chose. With tracking off Qt emits once, on
    // Return or focus-out, and once per arrow click. No debounce timer,
    // because Qt already has the concept.
    if (auto* sp = findChild<QSpinBox*>(QStringLiteral("ppg_align_percent"))) {
        sp->setKeyboardTracking(false);
        sp->setRange(0, 100);
        sp->setSuffix(QStringLiteral(" %"));
        sp->setValue(m_ppgAlignPercent);
        connect(sp, qOverload<int>(&QSpinBox::valueChanged), this,
            [this](int v) {
                // TYPING A PERCENTAGE IS CHOOSING THE PERCENT ALIGNMENT. It
                // was otherwise possible for the box to read 10 while Foot was
                // checked, which is a number on screen describing nothing;
                // syncPpgAlignControls moves the radio to match.
                applyPpgAlignSelection(PpgAlign::Percent, v);
            });
    }
    else {
        fprintf(stderr, "[ppg-align] NOT WIRED: ppg_align_percent\n");
    }

    struct Btn { const char* name; PpgAlign mode; };
    static const Btn kBtns[] = {
        { "ppg_auto_align",    PpgAlign::Auto },
        { "ppg_foot_align",    PpgAlign::Foot },
        { "ppg_specify_align", PpgAlign::Percent },
    };
    for (const Btn& b : kBtns) {
        QRadioButton* rb =
            findChild<QRadioButton*>(QString::fromLatin1(b.name));
        if (!rb) {
            fprintf(stderr, "[ppg-align] NOT WIRED: %s\n", b.name);
            continue;
        }
        const PpgAlign mode = b.mode;
        connect(rb, &QRadioButton::toggled, this, [this, mode](bool on) {
            if (!on) return;              // only the newly-checked one acts
            applyPpgAlignSelection(mode, m_ppgAlignPercent);
            });
        // Sync the member to whichever button Designer has checked, WITHOUT
        // re-stacking: this runs during construction, before any page or panel
        // exists. (ppg_auto_align carries checked="true" in the .ui, so the
        // startup state is "as built" -- the honest default, since the window
        // has not re-aligned anything yet.)
        if (rb->isChecked()) m_ppgAlignMode = mode;
    }
}

// Blocked signals, not a re-entrancy flag: the handlers above call
// applyPpgAlignSelection, which calls this, which would re-fire them. Qt has
// the switch, so use it.
void TemplateViewerWindow::syncPpgAlignControls() {
    struct Btn { const char* name; PpgAlign mode; };
    static const Btn kBtns[] = {
        { "ppg_auto_align",    PpgAlign::Auto },
        { "ppg_foot_align",    PpgAlign::Foot },
        { "ppg_specify_align", PpgAlign::Percent },
    };
    for (const Btn& b : kBtns) {
        QRadioButton* rb =
            findChild<QRadioButton*>(QString::fromLatin1(b.name));
        if (!rb) continue;
        const bool was = rb->blockSignals(true);
        rb->setChecked(b.mode == m_ppgAlignMode);
        rb->blockSignals(was);
    }
    if (auto* sp = findChild<QSpinBox*>(QStringLiteral("ppg_align_percent"))) {
        const bool was = sp->blockSignals(true);
        sp->setValue(m_ppgAlignPercent);
        sp->blockSignals(was);
    }
}