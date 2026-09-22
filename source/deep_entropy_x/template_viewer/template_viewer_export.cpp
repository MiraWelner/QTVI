// ========================================================================
// template_viewer_export.cpp
// Normalization references and every CSV / bin the Save path writes.
//
// One of five units; TemplateViewerWindow is declared in template_viewer.hpp.
// ========================================================================

#include "template_viewer.hpp"


void TemplateViewerWindow::compute_global_refs() {
    /*compute the ecg global reference value for QRS complex height(abs(R) + abs(S))) and pulse global ref for
    PPG / ART / ART_PULM(abs(peak) - abs(foot))*/
    for (int c = 0; c < 3; ++c) {
        m_ecgGlobalRef[c] = normalize_features::compute_ecg_global_ref(m_bins, c, m_sampleRate);
    }
    for (int c = 0; c < 4; ++c) {
        m_pulseGlobalRef[c] = normalize_features::compute_pulse_global_ref(m_bins, c);
    }
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
            // ecgRef to R even inside a non-R part.
            const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
            const auto& ecg = chs[ch]->ecgTemplate_raw;
            if (ecg.empty()) continue;
            // slotMarks selects the lead, so the per-lead subscripts below
            const tbank::BankMarkerSet& rmk =
                b.slotMarks(ch, 0, AnchorType::R_PEAK);
            const FeatureMarks::ReactiveEcg rx = FeatureMarks::reactive_ecg(
                ecg, rmk.p_begin, rmk.q_onset, rmk.s_end, rmk.t_end, m_sampleRate);
            EcgFeatures f = computeEcgFeatures(ecg, rx.p_peak, rmk.q_onset,
                b.r_peak_ch[ch], rmk.s_end, rmk.t_end, m_sampleRate);
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

std::vector<double> TemplateViewerWindow::normalize_ppg_or_similar(const std::vector<double>& raw, double footIdx, int pulseChan) const {
    const double ref = (pulseChan >= 0 && pulseChan < 4) ? m_pulseGlobalRef[pulseChan] : std::nan("");
    return normalize_features::normalize_pulse_trace(raw, footIdx, ref);
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

// Merge an ordered list of in-memory CSV parts into one canonical file.
// parts[0] is the base (kept whole, with its keys); each subsequent part
// contributes only its value columns (first 3 key columns stripped), appended
// to every row. Done ONCE at the end over N parts, instead of the growing
// read-modify-write per pass this replaced. Row counts must match across all
// parts. Returns true on success.
//
// Each part is (label, content): the label names the alignment for diagnostics
// only, since the parts are no longer files on disk. Previously each part was
// staged as a temp file, re-read, rewritten as a sidecar, and read a third time
// here; the whole save now happens in memory and only the canonical file is
// written.
struct CsvPart { std::string label; std::string content; };

static bool mergeCsvParts(const std::string& canonicalPath,
    const std::vector<CsvPart>& parts)
{
    if (parts.empty()) return false;

    auto splitLines = [](const std::string& s) {
        std::vector<std::string> lines; std::string ln;
        std::istringstream in(s);
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

    std::vector<std::string> merged = splitLines(parts[0].content);
    if (merged.empty()) {
        fprintf(stderr, "[tmplcsv] merge: base part empty: %s\n", parts[0].label.c_str());
        return false;
    }
    for (size_t s = 1; s < parts.size(); ++s) {
        std::vector<std::string> next = splitLines(parts[s].content);
        if (next.size() != merged.size()) {
            fprintf(stderr, "[tmplcsv] merge: row mismatch (%zu vs %zu) for %s -- skipping\n",
                next.size(), merged.size(), parts[s].label.c_str());
            continue;   // skip a bad part rather than abort the whole merge
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

// Companion to <id>_bins.csv: the fitted curve for every ECG landmark, so the
// exact model and parameters that placed each anchor are on disk next to the
// templates. PEAKS (p/q/r-peak) get the weighted quadratic from
// symmetricExtremumFit -- coeff ascending in (t - seed): value = p0 + p1*(t-seed)
// + p2*(t-seed)^2. ONSETS/OFFSETS (p_begin, q_onset, s_end, t_end) get the
// curve_fit model BIC selected -- params meaning depends on curve_type
// (LINEAR {m1,c1,m2,c2,brk} / SIGMOID {a,k,t0,c} / FRACTIONAL {c0,c1,c2,p1,p2,lo}
// / FLAT {mean}). Recomputed here from the R-aligned average over a +-100 ms
// window (same convention as the boundary log), not stored during detection.
// Position of `a` in anchor_view::kAllAnchors, or -1. The export cache is
// indexed by this rather than by AnchorType's numeric value so its storage is
// a fixed four-element array that exists in full before any worker runs.
int TemplateViewerWindow::anchorSlot(AnchorType a)
{
    for (std::size_t i = 0; i < anchor_view::anchor_array.size(); ++i)
        if (anchor_view::anchor_array[i] == a) return static_cast<int>(i);
    return -1;
}

// Fill the export landmark cache for every (bin, lead, alignment).
//
// Called once at the start of the Finish export, before the CSV parts. Same
// per-bin concurrency argument as seedOneBin: detect_template_landmarks and
// everything under curve_fit are pure functions of a trace plus scalars, each
// task writes only its own [bin] row, and every row exists before the first
// task starts. Nothing here reads another bin.
//
// NOT primed at load. The export is the only caller that wants all four
// alignments for all bins at once; the grid only ever needs the page it is
// drawing, and reactiveGlyphs already caches that per panel.
void TemplateViewerWindow::primeExportLandmarks()
{
    const std::size_t n = m_bins.size();
    for (auto& per : m_exportLm) per.assign(n, {});
    if (n == 0 || !(m_sampleRate > 0.0)) return;

    std::vector<std::size_t> idx(n);
    for (std::size_t i = 0; i < n; ++i) idx[i] = i;

    const double fs = m_sampleRate;
    QtConcurrent::blockingMap(idx, [this, fs](std::size_t bi) {
        const TemplateBin& bin = m_bins[bi];
        for (std::size_t s = 0; s < anchor_view::anchor_array.size(); ++s) {
            const AnchorType a = anchor_view::anchor_array[s];
            for (int c = 0; c < 3; ++c)
                m_exportLm[s][bi][static_cast<std::size_t>(c)] =
                alignedLandmarks(bin, c, a, fs);
        }
        });
}

// Cache read. Returns an invalid (valid == false) set for anything out of
// range, which is exactly what alignedLandmarks itself returned for a missing
// lead or an empty template, so callers need no new branch.
const FeatureMarks::TemplateLandmarks&
TemplateViewerWindow::exportLandmarks(std::size_t bi, int lead,
    AnchorType a) const
{
    static const FeatureMarks::TemplateLandmarks kNone;
    const int s = anchorSlot(a);
    if (s < 0 || lead < 0 || lead > 2) return kNone;
    const auto& per = m_exportLm[static_cast<std::size_t>(s)];
    if (bi >= per.size()) return kNone;
    return per[bi][static_cast<std::size_t>(lead)];
}

void TemplateViewerWindow::writeLandmarkFitsCsv(const std::string& dir) {
    const std::string subj = m_subjectId.toStdString();
    std::ofstream f(dir + "/" + subj + "_landmark_fits.csv", std::ios::trunc);
    if (!f) {
        std::cerr << "[landmark-fits] could not write "
            << dir << "/" << subj << "_landmark_fits.csv\n";
        return;
    }
    f << "file_id,bin,channel,landmark,curve_type,seed_or_lo,rss,"
        "p0,p1,p2,p3,p4,p5\n";
    if (!(m_sampleRate > 0.0)) { std::cout << "Wrote landmark fits CSV (empty)\n"; return; }

    const int half = std::max(2, static_cast<int>(std::lround(0.100 * m_sampleRate)));
    const double NaNv = std::numeric_limits<double>::quiet_NaN();

    auto peakTypeName = [](subsample_refine::PeakCurveType t) -> const char* {
        switch (t) {
        case subsample_refine::PeakCurveType::QUADRATIC:  return "QUADRATIC";
        case subsample_refine::PeakCurveType::CUBIC:      return "CUBIC";
        case subsample_refine::PeakCurveType::FIVE_POINT: return "FIVE_POINT";
        default:                                      return "SEED";
        }
        };
    auto transTypeName = [](curve_fit::FitType t) -> const char* {
        switch (t) {
        case curve_fit::FitType::LINEAR:       return "LINEAR";
        case curve_fit::FitType::SIGMOID:      return "SIGMOID";
        case curve_fit::FitType::FRACTIONAL:   return "FRACTIONAL";
        case curve_fit::FitType::CUBIC_SPLINE: return "CUBIC_SPLINE";
        case curve_fit::FitType::CUBIC:        return "CUBIC";
        default:                               return "FLAT";
        }
        };

    for (size_t bi = 0; bi < m_bins.size(); ++bi) {
        const TemplateBin& b = m_bins[bi];
        for (int c = 0; c < 3; ++c) {
            const std::vector<double>& ecg = b.chFor(c, AnchorType::R_PEAK).ecgTemplate_raw;
            if (ecg.empty()) continue;
            const FeatureMarks::TemplateLandmarks& aa =
                exportLandmarks(bi, c, AnchorType::R_PEAK);
            const int N = static_cast<int>(ecg.size());

            auto emitRow = [&](const char* name, const char* ctype,
                double seedOrLo, double rss, const std::vector<double>& ps) {
                    f << subj << ',' << bi << ',' << (c + 1) << ',' << name << ','
                        << ctype << ',';
                    if (seedOrLo >= 0.0) f << seedOrLo; f << ',';
                    if (!std::isnan(rss)) f << rss; f << ',';
                    for (int k = 0; k < 6; ++k) {
                        if (k < static_cast<int>(ps.size())) f << ps[k];
                        if (k < 5) f << ',';
                    }
                    f << '\n';
                };

            // PEAKS: the selected peak model (Auto = BIC quad-vs-cubic),
            // per-marker sigma -- so the CSV reports what Save placed.
            // halfWidth IS A PARAMETER NOW, not the global constant. This CSV
            // is the record of what placed each mark, so it has to fit the
            // same span the detector fitted -- a hardcoded +-7 here would
            // report a different polynomial from the one on screen for P.
            // Returns the fit whole, so the reported model and the position
            // written below come from one computation. sigma and halfWidth are
            // resolved inside placeEcgPeak from the shared table.
            auto peak = [&](const char* name, EcgPeak which, double pos) {
                const subsample_refine::peak_fit fit =
                    placeEcgPeak(ecg, which, pos, m_peakFitMode);
                if (fit.position < 0.0) {
                    emitRow(name, "NONE", -1.0, NaNv, {});
                    return fit;
                }
                emitRow(name, peakTypeName(fit.type),
                    static_cast<double>(fit.seed), fit.rss,
                    { fit.coeff[0], fit.coeff[1], fit.coeff[2], fit.coeff[3] });
                return fit;
                };
            // ONSETS/OFFSETS: the selected transition model (Auto = BIC) over
            // +-100 ms.
            auto trans = [&](const char* name, double pos) {
                if (pos < 0.0 || pos > N - 1) { emitRow(name, "NONE", -1.0, NaNv, {}); return; }
                const int lo = std::max(0, static_cast<int>(std::lround(pos)) - half);
                const int hi = std::min(N - 1, static_cast<int>(std::lround(pos)) + half);
                if (hi - lo < 5) { emitRow(name, "NONE", (double)lo, NaNv, {}); return; }
                const curve_fit::FitResult fit = curve_fit::selectBestFit(ecg, lo, hi, m_onOffsetFitMode);
                emitRow(name, transTypeName(fit.type), (double)lo, fit.rss, fit.params);
                };

            // Kept: the position columns below read these fits, not aa.*.
            const auto plP = peak("p_peak", EcgPeak::P, aa.p_peak);
            const auto plQ = peak("q_peak", EcgPeak::Q, aa.q_peak);
            const auto plR = peak("r_peak", EcgPeak::R, aa.r_peak);
            (void)plP; (void)plQ; (void)plR;   // read by the emitters below
            trans("p_begin", aa.p_begin);
            trans("q_onset", aa.q_onset);
            trans("s_end", aa.s_end);
            trans("t_end", aa.t_end);
        }
    }
    std::cout << "Wrote landmark fits CSV: " << dir << "/"
        << subj << "_landmark_fits.csv\n";
}

std::string TemplateViewerWindow::buildAlignedTemplateCsv(AnchorType anchor) {
    if (m_bins.empty()) return {};
    std::ostringstream f;

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
        "p_begin", "p_peak", "q_onset", "q_peak", "r_peak", "s_peak",
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
            const FeatureMarks::TemplateLandmarks& aaF =
                exportLandmarks(bi, c, anchor);
            // NO ROUNDING: computeEcgFeatures takes doubles, AnchorAuto is
            // double, and the qrs/qt milliseconds this feeds are sub-sample.
            ftAuto[c] = computeEcgFeatures(ecg,
                aaF.p_peak, aaF.q_onset, aaF.r_peak,
                aaF.s_end, aaF.t_end, m_sampleRate);
            // Per lead, because slotMarks selects the lead -- the old bin-wide
            // MarkerSet held all three leads in one object and was fetched once
            // per bin.
            // THE ONE BAR SET, TRANSLATED INTO THIS BLOCK'S FRAME. There are
            // four bars in the record and they live in two cells: p_begin in
            // P_ONSET, q_onset / s_end / t_end in Q_ONSET. userMarks fetches
            // them and adds frameShift(owner -> anchor), so every block reports
            // the same marks expressed in its own alignment's columns.
            //
            // WHAT EACH BLOCK CARRIES is still whatever showsBar admits, which
            // is the point: R reports all four, Q reports q_onset / s_end /
            // t_end, P reports p_begin, J reports t_end. The intervals below
            // therefore resolve under R and Q, which hold q_onset, s_end and
            // t_end, and come back absent under P and J.
            // R AND J REPORT THEIR OWN CELLS, P AND Q THE CANONICAL SET.
            // Those two alignments carry their own measurement of a landmark
            // (anchor_view::ownsCanonicalBar), so a _user column in their block
            // means "placed against THIS waveform" and must not be the shared
            // bar translated in. P and Q own the shared bars, so userMarks --
            // which reads the owner cells and converts -- is the right answer
            // for them and for the R-block interval columns.
            tbank::BankMarkerSet umk = anchor_view::hasOwnBars(anchor)
                ? b.slotMarks(c, 0, anchor)
                : b.userMarks(c, 0, anchor);
            // userMarks returns BARS ONLY, and BankMarkerSet no longer has a
            // p_peak field at all: P peak is a glyph. It is recomputed from the
            // two bars that bracket it, on this alignment's own waveform.
            const FeatureMarks::ReactiveEcg rxF = FeatureMarks::reactive_ecg(
                ecg, umk.p_begin, umk.q_onset, umk.s_end, umk.t_end, m_sampleRate);
            ftUser[c] = computeEcgFeatures(ecg,
                rxF.p_peak, umk.q_onset, b.r_peak_ch[c],
                umk.s_end, umk.t_end, m_sampleRate);
        }

        // ECG marker positions, indices aligned with ECG_MARKERS.
        double ecgAuto[3][kNumEcgMarkers], ecgUser[3][kNumEcgMarkers];
        for (int c = 0; c < 3; ++c) {
            // Glyphs AND bars from THIS alignment: the glyphs from its own
            // detection, the bars from its own cells. No frame conversion --
            // both were measured on this alignment's average.
            // SECOND READ OF THE SAME THING. This used to be a second
            // alignedLandmarks() call with identical arguments to aaF above,
            // in a second loop over the same leads of the same bin -- so every
            // bin paid for its detections twice. Now both are cache reads.
            const FeatureMarks::TemplateLandmarks& aa =
                exportLandmarks(bi, c, anchor);
            // Same owner-cell read as the block above: one bar set, this
            // alignment's columns. See the note there.
            // R AND J REPORT THEIR OWN CELLS, P AND Q THE CANONICAL SET.
            // Those two alignments carry their own measurement of a landmark
            // (anchor_view::ownsCanonicalBar), so a _user column in their block
            // means "placed against THIS waveform" and must not be the shared
            // bar translated in. P and Q own the shared bars, so userMarks --
            // which reads the owner cells and converts -- is the right answer
            // for them and for the R-block interval columns.
            tbank::BankMarkerSet umk = anchor_view::hasOwnBars(anchor)
                ? b.slotMarks(c, 0, anchor)
                : b.userMarks(c, 0, anchor);
            const std::vector<double>& ecgA = b.chFor(c, anchor).ecgTemplate_raw;
            const FeatureMarks::ReactiveEcg rxA = FeatureMarks::reactive_ecg(
                ecgA, (int)std::lround(aa.p_begin), (int)std::lround(aa.q_onset),
                (int)std::lround(aa.s_end), (int)std::lround(aa.t_end), m_sampleRate);
            const FeatureMarks::ReactiveEcg rxU = FeatureMarks::reactive_ecg(ecgA, umk.p_begin, umk.q_onset, umk.s_end, umk.t_end, m_sampleRate);

            // Same placement the glyph and the focus panel draw: aa.r_peak is
            // a seed, placeEcgPeak turns it into a position. On a failed fit
            // the seed stands, so a found landmark is never blanked.
            //
            // ecgA, not ecg: this block's own alignment, whose columns these
            // numbers are expressed in.
            const auto plRa = placeEcgPeak(ecgA, EcgPeak::R, aa.r_peak,
                m_peakFitMode);
            const double rPeakOut = (plRa.position >= 0.0) ? plRa.position
                : aa.r_peak;

            ecgAuto[c][0] = aa.p_begin;
            ecgAuto[c][1] = rxA.p_peak;          // reactive glyph, detector brackets
            ecgAuto[c][2] = aa.q_onset;
            ecgAuto[c][3] = ftAuto[c].q_idx;
            ecgAuto[c][4] = rPeakOut;
            ecgAuto[c][5] = ftAuto[c].s_idx;
            ecgAuto[c][6] = aa.s_end;
            ecgAuto[c][7] = aa.t_end;
            ecgUser[c][0] = umk.p_begin;
            ecgUser[c][1] = rxU.p_peak;          // reactive glyph, operator brackets
            ecgUser[c][2] = umk.q_onset;
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
        // ALL DOUBLE. The (double) casts on the pulse bars are gone with
        // TemplateBin's int fields, and the arterial arrays were narrowing
        // fifteen sub-sample positions apiece.
        const double ppgUser[kNumPpgMarkers] = {
            b.ppg_onset, rxPpgUser.t50, b.ppg_peak,
            b.ppg_dicrotic, b.ppg_peak2, rxPpgUser.t80,
            b.ppg_end };
        const double abpAuto[kNumArterialMarkers] = { b.abp_onset_auto, b.abp_peak_auto,
            b.abp_dicrotic_auto, b.abp_peak2_auto, b.abp_end_auto };
        const double abpUser[kNumArterialMarkers] = { b.abp_onset, b.abp_peak,
            b.abp_dicrotic, b.abp_peak2, b.abp_end };
        const double artAuto[kNumArterialMarkers] = { b.art_onset_auto, b.art_peak_auto,
            b.art_dicrotic_auto, b.art_peak2_auto, b.art_end_auto };
        const double artUser[kNumArterialMarkers] = { b.art_onset, b.art_peak,
            b.art_dicrotic, b.art_peak2, b.art_end };
        const double artpAuto[kNumArterialMarkers] = { b.art_pulm_onset_auto, b.art_pulm_peak_auto,
            b.art_pulm_dicrotic_auto, b.art_pulm_peak2_auto, b.art_pulm_end_auto };
        const double artpUser[kNumArterialMarkers] = { b.art_pulm_onset, b.art_pulm_peak,
            b.art_pulm_dicrotic, b.art_pulm_peak2, b.art_pulm_end };

        // Derived T peak for the autodetect glyph group, bracketed by the AUTO J-point and T-end. The J-point is the left bracket because no T-onset
        double tPeakAutoGlyph[3];
        for (int gc = 0; gc < 3; ++gc)
            tPeakAutoGlyph[gc] = FeatureMarks::find_t_peak(
                chs[gc]->ecgTemplate_raw,
                b.s_end_auto_ch[gc], b.t_end_auto_ch[gc]);

        auto emitVal = [&](const std::vector<double>& v, int j) {
            f << ',';
            if (j >= 0 && j < (int)v.size() && !std::isnan(v[j])) f << v[j];
            };
        // Emit ",1" if the row is this marker's row, ",<blank>" otherwise.
        //normally you don't want to round the marker, but in this case it is a 1 hot encoding so you have to 
        auto emitLoc = [&](double markerIdx, int row) {
            f << ',';
            if (markerIdx >= 0.0 && (int)std::lround(markerIdx) == row) f << '1';
            };
        const double toMs = 1000.0 / m_sampleRate;
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
                emitLoc(b.q_onset_auto_ch[gc], row);
                emitLoc(b.r_peak_auto_ch[gc], row);
                emitLoc(tPeakAutoGlyph[gc], row);
            }
            for (const auto& gl : ppg_and_artpulse_automated_markers)
                emitLoc(b.*gl.idx, row);
            f << '\n';
        }
    }

    // Serialize and suffix the value columns with the alignment tag. Row-key
    // columns (file_id/bin_num/x_ms) are never suffixed, so the parts line up
    // on the merge in save_bin_and_csv. Returned as a string -- the caller
    // merges it directly instead of staging a sidecar file.
    std::string content = f.str();
    const size_t nl = content.find('\n');
    const std::string header = content.substr(0, nl);
    const std::string body = content.substr(nl);   // includes the leading '\n'
    // From the ARGUMENT, not from window state: this function is called once
    // per alignment inside a single save.
    const std::string suffix = std::string("_") + anchor_view::label(anchor);
    return suffixValueColumns(header, suffix) + body;
}

// AnchorType -> short name for the boundary log's `anchor` column.
static const char* anchorName_boundary(AnchorType a) {
    switch (a) {
    case AnchorType::P_PEAK: return "P_ONSET";
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
        case Landmark::Q_ONSET:  return (int)std::lround(tb.q_onset_auto_ch[lead]);
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
            // QRS duration (ms) = distance between q_onset and s_end glyphs.
            double qrsMs = 0.0;
            {
                const int q = (int)std::lround(b.q_onset_auto_ch[lead]);
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
                // ROUNDED HERE, DELIBERATELY: the record's confirmedIndex is a
                // sample offset into rec.segment, which has no fractional
                // representation. The touched store keeps the sub-sample click
                // position; only this one field quantises.
                int confirmed = -1;
                const int mid = markerIdOf(tgt.lm);
                // THE CANONICAL CELL. This record has no alignment of its own
                // -- it logs a segment of the R-framed signal -- so it reports
                // the confirmation of the bar's canonical copy, which is the
                // one the Automatic view draws. A bar confirmed only on a
                // non-canonical alignment reads as unconfirmed here.
                auto it = m_touchedMarks.find(
                    touchKey(i, lead, mid, anchor_view::anchorFor(mid)));
                if (it != m_touchedMarks.end())
                    confirmed = static_cast<int>(std::lround(it->second)) - lo;

                boundary_training::BoundaryTrainingRecord rec;
                rec.segment.assign(sig.begin() + lo, sig.begin() + hi);
                rec.confirmedIndex = confirmed;
                // fit from curve_fit: fit-and-select on this landmark's window.
                const curve_fit::FitResult fit = curve_fit::selectBestFit(sig, lo, hi - 1);
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

    // NO RE-SEED BEFORE SAVE. The writers call alignedLandmarks /
    // detect_template_landmarks with m_onOffsetFitMode and m_peakFitMode
    // directly, so the selected models are already honoured; the pass that
    // was here wrote the same numbers into the cells first.

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

        // ---- ONE WRITE, NO MERGE ------------------------------------
        //
        // The writer emits all four alignment blocks and the pulse block into
        // one row itself. This replaced four suffixed EcgOnly parts stitched by
        // mergeCsvParts, which compared row COUNTS and concatenated line i of
        // each -- so what a row meant was a contract between five files.
        const QString csvPath = csvDir.absolutePath() + "/" + m_subjectId + "_template_markings.csv";

        {
            std::ofstream mf(csvPath.toStdString(), std::ios::trunc);
            if (!mf)
                throw std::runtime_error("cannot open for write: " + csvPath.toStdString());
            // `anchor` is unused for EcgAndPulse: the writer walks
            // kAllAnchors itself. R_PEAK only satisfies the signature.
            writeTemplateMarkingsCsv(mf, m_bins,
                m_subjectId.toStdString(), m_sampleRate,
                AnchorType::R_PEAK, MarkingsCsvSection::EcgAndPulse,
                m_peakFitMode, m_onOffsetFitMode);
            if (!mf.good())
                throw std::runtime_error("failed writing " + csvPath.toStdString());
            std::cout << "Wrote markings CSV: " << csvPath.toStdString() << "\n";
        }

        // ---- WHAT THE OPERATOR ACTUALLY RULED ON ---------------------
        //
        // ITS OWN FILE, because templates.csv cannot carry it. That one is
        // written by GenerateTemplatesFast inside prepareViewerJob, BEFORE
        // this window exists, so its `confirmed` column could only ever read
        // "presumed" -- for every template in every record, however much
        // marking followed. It could not be fixed by rewriting either: the
        // ChannelBlocks that writer takes hold raw pointers into
        // GenerateTemplatesFast's own frame, which is gone by the time the
        // operator finishes.
        //
        // This is written from m_bins, after the session, which is the only
        // place and time the answer exists. Joined to templates.csv on
        // (bin, channel, template) -- three rows that file already carries.
        //
        // `template` IS THE SLOT INDEX, not the letter. templates.csv's
        // `template` row holds the NAME (PQRST_A), and the letter comes from
        // tbank::letterRanks over the surviving templates -- so joining on the
        // name would mean recomputing those ranks here and keeping the two in
        // step. The slot index is the same key on both sides with nothing to
        // recompute.
        {
            const QString cPath = csvDir.absolutePath() + "/" + m_subjectId + "_template_confirmations.csv";
            std::ofstream cf(cPath.toStdString(), std::ios::trunc);

            cf << "file_id,bin,channel,template,state,n_members\n";
            static const char* kChan[4] = { "CH1", "CH2", "CH3", "PPG" };
            for (size_t i = 0; i < m_bins.size(); ++i) {
                const TemplateBin& b = m_bins[i];
                for (int c = 0; c < 4; ++c) {
                    const tbank::TemplateBank& bk =
                        (c < 3) ? b.ecg_bank[c] : b.ppg_bank;
                    for (int t = 0; t < bk.size(); ++t) {
                        const tbank::BankTemplate& tp = bk.templates[t];
                        // CROSSED-OUT IS TESTED FIRST. You have to view a
                        // template to cross it out, so both flags are set
                        // and the rejection is the later and more specific
                        // statement; the other order reports every
                        // rejected template as confirmed.
                        const char* st =
                            (tp.marked_invalid_template != 0) ? "unconfirmed"
                            : tp.confirmed() ? "confirmed"
                            : "presumed";
                        cf << m_subjectId.toStdString() << ',' << i << ','
                            << kChan[c] << ',' << t << ',' << st << ','
                            << tp.memberCount() << '\n';
                    }
                }
            }
        }
    }
    catch (const std::exception& e) {
        QMessageBox::critical(this, "Save failed",
            QString("Could not write markings for %1:\n\n%2\n\n"
                "If the CSV is open in Excel, close it and try again.")
            .arg(m_subjectId, e.what()));
        return;   // don't emit finished(); let the user retry
    }

    // Aligned-template CSV: one part per alignment, holding that alignment's
    // own averages, merged into the canonical <id>_template.csv in one write.
    // Same restructuring as the markings parts above, same reason -- the
    // sidecars only existed to survive window teardowns between passes.
    {
        QDir alignedDir(m_templateDir);
        if (!alignedDir.exists()) alignedDir.mkpath(".");
        const QString canonical = alignedDir.filePath(m_subjectId + "_bins.csv");


        // ONCE, ACROSS CORES, BEFORE ANY PART IS BUILT. Every landmark the
        // five CSV passes below report is a pure function of a stored average
        // and its R column, so they all read this instead of re-detecting.
        QGuiApplication::setOverrideCursor(Qt::WaitCursor);
        primeExportLandmarks();
        QGuiApplication::restoreOverrideCursor();

        std::vector<CsvPart> parts;
        for (AnchorType a : anchor_view::anchor_array) {
            std::string content = buildAlignedTemplateCsv(a);
            if (content.empty()) continue;
            parts.push_back(CsvPart{ anchor_view::label(a), std::move(content) });
        }
        if (!parts.empty() && mergeCsvParts(canonical.toStdString(), parts)) {
            std::cout << "Wrote bins CSV: " << canonical.toStdString() << "\n";
        }

        writeLandmarkFitsCsv(alignedDir.absolutePath().toStdString());
    }
    emit finished();
}