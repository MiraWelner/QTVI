#pragma once
/**
 * @file   record_summary.hpp
 * @brief  Task H (Section 4.8): the per-record diagnostic summary -- the
 *         arrays behind the single-page QA output, built from the outputs of
 *         the earlier tasks and exported as CSV. record_summary_plot.hpp draws
 *         the page from the same struct.
 *
 * ---------------------------------------------------------------------------
 * PANELS AND WHERE EACH COMES FROM
 * ---------------------------------------------------------------------------
 *
 *   chiSq0_vs_beat, chiSqAbs_vs_beat
 *       Task A's ECG SQI (sqi_ecg.hpp, scoreEcgSQI) -- the SAME scoring pass,
 *       not a re-derivation -- on the summary lead. Indexed by record beat
 *       number; NaN where that lead scored no beat (a row no template claims,
 *       or a beat the lead pruned). The other two leads are kept as well, in
 *       chiSq*_by_lead, for the CSV.
 *
 *   rr_vs_pp
 *       RR(n) against PP(n), per beat, in ms: ECG P wave to P wave. RR is
 *       CH1's R-pair interval in recording time. PP is the time from beat n's
 *       P peak to beat n+1's, each found ON THAT BEAT'S OWN WAVEFORM in the
 *       summary lead -- not read off its template, which would make PP just
 *       RR plus a constant. Each P is measured from its own R, so
 *           PP(n) = RR(n) + PR(n+1) - PR(n)
 *       and a changing PR interval moves the point off the diagonal.
 *
 *       Each beat's P peak is FeatureMarks::find_p_peak on that beat,
 *       between 300 ms before R and its template's Q onset. Ectopic templates
 *       are skipped (a PVC has no P). A pair needs both beats' P and a valid
 *       RR, so it never crosses a bin or a splice gap.
 *
 *   rr_vs_qt
 *       PER TEMPLATE, not per beat: this pipeline measures QT on a template,
 *       never on a single beat, so a per-beat QT would be the template's QT
 *       repeated over every member and the panel would be horizontal streaks.
 *       One point per summary-lead template: x = the median RR over the
 *       template's averaged members, y = the template's QT. The detail behind
 *       each point is in qt_templates.
 *
 *       QT is Q onset -> T end on the Q_ONSET alignment, where both bars live
 *       (anchor_view::anchorFor). In order: the operator's bar, else the
 *       viewer's stored detection (*_auto), else the landmarks Task A's
 *       segments were built from on the template's own median. Which one was
 *       used is recorded per template (qt_source).
 *
 *   poincare
 *       RR(n) against RR(n+1), consecutive beats only: a pair is formed only
 *       inside one bin and when both intervals are valid, which also means no
 *       splice gap sits between them.
 *
 *   normalBeats, abnormalBeats
 *       Individual beat tracings on the summary lead, cropped to a common
 *       window around R (kTraceBeforeR_s .. kTraceAfterR_s) so beats from bins
 *       with different R columns superimpose. Taken from the scored beat
 *       matrix, so they are the waveforms Task A scored.
 *         normal   = averaged members of a REGULAR template, rhythm NORMAL
 *         abnormal = members of an ECTOPIC template, or any member whose
 *                    rhythm verdict is PVC / VOTED_PVC
 *       NOISE templates and templates the operator marked bad ECG are in
 *       neither. At most kMaxTracings of each, spread evenly over the record.
 *
 *   stratified
 *       Task J. A stub today (state_stratified.hpp).
 *
 * BEAT NUMBER is the record-wide R-pair ordinal: bins in order, slices in
 * order within a bin. It is the same beat on every lead, so a beat number in
 * this summary names the same heartbeat in every panel.
 *
 * BUILT AT COMMIT, from the operator's final banks: QT needs the operator's
 * bars, and templates are what the operator edited. scoreEcgSQI is run there
 * against those banks, so chi-sq and QT describe the same templates.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "config_file_handling/config.hpp"
#include "template_generation/template_structs.hpp"
#include "peak_finding/peakfinding_structs.hpp"     // output_binfile_data
#include "prep_for_peakfinding/beat_times.hpp"      // Splice, BeatTimes
#include "prep_for_peakfinding/record_sleep.hpp"
#include "fiducial_marker_finding/anchor_view.hpp"  // AnchorType
#include "logging/sqi_ecg.hpp"                      // EcgSQIResult (Task A)
#include "state/state_stratified.hpp"             // Task J (stub)

namespace record_summary {

    // The spec's name for Task J's type (Section 4.8 struct).
    using StateStratifiedFeatures = state_stratified::StateStratifiedFeatures;

    inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    inline constexpr double kTraceBeforeR_s = 0.30;
    inline constexpr double kTraceAfterR_s = 0.60;
    inline constexpr std::size_t kMaxTracings = 200;

    enum class QtSource : uint8_t { NONE = 0, OPERATOR = 1, VIEWER_AUTO = 2, DETECTED = 3 };
    inline const char* qtSourceName(QtSource s) {
        switch (s) {
        case QtSource::OPERATOR:    return "OPERATOR";
        case QtSource::VIEWER_AUTO: return "VIEWER_AUTO";
        case QtSource::DETECTED:    return "DETECTED";
        default:                    return "NONE";
        }
    }

    // One point of the RR-QT panel, with what produced it.
    struct TemplateQT {
        std::size_t bin = 0;
        int slot = 0;
        std::string name;
        int n_beats = 0;          // averaged members
        double rr_ms = kNaN;      // median RR over the averaged members
        double qt_ms = kNaN;
        QtSource qt_source = QtSource::NONE;
        bool abnormal = false;    // ECTOPIC template
    };

    struct DiagnosticSummary {
        std::string recordID;
        int lead = 0;                                    // the summary lead, 0..2

        std::vector<double> chiSq0_vs_beat;              // summary lead, per beat number
        std::vector<double> chiSqAbs_vs_beat;
        std::array<std::vector<double>, 3> chiSq0_by_lead;
        std::array<std::vector<double>, 3> chiSqAbs_by_lead;

        std::vector<std::pair<double, double>> rr_vs_pp;   // ms, per beat: ECG P to P
        std::vector<std::pair<double, double>> rr_vs_qt;   // ms, per template
        std::vector<std::pair<double, double>> poincare;   // ms, RR(n) vs RR(n+1)
        std::vector<TemplateQT> qt_templates;              // behind rr_vs_qt, same order

        std::vector<std::vector<double>> normalBeats;      // superimposed tracings
        std::vector<std::vector<double>> abnormalBeats;
        std::vector<std::size_t> normalBeatNumbers;        // beat number of each tracing
        std::vector<std::size_t> abnormalBeatNumbers;
        double trace_t0_ms = kNaN;                         // time of sample 0, re R
        double trace_dt_ms = kNaN;

        StateStratifiedFeatures stratified;   // from Task J

        std::size_t n_beats = 0;
        std::size_t n_p_waves = 0;                        // beats whose own P was found
        std::vector<double> p_re_r_ms;                    // per beat number: its own P peak, ms re its R (NaN = none)
    };

    // Everything buildSummary reads. Pointers are borrowed; nullptr / empty
    // means that input is absent and its panel stays empty.
    struct SummaryInputs {
        std::string recordID;
        const template_structs::TemplateFile* tmpl = nullptr;
        const template_structs::BeatsFile* beats = nullptr;
        const std::vector<output_binfile_data>* peakResults = nullptr;
        // [bin][lead][row] -> R-pair ordinal (TemplateInfo::ecg_slice_of_row).
        const std::vector<std::array<std::vector<std::size_t>, 3>>* ecgSliceOfRow = nullptr;
        const EcgSQIResult* sqi = nullptr;                 // Task A
        const beat_times::BeatTimes* times = nullptr;
        const record_sleep::RecordSleep* sleep = nullptr;
        double ecgFs = 0.0;          // rate of the beat matrices / templates
        double ecgRecRateHz = 0.0;   // rate the R positions are in
        int lead = 0;                // preferred summary lead
    };

    namespace detail {
        inline double median(std::vector<double> v) {
            v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
            if (v.empty()) return kNaN;
            const std::size_t m = v.size() / 2;
            std::nth_element(v.begin(), v.begin() + m, v.end());
            const double hi = v[m];
            return (v.size() % 2) ? hi : 0.5 * (*std::max_element(v.begin(), v.begin() + m) + hi);
        }

        // Up to `cap` of `n` items, evenly spaced, in order.
        inline std::vector<std::size_t> evenPick(std::size_t n, std::size_t cap) {
            std::vector<std::size_t> out;
            if (n == 0 || cap == 0) return out;
            if (n <= cap) { for (std::size_t i = 0; i < n; ++i) out.push_back(i); return out; }
            for (std::size_t k = 0; k < cap; ++k)
                out.push_back(static_cast<std::size_t>((static_cast<double>(k) + 0.5) * n / cap));
            return out;
        }
    }

    inline DiagnosticSummary buildSummary(const SummaryInputs& in)
    {
        DiagnosticSummary ds;
        ds.recordID = in.recordID;
        if (!in.tmpl || !in.beats || !in.peakResults) return ds;
        const auto& bins = *in.peakResults;
        const std::size_t nBins = std::min(bins.size(), in.tmpl->bins.size());

        // ---- record beat numbering, and RR per beat -------------------------
        // offset[b] + s is beat number of R-pair s in bin b.
        std::vector<std::size_t> offset(nBins + 1, 0);
        for (std::size_t b = 0; b < nBins; ++b) {
            const std::size_t nR = bins[b].ch1.raw.size();
            offset[b + 1] = offset[b] + (nR >= 2 ? nR - 1 : 0);
        }
        ds.n_beats = offset[nBins];
        std::vector<double> rrMs(ds.n_beats, kNaN);
        for (std::size_t b = 0; b < nBins; ++b) {
            const auto& rp = bins[b].ch1.raw;
            if (rp.size() < 2 || !(in.ecgRecRateHz > 0.0)) continue;
            const beat_times::Splice sp(bins[b].ecg_bin_indexs, std::numeric_limits<uint64_t>::max());
            for (std::size_t s = 0; s + 1 < rp.size(); ++s) {
                uint64_t o0 = 0, o1 = 0;
                if (!sp.original(rp[s], o0) || !sp.original(rp[s + 1], o1)) continue;
                // A pair straddling a splice junction is not an interval: its
                // original-record distance differs from its spliced one.
                if (o1 <= o0 || o1 - o0 != static_cast<uint64_t>(rp[s + 1] - rp[s])) continue;
                rrMs[offset[b] + s] = 1000.0 * static_cast<double>(o1 - o0) / in.ecgRecRateHz;
            }
        }

        // ---- Poincare: consecutive valid intervals, WITHIN A BIN -----------
        // Beat numbers run on across a bin boundary, but the last beat of one
        // bin and the first of the next are not neighbours in time, so the
        // pair is formed only inside a bin.
        for (std::size_t b = 0; b < nBins; ++b)
            for (std::size_t n = offset[b]; n + 1 < offset[b + 1]; ++n)
                if (std::isfinite(rrMs[n]) && std::isfinite(rrMs[n + 1]))
                    ds.poincare.push_back({ rrMs[n], rrMs[n + 1] });

        // ---- the summary lead: the preferred one if it scored anything ------
        int lead = std::clamp(in.lead, 0, 2);
        if (in.sqi) {
            std::array<bool, 3> has{ false, false, false };
            for (const EcgSQIBeat& r : in.sqi->beats) has[r.lead] = true;
            if (!has[lead]) for (int c = 0; c < 3; ++c) if (has[c]) { lead = c; break; }
        }
        ds.lead = lead;

        auto sliceOf = [&](std::size_t b, int c, uint32_t row) -> std::size_t {
            if (!in.ecgSliceOfRow || b >= in.ecgSliceOfRow->size()) return SIZE_MAX;
            const auto& v = (*in.ecgSliceOfRow)[b][c];
            return (row < v.size()) ? v[row] : SIZE_MAX;
            };
        auto beatNumber = [&](std::size_t b, int c, uint32_t row) -> std::size_t {
            const std::size_t s = sliceOf(b, c, row);
            if (s == SIZE_MAX || b >= nBins || offset[b] + s >= offset[b + 1]) return SIZE_MAX;
            return offset[b] + s;
            };

        // ---- chi-sq vs beat (Task A) -----------------------------------------
        for (int c = 0; c < 3; ++c) {
            ds.chiSq0_by_lead[c].assign(ds.n_beats, kNaN);
            ds.chiSqAbs_by_lead[c].assign(ds.n_beats, kNaN);
        }
        if (in.sqi) {
            for (const EcgSQIBeat& r : in.sqi->beats) {
                const std::size_t n = beatNumber(r.bin, r.lead, r.row);
                if (n == SIZE_MAX) continue;
                ds.chiSq0_by_lead[r.lead][n] = r.q.chiSq0;
                ds.chiSqAbs_by_lead[r.lead][n] = r.q.chiSqAbs;
            }
        }
        ds.chiSq0_vs_beat = ds.chiSq0_by_lead[lead];
        ds.chiSqAbs_vs_beat = ds.chiSqAbs_by_lead[lead];

        // ---- per template: RR vs QT, and the tracing pools -------------------
        static const char* const keys[] = { "CH1", "CH2", "CH3" };
        const auto bIt = in.beats->per_channel_beats.find(keys[lead]);
        const auto rhIt = in.beats->per_channel_rhythm.find(keys[lead]);
        struct Cand { std::size_t beat, bin; uint32_t row; int rCol; };
        std::vector<double> pRelMs(ds.n_beats, kNaN);   // each beat's own P, ms re its R
        std::vector<Cand> normalPool, abnormalPool;
        const double msPer = (in.ecgFs > 0.0) ? 1000.0 / in.ecgFs : kNaN;

        if (in.sqi) {
            for (const EcgSQITemplate& et : in.sqi->templates) {
                if (et.lead != lead || et.bin >= nBins) continue;
                const tbank::TemplateBank& bank = in.tmpl->bins[et.bin].ecg_bank[lead];
                if (et.slot < 0 || et.slot >= bank.size()) continue;
                const tbank::template_of_all_signals& tp = bank.templates[et.slot];
                const tbank::Category cat = tp.presumedCategory();
                const bool bad = tp.badEcgMarked() || tp.marked_invalid_template;
                const std::vector<uint32_t>& avg = tbank::averagedRows(tp);

                // RR-QT point.
                TemplateQT pt;
                pt.bin = et.bin; pt.slot = et.slot; pt.name = et.name;
                pt.n_beats = static_cast<int>(avg.size());
                pt.abnormal = (cat == tbank::Category::ECTOPIC);
                {
                    std::vector<double> rr;
                    for (const uint32_t r : avg) {
                        const std::size_t n = beatNumber(et.bin, lead, r);
                        if (n != SIZE_MAX) rr.push_back(rrMs[n]);
                    }
                    pt.rr_ms = detail::median(std::move(rr));
                }
                const tbank::BankMarkerSet& m = tp.marks(static_cast<int32_t>(AnchorType::Q_ONSET));
                auto qtFrom = [&](double q, double t) {
                    return (q >= 0.0 && t > q && std::isfinite(msPer)) ? (t - q) * msPer : kNaN;
                    };
                if (m.q_onset >= 0.0 || m.t_end >= 0.0) {
                    // An operator edit on either bar; the other falls back to
                    // its stored detection, as the viewer draws it.
                    pt.qt_ms = qtFrom(m.q_onset >= 0.0 ? m.q_onset : m.q_onset_auto,
                        m.t_end >= 0.0 ? m.t_end : m.t_end_auto);
                    if (std::isfinite(pt.qt_ms)) pt.qt_source = QtSource::OPERATOR;
                }
                if (!std::isfinite(pt.qt_ms)) {
                    pt.qt_ms = qtFrom(m.q_onset_auto, m.t_end_auto);
                    if (std::isfinite(pt.qt_ms)) pt.qt_source = QtSource::VIEWER_AUTO;
                }
                if (!std::isfinite(pt.qt_ms)) {
                    pt.qt_ms = qtFrom(et.qOnset, et.tEnd);
                    if (std::isfinite(pt.qt_ms)) pt.qt_source = QtSource::DETECTED;
                }
                if (!bad && cat != tbank::Category::NOISE) {
                    ds.qt_templates.push_back(pt);
                    if (std::isfinite(pt.rr_ms) && std::isfinite(pt.qt_ms))
                        ds.rr_vs_qt.push_back({ pt.rr_ms, pt.qt_ms });
                }

                // Each beat's own P peak, for RR vs PP: find_p_peak between
                // R - 300 ms and the template's Q onset. A PVC has no P.
                if (cat != tbank::Category::ECTOPIC && bIt != in.beats->per_channel_beats.end()
                    && et.bin < bIt->second.size() && in.ecgFs > 0.0 && et.qOnset > 0) {
                    const auto& mat = bIt->second[et.bin];
                    const double lo = std::max(0.0, et.rCol - 0.300 * in.ecgFs);
                    for (const uint32_t r : tp.members) {
                        if (r >= mat.size()) continue;
                        const std::size_t n = beatNumber(et.bin, lead, r);
                        if (n == SIZE_MAX) continue;
                        const double p = FeatureMarks::find_p_peak(mat[r], lo, et.qOnset, in.ecgFs);
                        if (std::isfinite(p) && p >= 0.0) pRelMs[n] = (p - et.rCol) * 1000.0 / in.ecgFs;
                    }
                }

                // Tracing pools.
                if (bad || cat == tbank::Category::NOISE) continue;
                const std::vector<uint8_t>* rhythm =
                    (rhIt != in.beats->per_channel_rhythm.end() && et.bin < rhIt->second.size())
                    ? &rhIt->second[et.bin] : nullptr;
                std::vector<char> inAvg;
                for (const uint32_t r : avg) { if (r >= inAvg.size()) inAvg.resize(r + 1, 0); inAvg[r] = 1; }
                for (const uint32_t r : tp.members) {
                    const std::size_t n = beatNumber(et.bin, lead, r);
                    if (n == SIZE_MAX) continue;
                    const bool ectopicRhythm = rhythm && r < rhythm->size() && (*rhythm)[r] != 0;
                    const Cand c{ n, et.bin, r, et.rCol };
                    if (cat == tbank::Category::ECTOPIC || ectopicRhythm) abnormalPool.push_back(c);
                    else if (r < inAvg.size() && inAvg[r]) normalPool.push_back(c);
                }
            }
        }

        // ---- RR vs PP: P to P over each R-pair -------------------------------
        for (std::size_t n = 0; n < ds.n_beats; ++n) if (std::isfinite(pRelMs[n])) ++ds.n_p_waves;
        ds.p_re_r_ms = pRelMs;
        for (std::size_t b = 0; b < nBins; ++b)
            for (std::size_t n = offset[b]; n + 1 < offset[b + 1]; ++n) {
                const double rr = rrMs[n];
                if (!std::isfinite(rr) || !std::isfinite(pRelMs[n]) || !std::isfinite(pRelMs[n + 1]))
                    continue;
                const double pp = rr + pRelMs[n + 1] - pRelMs[n];
                if (pp > 0.0) ds.rr_vs_pp.push_back({ rr, pp });
            }

        // ---- tracings: crop around R, evenly over the record ---------------
        if (bIt != in.beats->per_channel_beats.end() && in.ecgFs > 0.0) {
            const int pre = static_cast<int>(std::lround(kTraceBeforeR_s * in.ecgFs));
            const int post = static_cast<int>(std::lround(kTraceAfterR_s * in.ecgFs));
            ds.trace_t0_ms = -1000.0 * pre / in.ecgFs;
            ds.trace_dt_ms = 1000.0 / in.ecgFs;
            auto take = [&](std::vector<Cand>& pool, std::vector<std::vector<double>>& dst,
                std::vector<std::size_t>& nums) {
                    std::sort(pool.begin(), pool.end(), [](const Cand& a, const Cand& b) { return a.beat < b.beat; });
                    pool.erase(std::unique(pool.begin(), pool.end(),
                        [](const Cand& a, const Cand& b) { return a.beat == b.beat; }), pool.end());
                    for (const std::size_t k : detail::evenPick(pool.size(), kMaxTracings)) {
                        const Cand& c = pool[k];
                        if (c.bin >= bIt->second.size() || c.row >= bIt->second[c.bin].size()) continue;
                        const std::vector<double>& beat = bIt->second[c.bin][c.row];
                        std::vector<double> tr(static_cast<std::size_t>(pre + post + 1), kNaN);
                        for (int i = -pre; i <= post; ++i) {
                            const int j = c.rCol + i;
                            if (j >= 0 && j < static_cast<int>(beat.size())) tr[i + pre] = beat[j];
                        }
                        dst.push_back(std::move(tr));
                        nums.push_back(c.beat);
                    }
                };
            take(normalPool, ds.normalBeats, ds.normalBeatNumbers);
            take(abnormalPool, ds.abnormalBeats, ds.abnormalBeatNumbers);
        }

        // ---- Task J ----------------------------------------------------------
        ds.stratified = state_stratified::stratify(*in.tmpl, in.times, in.sleep);
        return ds;
    }

    // ---- THE SPEC'S SIGNATURE (Section 4.8), AS WRITTEN ---------------------
    //
    // buildSummary(bins, sqis, rr, qt): one BeatSQI per beat, and rr / qt as
    // parallel per-beat vectors, the spec's code line for line. `Bin` is
    // BinTemplates, the per-bin record in this codebase; the spec's body does
    // not read it, so neither does this one.
    //
    // TWO DEPARTURES, both so the arrays stay correct:
    //   * a pair with a non-finite member is skipped rather than pushed, so a
    //     beat with no RR or no QT does not plot as NaN;
    //   * the Poincare pair needs both rr[i] and rr[i+1] finite, so a caller
    //     marks a bin boundary or splice gap with NaN and no pair crosses it.
    //
    // The record pipeline (commit) uses buildSummary(SummaryInputs) above,
    // which also fills RR vs PP, the tracings, and per-template QT.
    inline DiagnosticSummary buildSummary(const std::vector<template_structs::BinTemplates>& /*bins*/,
        const std::vector<BeatSQI>& sqis,
        const std::vector<double>& rr,
        const std::vector<double>& qt)
    {
        DiagnosticSummary ds;
        for (std::size_t i = 0; i < sqis.size(); ++i) {
            ds.chiSq0_vs_beat.push_back(sqis[i].chiSq0);
            ds.chiSqAbs_vs_beat.push_back(sqis[i].chiSqAbs);
        }
        for (std::size_t i = 0; i < rr.size(); ++i) {
            if (i < qt.size() && std::isfinite(rr[i]) && std::isfinite(qt[i]))
                ds.rr_vs_qt.push_back({ rr[i], qt[i] });
            if (i + 1 < rr.size() && std::isfinite(rr[i]) && std::isfinite(rr[i + 1]))
                ds.poincare.push_back({ rr[i], rr[i + 1] });
        }
        ds.n_beats = std::max(sqis.size(), rr.size());
        return ds;
    }

    // ---- CSV EXPORT ---------------------------------------------------------
    //
    // <dir>/<stem>_diagnostic_summary.csv, long format, every panel's arrays:
    //   panel,series,index,x,y
    //     chisq0 / chisqabs   series CHn, index = x = beat number
    //     rr_vs_pp            index = point, x = RR ms, y = PP ms
    //     poincare            index = point, x = RR(n) ms, y = RR(n+1) ms
    //     rr_vs_qt            series = template, index = bin, x = RR ms, y = QT ms
    //     normal_beats / abnormal_beats
    //                         series = beat number, index = sample, x = ms re R
    // <dir>/<stem>_diagnostic_summary_qt.csv: one row per template behind rr_vs_qt.
    inline bool writeSummaryCsv(const std::string& dir, const std::string& stem,
        const DiagnosticSummary& ds)
    {
        const std::string path = dir + "/" + stem + "_diagnostic_summary.csv";
        std::ofstream f(path);
        if (!f.is_open()) {
            std::cerr << "  WARNING: could not open " << path << " for the record summary\n";
            return false;
        }
        f.precision(10);
        f << "panel,series,index,x,y\n";
        static const char* const keys[] = { "CH1", "CH2", "CH3" };
        for (int c = 0; c < 3; ++c) {
            for (std::size_t n = 0; n < ds.chiSq0_by_lead[c].size(); ++n)
                if (std::isfinite(ds.chiSq0_by_lead[c][n]))
                    f << "chisq0," << keys[c] << ',' << n << ',' << n << ',' << ds.chiSq0_by_lead[c][n] << '\n';
            for (std::size_t n = 0; n < ds.chiSqAbs_by_lead[c].size(); ++n)
                if (std::isfinite(ds.chiSqAbs_by_lead[c][n]))
                    f << "chisqabs," << keys[c] << ',' << n << ',' << n << ',' << ds.chiSqAbs_by_lead[c][n] << '\n';
        }
        for (std::size_t i = 0; i < ds.rr_vs_pp.size(); ++i)
            f << "rr_vs_pp,," << i << ',' << ds.rr_vs_pp[i].first << ',' << ds.rr_vs_pp[i].second << '\n';
        for (std::size_t n = 0; n < ds.p_re_r_ms.size(); ++n)
            if (std::isfinite(ds.p_re_r_ms[n]))
                f << "p_wave,," << n << ',' << n << ',' << ds.p_re_r_ms[n] << '\n';
        for (std::size_t i = 0; i < ds.poincare.size(); ++i)
            f << "poincare,," << i << ',' << ds.poincare[i].first << ',' << ds.poincare[i].second << '\n';
        for (const TemplateQT& t : ds.qt_templates)
            if (std::isfinite(t.rr_ms) && std::isfinite(t.qt_ms))
                f << "rr_vs_qt," << t.name << ',' << t.bin << ',' << t.rr_ms << ',' << t.qt_ms << '\n';
        auto traces = [&](const char* panel, const std::vector<std::vector<double>>& tr,
            const std::vector<std::size_t>& nums) {
                for (std::size_t k = 0; k < tr.size(); ++k)
                    for (std::size_t i = 0; i < tr[k].size(); ++i)
                        if (std::isfinite(tr[k][i]))
                            f << panel << ',' << nums[k] << ',' << i << ','
                            << ds.trace_t0_ms + ds.trace_dt_ms * static_cast<double>(i) << ',' << tr[k][i] << '\n';
            };
        traces("normal_beats", ds.normalBeats, ds.normalBeatNumbers);
        traces("abnormal_beats", ds.abnormalBeats, ds.abnormalBeatNumbers);

        const std::string qpath = dir + "/" + stem + "_diagnostic_summary_qt.csv";
        std::ofstream q(qpath);
        if (!q.is_open()) {
            std::cerr << "  WARNING: could not open " << qpath << "\n";
            return false;
        }
        q << "bin,channel,template,slot,n_beats,rr_ms,qt_ms,qt_source,abnormal\n";
        for (const TemplateQT& t : ds.qt_templates)
            q << t.bin << ',' << keys[ds.lead] << ',' << t.name << ',' << t.slot << ',' << t.n_beats << ','
            << t.rr_ms << ',' << t.qt_ms << ',' << qtSourceName(t.qt_source) << ',' << int(t.abnormal) << '\n';
        return true;
    }

}  // namespace record_summary