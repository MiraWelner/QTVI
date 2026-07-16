#pragma once
//
// Per-subject amplitude normalization for template features.
//
// ECG:
//   RS_peak(bin, ch) = |R_peak_y| + |S_peak_y|
//   Global_Ref_ecg(ch) = median over bins of RS_peak(bin, ch)
//   For every ECG amplitude feature F on that channel:
//     F_norm_abs(bin, ch) = F(bin, ch) / Global_Ref_ecg(ch)
//
// PPG / arterial (ABP / ART / ART_PULM) -- same shape, PI-based:
//   PI(bin, chan) = 100 * (systolic_peak_y - onset_y) / onset_y
//     where onset_y is the foot (the "diastolic trough" of the anchor pulse).
//   Global_Ref_pulse(chan) = median over bins of PI(bin, chan).
//   For every pulse amplitude feature F on that channel:
//     Feature_Local_Ratio(bin, chan) = 100 * (F(bin, chan) - onset_y(bin, chan))
//                                    / onset_y(bin, chan)
//     F_norm_abs(bin, chan) = Feature_Local_Ratio(bin, chan) / Global_Ref_pulse(chan)
//
// Bins are excluded from the median AND from the output whenever they
// aren't usable for that channel:
//   ECG: bad_segment OR bad_r_ch[ch] OR channel template empty
//   pulse (PPG): bad_segment OR ppg_issue != 0 OR ppgTemplate empty
//   pulse (arterial): bad_segment OR <chan>_issue != 0 OR template empty
//
// Output: a per-subject CSV `<subject>_template_markings_normalized.csv`
// alongside the raw `_template_markings.csv`. One row per bin, same file_id
// / bin_index columns, then the normalized amplitude columns. Cells left
// blank when the bin is excluded from that channel.
//

#include "TemplateBinIO.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace normalize_features {

    // Median of a copy (nth_element based). NaN entries are skipped.
    // Returns NaN if no finite entries.
    inline double median_finite(std::vector<double> v) {
        v.erase(std::remove_if(v.begin(), v.end(),
            [](double x) { return !std::isfinite(x); }), v.end());
        if (v.empty()) return std::nan("");
        const size_t mid = v.size() / 2;
        std::nth_element(v.begin(), v.begin() + mid, v.end());
        const double a = v[mid];
        if (v.size() % 2 == 1) return a;
        // Even count: pair with max of the lower half.
        double b = -std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < mid; ++i) if (v[i] > b) b = v[i];
        return 0.5 * (a + b);
    }

    // Look up sample y-value at marker index, NaN if unavailable.
    inline double sample_y(const std::vector<double>& v, int idx) {
        if (idx < 0 || idx >= static_cast<int>(v.size())) return std::nan("");
        const double y = v[idx];
        return std::isnan(y) ? std::nan("") : y;
    }

    // ------------------------------------------------------------------
    //This finds Global_Ref_person# = median [abs(R_peak) + abs(S_peak)]
    // ------------------------------------------------------------------
    inline double compute_ecg_global_ref(const std::vector<TemplateBin>& bins,
        int ch, double sampleRateHz)
    {
        std::vector<double> vals;
        vals.reserve(bins.size());
        for (const auto& b : bins) {
            if (b.bad_segment) continue;
            if (b.bad_r_ch[ch]) continue;
            const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
            const auto& ecg = chs[ch]->ecgTemplate_raw;
            if (ecg.empty()) continue;

            EcgFeatures f = computeEcgFeatures(ecg,
                b.p_peak_ch[ch], b.q_begin_ch[ch], b.r_peak_ch[ch],
                b.s_end_ch[ch], b.t_begin_ch[ch], b.t_end_ch[ch], sampleRateHz);
            const double ry = sample_y(ecg, f.r_idx);
            const double sy = sample_y(ecg, f.s_idx);
            if (std::isnan(ry) || std::isnan(sy)) continue;
            vals.push_back(std::abs(ry) + std::abs(sy));
        }
        return median_finite(std::move(vals));
    }

    // ------------------------------------------------------------------
    // Pulse per-channel Global_Ref (median PI).
    // PI(bin) = (peak_y - foot_y) / |foot_y| * 100
    // ------------------------------------------------------------------
    struct PulseChannel {
        const std::vector<double>* trace;
        int   foot_idx;
        int   peak_idx;
        uint8_t issue;   // 0 = ok, 1 = user-bad, 2 = absent
    };

    inline PulseChannel pulseChan(const TemplateBin& b, int which) {
        // which: 0=PPG, 1=ABP, 2=ART, 3=ART_PULM
        switch (which) {
        case 0: return { &b.ppgTemplate,     b.ppg_onset,    b.ppg_peak,    b.ppg_issue };
        case 1: return { &b.abpTemplate,     b.abp_onset,    b.abp_peak,    b.abp_issue };
        case 2: return { &b.artTemplate,     b.art_onset,    b.art_peak,    b.art_issue };
        default: return { &b.artPulmTemplate, b.art_pulm_onset, b.art_pulm_peak, b.art_pulm_issue };
        }
    }

    inline double compute_pulse_global_ref(const std::vector<TemplateBin>& bins,
        int which)
    {
        std::vector<double> vals;
        vals.reserve(bins.size());
        for (const auto& b : bins) {
            if (b.bad_segment) continue;
            const PulseChannel pc = pulseChan(b, which);
            if (pc.issue != 0) continue;
            if (pc.trace->empty()) continue;
            const double foot_y = sample_y(*pc.trace, pc.foot_idx);
            const double peak_y = sample_y(*pc.trace, pc.peak_idx);
            if (std::isnan(foot_y) || std::isnan(peak_y)) continue;
            if (std::abs(foot_y) < 1e-12) continue;   // divide-by-zero guard
            vals.push_back(100.0 * (peak_y - foot_y) / foot_y);
        }
        return median_finite(std::move(vals));
    }

    // Emit ",value" or "," (blank) for a possibly-NaN value.
    inline void emit_csv_value(std::ofstream& f, double v) {
        f << ',';
        if (std::isfinite(v)) f << v;
    }

    // ------------------------------------------------------------------
    // CSV writer.
    // Columns per bin:
    //   file_id, bin_index,
    //   -- ECG (per channel, using RS-median normalization) --
    //   p_peak_ch{N}_y_norm_abs, q_begin_ch{N}_y_norm_abs, q_peak_ch{N}_y_norm_abs,
    //   r_peak_ch{N}_y_norm_abs, s_peak_ch{N}_y_norm_abs, s_end_ch{N}_y_norm_abs,
    //   t_begin_ch{N}_y_norm_abs, t_end_ch{N}_y_norm_abs
    //   -- pulse (PI-based) --
    //   {chan}_onset_y_norm_abs, {chan}_peak_y_norm_abs,
    //   {chan}_dicrotic_y_norm_abs, {chan}_peak2_y_norm_abs, {chan}_end_y_norm_abs
    //     for chan in {ppg, abp, art, art_pulm}
    // Blanks where the bin was excluded for that channel.
    // ------------------------------------------------------------------
    inline void writeNormalizedCsv(const std::string& path,
        const std::vector<TemplateBin>& bins,
        const std::string& fileID,
        double sampleRateHz)
    {
        std::ofstream f(path);
        if (!f.is_open())
            throw std::runtime_error("cannot open for write: " + path);
        f << std::setprecision(6);

        // Per-channel global refs.
        double ecgRef[3];
        for (int c = 0; c < 3; ++c) ecgRef[c] = compute_ecg_global_ref(bins, c, sampleRateHz);
        double pulseRef[4];
        for (int c = 0; c < 4; ++c) pulseRef[c] = compute_pulse_global_ref(bins, c);

        // Header.
        f << "file_id,bin_index";
        // 8 ECG columns per channel: p_peak, q_begin, q_peak (computed),
        // r_peak, s_peak (computed), s_end, t_peak, t_end.
        const char* ecgNames[] = {
            "p_peak","q_begin","q_peak","r_peak","s_peak",
            "s_end","t_begin","t_end"
        };
        for (int c = 1; c <= 3; ++c)
            for (const char* n : ecgNames)
                f << ',' << n << "_ch" << c << "_y_norm_abs";
        const char* pulseNames[] = { "onset","peak","dicrotic","peak2","end" };
        const char* pulseChans[] = { "ppg","abp","art","art_pulm" };
        for (const char* chan : pulseChans)
            for (const char* n : pulseNames)
                f << ',' << chan << '_' << n << "_y_norm_abs";
        f << ",ecg_ch1_global_ref,ecg_ch2_global_ref,ecg_ch3_global_ref,"
            "ppg_global_ref,abp_global_ref,art_global_ref,art_pulm_global_ref";
        f << '\n';

        for (const auto& b : bins) {
            f << fileID << ',' << b.index;

            // ECG per-channel amplitudes.
            const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
            for (int c = 0; c < 3; ++c) {
                const auto& ecg = chs[c]->ecgTemplate_raw;
                const double ref = ecgRef[c];
                const bool usable = !b.bad_segment && !b.bad_r_ch[c]
                    && !ecg.empty() && std::isfinite(ref) && ref > 0.0;

                double vals[8] = {
                    std::nan(""), std::nan(""), std::nan(""), std::nan(""),
                    std::nan(""), std::nan(""), std::nan(""), std::nan("")
                };
                if (usable) {
                    EcgFeatures ft = computeEcgFeatures(ecg,
                        b.p_peak_ch[c], b.q_begin_ch[c], b.r_peak_ch[c],
                        b.s_end_ch[c], b.t_begin_ch[c], b.t_end_ch[c], sampleRateHz);
                    const int idxs[8] = {
                        b.p_peak_ch[c], b.q_begin_ch[c], ft.q_idx, b.r_peak_ch[c], ft.s_idx,
                        b.s_end_ch[c],  b.t_begin_ch[c], b.t_end_ch[c]
                    };
                    for (int k = 0; k < 8; ++k) {
                        const double y = sample_y(ecg, idxs[k]);
                        if (std::isfinite(y)) vals[k] = y / ref;
                    }
                }
                for (int k = 0; k < 8; ++k) emit_csv_value(f, vals[k]);
            }

            // Pulse channels (PPG, ABP, ART, ART_PULM).
            for (int pc = 0; pc < 4; ++pc) {
                const PulseChannel p = pulseChan(b, pc);
                const double ref = pulseRef[pc];
                const bool usable = !b.bad_segment && p.issue == 0
                    && !p.trace->empty() && std::isfinite(ref) && ref != 0.0;

                double vals[5] = {
                    std::nan(""), std::nan(""), std::nan(""), std::nan(""), std::nan("")
                };
                if (usable) {
                    const double foot_y = sample_y(*p.trace, p.foot_idx);
                    if (std::isfinite(foot_y) && std::abs(foot_y) >= 1e-12) {
                        const int idxs[5] = {
                            p.foot_idx, p.peak_idx,
                            (pc == 0 ? b.ppg_dicrotic
                             : pc == 1 ? b.abp_dicrotic
                             : pc == 2 ? b.art_dicrotic
                             : b.art_pulm_dicrotic),
                            (pc == 0 ? b.ppg_peak2
                             : pc == 1 ? b.abp_peak2
                             : pc == 2 ? b.art_peak2
                             : b.art_pulm_peak2),
                            (pc == 0 ? b.ppg_end
                             : pc == 1 ? b.abp_end
                             : pc == 2 ? b.art_end
                             : b.art_pulm_end)
                        };
                        for (int k = 0; k < 5; ++k) {
                            const double y = sample_y(*p.trace, idxs[k]);
                            if (!std::isfinite(y)) continue;
                            const double local_ratio =
                                100.0 * (y - foot_y) / foot_y;
                            vals[k] = local_ratio / ref;
                        }
                    }
                }
                for (int k = 0; k < 5; ++k) emit_csv_value(f, vals[k]);
            }

            // Global reference values (constant per subject, repeated per row
            // so the file stays a simple flat table).
            for (int c = 0; c < 3; ++c) emit_csv_value(f, ecgRef[c]);
            for (int c = 0; c < 4; ++c) emit_csv_value(f, pulseRef[c]);
            f << '\n';
        }
    }

}   // namespace normalize_features