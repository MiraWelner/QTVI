/**
 * @file   EnsembleTemplate.hpp
 * @brief  Build an ensemble-averaged template from segmented waveform data.
 *         Port of EnsembleTemplate.m
 *
 *         Steps: segment -> remove bad lengths -> remove bad peak amplitudes ->
 *         remove bad peak positions -> align -> wave-score pruning -> median.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-26
 */
#pragma once

#include "TemplateTypes.hpp"
#include "AlignWaves.hpp"

 // Simple local-maxima finder (no min-distance constraint).
 // Only used here for ECG alignment — not the full R-peak findpeaks.
static inline void et_findpeaks_simple(const vector<double>& data,
    vector<double>& pks,
    vector<size_t>& locs) {
    pks.clear();
    locs.clear();
    if (data.size() < 3) return;
    for (size_t i = 1; i < data.size() - 1; ++i) {
        if (!std::isnan(data[i]) && data[i] > data[i - 1] && data[i] >= data[i + 1]) {
            // Walk through plateau
            size_t j = i;
            while (j < data.size() - 1 && data[j] == data[j + 1]) ++j;
            if (j < data.size() - 1 && data[j] > data[j + 1]) {
                pks.push_back(data[i]);
                locs.push_back(i + (j - i) / 2);
                i = j;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Local helpers (static to avoid ODR issues in header-only)
// ---------------------------------------------------------------------------
static inline double et_nanmean(const vector<double>& v) {
    double s = 0; size_t c = 0;
    for (double x : v) if (!std::isnan(x)) { s += x; ++c; }
    return c > 0 ? s / c : NaN;
}

static inline double et_nanstd(const vector<double>& v) {
    double m = et_nanmean(v);
    if (std::isnan(m)) return NaN;
    double ss = 0; size_t c = 0;
    for (double x : v) if (!std::isnan(x)) { ss += (x - m) * (x - m); ++c; }
    return c > 1 ? std::sqrt(ss / (c - 1)) : 0.0;
}

static inline double et_nanmedian(const vector<double>& v) {
    vector<double> tmp;
    for (double x : v) if (!std::isnan(x)) tmp.push_back(x);
    if (tmp.empty()) return NaN;
    std::sort(tmp.begin(), tmp.end());
    size_t n = tmp.size();
    return (n % 2 == 0) ? (tmp[n / 2 - 1] + tmp[n / 2]) / 2.0 : tmp[n / 2];
}

static inline vector<double> et_col_nanmedian(const vector<vector<double>>& mat, size_t cols) {
    vector<double> result(cols, NaN);
    for (size_t c = 0; c < cols; ++c) {
        vector<double> cv;
        for (const auto& row : mat) {
            if (c < row.size() && !std::isnan(row[c])) cv.push_back(row[c]);
        }
        if (!cv.empty()) {
            std::sort(cv.begin(), cv.end());
            size_t n = cv.size();
            result[c] = (n % 2 == 0) ? (cv[n / 2 - 1] + cv[n / 2]) / 2.0 : cv[n / 2];
        }
    }
    return result;
}

static inline void et_col_stats(const vector<vector<double>>& mat, size_t cols,
    vector<double>& cmean, vector<double>& cstd) {
    cmean.assign(cols, NaN);
    cstd.assign(cols, NaN);
    for (size_t c = 0; c < cols; ++c) {
        vector<double> v;
        for (const auto& row : mat) {
            if (c < row.size() && !std::isnan(row[c])) v.push_back(row[c]);
        }
        if (!v.empty()) {
            cmean[c] = et_nanmean(v);
            cstd[c] = et_nanstd(v);
        }
    }
}

static inline size_t et_valid_len(const vector<double>& v) {
    size_t l = v.size();
    while (l > 0 && std::isnan(v[l - 1])) --l;
    return l;
}

static inline void et_compute_peaks(const vector<vector<double>>& segs,
    vector<double>& peak_vals,
    vector<size_t>& peak_pos) {
    peak_vals.resize(segs.size());
    peak_pos.resize(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        double mx = -Inf; size_t mi = 0;
        for (size_t j = 0; j < segs[i].size(); ++j) {
            if (!std::isnan(segs[i][j]) && segs[i][j] > mx) {
                mx = segs[i][j]; mi = j;
            }
        }
        peak_vals[i] = mx;
        peak_pos[i] = mi;
    }
}

// Remove rows where mask[i]==true from all parallel vectors
template <typename... Vecs>
static inline void et_remove_masked(const vector<bool>& mask, Vecs&... vecs) {
    auto filter = [&](auto& v) {
        size_t write = 0;
        for (size_t i = 0; i < v.size(); ++i) {
            if (!mask[i]) v[write++] = std::move(v[i]);
        }
        v.resize(write);
        };
    (filter(vecs), ...);
}

// ---------------------------------------------------------------------------
// EnsembleTemplate
// ---------------------------------------------------------------------------
inline vector<double> EnsembleTemplate(
    const vector<double>& wave,
    const vector<size_t>& segment_idxs,
    double std_multiplier,
    const string& type,
    const vector<size_t>& expand = {})
{
    if (segment_idxs.size() < 2) return {};

    size_t n_segs = segment_idxs.size() - 1;
    vector<std::pair<size_t, size_t>> ranges(n_segs);

    for (size_t i = 0; i < n_segs; ++i) {
        size_t s = segment_idxs[i];
        size_t e = segment_idxs[i + 1];
        if (!expand.empty() && i < expand.size()) {
            s = (s > expand[i]) ? s - expand[i] : 0;
            e = std::min(e + expand[i], wave.size() - 1);
        }
        ranges[i] = { s, e };
    }

    vector<size_t> seg_lengths(n_segs);
    for (size_t i = 0; i < n_segs; ++i) {
        seg_lengths[i] = ranges[i].second - ranges[i].first + 1;
    }

    // ---- 1. Remove bad lengths ----
    vector<double> sl_d(seg_lengths.begin(), seg_lengths.end());
    double med_len = et_nanmedian(sl_d);
    double std_len = et_nanstd(sl_d);
    double len_min = med_len - std_multiplier * std_len;
    double len_max = med_len + std_multiplier * std_len;

    size_t max_seg_len = 0;
    for (size_t i = 0; i < n_segs; ++i) {
        double sl = static_cast<double>(seg_lengths[i]);
        if (sl >= len_min && sl <= len_max && seg_lengths[i] > max_seg_len)
            max_seg_len = seg_lengths[i];
    }
    if (max_seg_len == 0) return {};

    vector<vector<double>> good_segs;
    vector<size_t> good_lens;

    for (size_t i = 0; i < n_segs; ++i) {
        double sl = static_cast<double>(seg_lengths[i]);
        if (sl < len_min || sl > len_max) continue;

        vector<double> seg(max_seg_len + 1, NaN);
        for (size_t j = 0; j < seg_lengths[i] && j <= max_seg_len; ++j) {
            size_t idx = ranges[i].first + j;
            if (idx < wave.size()) seg[j] = wave[idx];
        }
        good_segs.push_back(std::move(seg));
        good_lens.push_back(seg_lengths[i]);
    }
    if (good_segs.empty()) return {};

    // ---- 2. Remove bad peak amplitudes ----
    {
        vector<double> pv; vector<size_t> pp;
        et_compute_peaks(good_segs, pv, pp);
        double avg = et_nanmean(pv);
        double sd = et_nanstd(pv);
        double lo = avg - std_multiplier * 2.0 * sd;
        double hi = avg + std_multiplier * 2.0 * sd;

        vector<bool> bad(good_segs.size(), false);
        for (size_t i = 0; i < good_segs.size(); ++i)
            bad[i] = (pv[i] < lo || pv[i] > hi);
        et_remove_masked(bad, good_segs, good_lens);
    }
    if (good_segs.empty()) return {};

    // ---- 3. Remove bad peak positions ----
    {
        vector<double> pv; vector<size_t> pp;
        et_compute_peaks(good_segs, pv, pp);
        vector<double> pp_d(pp.begin(), pp.end());
        double avg = et_nanmean(pp_d);
        double sd = et_nanstd(pp_d);
        double lo = avg - std_multiplier * 3.0 * sd;
        double hi = avg + std_multiplier * 3.0 * sd;

        vector<bool> bad(good_segs.size(), false);
        for (size_t i = 0; i < good_segs.size(); ++i)
            bad[i] = (static_cast<double>(pp[i]) < lo || static_cast<double>(pp[i]) > hi);
        et_remove_masked(bad, good_segs, good_lens);
    }
    if (good_segs.empty()) return {};

    // ---- 4. Compute alignment points ----
    vector<size_t> align_pts(good_segs.size());

    if (type == "ppg") {
        for (size_t i = 0; i < good_segs.size(); ++i) {
            size_t vl = et_valid_len(good_segs[i]);
            double mx = -Inf; size_t mi = 0;
            for (size_t j = 0; j + 1 < vl; ++j) {
                double dv = good_segs[i][j + 1] - good_segs[i][j];
                if (dv > mx) { mx = dv; mi = j; }
            }
            align_pts[i] = mi;
        }
    }
    else {
        for (size_t i = 0; i < good_segs.size(); ++i) {
            double seg_med = et_nanmedian(good_segs[i]);
            double seg_std = et_nanstd(good_segs[i]);
            double min_h = seg_med + seg_std * 2.5;

            vector<double> pks;
            vector<size_t> locs;
            et_findpeaks_simple(good_segs[i], pks, locs);

            bool found = false;
            for (size_t k = 0; k < pks.size(); ++k) {
                if (pks[k] >= min_h) {
                    align_pts[i] = locs[k];
                    found = true;
                    break;
                }
            }
            if (!found) {
                vector<double> pv; vector<size_t> pp;
                et_compute_peaks({ good_segs[i] }, pv, pp);
                align_pts[i] = pp[0];
            }
        }
    }

    // ---- 5. Align ----
    AlignWavesResult aligned = AlignWaves(good_segs, align_pts);
    auto& aw = aligned.alignedWaves;
    if (aw.empty()) return {};
    size_t ncols = aw[0].size();

    // ---- 6. Wave-score pruning ----
    vector<double> cmean, cstd;
    et_col_stats(aw, ncols, cmean, cstd);

    vector<double> waveScore(aw.size(), 0.0);
    for (size_t i = 0; i < aw.size(); ++i) {
        for (size_t c = 0; c < ncols; ++c) {
            double v = aw[i][c];
            if (std::isnan(v) || std::isnan(cmean[c]) || std::isnan(cstd[c])) continue;
            double dev = std::abs(v - cmean[c]);
            if (dev > 4.0 * cstd[c]) waveScore[i] += 9.0;
            else if (dev > 3.0 * cstd[c]) waveScore[i] += 3.0;
            else if (dev > 2.5 * cstd[c]) waveScore[i] += 1.0;
        }
    }

    for (size_t i = 0; i < aw.size(); ++i) {
        if (i < good_lens.size() && good_lens[i] > 0)
            waveScore[i] /= static_cast<double>(good_lens[i]);
    }

    vector<vector<double>> final_segs;
    vector<size_t> final_lens;
    vector<size_t> final_align;
    for (size_t i = 0; i < aw.size(); ++i) {
        if (waveScore[i] <= 0.30) {
            final_segs.push_back(aw[i]);
            if (i < good_lens.size()) final_lens.push_back(good_lens[i]);
            if (i < align_pts.size()) final_align.push_back(align_pts[i]);
        }
    }
    if (final_segs.empty()) return {};

    // ---- 7. Build template ----
    size_t alignment_point = 0;
    for (auto a : align_pts) if (a > alignment_point) alignment_point = a;

    vector<double> fa_d(final_align.begin(), final_align.end());
    vector<double> fl_d(final_lens.begin(), final_lens.end());
    double med_align = et_nanmedian(fa_d);
    double med_flen = et_nanmedian(fl_d);

    size_t beginpos = 0;
    if (alignment_point > static_cast<size_t>(std::round(med_align)))
        beginpos = alignment_point - static_cast<size_t>(std::round(med_align));

    size_t endpos = alignment_point +
        static_cast<size_t>(std::round(med_flen - med_align));
    size_t total_cols = final_segs.empty() ? 0 : final_segs[0].size();
    if (endpos >= total_cols) endpos = total_cols - 1;
    if (beginpos > endpos) return {};

    size_t tmpl_len = endpos - beginpos + 1;
    vector<vector<double>> sub(final_segs.size());
    for (size_t i = 0; i < final_segs.size(); ++i) {
        sub[i].resize(tmpl_len);
        for (size_t c = 0; c < tmpl_len; ++c)
            sub[i][c] = final_segs[i][beginpos + c];
    }

    return et_col_nanmedian(sub, tmpl_len);
}