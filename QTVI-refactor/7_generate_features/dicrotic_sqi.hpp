#pragma once
// ═══════════════════════════════════════════════════════════════════════════════
// dicrotic_sqi.hpp — Exact MATLAB match of dumbDicrotic.m and PPG_SQI.m
//
// Key fixes vs old C++:
//  - dumbDicrotic SP search: MATLAB increments x AFTER setting found=1,
//    so the break happens correctly. The C++ now matches this.
//  - EP narrowing: uses MATLAB-style normalization for orthoginal_dist_thresh.
//  - PPG_SQI: OpenMP parallelization on the inner beat loop (matching parfor).
// ═══════════════════════════════════════════════════════════════════════════════
#include "ppg_features.hpp"

namespace ppg {

    // ─── Dicrotic Notch Detection (exact dumbDicrotic.m) ────────────────────────
    inline int dumbDicrotic(const std::vector<double>& beat_in, double sp_ratio) {
        if (beat_in.size() < 4) return -1;

        // MATLAB: beat = nanfastsmooth(beat, 15)
        std::vector<double> beat = nanfastsmooth(beat_in, 15);
        if (beat.empty()) return -1;

        // MATLAB: [~,mx] = max(beat)
        int mx = 0;
        for (int i = 1; i < (int)beat.size(); i++)
            if (beat[i] > beat[mx]) mx = i;

        int x1 = mx;
        int x2 = (int)beat.size() - 1; // MATLAB: length(beat)

        // MATLAB: if x1 == x2 | x2-x1 < 2
        if (x1 == x2 || x2 - x1 < 2) return -1;

        int pmax = x1, pmin = x2;
        int p_min_dpdt_region, init_EP;

        if (!std::isnan(sp_ratio)) {
            // MATLAB: endlen = length(beat)-pmax
            int endlen = (int)beat.size() - pmax;
            int notch_est = pmax + (int)std::round(endlen * sp_ratio);
            p_min_dpdt_region = pmax + (int)std::round((notch_est - pmax) * (2.0 / 3.0));
            if (p_min_dpdt_region < pmax + 2)
                p_min_dpdt_region = pmax + (int)std::round((pmin - pmax) / 3.0);
            init_EP = notch_est + (int)std::round(((int)beat.size() - notch_est) * 0.5);
            if (init_EP <= p_min_dpdt_region)
                init_EP = pmax + (int)std::round((pmin - pmax) * 0.75);
        }
        else {
            p_min_dpdt_region = pmax + (int)std::round((pmin - pmax) / 3.0);
            init_EP = pmax + (int)std::round((pmin - pmax) * 0.75);
        }

        p_min_dpdt_region = std::min(p_min_dpdt_region, (int)beat.size() - 1);
        init_EP = std::min(init_EP, (int)beat.size() - 1);

        // MATLAB: slice = beat(pmax:p_min_dpdt_region)
        if (p_min_dpdt_region <= pmax) return -1;
        std::vector<double> slice(beat.begin() + pmax, beat.begin() + p_min_dpdt_region + 1);
        if (slice.size() < 2) return -1;

        // MATLAB: [~,p_min_dpdt] = min(diff(slice))
        int p_min_dpdt_local = 0;
        double min_diff = slice[1] - slice[0];
        for (int i = 1; i < (int)slice.size() - 1; i++) {
            double d = slice[i + 1] - slice[i];
            if (d < min_diff) { min_diff = d; p_min_dpdt_local = i; }
        }
        int p_min_dpdt = p_min_dpdt_local + pmax;

        // MATLAB: p_half = beat(pmax) - (beat(pmax) - beat(p_min_dpdt))/2
        double p_half = beat[pmax] - (beat[pmax] - beat[p_min_dpdt]) / 2.0;

        // MATLAB: potential_sp — sorted by abs(slice - p_half)
        std::vector<std::pair<double, int>> potential_sp;
        for (int i = 0; i < (int)slice.size(); i++)
            potential_sp.push_back({ std::abs(slice[i] - p_half), i + pmax });
        std::sort(potential_sp.begin(), potential_sp.end());

        // MATLAB SP search:
        // x = 1; found = 0;
        // while x < length(potential_sp) && found == 0
        //     SP = potential_sp(x,1);
        //     transform = shear_transform(SP:init_EP, beat(SP:init_EP));
        //     if sum(transform < beat(SP:init_EP))/length(beat(SP:init_EP)) < .5
        //         found = 1;
        //     end
        //     x = x+1;
        // end
        int SP = pmax;
        {
            int x = 0;
            bool found = false;
            while (x < (int)potential_sp.size() && !found) {
                int candidate = potential_sp[x].second;
                if (candidate < pmax || candidate >= init_EP || candidate >= (int)beat.size()) {
                    x++;
                    continue;
                }
                std::vector<int> xr;
                std::vector<double> yv;
                for (int k = candidate; k <= init_EP && k < (int)beat.size(); k++) {
                    xr.push_back(k);
                    yv.push_back(beat[k]);
                }
                if (yv.size() >= 2) {
                    auto tf = shear_transform(xr, yv);
                    int below = 0;
                    for (size_t k = 0; k < tf.size(); k++)
                        if (tf[k] < yv[k]) below++;
                    if ((double)below / yv.size() < 0.5) {
                        SP = candidate;
                        found = true;
                    }
                }
                x++;
            }
        }

        // MATLAB EP narrowing:
        // EP = init_EP;
        // while EP > SP+1
        //     shearpressure = beat(SP:EP); time = SP:EP;
        //     c = polyfit([time(end),time(1)], [shearpressure(1), shearpressure(end)], 1);
        //     ... normalize and check orthoginal_dist_thresh
        int EP = init_EP;
        while (EP > SP + 1) {
            int seg_len = EP - SP + 1;
            if (seg_len < 2) break;

            // MATLAB: shearpressure = beat(SP:EP), time = SP:EP
            double sp_first = beat[SP], sp_last = beat[EP];
            double t_first = (double)SP, t_last = (double)EP;

            // MATLAB normalization
            double mn = sp_first, mx_sp = sp_first;
            for (int k = SP; k <= EP; k++) {
                if (beat[k] < mn) mn = beat[k];
                if (beat[k] > mx_sp) mx_sp = beat[k];
            }
            double range_p = mx_sp - mn;
            if (range_p == 0) break;

            double time_range = t_last - t_first;
            if (time_range == 0) break;

            // MATLAB: c = polyfit([time(end),time(1)], [shearpressure(1), shearpressure(end)], 1);
            // m = c(2)/(c(1))
            // shearline: m * x + b for x=1:length
            // Then normalize all three: norm_line, norm_press, norm_time
            // Actually in the MATLAB code it does:
            //   c = polyfit([time(end),time(1)], [shearpressure(1), shearpressure(end)], 1);
            //   m = c(2)/c(1); b = 0;
            //   for x=1:length(shearpressure), shearline(x) = m*x + b; end
            //   norm_line = (shearline-min)/(max-min)
            //   norm_press = (shearpressure-min)/(max-min)
            //   norm_time = (time-min)/(max-min)
            // The polyfit with 2 points gives a line: y = c(1)*x + c(2)
            // where x=time(end) → y=shearpressure(1), x=time(1) → y=shearpressure(end)
            // So: slope_c1 = (sp(1)-sp(end))/(time(end)-time(1))
            //     intercept_c2 = sp(end) - slope_c1 * time(1)
            // Then m = c(2)/c(1), and shearline(x) = m*x (for x=1..len)

            double slope_c1 = (sp_first - sp_last) / (t_last - t_first);
            double intercept_c2 = sp_last - slope_c1 * t_first;
            // MATLAB: m = c(2)/c(1) — this is the shear ratio
            double m_val = (slope_c1 != 0) ? intercept_c2 / slope_c1 : 0;

            // MATLAB: shearline(x) = m * x + 0 for x = 1:length(shearpressure)
            int nn = seg_len;
            std::vector<double> shearline(nn), norm_line(nn), norm_press(nn), norm_time(nn);

            double sl_min = 1e18, sl_max = -1e18;
            for (int k = 0; k < nn; k++) {
                shearline[k] = m_val * (k + 1);  // x = 1-based
                if (shearline[k] < sl_min) sl_min = shearline[k];
                if (shearline[k] > sl_max) sl_max = shearline[k];
            }
            double sl_range = sl_max - sl_min;

            for (int k = 0; k < nn; k++) {
                norm_press[k] = (beat[SP + k] - mn) / range_p;
                norm_time[k] = (double)k / (nn - 1);
                norm_line[k] = (sl_range > 0) ? (shearline[k] - sl_min) / sl_range : 0;
            }

            if (orthoginal_dist_thresh(norm_time, norm_line, norm_press, 0.3))
                EP--;
            else
                break;
        }

        if (SP >= EP || SP >= (int)beat.size() || EP >= (int)beat.size()) return -1;

        // MATLAB: transform = shear_transform(SP:EP, beat(SP:EP))
        std::vector<int> fr;
        std::vector<double> fv;
        for (int k = SP; k <= EP; k++) { fr.push_back(k); fv.push_back(beat[k]); }
        auto tf = shear_transform(fr, fv);

        // MATLAB: [~, min_shear] = min(transform); min_shear = min_shear+SP-1;
        int min_shear_local = 0;
        for (int i = 1; i < (int)tf.size(); i++)
            if (tf[i] < tf[min_shear_local]) min_shear_local = i;
        int min_shear = min_shear_local + SP;

        // MATLAB: [~,start_diastolic_relax] = max(beat(min_shear:pmin))
        int start_dr = min_shear;
        for (int i = min_shear; i <= pmin && i < (int)beat.size(); i++)
            if (beat[i] > beat[start_dr]) start_dr = i;

        // MATLAB: start_diastolic_relax = start_diastolic_relax + min_shear-1
        // Already computed correctly above.

        if (min_shear >= start_dr) return -1;

        // MATLAB: [~,dicrotic_notch_idx] = min(beat(min_shear:start_diastolic_relax))
        int dn = min_shear;
        for (int i = min_shear; i <= start_dr && i < (int)beat.size(); i++)
            if (beat[i] < beat[dn]) dn = i;

        // MATLAB: dicrotic_notch_idx = dicrotic_notch_idx + min_shear-1
        // Already computed as absolute index.

        return dn;
    }

    // ─── PPG Signal Quality Index (exact PPG_SQI.m) ────────────────────────────
    // Uses pre-allocated per-thread SqiWorkspace passed in from caller.
    // OpenMP parallelization on inner beat loop (matching MATLAB parfor).
    inline SqiResult PPG_SQI(const std::vector<double>& wave,
        const std::vector<int>& anntime,
        const std::vector<double>& tmplate,
        int /*windowlen*/, double Fs,
        std::vector<std::unique_ptr<SqiWorkspace>>& workspaces) {
        SqiResult result;
        result.labels = { "mean_corr_dtw", "corrcoff_direct", "corrcoff_interp", "dtw", "frechet" };
        int num_beats = (int)anntime.size() - 1;
        result.typeMnemonic.resize(num_beats, std::vector<double>(5, 0.0));

        if (tmplate.empty()) return result;

        // MATLAB: wave = PPGmedianfilter(wave, Fs, Fs)
        std::vector<double> filtered = PPGmedianfilter(wave, (int)Fs, Fs);

        const auto& t = tmplate;
        int tlen = (int)t.size();

        // MATLAB: d1 = (t - min(t)) / (max(t) - min(t)) * 100
        double t_min = *std::min_element(t.begin(), t.end());
        double t_max = *std::max_element(t.begin(), t.end());
        double t_range = t_max - t_min;
        std::vector<double> d1(tlen);
        for (int i = 0; i < tlen; i++)
            d1[i] = (t_range > 0) ? (t[i] - t_min) / t_range * 100.0 : 0.0;

        // MATLAB: [tmp, pla1] = PLA(d1, 1, 1)
        auto [y1, pla1] = PLA(d1, 1, 1.0);
        int n1_pla = (int)pla1.size();
        int n1_sig = (int)y1.size();

        // Pre-compute y1 normalized for Frechet (shared across all beats)
        std::vector<double> y1n_shared(n1_sig);
        for (int k = 0; k < n1_sig; k++) y1n_shared[k] = y1[k] / 100.0;

        // Parallel loop — workspaces already allocated by caller
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto& ws = *workspaces[tid];

#pragma omp for schedule(dynamic, 16)
            for (int j = 0; j < num_beats; j++) {
                try {
                    int bb = anntime[j], be = anntime[j + 1];
                    if (be - bb > 3 * (int)Fs) be = bb + 3 * (int)Fs;
                    int complength = std::min(tlen, be - bb - 1);

                    if (bb + complength - 1 >= (int)filtered.size() ||
                        be > (int)filtered.size() || bb < 0 || complength < 1) {
                        result.typeMnemonic[j] = { 0,0,0,0,0 };
                        continue;
                    }

                    int blen = be - bb;
                    if (blen > SqiWorkspace::MAX_BEAT || blen < 1) {
                        result.typeMnemonic[j] = { 0,0,0,0,0 };
                        continue;
                    }

                    // ── SQI1: Direct compare ──
                    double c1 = corrcoef(t.data(), filtered.data() + bb, complength);
                    if (c1 < 0) c1 = 0;

                    // ── SQI2: Linear resampling ──
                    for (int k = 0; k < blen; k++) {
                        ws.bx[k] = (double)(k + 1);
                        ws.by[k] = filtered[bb + k];
                    }
                    for (int k = 0; k < tlen; k++)
                        ws.ix_buf[k] = 1.0 + (double)k * (blen - 1.0) / (tlen - 1.0);
                    interp1_linear_ws(ws.bx, ws.by, blen, ws.ix_buf, ws.yi, tlen);
                    for (int k = 0; k < tlen; k++)
                        if (std::isnan(ws.yi[k])) ws.yi[k] = 0;
                    double c2 = corrcoef(t.data(), ws.yi, tlen);
                    if (c2 < 0) c2 = 0;

                    // ── SQI3: DTW ──
                    double c3 = 0;
                    // Copy beat into workspace
                    for (int k = 0; k < blen; k++) ws.d2[k] = filtered[bb + k];

                    if (blen <= (int)d1.size() * 10) {
                        double d2mn = ws.d2[0], d2mx = ws.d2[0];
                        for (int k = 1; k < blen; k++) {
                            if (ws.d2[k] < d2mn) d2mn = ws.d2[k];
                            if (ws.d2[k] > d2mx) d2mx = ws.d2[k];
                        }
                        double d2r = d2mx - d2mn;
                        for (int k = 0; k < blen; k++)
                            ws.d2n[k] = (d2r > 0) ? (ws.d2[k] - d2mn) / d2r * 100.0 : 0.0;

                        // PLA on beat
                        int n2_pla = PLA_ws(ws.d2n, blen, 1, 1.0, ws.pla2, SqiWorkspace::MAX_PLA);

                        // Check DTW matrix fits
                        if ((size_t)n1_pla * n2_pla <= SqiWorkspace::MAX_DTW) {
                            simmx_dtw_ws(y1.data(), pla1.data(), n1_pla,
                                ws.d2n, ws.pla2, n2_pla, ws.w);
                            int path_len = dp_dtw2_ws(ws.w, n1_pla, n2_pla,
                                ws.D, ws.tbi, ws.tbj,
                                ws.path_p, ws.path_q);
                            draw_dtw_ws(y1.data(), pla1.data(), n1_pla,
                                ws.path_p, ws.d2n, ws.pla2, n2_pla,
                                ws.path_q, path_len,
                                ws.ym2, n1_sig);
                            c3 = corrcoef(y1.data(), ws.ym2, n1_sig);
                            if (c3 < 0) c3 = 0;
                        }
                    }

                    // ── SQI4: Frechet ──
                    // Normalize d2 for Frechet
                    {
                        double mn2 = ws.d2[0], mx2 = ws.d2[0];
                        for (int k = 1; k < blen; k++) {
                            if (ws.d2[k] < mn2) mn2 = ws.d2[k];
                            if (ws.d2[k] > mx2) mx2 = ws.d2[k];
                        }
                        double r2 = mx2 - mn2;
                        for (int k = 0; k < blen; k++)
                            ws.y2n[k] = (r2 > 0) ? (ws.d2[k] - mn2) / r2 : 0.0;
                    }

                    double cm;
                    int flen_a, flen_b;
                    const double* fptr_a;
                    const double* fptr_b;

                    if (tlen < blen) {
                        // Interp y1n to blen
                        for (int k = 0; k < tlen; k++) ws.xo[k] = (double)k / (tlen - 1);
                        for (int k = 0; k < blen; k++) ws.xf[k] = (double)k / (blen - 1);
                        interp1_linear_ws(ws.xo, y1n_shared.data(), tlen, ws.xf, ws.interp_out, blen);
                        fptr_a = ws.interp_out; flen_a = blen;
                        fptr_b = ws.y2n;        flen_b = blen;
                    }
                    else if (tlen > blen) {
                        // Interp y2n to tlen
                        for (int k = 0; k < blen; k++) ws.xo[k] = (double)k / (blen - 1);
                        for (int k = 0; k < tlen; k++) ws.xf[k] = (double)k / (tlen - 1);
                        interp1_linear_ws(ws.xo, ws.y2n, blen, ws.xf, ws.interp_out, tlen);
                        fptr_a = y1n_shared.data(); flen_a = tlen;
                        fptr_b = ws.interp_out;     flen_b = tlen;
                    }
                    else {
                        fptr_a = y1n_shared.data(); flen_a = tlen;
                        fptr_b = ws.y2n;            flen_b = blen;
                    }

                    double* ca_buf = ws.get_frechet_buf(flen_a, flen_b);
                    cm = DiscreteFrechetDist_ws(fptr_a, flen_a, fptr_b, flen_b, ca_buf);

                    double Frechet = (std::sqrt(2.0) - cm / 2.0) / std::sqrt(2.0);
                    result.typeMnemonic[j] = { (c1 + c2 + c3) / 3.0, c1, c2, c3, Frechet };
                }
                catch (...) {
                    result.typeMnemonic[j] = { 0,0,0,0,0 };
                }
            }
        } // end parallel

        return result;
    }

    // ─── PPG_SQI_serial — single-thread version using caller's workspace ────────
    // Called when the OUTER segment loop is parallelized (one segment per thread).
    inline SqiResult PPG_SQI_serial(const std::vector<double>& wave,
        const std::vector<int>& anntime,
        const std::vector<double>& tmplate,
        int /*windowlen*/, double Fs,
        SqiWorkspace& ws) {
        SqiResult result;
        result.labels = { "mean_corr_dtw", "corrcoff_direct", "corrcoff_interp", "dtw", "frechet" };
        int num_beats = (int)anntime.size() - 1;
        result.typeMnemonic.resize(num_beats, std::vector<double>(5, 0.0));
        if (tmplate.empty()) return result;

        std::vector<double> filtered = PPGmedianfilter(wave, (int)Fs, Fs);
        const auto& t = tmplate;
        int tlen = (int)t.size();

        double t_min = *std::min_element(t.begin(), t.end());
        double t_max = *std::max_element(t.begin(), t.end());
        double t_range = t_max - t_min;
        std::vector<double> d1(tlen);
        for (int i = 0; i < tlen; i++)
            d1[i] = (t_range > 0) ? (t[i] - t_min) / t_range * 100.0 : 0.0;

        auto [y1, pla1] = PLA(d1, 1, 1.0);
        int n1_pla = (int)pla1.size();
        int n1_sig = (int)y1.size();

        std::vector<double> y1n_shared(n1_sig);
        for (int k = 0; k < n1_sig; k++) y1n_shared[k] = y1[k] / 100.0;

        for (int j = 0; j < num_beats; j++) {
            try {
                int bb = anntime[j], be = anntime[j + 1];
                if (be - bb > 3 * (int)Fs) be = bb + 3 * (int)Fs;
                int complength = std::min(tlen, be - bb - 1);
                if (bb + complength - 1 >= (int)filtered.size() ||
                    be > (int)filtered.size() || bb < 0 || complength < 1) {
                    result.typeMnemonic[j] = { 0,0,0,0,0 }; continue;
                }
                int blen = be - bb;
                if (blen > SqiWorkspace::MAX_BEAT || blen < 1) {
                    result.typeMnemonic[j] = { 0,0,0,0,0 }; continue;
                }

                double c1 = corrcoef(t.data(), filtered.data() + bb, complength);
                if (c1 < 0) c1 = 0;

                for (int k = 0; k < blen; k++) { ws.bx[k] = (double)(k + 1); ws.by[k] = filtered[bb + k]; }
                for (int k = 0; k < tlen; k++) ws.ix_buf[k] = 1.0 + (double)k * (blen - 1.0) / (tlen - 1.0);
                interp1_linear_ws(ws.bx, ws.by, blen, ws.ix_buf, ws.yi, tlen);
                for (int k = 0; k < tlen; k++) if (std::isnan(ws.yi[k])) ws.yi[k] = 0;
                double c2 = corrcoef(t.data(), ws.yi, tlen);
                if (c2 < 0) c2 = 0;

                double c3 = 0;
                for (int k = 0; k < blen; k++) ws.d2[k] = filtered[bb + k];
                if (blen <= (int)d1.size() * 10) {
                    double d2mn = ws.d2[0], d2mx = ws.d2[0];
                    for (int k = 1; k < blen; k++) { if (ws.d2[k] < d2mn) d2mn = ws.d2[k]; if (ws.d2[k] > d2mx) d2mx = ws.d2[k]; }
                    double d2r = d2mx - d2mn;
                    for (int k = 0; k < blen; k++) ws.d2n[k] = (d2r > 0) ? (ws.d2[k] - d2mn) / d2r * 100.0 : 0.0;
                    int n2_pla = PLA_ws(ws.d2n, blen, 1, 1.0, ws.pla2, SqiWorkspace::MAX_PLA);
                    if ((size_t)n1_pla * n2_pla <= SqiWorkspace::MAX_DTW) {
                        simmx_dtw_ws(y1.data(), pla1.data(), n1_pla, ws.d2n, ws.pla2, n2_pla, ws.w);
                        int path_len = dp_dtw2_ws(ws.w, n1_pla, n2_pla, ws.D, ws.tbi, ws.tbj, ws.path_p, ws.path_q);
                        draw_dtw_ws(y1.data(), pla1.data(), n1_pla, ws.path_p, ws.d2n, ws.pla2, n2_pla, ws.path_q, path_len, ws.ym2, n1_sig);
                        c3 = corrcoef(y1.data(), ws.ym2, n1_sig);
                        if (c3 < 0) c3 = 0;
                    }
                }

                {
                    double mn2 = ws.d2[0], mx2 = ws.d2[0];
                    for (int k = 1; k < blen; k++) { if (ws.d2[k] < mn2)mn2 = ws.d2[k]; if (ws.d2[k] > mx2)mx2 = ws.d2[k]; }
                    double r2 = mx2 - mn2;
                    for (int k = 0; k < blen; k++) ws.y2n[k] = (r2 > 0) ? (ws.d2[k] - mn2) / r2 : 0.0;
                }

                double cm;
                int flen_a, flen_b;
                const double* fptr_a, * fptr_b;
                if (tlen < blen) {
                    for (int k = 0; k < tlen; k++) ws.xo[k] = (double)k / (tlen - 1);
                    for (int k = 0; k < blen; k++) ws.xf[k] = (double)k / (blen - 1);
                    interp1_linear_ws(ws.xo, y1n_shared.data(), tlen, ws.xf, ws.interp_out, blen);
                    fptr_a = ws.interp_out; flen_a = blen; fptr_b = ws.y2n; flen_b = blen;
                }
                else if (tlen > blen) {
                    for (int k = 0; k < blen; k++) ws.xo[k] = (double)k / (blen - 1);
                    for (int k = 0; k < tlen; k++) ws.xf[k] = (double)k / (tlen - 1);
                    interp1_linear_ws(ws.xo, ws.y2n, blen, ws.xf, ws.interp_out, tlen);
                    fptr_a = y1n_shared.data(); flen_a = tlen; fptr_b = ws.interp_out; flen_b = tlen;
                }
                else {
                    fptr_a = y1n_shared.data(); flen_a = tlen; fptr_b = ws.y2n; flen_b = blen;
                }
                double* ca_buf = ws.get_frechet_buf(flen_a, flen_b);
                cm = DiscreteFrechetDist_ws(fptr_a, flen_a, fptr_b, flen_b, ca_buf);
                double Frechet = (std::sqrt(2.0) - cm / 2.0) / std::sqrt(2.0);
                result.typeMnemonic[j] = { (c1 + c2 + c3) / 3.0, c1, c2, c3, Frechet };
            }
            catch (...) {
                result.typeMnemonic[j] = { 0,0,0,0,0 };
            }
        }
        return result;
    }

} // namespace ppg