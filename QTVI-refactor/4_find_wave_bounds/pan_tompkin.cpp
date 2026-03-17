#include "pan_tompkin.h"
#include "FilterUtils.h"
#include "StatsUtils.h"
#include "PeakFinder.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <iomanip>

using namespace std;

// Simple linear interpolation replicating MATLAB's interp1(x, v, xq, 'linear')
static vector<double> interp1_linear(const vector<double>& x, const vector<double>& v,
    const vector<double>& xq) {
    vector<double> result(xq.size());
    for (size_t i = 0; i < xq.size(); ++i) {
        double xi = xq[i];
        if (xi <= x.front()) { result[i] = v.front(); continue; }
        if (xi >= x.back()) { result[i] = v.back();  continue; }
        size_t k = 0;
        for (k = 0; k + 1 < x.size(); ++k) {
            if (x[k + 1] >= xi) break;
        }
        double t = (xi - x[k]) / (x[k + 1] - x[k]);
        result[i] = v[k] + t * (v[k + 1] - v[k]);
    }
    return result;
}

// Full causal convolution replicating MATLAB's conv(a, b) with output length = len(a)+len(b)-1
static vector<double> conv_full(const vector<double>& a, const vector<double>& b) {
    size_t na = a.size(), nb = b.size();
    vector<double> out(na + nb - 1, 0.0);
    for (size_t i = 0; i < na; ++i) {
        for (size_t j = 0; j < nb; ++j) {
            out[i + j] += a[i] * b[j];
        }
    }
    return out;
}

PanTompkinResult pan_tompkin(const vector<double>& ecg_input, double fs, int gr, const std::string& fileID) {
    PanTompkinResult result;
    result.delay = 0;
    if (ecg_input.empty()) return result;

    int N_sig = (int)ecg_input.size();
    vector<double> ecg(ecg_input.begin(), ecg_input.end());

    // ========================================================================
    // 1. Bandpass Filter  (5-15 Hz)
    // ========================================================================
    vector<double> ecg_h;
    int delay = 0;

    if ((int)fs == 200) {
        // Remove mean
        double avg = mean(ecg);
        for (auto& v : ecg) v -= avg;

        // Low-pass: butter(3, 12*2/fs, 'low') + filtfilt
        double Wn_lp = 12.0 * 2.0 / fs;
        vector<double> b_lp, a_lp;
        butter(3, Wn_lp, "low", b_lp, a_lp);
        vector<double> ecg_l = filtfilt(b_lp, a_lp, ecg);
        double mx_l = *max_element(ecg_l.begin(), ecg_l.end(),
            [](double a, double b) { return fabs(a) < fabs(b); });
        if (mx_l != 0) for (auto& v : ecg_l) v /= fabs(mx_l);

        // High-pass: butter(3, 5*2/fs, 'high') + filtfilt
        double Wn_hp = 5.0 * 2.0 / fs;
        vector<double> b_hp, a_hp;
        butter(3, Wn_hp, "high", b_hp, a_hp);
        ecg_h = filtfilt(b_hp, a_hp, ecg_l);
        double mx_h = *max_element(ecg_h.begin(), ecg_h.end(),
            [](double a, double b) { return fabs(a) < fabs(b); });
        if (mx_h != 0) for (auto& v : ecg_h) v /= fabs(mx_h);
    }
    else {
        // Remove mean
        //double avg = mean(ecg);
       // for (auto& v : ecg) v -= avg;

        // Bandpass: butter(3, [5 15]*2/fs) + filtfilt
        vector<double> b_bp, a_bp;

        if (std::abs(fs - 256.0) < 0.5) {
            // Precomputed from MATLAB: [b,a] = butter(3, [5 15]*2/256)
            b_bp = {
                 0.0014670007580121303,
                 0.0,
                -0.004401002274036391,
                 0.0,
                 0.004401002274036391,
                 0.0,
                -0.0014670007580121303
            };
            a_bp = {
                 1.0,
                -5.3856706261642193,
                12.210142717095941,
               -14.917369382919638,
                10.359253074547178,
                -3.8776092362023347,
                 0.61132583250308625
            };
        }
        else if (std::abs(fs - 128.0) < 0.5) {
            // If you ever need 128 Hz, run in MATLAB:
            //   [b,a] = butter(3, [5 15]*2/128);
            //   fprintf('%.17g ', b); fprintf('\n');
            //   fprintf('%.17g ', a); fprintf('\n');
            // and paste the results here.
            //
            // For now, fall back to runtime computation:
            double f1 = 5.0, f2 = 15.0;
            vector<double> Wn = { f1 * 2.0 / fs, f2 * 2.0 / fs };
            butter(3, Wn, b_bp, a_bp);
        }
        else {
            // Fallback for other sample rates
            double f1 = 5.0, f2 = 15.0;
            vector<double> Wn = { f1 * 2.0 / fs, f2 * 2.0 / fs };
            butter(3, Wn, b_bp, a_bp);
        }

        ecg_h = filtfilt(b_bp, a_bp, ecg);
        double mx = *max_element(ecg_h.begin(), ecg_h.end(),
            [](double a, double b) { return fabs(a) < fabs(b); });
        if (mx != 0) for (auto& v : ecg_h) v /= fabs(mx);
    }

    if (ecg_h.empty()) return result;
    int len = (int)ecg_h.size();


    // ========================================================================
    // 2. Derivative Filter
    //    H(z) = (1/8T)(-z^-2 - 2z^-1 + 2z + z^2)
    // ========================================================================
    double scale_d = (1.0 / 8.0) * fs;
    vector<double> base_coeffs = { 1.0 * scale_d, 2.0 * scale_d, 0.0, -2.0 * scale_d, -1.0 * scale_d };
    vector<double> b_d;

    if ((int)fs != 200) {
        double int_c = (5.0 - 1.0) / (fs * 1.0 / 40.0);
        vector<double> x_orig = { 1.0, 2.0, 3.0, 4.0, 5.0 };
        vector<double> xq;
        for (double v = 1.0; v <= 5.0 + 1e-9; v += int_c) xq.push_back(v);
        b_d = interp1_linear(x_orig, base_coeffs, xq);
    }
    else {
        b_d = base_coeffs;
    }

    vector<double> a_one = { 1.0 };
    vector<double> ecg_d = filtfilt(b_d, a_one, ecg_h);
    double mx_d = *max_element(ecg_d.begin(), ecg_d.end());
    if (mx_d != 0) for (auto& v : ecg_d) v /= mx_d;

    // ========================================================================
    // 3. Squaring
    // ========================================================================
    vector<double> ecg_s(ecg_d.size());
    for (size_t i = 0; i < ecg_d.size(); ++i) ecg_s[i] = ecg_d[i] * ecg_d[i];

    // ========================================================================
    // 4. Moving Average (Causal convolution, MATLAB-style)
    //    ecg_m = conv(ecg_s, ones(1,W)/W)  where W = round(0.150*fs)
    // ========================================================================
    int W = (int)round(0.150 * fs);
    if (W < 1) W = 1;
    vector<double> kernel(W, 1.0 / (double)W);
    vector<double> ecg_m = conv_full(ecg_s, kernel);
    delay += W / 2;

    // ========================================================================
    // 5. Find Peaks on ecg_m
    // ========================================================================
    vector<double> pks;
    vector<size_t> locs;
    findpeaks(ecg_m, pks, locs, round(0.2 * fs));
    int LLp = (int)pks.size();
    if (LLp == 0) return result;

    // ========================================================================
    // 6. Initialize thresholds from first 2 seconds
    // ========================================================================
    int skip_edge = W + 2;
    int init_start = std::min(skip_edge, (int)ecg_m.size());
    int init_end = std::min((int)ecg_m.size(), (int)(2.0 * fs));
    if (init_end - init_start < (int)(0.5 * fs)) {
        init_start = 0;
    }

    double max_m = 0, sum_m = 0;
    int count_m = 0;
    for (int k = init_start; k < init_end; ++k) {
        if (ecg_m[k] > max_m) max_m = ecg_m[k];
        sum_m += ecg_m[k];
        count_m++;
    }
    double THR_SIG = max_m * (1.0 / 3.0);
    double THR_NOISE = (count_m > 0) ? (sum_m / count_m) * 0.5 : 0;
    double SIG_LEV = THR_SIG;
    double NOISE_LEV = THR_NOISE;

    int init_start_h = std::min(skip_edge, (int)ecg_h.size());
    int init_end_h = std::min((int)ecg_h.size(), (int)(2.0 * fs));
    if (init_end_h - init_start_h < (int)(0.5 * fs)) {
        init_start_h = 0;
    }

    double max_h = 0, sum_h = 0;
    int count_h = 0;
    for (int k = init_start_h; k < init_end_h; ++k) {
        if (ecg_h[k] > max_h) max_h = ecg_h[k];
        sum_h += ecg_h[k];
        count_h++;
    }
    double THR_SIG1 = max_h * (1.0 / 3.0);
    double THR_NOISE1 = (count_h > 0) ? (sum_h / count_h) * 0.5 : 0;
    double SIG_LEV1 = THR_SIG1;
    double NOISE_LEV1 = THR_NOISE1;

    // ========================================================================
    // 7. Detection buffers
    // ========================================================================
    vector<double> qrs_c(LLp, 0), qrs_i(LLp, 0);
    vector<double> qrs_i_raw_buf(LLp, 0), qrs_amp_raw_buf(LLp, 0);
    vector<double> nois_c(LLp, 0), nois_i(LLp, 0);

    int Beat_C = 0;
    int Beat_C1 = 0;
    int Noise_Count = 0;
    double m_selected_RR = 0, mean_RR = 0;
    int skip = 0, not_nois = 0, ser_back = 0;
    int round_150 = (int)round(0.150 * fs);

    // ========================================================================
    // 8. Main Detection Loop
    // ========================================================================
    for (int i = 0; i < LLp; ++i) {
        // --- Locate corresponding peak in bandpass signal ecg_h ---
        double y_i = 0;
        int x_i = 0;

        if ((int)locs[i] - round_150 >= 0 && (int)locs[i] < (int)ecg_h.size()) {
            int search_start = (int)locs[i] - round_150;
            int search_end = (int)locs[i];
            if (search_start < 0) search_start = 0;
            if (search_end >= (int)ecg_h.size()) search_end = (int)ecg_h.size() - 1;

            y_i = ecg_h[search_start];
            x_i = 0;
            for (int k = search_start; k <= search_end; ++k) {
                if (ecg_h[k] > y_i) {
                    y_i = ecg_h[k];
                    x_i = k - search_start;
                }
            }
        }
        else {
            if (i == 0) {
                int search_end = min((int)locs[0], (int)ecg_h.size() - 1);
                y_i = ecg_h[0]; x_i = 0;
                for (int k = 0; k <= search_end; ++k) {
                    if (ecg_h[k] > y_i) { y_i = ecg_h[k]; x_i = k; }
                }
                ser_back = 1;
            }
            else if ((int)locs[i] >= (int)ecg_h.size()) {
                int search_start = (int)locs[i] - round_150;
                if (search_start < 0) search_start = 0;
                y_i = ecg_h[search_start]; x_i = 0;
                for (int k = search_start; k < (int)ecg_h.size(); ++k) {
                    if (ecg_h[k] > y_i) { y_i = ecg_h[k]; x_i = k - search_start; }
                }
            }
        }

        // --- Update heart rate estimate (after 9 beats) ---
        if (Beat_C >= 9) {
            double diffRR_sum = 0;
            int n_rr = 0;
            for (int d = Beat_C - 8; d < Beat_C; ++d) {
                diffRR_sum += (qrs_i[d] - qrs_i[d - 1]);
                n_rr++;
            }
            mean_RR = (n_rr > 0) ? diffRR_sum / n_rr : 0;
            double comp = qrs_i[Beat_C - 1] - qrs_i[Beat_C - 2];

            if (comp <= 0.92 * mean_RR || comp >= 1.16 * mean_RR) {
                THR_SIG *= 0.5;
                THR_SIG1 *= 0.5;
            }
            else {
                m_selected_RR = mean_RR;
            }
        }

        // --- Searchback: check if a QRS was missed ---
        double test_m = 0;
        if (m_selected_RR)                    test_m = m_selected_RR;
        else if (mean_RR && !m_selected_RR)   test_m = mean_RR;

        if (test_m && Beat_C >= 1) {
            double last_qrs = qrs_i[Beat_C - 1];
            if (((double)locs[i] - last_qrs) >= round(1.66 * test_m)) {
                int sb_start = (int)last_qrs + (int)round(0.200 * fs);
                int sb_end = (int)locs[i] - (int)round(0.200 * fs);
                if (sb_start < 0) sb_start = 0;
                if (sb_end >= (int)ecg_m.size()) sb_end = (int)ecg_m.size() - 1;

                if (sb_start < sb_end) {
                    double pks_temp = ecg_m[sb_start];
                    int locs_temp = sb_start;
                    for (int k = sb_start; k <= sb_end; ++k) {
                        if (ecg_m[k] > pks_temp) {
                            pks_temp = ecg_m[k];
                            locs_temp = k;
                        }
                    }

                    if (pks_temp > THR_NOISE) {
                        if (Beat_C >= (int)qrs_c.size()) { qrs_c.resize(Beat_C + 1); qrs_i.resize(Beat_C + 1); }
                        qrs_c[Beat_C] = pks_temp;
                        qrs_i[Beat_C] = (double)locs_temp;
                        Beat_C++;

                        // Locate in bandpass signal
                        double y_i_t = 0;
                        int x_i_t = 0;
                        int sb_h_start = locs_temp - round_150;
                        int sb_h_end = locs_temp;
                        if (sb_h_start < 0) sb_h_start = 0;
                        if (sb_h_end >= (int)ecg_h.size()) sb_h_end = (int)ecg_h.size() - 1;

                        for (int k = sb_h_start; k <= sb_h_end; ++k) {
                            if (ecg_h[k] > y_i_t) {
                                y_i_t = ecg_h[k];
                                x_i_t = k - sb_h_start;
                            }
                        }

                        if (y_i_t > THR_NOISE1) {
                            if (Beat_C1 >= (int)qrs_i_raw_buf.size()) {
                                qrs_i_raw_buf.resize(Beat_C1 + 1);
                                qrs_amp_raw_buf.resize(Beat_C1 + 1);
                            }
                            qrs_i_raw_buf[Beat_C1] = (double)(sb_h_start + x_i_t);
                            qrs_amp_raw_buf[Beat_C1] = y_i_t;
                            Beat_C1++;
                            SIG_LEV1 = 0.25 * y_i_t + 0.75 * SIG_LEV1;
                        }

                        not_nois = 1;
                        SIG_LEV = 0.25 * pks_temp + 0.75 * SIG_LEV;
                    }
                }
            }
            else {
                not_nois = 0;
            }
        }

        // --- Main classification: Signal, Noise, or T-wave ---
        if (pks[i] >= THR_SIG) {
            // T-wave check: if within 360ms of previous QRS
            if (Beat_C >= 3) {
                double last_qrs_loc = qrs_i[Beat_C - 1];
                if (((double)locs[i] - last_qrs_loc) <= round(0.3600 * fs)) {
                    int slope1_start = (int)locs[i] - (int)round(0.075 * fs);
                    int slope1_end = (int)locs[i];
                    int slope2_start = (int)last_qrs_loc - (int)round(0.075 * fs);
                    int slope2_end = (int)last_qrs_loc;

                    slope1_start = max(0, slope1_start);
                    slope1_end = min((int)ecg_m.size() - 1, slope1_end);
                    slope2_start = max(0, slope2_start);
                    slope2_end = min((int)ecg_m.size() - 1, slope2_end);

                    double s1 = 0, s2 = 0;
                    int n1 = 0, n2 = 0;
                    for (int k = slope1_start; k < slope1_end; ++k) {
                        s1 += (ecg_m[k + 1] - ecg_m[k]); n1++;
                    }
                    for (int k = slope2_start; k < slope2_end; ++k) {
                        s2 += (ecg_m[k + 1] - ecg_m[k]); n2++;
                    }
                    if (n1 > 0) s1 /= n1;
                    if (n2 > 0) s2 /= n2;

                    if (fabs(s1) <= fabs(0.5 * s2)) {
                        Noise_Count++;
                        if (Noise_Count > (int)nois_c.size()) {
                            nois_c.resize(Noise_Count);
                            nois_i.resize(Noise_Count);
                        }
                        nois_c[Noise_Count - 1] = pks[i];
                        nois_i[Noise_Count - 1] = (double)locs[i];
                        skip = 1;
                        NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1;
                        NOISE_LEV = 0.125 * pks[i] + 0.875 * NOISE_LEV;
                    }
                    else {
                        skip = 0;
                    }
                }
            }

            if (skip == 0) {
                if (Beat_C >= (int)qrs_c.size()) { qrs_c.resize(Beat_C + 1); qrs_i.resize(Beat_C + 1); }
                qrs_c[Beat_C] = pks[i];
                qrs_i[Beat_C] = (double)locs[i];
                Beat_C++;

                if (y_i >= THR_SIG1) {
                    if (Beat_C1 >= (int)qrs_i_raw_buf.size()) {
                        qrs_i_raw_buf.resize(Beat_C1 + 1);
                        qrs_amp_raw_buf.resize(Beat_C1 + 1);
                    }
                    if (ser_back) {
                        qrs_i_raw_buf[Beat_C1] = (double)x_i;
                    }
                    else {
                        qrs_i_raw_buf[Beat_C1] = (double)((int)locs[i] - round_150 + x_i);
                    }
                    qrs_amp_raw_buf[Beat_C1] = y_i;
                    Beat_C1++;
                    SIG_LEV1 = 0.125 * y_i + 0.875 * SIG_LEV1;
                }
                SIG_LEV = 0.125 * pks[i] + 0.875 * SIG_LEV;
            }
        }
        else if (pks[i] >= THR_NOISE && pks[i] < THR_SIG) {
            NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1;
            NOISE_LEV = 0.125 * pks[i] + 0.875 * NOISE_LEV;
        }
        else {
            Noise_Count++;
            if (Noise_Count > (int)nois_c.size()) {
                nois_c.resize(Noise_Count);
                nois_i.resize(Noise_Count);
            }
            nois_c[Noise_Count - 1] = pks[i];
            nois_i[Noise_Count - 1] = (double)locs[i];
            NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1;
            NOISE_LEV = 0.125 * pks[i] + 0.875 * NOISE_LEV;
        }

        // --- Update adaptive thresholds ---
        if (NOISE_LEV != 0 || SIG_LEV != 0) {
            THR_SIG = NOISE_LEV + 0.25 * fabs(SIG_LEV - NOISE_LEV);
            THR_NOISE = 0.5 * THR_SIG;
        }
        if (NOISE_LEV1 != 0 || SIG_LEV1 != 0) {
            THR_SIG1 = NOISE_LEV1 + 0.25 * fabs(SIG_LEV1 - NOISE_LEV1);
            THR_NOISE1 = 0.5 * THR_SIG1;
        }

        // Reset per-iteration flags
        skip = 0;
        not_nois = 0;
        ser_back = 0;
    }

    // ========================================================================
    // 9. Collect output
    // ========================================================================
    result.qrs_i_raw.resize(Beat_C1);
    result.qrs_amp_raw.resize(Beat_C1);
    for (int i = 0; i < Beat_C1; ++i) {
        result.qrs_i_raw[i] = (size_t)qrs_i_raw_buf[i];
        result.qrs_amp_raw[i] = qrs_amp_raw_buf[i];
    }
    result.delay = delay;

    return result;
}