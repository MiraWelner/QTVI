// ============================================================================
// File: pan_tompkin.cpp   `  
// Complete implementation of Pan-Tompkins algorithm
// References:
// [1] Sedghamiz. H, "Matlab Implementation of Pan Tompkins ECG QRS detector.", 2014
// [2] PAN.J, TOMPKINS. W.J, "A Real-Time QRS Detection Algorithm" IEEE
//     TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. BME-32, NO. 3, MARCH 1985
// ============================================================================
#include "pan_tompkin.h"
#include "FilterUtils.h"
#include "StatsUtils.h"
#include "PeakFinder.h"

PanTompkinResult pan_tompkin(const vector<double>& ecg, double fs, int gr) {
    PanTompkinResult result;
    result.delay = 0;

    if (ecg.empty()) {
        return result;
    }

    int skip = 0;
    double m_selected_RR = 0;
    double mean_RR = 0;
    int ser_back = 0;

    // ============ Noise cancelation (Filtering) (5-15 Hz) ===============
    vector<double> ecg_h;

    if (fs == 200) {
        // Remove the mean of signal
        double meanVal = mean(ecg);
        vector<double> ecg_centered = ecg;
        for (auto& val : ecg_centered) {
            val -= meanVal;
        }

        // Low Pass Filter
        double Wn = 12.0 * 2.0 / fs;
        int N = 3;
        vector<double> a, b;
        butter(N, Wn, "low", b, a);
        vector<double> ecg_l = filtfilt(a, b, ecg_centered);

        // Normalize
        double maxAbs = 0.0;
        for (const auto& val : ecg_l) {
            maxAbs = std::max(maxAbs, std::abs(val));
        }
        if (maxAbs > 0) {
            for (auto& val : ecg_l) {
                val /= maxAbs;
            }
        }

        // High Pass filter
        Wn = 5.0 * 2.0 / fs;
        butter(N, Wn, "high", b, a);
        ecg_h = filtfilt(a, b, ecg_l);

        // Normalize
        maxAbs = 0.0;
        for (const auto& val : ecg_h) {
            maxAbs = std::max(maxAbs, std::abs(val));
        }
        if (maxAbs > 0) {
            for (auto& val : ecg_h) {
                val /= maxAbs;
            }
        }
    }
    else {
        // Bandpass filter for other sampling frequencies
        double f1 = 5.0;
        double f2 = 15.0;
        vector<double> Wn = { f1 * 2.0 / fs, f2 * 2.0 / fs };
        int N = 3;
        vector<double> a, b;
        butter(N, Wn, b, a);
        ecg_h = filtfilt(a, b, ecg);

        // Normalize
        double maxAbs = 0.0;
        for (const auto& val : ecg_h) {
            maxAbs = std::max(maxAbs, std::abs(val));
        }
        if (maxAbs > 0) {
            for (auto& val : ecg_h) {
                val /= maxAbs;
            }
        }
    }

    // ==================== derivative filter ==========================
    vector<double> b;
    if (fs != 200) {
        double int_c = (5.0 - 1.0) / (fs * 1.0 / 40.0);
        vector<double> template_b = { 1, 2, 0, -2, -1 };
        for (auto& val : template_b) {
            val *= (1.0 / 8.0) * fs;
        }

        // Interpolate
        for (double i = 1.0; i <= 5.0; i += int_c) {
            size_t idx = static_cast<size_t>(i) - 1;
            if (idx < template_b.size()) {
                b.push_back(template_b[idx]);
            }
        }
    }
    else {
        b = { 1, 2, 0, -2, -1 };
        for (auto& val : b) {
            val *= (1.0 / 8.0) * fs;
        }
    }

    vector<double> ecg_d = filtfilt(b, { 1.0 }, ecg_h);

    // Normalize
    double maxAbs = 0.0;
    for (const auto& val : ecg_d) {
        maxAbs = std::max(maxAbs, std::abs(val));
    }
    if (maxAbs > 0) {
        for (auto& val : ecg_d) {
            val /= maxAbs;
        }
    }

    // ========== Squaring nonlinearly enhance the dominant peaks ==========
    vector<double> ecg_s(ecg_d.size());
    for (size_t i = 0; i < ecg_d.size(); ++i) {
        ecg_s[i] = ecg_d[i] * ecg_d[i];
    }

    // ============ Moving average ==================
    int window_size = static_cast<int>(std::round(0.150 * fs));
    vector<double> ones(window_size, 1.0 / window_size);
    vector<double> ecg_m = conv(ecg_s, ones);
    result.delay += window_size / 2;

    // ===================== Fiducial Marks ==============================
    double minPeakDistance = std::round(0.2 * fs);
    vector<double> pks;
    vector<size_t> locs;
    findpeaks(ecg_m, pks, locs, minPeakDistance);

    size_t LLp = pks.size();

    // =================== Initialize Some Other Parameters ===============
    vector<double> qrs_c(LLp, 0.0);
    vector<size_t> qrs_i(LLp, 0);
    vector<size_t> qrs_i_raw(LLp, 0);
    vector<double> qrs_amp_raw(LLp, 0.0);

    vector<double> nois_c(LLp, 0.0);
    vector<size_t> nois_i(LLp, 0);

    vector<double> SIGL_buf(LLp, 0.0);
    vector<double> NOISL_buf(LLp, 0.0);
    vector<double> THRS_buf(LLp, 0.0);

    vector<double> SIGL_buf1(LLp, 0.0);
    vector<double> NOISL_buf1(LLp, 0.0);
    vector<double> THRS_buf1(LLp, 0.0);

    // Initialize the training phase (2 seconds of the signal)
    size_t training_size = static_cast<size_t>(2 * fs);
    if (training_size > ecg_m.size()) {
        training_size = ecg_m.size();
    }

    double max_train = -Inf;
    double sum_train = 0.0;
    for (size_t i = 0; i < training_size; ++i) {
        if (ecg_m[i] > max_train) {
            max_train = ecg_m[i];
        }
        sum_train += ecg_m[i];
    }

    double THR_SIG = max_train / 3.0;
    double THR_NOISE = (sum_train / training_size) / 2.0;
    double SIG_LEV = THR_SIG;
    double NOISE_LEV = THR_NOISE;

    // Initialize bandpass filter threshold (2 seconds)
    max_train = -Inf;
    sum_train = 0.0;
    training_size = std::min(training_size, ecg_h.size());
    for (size_t i = 0; i < training_size; ++i) {
        if (ecg_h[i] > max_train) {
            max_train = ecg_h[i];
        }
        sum_train += ecg_h[i];
    }

    double THR_SIG1 = max_train / 3.0;
    double THR_NOISE1 = (sum_train / training_size) / 2.0;
    double SIG_LEV1 = THR_SIG1;
    double NOISE_LEV1 = THR_NOISE1;

    // ============ Thresholding and decision rule =============
    size_t Beat_C = 0;
    size_t Beat_C1 = 0;
    size_t Noise_Count = 0;

    for (size_t i = 0; i < LLp; ++i) {
        // ===== locate the corresponding peak in the filtered signal =====
        double y_i = 0;
        size_t x_i = 0;

        int search_start = static_cast<int>(locs[i]) - static_cast<int>(std::round(0.150 * fs));
        int search_end = static_cast<int>(locs[i]);

        if (search_start >= 0 && locs[i] < ecg_h.size()) {
            auto maxResult = max_element_index(ecg_h, search_start, search_end + 1);
            y_i = maxResult.first;
            x_i = maxResult.second;
        }
        else {
            if (i == 0) {
                auto maxResult = max_element_index(ecg_h, 0, locs[i] + 1);
                y_i = maxResult.first;
                x_i = maxResult.second;
                ser_back = 1;
            }
            else if (locs[i] >= ecg_h.size()) {
                size_t start = locs[i] >= static_cast<size_t>(std::round(0.150 * fs)) ?
                    locs[i] - static_cast<size_t>(std::round(0.150 * fs)) : 0;
                auto maxResult = max_element_index(ecg_h, start, ecg_h.size());
                y_i = maxResult.first;
                x_i = maxResult.second;
            }
        }

        // ================= update the heart_rate ====================
        if (Beat_C >= 9) {
            // Calculate RR interval
            double sum_diff = 0.0;
            for (size_t j = Beat_C - 8; j < Beat_C; ++j) {
                sum_diff += (qrs_i[j + 1] - qrs_i[j]);
            }
            mean_RR = sum_diff / 8.0;

            double comp = qrs_i[Beat_C] - qrs_i[Beat_C - 1];

            if (comp <= 0.92 * mean_RR || comp >= 1.16 * mean_RR) {
                THR_SIG = 0.5 * THR_SIG;
                THR_SIG1 = 0.5 * THR_SIG1;
            }
            else {
                m_selected_RR = mean_RR;
            }
        }

        // == calculate the mean last 8 R waves to ensure that QRS is not ====
        double test_m = 0;
        if (m_selected_RR > 0) {
            test_m = m_selected_RR;
        }
        else if (mean_RR > 0) {
            test_m = mean_RR;
        }

        int not_nois = 0;

        if (test_m > 0) {
            if (Beat_C > 0 && (locs[i] - qrs_i[Beat_C - 1]) >= std::round(1.66 * test_m)) {
                // Search back and locate the max in this interval
                size_t search_start = qrs_i[Beat_C - 1] + static_cast<size_t>(std::round(0.200 * fs));
                size_t search_end = locs[i] - static_cast<size_t>(std::round(0.200 * fs));

                if (search_end > ecg_m.size()) search_end = ecg_m.size();
                if (search_start < search_end) {
                    auto maxResult = max_element_index(ecg_m, search_start, search_end);
                    double pks_temp = maxResult.first;
                    size_t locs_temp = search_start + maxResult.second;

                    if (pks_temp > THR_NOISE) {
                        qrs_c[Beat_C] = pks_temp;
                        qrs_i[Beat_C] = locs_temp;
                        Beat_C++;

                        // Locate in filtered signal
                        size_t filter_start = locs_temp >= static_cast<size_t>(std::round(0.150 * fs)) ?
                            locs_temp - static_cast<size_t>(std::round(0.150 * fs)) : 0;
                        size_t filter_end = std::min(locs_temp + 1, ecg_h.size());

                        auto maxResult2 = max_element_index(ecg_h, filter_start, filter_end);
                        double y_i_t = maxResult2.first;
                        size_t x_i_t = maxResult2.second;

                        if (y_i_t > THR_NOISE1) {
                            qrs_i_raw[Beat_C1] = filter_start + x_i_t;
                            qrs_amp_raw[Beat_C1] = y_i_t;
                            Beat_C1++;
                            SIG_LEV1 = 0.25 * y_i_t + 0.75 * SIG_LEV1;
                        }

                        not_nois = 1;
                        SIG_LEV = 0.25 * pks_temp + 0.75 * SIG_LEV;
                    }
                }
            }
        }

        // ===================  find noise and QRS peaks ==================
        if (pks[i] >= THR_SIG) {
            // if No QRS in 360ms of the previous QRS See if T wave
            if (Beat_C >= 3) {
                if (Beat_C > 0 && (locs[i] - qrs_i[Beat_C - 1]) <= std::round(0.3600 * fs)) {
                    // Calculate slopes
                    double sum1 = 0.0, sum2 = 0.0;
                    int count1 = 0, count2 = 0;

                    size_t slope_start = locs[i] >= static_cast<size_t>(std::round(0.075 * fs)) ?
                        locs[i] - static_cast<size_t>(std::round(0.075 * fs)) : 0;

                    for (size_t j = slope_start; j < locs[i] && j + 1 < ecg_m.size(); ++j) {
                        sum1 += (ecg_m[j + 1] - ecg_m[j]);
                        count1++;
                    }
                    double Slope1 = count1 > 0 ? sum1 / count1 : 0;

                    slope_start = qrs_i[Beat_C - 1] >= static_cast<size_t>(std::round(0.075 * fs)) ?
                        qrs_i[Beat_C - 1] - static_cast<size_t>(std::round(0.075 * fs)) : 0;

                    for (size_t j = slope_start; j < qrs_i[Beat_C - 1] && j + 1 < ecg_m.size(); ++j) {
                        sum2 += (ecg_m[j + 1] - ecg_m[j]);
                        count2++;
                    }
                    double Slope2 = count2 > 0 ? sum2 / count2 : 0;

                    if (std::abs(Slope1) <= std::abs(0.5 * Slope2)) {
                        nois_c[Noise_Count] = pks[i];
                        nois_i[Noise_Count] = locs[i];
                        Noise_Count++;
                        skip = 1;
                        NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1;
                        NOISE_LEV = 0.125 * pks[i] + 0.875 * NOISE_LEV;
                    }
                    else {
                        skip = 0;
                    }
                }
            }

            // skip is 1 when a T wave is detected
            if (skip == 0) {
                qrs_c[Beat_C] = pks[i];
                qrs_i[Beat_C] = locs[i];
                Beat_C++;

                // bandpass filter check threshold
                if (y_i >= THR_SIG1) {
                    if (ser_back) {
                        qrs_i_raw[Beat_C1] = x_i;
                    }
                    else {
                        qrs_i_raw[Beat_C1] = locs[i] >= static_cast<size_t>(std::round(0.150 * fs)) ?
                            locs[i] - static_cast<size_t>(std::round(0.150 * fs)) + x_i : x_i;
                    }
                    qrs_amp_raw[Beat_C1] = y_i;
                    Beat_C1++;
                    SIG_LEV1 = 0.125 * y_i + 0.875 * SIG_LEV1;
                }
                SIG_LEV = 0.125 * pks[i] + 0.875 * SIG_LEV;
            }
        }
        else if (THR_NOISE <= pks[i] && pks[i] < THR_SIG) {
            NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1;
            NOISE_LEV = 0.125 * pks[i] + 0.875 * NOISE_LEV;
        }
        else if (pks[i] < THR_NOISE) {
            nois_c[Noise_Count] = pks[i];
            nois_i[Noise_Count] = locs[i];
            Noise_Count++;
            NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1;
            NOISE_LEV = 0.125 * pks[i] + 0.875 * NOISE_LEV;
        }

        // ================== adjust the threshold with SNR =============
        if (NOISE_LEV != 0 || SIG_LEV != 0) {
            THR_SIG = NOISE_LEV + 0.25 * std::abs(SIG_LEV - NOISE_LEV);
            THR_NOISE = 0.5 * THR_SIG;
        }

        // adjust the threshold with SNR for bandpassed signal
        if (NOISE_LEV1 != 0 || SIG_LEV1 != 0) {
            THR_SIG1 = NOISE_LEV1 + 0.25 * std::abs(SIG_LEV1 - NOISE_LEV1);
            THR_NOISE1 = 0.5 * THR_SIG1;
        }

        // take a track of thresholds
        SIGL_buf[i] = SIG_LEV;
        NOISL_buf[i] = NOISE_LEV;
        THRS_buf[i] = THR_SIG;

        SIGL_buf1[i] = SIG_LEV1;
        NOISL_buf1[i] = NOISE_LEV1;
        THRS_buf1[i] = THR_SIG1;

        // reset parameters
        skip = 0;
        not_nois = 0;
        ser_back = 0;
    }

    // ======================= Adjust Lengths ============================
    result.qrs_i_raw.assign(qrs_i_raw.begin(), qrs_i_raw.begin() + Beat_C1);
    result.qrs_amp_raw.assign(qrs_amp_raw.begin(), qrs_amp_raw.begin() + Beat_C1);

    return result;
}