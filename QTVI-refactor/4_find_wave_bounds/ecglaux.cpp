// ============================================================================
// File: ecglaux.cpp
// ============================================================================
#include "ecglaux.h"
#include "StatsUtils.h"

tuple<vector<size_t>, double, double> ecglaux(
    const vector<double>& ecg,
    const vector<double>& mwisignal,
    int sampling,
    double mwithold,
    double mvimaxval,
    int mwiwidthpts,
    int refractpts,
    double mwitholdfract,
    double mwitholdff
) {
    const double sl = sampling / 1000.0;
    const size_t n_mwi = mwisignal.size();
    const size_t n_ecg = ecg.size();

    // Pre-calculate all integer offsets outside the loop
    const size_t examwindow = static_cast<size_t>(std::round(200 * sl));
    const size_t sub1 = static_cast<size_t>(std::round(275 * sl));
    const size_t lookmorepts = static_cast<size_t>(std::round(0 * sl));
    const size_t ifno = static_cast<size_t>(std::round(25 * sl));
    const size_t bufindA = static_cast<size_t>(std::round(6 * sampling / 120.0));
    const size_t bufindC = static_cast<size_t>(std::round(4 * sampling / 120.0));
    const size_t sl10 = static_cast<size_t>(std::round(10 * sl));
    const size_t sl20 = static_cast<size_t>(std::round(20 * sl));
    const size_t mwiwidth_ext = static_cast<size_t>(std::round(1.25 * mwiwidthpts));
    const double divbufindC = static_cast<double>(bufindC + 1);

    size_t estimated_size = static_cast<size_t>(100 * n_ecg / sampling / 60.0);
    if (estimated_size < 10) estimated_size = 100; // sensible minimum

    vector<double> Rpickval; Rpickval.reserve(estimated_size);
    vector<size_t> Rpickind; Rpickind.reserve(estimated_size);

    double prevslopeup = 0.0;
    size_t perpt = sub1 + 1 + bufindA;
    // Safety boundary
    const size_t stop_limit = (n_mwi > (examwindow + mwiwidth_ext + bufindA + 1)) ?
        n_mwi - (examwindow + mwiwidth_ext + bufindA + 1) : 0;

    while (perpt < stop_limit) {
        // 1. Find max in current window
        size_t windowEnd = perpt + examwindow;
        auto maxResult = max_element_index(mwisignal, perpt, windowEnd);
        double val = maxResult.first;
        size_t absind = maxResult.second + perpt;

        // 2. Check threshold and upstroke (fast check)
        bool possible = false;
        if (val > mwithold) {
            size_t idx1 = absind + sl10;
            size_t idx2 = absind + sl20;
            if ((idx1 < n_mwi && mwisignal[idx1] > val) || (idx2 < n_mwi && mwisignal[idx2] > val)) {
                possible = true;
            }
        }

        if (possible) {
            // 3. Find MWI peak (Only call this ONCE)
            size_t mwiEnd = std::min(absind + mwiwidth_ext, n_mwi);
            auto mwiPeak = max_element_index(mwisignal, absind, mwiEnd);
            double val2 = mwiPeak.first;
            size_t Aind2 = mwiPeak.second + absind;

            if (val2 < 3.0 * mvimaxval) {
                // 4. Find ECG peak (Only call this ONCE)
                size_t pt1 = (Aind2 > sub1) ? Aind2 - sub1 : 0;
                size_t pt2 = std::min(Aind2 + lookmorepts, n_ecg);

                if (pt1 < pt2) {
                    auto ecgPeak = max_element_index(ecg, pt1, pt2);
                    size_t RpickindUP = ecgPeak.second + pt1;

                    // 5. Conditions check
                    if (RpickindUP >= bufindA && RpickindUP + bufindA < n_ecg) {
                        double nowslopeup = (RpickindUP >= bufindC) ?
                            std::abs(ecg[RpickindUP] - ecg[RpickindUP - bufindC]) / divbufindC : 0;

                        bool c1 = ecg[RpickindUP - bufindA] < ecg[RpickindUP];
                        bool c2 = ecg[RpickindUP + bufindA] < ecg[RpickindUP];
                        bool c3 = Rpickind.empty() || (RpickindUP > Rpickind.back() + refractpts);
                        bool c4 = nowslopeup > 0.33 * prevslopeup;

                        if (c1 && c2 && c3 && c4) {
                            // ACCEPTED - Use already calculated values
                            Rpickval.push_back(ecgPeak.first);
                            Rpickind.push_back(RpickindUP);

                            // Update adaptive thresholds
                            mwithold = mwitholdfract * val2 + mwitholdff * (mwithold - mwitholdfract * val2);
                            mvimaxval = val2 + mwitholdff * (mvimaxval - val2);
                            prevslopeup = nowslopeup + mwitholdff * (prevslopeup - nowslopeup);

                            perpt = RpickindUP + refractpts;
                            continue; // Move to next peak
                        }
                    }
                }
            }
        }
        perpt += ifno;
    }

    return std::make_tuple(Rpickind, mwithold, mvimaxval);
}
