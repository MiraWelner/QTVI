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
    double sl = sampling / 1000.0;

    // now i will examine the mwisignal for threshold detection
    int examwindow = static_cast<int>(std::round(200 * sl));  // should be <50% of a typical RRint
    int sub1 = static_cast<int>(std::round(275 * sl));  // subtract val 1
    int sub2 = static_cast<int>(std::round(150 * sl));  // subtract val 2
    int lookmorepts = static_cast<int>(std::round(0 * sl));
    int ifno = static_cast<int>(std::round(25 * sl));  // 25 msec ahead if not picked

    int bufindA = static_cast<int>(std::round(6 * sampling / 120.0));
    int bufindB = static_cast<int>(std::round(2 * sampling / 120.0));
    int bufindC = static_cast<int>(std::round(4 * sampling / 120.0));
    int divbufindC = bufindC + 1;

    // pre-make with zeros
    size_t estimated_size = static_cast<size_t>(std::round(100 * ecg.size() / sampling / 60.0));
    vector<double> Rpickval(estimated_size, 0.0);
    vector<size_t> Rpickind(estimated_size, 0);
    Rpickind[0] = static_cast<size_t>(-9999 * refractpts);
    double prevslopeup = 0.0;

    size_t limitofwhile = mwisignal.size() - examwindow - static_cast<size_t>(1.25 * mwiwidthpts) - 1 - lookmorepts - bufindA;
    double mwiwidthpts1pt25 = 1.25 * mwiwidthpts;
    double sl10 = 10 * sl;
    double sl20 = 20 * sl;

    size_t perpt = sub1 + 1 + bufindA;
    size_t Rcount = 0;

    while (perpt < limitofwhile) {
        char taken = 'N';

        // Find max in window
        auto maxResult = max_element_index(mwisignal, perpt, perpt + examwindow);
        double val = maxResult.first;
        size_t ind = maxResult.second;
        size_t absind = ind + perpt;

        // i will catch the upstroke of the MWI, at the extreme right of this window
        size_t idx1 = static_cast<size_t>(std::round(absind + sl10));
        size_t idx2 = static_cast<size_t>(std::round(absind + sl20));

        if (val > mwithold && (idx1 < mwisignal.size() && mwisignal[idx1] > val) ||
            (idx2 < mwisignal.size() && mwisignal[idx2] > val)) {
            // i caught the upstroke in the right window, find the max of this bump
            size_t endIdx = static_cast<size_t>(std::round(absind + mwiwidthpts1pt25));
            if (endIdx > mwisignal.size()) endIdx = mwisignal.size();

            auto maxResult2 = max_element_index(mwisignal, absind, endIdx);
            double val2 = maxResult2.first;
            size_t ind2 = maxResult2.second;
            size_t Aind2 = ind2 + absind;

            // this is the max of the MWI bump
            if (val2 < 3.0 * mvimaxval) {
                size_t pt1 = Aind2 - sub1;
                size_t pt2 = Aind2 + lookmorepts;

                if (pt2 >= ecg.size()) pt2 = ecg.size();

                auto maxECG = max_element_index(ecg, pt1, pt2);
                double RpickvalUP = maxECG.first;
                size_t RpickindUP = maxECG.second + pt1;

                double nowslopeup = 0.0;
                if (RpickindUP >= bufindC) {
                    nowslopeup = std::abs(ecg[RpickindUP] - ecg[RpickindUP - bufindC]) / divbufindC;
                }

                // need to decide whether to keep slope crit or not
                bool condition1 = RpickindUP >= bufindA && ecg[RpickindUP - bufindA] < ecg[RpickindUP];
                bool condition2 = RpickindUP + bufindA < ecg.size() && ecg[RpickindUP + bufindA] < ecg[RpickindUP];
                bool condition3 = RpickindUP > Rpickind[std::max(static_cast<int>(Rcount), 1) - 1] + refractpts;
                bool condition4 = nowslopeup > 0.33 * prevslopeup;

                if (condition1 && condition2 && condition3 && condition4) {
                    taken = 'Y';
                }
            }
        }

        if (taken == 'Y') {
            size_t pt1 = static_cast<size_t>(std::max(0, static_cast<int>(absind) - sub1));
            size_t pt2 = absind + lookmorepts;
            if (pt2 >= ecg.size()) pt2 = ecg.size();

            auto maxECG = max_element_index(ecg, pt1, pt2);
            double RpickvalUP = maxECG.first;
            size_t RpickindUP = maxECG.second + pt1;

            auto maxResult2 = max_element_index(mwisignal, absind,
                static_cast<size_t>(std::round(absind + mwiwidthpts1pt25)));
            double val2 = maxResult2.first;

            double nowslopeup = 0.0;
            if (RpickindUP >= bufindC) {
                nowslopeup = std::abs(ecg[RpickindUP] - ecg[RpickindUP - bufindC]) / divbufindC;
            }

            if (Rcount < Rpickind.size()) {
                Rpickval[Rcount] = RpickvalUP;
                Rpickind[Rcount] = RpickindUP;
                Rcount++;
            }

            perpt = RpickindUP + refractpts;
            mwithold = mwitholdfract * val2 + mwitholdff * (mwithold - mwitholdfract * val2);
            mvimaxval = val2 + mwitholdff * (mvimaxval - val2);
            prevslopeup = nowslopeup + mwitholdff * (prevslopeup - nowslopeup);
        }
        else {
            perpt = perpt + ifno;
        }
    }

    // Resize to actual count
    vector<size_t> rwave(Rpickind.begin(), Rpickind.begin() + Rcount);

    return std::make_tuple(rwave, mwithold, mvimaxval);
}
