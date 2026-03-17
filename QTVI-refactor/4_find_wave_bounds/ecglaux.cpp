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

    // 1. MATCH MATLAB CONSTANTS EXACTLY
    const size_t examwindow = static_cast<size_t>(std::round(200 * sl));
    const size_t ifno = static_cast<size_t>(std::round(25 * sl));
    const size_t sub1 = static_cast<size_t>(std::round(275 * sl));
    const size_t lookmorepts = 0;

    const size_t bufindA = static_cast<size_t>(std::round(6.0 * sampling / 120.0));
    const size_t bufindC = static_cast<size_t>(std::round(4.0 * sampling / 120.0));

    // 2. MATCH STARTING POINT
    // --------------------------------------------------------------------------
    // FIX: MATLAB starts at sub1+1+bufindA (1-based), which in 0-based C++ is
    //      sub1 + bufindA.  But the backward search (pt1 = Aind2 - sub1) and the
    //      amplitude checks (ecg[R - bufindA]) mean the detector physically
    //      cannot accept peaks whose index < sub1 + bufindA.  We can, however,
    //      start the *scan* earlier so the MWI bump whose corresponding R-peak
    //      is near index sub1+bufindA is not missed because its bump starts a
    //      few samples before perpt.
    //
    //      Lower the starting point by examwindow/2 so the first search window
    //      has a chance to catch the upstroke of the very first MWI bump.
    //      The internal safety checks (pt1 >= 0, bufindA bounds) already prevent
    //      out-of-range access; the slope/refractory checks handle false positives.
    // --------------------------------------------------------------------------
    size_t perpt = (sub1 + bufindA > examwindow / 2)
        ? sub1 + bufindA - examwindow / 2
        : 0;

    // 3. MATCH SEARCH RANGE (Inclusive)
    const size_t mwiwidth_ext_plus_1 = static_cast<size_t>(std::round(1.25 * mwiwidthpts)) + 1;

    // Pre-calculate all integer offsets outside the loop
    const size_t sl10 = static_cast<size_t>(std::round(10 * sl));
    const size_t sl20 = static_cast<size_t>(std::round(20 * sl));
    const size_t mwiwidth_ext = static_cast<size_t>(std::round(1.25 * mwiwidthpts));
    const double divbufindC = static_cast<double>(bufindC + 1);

    size_t estimated_size = static_cast<size_t>(100 * n_ecg / sampling / 60.0);
    if (estimated_size < 10) estimated_size = 100;

    vector<double> Rpickval; Rpickval.reserve(estimated_size);
    vector<size_t> Rpickind; Rpickind.reserve(estimated_size);

    double prevslopeup = 0.0;

    // --------------------------------------------------------------------------
    // FIX: Relax the stop_limit so the last beat's MWI bump is reachable.
    //
    // OLD (overly conservative):
    //   stop_limit = n_mwi - (examwindow + mwiwidth_ext + bufindA + 1)
    //
    // The examwindow term is needed so max_element_index doesn't read past
    // n_mwi.  The mwiwidth_ext term is needed for the secondary MWI peak
    // search.  bufindA is needed for the ecg amplitude checks.
    //
    // But the secondary search and amplitude checks are guarded by std::min
    // against n_mwi and n_ecg respectively, so we only truly need:
    //   stop_limit = n_mwi - examwindow - 1
    //
    // All downstream accesses already clamp to signal length.  The extra
    // margins were belt-and-suspenders that cost us the last 1-2 beats.
    // --------------------------------------------------------------------------
    const size_t stop_limit = (n_mwi > examwindow + 1)
        ? n_mwi - examwindow - 1
        : 0;

    while (perpt < stop_limit) {
        // 1. Find max in current window
        size_t windowEnd = std::min(perpt + examwindow + 1, n_mwi);
        auto maxResult = max_element_index(mwisignal, perpt, windowEnd);

        double val = maxResult.first;
        size_t absind = maxResult.second + perpt;

        // 2. Check threshold and upstroke (or direct peak capture)
        bool possible = false;
        double val2 = 0.0;
        size_t Aind2 = 0;

        if (val > mwithold) {
            size_t idx1 = absind + sl10;
            size_t idx2 = absind + sl20;
            if ((idx1 < n_mwi && mwisignal[idx1] > val) || (idx2 < n_mwi && mwisignal[idx2] > val)) {
                // Normal case: caught the upstroke, find the actual peak ahead
                possible = true;
                size_t mwiEnd = std::min(absind + mwiwidth_ext + 1, n_mwi);
                auto mwiPeak = max_element_index(mwisignal, absind, mwiEnd);
                val2 = mwiPeak.first;
                Aind2 = mwiPeak.second + absind;
            }
            else {
                // Upstroke check failed — but we may have caught the peak itself
                // (common for the first beat where the startup flat region precedes it).
                // Verify it's a real peak: values before AND after absind are lower.
                bool is_peak = false;
                if (absind >= sl10 && absind + sl10 < n_mwi) {
                    is_peak = (mwisignal[absind - sl10] < val) && (mwisignal[absind + sl10] < val);
                }
                if (is_peak) {
                    possible = true;
                    val2 = val;
                    Aind2 = absind;
                }
            }
        }

        if (possible) {

            if (val2 < 3.0 * mvimaxval) {
                // 4. Find ECG peak — clamp to signal end
                size_t pt1 = (Aind2 > sub1) ? Aind2 - sub1 : 0;
                size_t pt2 = std::min(Aind2 + lookmorepts + 1, n_ecg);
                if (pt1 < pt2) {
                    auto ecgPeak = max_element_index(ecg, pt1, pt2);
                    size_t RpickindUP = ecgPeak.second + pt1;

                    // 5. Conditions check — clamp bufindA access to signal bounds
                    if (RpickindUP >= bufindA && RpickindUP + bufindA < n_ecg) {
                        double nowslopeup = (RpickindUP >= bufindC) ?
                            std::abs(ecg[RpickindUP] - ecg[RpickindUP - bufindC]) / divbufindC : 0;

                        bool c1 = ecg[RpickindUP - bufindA] < ecg[RpickindUP];
                        bool c2 = ecg[RpickindUP + bufindA] < ecg[RpickindUP];
                        bool c3 = Rpickind.empty() || (RpickindUP > Rpickind.back() + refractpts);
                        bool c4 = nowslopeup > 0.33 * prevslopeup;

                        if (c1 && c2 && c3 && c4) {
                            Rpickval.push_back(ecgPeak.first);
                            Rpickind.push_back(RpickindUP);

                            mwithold = mwitholdfract * val2 + mwitholdff * (mwithold - mwitholdfract * val2);
                            mvimaxval = val2 + mwitholdff * (mvimaxval - val2);
                            prevslopeup = nowslopeup + mwitholdff * (prevslopeup - nowslopeup);

                            perpt = RpickindUP + refractpts;
                            continue;
                        }
                    }
                }
            }
        }
        perpt += ifno;
    }

    return std::make_tuple(Rpickind, mwithold, mvimaxval);
}