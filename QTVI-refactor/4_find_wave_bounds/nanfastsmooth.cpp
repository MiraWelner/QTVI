// ============================================================================
// File: nanfastsmooth.cpp
// ============================================================================
#include "nanfastsmooth.h"

// Helper function for sliding average
vector<double> sa(const vector<double>& Y, int smoothwidth, double tol) {
    if (smoothwidth == 1) {
        return Y;
    }

    // Bound tolerance
    if (tol < 0) tol = 0;
    if (tol > 1) tol = 1;

    int w = smoothwidth;
    int halfw = w / 2;
    size_t L = Y.size();

    vector<double> s(L, 0.0);
    vector<double> np(L, 0.0);

    if (w % 2 == 1) {  // Odd window
        // Initialize sums and counts
        double SumPoints = 0.0;
        double NumPoints = 0.0;

        for (int i = 0; i <= halfw && i < static_cast<int>(L); ++i) {
            if (!std::isnan(Y[i])) {
                SumPoints += Y[i];
                NumPoints += 1.0;
            }
        }

        s[0] = SumPoints;
        np[0] = NumPoints;

        for (size_t k = 1; k < L; ++k) {
            int removeIdx = static_cast<int>(k) - halfw - 1;
            int addIdx = static_cast<int>(k) + halfw;

            if (removeIdx >= 0 && !std::isnan(Y[removeIdx])) {
                SumPoints -= Y[removeIdx];
                NumPoints -= 1.0;
            }

            if (addIdx < static_cast<int>(L) && !std::isnan(Y[addIdx])) {
                SumPoints += Y[addIdx];
                NumPoints += 1.0;
            }

            s[k] = SumPoints;
            np[k] = NumPoints;
        }
    }
    else {  // Even window
        // Initialize sums and counts
        double SumPoints = 0.0;
        double NumPoints = 0.0;

        for (int i = 0; i < halfw && i < static_cast<int>(L); ++i) {
            if (!std::isnan(Y[i])) {
                SumPoints += Y[i];
                NumPoints += 1.0;
            }
        }

        if (halfw < static_cast<int>(L) && !std::isnan(Y[halfw])) {
            SumPoints += 0.5 * Y[halfw];
            NumPoints += 0.5;
        }

        s[0] = SumPoints;
        np[0] = NumPoints;

        for (size_t k = 1; k < L; ++k) {
            int removeIdx1 = static_cast<int>(k) - halfw - 1;
            int removeIdx2 = static_cast<int>(k) - halfw;
            int addIdx1 = static_cast<int>(k) + halfw - 1;
            int addIdx2 = static_cast<int>(k) + halfw;

            if (removeIdx1 >= 0 && !std::isnan(Y[removeIdx1])) {
                SumPoints -= 0.5 * Y[removeIdx1];
                NumPoints -= 0.5;
            }

            if (removeIdx2 >= 0 && !std::isnan(Y[removeIdx2])) {
                SumPoints -= 0.5 * Y[removeIdx2];
                NumPoints -= 0.5;
            }

            if (addIdx1 < static_cast<int>(L) && !std::isnan(Y[addIdx1])) {
                SumPoints += 0.5 * Y[addIdx1];
                NumPoints += 0.5;
            }

            if (addIdx2 < static_cast<int>(L) && !std::isnan(Y[addIdx2])) {
                SumPoints += 0.5 * Y[addIdx2];
                NumPoints += 0.5;
            }

            s[k] = SumPoints;
            np[k] = NumPoints;
        }
    }

    // Remove the amount of interpolated datapoints desired
    double minPoints = std::max(w * (1 - tol), 1.0);
    for (size_t i = 0; i < L; ++i) {
        if (np[i] < minPoints) {
            np[i] = NaN;
        }
    }

    // Calculate smoothed signal
    vector<double> SmoothY(L);
    for (size_t i = 0; i < L; ++i) {
        SmoothY[i] = s[i] / np[i];
    }

    return SmoothY;
}

vector<double> nanfastsmooth(const vector<double>& Y, double w, int type, double tol) {
    int smoothwidth = static_cast<int>(std::round(w));

    switch (type) {
    case 1:
        return sa(Y, smoothwidth, tol);
    case 2:
        return sa(sa(Y, smoothwidth, tol), smoothwidth, tol);
    case 3:
        return sa(sa(sa(Y, smoothwidth, tol), smoothwidth, tol), smoothwidth, tol);
    default:
        return sa(Y, smoothwidth, tol);
    }
}