
// ============================================================================
// File: stdoutlier.cpp
// ============================================================================
#include "stdoutlier.h"
#include "StatsUtils.h"

vector<bool> stdoutlier(const vector<double>& data,
    double multiplier,
    size_t mean_window,
    const string& direction,
    bool debug_plot) {
    // data == nx1 or 1xn matrix, Multiplier == std multipler (normally ~2.5 or 3)
    // mean_window == length(data) * .02 normally good value

    vector<double> d = diff(data);
    vector<double> x = movmean(d, mean_window);
    double s = std_dev(d);

    vector<double> upper_bound(d.size());
    vector<double> lower_bound(d.size());

    for (size_t i = 0; i < d.size(); ++i) {
        upper_bound[i] = x[i] + s * multiplier;
        lower_bound[i] = x[i] - s * multiplier;
    }

    vector<bool> weridOnes(d.size(), false);

    if (direction == "lower") {
        for (size_t i = 0; i < d.size(); ++i) {
            weridOnes[i] = (d[i] < lower_bound[i]);
        }
    }
    else if (direction == "upper") {
        for (size_t i = 0; i < d.size(); ++i) {
            weridOnes[i] = (d[i] > upper_bound[i]);
        }
    }
    else {
        for (size_t i = 0; i < d.size(); ++i) {
            weridOnes[i] = (d[i] > upper_bound[i] || d[i] < lower_bound[i]);
        }
    }

    vector<bool> outliers(data.size(), false);

    for (size_t i = 0; (i + 1) < data.size(); ++i) {
        if (i < weridOnes.size() && weridOnes[i]) { // Also safety check weridOnes
            outliers[i] = true;
            outliers[i + 1] = true;
        }
    }


    return outliers;
}