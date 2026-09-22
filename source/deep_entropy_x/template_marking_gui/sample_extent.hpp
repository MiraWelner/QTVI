#pragma once
/**
 * @file   sample_extent.hpp
 * @brief  "Where does this array's usable data start and end" -- once.
 *
 *         This existed in ten places: three global scans at the head of
 *         compute_t_peak / compute_p_peak / compute_p_begin, four inline
 *         window trims, and two in BinPlotWidget. They were not all the same,
 *         and the difference was invisible because no two of them sat in the
 *         same file.
 *
 *         THE FINITE EXTENT AND THE DRAWN EXTENT ARE DIFFERENT THINGS, and
 *         conflating them is a real bug rather than a tidiness matter:
 *
 *           - firstFinite / lastFinite stop at NaN. Every detector in
 *             feature_marks used this.
 *           - firstDrawn / lastDrawn stop at NaN AND then skip columns whose
 *             IQR is 0.0. align_beat_matrix leaves 0.0 wherever fewer than two
 *             beats contributed, and BinPlotWidget::recomputeFrame trims those
 *             off the drawn extent, so the panel will not paint them.
 *
 *         A detector that used the finite extent could therefore place a
 *         landmark in a band the painter refuses to draw -- non-NaN, real
 *         signal, one contributing beat -- and the marker loop's LEFT wall
 *         clamps such a column to the wall instead of skipping it. That is how
 *         a P onset came to sit at the very start of the ECG, hit-tested at a
 *         column it was not drawn at, and impossible to drag.
 *
 *         Naming the two extents separately is the point. Call the one you
 *         mean.
 */

#include <cmath>
#include <vector>

namespace sample_extent {

    double FeatureMarks::sample_at(const std::vector<double>& v, double p) {
        const int n = static_cast<int>(v.size());
        const double NaND = std::numeric_limits<double>::quiet_NaN();
        if (n == 0 || !std::isfinite(p) || p < 0.0 || p > n - 1) return NaND;
        const int i = static_cast<int>(std::floor(p));
        const double f = p - static_cast<double>(i);
        if (f == 0.0) return v[i];
        const double a = v[i], b = v[std::min(n - 1, i + 1)];
        if (std::isnan(a) || std::isnan(b)) return NaND;   // a gap stays a gap
        return a + f * (b - a);
    }

    // First / last non-NaN column; -1 when the array is empty or all NaN.
    inline int firstFinite(const std::vector<double>& v) {
        const int n = static_cast<int>(v.size());
        for (int i = 0; i < n; ++i) if (!std::isnan(v[i])) return i;
        return -1;
    }
    inline int lastFinite(const std::vector<double>& v) {
        for (int i = static_cast<int>(v.size()) - 1; i >= 0; --i)
            if (!std::isnan(v[i])) return i;
        return -1;
    }

    // Both ends at once. False when nothing is finite, in which case `first`
    // and `last` are left at -1 -- the "entirely NaN: no data" case every
    // caller already tested for by hand.
    inline bool finiteExtent(const std::vector<double>& v, int& first, int& last) {
        first = firstFinite(v);
        last = (first < 0) ? -1 : lastFinite(v);
        return first >= 0;
    }

    // The extent BinPlotWidget will actually paint: finite, then past the
    // zero-IQR shoulders. `iqr` shorter than `v` means the trim cannot be
    // applied, so the finite extent is returned unchanged -- the same
    // behaviour as the size check the widget already made.
    inline int firstDrawn(const std::vector<double>& v,
        const std::vector<double>& iqr) {
        int first = firstFinite(v);
        if (first < 0 || iqr.size() != v.size()) return first;
        const int n = static_cast<int>(v.size());
        while (first < n - 1 && iqr[first] == 0.0) ++first;
        return first;
    }
    inline int lastDrawn(const std::vector<double>& v,
        const std::vector<double>& iqr) {
        int last = lastFinite(v);
        if (last < 0 || iqr.size() != v.size()) return last;
        while (last > 0 && iqr[last] == 0.0) --last;
        return last;
    }

    // Shrink [a, b] inward past NaN. False when the window held nothing
    // finite, i.e. the "window was a NaN gap" case; a and b are then crossed
    // and must not be used.
    inline bool trimToFinite(const std::vector<double>& v, int& a, int& b) {
        const int n = static_cast<int>(v.size());
        if (a < 0) a = 0;
        if (b > n - 1) b = n - 1;
        while (a <= b && std::isnan(v[a])) ++a;
        while (b >= a && std::isnan(v[b])) --b;
        return b >= a;
    }

    // Trailing-only trim, for callers that carry their own length rather than
    // an end index (find_foot_pulseox works on a row length). Returns the new
    // length, 0 when the row is entirely NaN.
    inline int finiteLength(const std::vector<double>& v) {
        int len = static_cast<int>(v.size());
        while (len > 0 && std::isnan(v[len - 1])) --len;
        return len;
    }

}   // namespace sample_extent