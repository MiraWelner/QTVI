#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <set>

/**
 * @file interx.hpp
 * @brief Find intersection points of two 2D polyline curves.
 *
 * Port of NS's MATLAB InterX.m (v3.0, 2010).
 * Given two curves as sequences of (x, y) vertices, returns all intersection points.
 */

namespace ppg {

    /** @brief A 2D point. */
    struct Point2D {
        double x = 0.0;
        double y = 0.0;

        bool operator<(const Point2D& o) const {
            if (x != o.x) return x < o.x;
            return y < o.y;
        }

        bool operator==(const Point2D& o) const {
            return std::abs(x - o.x) < 1e-12 && std::abs(y - o.y) < 1e-12;
        }
    };

    /**
     * @brief Find all intersection points between two polyline curves.
     * @param x1  X coordinates of curve 1.
     * @param y1  Y coordinates of curve 1.
     * @param x2  X coordinates of curve 2.
     * @param y2  Y coordinates of curve 2.
     * @return Vector of intersection points (deduplicated).
     **/
    inline std::vector<Point2D> interx(
        const std::vector<double>& x1, const std::vector<double>& y1,
        const std::vector<double>& x2, const std::vector<double>& y2)
    {
        const int n1 = static_cast<int>(x1.size());
        const int n2 = static_cast<int>(x2.size());
        if (n1 < 2 || n2 < 2) return {};

        std::vector<Point2D> result;

        for (int i = 0; i < n1 - 1; ++i) {
            double ax = x1[i], ay = y1[i];
            double bx = x1[i + 1], by = y1[i + 1];
            double dx1 = bx - ax, dy1 = by - ay;

            for (int j = 0; j < n2 - 1; ++j) {
                double cx = x2[j], cy = y2[j];
                double dx = x2[j + 1], dy = y2[j + 1];
                double dx2 = dx - cx, dy2 = dy - cy;

                double denom = dx1 * dy2 - dy1 * dx2;
                if (std::abs(denom) < 1e-15) continue; // parallel

                double t = ((cx - ax) * dy2 - (cy - ay) * dx2) / denom;
                double u = ((cx - ax) * dy1 - (cy - ay) * dx1) / denom;

                if (t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0) {
                    Point2D p;
                    p.x = ax + t * dx1;
                    p.y = ay + t * dy1;
                    result.push_back(p);
                }
            }
        }

        // Deduplicate
        std::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
        return result;
    }

} // namespace ppg