#pragma once
//
// vcg.hpp
//
// VCG reconstruction from the three ECG leads, and the spatial QRS-T angle.
//
// The transform is a 3x3 matrix applied per sample: each of X, Y, Z is a
// weighted sum of the three leads AT THE SAME INSTANT. That per-instant
// requirement is the one real constraint -- every caller has to get its three
// inputs onto a common time base before calling in here, or the combination
// mixes instants that are milliseconds apart.
//
// Conventions:
//   NaN in any lead at instant i => NaN in X, Y and Z at instant i. A partial
//   instant is not a vector; substituting 0 would pull the loop to the origin
//   and read downstream as a real low-voltage segment.
//
//   Amplitudes must be in the SAME units across the three leads (mV). A
//   per-lead normalization applied before this point rescales each axis
//   independently, which rotates and stretches the loop into something that is
//   not a VCG.
//

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace vcg {

    // M_PI is not defined by <cmath> on MSVC without _USE_MATH_DEFINES, so
    // `* 180.0 / M_PI` does not compile in this project. Local constant instead.
    inline constexpr double kPi = 3.14159265358979323846;
    inline constexpr int    kNumLeads = 3;

    struct VCGSample { double x = 0.0, y = 0.0, z = 0.0; };

    /// Rows are X, Y, Z; columns are lead 0, 1, 2 (ECG1, ECG2, ECG3).
    ///   X = m[0][0]*ecg1 + m[0][1]*ecg2 + m[0][2]*ecg3
    struct VcgMatrix {
        double m[3][kNumLeads];
        const char* name = "";

        /// One sample's worth. `v` is indexed by lead.
        VCGSample apply(const double v[kNumLeads]) const {
            VCGSample s;
            for (int k = 0; k < kNumLeads; ++k) {
                s.x += m[0][k] * v[k];
                s.y += m[1][k] * v[k];
                s.z += m[2][k] * v[k];
            }
            return s;
        }

        /// Signature-compatible with the draft: reads instant `i` out of the
        /// three lead signals and transforms it. Returns NaN for all three
        /// components if a lead is missing or short at `i`, so a ragged map
        /// cannot silently contribute a zero.
        VCGSample transform(const std::map<int, std::vector<double>>& leadSignals,
            std::size_t i) const {
            const double kNaN = std::numeric_limits<double>::quiet_NaN();
            double v[kNumLeads];
            for (int k = 0; k < kNumLeads; ++k) {
                auto it = leadSignals.find(k);
                if (it == leadSignals.end() || i >= it->second.size())
                    return VCGSample{ kNaN, kNaN, kNaN };
                v[k] = it->second[i];
                if (std::isnan(v[k])) return VCGSample{ kNaN, kNaN, kNaN };
            }
            return apply(v);
        }
    };

    /// ECG1 -> X, ECG2 -> Y, ECG3 -> Z.
    inline constexpr VcgMatrix kIdentity = { {{1,0,0},{0,1,0},{0,0,1}},
                                             "ECG1/ECG2/ECG3 -> X/Y/Z" };

    struct VcgResult {
        std::vector<VCGSample> samples;
        bool        valid = false;
        std::string why;        ///< populated when !valid, for the log / status bar
        const char* basisName = "";
    };

    /**
     * @brief Per-sample reconstruction over the three leads.
     *
     * @param leadSignals  Keys 0, 1, 2 = ECG1, ECG2, ECG3. All three required,
     *                     all on one common sample grid.
     *
     * @return `valid == false` with `why` set when the inputs cannot support a
     *         reconstruction: fewer than three leads, or unequal lengths. Both
     *         are refused rather than truncated -- unequal length means the
     *         leads are not on one grid, and truncating to the shortest would
     *         silently misalign whatever remains.
     */
    inline VcgResult reconstructVCG(const std::map<int, std::vector<double>>& leadSignals,
        const VcgMatrix& mat = kIdentity) {
        VcgResult out;
        out.basisName = mat.name;

        std::size_t n = 0;
        for (int k = 0; k < kNumLeads; ++k) {
            auto it = leadSignals.find(k);
            if (it == leadSignals.end() || it->second.empty()) {
                out.why = "lead " + std::to_string(k) + " (ECG" + std::to_string(k + 1)
                    + ") is absent";
                return out;
            }
            if (n == 0) n = it->second.size();
            else if (it->second.size() != n) {
                out.why = "lead " + std::to_string(k) + " has "
                    + std::to_string(it->second.size()) + " samples, expected "
                    + std::to_string(n);
                return out;
            }
        }
        if (n < 2) { out.why = "fewer than 2 samples"; return out; }

        out.samples.resize(n);
        for (std::size_t i = 0; i < n; ++i)
            out.samples[i] = mat.transform(leadSignals, i);
        out.valid = true;
        return out;
    }

    /// Convenience form taking the three signals directly.
    inline VcgResult reconstructVCG(const std::vector<const std::vector<double>*>& leads,
        const VcgMatrix& mat = kIdentity) {
        VcgResult out;
        out.basisName = mat.name;
        if (leads.size() != kNumLeads) { out.why = "needs exactly 3 leads"; return out; }
        for (int k = 0; k < kNumLeads; ++k)
            if (!leads[k] || leads[k]->size() < 2) {
                out.why = "ECG" + std::to_string(k + 1) + " is empty"; return out;
            }
        const std::size_t n = leads[0]->size();
        for (int k = 1; k < kNumLeads; ++k)
            if (leads[k]->size() != n) {
                out.why = "ECG" + std::to_string(k + 1) + " has "
                    + std::to_string(leads[k]->size()) + " samples, expected "
                    + std::to_string(n);
                return out;
            }

        const double kNaN = std::numeric_limits<double>::quiet_NaN();
        out.samples.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            const double v[kNumLeads] = { (*leads[0])[i], (*leads[1])[i], (*leads[2])[i] };
            out.samples[i] = (std::isnan(v[0]) || std::isnan(v[1]) || std::isnan(v[2]))
                ? VCGSample{ kNaN, kNaN, kNaN }
            : mat.apply(v);
        }
        out.valid = true;
        return out;
    }

    // -----------------------------------------------------------------------
    // The derived scalar lead
    // -----------------------------------------------------------------------
    /// Which scalar the kors_matrix panel shows and the GUI marks.
    enum class DerivedLead {
        VectorMagnitude,  ///< sqrt(x^2+y^2+z^2). Always >= 0: no negative
        ///< deflections, so QRS polarity logic elsewhere
        ///< must not be applied to it blindly.
        X, Y, Z           ///< a single axis, keeps its sign
    };

    inline const char* derivedLeadName(DerivedLead d) {
        switch (d) {
        case DerivedLead::VectorMagnitude: return "VCG |XYZ|";
        case DerivedLead::X:               return "VCG X";
        case DerivedLead::Y:               return "VCG Y";
        case DerivedLead::Z:               return "VCG Z";
        }
        return "VCG";
    }

    /// Flatten X/Y/Z to the single trace the GUI marks and renders. Same length
    /// and sample grid as the inputs; NaNs pass through as NaNs.
    inline std::vector<double> derivedLeadTrace(const VcgResult& r, DerivedLead which) {
        std::vector<double> out;
        if (!r.valid) return out;
        out.resize(r.samples.size());
        for (std::size_t i = 0; i < r.samples.size(); ++i) {
            const VCGSample& s = r.samples[i];
            switch (which) {
            case DerivedLead::X: out[i] = s.x; break;
            case DerivedLead::Y: out[i] = s.y; break;
            case DerivedLead::Z: out[i] = s.z; break;
            case DerivedLead::VectorMagnitude:
                out[i] = std::sqrt(s.x * s.x + s.y * s.y + s.z * s.z);
                break;
            }
        }
        return out;
    }

    // -----------------------------------------------------------------------
    // Spatial QRS-T angle
    // -----------------------------------------------------------------------
    /**
     * @brief Angle in degrees between the mean QRS vector and the mean T
     *        vector, over [qrsOnset, qrsOffset) and [qrsOffset, tEnd).
     *
     * @return NaN when the windows are empty, out of range, contain no usable
     *         samples, or either mean vector is degenerate. The draft's
     *         `dot / (magQ * magT + 1e-12)` returns 0 -> 90 degrees for a null
     *         vector, which is indistinguishable from a real perpendicular
     *         result; a NaN cannot be mistaken for a measurement.
     */
    inline double spatialQRSTAngle(const std::vector<VCGSample>& v,
        int qrsOnset, int qrsOffset, int tEnd) {
        const double kNaN = std::numeric_limits<double>::quiet_NaN();
        const int N = static_cast<int>(v.size());
        if (qrsOnset < 0 || qrsOffset <= qrsOnset || tEnd <= qrsOffset || tEnd > N)
            return kNaN;

        // NaN-skipping means the divisor is the count of USABLE samples, not
        // the window width -- dividing by the width when part of the window is
        // a dropout shrinks the mean vector toward zero and rotates the angle.
        auto meanOver = [&](int lo, int hi, VCGSample& m) -> bool {
            double sx = 0, sy = 0, sz = 0; int cnt = 0;
            for (int i = lo; i < hi; ++i) {
                const VCGSample& s = v[i];
                if (std::isnan(s.x) || std::isnan(s.y) || std::isnan(s.z)) continue;
                sx += s.x; sy += s.y; sz += s.z; ++cnt;
            }
            if (cnt == 0) return false;
            m.x = sx / cnt; m.y = sy / cnt; m.z = sz / cnt;
            return true;
            };

        VCGSample q, t;
        if (!meanOver(qrsOnset, qrsOffset, q)) return kNaN;
        if (!meanOver(qrsOffset, tEnd, t))     return kNaN;

        const double magQ = std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z);
        const double magT = std::sqrt(t.x * t.x + t.y * t.y + t.z * t.z);
        if (!(magQ > 0.0) || !(magT > 0.0)) return kNaN;   // undefined, not 90

        const double cosA = std::clamp((q.x * t.x + q.y * t.y + q.z * t.z) / (magQ * magT),
            -1.0, 1.0);
        return std::acos(cosA) * 180.0 / kPi;
    }

}  // namespace vcg#pragma once
