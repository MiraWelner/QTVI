#pragma once
//
// vcg.hpp
//
// VCG reconstruction from the three ECG leads, the SVD-derived orthogonal
// basis that supplies the transform, and the spatial QRS-T angle.
//
// NOT A PUBLISHED REGRESSION. Kors, Dower and inverse-Dower are least-squares
// fits FROM the standard 12-lead TO Frank XYZ, reading eight independent leads
// (V1..V6, I, II). Three arbitrary channels are not a subset of that input, so
// there is nothing to substitute in, and a 3x3 slice of an 8-column regression
// is a different, unvalidated transform rather than the regression restricted
// to three leads. The basis used here is derived from the recording itself
// (see the SVD section below) and makes no claim to anatomical X/Y/Z.
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
#include <fstream>
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

    /// The three recorded channels passed through unchanged. NOT an orthogonal
    /// basis: ECG1/2/3 are correlated projections, so calling their components
    /// X/Y/Z asserts an orthogonality the electrodes do not have. Kept for the
    /// Basis::Fixed debug path, and for callers that want only a
    /// rotation-invariant quantity out of it (see derivedLeadTrace's note on
    /// VectorMagnitude), where the choice of basis cannot matter.
    inline constexpr VcgMatrix kIdentity = { {{1,0,0},{0,1,0},{0,0,1}},
                                             "ECG1/ECG2/ECG3 as-is (not orthogonal)" };

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
    /// Which scalar the derived-lead panel shows and the GUI marks. Under an
    /// SVD basis the axis names are positional, not anatomical: X/Y/Z are
    /// PC1/PC2/PC3, ordered by variance. Use orthoAxisName() to label them.
    enum class DerivedLead {
        VectorMagnitude,  ///< sqrt(x^2+y^2+z^2). Always >= 0: no negative
        ///< deflections, so QRS polarity logic elsewhere
        ///< must not be applied to it blindly. INVARIANT
        ///< under any orthonormal basis, so choosing this
        ///< makes the SVD transform a no-op -- see the
        ///< SVD section below before selecting it.
        X, Y, Z           ///< a single axis, keeps its sign. X is the
        ///< max-variance axis (PC1) under an SVD basis.
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

    /// Flatten the three components to the single trace the GUI marks and
    /// renders. Same length and sample grid as the inputs; NaNs pass through as
    /// NaNs. VectorMagnitude is rotation-invariant: it returns the same numbers
    /// whichever orthonormal basis produced `r`, kIdentity included.
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

    // =======================================================================
    // THE SVD-ORTHOGONAL BASIS
    // =======================================================================
    //
    // A DATA-DRIVEN orthogonal basis for the three recorded ECG channels,
    // supplying the VcgMatrix that everything above applies.
    //
    // ---------------------------------------------------------------------------
    // WHAT THE SVD ACTUALLY BUYS, AND WHAT IT DOES NOT
    // ---------------------------------------------------------------------------
    // The basis here is the right singular vectors of the centred n-by-3 lead
    // matrix: three mutually orthogonal unit directions in lead space, ordered by
    // how much signal variance projects onto each. Equivalently the eigenvectors
    // of the 3x3 lead covariance, which is how it is computed (n is in the
    // millions; the 3x3 is not).
    //
    // READ THIS BEFORE CHANGING THE PANEL'S DerivedLead:
    //
    //   An orthonormal change of basis does NOT change the vector magnitude
    //   sqrt(x^2+y^2+z^2), and does NOT change the angle between two vectors.
    //
    // So DerivedLead::VectorMagnitude and vcg::spatialQRSTAngle produce bitwise
    // the same numbers under this basis as under kIdentity. If the panel shows
    // |XYZ|, switching to SVD changes nothing at all. What the basis does buy:
    //
    //   PC1 (row 0)  The unit-norm linear combination of the three channels with
    //                maximum variance -- the best-SNR single scalar lead the three
    //                channels can produce, which is the useful thing to run a peak
    //                detector on and the recommended DerivedLead here.
    //   PC2, PC3     Decorrelated residual content. Under kIdentity, ECG1/2/3 are
    //                correlated projections and the "axes" duplicate each other's
    //                information; these do not.
    //   sigma[]      A dimensionality measure that kIdentity cannot give:
    //                sigma[2]/sigma[0] near zero means the loop is essentially
    //                planar, i.e. the three channels span only two useful
    //                directions and the third "axis" is noise.
    //
    // ---------------------------------------------------------------------------
    // THE TWO INSTABILITIES, AND WHY THE BASIS IS BUILT ONCE PER FILE
    // ---------------------------------------------------------------------------
    // An SVD basis is only defined up to (a) the sign of each vector and (b) the
    // ordering of vectors with equal singular values. Both matter here because the
    // derived lead is MARKABLE: an annotation placed on it in one chunk has to
    // mean the same linear combination in the next chunk, or the mark refers to a
    // signal that no longer exists.
    //
    //   Sign: fixed by convention -- each axis is flipped so its dominant
    //         deflection is positive, which also keeps R upright for the detector.
    //         See fixSigns().
    //   Order: unstable when two singular values are close. degeneracy() reports
    //         the separation; alignTo() locks a new basis onto a previous one by
    //         matching directions instead of re-sorting by variance.
    //
    // The intended use is therefore OrthoAccumulator over the WHOLE file (or a long
    // representative stretch), finish() once, then that one basis applied to every
    // chunk. Rebuilding per chunk gives a channel whose meaning drifts as the
    // operator scrolls. Accumulator is incremental precisely so the basis can be
    // gathered across chunk loads before any marking begins.
    //
    // Amplitude convention: the basis is orthonormal, so the output is in the same
    // units as the input (mV) and no channel is rescaled relative to another.
    // Whitening (dividing each axis by its singular value) is deliberately NOT
    // offered: it rescales the axes independently, which is exactly the operation
    // this file's header warns turns a VCG into something that is not a VCG.
    //


    // -----------------------------------------------------------------------
    // Symmetric 3x3 eigendecomposition (cyclic Jacobi)
    // -----------------------------------------------------------------------
    /**
     * @brief Eigendecomposition of a symmetric 3x3, in place on a copy.
     *
     * @param in    Symmetric input. Only the upper triangle needs to be
     *              correct; the lower is assumed to mirror it.
     * @param eval  Eigenvalues, in the order the rotations happened to leave
     *              them (NOT sorted -- the caller sorts).
     * @param evec  evec[k] is the unit eigenvector for eval[k], as a ROW.
     *
     * Jacobi rather than the closed-form cubic: the cubic loses most of its
     * precision on the near-degenerate covariances that correlated ECG
     * channels routinely produce, and near-degeneracy is exactly the case
     * whose eigenvectors this code has to report honestly.
     */
    inline void symmetricEigen3(const double in[3][3],
        double eval[3], double evec[3][3]) {
        double a[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) a[i][j] = in[i][j];

        double v[3][3] = { {1,0,0},{0,1,0},{0,0,1} };   // columns are eigenvectors

        for (int sweep = 0; sweep < 32; ++sweep) {
            const double off = std::fabs(a[0][1]) + std::fabs(a[0][2]) + std::fabs(a[1][2]);
            const double diag = std::fabs(a[0][0]) + std::fabs(a[1][1]) + std::fabs(a[2][2]);
            if (off <= 1e-18 * diag) break;             // converged (also handles all-zero)

            for (int p = 0; p < 2; ++p) {
                for (int q = p + 1; q < 3; ++q) {
                    if (std::fabs(a[p][q]) <= 1e-300) continue;
                    const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                    const double t = std::copysign(1.0, theta)
                        / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                    const double c = 1.0 / std::sqrt(t * t + 1.0);
                    const double s = t * c;

                    for (int k = 0; k < 3; ++k) {       // A <- A * J
                        const double kp = a[k][p], kq = a[k][q];
                        a[k][p] = c * kp - s * kq;
                        a[k][q] = s * kp + c * kq;
                    }
                    for (int k = 0; k < 3; ++k) {       // A <- J^T * A
                        const double pk = a[p][k], qk = a[q][k];
                        a[p][k] = c * pk - s * qk;
                        a[q][k] = s * pk + c * qk;
                    }
                    for (int k = 0; k < 3; ++k) {       // V <- V * J
                        const double kp = v[k][p], kq = v[k][q];
                        v[k][p] = c * kp - s * kq;
                        v[k][q] = s * kp + c * kq;
                    }
                }
            }
        }

        for (int k = 0; k < 3; ++k) {
            eval[k] = a[k][k];
            for (int i = 0; i < 3; ++i) evec[k][i] = v[i][k];   // column k -> row k
        }
    }

    // -----------------------------------------------------------------------
    // The basis
    // -----------------------------------------------------------------------
    /**
     * @brief An orthonormal 3x3 basis derived from the data, plus the evidence
     *        needed to decide whether to trust it.
     *
     *        `mat` is directly usable as vcg::VcgMatrix: row k is axis k, so
     *        VcgMatrix::apply() projects a lead triple onto the basis. Rows are
     *        ordered by descending singular value, so row 0 is PC1.
     */
    struct OrthoBasis {
        VcgMatrix   mat{};                  ///< rows = orthonormal axes, PC1 first
        double      sigma[3] = { 0,0,0 };   ///< singular values of the centred matrix, descending
        double      mean[3] = { 0,0,0 };    ///< per-channel mean that was removed
        long long   nUsed = 0;              ///< instants with all three channels finite
        bool        valid = false;
        std::string why;                    ///< populated when !valid
        std::string label;                  ///< human-readable, for the panel title

        /// sigma[2]/sigma[0]. Near 0 => the loop is planar and axis 3 is noise.
        double planarity() const {
            return (sigma[0] > 0.0) ? sigma[2] / sigma[0]
                : std::numeric_limits<double>::quiet_NaN();
        }
        /// Smallest relative gap between consecutive singular values. Below a
        /// few percent, the axis ORDER is not reliably reproducible and a
        /// rebuild may permute the axes; use alignTo() rather than re-sorting.
        double degeneracy() const {
            if (!(sigma[0] > 0.0)) return std::numeric_limits<double>::quiet_NaN();
            return std::min(sigma[0] - sigma[1], sigma[1] - sigma[2]) / sigma[0];
        }
    };

    /**
     * @brief Write a basis to CSV. Returns false and fills `why` when it
     *        cannot.
     *
     *        Lives here because the format belongs with the struct it
     *        serializes: anything that adds a field to OrthoBasis sees this
     *        function in the same file. Takes a full PATH rather than a
     *        directory or a config object, so this header keeps depending on
     *        nothing outside the standard library.
     */
    inline bool writeBasisCsv(const OrthoBasis& b,
        const std::string& path,
        const std::string& subject,
        std::string& why) {
        if (!b.valid) {
            why = b.why.empty() ? "basis not computed" : b.why;
            return false;
        }
        std::ofstream f(path, std::ios::trunc);
        if (!f) { why = "cannot open " + path; return false; }

        f << "subject,axis,ecg1_coeff,ecg2_coeff,ecg3_coeff,"
            "singular_value,rms_per_sample,mean_removed,"
            "n_instants,planarity_s3_over_s1\n";
        f.setf(std::ios::fixed);
        f.precision(9);
        // sigma is the singular value of the CENTRED matrix, so it scales with
        // sqrt(n) -- at millions of samples it reads in the hundreds of
        // thousands. rms_per_sample is the same quantity in the input's own
        // units, and the only one of the two comparable across recordings.
        const double rootN = (b.nUsed > 1)
            ? std::sqrt(static_cast<double>(b.nUsed - 1)) : 1.0;
        for (int k = 0; k < 3; ++k)
            f << subject << ",PC" << (k + 1) << ','
            << b.mat.m[k][0] << ',' << b.mat.m[k][1] << ',' << b.mat.m[k][2] << ','
            << b.sigma[k] << ',' << (b.sigma[k] / rootN) << ','
            << b.mean[k] << ',' << b.nUsed << ',' << b.planarity() << '\n';
        return static_cast<bool>(f);
    }

    /**
     * @brief Incremental covariance accumulator, so one basis can be built
     *        across many chunk loads without holding the whole file.
     *
     *        Uses raw sums of products rather than Welford. ECG in mV over even
     *        a long recording keeps the sums far inside double's range, and the
     *        means are small relative to the deflections, so the cancellation
     *        that makes the raw form dangerous in general does not bite here.
     *        add() skips any instant where a channel is NaN: a partial instant
     *        is not a vector, and letting it through as 0 would drag the
     *        covariance toward the origin and rotate the basis.
     */
    struct OrthoAccumulator {
        double    s1[3] = { 0,0,0 };        ///< sum of each channel
        double    s2[3][3] = { {0,0,0},{0,0,0},{0,0,0} };  ///< sum of products
        long long n = 0;

        void reset() { *this = OrthoAccumulator{}; }

        /// One instant. `v` is indexed by channel (0,1,2).
        void addSample(const double v[3]) {
            if (std::isnan(v[0]) || std::isnan(v[1]) || std::isnan(v[2])) return;
            for (int i = 0; i < 3; ++i) {
                s1[i] += v[i];
                for (int j = i; j < 3; ++j) s2[i][j] += v[i] * v[j];
            }
            ++n;
        }

        /// A whole chunk. Lanes must be equal length and on ONE time grid --
        /// see vcg_lead.hpp, which resamples before it gets here.
        void addLanes(const std::vector<const std::vector<double>*>& lanes) {
            if (lanes.size() != 3 || !lanes[0] || !lanes[1] || !lanes[2]) return;
            const std::size_t len = std::min({ lanes[0]->size(),
                                               lanes[1]->size(),
                                               lanes[2]->size() });
            for (std::size_t k = 0; k < len; ++k) {
                const double v[3] = { (*lanes[0])[k], (*lanes[1])[k], (*lanes[2])[k] };
                addSample(v);
            }
        }

        /**
         * @brief Reduce to a basis.
         *
         * @param minSamples  Refuse below this many usable instants. The
         *                    default is deliberately not 3 (the algebraic
         *                    minimum): a covariance from a handful of samples
         *                    is dominated by whatever noise happened to be in
         *                    them, and the resulting "axes" are arbitrary.
         *
         *        Signs are NOT set here -- finish() has no access to the
         *        samples, and the convention needs the projected traces. Call
         *        fixSigns() with the data, or alignTo() with a previous basis.
         */
        OrthoBasis finish(long long minSamples = 1000) const {
            OrthoBasis out;
            out.nUsed = n;
            if (n < minSamples) {
                out.why = "only " + std::to_string(n)
                    + " instants with all three channels present; need "
                    + std::to_string(minSamples);
                return out;
            }

            const double inv = 1.0 / static_cast<double>(n);
            for (int i = 0; i < 3; ++i) out.mean[i] = s1[i] * inv;

            // Covariance, normalised by n-1.
            const double invDof = 1.0 / static_cast<double>(n - 1);
            double cov[3][3];
            for (int i = 0; i < 3; ++i)
                for (int j = i; j < 3; ++j) {
                    const double c = (s2[i][j] - static_cast<double>(n) * out.mean[i] * out.mean[j]) * invDof;
                    cov[i][j] = c;
                    cov[j][i] = c;
                }

            double eval[3], evec[3][3];
            symmetricEigen3(cov, eval, evec);

            int order[3] = { 0,1,2 };
            std::sort(order, order + 3, [&](int a, int b) { return eval[a] > eval[b]; });

            if (!(eval[order[0]] > 0.0)) {
                out.why = "lead covariance is degenerate (all three channels flat)";
                return out;
            }

            for (int k = 0; k < 3; ++k) {
                const int e = order[k];
                // Rounding can leave a near-zero eigenvalue slightly negative;
                // clamp rather than emit NaN from the sqrt.
                const double lam = std::max(0.0, eval[e]);
                out.sigma[k] = std::sqrt(lam * static_cast<double>(n - 1));
                double norm = 0.0;
                for (int i = 0; i < 3; ++i) norm += evec[e][i] * evec[e][i];
                norm = (norm > 0.0) ? std::sqrt(norm) : 1.0;
                for (int i = 0; i < 3; ++i) out.mat.m[k][i] = evec[e][i] / norm;
            }

            out.mat.name = "SVD-orthogonal (PC1/PC2/PC3)";
            out.label = "SVD-orthogonal (PC1/PC2/PC3)";
            out.valid = true;
            return out;
        }
    };

    /**
     * @brief Apply the sign convention: each axis is flipped, if needed, so
     *        that its largest excursion from the projected median is POSITIVE.
     *
     *        Rationale: the largest excursion in an ECG lead is the QRS, so
     *        this puts R upright, which is what the existing peak detector and
     *        the operator both expect. It is a convention, not a measurement --
     *        the SVD has no opinion about which end of an axis is "up".
     *
     *        It is also morphology-dependent: a chunk where the T wave
     *        outgrows the R (or a lead whose QRS is genuinely biphasic) can
     *        flip relative to another chunk. That is the reason to build ONE
     *        basis per file rather than one per chunk, and the reason alignTo()
     *        exists for the case where a rebuild is unavoidable.
     *
     * @param lanes  The same three lanes the basis was accumulated from.
     */
    inline void fixSigns(OrthoBasis& b,
        const std::vector<const std::vector<double>*>& lanes) {
        if (!b.valid || lanes.size() != 3) return;
        if (!lanes[0] || !lanes[1] || !lanes[2]) return;
        const std::size_t len = std::min({ lanes[0]->size(),
                                           lanes[1]->size(),
                                           lanes[2]->size() });
        if (len == 0) return;

        // Project, collect per-axis values, take the median as baseline, then
        // the signed extreme relative to it.
        std::vector<double> proj[3];
        for (int k = 0; k < 3; ++k) proj[k].reserve(len);
        for (std::size_t i = 0; i < len; ++i) {
            const double v[3] = { (*lanes[0])[i], (*lanes[1])[i], (*lanes[2])[i] };
            if (std::isnan(v[0]) || std::isnan(v[1]) || std::isnan(v[2])) continue;
            const double c[3] = { v[0] - b.mean[0], v[1] - b.mean[1], v[2] - b.mean[2] };
            for (int k = 0; k < 3; ++k)
                proj[k].push_back(b.mat.m[k][0] * c[0]
                    + b.mat.m[k][1] * c[1]
                    + b.mat.m[k][2] * c[2]);
        }

        for (int k = 0; k < 3; ++k) {
            if (proj[k].empty()) continue;
            std::vector<double> tmp = proj[k];
            const auto mid = tmp.begin() + tmp.size() / 2;
            std::nth_element(tmp.begin(), mid, tmp.end());
            const double med = *mid;

            double lo = 0.0, hi = 0.0;
            for (double x : proj[k]) {
                const double d = x - med;
                if (d > hi) hi = d;
                if (d < lo) lo = d;
            }
            if (-lo > hi)                                  // dominant deflection is negative
                for (int i = 0; i < 3; ++i) b.mat.m[k][i] = -b.mat.m[k][i];
        }
    }

    /**
     * @brief Lock a freshly computed basis onto a previous one: permute and
     *        flip `b`'s axes so each matches the corresponding axis of `prev`
     *        as closely as possible, instead of trusting variance order.
     *
     *        Use when the basis must be recomputed on data that has changed
     *        (new file, or a rebuild after a long absent stretch) but the
     *        derived channel's existing annotations still have to mean the same
     *        thing. Greedy: strongest |dot| pairing first. sigma[] is permuted
     *        with the rows so it keeps describing its own axis.
     *
     * @return Smallest |dot| among the three pairings, i.e. how well the new
     *         basis could be matched to the old one. Near 1 means the axes
     *         barely moved; near 0 means the recording's geometry genuinely
     *         changed and marks placed on the old axes should be treated as
     *         suspect rather than silently carried forward.
     */
    inline double alignTo(OrthoBasis& b, const OrthoBasis& prev) {
        if (!b.valid || !prev.valid) return std::numeric_limits<double>::quiet_NaN();

        double dot[3][3];
        for (int p = 0; p < 3; ++p)
            for (int k = 0; k < 3; ++k)
                dot[p][k] = prev.mat.m[p][0] * b.mat.m[k][0]
                + prev.mat.m[p][1] * b.mat.m[k][1]
                + prev.mat.m[p][2] * b.mat.m[k][2];

        int pick[3] = { -1,-1,-1 };
        bool takenPrev[3] = { false,false,false }, takenNew[3] = { false,false,false };
        double worst = 1.0;
        for (int step = 0; step < 3; ++step) {
            int bp = -1, bk = -1; double best = -1.0;
            for (int p = 0; p < 3; ++p) {
                if (takenPrev[p]) continue;
                for (int k = 0; k < 3; ++k) {
                    if (takenNew[k]) continue;
                    const double d = std::fabs(dot[p][k]);
                    if (d > best) { best = d; bp = p; bk = k; }
                }
            }
            if (bp < 0) break;
            takenPrev[bp] = takenNew[bk] = true;
            pick[bp] = bk;
            worst = std::min(worst, best);
        }

        VcgMatrix out{};
        out.name = b.mat.name;
        double sig[3] = { 0,0,0 };
        for (int p = 0; p < 3; ++p) {
            const int k = (pick[p] >= 0) ? pick[p] : p;
            const double sign = (dot[p][k] < 0.0) ? -1.0 : 1.0;
            for (int i = 0; i < 3; ++i) out.m[p][i] = sign * b.mat.m[k][i];
            sig[p] = b.sigma[k];
        }
        b.mat = out;
        for (int p = 0; p < 3; ++p) b.sigma[p] = sig[p];
        return worst;
    }

    /**
     * @brief One-shot: accumulate, reduce, apply the sign convention.
     *        Convenient for a single chunk; prefer OrthoAccumulator across the
     *        whole file for a channel that will be marked.
     */
    inline OrthoBasis computeOrthoBasis(
        const std::vector<const std::vector<double>*>& lanes,
        long long minSamples = 1000) {
        OrthoAccumulator acc;
        acc.addLanes(lanes);
        OrthoBasis b = acc.finish(minSamples);
        fixSigns(b, lanes);
        return b;
    }

    /// Axis names under this basis. The X/Y/Z names in DerivedLead are
    /// anatomical and do not apply; row k is PC(k+1).
    inline const char* orthoAxisName(DerivedLead d) {
        switch (d) {
        case DerivedLead::X: return "PC1 (max-variance)";
        case DerivedLead::Y: return "PC2";
        case DerivedLead::Z: return "PC3";
        case DerivedLead::VectorMagnitude: return "|PC1,PC2,PC3| (= |ECG1,ECG2,ECG3|)";
        }
        return "PC?";
    }

}  // namespace vcg