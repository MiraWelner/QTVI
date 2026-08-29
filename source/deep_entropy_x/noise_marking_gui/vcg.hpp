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
    //   Sign: fixed by convention -- each axis is flipped so it correlates
    //         positively with the sum of the three leads over the whole
    //         trace, which also keeps R generally upright for the detector.
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

        // Einthoven's triangle closure check: for genuine limb leads
        // (ECG1=I, ECG2=II, ECG3=III), III is defined as II - I, so
        // L3 - (L2 - L1) should be ~0 at the DC level. Computed from the
        // SAME per-channel means (b.mean[]) the basis itself was built from
        // (OrthoAccumulator::finish()'s mean over every instant with all
        // three channels present).
        //
        // The common-mode pedestal (the three channels' shared average) is
        // removed FIRST. Without that, a raw ADC-count offset that is nearly
        // IDENTICAL across all three channels (e.g. an unsigned ADC's
        // half-scale zero, un-converted to calibrated mV) mostly cancels
        // between L2 and L1 but leaves one uncancelled copy behind in L3 --
        // so the "closure" reads as roughly the pedestal itself (thousands)
        // instead of the genuine cross-lead residual (a fraction of a unit).
        // Subtracting the common mode from each channel first cancels that
        // shared pedestal exactly, leaving only the actual Einthoven-law
        // violation, whatever preprocessing produced it.
        //
        // A single basis-wide scalar, so it is repeated on every row, the
        // same convention n_instants and planarity_s3_over_s1 already use.
        const double commonMode = (b.mean[0] + b.mean[1] + b.mean[2]) / 3.0;
        const double l1 = b.mean[0] - commonMode;
        const double l2 = b.mean[1] - commonMode;
        const double l3 = b.mean[2] - commonMode;
        const double einthovenClosure = l3 - (l2 - l1);

        f << "subject,axis,ecg1_coeff,ecg2_coeff,ecg3_coeff,"
            "singular_value,rms_per_sample,mean_removed,"
            "n_instants,planarity_s3_over_s1,l3_minus_l2_minus_l1\n";
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
            << b.mean[k] << ',' << b.nUsed << ',' << b.planarity() << ','
            << einthovenClosure << '\n';
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
     * @brief Measures, rather than asks, whether the three recorded leads
     *        are consistent with limb leads (ECG1=I, ECG2=II, ECG3=III) and,
     *        if so, which (if any) is polarity-inverted.
     *
     *        Kirchhoff's voltage law gives II = I + III exactly for genuine
     *        limb leads, i.e. L3 - (L2 - L1) = 0, REGARDLESS of any shared DC
     *        offset -- it is already a property of differences of electrode
     *        potentials. A single inverted lead breaks that identity. Trying
     *        every one of the 2^3 = 8 sign combinations and keeping the one
     *        with the smallest residual finds exactly which lead(s), if any,
     *        need flipping -- a measurement, not a checkbox.
     *
     *        If even the BEST of the 8 combinations has a large residual,
     *        the three channels are not limb leads at all (e.g. one is a
     *        precordial lead), or have a magnitude/calibration mismatch a
     *        sign flip cannot fix -- no sign flip corrects a wrong lead SET
     *        or a wrong per-channel GAIN, only a wrong lead SIGN.
     *
     *        A global flip of all three signs leaves the residual unchanged
     *        (it cancels algebraically), so this cannot detect "all three
     *        inverted together" -- that is a separate, harmless convention
     *        question fixSigns() already handles.
     *
     *        Uses only the covariance the SVD basis already accumulates (no
     *        second pass over the raw samples, and callable on the SAME
     *        OrthoAccumulator mid-accumulation): the RMS residual for signs
     *        (s1,s2,s3) is sqrt(w^T * Cov * w) with w = (s1, -s2, s3), since
     *        residual(t) = s1*L1(t) - s2*L2(t) + s3*L3(t). This is the AC
     *        content only (the covariance is already mean-removed) -- the
     *        DC-level equivalent of this same check is the
     *        l3_minus_l2_minus_l1 column writeBasisCsv
     *        already writes.
     *
     * @param thresh  Normalized-residual cutoff below which the leads are
     *                declared consistent with limb leads. 0.1 is a starting
     *                point, not a validated clinical threshold.
     */
    struct PolarityCheckResult {
        int    sign[3] = { 1, 1, 1 };   ///< best combination found: {ECG1, ECG2, ECG3}; -1 = should be inverted
        double normalizedResidual = std::numeric_limits<double>::quiet_NaN();  ///< for the best combination
        bool   consistentWithLimbLeads = false;
        std::string why;                ///< populated when the check could not run at all
    };

    inline PolarityCheckResult checkLimbLeadPolarity(const OrthoAccumulator& acc, double thresh = 0.10) {
        PolarityCheckResult out;
        if (acc.n < 2) { out.why = "fewer than 2 accumulated instants"; return out; }

        const double inv = 1.0 / static_cast<double>(acc.n);
        double mean[3];
        for (int i = 0; i < 3; ++i) mean[i] = acc.s1[i] * inv;

        const double invDof = 1.0 / static_cast<double>(acc.n - 1);
        double cov[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = i; j < 3; ++j) {
                const double c = (acc.s2[i][j] - static_cast<double>(acc.n) * mean[i] * mean[j]) * invDof;
                cov[i][j] = c;
                cov[j][i] = c;
            }

        const double totalVar = cov[0][0] + cov[1][1] + cov[2][2];
        if (!(totalVar > 0.0)) { out.why = "channels are flat (zero variance)"; return out; }

        out.normalizedResidual = std::numeric_limits<double>::infinity();
        for (int s1 = -1; s1 <= 1; s1 += 2)
            for (int s2 = -1; s2 <= 1; s2 += 2)
                for (int s3 = -1; s3 <= 1; s3 += 2) {
                    const double w[3] = { static_cast<double>(s1), static_cast<double>(-s2), static_cast<double>(s3) };
                    double var = 0.0;
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j)
                            var += w[i] * w[j] * cov[i][j];
                    const double normalized = std::sqrt(std::max(0.0, var) / totalVar);
                    if (normalized < out.normalizedResidual) {
                        out.normalizedResidual = normalized;
                        out.sign[0] = s1; out.sign[1] = s2; out.sign[2] = s3;
                    }
                }
        out.consistentWithLimbLeads = out.normalizedResidual < thresh;

        // Global negation of all three signs leaves the residual unchanged
        // (it cancels algebraically), so the winning combination is only
        // ever meaningful up to that symmetry. Report whichever of the two
        // equivalent representatives has fewer negative signs (always 0 or
        // 1 -- every combo and its global flip have negative-counts that
        // sum to 3, so the smaller of the pair is always <= 1), so the
        // result reads as "at most one lead is inverted" rather than an
        // arbitrary tie-break between two equally-valid descriptions of the
        // same finding.
        const int nNeg = (out.sign[0] < 0) + (out.sign[1] < 0) + (out.sign[2] < 0);
        if (nNeg > 1)
            for (int i = 0; i < 3; ++i) out.sign[i] = -out.sign[i];

        return out;
    }

    /**
     * @brief Apply the sign convention: each axis is flipped, if needed, so
     *        that it moves in the SAME direction as the leads it is built
     *        from, measured over the whole trace rather than guessed from a
     *        single extreme sample.
     *
     *        There is no "correct" sign for a PCA axis the way there is a
     *        correct sign for a limb lead: Einthoven's law gives limb leads
     *        a real right/wrong answer (see checkLimbLeadPolarity()); an
     *        orthonormal rotation has no such constraint, so this remains a
     *        CONVENTION, not a correctness measurement. It is, however, a
     *        more ROBUST convention than "the single biggest sample is
     *        positive": correlating against the sum of all three leads over
     *        every sample is far less sensitive to one large T-wave, noise
     *        spike, or biphasic QRS than looking at one extreme point, and it
     *        ties the derived axis's orientation to leads whose own polarity
     *        has (when possible) already been measured and corrected by
     *        checkLimbLeadPolarity(), rather than to an arbitrary
     *        "which sample happens to be biggest" rule.
     *
     *        Still morphology-dependent in the sense any whole-trace summary
     *        is: build ONE basis per file rather than one per chunk, and use
     *        alignTo() where a rebuild is unavoidable, same as before.
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

        // Project each axis AND build the reference trace (the mean-removed
        // SUM of the three leads) over the same samples, so both line up
        // index-for-index for the dot product below.
        std::vector<double> proj[3];
        std::vector<double> ref;
        for (int k = 0; k < 3; ++k) proj[k].reserve(len);
        ref.reserve(len);
        for (std::size_t i = 0; i < len; ++i) {
            const double v[3] = { (*lanes[0])[i], (*lanes[1])[i], (*lanes[2])[i] };
            if (std::isnan(v[0]) || std::isnan(v[1]) || std::isnan(v[2])) continue;
            const double c[3] = { v[0] - b.mean[0], v[1] - b.mean[1], v[2] - b.mean[2] };
            for (int k = 0; k < 3; ++k)
                proj[k].push_back(b.mat.m[k][0] * c[0]
                    + b.mat.m[k][1] * c[1]
                    + b.mat.m[k][2] * c[2]);
            ref.push_back(c[0] + c[1] + c[2]);
        }

        for (int k = 0; k < 3; ++k) {
            if (proj[k].size() != ref.size() || proj[k].empty()) continue;
            // Sign of the dot product over the WHOLE trace -- equivalent to
            // the sign of the correlation (the normalizing denominator is
            // always positive), so this measures "does this axis move with
            // or against the combined leads, on average", integrated over
            // every sample rather than decided by whichever single sample
            // happens to be largest.
            double dot = 0.0;
            for (std::size_t i = 0; i < proj[k].size(); ++i) dot += proj[k][i] * ref[i];
            if (dot < 0.0)
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