#pragma once
/**
 * @file   envelopes.hpp
 * @brief  Section 4.7 DYNAMIC envelopes: moving-average, power-spectral and
 *         wavelet envelopes of the P, QRS and T waves over short (5-10 beat)
 *         and long (30 beat) windows.
 *
 *         DISTINCT FROM morphology_envelope.hpp. That builds ONE STATIC
 *         corridor per bin -- a 2.5/97.5 band over the bin's clean beats, used
 *         to score a beat against its own bin. This builds ROLLING statistics
 *         that move with the record, which is what "dynamic" means and what
 *         makes them usable by the dynamic master template selection in
 *         Section 9.5. A static corridor cannot detect drift: it is computed
 *         from the very beats that drifted.
 *
 *         WHY TWO WINDOW LENGTHS. Short windows (5-10 beats) catch PVC onset
 *         and conduction change -- events that appear and are gone inside a few
 *         beats, and which a 30-beat window averages away. Long windows (30
 *         beats) catch ischemic ST drift, which is monotone and small per beat
 *         and which a short window cannot distinguish from noise. Neither
 *         length subsumes the other, so both are computed for every measure and
 *         both are reported.
 *
 *         WHY THREE FAMILIES. They respond to different kinds of change:
 *           - MOVING AVERAGE tracks amplitude and shape in the time domain. It
 *             sees a beat that got taller or a segment that shifted level.
 *           - POWER SPECTRAL tracks where the segment's energy sits in
 *             frequency. A widened QRS moves energy downward without
 *             necessarily changing its amplitude, so the moving average can
 *             miss entirely what this reports immediately.
 *           - WAVELET tracks energy per scale with time localisation, which is
 *             what neither of the other two has: a spectrum is blind to WHERE
 *             in the segment the change happened, and a moving average is blind
 *             to its scale.
 *
 *         The output is deliberately per beat per window, not per bin: 9.5
 *         selects a master template at a point in the record, so it needs the
 *         value AT that point plus the local mean and spread to judge it
 *         against.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <vector>

namespace envelopes {

    // Section 4.7 window lengths. The spec gives the short window as a range
    // (5 to 10); 8 sits in the middle and is a power of two, which matters for
    // the wavelet decomposition below.
    inline constexpr int kShortWindow = 8;
    inline constexpr int kLongWindow = 30;

    // Wavelet decomposition depth. Three Haar levels on a segment of a few tens
    // of samples gives scales of 2, 4 and 8 samples -- fine enough to separate a
    // notch from a slur, coarse enough that the top level still has support.
    inline constexpr int kWaveletLevels = 3;

    // Spectral bands as fractions of Nyquist. Deliberately coarse: the point is
    // to detect energy MOVING between low and high, not to resolve a spectrum
    // whose frequency step is set by however many samples the segment happens
    // to span.
    inline constexpr int kSpectralBands = 3;

    enum class Segment : uint8_t { P = 0, QRS = 1, T = 2 };
    inline const char* segmentName(Segment s) {
        switch (s) {
        case Segment::P:   return "P";
        case Segment::QRS: return "QRS";
        default:           return "T";
        }
    }

    // Sample span of each wave within a beat, from the template's landmarks.
    // Half-open [begin, end). A segment whose landmarks are missing (-1) is
    // absent rather than empty: a ventricular beat legitimately has no P wave,
    // and every P-dependent measure must come out NaN rather than 0.
    struct SegmentSpans {
        int p_begin = -1, p_end = -1;
        int qrs_begin = -1, qrs_end = -1;
        int t_begin = -1, t_end = -1;

        bool has(Segment s) const {
            switch (s) {
            case Segment::P:   return p_begin >= 0 && p_end > p_begin;
            case Segment::QRS: return qrs_begin >= 0 && qrs_end > qrs_begin;
            default:           return t_begin >= 0 && t_end > t_begin;
            }
        }
        int begin(Segment s) const {
            switch (s) {
            case Segment::P:   return p_begin;
            case Segment::QRS: return qrs_begin;
            default:           return t_begin;
            }
        }
        int end(Segment s) const {
            switch (s) {
            case Segment::P:   return p_end;
            case Segment::QRS: return qrs_end;
            default:           return t_end;
            }
        }
    };

    // ---------------------------------------------------------------------
    // Per-beat, per-segment measures
    // ---------------------------------------------------------------------

    struct BeatMeasures {
        // Moving-average family: time-domain amplitude and level.
        double mean = std::numeric_limits<double>::quiet_NaN();
        double amplitude = std::numeric_limits<double>::quiet_NaN();  // max - min
        double area = std::numeric_limits<double>::quiet_NaN();  // sum|x| per sample

        // Power-spectral family: band powers as FRACTIONS of total, plus the
        // spectral centroid. Fractions rather than absolute powers because a
        // beat that simply got taller raises every band equally -- that is the
        // moving average's job to report, and mixing it in here would make the
        // spectral measure a worse amplitude detector instead of a shape one.
        double band_frac[kSpectralBands] = {
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN() };
        double centroid = std::numeric_limits<double>::quiet_NaN();

        // Wavelet family: Haar detail energy per level, as a fraction of total
        // detail energy, for the same reason the spectral bands are fractions.
        double wave_frac[kWaveletLevels] = {
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN() };

        // ---- INTERMEDIATES ------------------------------------------------
        // The quantities the fractions above are fractions OF, kept so a
        // reported value can be audited without re-running the math. Every
        // fraction hides its denominator, and two beats with identical
        // band_frac can have total powers an order of magnitude apart -- which
        // is the difference between a shape change and a quiet stretch. Filled
        // by measure(); nothing routes on them.
        //
        // Recorded rather than recomputed on purpose: a second copy of the
        // spectral/wavelet math in a reporting path would drift from this one.
        int    n_samples = 0;     // non-NaN samples the segment contributed
        double spectral_total_power = std::numeric_limits<double>::quiet_NaN();
        double band_power[kSpectralBands] = {
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN() };
        double wave_energy[kWaveletLevels] = {
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN() };
        double wave_total_energy = std::numeric_limits<double>::quiet_NaN();

        bool valid = false;   // false when the segment was absent or too short
    };

    namespace detail {

        // Collect a segment's non-NaN samples. The shared axis is NaN-padded, so
        // a segment near a beat's edge can be partly padding.
        inline std::vector<double> gather(const std::vector<double>& beat,
            int begin, int end)
        {
            std::vector<double> v;
            if (begin < 0 || end <= begin) return v;
            const int hi = std::min<int>(end, static_cast<int>(beat.size()));
            v.reserve(std::max(0, hi - begin));
            for (int i = std::max(0, begin); i < hi; ++i)
                if (!std::isnan(beat[i])) v.push_back(beat[i]);
            return v;
        }

        // Real DFT magnitude spectrum, bins 1..N/2. Direct evaluation rather
        // than an FFT: a P/QRS/T segment is a few tens of samples, so O(n^2) on
        // n = 40 is under two thousand multiply-adds and an FFT's setup would
        // cost more than it saves. Bin 0 (DC) is skipped -- it is the segment's
        // mean, which the moving-average family already reports, and leaving it
        // in would make every band fraction a function of baseline level.
        inline std::vector<double> magnitudeSpectrum(const std::vector<double>& x)
        {
            const int n = static_cast<int>(x.size());
            std::vector<double> mag;
            if (n < 4) return mag;

            double mean = 0.0;
            for (double v : x) mean += v;
            mean /= n;

            const int nb = n / 2;
            mag.resize(nb, 0.0);

            // ---- TWIDDLE TABLE, CACHED BY LENGTH ---------------------------
            // This loop is still O(n^2), but it no longer calls std::cos and
            // std::sin inside it. That mattered enormously: measure() runs per
            // segment per channel per beat, so a night at 40 bins is on the
            // order of a million calls here, each previously doing n^2/2
            // transcendental PAIRS. At n ~ 200 that is ~10^10 trig evaluations
            // for one subject -- tens of minutes, indistinguishable from a hang.
            //
            // cos(w*k*t) depends only on (k*t) mod n, so one table of n entries
            // serves every k and t. Segment lengths repeat constantly (beats in
            // a bin are near-identical in span), so the table is rebuilt only
            // when the length actually changes.
            //
            // thread_local, not static: measure() is callable from parallel
            // per-bin loops, and a shared mutable table would be a race that
            // corrupts spectra rather than crashing -- the failure mode that
            // produces plausible numbers.
            thread_local std::vector<double> tw_cos, tw_sin;
            thread_local int tw_n = -1;
            if (tw_n != n) {
                tw_cos.resize(n);
                tw_sin.resize(n);
                const double base = -2.0 * 3.14159265358979323846 / n;
                for (int j = 0; j < n; ++j) {
                    tw_cos[j] = std::cos(base * j);
                    tw_sin[j] = std::sin(base * j);
                }
                tw_n = n;
            }

            for (int k = 1; k <= nb; ++k) {
                double re = 0.0, im = 0.0;
                int idx = 0;                  // (k*t) mod n, advanced by k
                for (int t = 0; t < n; ++t) {
                    const double xt = x[t] - mean;   // DC removed
                    re += xt * tw_cos[idx];
                    im += xt * tw_sin[idx];
                    idx += k;
                    if (idx >= n) idx -= n;   // cheaper than a modulo per sample
                }
                mag[k - 1] = std::sqrt(re * re + im * im);
            }
            return mag;
        }

        // Haar detail energy per level. Averaging and differencing in place:
        // each level halves the approximation and yields one detail band, so
        // level 0 is the finest scale (adjacent-sample differences) and higher
        // levels are progressively coarser.
        inline void haarDetailEnergy(std::vector<double> x,
            double out_energy[kWaveletLevels])
        {
            for (int L = 0; L < kWaveletLevels; ++L)
                out_energy[L] = std::numeric_limits<double>::quiet_NaN();

            for (int L = 0; L < kWaveletLevels; ++L) {
                const int n = static_cast<int>(x.size()) / 2;
                if (n < 2) return;   // no support left at this scale
                std::vector<double> approx(n), detail(n);
                for (int i = 0; i < n; ++i) {
                    const double a = x[2 * i], b = x[2 * i + 1];
                    approx[i] = 0.5 * (a + b);
                    detail[i] = 0.5 * (a - b);
                }
                double e = 0.0;
                for (double d : detail) e += d * d;
                out_energy[L] = e;
                x = std::move(approx);
            }
        }

    }  // namespace detail

    // All three families for one segment of one beat.
    inline BeatMeasures measure(const std::vector<double>& beat,
        const SegmentSpans& spans, Segment seg)
    {
        BeatMeasures m;
        if (!spans.has(seg)) return m;   // absent landmark: stays NaN, not zero

        const std::vector<double> x =
            detail::gather(beat, spans.begin(seg), spans.end(seg));
        if (x.size() < 8) return m;      // too short for a 3-level decomposition

        // ---- moving-average family --------------------------------------
        double sum = 0.0, absum = 0.0;
        double lo = x[0], hi = x[0];
        for (double v : x) {
            sum += v; absum += std::fabs(v);
            lo = std::min(lo, v); hi = std::max(hi, v);
        }
        const double n = static_cast<double>(x.size());
        m.n_samples = static_cast<int>(x.size());
        m.mean = sum / n;
        m.amplitude = hi - lo;
        m.area = absum / n;

        // ---- power-spectral family --------------------------------------
        const std::vector<double> mag = detail::magnitudeSpectrum(x);
        if (!mag.empty()) {
            double total = 0.0, weighted = 0.0;
            for (size_t k = 0; k < mag.size(); ++k) {
                const double p = mag[k] * mag[k];
                total += p;
                weighted += p * static_cast<double>(k + 1);
            }
            m.spectral_total_power = total;
            if (total > 0.0) {
                double bandPow[kSpectralBands] = { 0.0, 0.0, 0.0 };
                for (size_t k = 0; k < mag.size(); ++k) {
                    int b = static_cast<int>((kSpectralBands * k) / mag.size());
                    b = std::min(b, kSpectralBands - 1);
                    bandPow[b] += mag[k] * mag[k];
                }
                for (int b = 0; b < kSpectralBands; ++b) {
                    m.band_power[b] = bandPow[b];
                    m.band_frac[b] = bandPow[b] / total;
                }
                // Centroid normalised to [0,1] of Nyquist so it is comparable
                // across segments of different lengths -- a QRS span and a T
                // span do not have the same frequency resolution.
                m.centroid = (weighted / total) / static_cast<double>(mag.size());
            }
        }

        // ---- wavelet family ---------------------------------------------
        double we[kWaveletLevels];
        detail::haarDetailEnergy(x, we);
        double wtot = 0.0;
        for (int L = 0; L < kWaveletLevels; ++L)
            if (!std::isnan(we[L])) wtot += we[L];
        for (int L = 0; L < kWaveletLevels; ++L) m.wave_energy[L] = we[L];
        m.wave_total_energy = wtot;
        if (wtot > 0.0)
            for (int L = 0; L < kWaveletLevels; ++L)
                if (!std::isnan(we[L])) m.wave_frac[L] = we[L] / wtot;

        m.valid = true;
        return m;
    }

    // ---------------------------------------------------------------------
    // Rolling envelope
    // ---------------------------------------------------------------------

    // The envelope proper: local mean and spread of one measure over a window
    // of preceding beats, plus the current beat's deviation from it in SDs.
    //
    // z is the output 9.5 actually consumes. A raw value cannot say whether a
    // beat is unusual without something to compare it to, and the comparison
    // has to be LOCAL or drift makes every late beat look abnormal against an
    // early baseline.
    struct EnvelopePoint {
        double value = std::numeric_limits<double>::quiet_NaN();
        double mean = std::numeric_limits<double>::quiet_NaN();
        double sd = std::numeric_limits<double>::quiet_NaN();
        double z = std::numeric_limits<double>::quiet_NaN();
        int    n = 0;      // beats actually contributing
        bool   ready = false;  // window filled
    };

    // Trailing window, excluding the current beat. Excluding it matters: a beat
    // included in its own baseline pulls the mean toward itself, which is
    // exactly wrong for onset detection -- the first beat of a run would partly
    // mask itself, and with a 5-beat window it would mask a fifth of the change.
    class RollingEnvelope {
    public:
        explicit RollingEnvelope(int window) : m_window(std::max(2, window)) {}

        EnvelopePoint push(double value) {
            EnvelopePoint p;
            p.value = value;

            const int n = static_cast<int>(m_buf.size());
            if (n >= 2) {
                double s = 0.0;
                for (double v : m_buf) s += v;
                const double mu = s / n;
                double ss = 0.0;
                for (double v : m_buf) ss += (v - mu) * (v - mu);
                const double sd = std::sqrt(ss / (n - 1));   // ddof = 1
                p.mean = mu;
                p.sd = sd;
                p.n = n;
                p.ready = (n >= m_window);
                if (!std::isnan(value) && sd > 0.0) p.z = (value - mu) / sd;
            }

            // Only real observations enter the baseline. Letting NaN through
            // would silently shorten the window; skipping it keeps the window a
            // window of BEATS THAT WERE MEASURED, which is what the SD means.
            if (!std::isnan(value)) {
                m_buf.push_back(value);
                while (static_cast<int>(m_buf.size()) > m_window) m_buf.pop_front();
            }
            return p;
        }

        void reset() { m_buf.clear(); }

    private:
        int m_window;
        std::deque<double> m_buf;
    };

    // ---------------------------------------------------------------------
    // Per-beat output for one segment, both window lengths
    // ---------------------------------------------------------------------

    struct SegmentEnvelopes {
        // Short window: PVC onset and conduction change.
        EnvelopePoint s_amplitude, s_centroid, s_wave0;
        // Long window: ischemic ST drift.
        EnvelopePoint l_mean, l_amplitude, l_centroid, l_wave0;
    };

    struct BeatEnvelopes {
        BeatMeasures     measures[3];   // indexed by Segment
        SegmentEnvelopes env[3];
    };

    // One tracker per segment, holding both window lengths for every measure it
    // reports. Feed beats in record order; each push returns that beat's point.
    class Tracker {
    public:
        Tracker(int shortWindow = kShortWindow, int longWindow = kLongWindow)
            : m_short(shortWindow), m_long(longWindow) {
            for (int s = 0; s < 3; ++s) {
                m_sAmp[s] = RollingEnvelope(shortWindow);
                m_sCent[s] = RollingEnvelope(shortWindow);
                m_sW0[s] = RollingEnvelope(shortWindow);
                m_lMean[s] = RollingEnvelope(longWindow);
                m_lAmp[s] = RollingEnvelope(longWindow);
                m_lCent[s] = RollingEnvelope(longWindow);
                m_lW0[s] = RollingEnvelope(longWindow);
            }
        }

        BeatEnvelopes push(const std::vector<double>& beat,
            const SegmentSpans& spans)
        {
            BeatEnvelopes out;
            for (int s = 0; s < 3; ++s) {
                const Segment seg = static_cast<Segment>(s);
                out.measures[s] = measure(beat, spans, seg);
                const BeatMeasures& m = out.measures[s];

                out.env[s].s_amplitude = m_sAmp[s].push(m.amplitude);
                out.env[s].s_centroid = m_sCent[s].push(m.centroid);
                out.env[s].s_wave0 = m_sW0[s].push(m.wave_frac[0]);

                out.env[s].l_mean = m_lMean[s].push(m.mean);
                out.env[s].l_amplitude = m_lAmp[s].push(m.amplitude);
                out.env[s].l_centroid = m_lCent[s].push(m.centroid);
                out.env[s].l_wave0 = m_lW0[s].push(m.wave_frac[0]);
            }
            return out;
        }

        int shortWindow() const { return m_short; }
        int longWindow()  const { return m_long; }

    private:
        int m_short, m_long;
        RollingEnvelope m_sAmp[3]{ RollingEnvelope(kShortWindow),
            RollingEnvelope(kShortWindow), RollingEnvelope(kShortWindow) };
        RollingEnvelope m_sCent[3]{ RollingEnvelope(kShortWindow),
            RollingEnvelope(kShortWindow), RollingEnvelope(kShortWindow) };
        RollingEnvelope m_sW0[3]{ RollingEnvelope(kShortWindow),
            RollingEnvelope(kShortWindow), RollingEnvelope(kShortWindow) };
        RollingEnvelope m_lMean[3]{ RollingEnvelope(kLongWindow),
            RollingEnvelope(kLongWindow), RollingEnvelope(kLongWindow) };
        RollingEnvelope m_lAmp[3]{ RollingEnvelope(kLongWindow),
            RollingEnvelope(kLongWindow), RollingEnvelope(kLongWindow) };
        RollingEnvelope m_lCent[3]{ RollingEnvelope(kLongWindow),
            RollingEnvelope(kLongWindow), RollingEnvelope(kLongWindow) };
        RollingEnvelope m_lW0[3]{ RollingEnvelope(kLongWindow),
            RollingEnvelope(kLongWindow), RollingEnvelope(kLongWindow) };
    };

    // ---------------------------------------------------------------------
    // Detectors, the two things the two window lengths exist to catch
    // ---------------------------------------------------------------------

    inline constexpr double kOnsetZ = 3.0;   // short-window departure
    inline constexpr double kDriftZ = 2.0;   // long-window departure

    // Short-window departure in QRS amplitude, spectral centroid or finest
    // wavelet scale. A PVC is wide and differently shaped, so its energy moves
    // down in frequency and its finest-scale detail drops -- the centroid and
    // wavelet terms catch a morphology change that left amplitude alone, which
    // is the case a purely amplitude-based detector misses.
    inline bool isOnsetEvent(const BeatEnvelopes& be, double zthr = kOnsetZ) {
        const SegmentEnvelopes& q = be.env[static_cast<int>(Segment::QRS)];
        auto hit = [&](const EnvelopePoint& p) {
            return p.ready && !std::isnan(p.z) && std::fabs(p.z) >= zthr;
            };
        return hit(q.s_amplitude) || hit(q.s_centroid) || hit(q.s_wave0);
    }

    // Long-window departure in T-wave level. ST drift is monotone and small per
    // beat, so it shows in the LONG mean and is invisible to the short window,
    // which follows it and treats it as the new normal within a few beats.
    inline bool isDriftEvent(const BeatEnvelopes& be, double zthr = kDriftZ) {
        const EnvelopePoint& t = be.env[static_cast<int>(Segment::T)].l_mean;
        return t.ready && !std::isnan(t.z) && std::fabs(t.z) >= zthr;
    }

}  // namespace envelopes
