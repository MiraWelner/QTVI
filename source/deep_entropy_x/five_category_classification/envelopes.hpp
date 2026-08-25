/**
 * @file   envelopes.hpp
 * @brief  Dynamic morphology envelopes over the P, QRS and T waves
 *         (Section 4.7): moving-average, power-spectral and wavelet envelopes
 *         on short (5 to 10 beat) and long (30 beat) windows. These feed the
 *         dynamic master template selection of Section 9.5.
 *
 *         SPEC FIDELITY. Section 4.7 states this task in prose and gives no
 *         code, no window arithmetic and no thresholds beyond the two window
 *         lengths. The choices below are therefore implementation, not
 *         transcription, and each is named as such where it is made. The short
 *         window defaults to 8 beats, the midpoint of the specified 5-to-10
 *         range.
 *
 *         RELATION TO morphology_envelope.hpp. Two different objects share the
 *         word "envelope" in Section 4.7 and they are not interchangeable:
 *
 *           morphology_envelope.hpp (4.7.1) is a corridor ACROSS SAMPLES --
 *           for one bin, the 2.5/97.5 percentile at each column of the beat.
 *           It answers "is this beat's shape inside the usual shape".
 *
 *           this file (4.7) is a corridor ACROSS BEATS -- for one scalar per
 *           beat (QRS area, P amplitude, T area), the trailing-window trend,
 *           dispersion, spectral content and wavelet detail energy. It answers
 *           "is this quantity drifting, alternating, or stepping".
 *
 *         Four implementation decisions:
 *
 *          1. WINDOWS ARE TRAILING, NOT CENTRED. A centred window cannot be
 *             evaluated until half of it is in the future, which makes the
 *             output unusable for the real-time progressive path in 4.7.3.
 *             Every value at beat i comes from beats (i-win, i]. Pass
 *             centred = true for the offline variant.
 *
 *          2. FREQUENCY IS IN CYCLES PER BEAT, NOT HERTZ. The series is
 *             indexed by beat and beat spacing is not uniform in time, so a
 *             periodogram of it has no Hz axis. Bands are drift [0, 0.1),
 *             mid [0.1, 0.4) and alternans [0.4, 0.5] cycles/beat. The
 *             alternans band is the useful one: 2:1 alternation -- bigeminy,
 *             QRS and T alternans -- lands at 0.5 cycles/beat and nowhere else.
 *
 *          3. THE WAVELET IS AN UNDECIMATED (a-trous) HAAR. Decimating halves
 *             the coefficient count per level, leaving three or four at level
 *             3 on a 30-beat window -- too few for an energy estimate. The
 *             undecimated form keeps one coefficient per beat at every level.
 *
 *          4. SHORT WINDOWS ARE MEAN-DETRENDED, LONG WINDOWS ARE
 *             LINE-DETRENDED. Over 30 beats an ischaemic ST ramp puts most of
 *             its energy in the low bins and would swamp the alternans band.
 *             The ramp is reported by the moving-average slope instead.
 *
 *         quantile_of is the quantile shared with morphology_envelope.hpp,
 *         which references it by name.
 *
 * @date   2026-08-24
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace envelopes {

    inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // Default window lengths. The spec gives the short window as a 5-to-10
    // beat range; 8 is the midpoint and is also a power of two, which the
    // Haar levels below want.
    inline constexpr int kShortWindow = 8;
    inline constexpr int kLongWindow = 30;
    inline constexpr int kMinForSpectrum = 8;   // below this a periodogram is noise
    inline constexpr int kWaveletLevels = 3;

    // ---------------------------------------------------------------------
    // Shared quantile
    // ---------------------------------------------------------------------
    // Linear-interpolated quantile. Takes the vector BY VALUE and sorts it,
    // dropping non-finite entries first, so callers can pass a raw column
    // without pre-cleaning it. Identical interpolation to
    // quantile_sorted()/quantile_nth() in morphology_envelope.hpp, so the
    // per-sample corridor and the per-beat corridor report comparable numbers
    // -- that agreement is the reason this function is shared rather than
    // reimplemented on each side.
    inline double quantile_of(std::vector<double> v, double q) {
        v.erase(std::remove_if(v.begin(), v.end(),
            [](double x) { return !std::isfinite(x); }), v.end());
        if (v.empty()) return kNaN;
        if (v.size() == 1) return v[0];
        q = std::clamp(q, 0.0, 1.0);
        std::sort(v.begin(), v.end());
        const double pos = q * static_cast<double>(v.size() - 1);
        const std::size_t i = static_cast<std::size_t>(std::floor(pos));
        const std::size_t k = std::min(i + 1, v.size() - 1);
        const double f = pos - static_cast<double>(i);
        return v[i] * (1.0 - f) + v[k] * f;
    }

    // ---------------------------------------------------------------------
    // Which per-beat scalar an envelope is tracking
    // ---------------------------------------------------------------------
    enum class Wave { P, QRS, T };

    inline const char* waveName(Wave w) {
        switch (w) {
        case Wave::P:   return "P";
        case Wave::QRS: return "QRS";
        default:        return "T";
        }
    }

    // One scalar per wave per beat. The caller decides what the scalar is
    // (peak amplitude, signed area, duration); the envelope machinery is
    // indifferent. NaN means "not measurable on this beat" and is carried
    // through as a gap rather than interpolated away -- an undetected P wave
    // is information, and filling it would manufacture a stable P series out
    // of a record that has none.
    struct BeatWaveScalars {
        double p = kNaN, qrs = kNaN, t = kNaN;
        double get(Wave w) const {
            return w == Wave::P ? p : (w == Wave::QRS ? qrs : t);
        }
    };

    inline std::vector<double> extractSeries(const std::vector<BeatWaveScalars>& s,
        Wave w)
    {
        std::vector<double> out(s.size(), kNaN);
        for (std::size_t i = 0; i < s.size(); ++i) out[i] = s[i].get(w);
        return out;
    }

    // ---------------------------------------------------------------------
    // Window extraction
    // ---------------------------------------------------------------------
    // Finite values of the trailing (or centred) window ending at beat i.
    // Returns fewer than `win` values near the record start, and callers
    // check `size()` against their own minimum rather than being handed a
    // silently short-sample statistic.
    inline std::vector<double> windowValues(const std::vector<double>& x, int i,
        int win, bool centred = false)
    {
        std::vector<double> w;
        if (win <= 0 || x.empty()) return w;
        const int n = static_cast<int>(x.size());
        int lo, hi;                                    // [lo, hi)
        if (centred) { lo = i - win / 2;  hi = lo + win; }
        else { hi = i + 1;      lo = hi - win; }
        lo = std::max(0, lo);
        hi = std::min(n, hi);
        w.reserve(static_cast<std::size_t>(std::max(0, hi - lo)));
        for (int k = lo; k < hi; ++k)
            if (std::isfinite(x[k])) w.push_back(x[k]);
        return w;
    }

    // ---------------------------------------------------------------------
    // 1. Moving-average envelope
    // ---------------------------------------------------------------------
    // Per beat: trailing mean, sd, and a 2.5/97.5 corridor over the window,
    // plus the z-score of the beat's own value against its window. The
    // z-score is the actionable output -- it is what a conduction change or a
    // sudden ectopic shows up in, and it is scale-free, so the same threshold
    // works on an amplitude series and an area series.
    struct MovingAverageEnvelope {
        std::vector<double> mean, sd, lo, hi, z;
        // Least-squares slope of the window, in series units per beat, and the
        // same slope expressed as total change across the window in units of
        // the window sd. `trend` is what a drift test should read: it is
        // scale-free and, unlike the spectral drift fraction, it survives the
        // line-detrending that the long-window spectrum applies.
        std::vector<double> slope, trend;
        std::vector<int>    nContrib;      // finite beats in each window
        int window = 0;
    };

    // HARD CEILING ON |z|. One outlier inside a window of n inflates that
    // window's own sd, and the arithmetic caps the achievable z-score at
    // (n-1)/sqrt(n) -- 2.47 for n = 8, 1.79 for n = 5. A step-detection
    // threshold above that ceiling can never fire, so the default zStep below
    // is set for the SHORT window, not for a textbook 3-sigma. Callers
    // changing the window length should re-check their threshold against this.
    inline double maxAttainableZ(int win) {
        return (win >= 2) ? (win - 1) / std::sqrt(static_cast<double>(win)) : kNaN;
    }

    inline MovingAverageEnvelope movingAverageEnvelope(const std::vector<double>& x,
        int win, bool centred = false,
        int minContrib = 3)
    {
        MovingAverageEnvelope e;
        e.window = win;
        const std::size_t n = x.size();
        e.mean.assign(n, kNaN); e.sd.assign(n, kNaN);
        e.lo.assign(n, kNaN);   e.hi.assign(n, kNaN);
        e.z.assign(n, kNaN);    e.nContrib.assign(n, 0);
        e.slope.assign(n, kNaN); e.trend.assign(n, kNaN);

        for (std::size_t i = 0; i < n; ++i) {
            const std::vector<double> w = windowValues(x, static_cast<int>(i), win, centred);
            e.nContrib[i] = static_cast<int>(w.size());
            if (static_cast<int>(w.size()) < minContrib) continue;

            double m = 0.0;
            for (double v : w) m += v;
            m /= static_cast<double>(w.size());
            double s2 = 0.0;
            for (double v : w) s2 += (v - m) * (v - m);
            const double s = std::sqrt(s2 / static_cast<double>(w.size() - 1));

            e.mean[i] = m;
            e.sd[i] = s;
            e.lo[i] = quantile_of(w, 0.025);
            e.hi[i] = quantile_of(w, 0.975);
            // A window of identical values has s == 0; z is then undefined
            // rather than infinite. Reporting inf would let one flat stretch
            // dominate every downstream max().
            if (s > 0.0 && std::isfinite(x[i])) e.z[i] = (x[i] - m) / s;

            // Least-squares slope over the surviving values, indexed by their
            // position within the window. Gaps shorten the window rather than
            // shifting the abscissa, which biases the slope slightly on a
            // gappy stretch -- acceptable, and visible through nContrib.
            const double xbar = 0.5 * (static_cast<double>(w.size()) - 1.0);
            double sxy = 0.0, sxx = 0.0;
            for (std::size_t k = 0; k < w.size(); ++k) {
                const double dx = static_cast<double>(k) - xbar;
                sxy += dx * (w[k] - m);
                sxx += dx * dx;
            }
            if (sxx > 0.0) {
                e.slope[i] = sxy / sxx;
                // Total change across the window, in sd units. A pure ramp has
                // sd proportional to its own span, so this saturates around
                // 3.4 for a noiseless ramp -- it is a drift DETECTOR, not a
                // drift magnitude. Read slope for magnitude.
                if (s > 0.0)
                    e.trend[i] = e.slope[i] * (static_cast<double>(w.size()) - 1.0) / s;
            }
        }
        return e;
    }

    // ---------------------------------------------------------------------
    // 2. Power-spectral envelope
    // ---------------------------------------------------------------------
    // Band powers of the beat-indexed series, in cycles/beat. `total` is the
    // detrended variance, so the three band fractions sum to 1 whenever the
    // window was long enough to transform.
    struct SpectralPoint {
        double total = kNaN;        // detrended power
        double drift = kNaN;        // [0, 0.1) cycles/beat
        double mid = kNaN;        // [0.1, 0.4)
        double alternans = kNaN;        // [0.4, 0.5] -- 2:1 alternation
        double driftFrac = kNaN, midFrac = kNaN, alternansFrac = kNaN;
        // Alternans power PER BIN divided by the mean power per bin everywhere
        // else. White noise sits at ~1 by construction; a real 2:1 alternation
        // is many times that. The FRACTION alone is not enough: on an 8-beat
        // window the alternans band is one of only four bins, so on pure noise
        // its share exceeds 0.35 about a third of the time. Measured on a
        // 70-beat noise stretch the fraction reached 0.91, which would read as
        // florid alternans. The SNR is what makes the test specific.
        double alternansSnr = kNaN;
        double dominantCpb = kNaN;  // frequency of the largest bin
    };

    struct SpectralEnvelope {
        std::vector<SpectralPoint> pt;
        int  window = 0;
        bool lineDetrended = false;
    };

    // Direct periodogram. N is at most 30 here, so an O(N^2) DFT is ~900
    // multiplies per beat per wave -- cheaper than setting up an FFT plan and
    // exact at every N, which matters because the window is short and N is
    // not a power of two.
    inline SpectralPoint spectrumOf(std::vector<double> w, bool lineDetrend) {
        SpectralPoint p;
        const int N = static_cast<int>(w.size());
        if (N < kMinForSpectrum) return p;

        double mean = 0.0;
        for (double v : w) mean += v;
        mean /= N;

        if (lineDetrend) {
            // Least-squares line on index, removed. Over a long window an
            // ischaemic ramp otherwise leaks across every low bin.
            const double xbar = 0.5 * (N - 1);
            double sxy = 0.0, sxx = 0.0;
            for (int i = 0; i < N; ++i) {
                sxy += (i - xbar) * (w[i] - mean);
                sxx += (i - xbar) * (i - xbar);
            }
            const double slope = (sxx > 0.0) ? sxy / sxx : 0.0;
            for (int i = 0; i < N; ++i) w[i] -= (mean + slope * (i - xbar));
        }
        else {
            for (int i = 0; i < N; ++i) w[i] -= mean;
        }

        const int kMax = N / 2;
        double drift = 0.0, mid = 0.0, alt = 0.0, best = -1.0, bestF = kNaN;
        int nAlt = 0, nOther = 0;
        for (int k = 1; k <= kMax; ++k) {
            double re = 0.0, im = 0.0;
            const double c = 2.0 * 3.14159265358979323846 * k / N;
            for (int i = 0; i < N; ++i) {
                re += w[i] * std::cos(c * i);
                im -= w[i] * std::sin(c * i);
            }
            const double pk = (re * re + im * im) / N;
            const double f = static_cast<double>(k) / N;   // cycles per beat
            if (f < 0.1) { drift += pk; ++nOther; }
            else if (f < 0.4) { mid += pk; ++nOther; }
            else { alt += pk;   ++nAlt; }
            if (pk > best) { best = pk; bestF = f; }
        }
        const double tot = drift + mid + alt;
        p.total = tot; p.drift = drift; p.mid = mid; p.alternans = alt;
        if (tot > 0.0) {
            p.driftFrac = drift / tot;
            p.midFrac = mid / tot;
            p.alternansFrac = alt / tot;
        }
        if (nAlt > 0 && nOther > 0) {
            const double otherPerBin = (drift + mid) / nOther;
            const double altPerBin = alt / nAlt;
            p.alternansSnr = (otherPerBin > 0.0) ? altPerBin / otherPerBin : kNaN;
        }
        p.dominantCpb = bestF;
        return p;
    }

    inline SpectralEnvelope spectralEnvelope(const std::vector<double>& x, int win,
        bool centred = false)
    {
        SpectralEnvelope e;
        e.window = win;
        e.lineDetrended = (win >= 16);       // long windows get the line removed
        e.pt.assign(x.size(), SpectralPoint{});
        for (std::size_t i = 0; i < x.size(); ++i)
            e.pt[i] = spectrumOf(windowValues(x, static_cast<int>(i), win, centred),
                e.lineDetrended);
        return e;
    }

    // ---------------------------------------------------------------------
    // 3. Wavelet envelope
    // ---------------------------------------------------------------------
    // Undecimated Haar. Level l detail responds to changes over 2^(l-1)
    // beats, so level 1 is beat-to-beat alternation, level 2 is a 2-to-4 beat
    // event (a couplet, a short run) and level 3 is a slower shift. Energy is
    // the mean square detail inside the window, which makes the three levels
    // directly comparable to each other.
    struct WaveletEnvelope {
        // energy[l][i]: level l+1 detail energy in the window ending at beat i
        std::vector<std::vector<double>> energy;
        std::vector<double> dominantLevel;   // 1-based level holding most energy
        int window = 0, levels = 0;
    };

    inline WaveletEnvelope waveletEnvelope(const std::vector<double>& x, int win,
        int levels = kWaveletLevels,
        bool centred = false)
    {
        WaveletEnvelope e;
        e.window = win;
        e.levels = std::max(1, levels);
        const int n = static_cast<int>(x.size());
        e.energy.assign(static_cast<std::size_t>(e.levels),
            std::vector<double>(static_cast<std::size_t>(n), kNaN));
        e.dominantLevel.assign(static_cast<std::size_t>(n), kNaN);
        if (n == 0) return e;

        // Gaps are held at the last finite value for the transform only. A
        // Haar difference across a NaN would propagate the NaN over 2^l beats
        // and blank the envelope well past the gap; holding produces a zero
        // detail there instead, which understates rather than erases. The
        // gaps themselves stay visible in MovingAverageEnvelope::nContrib.
        std::vector<double> s(x.begin(), x.end());
        double last = kNaN;
        for (int i = 0; i < n; ++i) {
            if (std::isfinite(s[i])) last = s[i];
            else                     s[i] = last;
        }
        for (int i = n - 1; i >= 0 && !std::isfinite(s[i]); --i) s[i] = 0.0;
        for (int i = 0; i < n; ++i) if (!std::isfinite(s[i])) s[i] = 0.0;

        std::vector<std::vector<double>> d(static_cast<std::size_t>(e.levels),
            std::vector<double>(static_cast<std::size_t>(n), 0.0));
        for (int l = 1; l <= e.levels; ++l) {
            const int lag = 1 << (l - 1);
            std::vector<double> sNext(static_cast<std::size_t>(n), 0.0);
            for (int i = 0; i < n; ++i) {
                const double prev = s[std::max(0, i - lag)];   // edge replication
                sNext[i] = 0.5 * (s[i] + prev);
                d[static_cast<std::size_t>(l - 1)][i] = s[i] - sNext[i];
            }
            s.swap(sNext);
        }

        for (int i = 0; i < n; ++i) {
            double bestE = -1.0; int bestL = -1;
            for (int l = 0; l < e.levels; ++l) {
                const std::vector<double> w = windowValues(d[static_cast<std::size_t>(l)],
                    i, win, centred);
                if (w.empty()) continue;
                double sum = 0.0;
                for (double v : w) sum += v * v;
                const double en = sum / static_cast<double>(w.size());
                e.energy[static_cast<std::size_t>(l)][static_cast<std::size_t>(i)] = en;
                if (en > bestE) { bestE = en; bestL = l + 1; }
            }
            if (bestL > 0) e.dominantLevel[static_cast<std::size_t>(i)] = bestL;
        }
        return e;
    }

    // ---------------------------------------------------------------------
    // Bundle: three envelope kinds x two window lengths, per wave
    // ---------------------------------------------------------------------
    struct WaveEnvelopes {
        Wave wave = Wave::QRS;
        MovingAverageEnvelope maShort, maLong;
        SpectralEnvelope      psShort, psLong;
        WaveletEnvelope       wvShort, wvLong;
    };

    struct DynamicEnvelopes {
        WaveEnvelopes p, qrs, t;
        int shortWindow = kShortWindow, longWindow = kLongWindow;
        int nBeats = 0;

        const WaveEnvelopes& of(Wave w) const {
            return w == Wave::P ? p : (w == Wave::QRS ? qrs : t);
        }
    };

    inline WaveEnvelopes buildWaveEnvelopes(const std::vector<double>& series, Wave w,
        int shortWin, int longWin, bool centred = false)
    {
        WaveEnvelopes e;
        e.wave = w;
        e.maShort = movingAverageEnvelope(series, shortWin, centred);
        e.maLong = movingAverageEnvelope(series, longWin, centred);
        e.psShort = spectralEnvelope(series, shortWin, centred);
        e.psLong = spectralEnvelope(series, longWin, centred);
        e.wvShort = waveletEnvelope(series, shortWin, kWaveletLevels, centred);
        e.wvLong = waveletEnvelope(series, longWin, kWaveletLevels, centred);
        return e;
    }

    inline DynamicEnvelopes buildDynamicEnvelopes(const std::vector<BeatWaveScalars>& s,
        int shortWin = kShortWindow,
        int longWin = kLongWindow,
        bool centred = false)
    {
        DynamicEnvelopes d;
        d.shortWindow = shortWin;
        d.longWindow = longWin;
        d.nBeats = static_cast<int>(s.size());
        d.p = buildWaveEnvelopes(extractSeries(s, Wave::P), Wave::P, shortWin, longWin, centred);
        d.qrs = buildWaveEnvelopes(extractSeries(s, Wave::QRS), Wave::QRS, shortWin, longWin, centred);
        d.t = buildWaveEnvelopes(extractSeries(s, Wave::T), Wave::T, shortWin, longWin, centred);
        return d;
    }

    // ---------------------------------------------------------------------
    // Readouts for Section 9.5 (dynamic master-template selection)
    // ---------------------------------------------------------------------
    // What 9.5 needs is not the whole envelope but the answer to three
    // questions per beat: is the morphology stable, is it stepping, is it
    // alternating. Those are the three fields below, and they are the only
    // thing the selector should have to read.
    struct EnvelopeVerdict {
        bool   stable = false;   // inside the long corridor, no step, no alternans
        bool   step = false;   // short-window excursion: PVC onset, conduction change
        bool   drifting = false;   // long-window ramp: ischaemic ST drift
        bool   alternating = false;   // 2:1 pattern in the short window
        double zShort = kNaN, zLong = kNaN, altFrac = kNaN, driftFrac = kNaN;
        double altSnr = kNaN, trendLong = kNaN;
    };

    struct VerdictThresholds {
        // 2.0, not 3.0: see maxAttainableZ. With the default 8-beat short
        // window nothing can exceed 2.47, so a 3-sigma rule would make `step`
        // dead code. 2.0 sits at 81% of the ceiling.
        double zStep = 2.0;
        // BOTH alternans gates must pass. See SpectralPoint::alternansSnr for
        // why the fraction on its own fires on white noise.
        double altFrac = 0.35;   // share of detrended power at 0.4-0.5 cyc/beat
        // 50, from the separation actually observed: on an 8-beat window a
        // white-noise stretch put this ratio at a median of 0.9 and a maximum
        // of 31, while an injected 2:1 alternation of 8 percent of the QRS
        // amplitude sat between 230 and 3600. 50 leaves margin on both sides.
        // A starting value: recalibrate against records with scored alternans
        // before reading anything clinical off it.
        double altSnr = 50.0;
        // Drift is read from the long window's TREND, not from its spectral
        // drift fraction. The long-window spectrum is line-detrended, which
        // removes a monotone ischaemic ramp by construction -- the ramp is
        // exactly what the detrending is there to take out so the alternans
        // band stays readable. driftFrac is still reported for diagnosis.
        double trendSd = 2.0;    // total change across the long window, in sd
        double driftFrac = 0.50;
        double minLongTotal = 0.0;    // ignore drift on a numerically flat series
    };

    inline std::vector<EnvelopeVerdict> verdicts(const WaveEnvelopes& e,
        const VerdictThresholds& th = {})
    {
        const std::size_t n = e.maShort.z.size();
        std::vector<EnvelopeVerdict> v(n);
        for (std::size_t i = 0; i < n; ++i) {
            EnvelopeVerdict& q = v[i];
            q.zShort = e.maShort.z[i];
            q.zLong = (i < e.maLong.z.size()) ? e.maLong.z[i] : kNaN;
            q.altFrac = (i < e.psShort.pt.size()) ? e.psShort.pt[i].alternansFrac : kNaN;
            q.altSnr = (i < e.psShort.pt.size()) ? e.psShort.pt[i].alternansSnr : kNaN;
            q.driftFrac = (i < e.psLong.pt.size()) ? e.psLong.pt[i].driftFrac : kNaN;

            q.trendLong = (i < e.maLong.trend.size()) ? e.maLong.trend[i] : kNaN;

            q.step = std::isfinite(q.zShort) && std::fabs(q.zShort) >= th.zStep;
            q.alternating = std::isfinite(q.altFrac) && q.altFrac >= th.altFrac
                && std::isfinite(q.altSnr) && q.altSnr >= th.altSnr;
            q.drifting = std::isfinite(q.trendLong)
                && std::fabs(q.trendLong) >= th.trendSd;
            // "Stable" requires positive evidence of stability, not merely the
            // absence of evidence: an unscorable beat is not a stable one.
            q.stable = std::isfinite(q.zShort) && !q.step && !q.alternating && !q.drifting;
        }
        return v;
    }

} // namespace envelopes