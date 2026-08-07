#pragma once
//
// beat_classify.hpp
//
// Per-beat classification into the five categories of spec 4.5, run BEFORE any
// template averaging. Only category 1 (bonafide PQRST) feeds feature extraction
// and the reference template; category 2 is retained with flags so arrhythmia
// burden can be reported; 3-5 are excluded.
//
// Three classifiers, all optional, combined by classify():
//
//   RULES      -- interpretable baseline. Prematurity, QRS correlation, depth
//                 and noise tests, each with a named threshold in Config.
//   GaussianNB -- diagonal-covariance naive Bayes over the same features, in
//                 log space. Returns a posterior per category, not a label.
//   OneClass   -- fits ONLY category 1 and scores novelty. Preferred when there
//                 is no annotation, because abundant normal beats make
//                 p(x | bonafide) estimable from a single record.
//
// The neural member of the spec's ensemble arrives through PosteriorFn, which is
// a std::function rather than a compiled-in dependency.
//
// ---------------------------------------------------------------------------
// ORDERING. Not the order the spec lists the categories in:
//
//   1. ARTIFACT_R first. An artifactual R-wave is a DETECTION error -- a beat
//      that was never there. It injects a spuriously short RR, which makes
//      is_premature() fire on both the false beat and its neighbour AND drags
//      the 10-beat rolling median down for the next ten beats. So these leave
//      the RR series before any rhythm test runs. Categories 4-5 are the
//      opposite case: real beats with corrupted signal, so their R positions
//      STAY in the RR series.
//
//   2. Noise evidence per beat, then the neighbourhood vote to confirm 4/5.
//
//   3. Rhythm (category 2): prematurity AND a morphology failure. Prematurity
//      alone cannot separate a PVC from a PAC -- see is_premature().
//
//   4. Correlation against the reference is circular (the reference is built
//      from category-1 beats, which needs the classification), so classify()
//      runs two passes -- see bootstrap_reference().
// ---------------------------------------------------------------------------

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <functional>
#include <memory>

#include "feature_marks.hpp"

namespace beat_classify {

    inline constexpr double kPi = 3.14159265358979323846;   // MSVC has no M_PI by default

    // ---- Categories (spec 4.5) --------------------------------------------
    enum class Category : uint8_t {
        UNCLASSIFIED = 0,
        BONAFIDE = 1,   // clean PQRST: the ONLY category that feeds templates
        ABNORMAL_RHYTHM = 2,   // ectopic / AF / PVC / VT -- excluded, retained + flagged
        ARTIFACT_R = 3,   // an R mark that isn't a beat (detection error)
        MACHINE_NOISE = 4,   // powerline, saturation, flatline, pacing spikes
        HUMAN_NOISE = 5,   // EMG, motion, wander, electrode pop
    };
    inline constexpr int kNumCategories = 6;   // includes UNCLASSIFIED at 0

    inline const char* category_name(Category c) {
        switch (c) {
        case Category::BONAFIDE:        return "bonafide";
        case Category::ABNORMAL_RHYTHM: return "abnormal_rhythm";
        case Category::ARTIFACT_R:      return "artifact_r";
        case Category::MACHINE_NOISE:   return "machine_noise";
        case Category::HUMAN_NOISE:     return "human_noise";
        default:                        return "unclassified";
        }
    }

    // Sub-reasons, so a rejection can be audited without re-running.
    enum ReasonBit : uint32_t {
        REASON_NONE = 0,
        REASON_PREMATURE = 1u << 0,
        REASON_LOW_ECG_CORR = 1u << 1,
        REASON_LOW_PPG_CORR = 1u << 2,
        REASON_WIDE_QRS = 1u << 3,
        REASON_NO_COMPENSATORY = 1u << 4,   // premature, no pause after -> PAC not PVC
        REASON_VOTE = 1u << 5,
        REASON_POWERLINE = 1u << 6,
        REASON_SATURATION = 1u << 7,
        REASON_FLATLINE = 1u << 8,
        REASON_RR_IMPLAUSIBLE = 1u << 9,
        REASON_EMG = 1u << 10,
        REASON_WANDER = 1u << 11,
        REASON_STEP = 1u << 12,
        REASON_NOVELTY = 1u << 13,   // one-class model
        REASON_AMPLITUDE = 1u << 14,
    };

    // ---- Per-beat result --------------------------------------------------
    struct BeatLabel {
        Category category = Category::UNCLASSIFIED;
        uint32_t reasons = REASON_NONE;
        double   ecg_corr = std::numeric_limits<double>::quiet_NaN();
        double   ppg_corr = std::numeric_limits<double>::quiet_NaN();
        double   rr_ms = std::numeric_limits<double>::quiet_NaN();
        double   novelty = std::numeric_limits<double>::quiet_NaN();  // mean z^2, ~1 = typical
        double   p_bonafide = std::numeric_limits<double>::quiet_NaN();
        bool accepted() const { return category == Category::BONAFIDE; }
    };

    struct CategoryCounts {
        int n[kNumCategories] = { 0, 0, 0, 0, 0, 0 };
        int total = 0;
        double percent(Category c) const {
            return (total > 0) ? 100.0 * n[static_cast<int>(c)] / total : 0.0;
        }
    };

    inline CategoryCounts tally(const std::vector<BeatLabel>& labels) {
        CategoryCounts out;
        for (const BeatLabel& l : labels) {
            const int i = static_cast<int>(l.category);
            if (i >= 0 && i < kNumCategories) ++out.n[i];
            ++out.total;
        }
        return out;
    }

    inline std::vector<char> accepted_mask(const std::vector<BeatLabel>& labels) {
        std::vector<char> m(labels.size(), 0);
        for (size_t i = 0; i < labels.size(); ++i) m[i] = labels[i].accepted() ? 1 : 0;
        return m;
    }

    // ---- Thresholds -------------------------------------------------------
    // Values from the spec are marked; the rest are defaults that need
    // calibration against annotated data.
    struct Config {
        // --- spec 4.6 ---
        double prematurity_frac = 0.80;   // RR(t) < 0.80 * rolling median
        int    prematurity_history = 10;
        double ecg_corr_min = 0.85;
        double ppg_corr_min = 0.80;
        int    vote_window = 8;
        int    vote_needed = 5;

        // --- not in the spec: placeholders ---
        double rr_min_ms = 250.0;   // shorter -> ARTIFACT_R, not a beat
        double rr_max_ms = 2500.0;
        double qrs_wide_ms = 120.0;   // premature + wide -> ventricular
        double compensatory_frac = 1.10;
        // Amplitude bounds for the artifact-R test. Deliberately wide: a beat
        // whose R amplitude is merely unusual is far more often noise-corrupted
        // than a false detection, and this test runs FIRST, so a tight bound
        // steals beats from categories 4-5. Only a grossly wrong amplitude, or
        // an impossible RR, declares "not a beat".
        double amp_ratio_min = 0.20;
        double amp_ratio_max = 4.00;
        double powerline_ratio_max = 0.25;   // 50/60 Hz power / total
        double flatline_sd_min = 1e-4;
        double saturation_frac_max = 0.02;   // fraction of samples at the rail
        double emg_ratio_max = 0.35;   // HF energy / total
        double wander_max = 0.50;   // |slope| across the beat, mV
        double step_max = 0.50;   // largest jump outside the QRS, mV
        double novelty_max = 12.0;   // mean z^2 above this -> not bonafide
        double p_bonafide_min = 0.90;   // posterior gate; see the note below
        bool   correlate_qrs_only = true;
        double powerline_hz = 60.0;
    };

    // =======================================================================
    // Small numeric helpers
    // =======================================================================

    inline double median_of(std::vector<double> v) {
        v.erase(std::remove_if(v.begin(), v.end(),
            [](double x) { return !std::isfinite(x); }), v.end());
        if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
        const size_t m = v.size() / 2;
        std::nth_element(v.begin(), v.begin() + m, v.end());
        const double hi = v[m];
        if (v.size() % 2) return hi;
        std::nth_element(v.begin(), v.begin() + m - 1, v.end());
        return 0.5 * (v[m - 1] + hi);
    }

    // Pearson r over [lo, hi] of two column-aligned beats. NaN if degenerate.
    inline double pearson(const std::vector<double>& a, const std::vector<double>& b,
        int lo, int hi)
    {
        const int n = static_cast<int>(std::min(a.size(), b.size()));
        lo = std::max(0, lo); hi = std::min(hi, n - 1);
        if (hi - lo < 3) return std::numeric_limits<double>::quiet_NaN();
        double sa = 0, sb = 0; int m = 0;
        for (int i = lo; i <= hi; ++i) {
            if (std::isnan(a[i]) || std::isnan(b[i])) continue;
            sa += a[i]; sb += b[i]; ++m;
        }
        if (m < 4) return std::numeric_limits<double>::quiet_NaN();
        const double ma = sa / m, mb = sb / m;
        double num = 0, da = 0, db = 0;
        for (int i = lo; i <= hi; ++i) {
            if (std::isnan(a[i]) || std::isnan(b[i])) continue;
            const double x = a[i] - ma, y = b[i] - mb;
            num += x * y; da += x * x; db += y * y;
        }
        if (da <= 0.0 || db <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        return num / std::sqrt(da * db);
    }

    inline double stddev_over(const std::vector<double>& v, int lo, int hi) {
        lo = std::max(0, lo); hi = std::min(hi, static_cast<int>(v.size()) - 1);
        if (hi - lo < 2) return 0.0;
        double s = 0.0; int n = 0;
        for (int i = lo; i <= hi; ++i) if (!std::isnan(v[i])) { s += v[i]; ++n; }
        if (n < 2) return 0.0;
        const double m = s / n;
        double q = 0.0;
        for (int i = lo; i <= hi; ++i) if (!std::isnan(v[i])) q += (v[i] - m) * (v[i] - m);
        return std::sqrt(q / (n - 1));
    }

    // Single-frequency power via Goertzel -- one bin of a DFT, O(n), no library.
    // Used for the powerline test; normalized by sample count so it compares
    // against variance directly.
    inline double goertzel_power(const std::vector<double>& v, int lo, int hi,
        double freq_hz, double fs)
    {
        lo = std::max(0, lo); hi = std::min(hi, static_cast<int>(v.size()) - 1);
        const int n = hi - lo + 1;
        if (n < 8 || fs <= 0.0 || freq_hz <= 0.0 || freq_hz >= 0.5 * fs) return 0.0;
        const double coeff = 2.0 * std::cos(2.0 * kPi * freq_hz / fs);
        double s1 = 0.0, s2 = 0.0;
        for (int i = lo; i <= hi; ++i) {
            const double x = std::isnan(v[i]) ? 0.0 : v[i];
            const double s0 = x + coeff * s1 - s2;
            s2 = s1; s1 = s0;
        }
        const double p = s1 * s1 + s2 * s2 - coeff * s1 * s2;
        return (p > 0.0) ? p / (static_cast<double>(n) * n) : 0.0;
    }

    // High-frequency energy fraction: mean squared SECOND difference over mean
    // squared deviation. Second differencing is a crude high-pass, so this rises
    // with EMG and quantization noise and is near zero on a smooth beat.
    inline double hf_energy_ratio(const std::vector<double>& v, int lo, int hi) {
        lo = std::max(1, lo); hi = std::min(hi, static_cast<int>(v.size()) - 2);
        if (hi - lo < 4) return 0.0;
        double hf = 0.0; int n = 0;
        for (int i = lo; i <= hi; ++i) {
            if (std::isnan(v[i - 1]) || std::isnan(v[i]) || std::isnan(v[i + 1])) continue;
            const double d2 = v[i + 1] - 2.0 * v[i] + v[i - 1];
            hf += d2 * d2; ++n;
        }
        if (n < 4) return 0.0;
        const double sd = stddev_over(v, lo, hi);
        if (sd <= 0.0) return 0.0;
        return std::sqrt(hf / n) / sd;
    }

    // =======================================================================
    // Rule primitives
    // =======================================================================

    // Prematurity, spec 4.6: RR(t) < 0.80 * median of the previous 10 RRs.
    //
    // Necessary but NOT sufficient for a PVC: a PAC is equally premature with a
    // normal narrow QRS. What separates them is morphology and the compensatory
    // pause, so classify() requires prematurity AND a morphology failure.
    inline bool is_premature(const std::vector<double>& rr_ms, int t, const Config& cfg = {})
    {
        const int h = cfg.prematurity_history;
        if (t < h || t >= static_cast<int>(rr_ms.size())) return false;
        const double med = median_of(std::vector<double>(rr_ms.begin() + (t - h), rr_ms.begin() + t));
        if (!std::isfinite(med) || med <= 0.0) return false;
        return rr_ms[t] < cfg.prematurity_frac * med;
    }

    // True if the beat AFTER t is lengthened -- the compensatory pause that
    // distinguishes a ventricular ectopic from an atrial one.
    inline bool has_compensatory_pause(const std::vector<double>& rr_ms, int t,
        const Config& cfg = {})
    {
        const int h = cfg.prematurity_history;
        const int n = static_cast<int>(rr_ms.size());
        if (t < h || t + 1 >= n) return false;
        const double med = median_of(std::vector<double>(rr_ms.begin() + (t - h), rr_ms.begin() + t));
        if (!std::isfinite(med) || med <= 0.0) return false;
        return rr_ms[t + 1] > cfg.compensatory_frac * med;
    }

    // Neighbourhood vote, spec 4.6. Three corrections against the snippet,
    // which read lo = max(0,t-4), hi = min(size,t+4), count over [lo,hi):
    //   -> asymmetric (t-4 .. t+3), counts t ITSELF, only 7 neighbours;
    //   -> at the edges the window shrinks below vote_needed, so the first and
    //      last four beats of a record can never be flagged whatever they hold.
    // This version is symmetric, excludes t, and scales the requirement when
    // the window is clipped so edge beats stay decidable.
    //
    // 5-of-8 is a SUSTAINED-corruption test. An isolated PVC has exactly one
    // flagged beat in its neighbourhood, counts 1, and passes as clean -- so
    // this gates categories 4-5, never category 2. The spec's acceptance test
    // asks the vote to flag PVCs, which it cannot do.
    inline bool vote_flag(const std::vector<char>& flag, int t, const Config& cfg = {})
    {
        const int n = static_cast<int>(flag.size());
        if (t < 0 || t >= n) return false;
        const int half = cfg.vote_window / 2;
        const int lo = std::max(0, t - half);
        const int hi = std::min(n - 1, t + half);
        int inspected = 0, count = 0;
        for (int i = lo; i <= hi; ++i) {
            if (i == t) continue;
            ++inspected;
            count += (flag[i] != 0) ? 1 : 0;
        }
        if (inspected <= 0) return false;
        const double need = static_cast<double>(cfg.vote_needed) * inspected
            / static_cast<double>(cfg.vote_window);
        return count >= static_cast<int>(std::ceil(need));
    }

    // =======================================================================
    // Reference template (breaks the correlation circularity)
    // =======================================================================

    // Per-sample median over the supplied beats -- the same aggregation
    // CreateEcgTemplates uses, so the reference is comparable to the template
    // the pipeline will eventually build. `mask` empty = use every beat.
    inline std::vector<double> median_template(
        const std::vector<std::vector<double>>& beats,
        const std::vector<char>& mask = {})
    {
        size_t len = 0, used = 0;
        for (size_t i = 0; i < beats.size(); ++i) {
            if (!mask.empty() && (i >= mask.size() || !mask[i])) continue;
            len = std::max(len, beats[i].size()); ++used;
        }
        if (used == 0 || len == 0) return {};
        std::vector<double> out(len, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> col; col.reserve(used);
        for (size_t j = 0; j < len; ++j) {
            col.clear();
            for (size_t i = 0; i < beats.size(); ++i) {
                if (!mask.empty() && (i >= mask.size() || !mask[i])) continue;
                if (j < beats[i].size() && !std::isnan(beats[i][j])) col.push_back(beats[i][j]);
            }
            if (!col.empty()) out[j] = median_of(col);
        }
        return out;
    }

    // First pass: no classification exists yet, so every beat contributes. The
    // median is order-statistic, so a minority of ectopics or noisy beats moves
    // it far less than a mean would -- good enough to correlate against, and the
    // second pass rebuilds from the survivors.
    inline std::vector<double> bootstrap_reference(
        const std::vector<std::vector<double>>& beats)
    {
        return median_template(beats);
    }

    inline std::vector<double> reference_from_accepted(
        const std::vector<std::vector<double>>& beats,
        const std::vector<BeatLabel>& labels)
    {
        return median_template(beats, accepted_mask(labels));
    }

    // Correlation window. Whole-beat r drops on ischaemic ST drift, which is a
    // category-1 beat we must KEEP -- detecting that drift is the point of the
    // long envelope in envelopes.hpp. So default to the QRS only, from the
    // reference's own landmarks.
    inline std::pair<int, int> correlation_window(const std::vector<double>& reference,
        int r_col, double fs, const Config& cfg = {})
    {
        const int n = static_cast<int>(reference.size());
        if (n < 8 || r_col < 0 || r_col >= n) return { 0, std::max(0, n - 1) };
        if (!cfg.correlate_qrs_only) return { 0, n - 1 };
        const double q = FeatureMarks::compute_q_onset(reference, fs, r_col);
        const double j = FeatureMarks::compute_j_point(reference, fs, r_col);
        int lo = (q >= 0.0) ? static_cast<int>(std::floor(q)) : r_col - static_cast<int>(0.05 * fs);
        int hi = (j >= 0.0) ? static_cast<int>(std::ceil(j)) : r_col + static_cast<int>(0.05 * fs);
        lo = std::clamp(lo, 0, n - 1); hi = std::clamp(hi, 0, n - 1);
        if (hi - lo < 4) { lo = std::max(0, r_col - 20); hi = std::min(n - 1, r_col + 20); }
        return { lo, hi };
    }

    // =======================================================================
    // Category 3: artifactual R waves. MUST run first -- see the ordering note.
    // =======================================================================
    // Returns the beat indices whose R mark isn't a beat, so the caller can drop
    // them from the RR series before any rhythm test sees it.
    inline std::vector<int> find_artifact_r(
        const std::vector<std::vector<double>>& beats,
        const std::vector<double>& rr_ms,
        int r_col, double /*fs*/, const Config& cfg = {})
    {
        std::vector<int> out;
        const int n = static_cast<int>(beats.size());

        // Running median R amplitude, for the amplitude test below.
        std::vector<double> amp(n, std::numeric_limits<double>::quiet_NaN());
        for (int i = 0; i < n; ++i)
            if (r_col >= 0 && r_col < static_cast<int>(beats[i].size()))
                amp[i] = std::abs(beats[i][r_col]);
        const double ampMed = median_of(amp);

        for (int i = 0; i < n; ++i) {
            bool bad = false;
            // (a) physiologically impossible interval -- decisive on its own
            if (i < static_cast<int>(rr_ms.size()) && std::isfinite(rr_ms[i])
                && (rr_ms[i] < cfg.rr_min_ms || rr_ms[i] > cfg.rr_max_ms)) bad = true;
            // (b) grossly wrong amplitude at the mark, i.e. no QRS-sized
            // deflection. Amplitude is only a PROXY for "is there a QRS here",
            // and motion can push a real beat outside any bound, so the window
            // is wide. A proper template-shape test would be the right fix.
            if (std::isfinite(ampMed) && ampMed > 0.0 && std::isfinite(amp[i])) {
                const double ratio = amp[i] / ampMed;
                if (ratio < cfg.amp_ratio_min || ratio > cfg.amp_ratio_max) bad = true;
            }
            // (c) degenerate beat
            if (beats[i].size() < 8) bad = true;
            if (bad) out.push_back(i);
        }
        return out;
    }

    // =======================================================================
    // Categories 4-5: per-beat noise evidence, BEFORE the vote is applied.
    // =======================================================================
    inline Category noise_category(const std::vector<double>& beat, double fs,
        uint32_t& reasons, const Config& cfg = {})
    {
        reasons = REASON_NONE;
        const int n = static_cast<int>(beat.size());
        if (n < 16 || fs <= 0.0) return Category::UNCLASSIFIED;

        const double sd = stddev_over(beat, 0, n - 1);

        // --- machine ---
        if (sd < cfg.flatline_sd_min) reasons |= REASON_FLATLINE;

        // Rail saturation: repeated samples at the extreme value.
        {
            double mx = -std::numeric_limits<double>::infinity();
            double mn = std::numeric_limits<double>::infinity();
            for (int i = 0; i < n; ++i)
                if (!std::isnan(beat[i])) { mx = std::max(mx, beat[i]); mn = std::min(mn, beat[i]); }
            if (std::isfinite(mx) && std::isfinite(mn)) {
                const double eps = 1e-9 + 1e-6 * (mx - mn);
                int atRail = 0;
                for (int i = 0; i < n; ++i)
                    if (!std::isnan(beat[i]) &&
                        (std::abs(beat[i] - mx) < eps || std::abs(beat[i] - mn) < eps)) ++atRail;
                if (static_cast<double>(atRail) / n > cfg.saturation_frac_max)
                    reasons |= REASON_SATURATION;
            }
        }

        // Powerline: single-bin power at 50/60 Hz against total variance.
        if (sd > 0.0) {
            const double pl = goertzel_power(beat, 0, n - 1, cfg.powerline_hz, fs);
            if (pl / (sd * sd) > cfg.powerline_ratio_max) reasons |= REASON_POWERLINE;
        }

        if (reasons & (REASON_FLATLINE | REASON_SATURATION | REASON_POWERLINE))
            return Category::MACHINE_NOISE;

        // --- human ---
        if (hf_energy_ratio(beat, 1, n - 2) > cfg.emg_ratio_max) reasons |= REASON_EMG;

        // Baseline wander: level difference between the two ends, measured on
        // short medians so a single spike can't trigger it.
        {
            const int w = std::max(3, n / 10);
            const double a = median_of(std::vector<double>(beat.begin(), beat.begin() + w));
            const double b = median_of(std::vector<double>(beat.end() - w, beat.end()));
            if (std::isfinite(a) && std::isfinite(b) && std::abs(b - a) > cfg.wander_max)
                reasons |= REASON_WANDER;
        }

        // Electrode pop: largest single-sample jump. The QRS itself is a large
        // legitimate jump, so this looks only at the outer thirds.
        {
            double mx = 0.0;
            const int a = 1, b = n / 3, c = 2 * n / 3, d = n - 1;
            for (int i = a; i <= b; ++i)
                if (!std::isnan(beat[i]) && !std::isnan(beat[i - 1]))
                    mx = std::max(mx, std::abs(beat[i] - beat[i - 1]));
            for (int i = c; i <= d; ++i)
                if (!std::isnan(beat[i]) && !std::isnan(beat[i - 1]))
                    mx = std::max(mx, std::abs(beat[i] - beat[i - 1]));
            if (mx > cfg.step_max) reasons |= REASON_STEP;
        }

        if (reasons & (REASON_EMG | REASON_WANDER | REASON_STEP))
            return Category::HUMAN_NOISE;
        return Category::UNCLASSIFIED;   // no noise evidence
    }

    // =======================================================================
    // Feature vector (shared by the rules and both statistical models)
    // =======================================================================
    enum FeatureIndex {
        F_RR_RATIO = 0,   // rr / median(prev 10)
        F_RR_NEXT_RATIO,      // rr(t+1) / median -- compensatory pause
        F_QRS_CORR,           // Pearson r vs reference, over the QRS
        F_QRS_WIDTH_MS,
        F_AMP_RATIO,          // |R| / running median |R|
        F_HF_RATIO,           // high-frequency energy fraction
        F_TP_SD,              // isoelectric-window sigma
        F_POWERLINE,          // 50/60 Hz power / variance
        F_COUNT
    };

    inline const char* feature_name(int j) {
        switch (j) {
        case F_RR_RATIO:      return "rr_ratio";
        case F_RR_NEXT_RATIO: return "rr_next_ratio";
        case F_QRS_CORR:      return "qrs_corr";
        case F_QRS_WIDTH_MS:  return "qrs_width_ms";
        case F_AMP_RATIO:     return "amp_ratio";
        case F_HF_RATIO:      return "hf_ratio";
        case F_TP_SD:         return "tp_sd";
        case F_POWERLINE:     return "powerline";
        default:              return "?";
        }
    }

    // Everything the extractor needs that isn't the beat itself.
    struct BeatContext {
        const std::vector<double>* reference = nullptr;
        const std::vector<double>* rr_ms = nullptr;
        int    index = 0;
        int    r_col = -1;
        double fs = 0.0;
        int    corr_lo = -1, corr_hi = -1;   // QRS window on the reference
        double amp_median = std::numeric_limits<double>::quiet_NaN();
        double qrs_width_ms = std::numeric_limits<double>::quiet_NaN();  // of the reference
    };

    inline std::vector<double> extract_features(const std::vector<double>& beat,
        const BeatContext& ctx, const Config& cfg = {})
    {
        std::vector<double> x(F_COUNT, std::numeric_limits<double>::quiet_NaN());
        const int n = static_cast<int>(beat.size());

        if (ctx.rr_ms) {
            const std::vector<double>& rr = *ctx.rr_ms;
            const int t = ctx.index, h = cfg.prematurity_history;
            if (t >= h && t < static_cast<int>(rr.size())) {
                const double med = median_of(std::vector<double>(rr.begin() + (t - h), rr.begin() + t));
                if (std::isfinite(med) && med > 0.0) {
                    x[F_RR_RATIO] = rr[t] / med;
                    if (t + 1 < static_cast<int>(rr.size())) x[F_RR_NEXT_RATIO] = rr[t + 1] / med;
                }
            }
        }

        if (ctx.reference && ctx.corr_lo >= 0)
            x[F_QRS_CORR] = pearson(beat, *ctx.reference, ctx.corr_lo, ctx.corr_hi);

        // Per-beat QRS width, from this beat's own landmarks.
        if (ctx.r_col >= 0 && ctx.fs > 0.0 && n > 8) {
            const double q = FeatureMarks::compute_q_onset(beat, ctx.fs, ctx.r_col);
            const double j = FeatureMarks::compute_j_point(beat, ctx.fs, ctx.r_col);
            if (q >= 0.0 && j > q) x[F_QRS_WIDTH_MS] = (j - q) * 1000.0 / ctx.fs;
        }

        if (ctx.r_col >= 0 && ctx.r_col < n && std::isfinite(ctx.amp_median) && ctx.amp_median > 0.0)
            x[F_AMP_RATIO] = std::abs(beat[ctx.r_col]) / ctx.amp_median;

        x[F_HF_RATIO] = hf_energy_ratio(beat, 1, n - 2);

        // TP window: from the reference's T-end to the end of the array.
        if (ctx.reference && ctx.r_col >= 0 && ctx.fs > 0.0) {
            const double j = FeatureMarks::compute_j_point(*ctx.reference, ctx.fs, ctx.r_col);
            const double tb = FeatureMarks::compute_t_begin(*ctx.reference, ctx.fs, ctx.r_col, j);
            const double te = FeatureMarks::compute_t_end(*ctx.reference, ctx.fs, ctx.r_col, tb);
            if (te >= 0.0) {
                const int lo = std::clamp(static_cast<int>(std::ceil(te)), 0, n - 1);
                x[F_TP_SD] = stddev_over(beat, lo, n - 1);
            }
        }

        const double sd = stddev_over(beat, 0, n - 1);
        if (sd > 0.0 && ctx.fs > 0.0)
            x[F_POWERLINE] = goertzel_power(beat, 0, n - 1, cfg.powerline_hz, ctx.fs) / (sd * sd);

        return x;
    }

    // =======================================================================
    // Gaussian naive Bayes, diagonal covariance, log space
    // =======================================================================
    // log P(c | x) = log P(c) + sum_j log p(x_j | c) - log Z
    //
    // "Naive" because the features are treated as conditionally independent
    // given the class. That is false here -- RR ratio and QRS width correlate in
    // a PVC -- but it mostly affects how SHARP the posteriors are rather than
    // their ranking, and a full covariance needs far more labelled data than a
    // single record provides.
    //
    // NaN features are skipped, so a beat missing its RR history still scores on
    // whatever else is available.
    struct GaussianNB {
        int nfeat = 0;
        std::vector<double> logPrior;   // [kNumCategories]
        std::vector<double> mean, var;  // [kNumCategories * nfeat]
        std::vector<int>    count;      // [kNumCategories]
        double varFloor = 1e-9;

        bool trained() const { return nfeat > 0 && !mean.empty(); }

        // One pass of sums and sums-of-squares per class.
        void fit(const std::vector<std::vector<double>>& X,
            const std::vector<Category>& y, int nfeatures)
        {
            nfeat = nfeatures;
            mean.assign(kNumCategories * nfeat, 0.0);
            var.assign(kNumCategories * nfeat, 0.0);
            count.assign(kNumCategories, 0);
            std::vector<int> per(kNumCategories * nfeat, 0);

            for (size_t i = 0; i < X.size() && i < y.size(); ++i) {
                const int c = static_cast<int>(y[i]);
                if (c <= 0 || c >= kNumCategories) continue;
                ++count[c];
                for (int j = 0; j < nfeat; ++j) {
                    const double v = X[i][j];
                    if (!std::isfinite(v)) continue;
                    mean[c * nfeat + j] += v;
                    var[c * nfeat + j] += v * v;
                    ++per[c * nfeat + j];
                }
            }
            for (int c = 0; c < kNumCategories; ++c) {
                for (int j = 0; j < nfeat; ++j) {
                    const int k = c * nfeat + j, m = per[k];
                    if (m < 2) { mean[k] = 0.0; var[k] = -1.0; continue; }   // -1 = unusable
                    const double mu = mean[k] / m;
                    mean[k] = mu;
                    var[k] = std::max(var[k] / m - mu * mu, varFloor);
                }
            }
            // Priors from the class frequencies.
            int total = 0;
            for (int c = 0; c < kNumCategories; ++c) total += count[c];
            logPrior.assign(kNumCategories, -std::numeric_limits<double>::infinity());
            for (int c = 1; c < kNumCategories; ++c)
                if (count[c] > 0 && total > 0)
                    logPrior[c] = std::log(static_cast<double>(count[c]) / total);
        }

        std::vector<double> posterior(const std::vector<double>& x) const {
            std::vector<double> lp(kNumCategories, -std::numeric_limits<double>::infinity());
            if (!trained()) return std::vector<double>(kNumCategories, 0.0);
            for (int c = 1; c < kNumCategories; ++c) {
                if (!std::isfinite(logPrior[c])) continue;
                double s = logPrior[c]; int used = 0;
                for (int j = 0; j < nfeat && j < static_cast<int>(x.size()); ++j) {
                    const double v = var[c * nfeat + j];
                    if (v < 0.0 || !std::isfinite(x[j])) continue;   // unusable or missing
                    const double d = x[j] - mean[c * nfeat + j];
                    s += -0.5 * (std::log(2.0 * kPi * v) + d * d / v);
                    ++used;
                }
                if (used > 0) lp[c] = s;
            }
            // log-sum-exp normalise
            double m = -std::numeric_limits<double>::infinity();
            for (double v : lp) m = std::max(m, v);
            std::vector<double> p(kNumCategories, 0.0);
            if (!std::isfinite(m)) return p;
            double z = 0.0;
            for (double v : lp) if (std::isfinite(v)) z += std::exp(v - m);
            if (z <= 0.0) return p;
            for (int c = 0; c < kNumCategories; ++c)
                if (std::isfinite(lp[c])) p[c] = std::exp(lp[c] - m) / z;
            return p;
        }

        // PVC burden varies 0-30% between subjects, so a fixed prior is wrong
        // for nearly everyone. Refit the priors to THIS record by EM: posteriors
        // with the current priors, then set P(c) to the mean posterior. Two or
        // three rounds is enough, and it is the single change that most improves
        // calibration.
        void refit_priors_em(const std::vector<std::vector<double>>& X, int rounds = 3) {
            if (!trained() || X.empty()) return;
            for (int r = 0; r < rounds; ++r) {
                std::vector<double> acc(kNumCategories, 0.0);
                int n = 0;
                for (const std::vector<double>& x : X) {
                    const std::vector<double> p = posterior(x);
                    for (int c = 0; c < kNumCategories; ++c) acc[c] += p[c];
                    ++n;
                }
                if (n == 0) return;
                for (int c = 1; c < kNumCategories; ++c) {
                    const double pi = acc[c] / n;
                    logPrior[c] = (pi > 1e-12) ? std::log(pi)
                        : -std::numeric_limits<double>::infinity();
                }
            }
        }
    };

    // =======================================================================
    // One-class model: fit ONLY category 1
    // =======================================================================
    // Preferred when there is no annotation. Normal beats are abundant, so
    // p(x | bonafide) is estimable from a single record with no labelling, and
    // it answers the question the pipeline actually asks -- may this beat enter
    // the reference template.
    //
    // novelty() is the mean squared z-score across usable features: ~1 for a
    // typical beat, growing with distance. It is NOT a probability; threshold it
    // directly (Config::novelty_max) rather than pretending otherwise.
    struct OneClassNormal {
        int nfeat = 0;
        std::vector<double> mean, var;
        std::vector<int>    per;
        double varFloor = 1e-9;

        bool trained() const { return nfeat > 0 && !mean.empty(); }

        // Relative floor. An ABSOLUTE floor is useless here: features like
        // rr_ratio sit at 1.0 +- 0.001 among accepted beats, so an unfloored
        // z-score divides by ~0 and novelty explodes to 1e8. Flooring at a
        // fraction of the feature's own scale keeps z interpretable.
        double relFloor = 1e-3;

        void fit(const std::vector<std::vector<double>>& X,
            const std::vector<char>& mask, int nfeatures)
        {
            nfeat = nfeatures;
            mean.assign(nfeat, 0.0); var.assign(nfeat, 0.0); per.assign(nfeat, 0);
            for (size_t i = 0; i < X.size(); ++i) {
                if (!mask.empty() && (i >= mask.size() || !mask[i])) continue;
                for (int j = 0; j < nfeat; ++j) {
                    const double v = X[i][j];
                    if (!std::isfinite(v)) continue;
                    mean[j] += v; var[j] += v * v; ++per[j];
                }
            }
            for (int j = 0; j < nfeat; ++j) {
                if (per[j] < 4) { mean[j] = 0.0; var[j] = -1.0; continue; }
                const double mu = mean[j] / per[j];
                mean[j] = mu;
                const double scale = relFloor * std::max(1.0, std::abs(mu));
                var[j] = std::max(var[j] / per[j] - mu * mu,
                    std::max(varFloor, scale * scale));
            }
        }

        // Per-feature z^2 is CAPPED, so novelty is bounded by zCap^2 and one
        // near-degenerate feature cannot swamp the score. Without this, a
        // feature that happens to be near-constant among accepted beats (rr
        // ratio on a metronomic record, say) drives novelty into the millions
        // and any threshold on it becomes meaningless.
        double zCap = 10.0;

        double novelty(const std::vector<double>& x) const {
            if (!trained()) return std::numeric_limits<double>::quiet_NaN();
            double s = 0.0; int used = 0;
            for (int j = 0; j < nfeat && j < static_cast<int>(x.size()); ++j) {
                if (var[j] < 0.0 || !std::isfinite(x[j])) continue;
                double z = (x[j] - mean[j]) / std::sqrt(var[j]);
                z = std::clamp(z, -zCap, zCap);
                s += z * z; ++used;
            }
            return (used > 0) ? s / used : std::numeric_limits<double>::quiet_NaN();
        }
    };

    // =======================================================================
    // Ensemble hook
    // =======================================================================
    // Spec 4.5 wants Bayesian + nonlinear + neural voting. GaussianNB above is
    // the Bayesian member; this injects a second, so a neural model can arrive
    // later without touching this header. Return an empty vector for "no
    // opinion" and classify() falls back to the rules.
    using PosteriorFn = std::function<std::vector<double>(
        const std::vector<double>& features)>;

    // =======================================================================
    // Top level
    // =======================================================================
    struct ClassifyResult {
        std::vector<BeatLabel>            labels;
        std::vector<std::vector<double>>  features;   // parallel to labels
        std::vector<double>               reference;  // final, from accepted beats
        CategoryCounts                    counts;
        OneClassNormal                    oneClass;   // fitted on the accepted set
    };

    // One label per beat, in the order documented at the top of this file.
    // `rr_ms` may be shorter than `beats`; missing entries just disable the
    // rhythm tests for those beats.
    inline ClassifyResult classify(
        const std::vector<std::vector<double>>& beats,
        const std::vector<double>& rr_ms_in,
        int r_col, double fs,
        const Config& cfg = {},
        PosteriorFn posterior = nullptr)
    {
        ClassifyResult res;
        const int n = static_cast<int>(beats.size());
        res.labels.assign(n, BeatLabel{});
        if (n == 0) return res;

        // ---- step 1: category 3, and remove those RRs -----------------------
        std::vector<double> rr = rr_ms_in;
        rr.resize(n, std::numeric_limits<double>::quiet_NaN());
        for (int i = 0; i < n; ++i) res.labels[i].rr_ms = rr[i];

        for (int i : find_artifact_r(beats, rr, r_col, fs, cfg)) {
            res.labels[i].category = Category::ARTIFACT_R;
            res.labels[i].reasons |= REASON_RR_IMPLAUSIBLE | REASON_AMPLITUDE;
            rr[i] = std::numeric_limits<double>::quiet_NaN();   // out of the RR series
        }

        // ---- step 2: noise evidence, then the vote --------------------------
        std::vector<Category> noiseCat(n, Category::UNCLASSIFIED);
        std::vector<uint32_t> noiseWhy(n, REASON_NONE);
        std::vector<char>     noisy(n, 0);
        for (int i = 0; i < n; ++i) {
            if (res.labels[i].category != Category::UNCLASSIFIED) continue;
            noiseCat[i] = noise_category(beats[i], fs, noiseWhy[i], cfg);
            noisy[i] = (noiseCat[i] != Category::UNCLASSIFIED) ? 1 : 0;
        }
        for (int i = 0; i < n; ++i) {
            if (res.labels[i].category != Category::UNCLASSIFIED) continue;
            if (!noisy[i]) continue;
            // A single noisy beat is accepted as noise on its own evidence; the
            // vote only ADDS beats whose neighbours are corrupted, catching the
            // middle of a burst whose own features happen to look clean.
            res.labels[i].category = noiseCat[i];
            res.labels[i].reasons |= noiseWhy[i];
        }
        for (int i = 0; i < n; ++i) {
            if (res.labels[i].category != Category::UNCLASSIFIED) continue;
            if (vote_flag(noisy, i, cfg)) {
                res.labels[i].category = Category::HUMAN_NOISE;
                res.labels[i].reasons |= REASON_VOTE;
            }
        }

        // ---- steps 3-4: two passes over reference + rhythm ------------------
        std::vector<double> reference = bootstrap_reference(beats);

        // Running median R amplitude, for F_AMP_RATIO.
        std::vector<double> amps;
        amps.reserve(n);
        for (int i = 0; i < n; ++i)
            if (r_col >= 0 && r_col < static_cast<int>(beats[i].size()))
                amps.push_back(std::abs(beats[i][r_col]));
        const double ampMed = median_of(amps);

        for (int pass = 0; pass < 2; ++pass) {
            const std::pair<int, int> cw = correlation_window(reference, r_col, fs, cfg);

            double refQrsMs = std::numeric_limits<double>::quiet_NaN();
            if (!reference.empty() && r_col >= 0 && fs > 0.0) {
                const double q = FeatureMarks::compute_q_onset(reference, fs, r_col);
                const double j = FeatureMarks::compute_j_point(reference, fs, r_col);
                if (q >= 0.0 && j > q) refQrsMs = (j - q) * 1000.0 / fs;
            }

            res.features.assign(n, std::vector<double>(F_COUNT,
                std::numeric_limits<double>::quiet_NaN()));

            for (int i = 0; i < n; ++i) {
                BeatContext ctx;
                ctx.reference = &reference; ctx.rr_ms = &rr; ctx.index = i;
                ctx.r_col = r_col; ctx.fs = fs;
                ctx.corr_lo = cw.first; ctx.corr_hi = cw.second;
                ctx.amp_median = ampMed; ctx.qrs_width_ms = refQrsMs;
                res.features[i] = extract_features(beats[i], ctx, cfg);

                BeatLabel& L = res.labels[i];
                L.ecg_corr = res.features[i][F_QRS_CORR];
                if (L.category == Category::ARTIFACT_R
                    || L.category == Category::MACHINE_NOISE
                    || L.category == Category::HUMAN_NOISE) continue;

                // Rhythm: premature AND a morphology failure. Prematurity alone
                // would flag every atrial ectopic as ventricular.
                const bool prem = is_premature(rr, i, cfg);
                const bool lowCorr = std::isfinite(L.ecg_corr) && L.ecg_corr < cfg.ecg_corr_min;
                const double wid = res.features[i][F_QRS_WIDTH_MS];
                const bool wide = std::isfinite(wid) && wid > cfg.qrs_wide_ms;

                uint32_t why = REASON_NONE;
                if (prem)     why |= REASON_PREMATURE;
                if (lowCorr)  why |= REASON_LOW_ECG_CORR;
                if (wide)     why |= REASON_WIDE_QRS;
                if (prem && !has_compensatory_pause(rr, i, cfg)) why |= REASON_NO_COMPENSATORY;

                if (prem && (lowCorr || wide)) {
                    L.category = Category::ABNORMAL_RHYTHM;
                    L.reasons |= why;
                }
                else if (lowCorr) {
                    // Not premature but morphologically wrong: keep it out of
                    // the reference without claiming it is ectopic.
                    L.category = Category::HUMAN_NOISE;
                    L.reasons |= REASON_LOW_ECG_CORR;
                }
                else {
                    L.category = Category::BONAFIDE;
                    L.reasons |= why;   // e.g. premature-but-normal (PAC) stays flagged
                }
            }

            if (pass == 0) {
                std::vector<double> rebuilt = reference_from_accepted(beats, res.labels);
                if (!rebuilt.empty()) reference = std::move(rebuilt);
                // Reset everything the rhythm stage decided; noise and artifact
                // labels from steps 1-2 stand.
                for (int i = 0; i < n; ++i)
                    if (res.labels[i].category == Category::BONAFIDE
                        || res.labels[i].category == Category::ABNORMAL_RHYTHM)
                        res.labels[i].category = Category::UNCLASSIFIED;
            }
        }

        // ---- step 5: one-class novelty on the accepted set ------------------
        res.oneClass.fit(res.features, accepted_mask(res.labels), F_COUNT);
        for (int i = 0; i < n; ++i) {
            res.labels[i].novelty = res.oneClass.novelty(res.features[i]);
            if (res.labels[i].category == Category::BONAFIDE
                && std::isfinite(res.labels[i].novelty)
                && res.labels[i].novelty > cfg.novelty_max) {
                res.labels[i].category = Category::HUMAN_NOISE;
                res.labels[i].reasons |= REASON_NOVELTY;
            }
        }

        // ---- step 6: optional injected posterior ---------------------------
        // The decision rule is NOT argmax. Costs are wildly asymmetric: there
        // are thousands of beats, so discarding a good one is nearly free, while
        // admitting one ectopic contaminates the reference every other beat is
        // compared against. So this only ever DEMOTES a beat the rules accepted.
        if (posterior) {
            for (int i = 0; i < n; ++i) {
                const std::vector<double> p = posterior(res.features[i]);
                if (p.size() != static_cast<size_t>(kNumCategories)) continue;
                res.labels[i].p_bonafide = p[static_cast<int>(Category::BONAFIDE)];
                if (res.labels[i].category == Category::BONAFIDE
                    && res.labels[i].p_bonafide < cfg.p_bonafide_min) {
                    // Demote to whichever abnormal class the model prefers.
                    int best = 0; double bv = -1.0;
                    for (int c = 2; c < kNumCategories; ++c) if (p[c] > bv) { bv = p[c]; best = c; }
                    if (best >= 2) res.labels[i].category = static_cast<Category>(best);
                }
            }
        }

        res.reference = reference_from_accepted(beats, res.labels);
        if (res.reference.empty()) res.reference = reference;
        res.counts = tally(res.labels);
        return res;
    }

    // Convenience: fit a GaussianNB from a completed rule pass and hand back a
    // PosteriorFn ready for a second classify() call. Self-training -- it
    // inherits every bias the rules have, so treat its output as a refinement of
    // the rules rather than an independent opinion.
    inline PosteriorFn make_nb_posterior(const ClassifyResult& res, int em_rounds = 3) {
        std::vector<Category> y;
        y.reserve(res.labels.size());
        for (const BeatLabel& L : res.labels) y.push_back(L.category);
        auto model = std::make_shared<GaussianNB>();
        model->fit(res.features, y, F_COUNT);
        model->refit_priors_em(res.features, em_rounds);
        return [model](const std::vector<double>& x) { return model->posterior(x); };
    }

}   // namespace beat_classify