#pragma once
//
// beat_classify.hpp
//
// Per-beat classification into the five categories of spec 4.5, run BEFORE
// any template averaging. Only category 1 (bonafide PQRST) feeds feature
// extraction and the reference template; category 2 is retained with flags so
// arrhythmia burden can be reported; 3-5 are excluded outright.
//
// STATUS: stubs. Every function below compiles and returns a defined
// "nothing detected" value. Bodies marked TODO.
//
// ---------------------------------------------------------------------------
// Ordering matters, and it is not the order the spec lists the categories in:
//
//   1. ARTIFACT_R first. An artifactual R-wave is a DETECTION error -- a beat
//      that was never there. It injects a spuriously short RR into the series,
//      which makes is_premature() fire on both the false beat and its
//      neighbour AND drags the 10-beat rolling median down for the next ten
//      beats. So these have to be pulled out of the RR series before any
//      rhythm test runs. Categories 4-5 are the opposite case: real beats with
//      corrupted signal, so their R positions STAY in the RR series.
//
//   2. Rhythm (category 2) next: prematurity AND morphology together. See
//      is_premature() for why prematurity alone cannot identify a PVC.
//
//   3. Noise (categories 4-5) last, via the neighbourhood vote. The vote is a
//      SUSTAINED-corruption test and belongs here, not on category 2 -- see
//      vote_flag().
//
//   4. Correlation against the reference template is circular (the template is
//      built from category-1 beats, which requires the classification), so
//      classify() runs two passes. See bootstrap_reference().
// ---------------------------------------------------------------------------

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <limits>
#include <functional>

namespace beat_classify {

    // ---- Categories (spec 4.5) --------------------------------------------
    enum class Category : uint8_t {
        UNCLASSIFIED = 0,   // not yet assigned
        BONAFIDE = 1,   // clean PQRST: the ONLY category that feeds templates
        ABNORMAL_RHYTHM = 2,   // ectopic / AF / PVC / VT -- excluded, but retained + flagged
        ARTIFACT_R = 3,   // R-wave that isn't a beat (detection error)
        MACHINE_NOISE = 4,   // powerline, rail saturation, lead-off, pacing spikes, quantization
        HUMAN_NOISE = 5,   // EMG, motion, baseline wander, electrode pop
    };

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

    // Sub-reasons, kept alongside the category so a rejection can be audited
    // without re-running the classifier. Bitmask: more than one can apply.
    enum ReasonBit : uint32_t {
        REASON_NONE = 0,
        REASON_PREMATURE = 1u << 0,   // RR(t) < prematurity_frac * rolling median
        REASON_LOW_ECG_CORR = 1u << 1,   // r < ecg_corr_min against the reference
        REASON_LOW_PPG_CORR = 1u << 2,   // r < ppg_corr_min
        REASON_WIDE_QRS = 1u << 3,   // QRS duration above qrs_wide_ms
        REASON_NO_COMPENSATORY = 1u << 4,   // premature but no following pause (-> PAC, not PVC)
        REASON_VOTE = 1u << 5,   // neighbourhood vote (sustained corruption)
        REASON_POWERLINE = 1u << 6,
        REASON_SATURATION = 1u << 7,
        REASON_FLATLINE = 1u << 8,
        REASON_RR_IMPLAUSIBLE = 1u << 9,   // RR too short to be physiological -> artifact R
    };

    // ---- Per-beat result --------------------------------------------------
    struct BeatLabel {
        Category category = Category::UNCLASSIFIED;
        uint32_t reasons = REASON_NONE;
        double   ecg_corr = std::numeric_limits<double>::quiet_NaN();   // vs reference
        double   ppg_corr = std::numeric_limits<double>::quiet_NaN();
        double   rr_ms = std::numeric_limits<double>::quiet_NaN();
        double   confidence = std::numeric_limits<double>::quiet_NaN();  // 0..1, for the ensemble hook
    };

    // Per-bin summary -- the acceptance test's "per-bin category percentages".
    struct CategoryCounts {
        int n[6] = { 0, 0, 0, 0, 0, 0 };   // indexed by Category
        int total = 0;
        double percent(Category c) const {
            return (total > 0) ? 100.0 * n[static_cast<int>(c)] / total : 0.0;
        }
    };

    inline CategoryCounts tally(const std::vector<BeatLabel>& labels) {
        CategoryCounts out;
        for (const BeatLabel& l : labels) {
            const int i = static_cast<int>(l.category);
            if (i >= 0 && i < 6) ++out.n[i];
            ++out.total;
        }
        return out;
    }

    // ---- Thresholds -------------------------------------------------------
    // All in one struct so a caller can tune without touching the algorithms,
    // and so the values that came from the spec are visibly separate from the
    // ones that had to be invented.
    struct Config {
        // --- from spec 4.6 ---
        double prematurity_frac = 0.80;   // RR(t) < 0.80 * rolling median
        int    prematurity_history = 10;     // beats in the rolling median
        double ecg_corr_min = 0.85;   // exclude below this
        double ppg_corr_min = 0.80;
        int    vote_window = 8;      // neighbours inspected
        int    vote_needed = 5;      // 5-of-8

        // --- NOT in the spec; defaults are placeholders ---
        double rr_min_ms = 250.0;   // shorter than this -> ARTIFACT_R, not a real beat
        double qrs_wide_ms = 120.0;   // above this, a premature beat is ventricular not atrial
        double compensatory_frac = 1.10;   // RR(t+1) > 1.10 * median => compensatory pause
        bool   correlate_qrs_only = true;    // see note in correlate()
    };

    // ---- Primitives -------------------------------------------------------

    // Pearson r over [lo, hi] of two equal-length, column-aligned beats.
    // NaN if either side is flat or the span is degenerate.
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

    // Prematurity, per spec 4.6: RR(t) < 0.80 * median of the previous 10 RRs.
    //
    // NOTE this is necessary but NOT sufficient to call something a PVC. A PAC
    // is equally premature with a perfectly normal narrow QRS; what separates
    // them is morphology (wide, aberrant) and the compensatory pause after.
    // classify() therefore requires prematurity AND a morphology failure --
    // this predicate on its own would flag every atrial ectopic as ventricular.
    inline bool is_premature(const std::vector<double>& rr_ms, int t, const Config& cfg = {})
    {
        const int h = cfg.prematurity_history;
        if (t < h || t >= static_cast<int>(rr_ms.size())) return false;
        std::vector<double> w(rr_ms.begin() + (t - h), rr_ms.begin() + t);
        w.erase(std::remove_if(w.begin(), w.end(),
            [](double x) { return !std::isfinite(x); }), w.end());
        if (w.empty()) return false;
        std::nth_element(w.begin(), w.begin() + w.size() / 2, w.end());
        const double med = w[w.size() / 2];
        return rr_ms[t] < cfg.prematurity_frac * med;
    }

    // Neighbourhood vote, per spec 4.6: flag beat t if >= vote_needed of the
    // surrounding vote_window beats are already flagged.
    //
    // Three corrections against the spec's snippet, which read
    //   lo = max(0,t-4), hi = min(size,t+4), count over [lo,hi)
    //   -> asymmetric (t-4 .. t+3), counts t ITSELF, and only 7 neighbours;
    //   -> at the edges the window shrinks below vote_needed, so the first and
    //      last four beats of a record can never be flagged whatever they hold.
    // This version is symmetric, excludes t, and scales the requirement when
    // the window is clipped so edge beats stay testable.
    //
    // Also: 5-of-8 is a SUSTAINED-corruption test. An isolated PVC has exactly
    // one flagged beat in its neighbourhood, counts 1, and passes as clean --
    // so this must gate categories 4-5 (motion/EMG lasting seconds), never
    // category 2. The spec's acceptance test asks the vote to flag PVCs, which
    // it cannot do.
    inline bool vote_flag(const std::vector<char>& flag, int t, const Config& cfg = {})
    {
        const int n = static_cast<int>(flag.size());
        if (t < 0 || t >= n) return false;
        const int half = cfg.vote_window / 2;
        const int lo = std::max(0, t - half);
        const int hi = std::min(n - 1, t + half);
        int inspected = 0, count = 0;
        for (int i = lo; i <= hi; ++i) {
            if (i == t) continue;                  // don't let t vote for itself
            ++inspected;
            count += (flag[i] != 0) ? 1 : 0;
        }
        if (inspected <= 0) return false;
        // Clipped window: keep the same required FRACTION rather than the raw
        // count, so beats near the record edges are still decidable.
        const double need = static_cast<double>(cfg.vote_needed) * inspected
            / static_cast<double>(cfg.vote_window);
        return count >= static_cast<int>(std::ceil(need));
    }

    // ---- Stubs ------------------------------------------------------------

    // Category 3. Must run FIRST: returns the indices of R marks that aren't
    // beats, so the caller can drop them from the RR series before any rhythm
    // test sees it.
    // TODO: RR below cfg.rr_min_ms; no QRS-shaped deflection at the mark;
    //       amplitude wildly off the running median; two marks inside one QRS.
    inline std::vector<int> find_artifact_r(const std::vector<std::vector<double>>& /*beats*/,
        const std::vector<double>& /*rr_ms*/,
        int /*r_col*/, double /*fs*/,
        const Config & /*cfg*/ = {})
    {
        return {};
    }

    // Categories 4-5. Per-beat noise evidence, before the vote is applied.
    // TODO machine: 50/60 Hz band power vs total, rail saturation (clipped
    //       runs at the ADC limit), flatline (zero variance), pacing spikes,
    //       quantization steps. notch_filter() and the existing noise-marking
    //       GUI in this repo already have prior art for the powerline part.
    // TODO human: EMG band power, |baseline wander| slope, step discontinuity
    //       (electrode pop), sample-to-sample jerk.
    // The spec names both categories but defines features for neither.
    inline Category noise_category(const std::vector<double>& /*beat*/,
        double /*fs*/, uint32_t& reasons,
        const Config & /*cfg*/ = {})
    {
        reasons = REASON_NONE;
        return Category::UNCLASSIFIED;   // "no noise evidence"
    }

    // Provisional reference template for the FIRST pass, built without a
    // classification -- this is the bootstrap that breaks the circularity of
    // "exclude if r < 0.85 (against a template built from the beats that pass)".
    // TODO: per-sample median over all beats, or better, over the largest
    //       tight correlation cluster; then classify against it, rebuild from
    //       the survivors, and optionally iterate once.
    inline std::vector<double> bootstrap_reference(
        const std::vector<std::vector<double>>& /*beats*/)
    {
        return {};
    }

    // Reference template from the accepted beats only (second pass onward).
    // TODO: per-sample median across beats whose label is BONAFIDE, matching
    //       how CreateEcgTemplates already aggregates.
    inline std::vector<double> reference_from_accepted(
        const std::vector<std::vector<double>>& /*beats*/,
        const std::vector<BeatLabel>& /*labels*/)
    {
        return {};
    }

    // Correlation window. Whole-beat r drops on ischaemic ST drift, which is a
    // category-1 beat we must KEEP (detecting that drift is the whole point of
    // the long envelope in envelopes.hpp). So default to the QRS only.
    // TODO: [q_onset, s_end] from FeatureMarks when cfg.correlate_qrs_only.
    inline std::pair<int, int> correlation_window(int /*r_col*/, double /*fs*/,
        const Config & /*cfg*/ = {})
    {
        return { -1, -1 };
    }

    // ---- Ensemble hook ----------------------------------------------------
    // Spec 4.5 wants Bayesian + nonlinear + neural voting. The rule baseline
    // above is the interpretable first member. A second member is injected
    // here rather than compiled in, so the neural model can arrive later
    // without touching this header: given a beat, return a per-category
    // posterior. Empty return = "no opinion", and classify() falls back to
    // the rules alone.
    using PosteriorFn = std::function<std::vector<double>(
        const std::vector<double>& beat, double fs)>;

    // ---- Top level --------------------------------------------------------
    // One label per beat, in the order documented at the top of this file.
    // TODO:
    //   1. find_artifact_r  -> mark 3, drop those RRs from a working copy
    //   2. per-beat noise evidence -> provisional flags
    //   3. vote_flag over those flags -> confirm 4/5
    //   4. bootstrap_reference, then per-beat correlation
    //   5. is_premature AND (low corr OR wide QRS) -> 2; note compensatory pause
    //   6. everything left -> 1
    //   7. rebuild reference from accepted, re-run 4-6 once
    //   8. optional: blend in posterior() where it is supplied
    inline std::vector<BeatLabel> classify(
        const std::vector<std::vector<double>>& beats,
        const std::vector<double>& rr_ms,
        int /*r_col*/, double /*fs*/,
        const Config & /*cfg*/ = {},
        PosteriorFn /*posterior*/ = nullptr)
    {
        std::vector<BeatLabel> out(beats.size());
        for (size_t i = 0; i < out.size(); ++i)
            out[i].rr_ms = (i < rr_ms.size()) ? rr_ms[i]
            : std::numeric_limits<double>::quiet_NaN();
        return out;
    }

}   // namespace beat_classify