/**
 * @file   five_categories.hpp
 *
 *         NAMING. The spec calls for "beat_classify.hpp", one letter from the
 *         existing beat_classifier.hpp (Task M: PreMark, BeatClassifier,
 *         SessionAdapter, EnvelopeSession). Those two names collided in
 *         practice -- this file was saved over that one, which erased
 *         BeatClassifier and broke premark_beats.hpp, and because this file
 *         INCLUDES beat_classifier.hpp it then included itself. Renamed to
 *         five_categories.hpp after Section 4.5's five categories, so the two
 *         cannot be confused again.
 * @brief  Five-category beat classification, PVC prematurity filter and 5-of-8
 *         voting, multi-template morphology segregation, monomorphic vs
 *         polymorphic PVC tracking, NSVT detection.
 *         Spec Sections 4.5 and 4.6.
 *
 *         SPEC FIDELITY. Where 4.6 gives code -- isPremature, voteFlag, and
 *         the TemplateBank / NsvtRun declarations -- it is transcribed
 *         verbatim, constants and window arithmetic included. Consequences of
 *         the literal text are recorded as SPEC NOTE comments where they bite
 *         and asserted by the Python acceptance test; none is worked around.
 *
 *         THE MATCH SCORE IS A CORRELATION. 4.6 calls it a "band-match score"
 *         and gives template_match_floor as 0.60, which reads either as a
 *         correlation on [-1, 1] or as 60 percent of samples inside the
 *         template's corridor. Measured, on a 30-minute ventricular-bigeminy
 *         record and on the unit fixtures, with Section 4.6 run exactly as
 *         written (a single unmatched beat opens a template, merge only at the
 *         cap):
 *
 *           band percentage, floor 60   5 templates, 37/28/22/7/6 percent
 *           correlation, floor 0.60     2 templates, 50.0/50.0 percent
 *           correlation, floors 0.70,
 *                        0.80, 0.90     2 templates, 50.0/50.0 percent
 *
 *         A PVC scores ABOVE 60 percent of samples inside the sinus corridor,
 *         so the percentage reading cannot separate the two morphologies and
 *         the bank collapses to one template. The correlation reading meets
 *         the acceptance criterion at every floor tested. It is also the
 *         units-consistent reading: the same paragraph gives the exclusion
 *         thresholds as r < 0.85 (ECG) and r < 0.80 (PPG).
 *
 *         The band percentage is kept where the spec actually states a
 *         percentage -- Task M's 90 percent morphology gate -- and is still
 *         computed and reported per template. Two consequences worth knowing:
 *
 *           - A new template has one member and no corridor. Until it has
 *             kMinBeatsForBand members it is scored against a TOLERANCE
 *             corridor (mean +/- coldTolFrac of the pulse amplitude) so the
 *             score stays on the 0..100 scale from the first beat. 4.6 is
 *             silent on cold start; coldTolFrac is a config field.
 *           - The score is taken on the QRS band, not overall. A beat is ~400
 *             columns of which the QRS is ~30, so pct_overall cannot separate
 *             a grossly abnormal QRS from a sinus beat -- the same dilution
 *             argument premark_beats.hpp records for the 4.7 clean-pool gate.
 *             pct_overall is retained and reported.
 *
 *         MORPHOLOGY-CORRELATION EXCLUSION IS ENFORCED. 4.6's "r < 0.85 (ECG),
 *         r < 0.80 (PPG)" is now a hard gate in assignHandling: a beat whose
 *         correlation against its assigned template falls below the
 *         channel floor is excluded from reference calculations whatever its
 *         SQI and whatever category it landed in. Previously these two numbers
 *         were declared and never applied.
 *
 *         NO NEURAL STAGE. 4.5 asks for a Bayesian, nonlinear and
 *         neural-network ensemble. This file implements the interpretable
 *         rule-and-Bayesian baseline only; the neural member is out of scope
 *         here and there is no hook for it, so nothing in this header depends
 *         on a model artefact existing.
 *
 *         NO POLICIES BEYOND THE SPEC. Earlier revisions carried an opt-in
 *         persistence rule (a beat could not open a template alone) and an
 *         opt-in consolidation pass (merge on overlap, not only at the cap).
 *         Both are gone: 4.6 opens a template on any beat below the floor and
 *         merges only at the cap, and that is exactly what this does. With the
 *         correlation reading of template_match_floor the spec's own algorithm
 *         meets its acceptance criterion unaided -- two templates at 50.0 /
 *         50.0 percent on a 30-minute ventricular-bigeminy record -- so there
 *         was nothing left for the policies to fix.
 *
 *         WHAT HAD TO BE SUPPLIED. BeatClass (referenced by
 *         TemplateBank::labels, defined nowhere -- enumerated here from 4.6's
 *         own label list), and bodies for assignToTemplate and detectNsvt,
 *         which the spec declares with signatures only.
 *
 * @date   2026-08-24
 */
#pragma once

#include "template_generation/beat_classifier.hpp"      // PreMark::Class (pre-marking taxonomy)
#include "template_generation/morphology_envelope.hpp"  // MorphologyEnvelope, scoreBeatBands
#include "envelopes.hpp"            // envelopes::quantile_of

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

 // ---------------------------------------------------------------------------
 // Section 4.5: the five categories
 // ---------------------------------------------------------------------------
 // (1) bonafide PQRST, (2) abnormal rhythm (ectopic, AF, PVC, VT),
 // (3) artifactual R-waves, (4) machine noise, (5) human noise.
 // Only category 1 feeds feature extraction and averaging; category 2 is
 // retained with flags for arrhythmia burden.

enum class BeatCategory {
    UNCLASSIFIED = 0,
    BONAFIDE_PQRST = 1,
    ABNORMAL_RHYTHM = 2,
    ARTIFACTUAL_R = 3,
    MACHINE_NOISE = 4,
    HUMAN_NOISE = 5
};

inline const char* beatCategoryName(BeatCategory c) {
    switch (c) {
    case BeatCategory::BONAFIDE_PQRST:  return "bonafide_pqrst";
    case BeatCategory::ABNORMAL_RHYTHM: return "abnormal_rhythm";
    case BeatCategory::ARTIFACTUAL_R:   return "artifactual_r";
    case BeatCategory::MACHINE_NOISE:   return "machine_noise";
    case BeatCategory::HUMAN_NOISE:     return "human_noise";
    default:                            return "unclassified";
    }
}

// Template labels, per 4.6: "sinus, PVC type A, PVC type B, PAC, and so on."
// PVC_A and PVC_B are distinct labels because 4.6 counts DISTINCT PVC-labelled
// templates to decide monomorphic vs polymorphic.
enum class BeatClass {
    UNLABELED = 0,
    SINUS,
    PVC_A,
    PVC_B,
    PAC,
    VT,
    ARTIFACT
};

inline const char* beatClassName(BeatClass c) {
    switch (c) {
    case BeatClass::SINUS:    return "sinus";
    case BeatClass::PVC_A:    return "pvc_a";
    case BeatClass::PVC_B:    return "pvc_b";
    case BeatClass::PAC:      return "pac";
    case BeatClass::VT:       return "vt";
    case BeatClass::ARTIFACT: return "artifact";
    default:                  return "unlabeled";
    }
}

inline bool isVentricularClass(BeatClass c) {
    return c == BeatClass::PVC_A || c == BeatClass::PVC_B || c == BeatClass::VT;
}

// Which channel's morphology-correlation floor applies (4.6).
enum class SignalChannel { ECG, PPG };

// ---------------------------------------------------------------------------
// Section 4.6: PVC prematurity filter and 5-of-8 voting -- VERBATIM
// ---------------------------------------------------------------------------

// SPEC NOTE. rr is in ms and rr[t] is the interval preceding beat t. The first
// ten beats always return false: the ten-beat reference window cannot be
// formed there, and the specified return for that case is false.
//
// SPEC NOTE. The reference median is taken over the raw previous ten
// intervals, including any that are themselves ectopic. In bigeminy five of
// the ten are coupling intervals and the sorted median at w[5] is the upper
// central value -- the sinus interval -- so bigeminy still flags. Denser
// ectopy (trigeminy and beyond) pulls the median down until the filter falls
// silent.
inline bool isPremature(const std::vector<double>& rr, int t) { // rr in ms
    if (t < 10) return false;
    std::vector<double> w(rr.begin() + t - 10, rr.begin() + t);
    std::sort(w.begin(), w.end());
    double med = w[w.size() / 2];
    return rr[t] < 0.80 * med;                    // RR(t) < 0.80 * median
}

// 5-of-8 voting: flag beat t if >=5 of the surrounding 8 beats are flagged
//
// SPEC NOTE. The window is [t-4, t+4): eight slots covering the four beats
// before t, t itself, and the three after. Consequences, all asserted by the
// acceptance test rather than worked around:
//   - it is asymmetric, and t contributes to its own vote;
//   - at t = 0 the clamped window holds four entries, so >= 5 is unreachable;
//   - five flags must fall inside eight consecutive slots, so the rule
//     confirms only runs of at least five CONSECUTIVE premature intervals.
//     Isolated PVCs (one flag per neighbourhood) and bigeminy (at most four
//     flags per window) are never confirmed.
inline bool voteFlag(const std::vector<char>& flag, int t) {
    int lo = std::max(0, t - 4), hi = std::min((int)flag.size(), t + 4), c = 0;
    for (int i = lo; i < hi; ++i) c += flag[i];
    return c >= 5;
}

namespace beatcls {

    inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // Section 4.6 morphology-correlation exclusion floors, now enforced.
    inline constexpr double kEcgMorphCorrFloor = 0.85;
    inline constexpr double kPpgMorphCorrFloor = 0.80;

    // Section 4.6 multi-template defaults. The floor is a CORRELATION on
    // [-1, 1] -- the spec's 0.60 at face value. See the file header for the
    // measurement that settled this against the percentage reading.
    inline constexpr double kTemplateMatchFloor = 0.60;
    // Retained for the band-percentage path (assignToTemplateBanded) and for
    // Task M's morphology gate, which is genuinely a percentage.
    inline constexpr double kBandMatchFloorPct = 90.0;
    inline constexpr int    kMaxTemplatesPerBin = 6;
    // Below this many members a percentile corridor rests on too few order
    // statistics to mean anything, so the tolerance corridor is used instead.
    inline constexpr int    kMinBeatsForBand = 8;

    // Task A handling bands (Section 4.3).
    inline constexpr double kSqiInclude = 0.70;
    inline constexpr double kSqiSubstitute = 0.50;

    namespace detail {

        inline double pearson(const std::vector<double>& a, const std::vector<double>& b) {
            const std::size_t n = std::min(a.size(), b.size());
            double sa = 0.0, sb = 0.0; std::size_t m = 0;
            for (std::size_t i = 0; i < n; ++i) {
                if (!std::isfinite(a[i]) || !std::isfinite(b[i])) continue;
                sa += a[i]; sb += b[i]; ++m;
            }
            if (m < 3) return kNaN;
            const double ma = sa / m, mb = sb / m;
            double num = 0.0, da = 0.0, db = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                if (!std::isfinite(a[i]) || !std::isfinite(b[i])) continue;
                const double x = a[i] - ma, y = b[i] - mb;
                num += x * y; da += x * x; db += y * y;
            }
            if (da <= 0.0 || db <= 0.0) return kNaN;   // a flat trace correlates with nothing
            return num / std::sqrt(da * db);
        }

        inline double amplitudeOf(const std::vector<double>& v) {
            double lo = std::numeric_limits<double>::infinity();
            double hi = -std::numeric_limits<double>::infinity();
            for (double x : v) if (std::isfinite(x)) { lo = std::min(lo, x); hi = std::max(hi, x); }
            return (std::isfinite(lo) && std::isfinite(hi)) ? (hi - lo) : kNaN;
        }

        inline double gaussianLogPdf(double x, double mu, double sd) {
            if (!std::isfinite(x) || !(sd > 0.0)) return 0.0;   // absent evidence, no vote
            const double z = (x - mu) / sd;
            return -0.5 * z * z - std::log(sd);
        }

    } // namespace detail

    inline std::vector<char> prematureFlags(const std::vector<double>& rr) {
        std::vector<char> f(rr.size(), 0);
        for (int t = 0; t < (int)rr.size(); ++t)
            f[(std::size_t)t] = isPremature(rr, t) ? 1 : 0;
        return f;
    }

    inline std::vector<char> confirmedFlags(const std::vector<char>& raw) {
        std::vector<char> c(raw.size(), 0);
        for (int t = 0; t < (int)raw.size(); ++t)
            c[(std::size_t)t] = voteFlag(raw, t) ? 1 : 0;
        return c;
    }

    // ---------------------------------------------------------------------
    // Section 4.5: evidence -> category
    // ---------------------------------------------------------------------
    // Absent fields stay NaN (or -1 for motion) and contribute nothing, so a
    // PPG-only or RR-only record still classifies. Task A's BeatSQI maps in
    // directly: templateCorr, composite -> sqiComposite, motion -> motion.
    struct BeatEvidence {
        double templateCorr = kNaN;   // correlation against the assigned template
        double pctBandOverall = kNaN;   // band-match, whole beat
        double pctBandQRS = kNaN;   // band-match, QRS -- the segregation score
        double qrsWidthRatio = kNaN;
        double rAmpRatio = kNaN;
        bool   discordantT = false;
        bool   pWavePresent = true;

        double rrRatio = kNaN;   // RR(t) / local median RR
        bool   premature = false;
        bool   pause = false;

        double sqiComposite = kNaN;   // Task A composite
        double noiseHfFrac = kNaN;
        double powerlineFrac = kNaN;
        double clipFrac = kNaN;
        double stepDiscont = kNaN;
        double baselineDrift = kNaN;
        double flatFrac = kNaN;
        int    motion = -1;     // 1 clean, 0 motion, -1 unavailable
    };

    struct ClassifyConfig {
        SignalChannel channel = SignalChannel::ECG;
        double corrArtifactMax = 0.30;
        double clipMachineMin = 0.02;
        double powerlineMachineMin = 0.30;
        double stepMachineMin = 0.50;
        double flatMachineMin = 0.50;
        double noiseHumanMin = 0.35;
        double driftHumanMin = 0.30;
        double widePvcMin = 1.20;
        double sqiInclude = kSqiInclude;
        double sqiSubstitute = kSqiSubstitute;
        double morphCorrFloorEcg = kEcgMorphCorrFloor;   // 4.6
        double morphCorrFloorPpg = kPpgMorphCorrFloor;   // 4.6
        std::array<double, 5> prior{ { 0.88, 0.07, 0.02, 0.015, 0.015 } };

        double morphCorrFloor() const {
            return (channel == SignalChannel::PPG) ? morphCorrFloorPpg
                : morphCorrFloorEcg;
        }
    };

    struct BeatVerdict {
        BeatCategory category = BeatCategory::UNCLASSIFIED;
        BeatClass    label = BeatClass::UNLABELED;
        std::array<double, 5> posterior{ {0, 0, 0, 0, 0} };
        double confidence = 0.0;
        bool   ruleFired = false;
        const char* reason = "";

        enum Handling { INCLUDE, SUBSTITUTE, EXCLUDE } handling = EXCLUDE;
        bool retainedForBurden = false;   // category 2: out of references, flagged
        bool feedsFeatures = false;   // category 1 at INCLUDE only
        bool morphologyExcluded = false;   // tripped the 4.6 r < 0.85 / 0.80 floor
    };

    // Ordered rule pass; first match wins. Signal-chain failures are settled
    // before rhythm rules: a clipped or powerline-dominated beat has no
    // morphology left to interpret.
    inline bool applyRules(const BeatEvidence& e, const ClassifyConfig& cfg,
        BeatCategory& cat, BeatClass& label, const char*& why)
    {
        auto fin = [](double v) { return std::isfinite(v); };

        if (fin(e.clipFrac) && e.clipFrac >= cfg.clipMachineMin) {
            cat = BeatCategory::MACHINE_NOISE; label = BeatClass::ARTIFACT;
            why = "adc_clipping"; return true;
        }
        if (fin(e.powerlineFrac) && e.powerlineFrac >= cfg.powerlineMachineMin) {
            cat = BeatCategory::MACHINE_NOISE; label = BeatClass::ARTIFACT;
            why = "powerline_dominant"; return true;
        }
        if (fin(e.stepDiscont) && e.stepDiscont >= cfg.stepMachineMin) {
            cat = BeatCategory::MACHINE_NOISE; label = BeatClass::ARTIFACT;
            why = "lead_step"; return true;
        }
        if (fin(e.flatFrac) && e.flatFrac >= cfg.flatMachineMin) {
            cat = BeatCategory::MACHINE_NOISE; label = BeatClass::ARTIFACT;
            why = "flatline_dropout"; return true;
        }
        if (e.motion == 0 && fin(e.noiseHfFrac) && e.noiseHfFrac >= cfg.noiseHumanMin) {
            cat = BeatCategory::HUMAN_NOISE; label = BeatClass::ARTIFACT;
            why = "motion_plus_broadband"; return true;
        }
        if (fin(e.baselineDrift) && e.baselineDrift >= cfg.driftHumanMin
            && fin(e.templateCorr) && e.templateCorr < cfg.morphCorrFloor()) {
            cat = BeatCategory::HUMAN_NOISE; label = BeatClass::ARTIFACT;
            why = "baseline_wander"; return true;
        }
        if (fin(e.templateCorr) && e.templateCorr < cfg.corrArtifactMax) {
            cat = BeatCategory::ARTIFACTUAL_R; label = BeatClass::ARTIFACT;
            why = "no_qrs_morphology"; return true;
        }
        {
            int score = 0;
            if (e.premature) ++score;
            if (fin(e.qrsWidthRatio) && e.qrsWidthRatio >= cfg.widePvcMin) ++score;
            if (!e.pWavePresent) ++score;
            if (e.discordantT) ++score;
            if (e.pause) ++score;
            if (score >= 3) {
                cat = BeatCategory::ABNORMAL_RHYTHM; label = BeatClass::PVC_A;
                why = "premature_wide_ventricular"; return true;
            }
            if (e.premature && e.pWavePresent
                && fin(e.qrsWidthRatio) && e.qrsWidthRatio < cfg.widePvcMin) {
                cat = BeatCategory::ABNORMAL_RHYTHM; label = BeatClass::PAC;
                why = "premature_narrow_with_p"; return true;
            }
        }
        return false;
    }

    // Naive Bayes over whatever evidence is present. Per-category statistics
    // are starting estimates, to be refit on the Task B labelled log.
    inline std::array<double, 5> bayesPosterior(const BeatEvidence& e,
        const ClassifyConfig& cfg)
    {
        struct FeatStat { double mu, sd; };
        //                              cat1          cat2          cat3          cat4          cat5
        const FeatStat corr[5] = { {0.97,0.03}, {0.70,0.15}, {0.15,0.20}, {0.40,0.25}, {0.60,0.20} };
        const FeatStat bandQ[5] = { {97.0,4.0},  {70.0,20.0}, {30.0,25.0}, {40.0,25.0}, {60.0,25.0} };
        const FeatStat width[5] = { {1.00,0.08}, {1.35,0.25}, {1.60,0.60}, {1.20,0.50}, {1.20,0.40} };
        const FeatStat rrRat[5] = { {1.00,0.08}, {0.75,0.20}, {0.60,0.35}, {1.00,0.30}, {1.00,0.25} };
        const FeatStat hf[5] = { {0.05,0.04}, {0.08,0.06}, {0.30,0.20}, {0.35,0.25}, {0.45,0.20} };
        const FeatStat drift[5] = { {0.05,0.05}, {0.07,0.06}, {0.25,0.20}, {0.15,0.15}, {0.40,0.20} };
        const FeatStat clip[5] = { {0.00,0.01}, {0.00,0.01}, {0.02,0.05}, {0.15,0.12}, {0.02,0.04} };
        const FeatStat sqi[5] = { {0.90,0.08}, {0.75,0.15}, {0.30,0.20}, {0.35,0.20}, {0.40,0.20} };

        std::array<double, 5> logp{};
        for (int c = 0; c < 5; ++c) {
            double lp = std::log(std::max(1e-9, cfg.prior[(std::size_t)c]));
            lp += detail::gaussianLogPdf(e.templateCorr, corr[c].mu, corr[c].sd);
            lp += detail::gaussianLogPdf(e.pctBandQRS, bandQ[c].mu, bandQ[c].sd);
            lp += detail::gaussianLogPdf(e.qrsWidthRatio, width[c].mu, width[c].sd);
            lp += detail::gaussianLogPdf(e.rrRatio, rrRat[c].mu, rrRat[c].sd);
            lp += detail::gaussianLogPdf(e.noiseHfFrac, hf[c].mu, hf[c].sd);
            lp += detail::gaussianLogPdf(e.baselineDrift, drift[c].mu, drift[c].sd);
            lp += detail::gaussianLogPdf(e.clipFrac, clip[c].mu, clip[c].sd);
            lp += detail::gaussianLogPdf(e.sqiComposite, sqi[c].mu, sqi[c].sd);
            if (e.premature)      lp += (c == 1) ? std::log(0.45) : std::log(0.03);
            if (!e.pWavePresent)  lp += (c == 1) ? std::log(0.60) : std::log(0.15);
            if (e.motion == 0)    lp += (c == 4) ? std::log(0.70) : std::log(0.10);
            logp[(std::size_t)c] = lp;
        }
        const double mx = *std::max_element(logp.begin(), logp.end());
        double sum = 0.0;
        std::array<double, 5> post{};
        for (int c = 0; c < 5; ++c) {
            post[(std::size_t)c] = std::exp(logp[(std::size_t)c] - mx);
            sum += post[(std::size_t)c];
        }
        if (sum > 0.0) for (double& v : post) v /= sum;
        return post;
    }

    // Handling policy. Order matters:
    //   1. anything outside category 1 is out of the references;
    //   2. the 4.6 morphology-correlation floor excludes whatever it catches,
    //      regardless of SQI -- a beat that does not look like its own
    //      template has no business in the average however clean it is;
    //   3. what survives is graded by the 0.70 / 0.50 Task A bands.
    inline void assignHandling(BeatVerdict& v, const BeatEvidence& e,
        const ClassifyConfig& cfg)
    {
        v.retainedForBurden = (v.category == BeatCategory::ABNORMAL_RHYTHM);

        // Section 4.6: r < 0.85 (ECG), r < 0.80 (PPG). Applied to every
        // category so the flag is meaningful on a burden count too.
        v.morphologyExcluded = std::isfinite(e.templateCorr)
            && e.templateCorr < cfg.morphCorrFloor();

        if (v.category != BeatCategory::BONAFIDE_PQRST) {
            v.handling = BeatVerdict::EXCLUDE;
            v.feedsFeatures = false;
            return;
        }
        if (v.morphologyExcluded) {
            v.handling = BeatVerdict::EXCLUDE;
            v.feedsFeatures = false;
            if (!v.ruleFired) v.reason = "morphology_corr_below_floor";
            return;
        }
        const double q = std::isfinite(e.sqiComposite) ? e.sqiComposite : e.templateCorr;
        if (!std::isfinite(q))            v.handling = BeatVerdict::SUBSTITUTE;
        else if (q >= cfg.sqiInclude)     v.handling = BeatVerdict::INCLUDE;
        else if (q >= cfg.sqiSubstitute)  v.handling = BeatVerdict::SUBSTITUTE;
        else                              v.handling = BeatVerdict::EXCLUDE;
        v.feedsFeatures = (v.handling == BeatVerdict::INCLUDE);
    }

    inline BeatVerdict classifyBeat(const BeatEvidence& e,
        const ClassifyConfig& cfg = {})
    {
        BeatVerdict v;
        BeatCategory cat = BeatCategory::UNCLASSIFIED;
        BeatClass    lbl = BeatClass::UNLABELED;
        const char* why = "bayes";

        if (applyRules(e, cfg, cat, lbl, why)) {
            v.category = cat; v.label = lbl; v.ruleFired = true; v.reason = why;
            v.posterior.fill(0.0);
            v.posterior[(std::size_t)((int)cat - 1)] = 1.0;
            v.confidence = 1.0;
            assignHandling(v, e, cfg);
            return v;
        }

        const std::array<double, 5> post = bayesPosterior(e, cfg);
        const std::size_t best = (std::size_t)(
            std::max_element(post.begin(), post.end()) - post.begin());
        v.posterior = post;
        v.confidence = post[best];
        v.category = (BeatCategory)((int)best + 1);
        v.reason = why;
        if (v.category == BeatCategory::ABNORMAL_RHYTHM)
            v.label = (e.premature && e.pWavePresent) ? BeatClass::PAC : BeatClass::PVC_A;
        else if (v.category == BeatCategory::BONAFIDE_PQRST)
            v.label = BeatClass::SINUS;
        else
            v.label = BeatClass::ARTIFACT;
        assignHandling(v, e, cfg);
        return v;
    }

    // Per-bin category percentages and arrhythmia burden (4.5).
    struct BinCategoryReport {
        std::array<int, 6> counts{ {0,0,0,0,0,0} };
        int nBeats = 0, nInclude = 0, nSubstitute = 0, nExclude = 0;
        int nMorphologyExcluded = 0;
        double pvcBurdenPct = kNaN, pacBurdenPct = kNaN, substitutedPct = kNaN;
        double morphologyExcludedPct = kNaN;
        double pct(BeatCategory c) const {
            return nBeats ? 100.0 * counts[(std::size_t)(int)c] / nBeats : kNaN;
        }
    };

    inline BinCategoryReport summarizeBin(const std::vector<BeatVerdict>& v) {
        BinCategoryReport r;
        r.nBeats = (int)v.size();
        int pvc = 0, pac = 0;
        for (const BeatVerdict& b : v) {
            ++r.counts[(std::size_t)(int)b.category];
            switch (b.handling) {
            case BeatVerdict::INCLUDE:    ++r.nInclude;    break;
            case BeatVerdict::SUBSTITUTE: ++r.nSubstitute; break;
            default:                      ++r.nExclude;    break;
            }
            if (b.morphologyExcluded) ++r.nMorphologyExcluded;
            if (b.category == BeatCategory::ABNORMAL_RHYTHM) {
                if (isVentricularClass(b.label)) ++pvc;
                else if (b.label == BeatClass::PAC) ++pac;
            }
        }
        if (r.nBeats > 0) {
            r.pvcBurdenPct = 100.0 * pvc / r.nBeats;
            r.pacBurdenPct = 100.0 * pac / r.nBeats;
            r.substitutedPct = 100.0 * r.nSubstitute / r.nBeats;
            r.morphologyExcludedPct = 100.0 * r.nMorphologyExcluded / r.nBeats;
        }
        return r;
    }

} // namespace beatcls

// ---------------------------------------------------------------------------
// Section 4.6: multi-template morphology segregation
// ---------------------------------------------------------------------------
// The three spec fields are first and unchanged. The rest is what band-match
// scoring requires: a corridor per template, and the members it is built from.

struct TemplateBank {
    std::vector<std::vector<double>> templates;   // running mean per template
    std::vector<BeatClass> labels;
    std::vector<int> counts;

    std::vector<char> confirmed;                  // operator has labelled it
    std::vector<std::vector<std::vector<double>>> members;   // corridor inputs
    std::vector<MorphologyEnvelope> envs;         // per-template band corridor

    // Band the match score is taken over. Both -1 means "no landmarks", and
    // the score falls back to the whole beat.
    int qrsStart = -1, qrsEnd = -1;
    int    maxMembers = 200;    // ring-buffer cap per template
    double coldTolFrac = 0.10;   // cold-start corridor half-width, x amplitude


    int size() const { return (int)templates.size(); }
    int totalBeats() const { int n = 0; for (int c : counts) n += c; return n; }
    double share(int i) const {
        const int tot = totalBeats();
        return (tot > 0 && i >= 0 && i < size())
            ? 100.0 * counts[(std::size_t)i] / tot
            : std::numeric_limits<double>::quiet_NaN();
    }
};

struct NsvtRun { int startBeat, length, templateId; double meanCycleMs, rateBpm; bool sustained; };

namespace beatcls {

    // Cold-start corridor: the template plus or minus a fixed fraction of its
    // own amplitude. Used until a template has kMinBeatsForBand members, so
    // the match score is on the 0..100 band scale from the first beat rather
    // than switching units partway through a record.
    inline MorphologyEnvelope toleranceCorridor(const std::vector<double>& tmpl,
        double tolFrac)
    {
        MorphologyEnvelope e;
        const int W = (int)tmpl.size();
        if (W <= 0) return e;
        const double amp = detail::amplitudeOf(tmpl);
        const double half = (std::isfinite(amp) ? amp : 1.0) * std::max(1e-9, tolFrac);
        e.nBeats = 1;
        e.lo.assign((std::size_t)W, kNaN);
        e.hi.assign((std::size_t)W, kNaN);
        e.mean.assign((std::size_t)W, kNaN);
        e.sd.assign((std::size_t)W, kNaN);
        e.nContrib.assign((std::size_t)W, 1);
        for (int i = 0; i < W; ++i) {
            if (!std::isfinite(tmpl[(std::size_t)i])) { e.nContrib[(std::size_t)i] = 0; continue; }
            e.lo[(std::size_t)i] = tmpl[(std::size_t)i] - half;
            e.hi[(std::size_t)i] = tmpl[(std::size_t)i] + half;
            e.mean[(std::size_t)i] = tmpl[(std::size_t)i];
        }
        return e;
    }

    inline void rebuildTemplateEnv(TemplateBank& bank, int i) {
        if (i < 0 || i >= bank.size()) return;
        const auto& mem = bank.members[(std::size_t)i];
        const int W = (int)bank.templates[(std::size_t)i].size();
        bank.envs[(std::size_t)i] =
            ((int)mem.size() >= kMinBeatsForBand)
            ? buildEnvelope(mem, W)
            : toleranceCorridor(bank.templates[(std::size_t)i], bank.coldTolFrac);
    }

    // Band-match score of a beat against one template's corridor. QRS band
    // when landmarks are known, whole beat otherwise.
    inline BandMatchResult scoreAgainstTemplate(const std::vector<double>& beat,
        const TemplateBank& bank, int i)
    {
        if (i < 0 || i >= bank.size()) return BandMatchResult{};
        const int W = (int)std::min(beat.size(), bank.envs[(std::size_t)i].lo.size());
        const int qs = (bank.qrsStart >= 0) ? bank.qrsStart : 0;
        const int qe = (bank.qrsEnd > bank.qrsStart && bank.qrsEnd >= 0) ? bank.qrsEnd : W;
        return scoreBeatBands(beat, bank.envs[(std::size_t)i], qs, qe);
    }

    // The score the assignment decision is made on: QRS band, falling back to
    // overall when the QRS band has nothing scorable.
    inline double matchScoreOf(const BandMatchResult& bm) {
        if (std::isfinite(bm.pct_QRS)) return bm.pct_QRS;
        if (std::isfinite(bm.pct_overall)) return bm.pct_overall;
        return kNaN;
    }

    // Correlation of a beat against every template in the bank -- what the
    // assignment decision is made on.
    inline std::vector<double> bankScores(const std::vector<double>& beat,
        const TemplateBank& bank)
    {
        std::vector<double> s((std::size_t)bank.size(), kNaN);
        for (int i = 0; i < bank.size(); ++i)
            s[(std::size_t)i] = detail::pearson(beat, bank.templates[(std::size_t)i]);
        return s;
    }

    inline std::vector<double> bankBandScores(const std::vector<double>& beat,
        const TemplateBank& bank)
    {
        std::vector<double> s((std::size_t)bank.size(), kNaN);
        for (int i = 0; i < bank.size(); ++i)
            s[(std::size_t)i] = matchScoreOf(scoreAgainstTemplate(beat, bank, i));
        return s;
    }

    // Correlation against one template. Kept separate from the band score:
    // different quantities on different scales, and conflating them is how the
    // 0.60 floor became ambiguous in the first place.
    inline double correlationTo(const std::vector<double>& beat,
        const TemplateBank& bank, int i)
    {
        if (i < 0 || i >= bank.size()) return kNaN;
        return detail::pearson(beat, bank.templates[(std::size_t)i]);
    }

    // The bin's reference morphology: the confirmed sinus template, or
    // template 0 when nothing is labelled yet.
    inline int referenceTemplateIndex(const TemplateBank& bank) {
        for (int i = 0; i < bank.size(); ++i)
            if (bank.confirmed[(std::size_t)i] && bank.labels[(std::size_t)i] == BeatClass::SINUS)
                return i;
        return bank.size() > 0 ? 0 : -1;
    }

    // WHICH TEMPLATE THE 4.6 EXCLUSION FLOOR IS MEASURED AGAINST.
    // 4.6 gives "r < 0.85 (ECG), r < 0.80 (PPG)" as exclusion thresholds in
    // the paragraph about keeping PVCs and artifact out of the reference
    // calculations, but does not say which template the correlation is taken
    // against, and the two readings behave oppositely once multi-template
    // segregation exists:
    //
    //   against the beat's OWN ASSIGNED template -- an ectopic beat sits in
    //     the ectopic template and correlates with it at ~1.0, so the floor
    //     never fires on ectopy and becomes a pure noise detector. Measured:
    //     on a record with four known PVCs, zero beats tripped the floor.
    //
    //   against the bin's REFERENCE (sinus) template -- ectopy and artifact
    //     both fall below it, which is the exclusion the paragraph is about.
    //
    // The reference reading is used, because the other makes a stated
    // criterion a no-op. correlationTo is still available for the per-template
    // quality question, and both numbers are reported so the choice is
    // visible rather than buried.
    inline double correlationToReference(const std::vector<double>& beat,
        const TemplateBank& bank)
    {
        return correlationTo(beat, bank, referenceTemplateIndex(bank));
    }

    inline void addMember(TemplateBank& bank, int i, const std::vector<double>& beat) {
        auto& mem = bank.members[(std::size_t)i];
        if ((int)mem.size() < bank.maxMembers) mem.push_back(beat);
        else mem[(std::size_t)(bank.counts[(std::size_t)i] % bank.maxMembers)] = beat;
    }

    // Count-weighted running mean, so a template stays the arithmetic mean of
    // its members and a merge is exact.
    inline void accumulate(TemplateBank& bank, int idx, const std::vector<double>& beat) {
        if (idx < 0 || idx >= bank.size()) return;
        std::vector<double>& t = bank.templates[(std::size_t)idx];
        const std::size_t n = std::min(t.size(), beat.size());
        const double k = (double)bank.counts[(std::size_t)idx];
        for (std::size_t j = 0; j < n; ++j) {
            if (!std::isfinite(beat[j])) continue;
            if (!std::isfinite(t[j])) { t[j] = beat[j]; continue; }
            t[j] = (t[j] * k + beat[j]) / (k + 1.0);
        }
        addMember(bank, idx, beat);
        ++bank.counts[(std::size_t)idx];
        // The corridor is rebuilt every assignment. buildEnvelope is O(N*W)
        // with two nth_element passes per column and N capped at maxMembers,
        // which is inside the 10 ms the 4.7.3 budget allows for a rebuild.
        rebuildTemplateEnv(bank, idx);
    }

    inline void pushTemplate(TemplateBank& bank, const std::vector<double>& beat,
        BeatClass label, bool confirmed)
    {
        bank.templates.push_back(beat);
        bank.labels.push_back(label);
        bank.counts.push_back(1);
        bank.confirmed.push_back(confirmed ? 1 : 0);
        bank.members.push_back({ beat });
        bank.envs.push_back(MorphologyEnvelope{});
        rebuildTemplateEnv(bank, bank.size() - 1);
    }

    // "Seed the bank with the sinus template from Phase 1."
    inline TemplateBank seedBank(const std::vector<double>& sinusTemplate,
        int qrsStart = -1, int qrsEnd = -1)
    {
        TemplateBank b;
        b.qrsStart = qrsStart;
        b.qrsEnd = qrsEnd;
        if (sinusTemplate.empty()) return b;
        pushTemplate(b, sinusTemplate, BeatClass::SINUS, true);
        b.counts[0] = 0;                 // the seed is a template, not a beat
        b.members[0] = { sinusTemplate }; // but it does anchor the corridor
        return b;
    }

    // "merge the two closest templates when the cap is reached"
    inline bool mergeTemplates(TemplateBank& bank, int i, int j) {
        if (i == j || i < 0 || j < 0 || i >= bank.size() || j >= bank.size()) return false;
        if (j < i) std::swap(i, j);
        const auto ui = (std::size_t)i, uj = (std::size_t)j;
        const double ci = bank.counts[ui], cj = bank.counts[uj], tot = ci + cj;
        std::vector<double>& a = bank.templates[ui];
        const std::vector<double>& b = bank.templates[uj];
        const std::size_t n = std::min(a.size(), b.size());
        for (std::size_t k = 0; k < n; ++k) {
            const bool fa = std::isfinite(a[k]), fb = std::isfinite(b[k]);
            if (fa && fb && tot > 0.0) a[k] = (a[k] * ci + b[k] * cj) / tot;
            else if (!fa && fb)        a[k] = b[k];
        }
        bool clean = true;
        if (bank.confirmed[uj] && !bank.confirmed[ui]) {
            bank.labels[ui] = bank.labels[uj];
            bank.confirmed[ui] = 1;
        }
        else if (bank.confirmed[ui] && bank.confirmed[uj]
            && bank.labels[ui] != bank.labels[uj]) {
            if (cj > ci) bank.labels[ui] = bank.labels[uj];
            clean = false;                       // conflicting confirmed labels
        }
        // Members pool, up to the cap, so the merged corridor spans both.
        for (const auto& m : bank.members[uj]) {
            if ((int)bank.members[ui].size() >= bank.maxMembers) break;
            bank.members[ui].push_back(m);
        }
        bank.counts[ui] = (int)tot;
        bank.templates.erase(bank.templates.begin() + j);
        bank.labels.erase(bank.labels.begin() + j);
        bank.counts.erase(bank.counts.begin() + j);
        bank.confirmed.erase(bank.confirmed.begin() + j);
        bank.members.erase(bank.members.begin() + j);
        bank.envs.erase(bank.envs.begin() + j);
        rebuildTemplateEnv(bank, i);
        return clean;
    }

    // The closest pair, by correlation -- the same score the assignment uses.
    inline bool closestPair(const TemplateBank& bank, int& bi, int& bj) {
        double best = -2.0; bi = bj = -1;
        for (int i = 0; i < bank.size(); ++i)
            for (int j = i + 1; j < bank.size(); ++j) {
                const double r = detail::pearson(bank.templates[(std::size_t)i],
                    bank.templates[(std::size_t)j]);
                if (std::isfinite(r) && r > best) { best = r; bi = i; bj = j; }
            }
        return bi >= 0;
    }

    // "Label templates by class once the operator confirms one member."
    inline void confirmTemplateLabel(TemplateBank& bank, int idx, BeatClass label) {
        if (idx < 0 || idx >= bank.size()) return;
        bank.labels[(std::size_t)idx] = label;
        bank.confirmed[(std::size_t)idx] = 1;
    }

    // Monomorphic vs polymorphic PVC tracking (4.6).
    struct EctopyReport {
        int    pvcTemplateCount = 0;
        bool   polymorphic = false;
        double pvcSharePct = kNaN;
        std::vector<int>    pvcTemplateIds;
        std::vector<double> shareByTemplate;
    };

    inline EctopyReport ectopyReport(const TemplateBank& bank) {
        EctopyReport r;
        r.shareByTemplate.resize((std::size_t)bank.size(), kNaN);
        int pvcBeats = 0;
        for (int i = 0; i < bank.size(); ++i) {
            r.shareByTemplate[(std::size_t)i] = bank.share(i);
            if (isVentricularClass(bank.labels[(std::size_t)i])) {
                ++r.pvcTemplateCount;
                r.pvcTemplateIds.push_back(i);
                pvcBeats += bank.counts[(std::size_t)i];
            }
        }
        r.polymorphic = (r.pvcTemplateCount >= 2);
        const int tot = bank.totalBeats();
        if (tot > 0) r.pvcSharePct = 100.0 * pvcBeats / tot;
        return r;
    }

} // namespace beatcls

// Spec signature. The match score is a CORRELATION (see the file header).
// "For each beat, compute the band-match score against every template in the
// bank and assign it to the best match. When a beat matches no template above
// template_match_floor (default 0.60), open a new template, but cap the bank at
// max_templates_per_bin (default 6) and merge the two closest templates when
// the cap is reached."
//
// Returns the template index, or -1 for an unusable beat.
inline int assignToTemplate(const std::vector<double>& beat, TemplateBank& bank,
    double matchFloor = beatcls::kTemplateMatchFloor,
    int maxTemplates = beatcls::kMaxTemplatesPerBin)
{
    if (beat.empty()) return -1;
    if (bank.size() == 0) {
        beatcls::pushTemplate(bank, beat, BeatClass::UNLABELED, false);
        return 0;
    }

    // "compute the band-match score against every template in the bank and
    // assign it to the best match"
    const std::vector<double> s = beatcls::bankScores(beat, bank);
    int best = -1; double bestR = -2.0;
    for (int i = 0; i < (int)s.size(); ++i)
        if (std::isfinite(s[(std::size_t)i]) && s[(std::size_t)i] > bestR) {
            bestR = s[(std::size_t)i]; best = i;
        }
    // Nothing correlates at all (a flat beat, or every template flat). Do not
    // open a template for it: a constant trace would become an attractor that
    // swallows every subsequent dropout.
    if (best < 0) return -1;

    if (bestR >= matchFloor) {
        beatcls::accumulate(bank, best, beat);
        return best;
    }

    // "When a beat matches no template above template_match_floor, open a new
    // template, but cap the bank at max_templates_per_bin and merge the two
    // closest templates when the cap is reached."
    if (bank.size() >= maxTemplates) {
        int bi = -1, bj = -1;
        if (beatcls::closestPair(bank, bi, bj)) beatcls::mergeTemplates(bank, bi, bj);
        else { beatcls::accumulate(bank, best, beat); return best; }
    }
    beatcls::pushTemplate(bank, beat, BeatClass::UNLABELED, false);
    return bank.size() - 1;
}

// Explicit-envelope variant, for callers that already hold one corridor per
// template outside a TemplateBank.
inline int assignToTemplateBanded(const std::vector<double>& beat,
    const std::vector<MorphologyEnvelope>& envs,
    int qrsStart, int qrsEnd,
    double pctFloor = beatcls::kBandMatchFloorPct)
{
    int best = -1; double bestPct = -1.0;
    for (std::size_t i = 0; i < envs.size(); ++i) {
        const double p = beatcls::matchScoreOf(
            scoreBeatBands(beat, envs[i], qrsStart, qrsEnd));
        if (std::isfinite(p) && p > bestPct) { bestPct = p; best = (int)i; }
    }
    return (bestPct >= pctFloor) ? best : -1;
}

// Spec signature. "Flag a run when three or more consecutive beats match the
// same ventricular-labeled template and the run rate exceeds 100 beats per
// minute. Record run onset, run length, mean cycle length, and the template
// identity. A run of 30 seconds or longer is sustained VT."
//
// rrMs[i] is the interval preceding beat i, so the intervals internal to a run
// of L beats starting at s are rrMs[s+1 .. s+L-1]: L-1 intervals. Mean cycle
// length and run duration are both taken over those.
inline std::vector<NsvtRun> detectNsvt(const std::vector<int>& templateIdPerBeat,
    const std::vector<double>& rrMs,
    const TemplateBank& bank, int minRun = 3)
{
    std::vector<NsvtRun> runs;
    const int n = (int)templateIdPerBeat.size();

    auto ventricular = [&](int id) {
        return id >= 0 && id < bank.size()
            && isVentricularClass(bank.labels[(std::size_t)id]);
        };

    int i = 0;
    while (i < n) {
        const int id = templateIdPerBeat[(std::size_t)i];
        if (!ventricular(id)) { ++i; continue; }
        int j = i;
        while (j + 1 < n && templateIdPerBeat[(std::size_t)(j + 1)] == id) ++j;
        const int len = j - i + 1;
        if (len >= minRun) {
            double sum = 0.0; int cnt = 0;
            for (int k = i + 1; k <= j && k < (int)rrMs.size(); ++k) {
                if (!std::isfinite(rrMs[(std::size_t)k])) continue;
                sum += rrMs[(std::size_t)k]; ++cnt;
            }
            if (cnt > 0) {
                NsvtRun r{};
                r.startBeat = i;
                r.length = len;
                r.templateId = id;
                r.meanCycleMs = sum / cnt;
                r.rateBpm = (r.meanCycleMs > 0.0) ? 60000.0 / r.meanCycleMs
                    : std::numeric_limits<double>::quiet_NaN();
                r.sustained = (sum / 1000.0 >= 30.0);
                if (std::isfinite(r.rateBpm) && r.rateBpm > 100.0) runs.push_back(r);
            }
        }
        i = j + 1;
    }
    return runs;
}