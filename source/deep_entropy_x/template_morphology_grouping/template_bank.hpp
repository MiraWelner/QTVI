#pragma once
/**
 * @file   template_bank.hpp
 * @brief  Multi-template morphology segregation (Spec Section 4.6).
 *
 *         One template per bin cannot represent a bin containing more than one
 *         beat morphology. Averaging sinus and ectopic beats produces a third
 *         shape matching neither and widens the per-sample corridor that every
 *         downstream feature is measured against -- so the damage is not
 *         confined to the ectopic beats, it degrades the sinus measurements
 *         too. This header holds the state for a small bank of templates per
 *         bin per channel, split by SHAPE and not by time: every template in a
 *         bank spans the bin's full 15 minutes, and a PVC at minute 2 shares a
 *         template with a PVC at minute 13.
 *
 *         Design decisions settled before this file was written, recorded here
 *         because none of them are recoverable from the spec text alone:
 *
 *          1. THE METRIC IS THE BAND-MATCH SCORE. Section 4.6 names the score
 *             and gives its thresholds -- 0.85 (ECG) and 0.80 (PPG), the same
 *             two numbers used for assignment and for spawning -- and those
 *             parts are spec.
 *
 *             THE BAND ITSELF IS NOT SPECIFIED. 4.6 says "band-match score"
 *             without saying what the band is. The 2.5/97.5 per-column
 *             corridor used here is a LOCAL CHOICE, made to agree with
 *             morphology_envelope.hpp, which builds 2.5/97.5 corridors for the
 *             envelope sections. Min/max, +/- 2 SD, or an IQR band would all
 *             satisfy the clause as written and would all move the score.
 *
 *             Flagged rather than left implicit because a threshold is
 *             meaningless without the quantity it compares, and this file used
 *             to present the corridor as though 4.6 stated it. Several
 *             constants below hang off the choice -- kMinMembersForCorridor is
 *             justified by the percentiles needing two distinct order
 *             statistics, and the corridor inflation and inheritance in
 *             template_assign.hpp exist because a small-N percentile band is a
 *             poor estimate -- so if the band is ever defined differently, all
 *             of that has to be rederived rather than retuned.
 *
 *             AS A FRACTION, NOT A PERCENTAGE. morphology_envelope.hpp scores
 *             on 0-100, which does not compare to 0.85. Read as a fraction the
 *             two are the same quantity and the threshold means "at least 85%
 *             of this beat's samples lie inside the corridor". That is the only
 *             self-consistent reading of the clause and it is what bandMatch()
 *             returns.
 *
 *             THE CORRIDOR OF A YOUNG TEMPLATE IS INHERITED, and this is the
 *             one thing the clause does not specify. A corridor is a spread
 *             estimate, and a template with one member has none:
 *             lo[c] == hi[c] at every column, so a second beat of the SAME
 *             morphology scores ~0 and spawns yet another template. Left alone
 *             that makes the spec's own metric unable to grow a template it
 *             just opened. So below kMinMembersForCorridor a template's
 *             corridor is widened to slot 0's spread at each column: "this
 *             morphology's variability is not yet known, and is assumed no
 *             tighter than sinus". The assumption is conservative in the right
 *             direction -- it can let a beat in, never keep one out -- and it
 *             introduces no second threshold, which the clause does forbid.
 *
 *             Pearson r is retained for ONE purpose only: reporting. It is
 *             written alongside the band-match score in the per-beat archive
 *             because the two disagree informatively -- r is blind to
 *             amplitude, the corridor is not, so a beat with high r and low
 *             band-match is an amplitude outlier of a known shape, which is a
 *             real finding and was invisible while r was the only number.
 *             Nothing routes on it.
 *
 *          2. ONE THRESHOLD, NOT TWO. The same floor decides "this beat
 *             belongs to that template" and "no new template is needed". The
 *             spec forbids a second looser assignment threshold, and the
 *             reason is worth keeping in view: a two-number scheme (assign at
 *             0.70, spawn below 0.85) files beats into templates they do not
 *             match, which is precisely the variance inflation this section
 *             exists to remove -- now hidden inside a cluster instead of
 *             visible as a bad median.
 *
 *          3. TEMPLATES UPDATE BY MEDIAN OVER MEMBERS. Column-wise
 *             NaN-skipping median, matching create_ecg_templates.hpp's
 *             medianOver(). NOT the alpha = 1/8 EWMA from the beat
 *             substitution section: an EWMA is a recursion whose result
 *             depends on arrival order, and pass 2 exists specifically to
 *             remove order dependence. Median-over-members also keeps the
 *             archive property that a template is exactly reconstructible
 *             from the per-beat flags -- anyone can recompute the median from
 *             the CSV and get the same numbers.
 *
 *          4. LABEL CODES ARE annotation_types CODES, NOT A PARALLEL ENUM.
 *             annotation_types.hpp is the stated single source of truth for
 *             marking codes; a BeatClass enum duplicating PVC/PAC/AF/SVT/VT
 *             here would be a second place for them to drift. label_code is
 *             the raw code, with 0 reserved for UNLABELED. Note that PVC is
 *             code 4 and PAC is code 5 in that table, while the addendum text
 *             says "PVC is annotation type 5" -- the table wins, or every PVC
 *             template gets labeled PAC.
 *
 *          5. SUBTYPE INDICES ARE STORED, NOT DERIVED FROM POSITION. Merging
 *             erases an element and shifts everything after it, so a subtype
 *             derived from bank position would silently renumber PVC-2 to
 *             PVC-1 and break comparability against previously archived
 *             features. Assigned once, on first confirmation, then immutable
 *             -- including across pass 2.
 *
 *          6. BANKS ARE PER CHANNEL, and channels are allowed to disagree on
 *             template count. The bin-level monomorphic/polymorphic verdict is
 *             the MAX over ECG channels of the confirmed-PVC-template count
 *             (see polymorphicVerdict). Because only CONFIRMED templates
 *             count, the algorithm can propose a second morphology but cannot
 *             declare polymorphy by itself -- it is an operator-gated finding.
 *
 *          7. LABELS PROPAGATE BY BEAT IDENTITY, not by template index.
 *             ch1's template 2 and ch2's template 2 are different beat sets,
 *             so confirming a beat must label whichever template contains that
 *             beat on EACH channel. propagateLabel() takes a beat index for
 *             this reason and not a template id.
 */

#include <algorithm>
#include <array>
#include "noise_marking_gui/annotation_types.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

 // NO PROJECT INCLUDES ON PURPOSE. Anchors are keyed as int32_t rather than
 // AnchorType so this header pulls in nothing from the tree: template_io.hpp
 // stores anchor tags as int for the same reason ("to keep this header free of
 // the feature_marks dependency"), and the bank has to be storable from there.
 // Cast with static_cast<int32_t>(AnchorType::...) at call sites that have it.

namespace tbank {

    // ---------------------------------------------------------------------
    // Constants
    // ---------------------------------------------------------------------

    // Section 4.6 morphology thresholds, as CORRELATIONS. Used both for
    // assignment and for spawning -- see design note 2.
    // ---- SECTION 4.6 MORPHOLOGY THRESHOLDS, RUNTIME-SETTABLE -------------
    //
    // The spec's defaults. Used for assignment AND for spawning -- see design
    // note 2; there is deliberately no second, looser assignment threshold.
    inline constexpr double kDefaultMatchFloorEcg = 0.85;
    inline constexpr double kDefaultMatchFloorPpg = 0.80;

    // The values actually in force. Set once from config.csv at startup, via
    // setMatchFloors(); read everywhere through the accessors below.
    //
    // WHY NOT constexpr ANY MORE. These were compile-time constants, which
    // meant the only way to try 0.90 was to rebuild. They are the two numbers
    // most likely to need tuning per dataset -- a pulse channel at 0.80 admits
    // far more on a clean arterial line than on a sleep-study pulse-ox -- and
    // the whole partition, every spawn and every merge, turns on them.
    //
    // READ THROUGH THE FUNCTIONS, never by capturing the variable. A caller
    // that copies the value into its own constant at static-init time gets
    // whatever was in force before the config was read, which is the default,
    // silently.
    namespace detail_floors {
        inline double g_ecg = kDefaultMatchFloorEcg;
        inline double g_ppg = kDefaultMatchFloorPpg;
    }

    inline double matchFloorEcg() { return detail_floors::g_ecg; }
    inline double matchFloorPpg() { return detail_floors::g_ppg; }

    // Returns false and changes NOTHING if either value is outside (0, 1].
    //
    // A blank config cell parses to 0.0 through the loader's stod_or_zero, and
    // a floor of 0.0 accepts every beat against every template -- one template
    // per bin, no ectopy ever separated, and no error anywhere to say why. So
    // an unusable value leaves the default in place rather than being applied.
    // A floor above 1.0 is the opposite failure: correlation cannot exceed 1,
    // so every beat spawns and the bank fills with singletons.
    // ---- MINIMUM BEATS FOR A TEMPLATE TO EXIST AT ALL --------------------
    //
    // A template built from too few beats is a median over too few
    // contributors to be a reference for anything. On a real record 44 of 162
    // columns held a single beat and 25 landmark columns held fewer than four,
    // so the operator was being asked to place fiducials on 2-beat waveforms.
    //
    // SEPARATE FROM kMinMembersForColumn, which is the junk/noise CATEGORY
    // gate and stays where it is. This one is set per dataset from config.csv
    // and decides whether the column is written and drawn at all.
    //
    // 0 MEANS NO MINIMUM, and it is the default, so a config without these
    // columns suppresses nothing.
    inline constexpr int kDefaultMinBeatsEcg = 0;
    inline constexpr int kDefaultMinBeatsPpg = 0;

    namespace detail_minbeats {
        inline int g_ecg = kDefaultMinBeatsEcg;
        inline int g_ppg = kDefaultMinBeatsPpg;
    }

    inline int minBeatsEcg() { return detail_minbeats::g_ecg; }
    inline int minBeatsPpg() { return detail_minbeats::g_ppg; }

    // Negative is refused and changes nothing.
    inline bool setMinBeats(int ecg, int ppg) {
        if (ecg < 0 || ppg < 0) return false;
        detail_minbeats::g_ecg = ecg;
        detail_minbeats::g_ppg = ppg;
        return true;
    }

    inline bool setMatchFloors(double ecg, double ppg) {
        if (!(ecg > 0.0 && ecg <= 1.0)) return false;
        if (!(ppg > 0.0 && ppg <= 1.0)) return false;
        detail_floors::g_ecg = ecg;
        detail_floors::g_ppg = ppg;
        return true;
    }

    // Default bank cap. A soft default, not an invariant: a bin whose two
    // closest templates are both confirmed raises its own cap rather than
    // merging them (see TemplateBank::effective_cap).
    inline constexpr int kDefaultMaxTemplatesPerBin = 6;

    // Minimum overlapping non-NaN columns for a correlation to mean anything.
    // alignment.hpp's pearson() returns 0.0 below its own floor, and 0.0 is
    // below every match floor, so an unscorable beat would spawn a template.
    // Beats below this get kUnscorable instead of an assignment.
    inline constexpr int kMinOverlapColumns = 8;

    // A template must hold at least this many members before it earns a
    // display column. Below it the template still exists, still counts, still
    // reaches the archive -- it is just not gridded, because a bad thirty
    // seconds produces many single-beat noise templates and they would drown
    // the record-wide left-to-right reading of drift.
    inline constexpr int kMinMembersForColumn = 2;

    // Merge-eligibility ceiling for the garbage tier, deliberately ABOVE
    // kMinMembersForColumn. Tying the two together is a trap: merging two
    // 1-member noise templates yields a 2-member template, which earns a column
    // and is therefore no longer garbage -- so junk gets promoted out of the
    // tier that exists to collect it, accumulates as 2-member zombies holding
    // slots, and total churn rises rather than falls. A ceiling of 3 lets a
    // merged pair remain collectable once more.
    inline constexpr int kMaxJunkMembers = 3;

    // Members below which a template's own 2.5/97.5 corridor is not an estimate
    // of anything and slot 0's spread is inherited instead (design note 1).
    // Four is the fewest that gives the 2.5 and 97.5 percentiles two distinct
    // values to interpolate between at all; below it they collapse onto the
    // min and max, which is a range, not a percentile.
    inline constexpr int kMinMembersForCorridor = 4;

    // Corridor half-width when even slot 0 has no usable spread at a column --
    // an all-NaN column, or a bin so sparse that slot 0 itself is degenerate.
    // Expressed as a fraction of the template's own peak-to-peak amplitude so
    // it carries no unit and survives normalization. Deliberately generous: a
    // corridor that is too wide admits a beat that should have spawned, which
    // pass 2 can still move; one that is too narrow spawns a template per beat,
    // which nothing downstream can undo.
    inline constexpr double kFallbackCorridorFrac = 0.15;

    // ESTIMATION INFLATION, sqrt(1 + 1/n). A corridor describes the spread of
    // the POPULATION around its median, but it is centred on a median estimated
    // from n members, and a beat is compared against that estimate, not the
    // population. The difference between a fresh draw and an n-member centre has
    // variance sigma^2 * (1 + 1/n), so a corridor sized for sigma alone is too
    // narrow by that factor, and most narrow at exactly n = 1.
    //
    // MEASURED, not derived and hoped for. Held-out beats from a homogeneous
    // population, scored against a 300-member corridor of that same population:
    // min 0.910, mean 0.942 -- comfortably over the 0.85 floor. The same beats
    // against a 1-member template of the same population: min 0.740, mean 0.832,
    // max 0.885. Almost every one FAILS. So without this correction a template
    // cannot acquire its second member no matter how right it is, every beat of a
    // real morphology spawns its own slot, and the bank saturates on one
    // morphology -- 80 spawns and 75 merges on a synthetic bin holding three.
    //
    // At n = 1 the factor is 1.414; by n = 30 it is 1.008 and stops mattering.
    // This is not a second threshold: the floor is untouched and the same
    // comparison is being made. It corrects the corridor's WIDTH for how well
    // its own centre is known, which is a property of the corridor, not of the
    // decision rule.
    inline double corridorInflation(int n_members) {
        if (n_members < 1) return 1.0;
        return std::sqrt(1.0 + 1.0 / static_cast<double>(n_members));
    }

    inline constexpr int kUnscorable = -2;   // returned by assignment
    inline constexpr int kNoMatch = -1;   // spawn required

    inline constexpr uint8_t kUnlabeled = 0;   // label_code sentinel

    // Class codes, READ FROM annotation_types.hpp rather than mirrored.
    //
    // These used to be literals duplicating that table, because the table
    // header pulled in QString and QColor and everything in this namespace --
    // bin_pipeline, seed_pool, nsvt_detect -- runs on the template-generation
    // side, which is Qt-free and should stay that way. The duplication was a
    // known drift risk and annotation_code_check.hpp existed to catch it at
    // startup. It was never called from anywhere, so for as long as the mirror
    // existed the drift was entirely unguarded.
    //
    // annotation_table.hpp is the table with the Qt helpers split off, so there
    // is nothing left to mirror: each constant is the table entry, resolved at
    // compile time by label. A renamed or removed label is a build failure on
    // the static_assert below rather than a 0 that silently becomes "unlabeled",
    // and there is no second copy for anyone to update out of step.
    //
    // Note that the Section 4.6 addendum asserts "PVC is annotation type 5",
    // while the table has 4) PVC at code 4 and 5) PAC at code 5. The table wins,
    // and now it wins structurally: taking the spec's number would mean writing
    // a literal here again.
    inline constexpr uint8_t kCodeMinorNoise = annotation_types::code_for_label("2) Minor Noise");
    inline constexpr uint8_t kCodePvc = annotation_types::code_for_label("4) PVC");
    inline constexpr uint8_t kCodePac = annotation_types::code_for_label("5) PAC");
    inline constexpr uint8_t kCodeCondDelay = annotation_types::code_for_label("6) Cond. Delay");
    inline constexpr uint8_t kCodeAf = annotation_types::code_for_label("7) AF");
    inline constexpr uint8_t kCodeSvt = annotation_types::code_for_label("8) SVT");
    inline constexpr uint8_t kCodeVt = annotation_types::code_for_label("9) VT");

    // The whole of what verifyAnnotationCodes() did, at compile time. 0 is the
    // "no such label" return and is also kUnlabeled, so a missing entry would
    // otherwise turn every template of that class into an unlabeled one.
    static_assert(kCodeMinorNoise != 0, "annotation_types: no \"2) Minor Noise\" entry");
    static_assert(kCodePvc != 0, "annotation_types: no \"4) PVC\" entry");
    static_assert(kCodePac != 0, "annotation_types: no \"5) PAC\" entry");
    static_assert(kCodeCondDelay != 0, "annotation_types: no \"6) Cond. Delay\" entry");
    static_assert(kCodeAf != 0, "annotation_types: no \"7) AF\" entry");
    static_assert(kCodeSvt != 0, "annotation_types: no \"8) SVT\" entry");
    static_assert(kCodeVt != 0, "annotation_types: no \"9) VT\" entry");
    static_assert(kCodePvc != kCodePac, "annotation_types: PVC and PAC share a code");

    // Ventricular origin, for NSVT run detection. PAC, AF and SVT are
    // supraventricular and must not qualify.
    inline bool isVentricular(uint8_t label_code) {
        return label_code == kCodePvc || label_code == kCodeVt;
    }

    // ---------------------------------------------------------------------
    // Per-beat flags. Three independent axes, written for EVERY beat.
    //
    // They are allowed to disagree, and each disagreement is informative:
    //
    //   REGULAR + pvc NONE  + rejected      subtle noise the classifier missed
    //   REGULAR + pvc set   + kept          a PAC the classifier missed
    //   ECTOPIC + pvc NONE  + not eligible  late ventricular beat, or a beat
    //                                       mid-run where the trailing median
    //                                       had already collapsed
    //
    // Because the PVC filter runs across ALL beats, not just category 1, the
    // timing-vs-shape agreement matrix falls out of the archive for free.
    // ---------------------------------------------------------------------

    // Derived from operator marks, so this is always populated: the absence of
    // a mark IS the regular verdict, and there is no moment at which a beat's
    // category is unknown. REGULAR is therefore the default, not a sentinel.
    //
    // Do not confuse this with beat_classifier.hpp's UNKNOWN, which reports
    // that the ONNX model is not wired. The category axis does not depend on
    // that model.
    //
    // NOISE here means the operator's 2) Minor Noise mark. The bank will find
    // additional noise the operator never marked -- artifact that slipped past
    // detection, correlating with nothing and sitting alone in its own
    // template. That is not an unset category; it is classification being
    // wrong in a way the bank surfaces, and the per-beat flags are what make
    // the disagreement visible.
    enum class Category : uint8_t {
        REGULAR = 1,   // unmarked; also 6) Cond. Delay, 7) AF, 8) SVT for now
        ECTOPIC = 2,   // marks 4) PVC, 5) PAC, 9) VT, and the postEligible
        //   beat following each of them
        NOISE = 3    // mark 2) Minor Noise
    };

    // An operator mark code -> the 4.5 category it implies. ONE definition, and
    // it is here because this is where the codes are named: bin_pipeline had its
    // own copy written as bare literals (case 4: case 5: case 9:), which is the
    // same table with nothing tying it to annotation_types, so a renumbered
    // annotation would silently reclassify beats there and not here.
    //
    // Absence of a mark is the REGULAR verdict, not an unknown one -- there is
    // no third state, which is why the default arm returns REGULAR rather than
    // failing.
    inline Category categoryForLabelCode(uint8_t code) {
        if (code == kCodeMinorNoise) return Category::NOISE;
        if (code == kCodePvc || code == kCodePac || code == kCodeVt)
            return Category::ECTOPIC;
        return Category::REGULAR;
    }

    // Per-beat categories from a mark vector, WITH THE POST-ECTOPIC RULE.
    //
    // The beat FOLLOWING an ectopic one inherits ECTOPIC when it carries no mark
    // of its own. That beat is the compensatory pause -- its RR is long, its
    // baseline is still recovering -- and admitting it to the sinus reference
    // contaminates the template with post-ectopic morphology. It is a separate
    // rule from categoryForLabelCode and cannot be derived from a single mark,
    // which is why this takes the whole vector.
    //
    // A MARKED BEAT IS NEVER OVERWRITTEN: the operator's verdict on beat i wins
    // over what beat i-1 implies about it.
    //
    // Moved here from bin_pipeline::categoryFromMarks, which died with
    // runChannel. Its postEligibleCode() also listed AF and SVT, but they only
    // mattered in conjunction with an ECTOPIC verdict on the same mark and both
    // map to REGULAR -- so the set that can actually trigger inheritance is
    // PVC / PAC / VT, and that is what this tests. Same behaviour, one fewer
    // table.
    inline std::vector<Category> categoriesFromMarks(
        const std::vector<uint8_t>& mark_code)
    {
        const size_t n = mark_code.size();
        std::vector<Category> out(n, Category::REGULAR);
        for (size_t i = 0; i < n; ++i) {
            out[i] = categoryForLabelCode(mark_code[i]);
            if (i > 0 && mark_code[i] == 0
                && categoryForLabelCode(mark_code[i - 1]) == Category::ECTOPIC)
                out[i] = Category::ECTOPIC;
        }
        return out;
    }



    // The two paths are recorded separately because their competence is
    // disjoint. The raw flag catches isolated ectopy. The 5-of-8 vote only
    // rescues beats inside CONSECUTIVE runs -- in perfect bigeminy exactly 4
    // of any 8 beats are ectopic, the count never reaches 5, and the vote
    // never fires. High raw / zero vote is bigeminy; high vote is a run.
    enum class PvcFilter : uint8_t {
        NONE = 0,
        PREMATURE = 1,   // RR(t) < 0.80 * median of trailing ten
        VOTE = 2    // 5-of-8 over the prematurity flags
    };

    // Tukey runs INSIDE the classified category-1 population, so most beats
    // were never eligible for it. Recording which pass rejected a beat costs
    // nothing and pays for itself: a REGULAR beat rejected on RR length is a
    // strong hint of ectopy the classifier missed, and that count rising while
    // total rejections fall is the high-burden failure announcing itself.
    enum class TukeyOutcome : uint8_t {
        NOT_ELIGIBLE = 0,   // not category 1
        KEPT = 1,
        REJ_RR_LENGTH = 2,
        REJ_AMPLITUDE = 3,
        REJ_R_LOCATION = 4,
        REJ_WAVE_SCORE = 5
    };

    struct BeatFlags {
        Category     category = Category::REGULAR;
        PvcFilter    pvc = PvcFilter::NONE;
        TukeyOutcome tukey = TukeyOutcome::NOT_ELIGIBLE;

        // Bank assignment per ECG channel, and PPG. -1 = unassigned,
        // kUnscorable = too little axis overlap to score.
        std::array<int32_t, 3> template_id_ecg = { -1, -1, -1 };
        int32_t                template_id_ppg = -1;

        // True when this beat's waveform is an EWMA substitution rather than
        // an observation (4.6 beat substitution, alpha = 1/8, run AFTER
        // assignment against the assigned template's own average). Variance
        // measures must not treat these as observations.
        bool substituted = false;
    };

    // ---------------------------------------------------------------------
    // Per-bin counts. This defect class produces output that looks fine and
    // shows up only in the counts, so they are first-class state, not
    // debug logging.
    //
    // Tukey assumes one dominant population. With ectopy in the minority both
    // quartiles sit in the sinus cluster and every ectopic beat falls outside
    // the fence. As ectopy approaches half the bin the IQR widens across both
    // clusters and Tukey stops rejecting anything, including real artifact.
    // Moving Tukey inside the classified pool does not remove that failure --
    // it makes Tukey's correctness CONDITIONAL ON THE CLASSIFIER'S RECALL,
    // one layer deeper and quieter, with no backstop behind it. The two
    // mechanisms also degrade together rather than independently: at high
    // burden the trailing-ten median is dragged down by ectopic RRs, so
    // RR(t) < 0.80*median fires less often at exactly the burden where the
    // fences have widened past usefulness.
    //
    // Hence fence_iqr alongside the rejection counts. Counts give the
    // outcome; the IQR gives the cause, and it is continuous where the counts
    // are discrete. The signature to watch for is LOW REJECTION RATE WITH
    // WIDE IQR -- indistinguishable from "the bin is genuinely clean" if you
    // only have the counts.
    // ---------------------------------------------------------------------

    struct TukeyPassCounts {
        uint32_t beats_in = 0;
        uint32_t beats_out = 0;
        uint32_t rejected = 0;
        double   q1 = std::numeric_limits<double>::quiet_NaN();
        double   q3 = std::numeric_limits<double>::quiet_NaN();
        double   fence_lo = std::numeric_limits<double>::quiet_NaN();
        double   fence_hi = std::numeric_limits<double>::quiet_NaN();
        double   iqr() const { return q3 - q1; }
    };

    struct BinCounts {
        uint32_t beats_detected = 0;

        // Classification, before anything is removed.
        uint32_t n_regular = 0;
        uint32_t n_ectopic = 0;
        uint32_t n_noise = 0;

        // PVC filter over ALL beats, paths kept separate.
        uint32_t n_premature = 0;
        uint32_t n_vote_only = 0;   // voted but not itself premature

        // Tukey inside category 1, one entry per pass.
        TukeyPassCounts tukey_rr, tukey_amplitude, tukey_r_location,
            tukey_wave_score;

        // Category-1 beats Tukey rejected on RR length: the cheapest estimate
        // of classifier recall you will get. A sinus-labeled beat rejected for
        // being short is very likely ectopy that was missed.
        uint32_t n_regular_rejected_on_rr = 0;

        // Bank churn. Merges are mostly garbage collection of lone noise
        // templates, so this doubles as a noise metric: heavy merge activity
        // means a high noise fraction in the bin.
        uint32_t n_spawns = 0;
        uint32_t n_merges = 0;
        // Merges split by tier. n_merges_garbage collected two templates that
        // both sat below kMinMembersForColumn; n_merges_real absorbed at least
        // one template with real membership, which is the case worth watching --
        // a bin with many real merges lost morphologies to cap pressure and its
        // template count understates what was there.
        uint32_t n_merges_garbage = 0;
        uint32_t n_merges_real = 0;
        uint32_t n_cap_raises = 0;
        uint32_t n_unscorable = 0;
        uint32_t n_reassigned_pass2 = 0;

        // --- how slot 0 was obtained ------------------------------------
        // Exactly one of the first two is 1. seed_fallback_used == 1 means no
        // Phase 1 template reached this bin and slot 0 is a locally rebuilt
        // seed-pool median, so this bin does NOT satisfy Section 4.6's first
        // bullet and its features are measured against a different reference
        // from Phase 1's. seed_verified_sinus == 0 with seed_from_phase1 == 1
        // means the template was supplied but the caller could not vouch that
        // it excludes ectopy -- which, given create_ecg_templates.hpp applies no
        // rhythm test, is the expected state today.
        uint32_t seed_from_phase1 = 0;
        uint32_t seed_fallback_used = 0;
        uint32_t seed_verified_sinus = 0;

        // NOT SERIALIZED, deliberately. writeCounts() lays BinCounts out
        // field-by-field with writeSeed() immediately after it inside the bank
        // block, so appending here shifts everything downstream and needs a
        // kVersionBank -> v2 bump plus a version branch in readCounts(). These
        // three are build-time diagnostics and surface in the per-bin stderr
        // log. Persist them with the next format bump, not before one.
    };

    // ---------------------------------------------------------------------
    // One template in a bank
    // ---------------------------------------------------------------------

    // Marker positions for ONE (channel, template) pair. TemplateBin's
    // MarkerSet carries arrays of 3 because it predates per-channel banks;
    // these are scalar, because a bank member belongs to exactly one channel.
    // Per-template marker sets are mandatory, not a refinement: a PVC's
    // Q-onset is at a different column than sinus's, so they cannot share.
    struct BankMarkerSet {
        // FOUR BARS AND ONE STORED GLYPH.
        //
        // p_begin / q_begin / s_end / t_end are the bars: the four positions
        // the operator drags, each belonging to exactly one alignment (see
        // anchor_view.hpp). p_peak is the one glyph kept here rather than
        // recomputed at every use, because it is written to both the markings
        // CSV and the markings bin and downstream consumers read it from the
        // file. It is REFRESHED FROM THE BARS, never detected independently --
        // FeatureMarks::reactive_ecg brackets it with p_begin and q_begin, and
        // TemplateBin::syncReactiveGlyphs is the only thing that assigns it. A
        // second, detector-sourced answer stored here is what made the
        // on-screen X and the CSV column disagree.
        //
        // t_begin REMOVED, RECORD AND ALL. It was a marker field nothing set
        // and nothing drew: maskFor had no entry for it, seed_all seeded
        // q_begin / s_end / t_end / p_begin only, and markerAtX never
        // hit-tested it -- so it sat at -1 for the life of every template while
        // the CSV's t_peak_*_user column bracketed T-peak against it and
        // therefore reported nothing, even though the on-screen X was bracketed
        // by s_end and t_end and sat in the right place.
        //
        // Both binary records drop the field rather than reserving its four
        // bytes, so FILES WRITTEN BEFORE THIS CHANGE DO NOT PARSE: neither
        // _template_markings.bin nor tbank_ser::detail::writeMarkerSet carries
        // a version, so there is nothing to branch on. Regenerate templates and
        // re-mark; do not attempt to read an old pair.
        //
        // The T-wave onset that morphology_envelope, premark_beats and
        // beat_classifier band on is a DIFFERENT quantity, computed from the
        // signal, and is untouched by this removal.
        int p_begin = -1;
        int p_peak = -1;
        int q_begin = -1;
        int s_end = -1;
        int t_end = -1;

        bool isUnset() const {
            return p_begin < 0 && p_peak < 0 && q_begin < 0
                && s_end < 0 && t_end < 0;
        }
    };
    struct BankPulseMarkerSet {
        int onset = -1, peak = -1, dicrotic = -1, peak2 = -1, end = -1;
        int t50 = -1, t80 = -1;
        double onset_auto = -1.0, peak_auto = -1.0, dicrotic_auto = -1.0;
        double peak2_auto = -1.0, end_auto = -1.0;
        bool notch_found = false;
        bool isUnset() const {
            return onset < 0 && peak < 0 && dicrotic < 0 && peak2 < 0 && end < 0;
        }
    };

    struct BankTemplate {
        // Column-wise NaN-skipping median over members, on the bin's shared
        // axis. Recomputed whenever membership changes (design note 3).
        std::vector<double> tmpl;
        std::vector<double> tmpl_iqr;      // per-sample spread
        int                 r_col = -1;

        // Bin-local beat indices. Needed for three things the spec requires
        // and the addendum's TemplateBank struct cannot express: recomputing
        // the median, propagating a label to "every other beat assigned to
        // the same template", and cross-channel label propagation by beat
        // identity.
        std::vector<uint32_t> members;

        // 0 = unlabeled. Unlabeled means NOT YET CONFIRMED, not unknown
        // class -- a distinct state, and the difference carries information.
        // Never infer this from morphological similarity to a labeled
        // template: whether two close morphologies are one class or two IS
        // the finding the operator is being asked to produce.
        uint8_t label_code = kUnlabeled;

        // CONFIRMATION IS A SEPARATE FIELD FROM THE CLASS. label_code may hold a
        // PRESUMED class -- an unmarked template presents as PQRST because that
        // is the sensible default, not because anyone said so. Only a confirmed
        // label is usable as training data later, and a training set that
        // silently included presumptions would be training on the algorithm's
        // own guesses with no way for a consumer to tell.
        //
        // So `confirmed` is never inferred from label_code being set. It is set
        // only by propagateLabel(), i.e. only when an operator confirmed a beat.
        bool confirmed_by_operator = false;

        // Appended by the bank in order of first appearance within the class:
        // PVC-1, PVC-2. -1 until the first confirmation. The operator never
        // types this and never sees it until the bank produces it.
        int32_t subtype = -1;

        // Spawn order within the bin, assigned in pass 1 and NOT recomputed
        // in pass 2. This is what "order of first appearance" resolves
        // against, so that reassignment cannot renumber subtypes between
        // runs.
        uint32_t spawn_seq = 0;

        // PPG members of THIS template, when the template is a projection of a
        // joint (four-channel) group. -1 = unknown, which is what a
        // per-channel-built bank leaves it as.
        //
        // Exists because the pulse cohort is a property of the GROUP, not of the
        // bin. Every column of a bin used to report the same bin-wide
        // ppg_n_beats, which made a split look like it had ignored PPG
        // entirely. A group's PPG members are the members it actually has on
        // that channel, and they differ between siblings.
        int32_t n_ppg_members = -1;

        // Members a Tukey pass rejected, counted where the flags are indexed
        // the same way the group's membership is -- in the joint projection,
        // over SLICES. It cannot be recounted downstream: BankTemplate::members
        // is in the channel's own local row space and the flag vector is in
        // slice space, so indexing one with the other silently reads another
        // beat's verdict. Which is exactly what the aggregate helpers in
        // morphology_csv used to do.
        int32_t n_tukey_members = 0;

        // ---- THE OPERATOR'S QUALITY VERDICT ON THIS PANEL -----------------
        //
        // 0 = good, 1 = bad R detection, 2 = bad pulse. Mirrors
        // BinPlotWidget::State, and it is PER TEMPLATE because a panel is a
        // (bin, template) pair.
        //
        // TemplateBin::bad_r_ch[3] and bad_ppg already existed and are per BIN.
        // That was right when a bin was one panel; a bin now occupies one panel
        // per markable morphology, so recording a right-click against the bin
        // crossed out every panel of that bin -- including morphologies the
        // operator had never looked at, and after a page rebuild rather than at
        // the moment of the click.
        //
        // THE BIN-LEVEL FLAGS STAY, AND SLOT 0 OWNS THEM. NormalizeFeatures
        // skips a whole bin's channel on bad_r_ch, template_marking_bin_io
        // serializes and exports it, and feature_marks sets it automatically --
        // so it has to keep meaning "this bin's lead is untrustworthy". Marking
        // the SEED panel sets it; marking a sub-template records the verdict
        // here and does not. Otherwise one bad 2-beat junk column would exclude
        // an entire bin from the feature reference, which is a much worse
        // outcome than the one being fixed.
        uint8_t operator_state = 0;

        // ---- TWO MEMBER LISTS, AND BOTH ARE REQUIRED --------------------
        //
        // `members` is EVERY beat the partition assigned to this template. It
        // is what the archive writes: a premature or Tukey-rejected beat is
        // still this template's beat, and the output has to say so and say why.
        // Dropping it from here would delete the record instead of marking it.
        //
        // `members_clean` is members minus the premature and minus the
        // Tukey-rejected. It is what the averaged waveform is built from, what
        // is drawn and marked on screen, and what "kept" means in the outputs.
        //
        // Empty means "not yet computed" -- the post-partition stage in
        // bin_pipeline fills it -- and consumers should fall back to `members`
        // in that case rather than treat it as an empty template.
        //
        // The exclusion REASON is not stored here. It is per beat, on
        // BeatFlags::pvc and BeatFlags::tukey, and duplicating it per template
        // would give two places to disagree about the same beat.
        std::vector<uint32_t> members_clean;

        int cleanCount() const {
            return static_cast<int>(members_clean.empty()
                ? members.size() : members_clean.size());
        }
        int excludedCount() const {
            return static_cast<int>(members_clean.empty()
                ? 0 : members.size() - members_clean.size());
        }

        // TOO FEW BEATS TO BE A TEMPLATE, against the configured minimum for
        // this channel kind (tbank::minBeatsEcg / minBeatsPpg). Reported in the
        // archive and used to suppress the column entirely.
        //
        // ON cleanCount(), not memberCount(): the question is how many beats
        // are actually behind the drawn waveform, and premature or
        // Tukey-rejected members are not.
        //
        // A zero minimum returns false for everything, which is the default --
        // nothing is suppressed unless the operator configured a threshold.
        bool tooFewBeats(bool is_ppg) const {
            const int lim = is_ppg ? minBeatsPpg() : minBeatsEcg();
            return lim > 0 && cleanCount() < lim;
        }

        // Landmarks per alignment anchor. A ventricular template has no P
        // wave, so p_begin/p_peak stay -1 legitimately and every P-dependent
        // feature must come out NaN rather than 0. Downstream extraction has
        // to treat -1 here as a valid state, not an error, or it will report
        // a PR interval measured from a P wave that does not exist.
        std::map<int32_t, BankMarkerSet> markers_by_anchor;   // key: AnchorType
        BankPulseMarkerSet pulse_marks;   // meaningful only on ppg_bank slots
        bool hasDetectedPulseMarks() const { return !pulse_marks.isUnset(); }

        // NOTE: INSERTS. markers_by_anchor[a] default-constructs an all -1 set
        // when the key is absent, so a read through this overload is a write.
        // Use hasDetectedMarks() to ask whether a set is populated.
        BankMarkerSet& marks(int32_t a) { return markers_by_anchor[a]; }

        // Non-inserting, and asks about CONTENT rather than key presence.
        bool hasDetectedMarks(int32_t a) const {
            auto it = markers_by_anchor.find(a);
            return it != markers_by_anchor.end() && !it->second.isUnset();
        }
        const BankMarkerSet& marks(int32_t a) const {
            static const BankMarkerSet kEmpty;
            auto it = markers_by_anchor.find(a);
            return (it == markers_by_anchor.end()) ? kEmpty : it->second;
        }

        bool confirmed() const { return confirmed_by_operator; }

        // Per-member census, filled after assignment settles. These are what
        // make a template's PRESUMED category computable before any operator
        // mark exists: prematurity is a timing verdict available at build time,
        // and non-reproducibility shows up as a member count that never grew.
        uint32_t n_premature_members = 0;
        uint32_t n_voted_members = 0;
        uint32_t n_noise_members = 0;

        // Presumed, NOT confirmed. Decides which templates go in front of an
        // operator for landmark marking: only category 1 templates get
        // landmarks, because only category 1 beats feed feature extraction -- a
        // P-onset on a PVC template has nothing downstream that consumes it,
        // and a PVC's QT is not comparable to a sinus QT.
        //
        // UNLABELED IS PQRST, NOT A THIRD STATE. Section 4.6: "A template with
        // no confirmed member stays unlabeled. Unlabeled means not yet
        // confirmed, NOT unknown class." That second sentence forbids exactly
        // what an UNCONFIRMED category would be -- an unknown bucket that gets
        // withheld from display. A template is displayed as PQRST until an
        // operator marks it otherwise; being unconfirmed is a statement about
        // the operator's progress, not about the morphology.
        //
        // What the clause DOES forbid is inferring a label from morphological
        // similarity to a labeled template, and that prohibition is enforced
        // where the inference actually happened: mergeTemplates() no longer
        // copies a confirmed label onto the template it absorbs, and
        // findMergePair() blocks any pair containing a confirmed template. That
        // is the operation the spec names. This function is not it.
        //
        // A confirmed label overrides the presumption outright; the operator's
        // verdict is not a hypothesis to be re-derived.
        Category presumedCategory() const {
            if (confirmed_by_operator) {
                if (label_code == kCodeMinorNoise) return Category::NOISE;
                if (label_code == kCodePvc || label_code == kCodePac
                    || label_code == kCodeVt) return Category::ECTOPIC;
                return Category::REGULAR;
            }
            // Never reproduced. A member-count test, not a class inference:
            // it says nothing about what the beat was, only that one beat is
            // not a morphology. This is the gate that keeps single-beat noise
            // fragments out of the marking grid.
            if (!earnsColumn()) return Category::NOISE;
            const uint32_t n = static_cast<uint32_t>(members.size());
            if (n == 0) return Category::NOISE;
            const uint32_t ect = n_premature_members + n_voted_members;
            // Majority-premature membership is the timing evidence for an
            // ectopic morphology. Deliberately a majority and not any: a sinus
            // template picks up the occasional premature beat, and one such beat
            // must not reclassify 800 others.
            if (ect * 2 > n) return Category::ECTOPIC;
            return Category::REGULAR;
        }

        // The Phase 1 sinus seed. spawn_seq 0 is issued once, to slot 0, by the
        // pipeline before pass 1 begins, and mergeTemplates() keeps the LOWER
        // spawn_seq of a merged pair, so this survives merge history without a
        // stored flag -- which also keeps the serialized layout byte-identical.
        bool isSeed() const { return spawn_seq == 0; }

        bool wantsLandmarkMarking() const {
            return presumedCategory() == Category::REGULAR;
        }
        // Per-column 2.5/97.5 corridor over the members, and the flag saying
        // whether it is the template's OWN spread or slot 0's inherited one.
        // Derived from members exactly as tmpl is, so neither is serialized --
        // template_bank_serialize.hpp stores nothing derived, on purpose.
        std::vector<double> band_lo;
        std::vector<double> band_hi;
        bool corridor_inherited = false;

        int  memberCount() const { return static_cast<int>(members.size()); }
        bool earnsColumn() const { return memberCount() >= kMinMembersForColumn; }

        // Eligible for the garbage-collection merge tier. Not the negation of
        // earnsColumn() -- see kMaxJunkMembers.
        bool isJunk() const { return memberCount() <= kMaxJunkMembers; }
    };

    // ---------------------------------------------------------------------
    // The bank: one per channel per bin
    // ---------------------------------------------------------------------

    struct TemplateBank {
        std::vector<BankTemplate> templates;

        // Per-bin, because the confirmed-member rule raises it. The addendum's
        // assignToTemplate() takes maxTemplates by value and returns a single
        // int, so it can neither persist a raise nor report one for logging;
        // the cap lives here instead.
        // max_templates_per_bin. A CONFIG VALUE with a default, which is what
        // the spec says it is: "Cap the bank at max_templates_per_bin (default
        // 6)". It was previously reachable only by editing
        // kDefaultMaxTemplatesPerBin and recompiling, so the "default" was the
        // only value the program had. Set from cfg at pipeline entry
        // (BinInput::max_templates_per_bin); the constant is now only the
        // fallback when nothing supplies one.
        //
        // configured_cap is the value the operator asked for and never changes.
        // effective_cap starts equal to it and only ever RISES, by the
        // confirmed-member rule. Keeping both means a raise is visible as a
        // difference rather than having to be reconstructed from the log.
        int32_t configured_cap = kDefaultMaxTemplatesPerBin;
        int32_t effective_cap = kDefaultMaxTemplatesPerBin;

        void setCap(int32_t cap) {
            if (cap < 1) cap = 1;   // slot 0 always exists
            configured_cap = cap;
            effective_cap = cap;
        }

        uint32_t next_spawn_seq = 0;

        // Total beats routed into this bank, for beat-share reporting.
        uint32_t assigned_beats = 0;

        double matchFloorFor(bool is_ppg) const {
            return is_ppg ? matchFloorPpg() : matchFloorEcg();
        }

        int  size() const { return static_cast<int>(templates.size()); }
        bool atCap() const { return size() >= effective_cap; }

        int nConfirmed() const {
            int c = 0;
            for (const auto& t : templates) if (t.confirmed()) ++c;
            return c;
        }

        // Distinct templates confirmed as `code`, i.e. the highest subtype
        // index reached for that class. One PVC template is monomorphic; two
        // or more is polymorphic.
        int countLabeled(uint8_t code) const {
            int c = 0;
            for (const auto& t : templates) if (t.label_code == code) ++c;
            return c;
        }

        // Next subtype index for a class, by order of first appearance. Reads
        // the highest index already issued rather than counting current
        // members, so a merge that removes a labeled template cannot cause
        // the next one to reuse a retired index.
        int32_t nextSubtypeFor(uint8_t code) const {
            int32_t hi = 0;
            for (const auto& t : templates)
                if (t.label_code == code) hi = std::max(hi, t.subtype);
            return hi + 1;
        }

        int findByBeat(uint32_t beat_idx) const {
            for (int i = 0; i < size(); ++i)
                for (uint32_t m : templates[i].members)
                    if (m == beat_idx) return i;
            return -1;
        }

        double beatShare(int template_idx) const {
            if (template_idx < 0 || template_idx >= size() || assigned_beats == 0)
                return std::numeric_limits<double>::quiet_NaN();
            return static_cast<double>(templates[template_idx].memberCount())
                / static_cast<double>(assigned_beats);
        }
    };

    // ---------------------------------------------------------------------
    // Cap-raise events. "Log the event" is a requirement, not bookkeeping: a
    // bin that raised its cap has three or more confirmed morphologies
    // competing for six slots, which makes it the most clinically interesting
    // bin in the record and the one most likely to expose a bug.
    // ---------------------------------------------------------------------

    struct CapRaiseEvent {
        uint64_t bin_index = 0;
        int      channel = 0;
        int32_t  old_cap = 0;
        int32_t  new_cap = 0;
        int32_t  template_a = -1;   // the two closest, both confirmed
        int32_t  template_b = -1;
        double   closeness = 0.0;   // r between them
        uint8_t  label_a = kUnlabeled;
        uint8_t  label_b = kUnlabeled;
    };

    // ---------------------------------------------------------------------
    // Bin-level verdict
    // ---------------------------------------------------------------------

    struct PolymorphicVerdict {
        int  count = 0;      // max over channels of confirmed PVC templates
        int  driving_channel = -1;     // which channel produced the max
        std::array<int, 3> per_channel = { 0, 0, 0 };  // stored for audit
        bool polymorphic() const { return count >= 2; }
        bool monomorphic() const { return count == 1; }
    };

    // Max across ECG channels. Per-channel counts are kept alongside the max
    // so the verdict is auditable and can be recomputed under a different
    // rule later without rerunning the marking. PPG is excluded on purpose: a
    // split in the PPG bank is far more likely perfusion or motion than a
    // second ventricular focus, and it is reported separately.
    inline PolymorphicVerdict polymorphicVerdict(const TemplateBank& ch1,
        const TemplateBank& ch2,
        const TemplateBank& ch3,
        uint8_t code = kCodePvc)
    {
        const TemplateBank* bk[3] = { &ch1, &ch2, &ch3 };
        PolymorphicVerdict v;
        for (int c = 0; c < 3; ++c) {
            v.per_channel[c] = bk[c]->countLabeled(code);
            if (v.per_channel[c] > v.count) {
                v.count = v.per_channel[c];
                v.driving_channel = c;
            }
        }
        return v;
    }

    inline PolymorphicVerdict polymorphicVerdict(
        const std::array<TemplateBank, 3>& ecg_banks, uint8_t code = kCodePvc)
    {
        return polymorphicVerdict(ecg_banks[0], ecg_banks[1], ecg_banks[2], code);
    }

    // ---------------------------------------------------------------------
    // Label propagation
    //
    // Operator confirms ONE beat -> the class label attaches to the template
    // that beat is assigned to -> and from there to every other beat in that
    // template. One click labels a whole morphology cluster.
    //
    // Keyed on the beat, not the template, so a single confirmation reaches
    // every channel's bank correctly even though the banks disagree on
    // template indices.
    // ---------------------------------------------------------------------

    struct PropagationResult {
        std::array<int32_t, 3> labeled_template = { -1, -1, -1 };
        std::array<int32_t, 3> subtype = { -1, -1, -1 };
        int beats_relabeled = 0;

        // The PPG bank's outcome, kept separate from the ECG array rather than
        // appended to it: a fourth entry in `labeled_template` would be read as
        // a fourth ECG lead by anything that iterates the array, and several
        // things do.
        int32_t ppg_labeled_template = -1;
        int32_t ppg_subtype = -1;
        int     ppg_beats_relabeled = 0;
    };

    // ---- SLOT-KEYED PROPAGATION, FOR THE JOINT PARTITION -----------------
    //
    // Use THIS one, not the beat-keyed pair below, on any bank that came out of
    // jbank::projectToChannel.
    //
    // The beat-keyed versions exist because three independently built banks
    // disagreed about which slot a morphology occupied, so a beat index was the
    // only shared handle. Under the joint partition template i IS group i on
    // every channel -- projectToChannel walks bank.groups in order, and
    // copyBanks preserves that order into BinTemplates -- so the slot is itself
    // the shared handle, and the beat lookup is not just unnecessary but wrong:
    // findByBeat searches `members`, which the projection fills with CHANNEL
    // LOCAL ROW indices. A representative beat's row number on CH1 means a
    // different beat on CH2 and a different one again on PPG, so keying on it
    // labels whichever unrelated morphology happens to hold that row number.
    //
    // AND IT REACHES PPG. The beat-keyed three-bank overload never touched the
    // pulse channel at all, so a confirmed PVC left its pulse cohort unlabeled.
    inline PropagationResult propagateLabelBySlot(
        std::array<TemplateBank, 3>& ecg_banks,
        TemplateBank& ppg_bank,
        int slot,
        uint8_t label_code)
    {
        PropagationResult out;
        if (label_code == kUnlabeled || slot < 0) return out;

        auto label = [&](TemplateBank& bk, int32_t* which, int32_t* sub,
            int* relabeled) {
                if (slot >= bk.size()) return;
                BankTemplate& t = bk.templates[slot];
                // Subtype is issued once and then immutable: a re-confirmation,
                // or a confirmation of another beat in an already-labeled
                // template, must not mint a second index for one class.
                if (t.subtype < 0 || t.label_code != label_code) {
                    if (t.label_code != label_code)
                        t.subtype = bk.nextSubtypeFor(label_code);
                    t.label_code = label_code;
                }
                t.confirmed_by_operator = true;
                if (which) *which = slot;
                if (sub) *sub = t.subtype;
                if (relabeled) *relabeled += t.memberCount();
            };

        for (int c = 0; c < 3; ++c)
            label(ecg_banks[c], &out.labeled_template[c], &out.subtype[c],
                &out.beats_relabeled);
        label(ppg_bank, &out.ppg_labeled_template, &out.ppg_subtype,
            &out.ppg_beats_relabeled);
        return out;
    }

    inline PropagationResult propagateLabel(std::array<TemplateBank, 3>& banks,
        uint32_t beat_idx,
        uint8_t label_code)
    {
        PropagationResult out;
        if (label_code == kUnlabeled) return out;

        for (int c = 0; c < 3; ++c) {
            TemplateBank& bk = banks[c];
            const int ti = bk.findByBeat(beat_idx);
            if (ti < 0) continue;

            BankTemplate& t = bk.templates[ti];

            // Subtype is issued once and then immutable. A re-confirmation, or
            // a confirmation of a different beat in an already-labeled
            // template, must not mint a new index.
            if (t.subtype < 0 || t.label_code != label_code) {
                if (t.label_code != label_code)
                    t.subtype = bk.nextSubtypeFor(label_code);
                t.label_code = label_code;
            }
            t.confirmed_by_operator = true;

            out.labeled_template[c] = ti;
            out.subtype[c] = t.subtype;
            out.beats_relabeled += t.memberCount();
        }
        return out;
    }

    // ---------------------------------------------------------------------
    // ONE BANK, for PPG. The rule -- "the class label propagates to the
    // template that beat is assigned to, and from there to every other beat
    // assigned to the same template" -- says nothing about ECG, and the PPG
    // bank obeys it or it has confirmed members that no operator ever
    // confirmed. Split out of the three-bank version rather than duplicated:
    // subtype issuance and the confirmation flag are the two things that must
    // not diverge between channel types.
    //
    // TAKES A PPG BEAT INDEX, WHICH IS NOT AN ECG BEAT INDEX. The PPG kept set
    // is filtered independently, so beat k of CH1 and beat k of PPG are
    // different beats. The caller has to map the operator's click into PPG
    // space before calling this; passing an ECG index would label a pulse
    // morphology chosen at random.
    // ---------------------------------------------------------------------
    inline void propagateLabelToBank(TemplateBank& bank,
        uint32_t beat_idx,
        uint8_t label_code,
        int32_t* out_template = nullptr,
        int32_t* out_subtype = nullptr,
        int* out_beats = nullptr)
    {
        if (label_code == kUnlabeled) return;
        const int ti = bank.findByBeat(beat_idx);
        if (ti < 0) return;

        BankTemplate& t = bank.templates[ti];

        // Subtype is issued once and then immutable, exactly as in the ECG
        // path: a re-confirmation, or a confirmation of a different beat in an
        // already-labeled template, must not mint a new index.
        if (t.subtype < 0 || t.label_code != label_code) {
            if (t.label_code != label_code)
                t.subtype = bank.nextSubtypeFor(label_code);
            t.label_code = label_code;
        }
        t.confirmed_by_operator = true;

        if (out_template) *out_template = ti;
        if (out_subtype)  *out_subtype = t.subtype;
        if (out_beats)    *out_beats = t.memberCount();
    }

    // Three ECG banks plus the PPG bank, for a confirmation the caller has
    // already mapped into both index spaces. ecg_beat_idx and ppg_beat_idx are
    // separate parameters precisely because they are separate index spaces --
    // see propagateLabelToBank.
    inline PropagationResult propagateLabel(std::array<TemplateBank, 3>& banks,
        TemplateBank& ppg_bank,
        uint32_t ecg_beat_idx,
        uint32_t ppg_beat_idx,
        uint8_t label_code)
    {
        PropagationResult out = propagateLabel(banks, ecg_beat_idx, label_code);
        propagateLabelToBank(ppg_bank, ppg_beat_idx, label_code,
            &out.ppg_labeled_template, &out.ppg_subtype,
            &out.ppg_beats_relabeled);
        return out;
    }

    // ---------------------------------------------------------------------
    // Merge candidate selection
    //
    // Never merge two templates that both have confirmed members, even at the
    // cap. Merging them does not lose a beat -- it changes a clinical
    // finding, silently, because the merged template looks perfectly normal
    // afterwards. It collapses polymorphic ectopy into monomorphic, which is
    // the one distinction this section exists to preserve. Given that the
    // bin-level verdict counts only CONFIRMED templates, such a merge would
    // destroy the sole evidence for polymorphy.
    //
    // TIER ORDER, AND WHY IT DEPARTS FROM "MERGE THE TWO CLOSEST".
    //
    // The spec says merge the two closest templates. Taken literally on real
    // data that rule inverts its own intent. Measured on one record, bin 0,
    // CH1: T0 was sinus (996 members), T5 was a real second morphology (152
    // members, r = 0.750 against sinus), and T1-T4 were single-beat noise
    // templates correlating with everything at 0.1-0.2. The CLOSEST pair in
    // that bank is T0 and T5 -- the two real morphologies -- because noise,
    // being non-reproducible, is by construction the FURTHEST from everything.
    // So "closest" merges the finding and preserves the junk, and it does so
    // every time the cap saturates: that bin logged 34 spawns and 29 merges,
    // and a neighbouring bin 72 and 67.
    //
    // Hence tier 1: templates that do not earn a column (below
    // kMinMembersForColumn) are merged first, closest pair among them. That is
    // garbage collection, and the product correlates with nothing so it will
    // not steal real beats afterwards. Only when no such pair exists does
    // tier 2 fall back to the spec's rule over the whole bank.
    //
    // This is the spec's own reasoning applied one level down. It forbids a
    // merge that would "silently collapse polymorphic ectopy into
    // monomorphic", and protects confirmed templates for that reason -- but
    // confirmation happens during review, long after the bank is built, so at
    // build time nothing is confirmed and nothing is protected. Member count
    // is the only evidence available at that point about which templates are
    // real.
    //
    // Unlabeled-to-unlabeled and unlabeled-to-confirmed merges remain allowed.
    // ---------------------------------------------------------------------

    struct MergeCandidate {
        int    a = -1;
        int    b = -1;
        double closeness = -std::numeric_limits<double>::infinity();
        bool   both_confirmed = false;   // true => raise the cap instead

        // MISNAMED, DELIBERATELY NOT RENAMED YET. This reads true whenever the
        // closest pair contained a confirmed template, not only when BOTH were
        // confirmed. findMergePair blocks any pair touching a confirmed
        // template -- stricter than 4.6's literal "both", because merging a
        // confirmed template into an unconfirmed one would absorb unlabelled
        // beats into a labelled class, which is the label-by-similarity
        // inference the section forbids elsewhere. The consequence is more cap
        // raises than the literal rule produces, which is the safe direction.
        //
        // Anyone reading a CapRaiseEvent should treat this as
        // "confirmed_involved". Renaming it touches template_assign.hpp,
        // bin_pipeline.hpp and joint_bank.hpp together and is left for a change
        // that can be compiled and run, not folded into a comment fix.
        bool   garbage_pair = false;   // both below the column threshold
        bool   valid() const { return a >= 0 && b >= 0; }
    };

    // `sim` returns the closeness of two templates. Returns the closest
    // MERGEABLE pair; if the closest pair is blocked by the confirmed-member
    // rule, it comes back with both_confirmed set so the caller raises the cap
    // and logs, rather than merging.
    //
    // THE SPEC'S RULE, SINGLE TIER: merge the two closest. The junk-first tier
    // that used to sit in front of it is gone. What it was for is on record and
    // is not hypothetical -- measured on one record, bin 0, CH1: T0 sinus (996
    // members), T5 a real second morphology (152 members, closeness 0.750
    // against sinus), T1-T4 single-beat noise templates scoring 0.1-0.2 against
    // everything. The closest pair in that bank is T0 and T5, the two REAL
    // morphologies, because noise is non-reproducible and therefore furthest
    // from everything. So "closest" merges the finding and preserves the junk,
    // and it does so every time the cap saturates: that bin logged 34 spawns
    // and 29 merges, a neighbouring one 72 and 67.
    //
    // That behaviour is restored deliberately, on instruction, because the rule
    // is the spec's. n_merges_garbage and n_merges_real still classify every
    // merge, so the damage is counted rather than argued about: a bin with many
    // REAL merges lost morphologies to cap pressure and its template count
    // understates what was there. Watch that counter.
    //
    // WHY ANY PAIR CONTAINING A CONFIRMED TEMPLATE IS BLOCKED, not only
    // confirmed-confirmed. The spec blocks confirmed-confirmed (rule 7) and
    // separately forbids inferring a label (rule 6). A confirmed-into-unlabeled
    // merge satisfies the first and breaks the second: labels are per template
    // and propagate to every beat in it, so absorbing an unlabeled template into
    // a confirmed PVC labels its beats PVC on the evidence of a correlation --
    // which is the operator's judgment, made by the algorithm. The two clauses
    // together therefore imply a stronger restriction than either states, and
    // this is it. Blocked pairs raise the cap, which rule 7 already establishes
    // as the correct response to an unmergeable bank.
    template <class SimFn>
    inline MergeCandidate findMergePair(const TemplateBank& bank, SimFn sim)
    {
        MergeCandidate best_any, best_blocked;
        const int n = bank.size();
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const double s = sim(bank.templates[i], bank.templates[j]);
                if (std::isnan(s)) continue;
                if (bank.templates[i].confirmed() || bank.templates[j].confirmed()) {
                    if (s > best_blocked.closeness) {
                        best_blocked.a = i; best_blocked.b = j;
                        best_blocked.closeness = s;
                        best_blocked.both_confirmed = true;
                    }
                    continue;
                }
                if (s > best_any.closeness) {
                    best_any.a = i; best_any.b = j;
                    best_any.closeness = s;
                    // Reported, never used for selection. The spec's rule is
                    // "the two closest" and this is only how the outcome is
                    // classified afterwards.
                    best_any.garbage_pair = bank.templates[i].isJunk()
                        && bank.templates[j].isJunk();
                }
            }
        }
        if (best_any.valid()) return best_any;   // the spec's rule
        return best_blocked;   // only confirmed-touching pairs remain
    }

    // ---------------------------------------------------------------------
    // Display letter per bank slot: A, B, C ... contiguous over the
    // templates that actually exist.
    //
    // THE ONE DEFINITION. The viewer and the morphology CSV/bin writers both
    // name templates, and they used to letter them independently -- the CSV
    // ranked non-empty templates, the viewer used the raw slot index. They
    // agree only when no lower slot is empty and spawn order matches slot
    // order, so the same template could appear as PQRST_C in one place and
    // PQRST_F in the other. Both now call this.
    //
    // ORDERED BY spawn_seq, NOT BY SLOT. A merge erases an element and shifts
    // everything after it, so slot position is not stable between runs while
    // spawn order is. Empty slots are skipped rather than consuming a letter,
    // which is what keeps the sequence contiguous.
    //
    // A CONFIRMED template takes its letter from the subtype the bank issued
    // (PVC_A, PVC_B are per-class), so its letter tracks the class rather than
    // the bank-wide rank. Returns one entry per slot, indexed by slot.
    // ---------------------------------------------------------------------
    inline std::vector<uint8_t> letterRanks(const TemplateBank& bank) {
        const int n = bank.size();
        std::vector<int> order;
        order.reserve(n);
        for (int i = 0; i < n; ++i)
            if (!bank.templates[i].tmpl.empty()) order.push_back(i);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return bank.templates[a].spawn_seq < bank.templates[b].spawn_seq;
            });

        std::vector<uint8_t> letter(n, 0);
        int rank = 0;
        for (int i : order) {
            const BankTemplate& t = bank.templates[i];
            const int idx = (t.label_code != kUnlabeled && t.subtype > 0)
                ? t.subtype - 1 : rank;
            letter[i] = static_cast<uint8_t>(idx % 26);
            ++rank;
        }
        return letter;
    }

}  // namespace tbank