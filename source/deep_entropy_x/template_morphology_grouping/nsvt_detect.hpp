#pragma once
/**
 * @file   nsvt_detect.hpp
 * @brief  Non-sustained ventricular tachycardia runs: three or more
 *         consecutive beats on the same ventricular-labeled template at a rate
 *         above 100 bpm.
 *
 *         THE PROBLEM THE SPEC'S SIGNATURE HIDES. detectNsvt() is declared over
 *         flat per-beat vectors -- templateIdPerBeat, rrMs -- and one
 *         TemplateBank. That is a record-level view. But banks are PER BIN, so
 *         template id 2 in bin 7 and template id 2 in bin 8 are unrelated
 *         morphologies that were never compared to each other. Concatenating
 *         the per-bin ids and scanning for "the same template" therefore
 *         MANUFACTURES RUNS across every bin boundary. And a real run that
 *         straddles a boundary is invisible, because its beats carry two
 *         different local ids.
 *
 *         So NSVT needs something the addendum never mentions: template
 *         identity that is global across bins. globalizeTemplates() below
 *         supplies it by matching each bin's templates against the running set
 *         of known morphologies, using the same correlation and the same floor
 *         the bank uses internally.
 *
 *         WHAT THIS DETECTOR CANNOT DO, BY CONSTRUCTION. The criterion is three
 *         or more consecutive beats on THE SAME template. Polymorphic VT and
 *         torsades change morphology beat to beat, so their beats scatter
 *         across templates and never form a run. That is the clinically more
 *         dangerous presentation, and it sits oddly beside a section whose
 *         stated purpose is preserving the polymorphic/monomorphic
 *         distinction. Implemented as specified; the gap is recorded here and
 *         countPolymorphicCandidates() below reports the runs it would have
 *         caught, so the blind spot is measurable rather than merely known.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "template_bank.hpp"
#include "template_assign.hpp"

namespace nsvt {

    // Section 4.6 addendum thresholds.
    inline constexpr int    kMinRun = 3;       // consecutive beats
    inline constexpr double kRateThresholdBpm = 100.0;   // exclusive
    inline constexpr double kSustainedMs = 30000.0; // 30 s

    // Floor for calling two templates from different bins the same morphology.
    // Deliberately the bank's own assignment floor: a beat and a template that
    // correlate at 0.85 are the same morphology inside a bin, so two templates
    // that correlate at 0.85 are the same morphology across bins. Using a
    // different number here would mean the record-level and bin-level notions
    // of "same shape" disagree.
    inline constexpr double kGlobalMatchFloor = tbank::kMatchFloorEcg;

    // ---------------------------------------------------------------------
    // Cross-bin template identity
    // ---------------------------------------------------------------------

    struct GlobalTemplate {
        std::vector<double> tmpl;        // the exemplar that defined this id
        uint8_t  label_code = tbank::kUnlabeled;
        int32_t  subtype = -1;
        uint32_t first_bin = 0;
        uint32_t n_bins_seen = 0;
        uint32_t n_beats_total = 0;
    };

    struct GlobalMap {
        std::vector<GlobalTemplate> morphologies;

        // key: (bin, channel, local template index) -> global id.
        // Flattened by the caller through globalIdOf().
        std::vector<std::vector<std::vector<int32_t>>> id;   // [bin][chan][local]

        int32_t globalIdOf(size_t bin, int chan, int local) const {
            if (bin >= id.size()) return -1;
            if (chan < 0 || chan >= static_cast<int>(id[bin].size())) return -1;
            if (local < 0 || local >= static_cast<int>(id[bin][chan].size()))
                return -1;
            return id[bin][chan][local];
        }
    };

    // Walks bins in order and assigns each local template either an existing
    // global id (best correlation at or above the floor) or a new one.
    //
    // Label handling here is strictly conservative. A global id inherits a
    // label only from a template that CARRIES one, and never propagates a
    // label backwards onto templates already mapped to that id. Inferring a
    // label from cross-bin morphological similarity would be the same
    // prohibited judgment as inferring one within a bin: whether two similar
    // morphologies are one class is the operator's call. The consequence is
    // that a global id can hold both labeled and unlabeled members, which is
    // correct -- it means the operator confirmed the morphology in one bin and
    // has not yet looked at another.
    inline GlobalMap globalizeTemplates(
        const std::vector<std::array<tbank::TemplateBank, 3>>& per_bin,
        double floor = kGlobalMatchFloor)
    {
        GlobalMap gm;
        gm.id.resize(per_bin.size());

        for (size_t b = 0; b < per_bin.size(); ++b) {
            gm.id[b].resize(3);
            for (int c = 0; c < 3; ++c) {
                const tbank::TemplateBank& bank = per_bin[b][c];
                gm.id[b][c].assign(bank.templates.size(), -1);

                for (int t = 0; t < bank.size(); ++t) {
                    const tbank::BankTemplate& lt = bank.templates[t];
                    if (lt.tmpl.empty()) continue;

                    int    best = -1;
                    double best_r = -std::numeric_limits<double>::infinity();
                    for (size_t g = 0; g < gm.morphologies.size(); ++g) {
                        const tbank::CorrResult cr =
                            tbank::correlate(lt.tmpl, gm.morphologies[g].tmpl);
                        if (!cr.scorable()) continue;
                        if (cr.r > best_r) { best_r = cr.r; best = static_cast<int>(g); }
                    }

                    int32_t gid;
                    if (best >= 0 && best_r >= floor) {
                        gid = best;
                        GlobalTemplate& g = gm.morphologies[gid];
                        ++g.n_bins_seen;
                        g.n_beats_total += lt.memberCount();
                        // Adopt a label only if this local template has one and
                        // the global id does not. Never overwrite: two
                        // different confirmed labels on one global id is a
                        // finding for the operator, not something to silently
                        // resolve.
                        if (g.label_code == tbank::kUnlabeled && lt.confirmed()) {
                            g.label_code = lt.label_code;
                            g.subtype = lt.subtype;
                        }
                    }
                    else {
                        GlobalTemplate g;
                        g.tmpl = lt.tmpl;
                        g.label_code = lt.label_code;
                        g.subtype = lt.subtype;
                        g.first_bin = static_cast<uint32_t>(b);
                        g.n_bins_seen = 1;
                        g.n_beats_total = static_cast<uint32_t>(lt.memberCount());
                        gid = static_cast<int32_t>(gm.morphologies.size());
                        gm.morphologies.push_back(std::move(g));
                    }
                    gm.id[b][c][t] = gid;
                }
            }
        }
        return gm;
    }

    // ---------------------------------------------------------------------
    // Runs
    // ---------------------------------------------------------------------

    struct NsvtRun {
        uint32_t start_beat = 0;
        uint32_t length = 0;        // beats
        int32_t  global_template = -1;
        int32_t  subtype = -1;       // e.g. the 2 in PVC-2
        uint8_t  label_code = tbank::kUnlabeled;

        double   mean_cycle_ms = 0.0;
        double   max_cycle_ms = 0.0;   // see the rate note below
        double   rate_bpm = 0.0;
        double   duration_ms = 0.0;

        bool     sustained = false;
        bool     crosses_bin = false;  // spanned a bin boundary

        uint32_t first_bin = 0;
        uint32_t last_bin = 0;
    };

    struct DetectInput {
        // Per beat, record-level, all indexed identically.
        std::vector<int32_t>  global_template;   // from GlobalMap, -1 unassigned
        std::vector<double>   rr_after;          // R[i+1] - R[i], ms
        std::vector<uint32_t> bin_of_beat;
        int                   min_run = kMinRun;
    };

    // A run's cycle lengths are the intervals BETWEEN its beats, so a run of L
    // beats starting at s uses rr_after[s .. s+L-2] -- L-1 intervals. Using L
    // intervals would fold in the interval leaving the run, which is the
    // compensatory pause and would drag the mean upward, understating the rate
    // of exactly the runs that matter.
    inline void summarizeRun(NsvtRun& r, const std::vector<double>& rr_after)
    {
        double sum = 0.0, mx = 0.0;
        int n = 0;
        for (uint32_t i = r.start_beat; i + 1 < r.start_beat + r.length; ++i) {
            if (i >= rr_after.size() || std::isnan(rr_after[i])) continue;
            sum += rr_after[i];
            mx = std::max(mx, rr_after[i]);
            ++n;
        }
        r.mean_cycle_ms = (n > 0) ? sum / n : 0.0;
        r.max_cycle_ms = mx;
        r.duration_ms = sum;
        r.rate_bpm = (r.mean_cycle_ms > 0.0) ? 60000.0 / r.mean_cycle_ms : 0.0;

        // Sustained is a DURATION, not a beat count. The spec says "a run of 30
        // seconds or longer", and NsvtRun carries a length in beats, so the two
        // are only equivalent at a fixed rate. Duration is the summed cycle
        // lengths.
        r.sustained = r.duration_ms >= kSustainedMs;
    }

    // Scans for runs of consecutive beats sharing one ventricular-labeled
    // global template.
    //
    // THE RATE CRITERION is applied to the MEAN cycle length, per "the run rate
    // exceeds 100 beats per minute". That admits a run containing one slow
    // beat, so max_cycle_ms is retained on every run: a run whose mean clears
    // 100 bpm but whose max cycle is 900 ms is a different object from one
    // where every interval is short, and only the max shows it. If the
    // criterion should instead be every-RR, the change is one line here and the
    // stored max makes previously archived runs re-filterable without a rerun.
    //
    // Runs are allowed to CROSS BIN BOUNDARIES, which is the entire reason
    // global ids exist. crosses_bin is recorded because such a run is assembled
    // from two banks' templates that globalizeTemplates() judged identical, and
    // that judgment is worth being able to audit.
    inline std::vector<NsvtRun> detectRuns(const DetectInput& in,
        const GlobalMap& gm)
    {
        std::vector<NsvtRun> runs;
        const size_t n = in.global_template.size();
        if (n == 0) return runs;

        size_t i = 0;
        while (i < n) {
            const int32_t gid = in.global_template[i];
            if (gid < 0) { ++i; continue; }

            // Only ventricular-labeled morphologies qualify. An UNLABELED
            // template is not eligible however ventricular it looks -- that
            // label is the operator's judgment and NSVT must not manufacture
            // it. So a genuine VT run in an unconfirmed template produces no
            // run until someone confirms a beat in it, which is deliberate:
            // the alternative is an algorithm reporting ventricular
            // tachycardia on its own authority.
            const bool eligible =
                gid < static_cast<int32_t>(gm.morphologies.size())
                && tbank::isVentricular(gm.morphologies[gid].label_code);
            if (!eligible) { ++i; continue; }

            size_t j = i;
            while (j + 1 < n && in.global_template[j + 1] == gid) ++j;
            const uint32_t len = static_cast<uint32_t>(j - i + 1);

            if (static_cast<int>(len) >= in.min_run) {
                NsvtRun r;
                r.start_beat = static_cast<uint32_t>(i);
                r.length = len;
                r.global_template = gid;
                r.label_code = gm.morphologies[gid].label_code;
                r.subtype = gm.morphologies[gid].subtype;
                summarizeRun(r, in.rr_after);

                if (r.rate_bpm > kRateThresholdBpm) {
                    if (i < in.bin_of_beat.size()) r.first_bin = in.bin_of_beat[i];
                    if (j < in.bin_of_beat.size()) r.last_bin = in.bin_of_beat[j];
                    r.crosses_bin = (r.first_bin != r.last_bin);
                    runs.push_back(r);
                }
            }
            i = j + 1;
        }
        return runs;
    }

    // The spec says a run of 30 s or more "is sustained VT and is escalated
    // rather than logged as NSVT", but NsvtRun carries a `sustained` flag,
    // which implies it IS logged, with the flag set. Both are honoured: every
    // run stays in the one list with its flag, and these accessors split it, so
    // an escalation path can consume one view while the NSVT archive consumes
    // the other. Nothing is dropped, because a run silently absent from both
    // lists would be the worst outcome available.
    inline std::vector<NsvtRun> nonSustained(const std::vector<NsvtRun>& runs) {
        std::vector<NsvtRun> out;
        for (const auto& r : runs) if (!r.sustained) out.push_back(r);
        return out;
    }
    inline std::vector<NsvtRun> escalate(const std::vector<NsvtRun>& runs) {
        std::vector<NsvtRun> out;
        for (const auto& r : runs) if (r.sustained) out.push_back(r);
        return out;
    }

    // Measures the blind spot rather than leaving it as a comment. Counts
    // stretches of kMinRun or more consecutive beats that are all on
    // ventricular-labeled templates but NOT all the same one -- i.e. exactly
    // what polymorphic VT looks like and what detectRuns() cannot report. A
    // non-zero count here on a record with no detected runs is the signal that
    // the same-template criterion is costing something real.
    inline uint32_t countPolymorphicCandidates(const DetectInput& in,
        const GlobalMap& gm)
    {
        const size_t n = in.global_template.size();
        uint32_t found = 0;
        size_t i = 0;
        auto ventricular = [&](size_t k) {
            const int32_t g = in.global_template[k];
            return g >= 0 && g < static_cast<int32_t>(gm.morphologies.size())
                && tbank::isVentricular(gm.morphologies[g].label_code);
            };
        while (i < n) {
            if (!ventricular(i)) { ++i; continue; }
            size_t j = i;
            bool mixed = false;
            while (j + 1 < n && ventricular(j + 1)) {
                if (in.global_template[j + 1] != in.global_template[i]) mixed = true;
                ++j;
            }
            if (mixed && static_cast<int>(j - i + 1) >= in.min_run) ++found;
            i = j + 1;
        }
        return found;
    }

}  // namespace nsvt
