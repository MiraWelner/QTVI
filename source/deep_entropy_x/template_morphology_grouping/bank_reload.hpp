#pragma once
//
// bank_reload.hpp
// if a user marked the templates, the template markings are reloaded. but maybe the config file or algorithm changes - then the bins get recalculated 
// and the markings are wrong. This file reloads the old templates which were split however they were when the user made the markings.
//

#include <array>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "template_generation\template_io.hpp"
#include "template_morphology_grouping/morphology_csv.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>

namespace bank_reload {




    // =======================================================================
    // THE SPLIT, RELOADED FROM THE MORPHOLOGY ARCHIVE
    // =======================================================================
    //
    // THE SPLIT COMES FROM THE MORPHOLOGY ARCHIVE, not from <stem>_bins.bin.
    // A previous reloader (reloadBanks, removed) took whole
    // tbank::TemplateBank objects out of that template_io file, but it holds
    // the per-bin averages and is written and never read, so the banks in it
    // are not the source of truth for the split. This reloader takes it from
    // <stem>_templates.bin instead: morphology_csv's archive, one record per
    // TEMPLATE, which since v4 carries `members` / `members_clean` and since v5
    // carries every remaining BankTemplate field and the bank's own scalars.
    //
    // THE EXISTENCE OF THE FILE IS THE SWITCH. No
    // validation, no config flag. Normally there is no prior archive and the
    // fresh split stands; when there is one, a config change is deliberately a
    // no-op on the partition.
    //
    // WHAT IS RESTORED -- everything the archive carries, which as of v5 is
    // every field of BankTemplate except the landmarks:
    //   members, members_clean      the partition itself
    //   tmpl, tmpl_iqr              the waveform and its spread
    //   r_col, label_code           the column and the class
    //   subtype, spawn_seq          PVC_2 / PAC_3 and the order of first
    //                               appearance that subtype letters resolve
    //                               against
    //   confirmed_by_operator       the operator's verdict
    //   marked_invalid_template,    the right-click quality marks
    //   operator_state
    //   mean_rr_ms                  mean R-R over member slices
    //   split_source                which channel's rejection spawned it
    //   the census counts           n_premature / n_voted / n_noise /
    //                               n_tukey / n_blended / n_ppg, which is what
    //                               presumedCategory() reads -- absent, an
    //                               ectopic template reads REGULAR and is
    //                               handed a landmark column it never had
    //   the bank's caps/counters    configured_cap, effective_cap,
    //                               next_spawn_seq, assigned_beats
    //
    // NOT RESTORED, DELIBERATELY: markers_by_anchor and pulse_marks. Landmarks
    // are operator judgement and live in _template_markings.bin, which nothing
    // but the marking session writes -- see template_bank_serialize.hpp. A
    // reload here must not touch them, and does not.
    //
    // THE SLICE-ORDINAL CAVEAT AT THE TOP OF THIS FILE APPLIES UNCHANGED.
    // members are R-pair slice ordinals, so a config change that altered the
    // slicing makes them name different heartbeats. Accepted for the same
    // reason: to force a repartition, delete the archive.

    struct SplitReport {
        bool prior_read = false;         // the file opened and parsed
        bool prior_present = false;      // the path exists
        bool too_old = false;            // parsed, but pre-v5: refused
        std::string prior_path;
        std::string error;
        size_t templates_restored = 0;   // (bin, channel, slot) triples
        size_t banks_restored = 0;       // (bin, channel) pairs
        size_t banks_skipped = 0;        // slot-count mismatch: left fresh
        size_t beats_restored = 0;       // member indices written back

        bool anythingReloaded() const { return banks_restored != 0; }
    };

    namespace detail {

        // The archive pads every column out to the block's widest template, so
        // a shorter waveform arrives with a NaN tail that is padding, not
        // signal. Trimmed here rather than carried, because tmpl.size() is what
        // every consumer treats as the template's length.
        inline std::vector<double> trimTrailingNaN(const double* p, uint32_t w) {
            if (!p || w == 0) return {};
            uint32_t n = w;
            while (n > 0 && std::isnan(p[n - 1])) --n;
            return std::vector<double>(p, p + n);
        }

        inline int channelIndexOf(const std::string& name) {
            if (name == "CH1") return 0;
            if (name == "CH2") return 1;
            if (name == "CH3") return 2;
            if (name == "PPG") return 3;   // the pulse bank
            return -1;
        }

    }  // namespace detail

    // ---- READ AND APPLY ARE SEPARATE, AND THE ORDER IS THE WHOLE POINT ----
    //
    // <stem>_templates.bin is REWRITTEN BY THIS RUN, inside
    // buildTemplatesAndBeatsFast -> morphology_csv::writeTemplatesBin. So a
    // reload that reads the file after the build reads THIS run's own output:
    // it finds a file every time, reports a successful restore, puts the fresh
    // split back over itself, and a config change is silently undone by the
    // thing it was supposed to survive. That is not a subtle failure mode -- it
    // is indistinguishable from working, because the numbers in the report are
    // real.
    //
    // So the archive is READ BEFORE the build and APPLIED AFTER it. Holding the
    // parsed blocks across the build costs one copy of the per-template records
    // and their waveforms, which is small beside the beat matrices the build
    // already holds.
    struct SplitArchive {
        SplitReport rep;
        std::vector<morphology_csv::BinBlock<morphology_csv::TemplateRecord>> blocks;

        // True when there is something to apply. A first run has no archive and
        // this is false, which is the normal case and not an error.
        bool usable() const { return rep.prior_read && !rep.too_old; }
    };

    // Call BEFORE buildTemplatesAndBeatsFast.
    inline SplitArchive readSplit(const std::string& priorPath) {
        SplitArchive out;
        SplitReport& rep = out.rep;
        rep.prior_path = priorPath;

        std::error_code ec;
        rep.prior_present = std::filesystem::exists(priorPath, ec) && !ec;
        if (!rep.prior_present) return out;

        if (!morphology_csv::readTemplatesBin(priorPath, out.blocks)) {
            rep.error = "readTemplatesBin failed (wrong magic, newer version,"
                " or truncated)";
            return out;
        }
        rep.prior_read = true;

        // Counted here rather than left to the caller: readSplit runs BEFORE
        // the build, so this is the last moment the archive on disk is the
        // PREVIOUS run's. Once buildTemplatesAndBeatsFast has run, the file has
        // been overwritten and there is nothing left to compare against.
        {
            size_t nrec = 0, ntrail = 0;
            for (const auto& blk : out.blocks) {
                nrec += blk.records.size();
                ntrail += blk.trailers.size();
            }
            std::fprintf(stderr, "  [split-reload] read %s: %zu block(s),"
                " %zu template record(s), %zu trailer(s)\n",
                priorPath.c_str(), out.blocks.size(), nrec, ntrail);
        }

        // ---- ALL OF IT OR NONE OF IT ----------------------------------
        //
        // The v5 trailer is what makes a reload exact. Without it there is no
        // census, so presumedCategory() would read an ectopic template as
        // REGULAR, and no subtype, so letters would come back from bank order.
        // Applying the partition anyway produces a bank that disagrees with
        // itself, which is worse than a clean repartition.
        for (const auto& blk : out.blocks) {
            if (blk.records.empty()) continue;
            if (blk.trailers.size() != blk.records.size()) {
                rep.too_old = true;
                rep.error = "archive predates v5 (no per-template trailer), so"
                    " a reload could not be exact";
                out.blocks.clear();
                return out;
            }
        }
        return out;
    }

    // Call AFTER the build, with what readSplit returned.
    inline SplitReport applySplit(SplitArchive& arch,
        template_io::TemplateFile& fresh)
    {
        SplitReport rep = arch.rep;
        if (!arch.usable()) return rep;
        const auto& blocks = arch.blocks;

        for (const auto& blk : blocks) {
            const int c = detail::channelIndexOf(blk.channel);
            if (c < 0) continue;

            // ---- WHICH BINS, AND HOW MANY SLOTS EACH HAD ---------------
            //
            // Records are flat: one per (bin, template), in bin order but with
            // no per-bin header. So the highest template_id per bin is what
            // says how many slots the prior partition had, and it has to be
            // known before anything is written -- a bin is reloaded whole or
            // not at all: a partially reloaded bank is a partition of no population.
            std::vector<int> needSlots(fresh.bins.size(), -1);
            for (const auto& rec : blk.records) {
                if (rec.bin >= fresh.bins.size()) continue;
                if (rec.template_id < 0) continue;
                needSlots[rec.bin] = std::max(needSlots[rec.bin],
                    static_cast<int>(rec.template_id));
            }

            std::vector<char> take(fresh.bins.size(), 0);
            for (size_t b = 0; b < fresh.bins.size(); ++b) {
                if (needSlots[b] < 0) continue;   // no prior templates here
                const tbank::TemplateBank& bank = (c < 3)
                    ? fresh.bins[b].ecg_bank[c] : fresh.bins[b].ppg_bank;
                if (static_cast<int>(bank.templates.size()) > needSlots[b]) {
                    take[b] = 1;
                    ++rep.banks_restored;
                }
                else {
                    ++rep.banks_skipped;
                    std::fprintf(stderr,
                        "  [split-reload] bin %zu %s: prior split names slot %d"
                        " but the fresh bank has %zu -- left freshly split\n",
                        b, blk.channel.c_str(), needSlots[b],
                        bank.templates.size());
                }
            }

            for (size_t k = 0; k < blk.records.size(); ++k) {
                const auto& rec = blk.records[k];
                if (rec.bin >= fresh.bins.size()) continue;
                if (!take[rec.bin]) continue;
                if (rec.template_id < 0) continue;

                tbank::TemplateBank& bank = (c < 3)
                    ? fresh.bins[rec.bin].ecg_bank[c]
                    : fresh.bins[rec.bin].ppg_bank;
                const size_t sl = static_cast<size_t>(rec.template_id);
                if (sl >= bank.templates.size()) continue;   // guarded by take[]
                tbank::BankTemplate& tp = bank.templates[sl];

                tp.members = blk.members[k];
                tp.members_clean = blk.members_clean[k];
                rep.beats_restored += tp.members.size();

                // BOTH TRIMMED TO THE SAME LENGTH, or neither. The waveform
                // was trimmed and the spread was not, so a reloaded slot came
                // back with tmpl_iqr LONGER than tmpl -- and every consumer
                // that pairs them by index and tests the two lengths for
                // equality then decided there was no spread at all. That is
                // why a reloaded pulse slot had no std band.
                tp.tmpl = detail::trimTrailingNaN(blk.column(k), blk.width);
                tp.tmpl_iqr = blk.tmpl_iqr[k];
                if (tp.tmpl_iqr.size() > tp.tmpl.size())
                    tp.tmpl_iqr.resize(tp.tmpl.size());
                tp.r_col = rec.r_col;
                tp.label_code = rec.label_code;

                // confirmed_by_operator comes from the TRAILER, not from
                // TemplateRecord::confirmed. That field is tri-state for the
                // CSV's benefit -- 0 presumed, 1 confirmed, 2 n/a -- so it
                // cannot round-trip a bool.
                const auto& tr = blk.trailers[k];
                tp.subtype = tr.subtype;
                tp.spawn_seq = tr.spawn_seq;
                tp.n_ppg_members = tr.n_ppg_members;
                tp.n_tukey_members = tr.n_tukey_members;
                tp.n_blended_members = tr.n_blended_members;
                tp.n_premature_members = tr.n_premature_members;
                tp.n_voted_members = tr.n_voted_members;
                tp.n_noise_members = tr.n_noise_members;
                tp.mean_rr_ms = tr.mean_rr_ms;
                tp.marked_invalid_template = (tr.marked_invalid_template != 0);
                tp.operator_state = tr.operator_state;
                tp.confirmed_by_operator = (tr.confirmed_by_operator != 0);
                tp.split_source = tr.split_source;

                // Bank scalars ride on every record of the bank -- the format
                // has no per-bank framing to hang them on. The values agree by
                // construction, so the last write equals the first.
                bank.configured_cap = tr.configured_cap;
                bank.effective_cap = tr.effective_cap;
                bank.next_spawn_seq = tr.next_spawn_seq;
                bank.assigned_beats = tr.assigned_beats;

                ++rep.templates_restored;
            }
        }

        return rep;
    }

    inline void printReport(const SplitReport& rep, std::FILE* out = stderr) {
        if (!rep.prior_present) {
            // SAID, not silent. This was silent for one build, on the grounds
            // that a first run has nothing to report -- and that made "no
            // archive on disk" indistinguishable from "archive found and
            // refused", which is the one distinction anybody debugging a
            // failed reload needs. The path is the useful part: it is composed
            // from cfg.template_path and the stem, and a reload that silently
            // does nothing is usually a reload looking in the wrong place.
            std::fprintf(out, "  [split-reload] no archive at %s"
                " -- fresh split stands\n", rep.prior_path.c_str());
            return;
        }
        if (rep.too_old) {
            std::fprintf(out,
                "  [split-reload] %s: %s -- repartitioning. Regenerate once and"
                " the next run reloads exactly.\n",
                rep.prior_path.c_str(), rep.error.c_str());
            return;
        }
        if (!rep.prior_read) {
            // PRESENT AND REJECTED. Shouted, because the record HAS a prior
            // partition and this run is about to replace it with a different
            // one. This is the case that used to print as "ignoring <path>:
            // vector too long" and be taken for a first run.
            std::fprintf(out,
                "  [split-reload] *** PRIOR SPLIT DISCARDED *** %s exists but"
                " would not parse: %s\n"
                "  [split-reload]     the fresh repartition will NOT match it.\n",
                rep.prior_path.c_str(), rep.error.c_str());
            return;
        }
        std::fprintf(out,
            "  [split-reload] %s: %zu template(s) over %zu (bin,channel) bank(s),"
            " %zu beat assignment(s) restored, %zu bank(s) left freshly split\n",
            rep.prior_path.c_str(), rep.templates_restored, rep.banks_restored,
            rep.beats_restored, rep.banks_skipped);
    }

}  // namespace bank_reload