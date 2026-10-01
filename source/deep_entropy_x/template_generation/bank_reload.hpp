#pragma once
//
// bank_reload.hpp
// if a user marked the templates, the template markings are reloaded. but maybe the config file or algorithm changes - then the bins get recalculated 
// and the markings are wrong. This file reloads the old templates which were split however they were when the user made the markings.
//
// ---- WHEN A PRIOR PARTITION IS RELOADED, AND WHEN IT IS NOT -----------------
//
// PER BIN, AND ONLY WHEN BOTH ARE TRUE:
//
//   1. THE OPERATOR RULED ON IT. At least one of the bin's prior templates, on
//      any channel, carries an operator verdict: confirmed, crossed out, or
//      marked bad. A bin nobody reviewed has nothing worth preserving, so it
//      takes the fresh partition -- which is what a config change is for.
//
//   2. THE SLICING IS UNCHANGED. Members are LOCAL ROWS of each channel's kept
//      beats, so they only name the same heartbeats if the bin was sliced from
//      the same R-peaks and kept the same beats. Checked against a fingerprint
//      of exactly that, written beside the archive (<stem>_slicing.bin) by the
//      run that produced it. No fingerprint, or a different one, is a changed
//      slicing: a prior partition applied to different rows is a partition of
//      a different population, and it looks identical to a correct one.
//
// IN EVERY OTHER CASE THE BIN IS REPARTITIONED -- including an archive from a
// build that predates the fingerprint, which therefore repartitions once and
// reloads exactly from then on.
//
// WHOLE BINS. The four channels of a bin are faces of ONE joint partition
// (template i of each channel is group i), so a bin is reloaded on every
// channel or on none: a bin whose ECG came back from the archive and whose
// pulse was freshly split would pair groups that are not the same heartbeats.
//
// APPLIED INSIDE THE BUILD, BEFORE <stem>_templates.bin IS WRITTEN. The build
// rewrites that archive, and the commit after the session patches only the
// operator's verdict into it. Applied after the write, the archive on disk
// would hold the FRESH members under the reloaded verdicts, and the next run
// would reload a partition the operator never saw. See GenerateTemplatesFast.
//

#include <array>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "template_generation/template_structs.hpp"
#include "template_generation/template_io.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>

namespace bank_reload {

    // =======================================================================
    // THE SLICING FINGERPRINT: <stem>_slicing.bin
    // =======================================================================
    //
    // One 64-bit hash per bin over everything that decides which heartbeat a
    // member index names: the bin's R-peaks (ch1.raw, which drives every
    // channel's slicer), the ECG and pulse rates, and each channel's kept-row
    // -> R-pair map. Equal hashes mean local row k of every channel is the same
    // heartbeat, sliced from the same samples, in both runs.
    //
    // 0 IS "NOT FINGERPRINTED" (a bad segment, or a bin with no R-peaks), and
    // never matches -- so such a bin always repartitions.
    //
    // Layout: 8-byte magic, uint32 version, uint64 bin count, then one uint64
    // per bin. Small, and read whole.
    namespace slicing {

        inline constexpr char kMagic[8] = { 'S','L','I','C','E','F','P','1' };
        inline constexpr uint32_t kVersion = 1;

        struct Hasher {
            uint64_t h = 1469598103934665603ull;   // FNV-1a 64 offset basis
            void bytes(const void* p, std::size_t n) {
                const unsigned char* c = static_cast<const unsigned char*>(p);
                for (std::size_t i = 0; i < n; ++i) { h ^= c[i]; h *= 1099511628211ull; }
            }
            void u64(uint64_t v) { bytes(&v, sizeof v); }
            void f64(double v) { bytes(&v, sizeof v); }
        };

        // ecgKept[c] / ppgKept: kept row -> R-pair ordinal for that channel,
        // or null when the channel has none in this bin. Lengths are hashed
        // with the contents, so an absent channel and an empty one agree and
        // neither collides with a present one.
        inline uint64_t hashBin(const std::vector<size_t>& rPeaks,
            double ecgRate, double ppgRate,
            const std::array<const std::vector<size_t>*, 3>& ecgKept,
            const std::vector<uint32_t>* ppgKept)
        {
            if (rPeaks.size() < 2) return 0;
            Hasher hs;
            hs.f64(ecgRate);
            hs.f64(ppgRate);
            hs.u64(rPeaks.size());
            for (const size_t r : rPeaks) hs.u64(r);
            for (int c = 0; c < 3; ++c) {
                hs.u64(0xEC60ull + static_cast<uint64_t>(c));   // channel tag
                const std::vector<size_t>* k = ecgKept[c];
                hs.u64(k ? k->size() : 0);
                if (k) for (const size_t v : *k) hs.u64(v);
            }
            hs.u64(0x9960ull);                                   // pulse tag
            hs.u64(ppgKept ? ppgKept->size() : 0);
            if (ppgKept) for (const uint32_t v : *ppgKept) hs.u64(v);
            return (hs.h == 0) ? 1 : hs.h;   // 0 is reserved for "none"
        }

        // <stem>_templates.bin -> <stem>_slicing.bin, so the reader derives
        // the fingerprint's path from the archive's and cannot look elsewhere.
        inline std::string pathFor(const std::string& templatesPath) {
            const std::string suf = "_templates.bin";
            if (templatesPath.size() >= suf.size()
                && templatesPath.compare(templatesPath.size() - suf.size(),
                    suf.size(), suf) == 0)
                return templatesPath.substr(0, templatesPath.size() - suf.size())
                + "_slicing.bin";
            return templatesPath + ".slicing";
        }

        inline bool write(const std::string& path, const std::vector<uint64_t>& perBin) {
            std::ofstream f(path, std::ios::binary | std::ios::trunc);
            if (!f) return false;
            const uint64_t n = perBin.size();
            f.write(kMagic, sizeof kMagic);
            f.write(reinterpret_cast<const char*>(&kVersion), sizeof kVersion);
            f.write(reinterpret_cast<const char*>(&n), sizeof n);
            if (n) f.write(reinterpret_cast<const char*>(perBin.data()),
                static_cast<std::streamsize>(n * sizeof(uint64_t)));
            return f.good();
        }

        inline bool read(const std::string& path, std::vector<uint64_t>& out) {
            out.clear();
            std::ifstream f(path, std::ios::binary);
            if (!f) return false;
            char magic[8] = {};
            uint32_t ver = 0;
            uint64_t n = 0;
            if (!f.read(magic, sizeof magic)
                || !std::equal(magic, magic + 8, kMagic)) return false;
            if (!f.read(reinterpret_cast<char*>(&ver), sizeof ver)
                || ver != kVersion) return false;
            if (!f.read(reinterpret_cast<char*>(&n), sizeof n)
                || n > (1ull << 32)) return false;
            out.resize(static_cast<size_t>(n));
            if (n && !f.read(reinterpret_cast<char*>(out.data()),
                static_cast<std::streamsize>(n * sizeof(uint64_t)))) {
                out.clear();
                return false;
            }
            return true;
        }

    }  // namespace slicing


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

        // ---- THE TWO GATES (see the top of this file), counted per bin ----
        bool fp_present = false;         // <stem>_slicing.bin exists
        bool fp_read = false;            // ...and parsed
        std::string fp_path;
        size_t bins_with_prior = 0;      // bins the archive has templates for
        size_t bins_reloaded = 0;
        size_t bins_no_verdict = 0;      // nobody ruled on it: repartitioned
        size_t bins_slicing_changed = 0; // fingerprint differs or is absent
        size_t bins_slot_mismatch = 0;   // fresh bank too small on a channel
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
        std::vector<templates_io::BinBlock<templates_io::TemplateRecord>> blocks;
        // The PREVIOUS run's slicing fingerprint, read with the archive and for
        // the same reason: the build rewrites both.
        std::vector<uint64_t> prior_fp;

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

        // The fingerprint FIRST, before any early return below: it is read now
        // or never, because the build overwrites it.
        rep.fp_path = slicing::pathFor(priorPath);
        rep.fp_present = std::filesystem::exists(rep.fp_path, ec) && !ec;
        if (rep.fp_present)
            rep.fp_read = slicing::read(rep.fp_path, out.prior_fp);

        if (!templates_io::readTemplatesBin(priorPath, out.blocks)) {
            // NO MAGIC TO BE WRONG. This file has never had one -- see the
            // "NO MAGIC" note in morphology_csv.hpp -- so the only header
            // check is the version, and the rest is a short read.
            rep.error = "readTemplatesBin failed (stale version or truncated)";
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
        // THE PER-TEMPLATE TRAILER IS WHAT MAKES A RELOAD EXACT. Without it
        // there is no census, so presumedCategory() would read an ectopic
        // template as REGULAR, and no subtype, so letters would come back from
        // bank order. Applying the partition anyway produces a bank that
        // disagrees with itself, which is worse than a clean repartition.
        //
        // Not a version test: one format, always a trailer per record, so a
        // mismatch means truncated or foreign.
        for (const auto& blk : out.blocks) {
            if (blk.records.empty()) continue;
            if (blk.trailers.size() != blk.records.size()) {
                rep.too_old = true;
                rep.error = "one trailer per template record is missing, so the"
                    " archive is truncated or foreign and a reload could not be"
                    " exact";
                out.blocks.clear();
                return out;
            }
        }
        return out;
    }

    // A verdict the operator made on one prior template: confirmed, crossed
    // out, or marked bad on either channel.
    inline bool hasOperatorVerdict(const templates_io::detail::TemplateTrailer& tr) {
        return tr.confirmed_by_operator != 0
            || tr.marked_invalid_template != 0
            || tr.operator_state != tbank::template_of_all_signals::kOperatorGood;
    }

    // Call INSIDE THE BUILD, after the fresh partition exists and BEFORE
    // <stem>_templates.bin is written -- see the top of this file.
    //
    // bankOf(bin, c) is the fresh bank for channel c (0-2 ECG, 3 pulse), or
    // null where the build has none. freshFp is this run's fingerprint per bin
    // (slicing::hashBin), in the same bin order as the archive.
    inline SplitReport applySplit(SplitArchive& arch, size_t nBins,
        const std::function<tbank::TemplateBank* (size_t, int)>& bankOf,
        const std::vector<uint64_t>& freshFp)
    {
        SplitReport rep = arch.rep;
        // Counts are this application's, not carried from an earlier one.
        rep.templates_restored = rep.banks_restored = rep.banks_skipped = 0;
        rep.beats_restored = 0;
        rep.bins_with_prior = rep.bins_reloaded = rep.bins_no_verdict = 0;
        rep.bins_slicing_changed = rep.bins_slot_mismatch = 0;
        if (!arch.usable()) { arch.rep = rep; return rep; }
        const auto& blocks = arch.blocks;

        // ---- WHICH BINS THE ARCHIVE SPEAKS FOR, AND WHAT IT SAYS ----------
        //
        // Records are flat: one per (bin, template), in bin order but with no
        // per-bin header. So the highest template_id per (bin, channel) is what
        // says how many slots the prior partition had.
        std::vector<char> hasPrior(nBins, 0), verdict(nBins, 0);
        std::vector<std::array<int, 4>> needSlots(nBins, { -1, -1, -1, -1 });
        for (const auto& blk : blocks) {
            const int c = detail::channelIndexOf(blk.channel);
            if (c < 0) continue;
            for (size_t k = 0; k < blk.records.size(); ++k) {
                const auto& rec = blk.records[k];
                if (rec.bin >= nBins || rec.template_id < 0) continue;
                hasPrior[rec.bin] = 1;
                needSlots[rec.bin][c] = std::max(needSlots[rec.bin][c],
                    static_cast<int>(rec.template_id));
                if (k < blk.trailers.size() && hasOperatorVerdict(blk.trailers[k]))
                    verdict[rec.bin] = 1;
            }
        }

        // ---- THE GATES, PER BIN, ALL CHANNELS OR NONE ---------------------
        std::vector<char> take(nBins, 0);
        for (size_t b = 0; b < nBins; ++b) {
            if (!hasPrior[b]) continue;
            ++rep.bins_with_prior;
            if (!verdict[b]) { ++rep.bins_no_verdict; continue; }

            const bool sameSlicing = rep.fp_read
                && b < arch.prior_fp.size() && b < freshFp.size()
                && arch.prior_fp[b] != 0 && arch.prior_fp[b] == freshFp[b];
            if (!sameSlicing) { ++rep.bins_slicing_changed; continue; }

            bool slotsOk = true;
            for (int c = 0; c < 4 && slotsOk; ++c) {
                if (needSlots[b][c] < 0) continue;
                const tbank::TemplateBank* bank = bankOf(b, c);
                const int have = bank ? static_cast<int>(bank->templates.size()) : 0;
                if (have <= needSlots[b][c]) {
                    slotsOk = false;
                    std::fprintf(stderr,
                        "  [split-reload] bin %zu channel %d: prior split names slot %d"
                        " but the fresh bank has %d -- bin left freshly split\n",
                        b, c, needSlots[b][c], have);
                }
            }
            if (!slotsOk) { ++rep.bins_slot_mismatch; continue; }

            take[b] = 1;
            ++rep.bins_reloaded;
            for (int c = 0; c < 4; ++c) if (needSlots[b][c] >= 0) ++rep.banks_restored;
        }

        for (const auto& blk : blocks) {
            const int c = detail::channelIndexOf(blk.channel);
            if (c < 0) continue;

            for (size_t k = 0; k < blk.records.size(); ++k) {
                const auto& rec = blk.records[k];
                if (rec.bin >= nBins) continue;
                if (!take[rec.bin]) continue;
                if (rec.template_id < 0) continue;

                tbank::TemplateBank* bankP = bankOf(rec.bin, c);
                if (!bankP) continue;                        // guarded by take[]
                tbank::TemplateBank& bank = *bankP;
                const size_t sl = static_cast<size_t>(rec.template_id);
                if (sl >= bank.templates.size()) continue;   // guarded by take[]
                tbank::template_of_all_signals& tp = bank.templates[sl];

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
                tp.tmpl_std = blk.tmpl_iqr[k];
                if (tp.tmpl_std.size() > tp.tmpl.size())
                    tp.tmpl_std.resize(tp.tmpl.size());
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

        arch.rep = rep;
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
        if (!rep.fp_read) {
            // THE ARCHIVE IS FINE AND THE SLICING CANNOT BE VOUCHED FOR, which
            // under the rule at the top of this file repartitions everything.
            // Said loudly: if the prior run had operator verdicts, they are
            // about to be dropped, and that has to be visible.
            std::fprintf(out,
                "  [split-reload] %s: no usable slicing fingerprint at %s (%s)"
                " -- every bin repartitioned. %zu bin(s) had operator verdicts"
                " that do NOT carry over. The next run reloads exactly.\n",
                rep.prior_path.c_str(), rep.fp_path.c_str(),
                rep.fp_present ? "unreadable" : "absent",
                rep.bins_with_prior - rep.bins_no_verdict);
            return;
        }
        std::fprintf(out,
            "  [split-reload] %s: %zu of %zu bin(s) reloaded (%zu template(s),"
            " %zu beat assignment(s)); repartitioned: %zu with no operator"
            " verdict, %zu with changed slicing, %zu with too few fresh slots\n",
            rep.prior_path.c_str(), rep.bins_reloaded, rep.bins_with_prior,
            rep.templates_restored, rep.beats_restored, rep.bins_no_verdict,
            rep.bins_slicing_changed, rep.bins_slot_mismatch);
    }

}  // namespace bank_reload