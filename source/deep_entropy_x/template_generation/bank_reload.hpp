#pragma once
//
// bank_reload.hpp
//
// A PRIOR templates.bin IS THE RECORD'S MORPHOLOGY SPLIT.
//
// GenerateTemplatesFast repartitions every bin from the beats on every run
// (~1 s per bin), and nothing ever read the previous answer back -- so editing
// a morphology threshold in config.csv silently discarded a split that may
// already have been reviewed and labelled by an operator, with no way to get
// it back except restoring the old config and hoping the rest of the pipeline
// was byte-identical.
//
// The split is already on disk. writeTemplateInfoBin persists it as four
// trailing sections: v3 (ECG banks, per bin per channel), v4 (PPG bank, per
// bin), v5 (per-template extras -- confirmed_by_operator, members_clean, the
// census), v6 (per-anchor bank slot averages). read_template_binfile reads all
// four back. What was missing was a caller that preferred them.
//
// THE RULE: if the file exists, its banks replace the freshly built ones.
// Unconditionally. No validation, no per-bin fallback, no config flag. THE
// EXISTENCE OF THE FILE IS THE SWITCH -- normally there is no prior file and
// the fresh split stands; when there is one, a config change is deliberately a
// no-op on the partition. Anything that made the reuse conditional would
// defeat the purpose.
//
// ---------------------------------------------------------------------------
// WHAT THIS MEANS IF THE BEAT SET MOVED
// ---------------------------------------------------------------------------
//
// A slot's membership (BankTemplate::members / members_clean) is a list of
// R-pair SLICE ORDINALS, meaningful only against the slicing that produced
// them. If a config change altered the slicing itself -- bin_size_minutes,
// ecg_upsample_rate, the drop rules, or anything upstream that re-runs R-peak
// detection -- the reloaded memberships name slices that either do not exist
// or are different heartbeats, and consumers that resolve a member to a beat
// resolve it wrongly.
//
// Accepted by design. A guard that rejected those bins would hand the
// partition back to config.csv for exactly the bins where an operator's review
// is most expensive to lose. TO FORCE A REPARTITION, DELETE templates.bin --
// that is the documented way, and the only one.
//
// ---------------------------------------------------------------------------
// WHAT IS AND IS NOT RELOADED
// ---------------------------------------------------------------------------
//
// RELOADED: the bank, and everything riding inside it -- tmpl, tmpl_iqr,
// r_col, label_code, subtype, confirmed_by_operator, operator_state, members,
// members_clean, spawn_seq, n_ppg_members, the per-template counts. All are
// fields of tbank::TemplateBank, so copying the bank carries them. That is why
// the bank is taken WHOLE rather than merged: an operator confirmation belongs
// to the split it was made against, and half of one split plus half of another
// describes no beat set that ever existed.
//
// NOT RELOADED, because templates.bin does not carry it: the per-slice
// bookkeeping in TemplateInfo::joint -- group_of_slice, flags, pvc,
// excluded_reason, rr_after_ms, and the counts/clean/subs census. Those are
// rebuilt fresh every run. So after a reload the templates on screen are the
// prior split while _bins.csv and the NSVT rows describe the newly computed
// one, and the two will disagree. Closing that needs a v7 section carrying the
// per-bin JointBinResult -- which is also the prerequisite for skipping
// buildBinBank altogether and recovering the ~1 s per bin.
//

#include <array>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "template_io.hpp"

namespace bank_reload {

    struct Report {
        bool   prior_read = false;      // the file opened and parsed
        std::string prior_path;
        std::string error;              // set when prior_read is false
        size_t prior_bins = 0;
        size_t fresh_bins = 0;
        size_t ecg_reloaded = 0;        // (bin, channel) banks replaced
        size_t ppg_reloaded = 0;
        size_t anchor_slot_sets = 0;    // (anchor, bin, channel) slot vectors carried over

        bool anythingReloaded() const { return ecg_reloaded || ppg_reloaded; }
    };

    // ------------------------------------------------------------------------
    // Replace `fresh`'s banks with `priorPath`'s, for every bin the prior file
    // covers. Call after the build and before write_template_binfile.
    //
    // A missing or unreadable prior file is not an error -- it is the first
    // run, which is the normal case. `fresh` is left untouched and
    // Report::prior_read stays false.
    //
    // BINS BEYOND THE PRIOR FILE'S COUNT keep their freshly built banks, since
    // there is nothing on disk to put there. That is the only case where a
    // fresh bank survives a reload.
    // ------------------------------------------------------------------------
    inline Report reloadBanks(const std::string& priorPath,
        template_io::TemplateFile& fresh) {
        Report rep;
        rep.prior_path = priorPath;
        rep.fresh_bins = fresh.bins.size();

        template_io::TemplateFile prior;
        try {
            prior = template_io::read_template_binfile(priorPath);
        }
        catch (const std::exception& e) {
            rep.error = e.what();
            return rep;
        }
        rep.prior_read = true;
        rep.prior_bins = prior.bins.size();

        const size_t n = (fresh.bins.size() < prior.bins.size())
            ? fresh.bins.size() : prior.bins.size();

        // ---- v3 / v4 / v5: the banks themselves ---------------------------
        //
        // An EMPTY prior bank is skipped rather than copied. That is not a
        // validation check -- it is the difference between "the prior file says
        // this bin has one group" and "the prior file has nothing for this bin"
        // (a pre-v3 file, or a bin its build skipped). Copying an empty bank
        // over a real one would delete the split rather than reload it.
        for (size_t i = 0; i < n; ++i) {
            for (int c = 0; c < 3; ++c) {
                if (prior.bins[i].ecg_bank[c].templates.empty()) continue;
                fresh.bins[i].ecg_bank[c] = prior.bins[i].ecg_bank[c];
                ++rep.ecg_reloaded;
            }
            if (!prior.bins[i].ppg_bank.templates.empty()) {
                fresh.bins[i].ppg_bank = prior.bins[i].ppg_bank;
                ++rep.ppg_reloaded;
            }
        }

        // ---- v6: per-anchor slot averages ---------------------------------
        //
        // Averages OF the slots, so they belong with the bank they came from
        // and are copied for the same (bin, channel) pairs. A fresh anchor
        // entry is created sized to fresh.bins.size() when the prior file
        // carries an anchor this build does not, so the indexing contract
        // (parallel to `bins`) holds either way.
        for (const auto& kv : prior.bank_anchors) {
            const int tag = kv.first;
            const auto& priorPerBin = kv.second;

            auto& freshPerBin = fresh.bank_anchors[tag];
            if (freshPerBin.size() < fresh.bins.size())
                freshPerBin.resize(fresh.bins.size());

            for (size_t i = 0; i < n && i < priorPerBin.size(); ++i) {
                for (int c = 0; c < 3; ++c) {
                    if (priorPerBin[i][c].empty()) continue;
                    freshPerBin[i][c] = priorPerBin[i][c];
                    ++rep.anchor_slot_sets;
                }
            }
        }

        return rep;
    }

    inline void printReport(const Report& rep, std::FILE* out = stderr) {
        if (!rep.prior_read) {
            std::fprintf(out, "[bank-reload] no prior split at %s%s%s\n",
                rep.prior_path.c_str(),
                rep.error.empty() ? " (absent -- fresh split stands)" : " -- ",
                rep.error.c_str());
            return;
        }
        std::fprintf(out,
            "[bank-reload] %s: %zu prior bins vs %zu fresh | ECG %zu bank(s) reloaded"
            " | PPG %zu | %zu anchor slot sets\n",
            rep.prior_path.c_str(), rep.prior_bins, rep.fresh_bins,
            rep.ecg_reloaded, rep.ppg_reloaded, rep.anchor_slot_sets);

        // A bin-count difference means the SLICING changed, not just a
        // threshold -- so the reloaded memberships refer to a beat set that is
        // no longer there. Said plainly, once, because the reload proceeds
        // anyway and this is the only warning of it.
        if (rep.fresh_bins > rep.prior_bins)
            std::fprintf(out, "[bank-reload] %zu bin(s) beyond the prior file keep"
                " their freshly built split\n",
                rep.fresh_bins - rep.prior_bins);
        if (rep.prior_bins > rep.fresh_bins)
            std::fprintf(out, "[bank-reload] WARNING: prior file has %zu MORE bin(s)"
                " than this build -- the slicing changed, so reloaded slot"
                " memberships name slices that no longer exist. Delete"
                " templates.bin to repartition.\n",
                rep.prior_bins - rep.fresh_bins);
    }

}  // namespace bank_reload
