#pragma once
/**
 * @file   template_bank_serialize.hpp
 * @brief  Per-bin bank state on disk, appended to the existing markings record
 *         rather than woven into it.
 *
 *         APPEND-ONLY, AND WHY THAT MATTERS HERE. The v1 per-bin record is the
 *         legacy record, byte for byte, followed by one bank block. Nothing in
 *         the legacy layout moves. That buys three things:
 *
 *          1. writeTemplateMarkingsBin/readTemplateMarkingsBin keep their
 *             existing field order untouched -- the riskiest edit in a format
 *             with no version field is reordering what is already there.
 *          2. A legacy file loads as a one-template-per-channel bank with no
 *             migration pass: absent bank block means slot 0 only, which is
 *             exactly what those files describe.
 *          3. If the bank block ever needs to grow again, it grows at the end
 *             of itself, and the block's own length prefix lets an older v1
 *             reader skip what it does not understand instead of walking off
 *             alignment.
 *
 *         That last point is the reason for block_bytes. A format that cannot
 *         skip an unknown tail can only ever be read by exactly the build that
 *         wrote it, which is how marker positions end up silently in the wrong
 *         bins -- the failure mode the versioned header in template_bank_io.hpp
 *         exists to prevent.
 *
 *         WHAT IS NOT STORED. Nothing derived. Beat shares are recomputed from
 *         member counts, the polymorphic verdict is recomputed from label
 *         codes, and template waveforms are recomputed from member indices plus
 *         the beats. Storing a median alongside the members that produce it
 *         invites the two to disagree, and the members are the smaller of the
 *         two. Templates ARE written, because the viewer needs them before the
 *         beats are loaded, but they are written as a cache that a reader may
 *         discard and rebuild -- see kTemplatesAreCache.
 */

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <vector>

#include "template_bank.hpp"
#include "template_bank_io.hpp"
#include "seed_pool.hpp"

namespace tbank_ser {

    // Templates on disk are a rebuildable cache, not the source of truth. A
    // reader that has the beats should prefer recomputeTemplate() over what it
    // read, and any disagreement is a bug in whoever wrote the file.
    inline constexpr bool kTemplatesAreCache = true;

    // ---------------------------------------------------------------------
    // Primitives. Explicit widths throughout: the legacy writer's mix of
    // 1/4/8-byte fields with no width markers is exactly what makes a stride
    // change unrecoverable, and there is no reason to repeat it.
    // ---------------------------------------------------------------------

    namespace detail {

        inline void w8(std::ofstream& f, uint8_t v) {
            f.write(reinterpret_cast<const char*>(&v), 1);
        }
        inline void w32(std::ofstream& f, int32_t v) {
            f.write(reinterpret_cast<const char*>(&v), 4);
        }
        inline void wu32(std::ofstream& f, uint32_t v) {
            f.write(reinterpret_cast<const char*>(&v), 4);
        }
        inline void wu64(std::ofstream& f, uint64_t v) {
            f.write(reinterpret_cast<const char*>(&v), 8);
        }
        inline void wf64(std::ofstream& f, double v) {
            f.write(reinterpret_cast<const char*>(&v), 8);
        }
        inline void wvecd(std::ofstream& f, const std::vector<double>& v) {
            wu32(f, static_cast<uint32_t>(v.size()));
            if (!v.empty())
                f.write(reinterpret_cast<const char*>(v.data()),
                    static_cast<std::streamsize>(v.size() * sizeof(double)));
        }
        inline void wvecu32(std::ofstream& f, const std::vector<uint32_t>& v) {
            wu32(f, static_cast<uint32_t>(v.size()));
            if (!v.empty())
                f.write(reinterpret_cast<const char*>(v.data()),
                    static_cast<std::streamsize>(v.size() * sizeof(uint32_t)));
        }

        inline uint8_t r8(std::ifstream& f) {
            uint8_t v = 0; f.read(reinterpret_cast<char*>(&v), 1); return v;
        }
        inline int32_t r32(std::ifstream& f) {
            int32_t v = 0; f.read(reinterpret_cast<char*>(&v), 4); return v;
        }
        inline uint32_t ru32(std::ifstream& f) {
            uint32_t v = 0; f.read(reinterpret_cast<char*>(&v), 4); return v;
        }
        inline uint64_t ru64(std::ifstream& f) {
            uint64_t v = 0; f.read(reinterpret_cast<char*>(&v), 8); return v;
        }
        inline double rf64(std::ifstream& f) {
            double v = 0.0; f.read(reinterpret_cast<char*>(&v), 8); return v;
        }

        // Length-prefixed reads are bounded so a corrupt or truncated file
        // cannot make the reader allocate gigabytes before it fails. The cap is
        // generous relative to any real bin -- a 15-minute bin holds low
        // thousands of beats and a template is low thousands of samples -- and
        // exists only to turn a wild length into an exception.
        inline constexpr uint32_t kMaxVectorLen = 4u * 1024u * 1024u;

        inline uint32_t readLen(std::ifstream& f) {
            const uint32_t nlen = ru32(f);
            if (nlen > kMaxVectorLen)
                throw std::runtime_error(
                    "template bank block: implausible vector length");
            return nlen;
        }
        inline std::vector<double> rvecd(std::ifstream& f) {
            const uint32_t nlen = readLen(f);
            std::vector<double> v(nlen);
            if (nlen)
                f.read(reinterpret_cast<char*>(v.data()),
                    static_cast<std::streamsize>(nlen * sizeof(double)));
            return v;
        }
        inline std::vector<uint32_t> rvecu32(std::ifstream& f) {
            const uint32_t nlen = readLen(f);
            std::vector<uint32_t> v(nlen);
            if (nlen)
                f.read(reinterpret_cast<char*>(v.data()),
                    static_cast<std::streamsize>(nlen * sizeof(uint32_t)));
            return v;
        }

    }  // namespace detail

    // ---------------------------------------------------------------------
    // The block: everything Section 4.6 adds to a bin
    // ---------------------------------------------------------------------

    struct BinBankState {
        std::array<tbank::TemplateBank, 3> ecg;
        tbank::TemplateBank                ppg;
        std::array<seed_pool::SeedSelection, 3> seed;
        std::vector<tbank::CapRaiseEvent>  cap_raises;
        tbank::BinCounts                   counts;
    };

    // ---------------------------------------------------------------------
    // Write
    // ---------------------------------------------------------------------

    namespace detail {

        inline void writeMarkerSet(std::ofstream& f, const tbank::BankMarkerSet& m) {
            w32(f, m.p_begin); w32(f, m.p_peak); w32(f, m.q_begin);
            w32(f, m.s_end);   w32(f, m.t_begin); w32(f, m.t_end);
        }
        inline tbank::BankMarkerSet readMarkerSet(std::ifstream& f) {
            tbank::BankMarkerSet m;
            m.p_begin = r32(f); m.p_peak = r32(f); m.q_begin = r32(f);
            m.s_end = r32(f); m.t_begin = r32(f); m.t_end = r32(f);
            return m;
        }

        inline void writeTemplate(std::ofstream& f, const tbank::BankTemplate& t) {
            wvecd(f, t.tmpl);
            wvecd(f, t.tmpl_iqr);
            w32(f, t.r_col);
            wvecu32(f, t.members);
            w8(f, t.label_code);
            w32(f, t.subtype);
            wu32(f, t.spawn_seq);

            wu32(f, static_cast<uint32_t>(t.markers_by_anchor.size()));
            for (const auto& kv : t.markers_by_anchor) {
                w32(f, kv.first);                        // AnchorType tag
                writeMarkerSet(f, kv.second);
            }
        }

        inline tbank::BankTemplate readTemplate(std::ifstream& f) {
            tbank::BankTemplate t;
            t.tmpl = rvecd(f);
            t.tmpl_iqr = rvecd(f);
            t.r_col = r32(f);
            t.members = rvecu32(f);
            t.label_code = r8(f);
            t.subtype = r32(f);
            t.spawn_seq = ru32(f);

            const uint32_t na = readLen(f);
            for (uint32_t k = 0; k < na; ++k) {
                const int32_t a = r32(f);
                t.markers_by_anchor[a] = readMarkerSet(f);
            }
            return t;
        }

        // ---- PER-TEMPLATE EXTRAS -----------------------------------------
        //
        // A SEPARATE BLOCK, appended after the banks, rather than fields added
        // to writeTemplate(). writeTemplate's record has no length prefix and
        // _templates.bin has no version field anywhere in it -- the whole format
        // rests on "old readers stop at end of file" -- so widening that record
        // makes new files misparse under the current reader and old files
        // misparse under the new one, in the middle of a section, with nothing
        // to detect it by. Appended, both directions degrade to defaults.
        //
        // WHAT IS IN HERE AND WHY EACH ONE MATTERS ON RELOAD:
        //
        //   confirmed_by_operator -- was persisted NOWHERE. A confirmed
        //     template reloaded with label_code set and confirmed() false, and
        //     confirmed() is what blocks a merge from collapsing two
        //     morphologies, what countLabeled counts, and what makes
        //     presumedCategory honour the operator's verdict over the
        //     algorithm's presumption. The operator's work survived the file
        //     only as a label nothing treated as authoritative.
        //
        //   members_clean -- which beats the waveform was averaged over. Absent,
        //     cleanCount() falls back to members.size(), excludedCount() reads
        //     zero, and a reader that honours kTemplatesAreCache and rebuilds
        //     averages the premature and Tukey-rejected beats back in.
        //
        //   operator_state -- the right-click quality verdict on this panel.
        //     Per template, unlike TemplateBin::bad_r_ch which is per bin. It
        //     goes in the extras rather than beside bad_r_ch because it is a new
        //     field on a record with no length prefix; see the header note.
        //
        //   the census counts -- n_premature_members and friends. Absent,
        //     presumedCategory() sees zero premature members and calls an
        //     ectopic template REGULAR, which hands it a landmark column it did
        //     not have when it was built.
        //
        // Nothing derived is stored: these are decisions and observations, not
        // recomputable from the templates.
        inline void writeTemplateExtras(std::ofstream& f,
            const tbank::BankTemplate& t) {
            w8(f, t.confirmed_by_operator ? 1u : 0u);
            wvecu32(f, t.members_clean);
            w32(f, t.n_premature_members);
            w32(f, t.n_voted_members);
            w32(f, t.n_noise_members);
            w32(f, t.n_tukey_members);
            w32(f, t.n_ppg_members);
            w8(f, t.operator_state);
        }

        inline void readTemplateExtras(std::ifstream& f,
            tbank::BankTemplate& t) {
            t.confirmed_by_operator = (r8(f) != 0);
            t.members_clean = rvecu32(f);
            t.n_premature_members = r32(f);
            t.n_voted_members = r32(f);
            t.n_noise_members = r32(f);
            t.n_tukey_members = r32(f);
            t.n_ppg_members = r32(f);
            t.operator_state = r8(f);
        }

        // Count-prefixed, and the count is checked against the bank rather than
        // trusted: an extras block whose length disagrees with the bank it
        // belongs to is a mismatched pair of sections, and applying it anyway
        // would attach one template's exclusions to another.
        inline void writeBankExtras(std::ofstream& f,
            const tbank::TemplateBank& b) {
            wu32(f, static_cast<uint32_t>(b.templates.size()));
            for (const auto& t : b.templates) writeTemplateExtras(f, t);
        }

        inline bool readBankExtras(std::ifstream& f, tbank::TemplateBank& b) {
            const uint32_t nt = readLen(f);
            if (nt != b.templates.size()) return false;
            for (uint32_t i = 0; i < nt; ++i)
                readTemplateExtras(f, b.templates[i]);
            return static_cast<bool>(f);
        }

        inline void writeBank(std::ofstream& f, const tbank::TemplateBank& b) {
            wu32(f, static_cast<uint32_t>(b.templates.size()));
            for (const auto& t : b.templates) writeTemplate(f, t);
            w32(f, b.effective_cap);
            wu32(f, b.next_spawn_seq);
            wu32(f, b.assigned_beats);
        }

        inline tbank::TemplateBank readBank(std::ifstream& f) {
            tbank::TemplateBank b;
            const uint32_t nt = readLen(f);
            b.templates.reserve(nt);
            for (uint32_t i = 0; i < nt; ++i) b.templates.push_back(readTemplate(f));
            b.effective_cap = r32(f);
            b.next_spawn_seq = ru32(f);
            b.assigned_beats = ru32(f);
            return b;
        }

        inline void writeTukeyCounts(std::ofstream& f,
            const tbank::TukeyPassCounts& c) {
            wu32(f, c.beats_in); wu32(f, c.beats_out); wu32(f, c.rejected);
            // Fence state is stored, not just the counts. A bin with a low
            // rejection rate and a wide IQR is contaminated; a bin with a low
            // rejection rate and a narrow IQR is clean. The counts alone cannot
            // tell those apart, and by the time anyone asks, the beats that
            // produced the fences may be long gone.
            wf64(f, c.q1); wf64(f, c.q3); wf64(f, c.fence_lo); wf64(f, c.fence_hi);
        }
        inline tbank::TukeyPassCounts readTukeyCounts(std::ifstream& f) {
            tbank::TukeyPassCounts c;
            c.beats_in = ru32(f); c.beats_out = ru32(f); c.rejected = ru32(f);
            c.q1 = rf64(f); c.q3 = rf64(f);
            c.fence_lo = rf64(f); c.fence_hi = rf64(f);
            return c;
        }

        inline void writeCounts(std::ofstream& f, const tbank::BinCounts& c) {
            wu32(f, c.beats_detected);
            wu32(f, c.n_regular); wu32(f, c.n_ectopic); wu32(f, c.n_noise);
            wu32(f, c.n_premature); wu32(f, c.n_vote_only);
            writeTukeyCounts(f, c.tukey_rr);
            writeTukeyCounts(f, c.tukey_amplitude);
            writeTukeyCounts(f, c.tukey_r_location);
            writeTukeyCounts(f, c.tukey_wave_score);
            wu32(f, c.n_regular_rejected_on_rr);
            wu32(f, c.n_spawns); wu32(f, c.n_merges); wu32(f, c.n_cap_raises);
            wu32(f, c.n_unscorable); wu32(f, c.n_reassigned_pass2);
        }
        inline tbank::BinCounts readCounts(std::ifstream& f) {
            tbank::BinCounts c;
            c.beats_detected = ru32(f);
            c.n_regular = ru32(f); c.n_ectopic = ru32(f); c.n_noise = ru32(f);
            c.n_premature = ru32(f); c.n_vote_only = ru32(f);
            c.tukey_rr = readTukeyCounts(f);
            c.tukey_amplitude = readTukeyCounts(f);
            c.tukey_r_location = readTukeyCounts(f);
            c.tukey_wave_score = readTukeyCounts(f);
            c.n_regular_rejected_on_rr = ru32(f);
            c.n_spawns = ru32(f); c.n_merges = ru32(f); c.n_cap_raises = ru32(f);
            c.n_unscorable = ru32(f); c.n_reassigned_pass2 = ru32(f);
            return c;
        }

        inline void writeSeed(std::ofstream& f, const seed_pool::SeedSelection& s) {
            // The basis is the point of storing this at all. A bin whose seed
            // fell back to ALL_USABLE has its features measured against a
            // reference that is not purely sinus, and that must be visible in
            // the archive rather than reconstructed later from the counts.
            w8(f, static_cast<uint8_t>(s.basis));
            wu32(f, s.n_candidates);
            wu32(f, s.n_excluded_premature);
            wu32(f, s.n_excluded_voted);
            wu32(f, s.n_excluded_category);
            wu32(f, s.n_selected);
            wf64(f, s.ectopic_fraction);
            w8(f, s.tukey_relocation_pending ? 1u : 0u);
            wvecu32(f, s.members);
        }
        inline seed_pool::SeedSelection readSeed(std::ifstream& f) {
            seed_pool::SeedSelection s;
            s.basis = static_cast<seed_pool::SeedBasis>(r8(f));
            s.n_candidates = ru32(f);
            s.n_excluded_premature = ru32(f);
            s.n_excluded_voted = ru32(f);
            s.n_excluded_category = ru32(f);
            s.n_selected = ru32(f);
            s.ectopic_fraction = rf64(f);
            s.tukey_relocation_pending = (r8(f) != 0u);
            s.members = rvecu32(f);
            return s;
        }

        inline void writeCapRaise(std::ofstream& f, const tbank::CapRaiseEvent& e) {
            wu64(f, e.bin_index);
            w32(f, e.channel);
            w32(f, e.old_cap); w32(f, e.new_cap);
            w32(f, e.template_a); w32(f, e.template_b);
            wf64(f, e.closeness);
            w8(f, e.label_a); w8(f, e.label_b);
        }
        inline tbank::CapRaiseEvent readCapRaise(std::ifstream& f) {
            tbank::CapRaiseEvent e;
            e.bin_index = ru64(f);
            e.channel = r32(f);
            e.old_cap = r32(f); e.new_cap = r32(f);
            e.template_a = r32(f); e.template_b = r32(f);
            e.closeness = rf64(f);
            e.label_a = r8(f); e.label_b = r8(f);
            return e;
        }

    }  // namespace detail

    // Writes the bank block for one bin, length-prefixed so that a reader of a
    // different v1 minor revision can skip the remainder instead of
    // misparsing. Call immediately after the legacy per-bin record.
    inline void writeBankBlock(std::ofstream& f, const BinBankState& st)
    {
        // Reserve the length, fill the body, then seek back. Streams here are
        // plain ofstreams on local files, so tellp/seekp are available; if this
        // ever needs to write to a non-seekable sink, build the body in a
        // std::string first and write the length ahead of it.
        const std::streampos len_pos = f.tellp();
        detail::wu32(f, 0u);
        const std::streampos body_start = f.tellp();

        for (int c = 0; c < 3; ++c) detail::writeBank(f, st.ecg[c]);
        detail::writeBank(f, st.ppg);
        for (int c = 0; c < 3; ++c) detail::writeSeed(f, st.seed[c]);

        detail::wu32(f, static_cast<uint32_t>(st.cap_raises.size()));
        for (const auto& e : st.cap_raises) detail::writeCapRaise(f, e);

        detail::writeCounts(f, st.counts);

        // ---- APPENDED, INSIDE THE LENGTH PREFIX --------------------------
        // No version bump for this. The block already carries body_bytes and
        // the reader below trusts it over its own field walk, which is exactly
        // the case this was built for: an OLD reader on a NEW block walks the
        // fields it knows, finds bytes left over, and seeks past them. Bumping
        // tbank_io::kVersionCurrent instead would make requireSupported() throw
        // on those files -- turning a format that degrades gracefully into one
        // that refuses to open.
        for (int c = 0; c < 3; ++c) detail::writeBankExtras(f, st.ecg[c]);
        detail::writeBankExtras(f, st.ppg);

        const std::streampos end = f.tellp();
        const uint32_t body_bytes =
            static_cast<uint32_t>(end - body_start);
        f.seekp(len_pos);
        detail::wu32(f, body_bytes);
        f.seekp(end);
    }

    // Reads one bank block. On a legacy file there is no block to read -- the
    // caller checks the header version and simply skips this, leaving slot 0 as
    // whatever the legacy per-channel template provided.
    inline BinBankState readBankBlock(std::ifstream& f)
    {
        BinBankState st;
        const uint32_t body_bytes = detail::ru32(f);
        const std::streampos body_start = f.tellg();

        for (int c = 0; c < 3; ++c) st.ecg[c] = detail::readBank(f);
        st.ppg = detail::readBank(f);
        for (int c = 0; c < 3; ++c) st.seed[c] = detail::readSeed(f);

        const uint32_t ne = detail::readLen(f);
        st.cap_raises.reserve(ne);
        for (uint32_t i = 0; i < ne; ++i)
            st.cap_raises.push_back(detail::readCapRaise(f));

        st.counts = detail::readCounts(f);

        // ---- EXTRAS, IF THIS FILE HAS THEM -------------------------------
        // Guarded on the length prefix rather than on a version, so a block
        // written before the extras existed reads correctly and leaves every
        // template's members_clean empty and confirmed_by_operator false -- the
        // same state such a file has always produced. A mismatched template
        // count aborts the extras and leaves the banks alone rather than
        // attaching one template's exclusions to another.
        if (f.tellg() - body_start
            < static_cast<std::streamoff>(body_bytes)) {
            bool ok = true;
            for (int c = 0; c < 3 && ok; ++c)
                ok = detail::readBankExtras(f, st.ecg[c]);
            if (ok) detail::readBankExtras(f, st.ppg);
        }

        // Trust the length prefix over the field-by-field walk. If a future
        // revision appended fields this reader does not know, the walk stops
        // short and the prefix carries it to the true end of the block; if the
        // walk overran, that is a format bug and must not be silently absorbed.
        const std::streamoff consumed = f.tellg() - body_start;
        if (consumed > static_cast<std::streamoff>(body_bytes))
            throw std::runtime_error(
                "template bank block: overran its declared length");
        f.seekg(body_start + static_cast<std::streamoff>(body_bytes));

        return st;
    }

    // ---------------------------------------------------------------------
    // Single-bank stream helpers, for embedding one TemplateBank inside another
    // file's format (template_io.cpp appends three per bin as a v3 trailing
    // section). Separate from writeBankBlock(): that one carries a whole bin's
    // worth of state with a length prefix, which is the wrong granularity when
    // the host file already has its own per-bin framing.
    //
    // readBankFromStream() returns false rather than throwing, because the
    // caller's contract is "a v1/v2 file simply ends here" -- a short read is
    // an expected outcome, not an error.
    // ---------------------------------------------------------------------

    // Public wrappers for the extras, for _templates.bin's trailing section.
    // That file has no version field anywhere and no length prefixes, so the
    // extras go in a section of their own after every existing one; see
    // template_io.cpp.
    inline void writeBankExtrasToStream(std::ofstream& f,
        const tbank::TemplateBank& b) {
        detail::writeBankExtras(f, b);
    }
    inline bool readBankExtrasFromStream(std::ifstream& f,
        tbank::TemplateBank& b) {
        try { return detail::readBankExtras(f, b); }
        catch (...) { return false; }
    }

    inline void writeBankToStream(std::ofstream& f, const tbank::TemplateBank& b) {
        detail::writeBank(f, b);
    }

    inline bool readBankFromStream(std::ifstream& f, tbank::TemplateBank& out) {
        try {
            out = detail::readBank(f);
        }
        catch (...) {
            return false;
        }
        return static_cast<bool>(f);
    }

    // A legacy bin: one template per channel, which is exactly a bank of size
    // one. Called when the header reports kVersionLegacy so that the rest of
    // the code never has to branch on file version again.
    inline tbank::TemplateBank bankFromLegacyTemplate(
        const std::vector<double>& tmpl,
        const std::vector<double>& tmpl_iqr,
        int r_col,
        uint32_t n_beats)
    {
        tbank::TemplateBank b;
        tbank::BankTemplate t;
        t.tmpl = tmpl;
        t.tmpl_iqr = tmpl_iqr;
        t.r_col = r_col;
        t.spawn_seq = 0;
        // Member indices are unrecoverable from a legacy file -- it stored the
        // median and a count, never which beats produced it. Left empty, with
        // assigned_beats carrying the count, so that anything needing members
        // fails loudly rather than reading an empty vector as "no beats".
        b.templates.push_back(std::move(t));
        b.next_spawn_seq = 1;
        b.assigned_beats = n_beats;
        return b;
    }

}  // namespace tbank_ser
