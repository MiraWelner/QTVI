#pragma once
/**
 * @file   template_bank_serialize.hpp
 * @brief  Per-bin bank state on disk: the Section 4.6 banks and their
 *         per-template extras, embedded in _templates.bin's trailing sections.
 *
 *         THE HOST FILE PROVIDES THE FRAMING. template_io.cpp appends three
 *         banks per bin as a v3 trailing section and the extras as a section of
 *         its own after that. This header owns the record layout inside those
 *         sections and nothing else -- no length prefixes, no version field of
 *         its own. _templates.bin has no version field anywhere in it; the
 *         whole format rests on "an older reader stops at end of file", which
 *         is why a new section is appended after every existing one rather than
 *         widening a record that is already there.
 *
 *         WHAT IS NOT STORED. Nothing derived. Beat shares are recomputed from
 *         member counts, the polymorphic verdict is recomputed from label
 *         codes, and template waveforms are recomputed from member indices plus
 *         the beats. Storing a median alongside the members that produce it
 *         invites the two to disagree, and the members are the smaller of the
 *         two. Templates ARE written, because the viewer needs them before the
 *         beats are loaded, but they are a rebuildable cache: a reader that has
 *         the beats should prefer recomputing over what it read, and any
 *         disagreement is a bug in whoever wrote the file.
 *
 *         NO MARKER POSITIONS. markers_by_anchor is deliberately absent from
 *         every record here -- see writeTemplate. Landmarks live in
 *         _template_markings.bin, which nothing but the marking session writes.
 */

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>   // std::to_string, for the extras mismatch message
#include <vector>

#include "template_bank.hpp"

namespace tbank_ser {

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
    // Records
    // ---------------------------------------------------------------------

    namespace detail {

        inline void writeTemplate(std::ofstream& f, const tbank::BankTemplate& t) {
            wvecd(f, t.tmpl);
            wvecd(f, t.tmpl_iqr);
            w32(f, t.r_col);
            wvecu32(f, t.members);
            w8(f, t.label_code);
            w32(f, t.subtype);
            wu32(f, t.spawn_seq);

            // markers_by_anchor IS DELIBERATELY NOT WRITTEN HERE.
            //
            // This file is _templates.bin, and the PIPELINE owns it: it is
            // rewritten by prepareViewerJob, again by finalizeViewerJob, and
            // once per anchor by regenerateWithAnchor. Operator landmarks stored
            // in it were therefore destroyed by the next "Finish and Next" --
            // silently, because the rewrite is a normal part of the anchor cycle.
            //
            // Landmarks are operator judgement and live in _template_markings.bin,
            // which nothing but the marking session writes. This file carries
            // morphology: the waveform, its spread, its members, its class.
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

            // No marker block to read -- see writeTemplate. markers_by_anchor is
            // left empty here and filled by readTemplateMarkingsBin, or seeded
            // fresh by seed_bank_template when no marking file exists.
            return t;
        }

        // ---- PER-TEMPLATE EXTRAS -----------------------------------------
        //
        // A SEPARATE SECTION, appended after the banks, rather than fields
        // added to writeTemplate(). writeTemplate's record has no length prefix
        // and _templates.bin has no version field anywhere in it -- the whole
        // format rests on "old readers stop at end of file" -- so widening that
        // record makes new files misparse under the current reader and old files
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
        //     zero, and a reader that treats the stored template as a cache and
        //     rebuilds averages the premature and Tukey-rejected beats back in.
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
            w8(f, t.marked_invalid_template);
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
            t.marked_invalid_template = r8(f);
        }

        // Count-prefixed, and the count is checked against the bank rather than
        // trusted: an extras section whose length disagrees with the bank it
        // belongs to is a mismatched pair of sections, and applying it anyway
        // would attach one template's exclusions to another.
        inline void writeBankExtras(std::ofstream& f,
            const tbank::TemplateBank& b) {
            wu32(f, static_cast<uint32_t>(b.templates.size()));
            for (const auto& t : b.templates) writeTemplateExtras(f, t);
        }

        inline bool readBankExtras(std::ifstream& f, tbank::TemplateBank& b) {
            const uint32_t nt = readLen(f);
            if (nt != b.templates.size())
                throw std::runtime_error(
                    "bank extras count " + std::to_string(nt)
                    + " does not match bank size "
                    + std::to_string(b.templates.size())
                    + " -- the bank and extras sections of template_io.cpp"
                    " disagree about which bins carry banks");
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

    }  // namespace detail

    // ---------------------------------------------------------------------
    // Single-bank stream helpers, for embedding one TemplateBank inside another
    // file's format. template_io.cpp appends three per bin as a v3 trailing
    // section, and the extras as a section after that.
    //
    // readBankFromStream() returns false rather than throwing, because the
    // caller's contract is "a v1/v2 file simply ends here" -- a short read is
    // an expected outcome, not an error.
    // ---------------------------------------------------------------------

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

    inline void writeBankExtrasToStream(std::ofstream& f,
        const tbank::TemplateBank& b) {
        detail::writeBankExtras(f, b);
    }

    inline bool readBankExtrasFromStream(std::ifstream& f,
        tbank::TemplateBank& b) {
        try { return detail::readBankExtras(f, b); }
        catch (...) { return false; }
    }

}  // namespace tbank_ser