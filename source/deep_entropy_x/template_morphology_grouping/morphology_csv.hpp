#pragma once
/**
 * @file   morphology_csv.hpp
 * @brief  Section 4.5-4.6 outputs, two granularities x two encodings:
 *
 *           <stem>_beats.csv      one column per BEAT       (text)
 *           <stem>_templates.csv  one column per TEMPLATE   (text)
 *           <stem>_beats.bin      the same, binary
 *           <stem>_templates.bin  the same, binary
 *
 *         The CSVs are descriptor rows on top, then waveform rows below,
 *         sharing one column layout so a column can be sliced and read whole.
 *         The .bin files carry identical content COLUMN-MAJOR: each column's
 *         descriptors immediately followed by its samples, so one seek reads one
 *         beat or one template whole.
 *
 *         WHY A BINARY FORM AT ALL. beats.csv is one column per beat, so a
 *         16-bin record runs to roughly 16,000 columns and every sample is
 *         written as decimal text -- the file is large, slow to parse, and past
 *         Excel's 16,384-column ceiling anyway. The .bin is a fifth the size,
 *         parses in one pass, and round-trips doubles exactly, which the CSV
 *         does not: text formatting loses low-order bits, so a CSV round-trip
 *         cannot reproduce a median bit-for-bit.
 *
 *         DESCRIPTORS ARE CODES, NOT STRINGS, in the binary form. The CSV spells
 *         out "premature" and "ectopic" for a human; the .bin stores the enum
 *         values those words came from. Nothing is lost -- the same enums
 *         generate both -- and it removes the parse step where a reader has to
 *         match strings that a later spelling change would silently break.
 *
 *         WHY TWO FILES AND NOT ONE. The descriptors mean different things at
 *         the two granularities and cannot be folded. Per BEAT, category /
 *         premature / tukey are the actual verdicts on that beat. Per TEMPLATE
 *         they can only be summaries of its members, which is a different claim:
 *         a template whose row 4 reads `premature` holds a majority of premature
 *         beats, not a premature beat. Putting both in one file would make row 4
 *         mean two things depending on which section you were in.
 *
 *         WHY THIS IS WRITTEN HERE AND NOT BY THE VIEWER. The viewer's
 *         writeAlignedTemplateCsv() previously produced this file, in long
 *         format: one row per SAMPLE per bin, columns being signals. It cannot
 *         produce a per-beat file, because the viewer has no beats -- it loads
 *         TemplateBin, which holds templates only, and there is no m_beats or
 *         BeatsFile anywhere in it. The beat matrices exist on the generation
 *         side (alignment's aligned.beats), so the writer lives here.
 *
 *         That makes the viewer's writeAlignedTemplateCsv() dead for this file
 *         and it must be retired, or the two writers will each clobber the
 *         other depending on run order.
 *
 *         Each column is one beat. Read top to bottom:
 *
 *           row 1  category    pqrst | ectopic | noise
 *           row 2  bin         bin index (bin length set by the config file)
 *           row 3  template    the template this beat was assigned to:
 *                              PQRST_A, PQRST_B ... or PVC_A, PVC_B ...
 *           row 4  premature   premature | vote | no
 *           row 5  tukey       removed | kept | not_eligible
 *
 *         WHY PER BEAT AND NOT PER TEMPLATE. Four of the five rows ARE per-beat
 *         verdicts. Category comes from the operator's mark on that beat,
 *         prematurity from the timing test on that beat's RR, and the Tukey
 *         outcome from that beat's position in the distribution. Only row 3 is a
 *         template property, and it is the assignment -- which template this
 *         beat landed in. Aggregating any of the other four to the template
 *         level would report a majority and discard the disagreements, and the
 *         disagreements are the informative part.
 *
 *         WHY ROWS 3 AND 4 MAY DISAGREE, AND MUST. 4.6 never reassigns a
 *         category; it only excludes. Row 3 is the CLASS -- what the beat is,
 *         which is the operator's judgment propagated through the template. Row
 *         4 is a REFERENCE-SET FLAG -- whether the timing filter excluded it
 *         from the sinus reference. A beat in PVC_A that reads `no` on row 4 is
 *         a late ventricular beat, or a beat inside a run where the trailing-ten
 *         median had already collapsed and the ratio stopped firing. A beat that
 *         reads `premature` but sits in a PQRST template is very likely a PAC:
 *         early, but conducted normally, so it looks like sinus. Reconciling the
 *         two rows would erase exactly these cases.
 *
 *         WHY ROW 1 IS THREE VALUES AND NOT FIVE. The five 4.5 categories come
 *         later; today the operator's markers support three. A `noise` column
 *         means the user marked that beat Minor Noise -- noise that does not
 *         disturb R-peak detection. Signal marked R Peak Noise cannot appear
 *         here at all: that mark suppresses detection, so no beats exist in
 *         those spans to have a column.
 *
 *         WHY ROW 3 SAYS PQRST WHEN NOTHING IS MARKED. Because it is not a
 *         label. "A template with no confirmed member stays unlabeled" -- PQRST
 *         is how an unlabeled template presents, and the letter after it is the
 *         algorithm's separation index, not a class. The `confirmed` row below
 *         keeps presumption and confirmation in separate fields, because only
 *         confirmed labels are usable as training data later and a training set
 *         that silently absorbed presumptions would be learning the algorithm's
 *         own guesses.
 *
 *         WIDTH WARNING. A 15-minute bin holds roughly a thousand beats, so a
 *         16-bin record is around 16,000 columns per channel. Excel stops at
 *         16,384 columns and will refuse or truncate; pandas, R and awk are
 *         fine. Channels are written as separate row BLOCKS rather than extra
 *         columns for this reason -- tripling the width would put every record
 *         past that limit.
 */

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "bin_pipeline.hpp"

namespace morphology_csv {

    inline std::string g_dir;
    inline std::string g_stem;

    inline void set(const std::string& dir, const std::string& stem,
        const std::string& subdir = "") {
        g_stem = stem;
        if (dir.empty()) { g_dir.clear(); return; }
        std::filesystem::path p(dir);
        if (!subdir.empty()) p /= subdir;
        std::error_code ec;
        std::filesystem::create_directories(p, ec);
        // ofstream on a missing directory fails silently and the writer would
        // return having produced nothing, so the destination is created here
        // rather than assumed.
        g_dir = ec ? dir : p.string();
    }

    namespace detail {

        // EVERY enumerator named explicitly, no catch-all. The `default:
        // return "pqrst"` this replaces mapped anything unrecognised to the
        // word for NORMAL, so when Category gained UNCONFIRMED every
        // unconfirmed template began writing itself into the archive as pqrst
        // -- the exact inference Section 4.6 forbids, arriving through a switch
        // default rather than through a decision anyone made. A missing case is
        // now a compiler warning instead of a plausible-looking column.
        inline const char* categoryWord(tbank::Category c) {
            switch (c) {
            case tbank::Category::REGULAR:     return "pqrst";
            case tbank::Category::ECTOPIC:     return "ectopic";
            case tbank::Category::NOISE:       return "noise";
            case tbank::Category::UNCONFIRMED: return "unconfirmed";
            }
            return "unconfirmed";   // unreachable; the safe direction if reached
        }

        // Row 3: the class, then an underscore and the letter the ALGORITHM
        // assigned when it separated the morphologies.
        //
        // For a confirmed template the letter follows the subtype index the bank
        // issued (PVC_A, PVC_B); for an unlabeled one it follows spawn order.
        // The two differ because subtypes are issued per class in order of first
        // confirmation, while spawn order is order of first appearance -- and
        // spawn_seq is used rather than the template's current position because
        // a merge erases an element and shifts everything after it, which would
        // otherwise renumber templates between runs.
        inline std::string templateName(const tbank::BankTemplate& t) {
            std::string cls = "PQRST";
            int letterIdx = -1;

            if (t.label_code != tbank::kUnlabeled) {
                switch (t.label_code) {
                case tbank::kCodePvc:        cls = "PVC";   break;
                case tbank::kCodePac:        cls = "PAC";   break;
                case tbank::kCodeVt:         cls = "VT";    break;
                case tbank::kCodeMinorNoise: cls = "NOISE"; break;
                default: cls = "CODE" + std::to_string(t.label_code); break;
                }
                if (t.subtype > 0) letterIdx = t.subtype - 1;
            }
            if (letterIdx < 0) letterIdx = static_cast<int>(t.spawn_seq);
            const char letter = static_cast<char>('A' + (letterIdx % 26));
            return cls + "_" + std::string(1, letter);
        }

        inline const char* prematureWord(tbank::PvcFilter p) {
            switch (p) {
            case tbank::PvcFilter::PREMATURE: return "premature";
            case tbank::PvcFilter::VOTE:      return "vote";
            default:                          return "no";
            }
        }

        // Row 5. Tukey runs only on beats NOT flagged premature: a premature
        // beat is already excluded from the reference set, so there is nothing
        // for Tukey to decide about it. `not_eligible` is therefore a real and
        // common state, distinct from "tested and kept" -- and the specific
        // rejection reason is kept because a beat classified pqrst and rejected
        // on RR LENGTH is very likely ectopy the classifier missed, which is the
        // cheapest estimate of classifier recall available.
        inline const char* tukeyWord(tbank::TukeyOutcome t) {
            switch (t) {
            case tbank::TukeyOutcome::KEPT:           return "kept";
            case tbank::TukeyOutcome::REJ_RR_LENGTH:  return "removed_rr";
            case tbank::TukeyOutcome::REJ_AMPLITUDE:  return "removed_amplitude";
            case tbank::TukeyOutcome::REJ_R_LOCATION: return "removed_r_location";
            case tbank::TukeyOutcome::REJ_WAVE_SCORE: return "removed_wave_score";
            default:                                  return "not_eligible";
            }
        }

    }  // namespace detail

    struct ChannelBlock {
        const char* channel = "CH1";
        const std::vector<bin_pipeline::ChannelOutput>* per_bin = nullptr;

        // Aligned beat matrices, [bin][beat][sample], on the shared axis. Same
        // indexing as per_bin[bin].flags, so column k of bin b is
        // beats[b][k] described by flags[k]. Null writes descriptors only.
        const std::vector<std::vector<std::vector<double>>>* beats = nullptr;

        // R column on the shared axis, per bin. Emitted as a descriptor row so
        // a reader can align columns to a fiducial without re-detecting it --
        // every beat's detected R sits at this column by construction.
        const std::vector<int>* r_col = nullptr;
    };

    // One column per beat. Channels are emitted as successive row blocks, each
    // preceded by a `channel` row naming it, so a reader can split on that row
    // rather than parsing three files.
    inline bool writeBeats(const std::vector<ChannelBlock>& blocks) {
        if (g_dir.empty() || g_stem.empty()) return false;
        std::ofstream f(g_dir + "/" + g_stem + "_beats.csv", std::ios::trunc);
        if (!f) return false;

        for (const auto& blk : blocks) {
            if (!blk.per_bin || blk.per_bin->empty()) continue;

            struct Col {
                std::string category, bin, name, premature, tukey, confirmed;
            };
            // Column -> (bin, beat) so the waveform rows below can find the
            // samples for the same column the descriptors describe.
            struct ColSrc { size_t bin; size_t beat; };
            std::vector<Col> cols;
            std::vector<ColSrc> colSrc;

            for (size_t b = 0; b < blk.per_bin->size(); ++b) {
                const bin_pipeline::ChannelOutput& out = (*blk.per_bin)[b];
                const size_t n = out.flags.size();
                for (size_t i = 0; i < n; ++i) {
                    Col c;
                    c.category = detail::categoryWord(out.flags[i].category);
                    c.bin = std::to_string(b);
                    c.premature = detail::prematureWord(out.flags[i].pvc);
                    c.tukey = detail::tukeyWord(out.flags[i].tukey);

                    // Assignment. kUnscorable (-2) means the beat had too little
                    // axis overlap for a correlation to mean anything -- it was
                    // deliberately NOT given a template, because a beat that
                    // cannot be compared is not evidence of a new morphology.
                    const int32_t t = (i < out.assignment.size())
                        ? out.assignment[i] : -1;
                    if (t == tbank::kUnscorable) {
                        c.name = "unscorable";
                        c.confirmed = "na";
                    }
                    else if (t >= 0 && t < out.bank.size()) {
                        const tbank::BankTemplate& tp = out.bank.templates[t];
                        c.name = detail::templateName(tp);
                        c.confirmed = tp.confirmed() ? "confirmed" : "presumed";
                    }
                    else {
                        c.name = "unassigned";
                        c.confirmed = "na";
                    }
                    cols.push_back(std::move(c));
                    colSrc.push_back({ b, i });
                }
            }
            if (cols.empty()) continue;

            auto row = [&](const char* label, std::string Col::* field) {
                f << label;
                for (auto& c : cols) f << ',' << (c.*field);
                f << '\n';
                };

            f << "channel";
            for (size_t k = 0; k < cols.size(); ++k) f << ',' << blk.channel;
            f << '\n';

            if (blk.r_col) {
                f << "r_col";
                for (const auto& src : colSrc)
                    f << ',' << ((src.bin < blk.r_col->size())
                        ? (*blk.r_col)[src.bin] : -1);
                f << '\n';
            }

            row("category", &Col::category);
            row("bin", &Col::bin);
            row("template", &Col::name);
            row("premature", &Col::premature);
            row("tukey", &Col::tukey);
            // Confirmation in its OWN row, separate from the class name in row
            // 3: PQRST_A is a presumed class, and only confirmed labels are
            // usable as training data.
            row("confirmed", &Col::confirmed);

            // ---- waveform: one row per sample index, one column per beat ----
            // Below the descriptors the same columns continue, so a reader can
            // slice a beat by column and get its descriptors and its samples
            // together. Row label is the sample index on the shared axis; the
            // r_col descriptor row above says where R sits, so the axis needs no
            // separate time column.
            //
            // NaN is written as an empty cell. The shared axis is NaN-padded at
            // both ends -- alignment lays every beat on one axis and fills the
            // overhang -- and writing "nan" or 0 would both be wrong: 0 is a
            // real amplitude, and a literal makes every consumer special-case a
            // string in a numeric column.
            if (blk.beats) {
                size_t width = 0;
                for (const auto& binBeats : *blk.beats)
                    for (const auto& bt : binBeats) width = std::max(width, bt.size());

                for (size_t sIdx = 0; sIdx < width; ++sIdx) {
                    f << sIdx;
                    for (const auto& col : colSrc) {
                        f << ',';
                        if (col.bin < blk.beats->size()) {
                            const auto& binBeats = (*blk.beats)[col.bin];
                            if (col.beat < binBeats.size()
                                && sIdx < binBeats[col.beat].size()) {
                                const double v = binBeats[col.beat][sIdx];
                                if (!std::isnan(v)) f << v;
                            }
                        }
                    }
                    f << '\n';
                }
            }
        }
        return true;
    }

    // ---------------------------------------------------------------------
    // <stem>_templates.csv -- one column per TEMPLATE
    //
    // Same rows, but every descriptor is now a SUMMARY over the template's
    // members rather than a verdict on one beat, and the distinction is real:
    // `premature` on a template row means a majority of its members were
    // flagged, not that the template is premature. Thresholded at a majority
    // rather than at any, because a sinus template picks up the occasional
    // premature beat and one such beat must not relabel eight hundred others.
    // ---------------------------------------------------------------------

    namespace detail {

        // Row 4, aggregated. `vote` and `mixed` are kept distinct from
        // `premature` because the two filter paths have disjoint competence: the
        // raw test catches isolated ectopy, the 5-of-8 vote only fires inside
        // runs, and in perfect bigeminy exactly 4 of any 8 beats are ectopic so
        // the vote never reaches 5 at all. A template reading `vote` was
        // therefore found by a mechanism that cannot see alternating rhythms.
        inline std::string prematureWordAgg(const tbank::BankTemplate& t) {
            const uint32_t n = static_cast<uint32_t>(t.members.size());
            if (n == 0) return "na";
            const uint32_t prem = t.n_premature_members;
            const uint32_t vote = t.n_voted_members;
            if (prem * 2 > n) return "premature";
            if ((prem + vote) * 2 > n) return "vote";
            if (prem + vote > 0) return "mixed";
            return "no";
        }

        // Row 5, aggregated. `not_eligible` dominates on ectopic templates by
        // construction: Tukey runs only on beats NOT flagged premature, since a
        // premature beat is already excluded from the reference set and there is
        // nothing left for Tukey to decide about it.
        inline std::string tukeyWordAgg(const tbank::BankTemplate& t,
            const bin_pipeline::ChannelOutput& out)
        {
            uint32_t eligible = 0, removed = 0;
            for (uint32_t m : t.members) {
                if (m >= out.flags.size()) continue;
                const tbank::TukeyOutcome o = out.flags[m].tukey;
                if (o == tbank::TukeyOutcome::NOT_ELIGIBLE) continue;
                ++eligible;
                if (o != tbank::TukeyOutcome::KEPT) ++removed;
            }
            if (eligible == 0) return "not_eligible";
            if (removed == 0) return "kept";
            if (removed == eligible) return "removed";
            return "partial";
        }

    }  // namespace detail

    inline bool writeTemplates(const std::vector<ChannelBlock>& blocks) {
        if (g_dir.empty() || g_stem.empty()) return false;
        std::ofstream f(g_dir + "/" + g_stem + "_templates.csv", std::ios::trunc);
        if (!f) return false;

        for (const auto& blk : blocks) {
            if (!blk.per_bin || blk.per_bin->empty()) continue;

            struct Col {
                std::string category, bin, name, premature, tukey, confirmed,
                    members, share, marking;
                const std::vector<double>* wave = nullptr;
                size_t binIdx = 0;
            };
            std::vector<Col> cols;

            for (size_t b = 0; b < blk.per_bin->size(); ++b) {
                const bin_pipeline::ChannelOutput& out = (*blk.per_bin)[b];
                for (int t = 0; t < out.bank.size(); ++t) {
                    const tbank::BankTemplate& tp = out.bank.templates[t];
                    if (tp.tmpl.empty()) continue;
                    Col c;
                    c.category = detail::categoryWord(tp.presumedCategory());
                    c.bin = std::to_string(b);
                    c.name = detail::templateName(tp);
                    c.premature = detail::prematureWordAgg(tp);
                    c.tukey = detail::tukeyWordAgg(tp, out);
                    c.confirmed = tp.confirmed() ? "confirmed" : "presumed";
                    c.members = std::to_string(tp.memberCount());
                    c.share = std::to_string(out.bank.beatShare(t));
                    // Only category 1 templates are landmark-marked: only
                    // category 1 beats feed feature extraction, so a P-onset on
                    // a PVC template has nothing downstream to consume it and a
                    // PVC's QT is not comparable to a sinus QT.
                    c.marking = tp.wantsLandmarkMarking() ? "landmark" : "class_only";
                    c.wave = &tp.tmpl;
                    c.binIdx = b;
                    cols.push_back(std::move(c));
                }
            }
            if (cols.empty()) continue;

            auto row = [&](const char* label, std::string Col::* field) {
                f << label;
                for (auto& c : cols) f << ',' << (c.*field);
                f << '\n';
                };

            f << "channel";
            for (size_t k = 0; k < cols.size(); ++k) f << ',' << blk.channel;
            f << '\n';
            if (blk.r_col) {
                f << "r_col";
                for (auto& c : cols)
                    f << ',' << ((c.binIdx < blk.r_col->size())
                        ? (*blk.r_col)[c.binIdx] : -1);
                f << '\n';
            }

            row("category", &Col::category);
            row("bin", &Col::bin);
            row("template", &Col::name);
            row("premature", &Col::premature);
            row("tukey", &Col::tukey);
            row("confirmed", &Col::confirmed);
            row("marking", &Col::marking);
            row("n_members", &Col::members);
            row("beat_share", &Col::share);

            // Waveform: the template median, one row per sample index. Empty
            // cell for NaN -- the shared axis is NaN-padded at both ends, and 0
            // is a real amplitude.
            size_t width = 0;
            for (auto& c : cols) width = std::max(width, c.wave->size());
            for (size_t sIdx = 0; sIdx < width; ++sIdx) {
                f << sIdx;
                for (auto& c : cols) {
                    f << ',';
                    if (sIdx < c.wave->size() && !std::isnan((*c.wave)[sIdx]))
                        f << (*c.wave)[sIdx];
                }
                f << '\n';
            }
        }
        return true;
    }

    // =====================================================================
    // Binary forms
    // =====================================================================
    //
    // Layout, both files:
    //
    //   [char[8] magic][uint32 version][uint32 nChannelBlocks]
    //   per block:
    //     [uint32 nameLen][char[nameLen] channel]
    //     [uint32 width]              samples per column, 0 = no waveform
    //     [uint64 nColumns]
    //     per column: a fixed descriptor record, then `width` doubles
    //
    // Magic first so a truncated or wrong-type file fails immediately rather
    // than being read as a plausible column count -- the failure mode the
    // markings .bin has, where a bare count opens the file and any 8 bytes look
    // like a valid header.
    //
    // NaN is written as NaN, not as a sentinel. The CSV writes an empty cell
    // because a literal would force consumers to special-case a string in a
    // numeric column; binary has no such problem, and IEEE NaN survives a
    // read/write pair exactly.

    inline constexpr char kBeatsMagic[9] = "DEXBEAT1";
    inline constexpr char kTemplatesMagic[9] = "DEXTMPL1";
    inline constexpr uint32_t kBinVersion = 1;

    // One record per beat column. Fixed size, so a reader can stride over
    // descriptors without parsing them.
    struct BeatRecord {
        uint32_t bin = 0;
        uint8_t  category = 0;   // tbank::Category
        uint8_t  premature = 0;   // tbank::PvcFilter
        uint8_t  tukey = 0;   // tbank::TukeyOutcome
        uint8_t  confirmed = 0;   // 0 presumed, 1 confirmed, 2 n/a
        uint8_t  label_code = 0;   // annotation code, 0 = unlabeled
        uint8_t  letter = 0;   // separation index: 0 = A
        int32_t  template_id = -1;  // bank slot, -2 unscorable, -1 unassigned
        int32_t  r_col = -1;
    };

    // One record per template column.
    struct TemplateRecord {
        uint32_t bin = 0;
        uint8_t  category = 0;   // presumed tbank::Category
        uint8_t  premature = 0;   // 0 no, 1 mixed, 2 vote, 3 premature
        uint8_t  tukey = 0;   // 0 not_eligible, 1 kept, 2 partial, 3 removed
        uint8_t  confirmed = 0;
        uint8_t  label_code = 0;
        uint8_t  letter = 0;
        uint8_t  landmark_marked = 0;   // 1 = category 1, gets landmarks
        int32_t  template_id = -1;
        int32_t  r_col = -1;
        uint32_t n_members = 0;
        double   beat_share = 0.0;
    };

    namespace detail {

        inline void w(std::ofstream& f, const void* p, size_t n) {
            f.write(reinterpret_cast<const char*>(p), static_cast<std::streamsize>(n));
        }
        inline bool r(std::ifstream& f, void* p, size_t n) {
            return static_cast<bool>(
                f.read(reinterpret_cast<char*>(p), static_cast<std::streamsize>(n)));
        }

        inline void writeBlockHeader(std::ofstream& f, const char* channel,
            uint32_t width, uint64_t nCols) {
            const uint32_t len = static_cast<uint32_t>(std::string(channel).size());
            w(f, &len, 4);
            w(f, channel, len);
            w(f, &width, 4);
            w(f, &nCols, 8);
        }

        // Letter index: the subtype for a confirmed template, spawn order for an
        // unlabeled one. They differ because subtypes are issued per class in
        // order of first confirmation while spawn order is order of first
        // appearance -- and spawn_seq is used rather than bank position because a
        // merge erases an element and shifts everything after it, which would
        // otherwise renumber templates between runs.
        inline uint8_t letterOf(const tbank::BankTemplate& t) {
            int idx = (t.label_code != tbank::kUnlabeled && t.subtype > 0)
                ? t.subtype - 1 : static_cast<int>(t.spawn_seq);
            return static_cast<uint8_t>(idx % 26);
        }

        inline uint8_t prematureCode(const std::string& word) {
            if (word == "premature") return 3;
            if (word == "vote")      return 2;
            if (word == "mixed")     return 1;
            return 0;
        }
        inline uint8_t tukeyAggCode(const std::string& word) {
            if (word == "removed") return 3;
            if (word == "partial") return 2;
            if (word == "kept")    return 1;
            return 0;
        }

    }  // namespace detail

    inline bool writeBeatsBin(const std::vector<ChannelBlock>& blocks) {
        if (g_dir.empty() || g_stem.empty()) return false;
        std::ofstream f(g_dir + "/" + g_stem + "_beats.bin",
            std::ios::binary | std::ios::trunc);
        if (!f) return false;

        uint32_t nBlocks = 0;
        for (const auto& b : blocks) if (b.per_bin && !b.per_bin->empty()) ++nBlocks;

        detail::w(f, kBeatsMagic, 8);
        detail::w(f, &kBinVersion, 4);
        detail::w(f, &nBlocks, 4);

        for (const auto& blk : blocks) {
            if (!blk.per_bin || blk.per_bin->empty()) continue;

            uint32_t width = 0;
            uint64_t nCols = 0;
            for (size_t b = 0; b < blk.per_bin->size(); ++b) {
                nCols += (*blk.per_bin)[b].flags.size();
                if (blk.beats && b < blk.beats->size())
                    for (const auto& bt : (*blk.beats)[b])
                        width = std::max(width, static_cast<uint32_t>(bt.size()));
            }
            detail::writeBlockHeader(f, blk.channel, width, nCols);

            const double nan = std::numeric_limits<double>::quiet_NaN();
            for (size_t b = 0; b < blk.per_bin->size(); ++b) {
                const bin_pipeline::ChannelOutput& out = (*blk.per_bin)[b];
                for (size_t i = 0; i < out.flags.size(); ++i) {
                    BeatRecord rec;
                    rec.bin = static_cast<uint32_t>(b);
                    rec.category = static_cast<uint8_t>(out.flags[i].category);
                    rec.premature = static_cast<uint8_t>(out.flags[i].pvc);
                    rec.tukey = static_cast<uint8_t>(out.flags[i].tukey);
                    rec.r_col = (blk.r_col && b < blk.r_col->size())
                        ? (*blk.r_col)[b] : -1;

                    const int32_t t = (i < out.assignment.size())
                        ? out.assignment[i] : -1;
                    rec.template_id = t;
                    if (t >= 0 && t < out.bank.size()) {
                        const tbank::BankTemplate& tp = out.bank.templates[t];
                        rec.label_code = tp.label_code;
                        rec.letter = detail::letterOf(tp);
                        rec.confirmed = tp.confirmed() ? 1 : 0;
                    }
                    else {
                        rec.confirmed = 2;   // n/a: unassigned or unscorable
                    }
                    detail::w(f, &rec, sizeof(rec));

                    for (uint32_t s = 0; s < width; ++s) {
                        double v = nan;
                        if (blk.beats && b < blk.beats->size()) {
                            const auto& binBeats = (*blk.beats)[b];
                            if (i < binBeats.size() && s < binBeats[i].size())
                                v = binBeats[i][s];
                        }
                        detail::w(f, &v, 8);
                    }
                }
            }
        }
        return static_cast<bool>(f);
    }

    inline bool writeTemplatesBin(const std::vector<ChannelBlock>& blocks) {
        if (g_dir.empty() || g_stem.empty()) return false;
        std::ofstream f(g_dir + "/" + g_stem + "_templates.bin",
            std::ios::binary | std::ios::trunc);
        if (!f) return false;

        uint32_t nBlocks = 0;
        for (const auto& b : blocks) if (b.per_bin && !b.per_bin->empty()) ++nBlocks;

        detail::w(f, kTemplatesMagic, 8);
        detail::w(f, &kBinVersion, 4);
        detail::w(f, &nBlocks, 4);

        for (const auto& blk : blocks) {
            if (!blk.per_bin || blk.per_bin->empty()) continue;

            uint32_t width = 0;
            uint64_t nCols = 0;
            for (const auto& out : *blk.per_bin)
                for (int t = 0; t < out.bank.size(); ++t) {
                    if (out.bank.templates[t].tmpl.empty()) continue;
                    ++nCols;
                    width = std::max(width,
                        static_cast<uint32_t>(out.bank.templates[t].tmpl.size()));
                }
            detail::writeBlockHeader(f, blk.channel, width, nCols);

            const double nan = std::numeric_limits<double>::quiet_NaN();
            for (size_t b = 0; b < blk.per_bin->size(); ++b) {
                const bin_pipeline::ChannelOutput& out = (*blk.per_bin)[b];
                for (int t = 0; t < out.bank.size(); ++t) {
                    const tbank::BankTemplate& tp = out.bank.templates[t];
                    if (tp.tmpl.empty()) continue;

                    TemplateRecord rec;
                    rec.bin = static_cast<uint32_t>(b);
                    rec.category = static_cast<uint8_t>(tp.presumedCategory());
                    rec.premature = detail::prematureCode(detail::prematureWordAgg(tp));
                    rec.tukey = detail::tukeyAggCode(detail::tukeyWordAgg(tp, out));
                    rec.confirmed = tp.confirmed() ? 1 : 0;
                    rec.label_code = tp.label_code;
                    rec.letter = detail::letterOf(tp);
                    rec.landmark_marked = tp.wantsLandmarkMarking() ? 1 : 0;
                    rec.template_id = t;
                    rec.r_col = (blk.r_col && b < blk.r_col->size())
                        ? (*blk.r_col)[b] : -1;
                    rec.n_members = static_cast<uint32_t>(tp.memberCount());
                    rec.beat_share = out.bank.beatShare(t);
                    detail::w(f, &rec, sizeof(rec));

                    for (uint32_t s = 0; s < width; ++s) {
                        const double v = (s < tp.tmpl.size()) ? tp.tmpl[s] : nan;
                        detail::w(f, &v, 8);
                    }
                }
            }
        }
        return static_cast<bool>(f);
    }

    // ---------------------------------------------------------------------
    // Readers. Present so the format has exactly one definition rather than
    // one here and another in whatever consumes it -- a binary layout with a
    // writer and no reader is a layout nobody can verify.
    // ---------------------------------------------------------------------

    template <class Rec>
    struct BinBlock {
        std::string channel;
        uint32_t width = 0;
        std::vector<Rec> records;
        std::vector<double> samples;   // records.size() * width, column-major

        const double* column(size_t k) const {
            return (k < records.size() && width) ? &samples[k * width] : nullptr;
        }
    };

    template <class Rec>
    inline bool readBin(const std::string& path, const char* magic,
        std::vector<BinBlock<Rec>>& out)
    {
        std::ifstream f(path, std::ios::binary);
        if (!f) return false;
        char m[8] = {};
        if (!detail::r(f, m, 8)) return false;
        if (std::memcmp(m, magic, 8) != 0) return false;   // wrong type or truncated
        uint32_t ver = 0, nBlocks = 0;
        if (!detail::r(f, &ver, 4) || !detail::r(f, &nBlocks, 4)) return false;
        if (ver > kBinVersion) return false;   // newer than this build understands

        out.clear();
        for (uint32_t bi = 0; bi < nBlocks; ++bi) {
            BinBlock<Rec> blk;
            uint32_t len = 0;
            if (!detail::r(f, &len, 4) || len > 64) return false;
            blk.channel.resize(len);
            if (len && !detail::r(f, blk.channel.data(), len)) return false;
            uint64_t nCols = 0;
            if (!detail::r(f, &blk.width, 4) || !detail::r(f, &nCols, 8)) return false;

            blk.records.resize(nCols);
            blk.samples.assign(static_cast<size_t>(nCols) * blk.width, 0.0);
            for (uint64_t k = 0; k < nCols; ++k) {
                if (!detail::r(f, &blk.records[k], sizeof(Rec))) return false;
                if (blk.width && !detail::r(f,
                    &blk.samples[static_cast<size_t>(k) * blk.width],
                    static_cast<size_t>(blk.width) * 8)) return false;
            }
            out.push_back(std::move(blk));
        }
        return true;
    }

    inline bool readBeatsBin(const std::string& path,
        std::vector<BinBlock<BeatRecord>>& out) {
        return readBin<BeatRecord>(path, kBeatsMagic, out);
    }
    inline bool readTemplatesBin(const std::string& path,
        std::vector<BinBlock<TemplateRecord>>& out) {
        return readBin<TemplateRecord>(path, kTemplatesMagic, out);
    }

}  // namespace morphology_csv