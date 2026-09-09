/**
 * @file   template_io.cpp
 * @brief  Implementation of the template-generation file I/O.
 */

#include "template_io.hpp"
#include "template_morphology_grouping/template_bank_serialize.hpp"

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <algorithm> 

namespace template_io {

    namespace {

        void writeVecD(std::ofstream& f, const std::vector<double>& v) {
            uint64_t sz = v.size();
            f.write(reinterpret_cast<const char*>(&sz), 8);
            if (sz > 0) f.write(reinterpret_cast<const char*>(v.data()), sz * 8);
        }

        void writeMethod(std::ofstream& f, const ChannelMethodTemplate& m) {
            writeVecD(f, m.ecgTemplate);
            // Empty for methods that don't compute std (sz=0, no payload).
            writeVecD(f, m.ecg_template_iqr);
            f.write(reinterpret_cast<const char*>(&m.alignment_point), 8);
            f.write(reinterpret_cast<const char*>(&m.r_col), 4);
        }

        bool read_template_vector(std::ifstream& f, std::vector<double>& v) {
            uint64_t sz;
            if (!f.read(reinterpret_cast<char*>(&sz), 8)) return false;
            if (sz > (1ull << 24)) return false;
            v.resize(sz);
            if (sz > 0) f.read(reinterpret_cast<char*>(v.data()), sz * 8);
            return static_cast<bool>(f);
        
        }

        bool readMethod(std::ifstream& f, ChannelMethodTemplate& m) {
            if (!read_template_vector(f, m.ecgTemplate)) return false;
            if (!read_template_vector(f, m.ecg_template_iqr)) return false;
            if (!f.read(reinterpret_cast<char*>(&m.alignment_point), 8)) return false;
            if (!f.read(reinterpret_cast<char*>(&m.r_col), 4)) return false;
            return true;
        }

    }  // anonymous namespace

    void write_template_binfile(const std::string& path, const TemplateFile& data) {
        std::ofstream f(path, std::ios::binary);
        if (!f) throw std::runtime_error("cannot create: " + path);

        char wbuf[1 << 16];
        f.rdbuf()->pubsetbuf(wbuf, sizeof(wbuf));

        uint64_t nBins = data.bins.size();
        f.write(reinterpret_cast<const char*>(&nBins), 8);

        for (const auto& b : data.bins) {
            writeMethod(f, b.ch1_raw); writeMethod(f, b.ch1_squared);
            writeMethod(f, b.ch1_absval); writeMethod(f, b.ch1_unfiltered);
            writeMethod(f, b.ch2_raw); writeMethod(f, b.ch2_squared);
            writeMethod(f, b.ch2_absval); writeMethod(f, b.ch2_unfiltered);
            writeMethod(f, b.ch3_raw); writeMethod(f, b.ch3_squared);
            writeMethod(f, b.ch3_absval); writeMethod(f, b.ch3_unfiltered);
            writeVecD(f, b.ppgTemplate);
            writeVecD(f, b.ppg_template_iqr);
            writeVecD(f, b.abpTemplate);
            writeVecD(f, b.abpTemplate_iqr);
            writeVecD(f, b.artTemplate);
            writeVecD(f, b.artTemplate_iqr);
            writeVecD(f, b.artPulmTemplate);
            writeVecD(f, b.artPulmTemplate_iqr);
            f.write(reinterpret_cast<const char*>(&b.ch1_n_beats_raw), 8);
            f.write(reinterpret_cast<const char*>(&b.ch2_n_beats_raw), 8);
            f.write(reinterpret_cast<const char*>(&b.ch3_n_beats_raw), 8);
            f.write(reinterpret_cast<const char*>(&b.ppg_n_beats), 8);
            f.write(reinterpret_cast<const char*>(&b.ppg_peak_col), 4);
            f.write(reinterpret_cast<const char*>(&b.ppg_onset_col), 4);
            uint8_t bad = b.bad_segment ? 1 : 0;
            f.write(reinterpret_cast<const char*>(&bad), 1);
        }

        // ---- v1 SECTION 2: per-anchor aligned templates ------------------
        // [uint64 nAnchors] then per anchor:
        //   [int32 anchorTag][uint64 nBinsForAnchor] then, per bin, 3x
        //   ChannelMethodTemplate (ch1_raw, ch2_raw, ch3_raw) via writeMethod.
        //
        {
            uint64_t nAnchors = data.raw_anchors.size();
            f.write(reinterpret_cast<const char*>(&nAnchors), 8);
            for (const auto& kv : data.raw_anchors) {
                int32_t tag = kv.first;
                f.write(reinterpret_cast<const char*>(&tag), 4);
                uint64_t nb = kv.second.size();
                f.write(reinterpret_cast<const char*>(&nb), 8);
                for (const auto& triplet : kv.second) {
                    writeMethod(f, triplet[0]);
                    writeMethod(f, triplet[1]);
                    writeMethod(f, triplet[2]);
                }
            }
        }

        // ---- v1 SECTION 3: Section 4.6 template bank, per bin, per lead ---
        // MUST come after the anchors section, unconditionally. The reader walks
        // these in order, so skipping or reordering one makes it parse this
        // count as another section's and produce plausible garbage.
        //
        // Layout: [uint64 nBinsWithBanks], then per such bin
        //         [uint64 binIndex][bank CH1][bank CH2][bank CH3].
        // Bins with no bank at all are skipped rather than written as zeros, so
        // a clean record costs almost nothing.
        {
            // THE SAME PREDICATE THE EXTRAS SECTION USES, and that is the fix
            // for the reload failing outright.
            //
            // This section listed a bin only when some ECG bank was non-empty,
            // while the extras section lists it when ANY of the four is -- PPG
            // included. So a bin with a pulse bank and no ECG bank got extras
            // written for ECG banks this section never wrote, the reader had
            // nothing to attach them to, readBankExtras found the counts
            // disagreed, and the stream desynchronised mid-section. On a record
            // where CH2 and CH3 contribute no beats that is not a corner case,
            // it is most of the file -- and the symptom was
            // length_error("vector too long") thrown from the LAST section,
            // dozens of kilobytes downstream of the cause.
            //
            // Listing the bin here writes three banks that may be empty --
            // eight bytes each -- and in exchange every extras block has a bank
            // to land on. Narrowing the extras section to ECG-only would also
            // resynchronise the two, but it would silently drop the PPG
            // confirmations for exactly those bins.
            auto hasAnyBank = [](const template_io::BinTemplates& b) {
                for (int c = 0; c < 3; ++c)
                    if (!b.ecg_bank[c].templates.empty()) return true;
                return !b.ppg_bank.templates.empty();
                };

            uint64_t nWithBanks = 0;
            for (const auto& b : data.bins) if (hasAnyBank(b)) ++nWithBanks;
            f.write(reinterpret_cast<const char*>(&nWithBanks), 8);

            for (uint64_t i = 0; i < data.bins.size(); ++i) {
                const auto& b = data.bins[i];
                if (!hasAnyBank(b)) continue;   // must match the count above
                f.write(reinterpret_cast<const char*>(&i), 8);
                for (int c = 0; c < 3; ++c)
                    tbank_ser::writeBankToStream(f, b.ecg_bank[c]);
            }
        }

        // ---- v1 SECTION 4: Section 4.6 PPG bank, per bin ------------------
        // MUST come after the ECG bank section, unconditionally, for the same
        // reason that one follows the anchors.
        //
        // Layout: [uint64 nBinsWithPpgBank], then per such bin
        //         [uint64 binIndex][bank].
        // Bins with no PPG bank are skipped rather than written as zeros, so a
        // record with no PPG costs 8 bytes total.
        {
            uint64_t nWithPpg = 0;
            for (const auto& b : data.bins)
                if (!b.ppg_bank.templates.empty()) ++nWithPpg;
            f.write(reinterpret_cast<const char*>(&nWithPpg), 8);

            for (uint64_t i = 0; i < data.bins.size(); ++i) {
                const auto& b = data.bins[i];
                if (b.ppg_bank.templates.empty()) continue;
                f.write(reinterpret_cast<const char*>(&i), 8);
                tbank_ser::writeBankToStream(f, b.ppg_bank);
            }
        }

        // ---- v1 SECTION 5: per-template extras, per bin -------------------
        // WHY A SEPARATE SECTION RATHER THAN WIDER TEMPLATE RECORDS. Nothing
        // here is length-prefixed, so a reader cannot skip a field it does not
        // know about. Keeping the extras in their own block means the bank
        // records stay one fixed shape that writer and reader agree on
        // literally, rather than two shapes that have to agree by convention.
        //
        // WHAT IS IN IT: confirmed_by_operator -- which was persisted nowhere,
        // so every operator confirmation was lost on reload and with it the
        // merge protection, the polymorphy count and the operator's override of
        // presumedCategory -- plus members_clean and the per-template census.
        //
        // MUST come after the PPG bank section, unconditionally. The reader walks
        // in order, so a conditional or reordered section makes it read one
        // section's count as another's.
        //
        // Layout: [uint64 nBinsWithExtras], then per such bin
        //         [uint64 binIndex][extras CH1][extras CH2][extras CH3][extras PPG].
        {
            auto hasBank = [](const template_io::BinTemplates& b) {
                for (int c = 0; c < 3; ++c)
                    if (!b.ecg_bank[c].templates.empty()) return true;
                return !b.ppg_bank.templates.empty();
                };

            uint64_t nWithExtras = 0;
            for (const auto& b : data.bins) if (hasBank(b)) ++nWithExtras;
            f.write(reinterpret_cast<const char*>(&nWithExtras), 8);

            for (uint64_t i = 0; i < data.bins.size(); ++i) {
                const auto& b = data.bins[i];
                if (!hasBank(b)) continue;
                f.write(reinterpret_cast<const char*>(&i), 8);
                for (int c = 0; c < 3; ++c)
                    tbank_ser::writeBankExtrasToStream(f, b.ecg_bank[c]);
                tbank_ser::writeBankExtrasToStream(f, b.ppg_bank);
            }
        }

        // ---- v1 SECTION 6: per-anchor bank slot templates -----------------
        // MUST come after the extras section, unconditionally: the reader walks
        // these in order, so a skipped section makes it parse the next one's
        // count as this one's.
        //
        // Layout: [uint64 nAnchors], then per anchor
        //         [int32 tag][uint64 nBins], then per bin, per 3 channels
        //         [uint32 nSlots], then per slot
        //         [vecd tmpl][vecd tmpl_iqr][uint32 n_members].
        {
            uint64_t nAnchors = data.bank_anchors.size();
            f.write(reinterpret_cast<const char*>(&nAnchors), 8);
            for (const auto& kv : data.bank_anchors) {
                int32_t tag = kv.first;
                f.write(reinterpret_cast<const char*>(&tag), 4);
                uint64_t nb = kv.second.size();
                f.write(reinterpret_cast<const char*>(&nb), 8);
                for (const auto& perBin : kv.second) {
                    for (int c = 0; c < 3; ++c) {
                        uint32_t nSlots = static_cast<uint32_t>(perBin[c].size());
                        f.write(reinterpret_cast<const char*>(&nSlots), 4);
                        for (const auto& st : perBin[c]) {
                            writeVecD(f, st.tmpl);
                            writeVecD(f, st.tmpl_iqr);
                            uint32_t nm = st.n_members;
                            f.write(reinterpret_cast<const char*>(&nm), 4);
                        }
                    }
                }
            }
        }
    }



    TemplateFile read_template_binfile(const std::string& path) {
        std::ifstream f(path, std::ios::binary);
        if (!f) throw std::runtime_error("cannot open: " + path);

        char rbuf[1 << 16];
        f.rdbuf()->pubsetbuf(rbuf, sizeof(rbuf));

        TemplateFile out;

        uint64_t nBins = 0;
        if (!f.read(reinterpret_cast<char*>(&nBins), 8))
            throw std::runtime_error("template file truncated: " + path);
        out.bins.resize(nBins);

        for (auto& b : out.bins) {
            if (!readMethod(f, b.ch1_raw) || !readMethod(f, b.ch1_squared) ||
                !readMethod(f, b.ch1_absval) || !readMethod(f, b.ch1_unfiltered) ||
                !readMethod(f, b.ch2_raw) || !readMethod(f, b.ch2_squared) ||
                !readMethod(f, b.ch2_absval) || !readMethod(f, b.ch2_unfiltered) ||
                !readMethod(f, b.ch3_raw) || !readMethod(f, b.ch3_squared) ||
                !readMethod(f, b.ch3_absval) || !readMethod(f, b.ch3_unfiltered) ||
                !read_template_vector(f, b.ppgTemplate) ||
                !read_template_vector(f, b.ppg_template_iqr) ||
                !read_template_vector(f, b.abpTemplate) ||
                !read_template_vector(f, b.abpTemplate_iqr) ||
                !read_template_vector(f, b.artTemplate) ||
                !read_template_vector(f, b.artTemplate_iqr) ||
                !read_template_vector(f, b.artPulmTemplate) ||
                !read_template_vector(f, b.artPulmTemplate_iqr))
                throw std::runtime_error("template file truncated mid-bin: " + path);
            if (!f.read(reinterpret_cast<char*>(&b.ch1_n_beats_raw), 8) ||
                !f.read(reinterpret_cast<char*>(&b.ch2_n_beats_raw), 8) ||
                !f.read(reinterpret_cast<char*>(&b.ch3_n_beats_raw), 8) ||
                !f.read(reinterpret_cast<char*>(&b.ppg_n_beats), 8))
                throw std::runtime_error("template file truncated (missing n_beats fields): " + path);
            if (!f.read(reinterpret_cast<char*>(&b.ppg_peak_col), 4) ||
                !f.read(reinterpret_cast<char*>(&b.ppg_onset_col), 4))
                throw std::runtime_error("template file truncated (missing ppg fiducials): " + path);
            uint8_t bad = 0;
            f.read(reinterpret_cast<char*>(&bad), 1);
            b.bad_segment = (bad != 0);
        }

        // ---- v1 SECTION 2: per-anchor aligned templates ------------------
        // FATAL ON A SHORT READ. There is one format version and this code has
        // not shipped, so every section the writer emits is present in every
        // file that exists. A count that will not read is a TRUNCATED OR
        // CORRUPT file, and the soft `if (f.read(...))` this replaced treated
        // it as "an older file that ends here" -- which is how a desync came to
        // present as length_error from a vector constructor five sections
        // later instead of as an error naming this one.
        uint64_t nAnchors = 0;
        if (!f.read(reinterpret_cast<char*>(&nAnchors), 8))
            throw std::runtime_error("template file truncated at anchor count: " + path);
        if (nAnchors > 8)
            throw std::runtime_error("template file has implausible anchor count: " + path);
        {
            for (uint64_t a = 0; a < nAnchors; ++a) {
                int32_t tag = 0;
                if (!f.read(reinterpret_cast<char*>(&tag), 4)) break;
                uint64_t nb = 0;
                if (!f.read(reinterpret_cast<char*>(&nb), 8)) break;
                // AN UNCHECKED COUNT IS NOT A COUNT. nb went straight from
                // disk into a vector constructor, so a stream one byte out
                // of position reached the allocator with garbage and threw
                // length_error("vector too long") -- reaching the caller
                // naming neither the cause nor the section. An anchor block
                // is per bin, so it cannot exceed the bin count this file
                // already declared.
                if (nb == 0 || nb > out.bins.size()) break;
                std::vector<std::array<ChannelMethodTemplate, 3>> perBin(nb);
                bool ok = true;
                for (uint64_t i = 0; i < nb && ok; ++i) {
                    ok = readMethod(f, perBin[i][0])
                        && readMethod(f, perBin[i][1])
                        && readMethod(f, perBin[i][2]);
                }
                if (!ok) break;   // truncated anchor section: keep what parsed cleanly
                out.raw_anchors[static_cast<int>(tag)] = std::move(perBin);
            }
        }

        // ---- v1 SECTION 3: template bank, per bin, per lead ---------------
        // Fatal on a short read, for the reason given at the anchors section.
        {
            uint64_t nWithBanks = 0;
            if (!f.read(reinterpret_cast<char*>(&nWithBanks), 8))
                throw std::runtime_error("template file truncated at bank count: " + path);
            if (nWithBanks > out.bins.size())
                throw std::runtime_error("bank section lists more bins than the file has: " + path);
            {
                for (uint64_t k = 0; k < nWithBanks; ++k) {
                    uint64_t bi = 0;
                    if (!f.read(reinterpret_cast<char*>(&bi), 8)) break;
                    bool ok = true;
                    for (int c = 0; c < 3 && ok; ++c) {
                        tbank::TemplateBank bank;
                        ok = tbank_ser::readBankFromStream(f, bank);
                        if (ok && bi < out.bins.size())
                            out.bins[bi].ecg_bank[c] = std::move(bank);
                    }
                    if (!ok) break;   // truncated: keep what parsed cleanly
                }
            }
        }

        // ---- v1 SECTION 4: PPG bank, per bin ------------------------------
        {
            uint64_t nWithPpg = 0;
            if (!f.read(reinterpret_cast<char*>(&nWithPpg), 8))
                throw std::runtime_error("template file truncated at PPG bank count: " + path);
            if (nWithPpg > out.bins.size())
                throw std::runtime_error("PPG bank section lists more bins than the file has: " + path);
            {
                for (uint64_t k = 0; k < nWithPpg; ++k) {
                    uint64_t bi = 0;
                    if (!f.read(reinterpret_cast<char*>(&bi), 8)) break;
                    tbank::TemplateBank bank;
                    if (!tbank_ser::readBankFromStream(f, bank)) break;
                    if (bi < out.bins.size())
                        out.bins[bi].ppg_bank = std::move(bank);
                }
            }
        }

        // ---- v1 SECTION 5: per-template extras, per bin -------------------
        //
        // The extras are applied ONTO the banks read above, so this must run
        // after both bank sections. A template-count mismatch inside
        // readBankExtrasFromStream stops the section rather than applying it:
        // extras that do not line up with their bank would attach one
        // template's exclusions and confirmation to another, which is worse
        // than not having them.
        {
            uint64_t nWithExtras = 0;
            if (f.read(reinterpret_cast<char*>(&nWithExtras), 8)) {
                for (uint64_t k = 0; k < nWithExtras; ++k) {
                    uint64_t bi = 0;
                    if (!f.read(reinterpret_cast<char*>(&bi), 8)) break;
                    if (bi >= out.bins.size()) break;
                    bool ok = true;
                    for (int c = 0; c < 3 && ok; ++c)
                        ok = tbank_ser::readBankExtrasFromStream(
                            f, out.bins[bi].ecg_bank[c]);
                    if (ok)
                        ok = tbank_ser::readBankExtrasFromStream(
                            f, out.bins[bi].ppg_bank);
                    if (!ok) break;
                }
            }
        }

        // ---- v1 SECTION 6: per-anchor bank slot templates -----------------
        {
            uint64_t nAnchors = 0;
            if (!f.read(reinterpret_cast<char*>(&nAnchors), 8))
                throw std::runtime_error("template file truncated at slot-anchor count: " + path);
            if (nAnchors > 8)
                throw std::runtime_error("slot-anchor section has implausible count: " + path);
            {
                for (uint64_t a = 0; a < nAnchors; ++a) {
                    int32_t tag = 0;
                    if (!f.read(reinterpret_cast<char*>(&tag), 4)) break;
                    uint64_t nb = 0;
                    if (!f.read(reinterpret_cast<char*>(&nb), 8)) break;
                    // AN UNCHECKED COUNT IS NOT A COUNT. nb went straight from
                    // disk into a vector constructor, so a stream one byte out
                    // of position reached the allocator with garbage and threw
                    // length_error("vector too long") -- reaching the caller
                    // naming neither the cause nor the section. An anchor block
                    // is per bin, so it cannot exceed the bin count this file
                    // already declared.
                    if (nb == 0 || nb > out.bins.size()) break;
                    std::vector<std::array<std::vector<TemplateFile::BankSlotTemplate>, 3>> perBin(nb);
                    bool ok = true;
                    for (uint64_t i = 0; i < nb && ok; ++i) {
                        for (int c = 0; c < 3 && ok; ++c) {
                            uint32_t nSlots = 0;
                            if (!f.read(reinterpret_cast<char*>(&nSlots), 4)) { ok = false; break; }
                            // Implausible-count guard: a corrupt length must
                            // stop the walk, not make the reader allocate on it.
                            if (nSlots > 4096u) { ok = false; break; }
                            perBin[i][c].resize(nSlots);
                            for (uint32_t sl = 0; sl < nSlots && ok; ++sl) {
                                if (!read_template_vector(f, perBin[i][c][sl].tmpl)) { ok = false; break; }
                                if (!read_template_vector(f, perBin[i][c][sl].tmpl_iqr)) { ok = false; break; }
                                uint32_t nm = 0;
                                if (!f.read(reinterpret_cast<char*>(&nm), 4)) { ok = false; break; }
                                perBin[i][c][sl].n_members = nm;
                            }
                        }
                    }
                    if (!ok) break;   // truncated: keep what parsed cleanly
                    out.bank_anchors[static_cast<int>(tag)] = std::move(perBin);
                }
            }
        }

        return out;
    }
}  // namespace template_io