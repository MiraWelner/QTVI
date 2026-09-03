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

        bool readVecD(std::ifstream& f, std::vector<double>& v) {
            uint64_t sz;
            if (!f.read(reinterpret_cast<char*>(&sz), 8)) return false;
            v.resize(sz);
            if (sz > 0) f.read(reinterpret_cast<char*>(v.data()), sz * 8);
            return static_cast<bool>(f);
        }

        bool readMethod(std::ifstream& f, ChannelMethodTemplate& m) {
            if (!readVecD(f, m.ecgTemplate)) return false;
            if (!readVecD(f, m.ecg_template_iqr)) return false;
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

        // ---- trailing per-anchor section (v2) --------------------------
        // [uint64 nAnchors] then per anchor:
        //   [int32 anchorTag][uint64 nBinsForAnchor] then, per bin, 3x
        //   ChannelMethodTemplate (ch1_raw, ch2_raw, ch3_raw) via writeMethod.
        // Old readers stop after the bins loop above and never see this, so
        // a v2 file still reads as R-only under the old reader.
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

        // ---- v3: Section 4.6 template bank, per bin, per ECG channel -------
        // Another TRAILING OPTIONAL SECTION, on exactly the contract the v2
        // anchors block above established: a v1/v2 reader stops at end-of-file,
        // and a v3 reader that finds nothing here leaves the banks empty. No
        // magic, no version byte, no migration -- an absent bank reads correctly
        // as "one template per channel", which is a bank of size one.
        //
        // MUST come after the anchors count, unconditionally. The reader walks
        // these sections in order, so skipping or reordering either one makes it
        // parse the bank count as an anchor count and produce plausible garbage.
        //
        // Layout: [uint64 nBinsWithBanks], then per such bin
        //         [uint64 binIndex][bank CH1][bank CH2][bank CH3].
        // Bins with no bank at all are skipped rather than written as zeros, so
        // a clean record costs almost nothing.
        {
            uint64_t nWithBanks = 0;
            for (const auto& b : data.bins)
                for (int c = 0; c < 3; ++c)
                    if (!b.ecg_bank[c].templates.empty()) { ++nWithBanks; break; }
            f.write(reinterpret_cast<const char*>(&nWithBanks), 8);

            for (uint64_t i = 0; i < data.bins.size(); ++i) {
                const auto& b = data.bins[i];
                bool any = false;
                for (int c = 0; c < 3; ++c)
                    if (!b.ecg_bank[c].templates.empty()) { any = true; break; }
                if (!any) continue;
                f.write(reinterpret_cast<const char*>(&i), 8);
                for (int c = 0; c < 3; ++c)
                    tbank_ser::writeBankToStream(f, b.ecg_bank[c]);
            }
        }

        // ---- v4: Section 4.6 PPG bank, per bin -----------------------------
        // A FOURTH TRAILING OPTIONAL SECTION, on the contract the v2 anchors and
        // v3 ECG banks established: a v1/v2/v3 reader stops at end-of-file, and
        // a v4 reader that finds nothing here leaves ppg_bank empty -- which
        // reads correctly as "one pulse template per bin", exactly what a v3
        // file holds.
        //
        // MUST come after the ECG bank section, unconditionally, for the same
        // reason that one had to follow the anchors: the reader walks these in
        // order, so reordering or conditionally skipping either makes it parse
        // one section's count as another's and produce plausible garbage.
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

        // ---- v5: per-template extras, per bin ------------------------------
        // A FIFTH TRAILING OPTIONAL SECTION, on the same contract as the v2
        // anchors, the v3 ECG banks and the v4 PPG bank: a v1-v4 reader stops at
        // end of file, and a v5 reader that finds nothing leaves every
        // template's extras at their defaults -- which is precisely the state a
        // v4 file has always loaded into.
        //
        // WHY A NEW SECTION RATHER THAN WIDER TEMPLATE RECORDS. This file has no
        // version field and no length prefixes; the entire format rests on old
        // readers hitting EOF. Adding fields to the bank records themselves
        // would make new files misparse under the current reader and old files
        // misparse under the new one, part-way through a section, with nothing
        // in either file able to detect it.
        //
        // WHAT IS IN IT: confirmed_by_operator -- which was persisted nowhere,
        // so every operator confirmation was lost on reload and with it the
        // merge protection, the polymorphy count and the operator's override of
        // presumedCategory -- plus members_clean and the per-template census.
        //
        // MUST come after the v4 count, unconditionally. The reader walks these
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
                !readVecD(f, b.ppgTemplate) ||
                !readVecD(f, b.ppg_template_iqr) ||
                !readVecD(f, b.abpTemplate) ||
                !readVecD(f, b.abpTemplate_iqr) ||
                !readVecD(f, b.artTemplate) ||
                !readVecD(f, b.artTemplate_iqr) ||
                !readVecD(f, b.artPulmTemplate) ||
                !readVecD(f, b.artPulmTemplate_iqr))
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

        // ---- trailing per-anchor section (v2, optional) ----------------
        // Old (v1) files end after the last bin, so a failed/short read here
        // just means "no extra anchors" -- leave raw_anchors empty and return
        // the R-only file. Only a CLEAN read of the count populates anchors.
        uint64_t nAnchors = 0;
        if (f.read(reinterpret_cast<char*>(&nAnchors), 8)) {
            for (uint64_t a = 0; a < nAnchors; ++a) {
                int32_t tag = 0;
                if (!f.read(reinterpret_cast<char*>(&tag), 4)) break;
                uint64_t nb = 0;
                if (!f.read(reinterpret_cast<char*>(&nb), 8)) break;
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

        // ---- v3: template bank (trailing, optional) ------------------------
        // A v1/v2 file simply ends here, so a failed read of the count means
        // "no banks" rather than an error -- the same contract as the anchors
        // section above. Anything already parsed stands.
        {
            uint64_t nWithBanks = 0;
            if (f.read(reinterpret_cast<char*>(&nWithBanks), 8)) {
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

        // ---- v4: PPG bank (trailing, optional) -----------------------------
        // A v1/v2/v3 file ends here, so a failed read of the count means "no
        // PPG bank" rather than an error -- the same contract as the two
        // sections above. Anything already parsed stands.
        {
            uint64_t nWithPpg = 0;
            if (f.read(reinterpret_cast<char*>(&nWithPpg), 8)) {
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

        // ---- v5: per-template extras (trailing, optional) ------------------
        // A v1-v4 file ends here, so a failed read of the count means "no
        // extras" rather than an error, on the same contract as the three
        // sections above. Anything already parsed stands.
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

        return out;
    }
}  // namespace template_io