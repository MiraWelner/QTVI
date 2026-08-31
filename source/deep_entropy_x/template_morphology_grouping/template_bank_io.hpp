#pragma once
/**
 * @file   template_bank_io.hpp
 * @brief  Format versioning for the template-markings binary, so Section 4.6
 *         bank state can be appended without silently misparsing files
 *         written before it existed.
 *
 *         THE PROBLEM. writeTemplateMarkingsBin() writes an 8-byte bin count
 *         and then bin records, with no magic and no version field. Adding a
 *         single field to the per-bin record therefore changes the stride with
 *         nothing to detect it by: an old reader on a new file, or a new
 *         reader on an old file, walks off alignment a few bins in and
 *         produces plausible-looking garbage rather than an error. Marker
 *         positions would land in the wrong bins. There is no version to bump
 *         because there was never one to begin with.
 *
 *         THE FIX, and why it is backward compatible. A legacy file opens with
 *         a bin count -- a small integer, tens to low thousands, so the top
 *         bytes of that uint64 are zero. A new file opens with 8 ASCII bytes,
 *         which read as a uint64 in the 5e18 range: an impossible bin count.
 *         So the first 8 bytes disambiguate the two formats with no ambiguity
 *         and no separate sidecar. Legacy files stay readable forever; new
 *         files are self-describing.
 *
 *         Version history is kept in this file rather than a comment in the
 *         writer, so that a reader handling several versions has one place to
 *         look.
 */

#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>

namespace tbank_io {

    // "TMKB" + version digits. Chosen so the little-endian uint64 read of these
    // bytes is astronomically larger than any real bin count.
    inline constexpr char kMagic[9] = "TMKBANK1";

    enum : uint32_t {
        // No magic, no version. Header is a bare uint64 bin count. Per-bin
        // record as of the pre-4.6 writer.
        kVersionLegacy = 0,

        // Adds, per bin: a per-channel template count, then for each bank
        // member its median, IQR, r_col, member beat indices, label code,
        // subtype, spawn_seq, and per-anchor BankMarkerSet. Adds per-bin
        // BinCounts and the CapRaiseEvent list. Template 0 of each channel
        // carries what ChannelTemplateData/MarkerSet carried in legacy files,
        // so a legacy file loads as a one-template-per-channel bank and needs
        // no migration pass.
        kVersionBank = 1,

        kVersionCurrent = kVersionBank
    };

    struct FormatHeader {
        uint32_t version = kVersionLegacy;
        uint64_t n_bins = 0;
    };

    // Reads whichever header the file has and leaves the stream positioned at
    // the first bin record. Throws only on truncation.
    inline FormatHeader readHeader(std::ifstream& f) {
        FormatHeader h;

        char probe[8] = {};
        f.read(probe, 8);
        if (f.gcount() != 8)
            throw std::runtime_error("template markings bin: truncated header");

        if (std::memcmp(probe, kMagic, 8) == 0) {
            uint32_t v = 0;
            f.read(reinterpret_cast<char*>(&v), 4);
            uint64_t n = 0;
            f.read(reinterpret_cast<char*>(&n), 8);
            if (!f)
                throw std::runtime_error("template markings bin: truncated versioned header");
            h.version = v;
            h.n_bins = n;
            return h;
        }

        // Legacy: the 8 bytes we just consumed WERE the bin count.
        uint64_t n = 0;
        std::memcpy(&n, probe, 8);
        h.version = kVersionLegacy;
        h.n_bins = n;
        return h;
    }

    inline void writeHeader(std::ofstream& f, uint64_t n_bins,
        uint32_t version = kVersionCurrent) {
        f.write(kMagic, 8);
        f.write(reinterpret_cast<const char*>(&version), 4);
        f.write(reinterpret_cast<const char*>(&n_bins), 8);
    }

    // A version the running build does not know how to parse is a hard error,
    // not a best-effort read: silently ignoring trailing fields is how marker
    // positions end up in the wrong bins.
    inline void requireSupported(const FormatHeader& h) {
        if (h.version > kVersionCurrent)
            throw std::runtime_error(
                "template markings bin: file version "
                + std::to_string(h.version)
                + " is newer than this build supports ("
                + std::to_string(static_cast<uint32_t>(kVersionCurrent)) + ")");
    }

}  // namespace tbank_io
