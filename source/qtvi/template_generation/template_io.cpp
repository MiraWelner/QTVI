/**
 * @file   template_io.cpp
 * @brief  Implementation of the template-generation file I/O.
 */

#include "template_io.hpp"

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
            writeVecD(f, m.ecgTemplate_std);
            f.write(reinterpret_cast<const char*>(&m.alignment_point), 8);
            f.write(reinterpret_cast<const char*>(&m.avg_r_expand), 8);
        }

        void writeAveraged(std::ofstream& f, const AveragedTemplate& a) {
            writeVecD(f, a.waveform);
            f.write(reinterpret_cast<const char*>(&a.n_contributing), 8);
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
            if (!readVecD(f, m.ecgTemplate_std)) return false;
            if (!f.read(reinterpret_cast<char*>(&m.alignment_point), 8)) return false;
            if (!f.read(reinterpret_cast<char*>(&m.avg_r_expand), 8)) return false;
            return true;
        }

        bool readAveraged(std::ifstream& f, AveragedTemplate& a) {
            if (!readVecD(f, a.waveform)) return false;
            if (!f.read(reinterpret_cast<char*>(&a.n_contributing), 8)) return false;
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
            writeVecD(f, b.ppgTemplate_std);
            writeVecD(f, b.abpTemplate);
            writeVecD(f, b.abpTemplate_std);
            writeVecD(f, b.artTemplate);
            writeVecD(f, b.artTemplate_std);
            writeVecD(f, b.artPulmTemplate);
            writeVecD(f, b.artPulmTemplate_std);
            uint8_t bad = b.bad_segment ? 1 : 0;
            f.write(reinterpret_cast<const char*>(&bad), 1);
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
                !readVecD(f, b.ppgTemplate_std) ||
                !readVecD(f, b.abpTemplate) ||
                !readVecD(f, b.abpTemplate_std) ||
                !readVecD(f, b.artTemplate) ||
                !readVecD(f, b.artTemplate_std) ||
                !readVecD(f, b.artPulmTemplate) ||
                !readVecD(f, b.artPulmTemplate_std))
                throw std::runtime_error("template file truncated mid-bin: " + path);
            uint8_t bad = 0;
            f.read(reinterpret_cast<char*>(&bad), 1);
            b.bad_segment = (bad != 0);
        }

        return out;
    }

    void write_beats_binfile(const std::string& path, const BeatsFile& data) {
        std::ofstream out(path, std::ios::binary);
        if (!out) {
            throw std::runtime_error("write_beats_binfile: cannot open " + path);
        }

        const uint64_t n_bins = static_cast<uint64_t>(data.per_bin_beats.size());
        out.write(reinterpret_cast<const char*>(&n_bins), sizeof(n_bins));

        for (uint64_t i = 0; i < n_bins; ++i) {
            const uint8_t bad = (i < data.bad_segment.size() && data.bad_segment[i]) ? 1 : 0;
            out.write(reinterpret_cast<const char*>(&bad), sizeof(bad));

            const auto& beats = data.per_bin_beats[i];
            const uint64_t n_beats = static_cast<uint64_t>(beats.size());
            out.write(reinterpret_cast<const char*>(&n_beats), sizeof(n_beats));

            for (const auto& beat : beats) {
                const uint64_t sz = static_cast<uint64_t>(beat.size());
                out.write(reinterpret_cast<const char*>(&sz), sizeof(sz));
                if (sz > 0) {
                    out.write(reinterpret_cast<const char*>(beat.data()),
                        sz * sizeof(double));
                }
            }
        }
    }

    void write_template_csvfile(const std::string& path, const TemplateFile& data,
        const std::string& fileID, double /*sampleRateHz*/) {
        std::ofstream f(path);
        if (!f) throw std::runtime_error("cannot create: " + path);

        // One row per sample. Templates within a bin can differ in length,
        // so each column is filled up to its own length and blank beyond it.
        // Only raw and ppg carry a _std column (the other methods never
        // compute std).
        f << "file_id,bin_num,"
            "ch1_raw_mv,ch1_raw_std,ch1_squared_mv,ch1_absval_mv,"
            "ch2_raw_mv,ch2_raw_std,ch2_squared_mv,ch2_absval_mv,"
            "ch3_raw_mv,ch3_raw_std,ch3_squared_mv,ch3_absval_mv,"
            "ppg_mv,ppg_std\n";

        f << std::setprecision(10);

        // Write v[i] if present and finite, else leave the cell blank.
        auto cell = [&](const std::vector<double>& v, size_t i) {
            if (i < v.size() && !std::isnan(v[i])) f << v[i];
            };

        for (size_t bi = 0; bi < data.bins.size(); ++bi) {
            const BinTemplates& b = data.bins[bi];

            struct Col { const std::vector<double>* mv; const std::vector<double>* sd; };
            const Col cols[] = {
                {&b.ch1_raw.ecgTemplate,     &b.ch1_raw.ecgTemplate_std},
                {&b.ch1_squared.ecgTemplate, nullptr},
                {&b.ch1_absval.ecgTemplate,  nullptr},
                {&b.ch2_raw.ecgTemplate,     &b.ch2_raw.ecgTemplate_std},
                {&b.ch2_squared.ecgTemplate, nullptr},
                {&b.ch2_absval.ecgTemplate,  nullptr},
                {&b.ch3_raw.ecgTemplate,     &b.ch3_raw.ecgTemplate_std},
                {&b.ch3_squared.ecgTemplate, nullptr},
                {&b.ch3_absval.ecgTemplate,  nullptr},
                {&b.ppgTemplate,             &b.ppgTemplate_std},
            };

            size_t maxLen = 0;
            for (const Col& c : cols) {
                maxLen = std::max(maxLen, c.mv->size());
                if (c.sd) maxLen = std::max(maxLen, c.sd->size());
            }

            for (size_t i = 0; i < maxLen; ++i) {
                f << fileID << ',' << bi;
                for (const Col& c : cols) {
                    f << ',';
                    cell(*c.mv, i);
                    if (c.sd) { f << ','; cell(*c.sd, i); }
                }
                f << '\n';
            }
        }
    }
}  // namespace template_io