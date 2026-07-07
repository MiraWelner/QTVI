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
            uint8_t bad = b.bad_segment ? 1 : 0;
            f.write(reinterpret_cast<const char*>(&bad), 1);
        }

        writeAveraged(f, data.saecg.ch1_raw);
        writeAveraged(f, data.saecg.ch1_squared);
        writeAveraged(f, data.saecg.ch1_absval);
        writeAveraged(f, data.saecg.ch1_unfiltered);
        writeAveraged(f, data.saecg.ch2_raw);
        writeAveraged(f, data.saecg.ch2_squared);
        writeAveraged(f, data.saecg.ch2_absval);
        writeAveraged(f, data.saecg.ch2_unfiltered);
        writeAveraged(f, data.saecg.ch3_raw);
        writeAveraged(f, data.saecg.ch3_squared);
        writeAveraged(f, data.saecg.ch3_absval);
        writeAveraged(f, data.saecg.ch3_unfiltered);
        writeAveraged(f, data.saecg.ppg);
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
                !readVecD(f, b.ppgTemplate_std))
                throw std::runtime_error("template file truncated mid-bin: " + path);
            uint8_t bad = 0;
            f.read(reinterpret_cast<char*>(&bad), 1);
            b.bad_segment = (bad != 0);
        }

        if (!readAveraged(f, out.saecg.ch1_raw) || !readAveraged(f, out.saecg.ch1_squared) ||
            !readAveraged(f, out.saecg.ch1_absval) || !readAveraged(f, out.saecg.ch1_unfiltered) ||
            !readAveraged(f, out.saecg.ch2_raw) || !readAveraged(f, out.saecg.ch2_squared) ||
            !readAveraged(f, out.saecg.ch2_absval) || !readAveraged(f, out.saecg.ch2_unfiltered) ||
            !readAveraged(f, out.saecg.ch3_raw) || !readAveraged(f, out.saecg.ch3_squared) ||
            !readAveraged(f, out.saecg.ch3_absval) || !readAveraged(f, out.saecg.ch3_unfiltered) ||
            !readAveraged(f, out.saecg.ppg))
            throw std::runtime_error("template file missing SAECG tail: " + path);

        return out;
    }

    void write_saecg_csvfile(const std::string& path, const TemplateFile& data) {
        std::ofstream f(path);
        if (!f) throw std::runtime_error("cannot create: " + path);

        // Column order matches the per-bin method order and the SAECG struct.
        struct Col { const char* name; const AveragedTemplate* avg; };
        const SAECG& s = data.saecg;
        const Col cols[] = {
            {"ch1_raw", &s.ch1_raw}, {"ch1_squared", &s.ch1_squared},
            {"ch1_absval", &s.ch1_absval}, {"ch1_unfiltered", &s.ch1_unfiltered},
            {"ch2_raw", &s.ch2_raw}, {"ch2_squared", &s.ch2_squared},
            {"ch2_absval", &s.ch2_absval}, {"ch2_unfiltered", &s.ch2_unfiltered},
            {"ch3_raw", &s.ch3_raw}, {"ch3_squared", &s.ch3_squared},
            {"ch3_absval", &s.ch3_absval}, {"ch3_unfiltered", &s.ch3_unfiltered},
            {"ppg", &s.ppg},
        };

        // Header.
        f << "sample_index";
        for (const auto& c : cols) f << ',' << c.name;
        f << '\n';

        // How many bins contributed to each average (0 => empty waveform).
        // Written as a labelled row so the count travels with the samples;
        // a strict typed parser can skip this one row.
        f << "n_contributing";
        for (const auto& c : cols) f << ',' << c.avg->n_contributing;
        f << '\n';

        // Longest waveform sets the row count; shorter columns get empty cells.
        size_t maxLen = 0;
        for (const auto& c : cols) maxLen = std::max(maxLen, c.avg->waveform.size());

        f << std::setprecision(9);
        for (size_t i = 0; i < maxLen; ++i) {
            f << i;
            for (const auto& c : cols) {
                f << ',';
                const auto& w = c.avg->waveform;
                if (i < w.size() && !std::isnan(w[i])) f << w[i];   // NaN / out-of-range -> empty
            }
            f << '\n';
        }
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

}  // namespace template_io