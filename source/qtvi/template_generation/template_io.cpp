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
            writeVecD(f, m.ecg_template_iqr);
            f.write(reinterpret_cast<const char*>(&m.alignment_point), 8);
            f.write(reinterpret_cast<const char*>(&m.r_col), 4);
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
            if (!readVecD(f, m.ecg_template_iqr)) return false;
            if (!f.read(reinterpret_cast<char*>(&m.alignment_point), 8)) return false;
            if (!f.read(reinterpret_cast<char*>(&m.r_col), 4)) return false;
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

        return out;
    }

    // Long/tidy snips CSV: for each channel, two columns -- the snip id
    // (0,1,2,... across all retained beats of that channel, repeated once per
    // sample) and the sample value. Bad-segment bins are skipped; NaN tails
    // are dropped. Channels with fewer total samples are blank-padded.
    void write_snips_csv(const std::string& path, const BeatsFile& beats) {
        std::ofstream f(path);
        if (!f) throw std::runtime_error("cannot create: " + path);
        f << std::setprecision(4);

        static const char* ORDER[] = { "CH1","CH2","CH3","PPG","ABP","ART","ART_PULM" };

        // One column per snip: header "{CHAN}_{n}", rows = sample index.
        struct Col { std::string name; const std::vector<double>* v; };
        std::vector<Col> cols;
        size_t maxLen = 0;
        for (const char* ch : ORDER) {
            auto it = beats.per_channel_beats.find(ch);
            if (it == beats.per_channel_beats.end()) continue;
            int snip = 0;
            for (size_t b = 0; b < it->second.size(); ++b) {
                if (b < beats.bad_segment.size() && beats.bad_segment[b]) continue;
                for (const auto& beat : it->second[b]) {
                    cols.push_back({ std::string(ch) + "_" + std::to_string(snip), &beat });
                    if (beat.size() > maxLen) maxLen = beat.size();
                    ++snip;
                }
            }
        }

        for (size_t c = 0; c < cols.size(); ++c) { if (c) f << ','; f << cols[c].name; }
        f << '\n';

        for (size_t s = 0; s < maxLen; ++s) {
            for (size_t c = 0; c < cols.size(); ++c) {
                if (c) f << ',';
                const std::vector<double>& v = *cols[c].v;
                if (s < v.size() && !std::isnan(v[s])) f << v[s];
            }
            f << '\n';
        }
    }
}  // namespace template_io