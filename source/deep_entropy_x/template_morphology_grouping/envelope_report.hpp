#pragma once
/**
 * @file   envelope_report.hpp
 * @brief  Section 4.7 dynamic-envelope REPORTING. Drives envelopes.hpp over a
 *         subject's kept beats and writes the result as
 *
 *           <stem>_envelopes.csv   one row per (channel, bin, beat, segment)
 *           <stem>_envelopes.bin   the same rows, fixed-width binary
 *
 *         WHY THIS FILE EXISTS SEPARATELY FROM envelopes.hpp. That header is
 *         the measurement: it knows how to score one beat and how to keep a
 *         rolling window. It has no opinion about where beats come from, what
 *         order they are in, or where results go. This file supplies all three
 *         and is the only place that does, so a second consumer of Section 4.7
 *         (the Section 9.5 selector, when it exists) can reuse the tracker
 *         without inheriting a file format.
 *
 *         EVERY INTERMEDIATE IS REPORTED, not just the envelope outputs. The
 *         columns run from the segment bounds used, through the raw per-beat
 *         measures, through the pre-normalisation denominators (total spectral
 *         power, per-band power, per-level wavelet energy), to the rolling
 *         mean/sd/z that 9.5 consumes. A z-score is three divisions away from
 *         the waveform and every one of them can be wrong in a way that still
 *         produces a plausible number; the point of the wide file is that any
 *         reported value can be recomputed by hand from its own row.
 *
 *         THE LAYOUT IS DEFINED TWICE, HERE AND IN read_envelope_bin.py, and
 *         the two must be changed together. writeBin stamps sizeof(record) and
 *         the four compile-time constants into the file header precisely so
 *         that a mismatch is refused by both readers rather than striding
 *         wrongly through the payload and producing plausible numbers.
 *
 *         ------------------------------------------------------------------
 *         THE THREE THINGS THIS FILE DECIDES, all of which the spec leaves open
 *         ------------------------------------------------------------------
 *
 *         1. WHICH BEATS. template_io::BeatsFile::per_channel_beats -- the
 *            KEPT, aligned beats on each bin's shared NaN-padded axis. They are
 *            already co-registered to the bin's r_col, which is what makes one
 *            set of segment bounds valid for every beat in the bin. Beats
 *            dropped upstream (Tukey, unscorable) are absent, so a rolling
 *            window is a window over MEASURED beats and not over wall-clock
 *            time. That is the same convention RollingEnvelope already applies
 *            to NaN values, so the two agree.
 *
 *         2. WHERE THE SEGMENT BOUNDS COME FROM. Auto-detected on the bin's own
 *            averaged template (chN_raw.ecgTemplate, anchored at its r_col) via
 *            the same FeatureMarks finders bin_archive.hpp uses, for the same
 *            reason: at this pipeline stage no operator marks exist. So these
 *            envelopes describe AUTOMATIC segment boundaries. The bounds are
 *            written into every row rather than left implicit, so a later
 *            re-run against operator-marked bounds is comparable rather than
 *            merely different.
 *
 *            P HAS NO OFFSET FINDER, so the P span ends at the Q-onset. That
 *            includes the PR segment, which is flat baseline and drags the P
 *            wave's mean toward it -- stated here because it makes P's
 *            amplitude the trustworthy P measure and P's mean the weak one.
 *
 *         3. ORDER. Beats are pushed in bin order, and within a bin in stored
 *            order, per channel. A Tracker is a sequential object: fed out of
 *            order it returns numbers that are arithmetically correct and
 *            physiologically meaningless. So the bin loop here is deliberately
 *            SERIAL even though bin_archive's equivalent is parallel, and each
 *            channel gets its own Tracker rather than sharing one.
 */

#include "envelopes.hpp"
#include "template_generation\template_io.hpp"
#include "template_marking_gui\feature_marks.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace envelope_report {

    inline constexpr int kNumEcgCh = 3;
    inline const double kNaN = std::numeric_limits<double>::quiet_NaN();

    // Keys into BeatsFile::per_channel_beats, same strings the rest of the
    // pipeline uses (see bin_archive.hpp's ChannelSpec).
    inline constexpr const char* kChannelKeys[kNumEcgCh] = { "CH1", "CH2", "CH3" };

    // ---------------------------------------------------------------------
    // Binary record layout
    // ---------------------------------------------------------------------
    //
    // WRITTEN RAW, so the padding is part of the format. Both structs are
    // pinned by static_assert below: a compiler that lays them out differently
    // fails the build rather than silently emitting a file no reader can parse.
    // Same contract as morphology_csv.hpp's BeatRecord/TemplateRecord.

    struct EnvPointRec {
        double  value, mean, sd, z;
        int32_t n;          // beats actually in the baseline
        uint8_t ready;      // 1 = window filled
        uint8_t pad[3];
    };

    struct EnvelopeRecord {
        uint32_t bin = 0;
        uint32_t beat_in_bin = 0;
        uint32_t beat_seq = 0;      // per channel, across the whole subject
        int32_t  seg_begin = -1;
        int32_t  seg_end = -1;
        int32_t  n_samples = 0;      // non-NaN samples the segment contributed

        uint8_t  channel = 0;      // 0..2 == CH1..CH3
        uint8_t  segment = 0;      // envelopes::Segment
        uint8_t  valid = 0;
        uint8_t  rhythm = 255;    // per_channel_rhythm, 255 = not supplied
        uint8_t  onset_event = 0;      // beat-level, repeated on each segment row
        uint8_t  drift_event = 0;      // beat-level, repeated on each segment row
        uint8_t  pad[2] = { 0, 0 };

        // Raw measures.
        double mean = 0, amplitude = 0, area = 0;
        double band_frac[envelopes::kSpectralBands] = { 0, 0, 0 };
        double centroid = 0;
        double wave_frac[envelopes::kWaveletLevels] = { 0, 0, 0 };

        // Intermediates: the denominators the fractions above divide by.
        double spectral_total_power = 0;
        double band_power[envelopes::kSpectralBands] = { 0, 0, 0 };
        double wave_energy[envelopes::kWaveletLevels] = { 0, 0, 0 };
        double wave_total_energy = 0;

        // Rolling envelopes, short then long, in SegmentEnvelopes order.
        EnvPointRec s_amplitude{}, s_centroid{}, s_wave0{};
        EnvPointRec l_mean{}, l_amplitude{}, l_centroid{}, l_wave0{};
    };

    static_assert(sizeof(EnvPointRec) == 40,
        "EnvPointRec layout changed; bump kVersion and update any reader");
    static_assert(sizeof(EnvelopeRecord) == 456,
        "EnvelopeRecord layout changed; bump kVersion and update any reader");

    inline constexpr char kMagic[9] = "DEXENVL1";
    inline constexpr uint32_t kVersion = 1;

    // ---------------------------------------------------------------------
    // Segment bounds from a bin's averaged template
    // ---------------------------------------------------------------------

    // Auto-detected, on the channel's own template, anchored at its own r_col.
    // Returns all -1 when the template or r_col is unusable, which
    // SegmentSpans::has() reads as "absent" and measure() turns into NaN rather
    // than zero -- the distinction that keeps a ventricular beat's missing P
    // wave out of the P statistics instead of in them as a zero.
    inline envelopes::SegmentSpans spansForChannel(
        const template_io::ChannelMethodTemplate& chRaw, double fs)
    {
        envelopes::SegmentSpans sp;
        const std::vector<double>& ecg = chRaw.ecgTemplate;
        const int rPeak = chRaw.r_col;
        if (ecg.empty() || rPeak < 0 || fs <= 0.0) return sp;

        // ONE DETECTOR. This used to inline the six finder calls, copied from
        // bin_archive.hpp -- which meant it inherited both of that copy's
        // discrepancies against the viewer: no R refinement, and unchained
        // T-onset/T-offset calls. Segment BOUNDS derived from landmarks that
        // disagree with the displayed ones make every envelope in the report
        // incomparable with what an operator sees.
        const FeatureMarks::TemplateLandmarks lm =
            FeatureMarks::detect_template_landmarks(ecg, rPeak, fs);
        if (!lm.valid) return sp;

        // -1 stays -1: SegmentSpans::has() reads it as absent and measure()
        // turns that into NaN rather than zero, which is what keeps a
        // ventricular beat's missing P wave out of the P statistics.
        auto idx = [](double d) {
            return (std::isnan(d) || d < 0.0) ? -1 : static_cast<int>(std::lround(d));
            };

        sp.p_begin = idx(lm.p_begin);
        sp.p_end = idx(lm.q_begin);  // no P-offset finder exists; see header note
        sp.qrs_begin = idx(lm.q_begin);
        sp.qrs_end = idx(lm.s_end);
        sp.t_end = idx(lm.t_end);
        return sp;
    }

    // ---------------------------------------------------------------------
    // Build
    // ---------------------------------------------------------------------

    namespace detail {

        inline EnvPointRec pack(const envelopes::EnvelopePoint& p) {
            EnvPointRec r{};
            r.value = p.value;
            r.mean = p.mean;
            r.sd = p.sd;
            r.z = p.z;
            r.n = p.n;
            r.ready = p.ready ? 1 : 0;
            return r;
        }

        // Rhythm verdict for one kept beat, or 255 when the file did not carry
        // one. Reported, never used: it is the single most useful column for
        // checking whether an onset_event landed on a beat the pipeline already
        // considered premature, and a report that cannot be cross-checked
        // against the rest of the pipeline is a report nobody trusts.
        inline uint8_t rhythmFor(const template_io::BeatsFile& beats,
            const char* key, size_t bin, size_t beat)
        {
            const auto it = beats.per_channel_rhythm.find(key);
            if (it == beats.per_channel_rhythm.end()) return 255;
            if (bin >= it->second.size()) return 255;
            const auto& v = it->second[bin];
            if (beat >= v.size()) return 255;
            return v[beat];
        }

    }  // namespace detail

    // One channel, every bin, in order. Appends to `out`.
    //
    // The Tracker is constructed here and lives for the whole channel, so the
    // rolling window carries ACROSS bin boundaries. That is deliberate: a bin
    // boundary is an artefact of how the recording was chunked, and resetting
    // at each one would blind the long window to drift for its first 30 beats
    // of every bin -- which is most of what the long window exists to see.
    inline void buildChannel(std::vector<EnvelopeRecord>& out,
        const std::vector<template_io::BinTemplates>& bins,
        const template_io::BeatsFile& beats,
        int channel, double fs)
    {
        const char* key = kChannelKeys[channel];
        const auto it = beats.per_channel_beats.find(key);
        if (it == beats.per_channel_beats.end()) return;
        const auto& perBin = it->second;

        envelopes::Tracker tracker;
        uint32_t beatSeq = 0;

        for (size_t b = 0; b < bins.size() && b < perBin.size(); ++b) {
            const template_io::BinTemplates& bt = bins[b];
            if (bt.bad_segment) continue;

            const template_io::ChannelMethodTemplate* chs[kNumEcgCh] = {
                &bt.ch1_raw, &bt.ch2_raw, &bt.ch3_raw };
            const envelopes::SegmentSpans sp = spansForChannel(*chs[channel], fs);

            const auto& binBeats = perBin[b];
            for (size_t k = 0; k < binBeats.size(); ++k) {
                const envelopes::BeatEnvelopes be = tracker.push(binBeats[k], sp);
                const bool onset = envelopes::isOnsetEvent(be);
                const bool drift = envelopes::isDriftEvent(be);

                for (int s = 0; s < 3; ++s) {
                    const envelopes::Segment seg = static_cast<envelopes::Segment>(s);
                    const envelopes::BeatMeasures& m = be.measures[s];
                    const envelopes::SegmentEnvelopes& e = be.env[s];

                    EnvelopeRecord r;
                    r.bin = static_cast<uint32_t>(b);
                    r.beat_in_bin = static_cast<uint32_t>(k);
                    r.beat_seq = beatSeq;
                    r.channel = static_cast<uint8_t>(channel);
                    r.segment = static_cast<uint8_t>(s);
                    r.rhythm = detail::rhythmFor(beats, key, b, k);
                    r.onset_event = onset ? 1 : 0;
                    r.drift_event = drift ? 1 : 0;

                    r.seg_begin = sp.begin(seg);
                    r.seg_end = sp.end(seg);
                    r.n_samples = m.n_samples;
                    r.valid = m.valid ? 1 : 0;

                    r.mean = m.mean;
                    r.amplitude = m.amplitude;
                    r.area = m.area;
                    r.centroid = m.centroid;
                    for (int i = 0; i < envelopes::kSpectralBands; ++i) {
                        r.band_frac[i] = m.band_frac[i];
                        r.band_power[i] = m.band_power[i];
                    }
                    for (int i = 0; i < envelopes::kWaveletLevels; ++i) {
                        r.wave_frac[i] = m.wave_frac[i];
                        r.wave_energy[i] = m.wave_energy[i];
                    }
                    r.spectral_total_power = m.spectral_total_power;
                    r.wave_total_energy = m.wave_total_energy;

                    r.s_amplitude = detail::pack(e.s_amplitude);
                    r.s_centroid = detail::pack(e.s_centroid);
                    r.s_wave0 = detail::pack(e.s_wave0);
                    r.l_mean = detail::pack(e.l_mean);
                    r.l_amplitude = detail::pack(e.l_amplitude);
                    r.l_centroid = detail::pack(e.l_centroid);
                    r.l_wave0 = detail::pack(e.l_wave0);

                    out.push_back(r);
                }
                ++beatSeq;
            }
        }
    }

    // PARALLEL OVER CHANNELS, SERIAL OVER BINS.
    //
    // Each channel has its own Tracker, so the three are genuinely independent
    // and can run at once. Bins cannot: a Tracker is a sequential object whose
    // rolling windows carry across bin boundaries, and feeding it out of order
    // returns numbers that are arithmetically correct and physiologically
    // meaningless. So the parallelism goes exactly one level up from where it
    // would be wrong.
    //
    // Each channel builds into its OWN vector and they are concatenated after,
    // rather than sharing `out` under a lock: the records are appended in
    // per-channel order and a shared push_back would interleave them
    // nondeterministically, making two runs of the same subject produce files
    // that differ only in row order. Reproducibility is worth one copy.
    inline std::vector<EnvelopeRecord> buildEnvelopeReport(
        const std::vector<template_io::BinTemplates>& bins,
        const template_io::BeatsFile& beats, double fs)
    {
        std::vector<std::vector<EnvelopeRecord>> per(kNumEcgCh);
#ifdef _OPENMP
#pragma omp parallel for num_threads(kNumEcgCh) schedule(static)
#endif
        for (int c = 0; c < kNumEcgCh; ++c)
            buildChannel(per[c], bins, beats, c, fs);

        size_t total = 0;
        for (const auto& v : per) total += v.size();
        std::vector<EnvelopeRecord> out;
        out.reserve(total);
        for (auto& v : per)
            out.insert(out.end(), v.begin(), v.end());
        return out;
    }

    // ---------------------------------------------------------------------
    // Writers
    // ---------------------------------------------------------------------

    namespace detail {

        // Buffer, stream, THEN open. pubsetbuf after open does nothing on
        // libstdc++ and is only dependable before any I/O on MSVC, and the
        // buffer has to outlive the flush in ~ofstream, hence the caller
        // declaring it first.
        inline void openBuffered(std::ofstream& f, std::vector<char>& iobuf,
            const std::string& path, std::ios::openmode mode)
        {
            f.rdbuf()->pubsetbuf(iobuf.data(),
                static_cast<std::streamsize>(iobuf.size()));
            f.open(path, mode);
        }

    }  // namespace detail

    // Layout:
    //   [char[8] magic]["DEXENVL1"]
    //   [uint32 version][uint32 record_size]
    //   [int32 short_window][int32 long_window]
    //   [int32 spectral_bands][int32 wavelet_levels]
    //   [double sample_rate_hz]
    //   [uint64 n_records]
    //   n_records x EnvelopeRecord
    //
    // record_size and the four constants are in the header so a reader can
    // refuse a file whose records it would misparse, instead of striding
    // wrongly through it and producing plausible numbers. That is the failure
    // mode the noise-markings .bin had before it grew a magic.
    inline bool writeBin(const std::string& path,
        const std::vector<EnvelopeRecord>& rows, double fs)
    {
        std::vector<char> iobuf(4u << 20);
        std::ofstream f;
        detail::openBuffered(f, iobuf, path, std::ios::binary | std::ios::trunc);
        if (!f) return false;

        const uint32_t recSize = static_cast<uint32_t>(sizeof(EnvelopeRecord));
        const int32_t sw = envelopes::kShortWindow;
        const int32_t lw = envelopes::kLongWindow;
        const int32_t sb = envelopes::kSpectralBands;
        const int32_t wl = envelopes::kWaveletLevels;
        const uint64_t n = rows.size();

        f.write(kMagic, 8);
        f.write(reinterpret_cast<const char*>(&kVersion), 4);
        f.write(reinterpret_cast<const char*>(&recSize), 4);
        f.write(reinterpret_cast<const char*>(&sw), 4);
        f.write(reinterpret_cast<const char*>(&lw), 4);
        f.write(reinterpret_cast<const char*>(&sb), 4);
        f.write(reinterpret_cast<const char*>(&wl), 4);
        f.write(reinterpret_cast<const char*>(&fs), 8);
        f.write(reinterpret_cast<const char*>(&n), 8);

        // One write for the whole block: the records are contiguous and
        // fixed-width, so there is nothing to serialise field by field.
        if (n > 0)
            f.write(reinterpret_cast<const char*>(rows.data()),
                static_cast<std::streamsize>(n * sizeof(EnvelopeRecord)));
        return static_cast<bool>(f);
    }

    // Reader, so the format has exactly one definition rather than one here and
    // another in whatever consumes it. Returns false on a magic, version or
    // record-size mismatch rather than attempting a best-effort parse.
    inline bool readBin(const std::string& path,
        std::vector<EnvelopeRecord>& out, double* fsOut = nullptr)
    {
        std::ifstream f(path, std::ios::binary);
        if (!f) return false;

        char magic[8] = {};
        uint32_t version = 0, recSize = 0;
        int32_t sw = 0, lw = 0, sb = 0, wl = 0;
        double fs = 0.0;
        uint64_t n = 0;

        if (!f.read(magic, 8)) return false;
        for (int i = 0; i < 8; ++i) if (magic[i] != kMagic[i]) return false;
        if (!f.read(reinterpret_cast<char*>(&version), 4)) return false;
        if (!f.read(reinterpret_cast<char*>(&recSize), 4)) return false;
        if (!f.read(reinterpret_cast<char*>(&sw), 4)) return false;
        if (!f.read(reinterpret_cast<char*>(&lw), 4)) return false;
        if (!f.read(reinterpret_cast<char*>(&sb), 4)) return false;
        if (!f.read(reinterpret_cast<char*>(&wl), 4)) return false;
        if (!f.read(reinterpret_cast<char*>(&fs), 8)) return false;
        if (!f.read(reinterpret_cast<char*>(&n), 8)) return false;

        if (version != kVersion) return false;
        if (recSize != sizeof(EnvelopeRecord)) return false;
        if (sb != envelopes::kSpectralBands || wl != envelopes::kWaveletLevels)
            return false;

        out.resize(static_cast<size_t>(n));
        if (n > 0 && !f.read(reinterpret_cast<char*>(out.data()),
            static_cast<std::streamsize>(n * sizeof(EnvelopeRecord)))) {
            out.clear();
            return false;
        }
        if (fsOut) *fsOut = fs;
        return true;
    }

    // ---------------------------------------------------------------------
    // One call for the pipeline
    // ---------------------------------------------------------------------
    /**
     * @brief Build every per-beat envelope record and write both encodings.
     *
     * @param dir        output directory (created by the caller).
     * @param subjectId  file stem; also the subject_id column.
     * @param bins       job.tmpl.bins -- supplies the segment bounds.
     * @param beats      job.beats -- supplies the beats, in order.
     * @param fs         ECG sampling rate, Hz.
     * @return false if the file could not be written.
     */
    inline bool writeEnvelopeReport(const std::string& dir, const std::string& subjectId,
        const std::vector<template_io::BinTemplates>& bins,
        const template_io::BeatsFile& beats, double fs)
    {
        const std::vector<EnvelopeRecord> rows =
            buildEnvelopeReport(bins, beats, fs);

        const std::string base = dir + "/" + subjectId + "_envelopes";
        const bool ok = writeBin(base + ".bin", rows, fs);

        std::fprintf(stderr,
            "  [envelopes] %s: %zu records over %zu bins (%d ch x 3 segments, "
            "short=%d long=%d)\n",
            subjectId.c_str(), rows.size(), bins.size(), kNumEcgCh,
            envelopes::kShortWindow, envelopes::kLongWindow);
        return ok;
    }

}  // namespace envelope_report
