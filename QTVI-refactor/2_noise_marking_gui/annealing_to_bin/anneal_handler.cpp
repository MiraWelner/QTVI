/**
 * @file   anneal_handler.cpp
 * @brief  Implementation of step 3. Reads the v2 data .bin, the noise
 *         markings .bin, runs AnnealSegments, and emits the annealed
 *         .bin with all 40 channels preserved.
 *
 *         All the I/O glue lives here; the algorithm itself is in
 *         AnnealSegments.hpp. The file_to_bin step is what produces the
 *         input .bin, so the on-disk format must match between
 *         file_to_bin's writer and read_data_bin() below.
 */

#include "anneal_handler.hpp"
#include "AnnealSegments.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

    // ============================================================================
    // Side-channel data carried alongside RawData.
    //
    // AnnealSegments only reads RawData (PPG/ECG/sleep + their rates).
    // Everything else needed to round-trip the .bin -- header scalars, all
    // 40 upsampled blocks, all 40 raw (t,v) blocks, native rates -- goes in
    // Extras so the writer can re-emit them.
    // ============================================================================

    constexpr int NUM_CHANNELS = 40;

    struct Extras {
        // NB: parenthesis-init via assignment to invoke the count constructor.
        // Brace-init {NUM_CHANNELS} on a vector<vector<double>> hits the
        // initializer_list ctor, producing a 1-element outer vector instead
        // of a NUM_CHANNELS-element one -- which silently corrupts memory the
        // moment we index past slot 0.
        std::vector<std::vector<double>> upsampled =
            std::vector<std::vector<double>>(NUM_CHANNELS);
        std::vector<std::vector<double>> rawFlat =
            std::vector<std::vector<double>>(NUM_CHANNELS);
        std::vector<float> nativeRates = std::vector<float>(NUM_CHANNELS, 0.0f);
        uint32_t signal_rate = 0, boolean_rate = 0, pacemaker_rate = 0, sleep_rate = 0;
    };

    // ============================================================================
    // Noise reader
    //
    // Format: [uint64 count] then count x [6 doubles per row]
    //   Row: startSample, endSample, startSec, endSec, labelId, typeId
    //   labelId: 0=unknown, 1=PPG, 2=ECG1, 3=ECG2, 4=ECG3, 5=ABP
    // (typeId is unused here -- we only care about which channel was marked.)
    // ============================================================================

    NoiseMarkings read_noise_bin(const std::filesystem::path& path) {
        NoiseMarkings m;
        std::ifstream f(path, std::ios::binary);
        if (!f.is_open()) return m;

        uint64_t count = 0;
        f.read(reinterpret_cast<char*>(&count), 8);

        for (uint64_t i = 0; i < count; ++i) {
            double row[6];
            f.read(reinterpret_cast<char*>(row), 48);
            std::pair<double, double> iv = { row[2], row[3] };
            switch (static_cast<int>(row[4])) {
            case 1: m.ppg.push_back(iv);  break;
            case 2: m.ecg1.push_back(iv); break;
            case 3: m.ecg2.push_back(iv); break;
            case 4: m.ecg3.push_back(iv); break;
                // labelId 5 (ABP) is silently dropped -- AnnealSegments only
                // operates on PPG/ECG. If you start excluding ABP in the future,
                // add a branch here and a field to NoiseMarkings.
            default: break;
            }
        }
        return m;
    }

    // ============================================================================
    // Data reader
    //
    // Reads the v2 data .bin produced by file_to_bin. The header layout MUST
    // match file_to_bin's writer: 4 uint32 scalars + 40 upsampled-sizes + 40
    // raw-sizes + 40 native-rate floats + 1 sleep-size = 125 fields = 500
    // bytes. (NUM_HEADER_FIELDS in file_to_bin.hpp.)
    // ============================================================================

    void read_data_bin(const std::filesystem::path& path,
        RawData& data, Extras& extras)
    {
        std::ifstream f(path, std::ios::binary);
        if (!f.is_open()) throw std::runtime_error("cannot open: " + path.string());

        f.seekg(0, std::ios::end);
        const uint64_t fileSize = static_cast<uint64_t>(f.tellg());
        f.seekg(0, std::ios::beg);

        constexpr size_t NHF = 4 + 3 * NUM_CHANNELS + 1;   // = 125 fields
        constexpr size_t HDR = NHF * 4;                    // = 500 bytes
        if (fileSize < HDR)
            throw std::runtime_error("bin too small: " + path.string());

        uint32_t hdr[NHF] = {};
        f.read(reinterpret_cast<char*>(hdr), HDR);

        extras.signal_rate = hdr[0];
        extras.boolean_rate = hdr[1];
        extras.pacemaker_rate = hdr[2];
        extras.sleep_rate = hdr[3];

        std::vector<uint32_t> sizes_up(NUM_CHANNELS), sizes_raw(NUM_CHANNELS);
        for (int i = 0; i < NUM_CHANNELS; ++i) {
            sizes_up[i] = hdr[4 + i];
            sizes_raw[i] = hdr[4 + NUM_CHANNELS + i];
            std::memcpy(&extras.nativeRates[i], &hdr[4 + 2 * NUM_CHANNELS + i], 4);
        }
        const uint32_t sleep_count = hdr[4 + 3 * NUM_CHANNELS];

        data.ecgSR = static_cast<double>(extras.signal_rate);
        data.ppgSR = static_cast<double>(extras.signal_rate);
        data.scoringEpochSec = static_cast<double>(extras.sleep_rate);

        // Sequential read: header is followed by 40 x {upsampled, raw} blocks,
        // then sleep stages. safeReadDoubles clamps to whatever's left in the
        // file so a truncated bin doesn't throw.
        auto safeReadDoubles = [&](std::vector<double>& dest, uint64_t count) {
            const uint64_t pos = static_cast<uint64_t>(f.tellg());
            const uint64_t avail = (fileSize > pos) ? (fileSize - pos) / 8 : 0;
            const uint64_t actual = std::min<uint64_t>(count, avail);
            dest.resize(actual);
            if (actual > 0)
                f.read(reinterpret_cast<char*>(dest.data()),
                    static_cast<std::streamsize>(actual * 8));
            };
        for (int i = 0; i < NUM_CHANNELS; ++i) {
            safeReadDoubles(extras.upsampled[i], sizes_up[i]);
            safeReadDoubles(extras.rawFlat[i], static_cast<uint64_t>(sizes_raw[i]) * 2);
        }
        safeReadDoubles(data.sleepStages, sleep_count);

        // Mirror algorithm-facing channels into RawData (slots 1..4).
        data.ecg1 = extras.upsampled[1];
        data.ecg2 = extras.upsampled[2];
        data.ecg3 = extras.upsampled[3];
        data.ppg = extras.upsampled[4];
    }

    // ============================================================================
    // Output writer
    //
    // Layout:
    //   Header: [uint64 nSegments][double ppgSR][double ecgSR][double epochSec]
    //           [uint32 nChannels=40][40 x float32 nativeRates]
    //   Per segment:
    //     ppg_bin_indexs, ecg_bin_indexs, ppg, ecg1, ecg2, ecg3, sleep,
    //     then 40 x {upsampled_slice, raw_slice}.
    //
    // The 40 trailing per-channel slices preserve every input channel sliced
    // to the segment's time window:
    //   - Upsampled slice: indices proportional to ecg_bin_indexs. ECG-index
    //     i in an N-sample ECG block maps to channel-X-index ceil(i*M/N) in
    //     an M-sample block of channel X. Avoids per-channel rate bookkeeping.
    //   - Raw slice: filter (t,v) pairs whose t falls inside the ECG time
    //     window. Timestamps stay in absolute seconds-from-recording-start.
    // ============================================================================

    void write_output_bin(const std::filesystem::path& path,
        const std::vector<FinalSegment>& segs,
        const Extras& extras)
    {
        std::ofstream out(path, std::ios::binary);
        if (!out.is_open())
            throw std::runtime_error("cannot create: " + path.string());

        uint64_t n = segs.size();
        out.write(reinterpret_cast<char*>(&n), 8);

        if (n > 0) {
            out.write(reinterpret_cast<const char*>(&segs[0].ppgSampleRate), 8);
            out.write(reinterpret_cast<const char*>(&segs[0].ecgSampleRate), 8);
            out.write(reinterpret_cast<const char*>(&segs[0].scoring_epoch_size_sec), 8);
        }
        const uint32_t nch = NUM_CHANNELS;
        out.write(reinterpret_cast<const char*>(&nch), 4);
        out.write(reinterpret_cast<const char*>(extras.nativeRates.data()),
            NUM_CHANNELS * sizeof(float));

        auto writePairs = [&](const std::vector<std::pair<uint64_t, uint64_t>>& v) {
            uint64_t sz = v.size();
            out.write(reinterpret_cast<char*>(&sz), 8);
            for (const auto& p : v) {
                out.write(reinterpret_cast<const char*>(&p.first), 8);
                out.write(reinterpret_cast<const char*>(&p.second), 8);
            }
            };
        auto writeVec = [&](const std::vector<double>& v) {
            uint64_t sz = v.size();
            out.write(reinterpret_cast<char*>(&sz), 8);
            if (!v.empty())
                out.write(reinterpret_cast<const char*>(v.data()), sz * 8);
            };

        const size_t ecgN = extras.upsampled[1].size();   // slot 1 = ECG1
        const double sr = (extras.signal_rate > 0)
            ? static_cast<double>(extras.signal_rate) : 1000.0;

        // Per-channel forward cursors into rwSrc, preserved across segments.
        // Raw pairs are time-sorted and segments are time-ordered, so we never
        // revisit pairs we've already passed. Without this, MESA files (~500
        // segments, dense raw blocks) re-scanned each raw vector from the
        // start for every segment -- quadratic in the number of segments.
        std::vector<size_t> rawCursor(NUM_CHANNELS, 0);

        for (const auto& s : segs) {
            writePairs(s.ppg_bin_indexs);
            writePairs(s.ecg_bin_indexs);
            writeVec(s.ppg);
            writeVec(s.ecg1);
            writeVec(s.ecg2);
            writeVec(s.ecg3);
            writeVec(s.sleep_stages);

            for (int ch = 0; ch < NUM_CHANNELS; ++ch) {
                const auto& upSrc = extras.upsampled[ch];
                const auto& rwSrc = extras.rawFlat[ch];

                // ---- Upsampled slice via proportional indexing ----
                // Sentinel passthrough: file_to_bin writes a single -1.0 for
                // missing channels. Pass it through so the schema is uniform.
                std::vector<double> upOut;
                if (upSrc.size() == 1 && upSrc[0] == -1.0) {
                    upOut = upSrc;
                }
                else if (!upSrc.empty() && ecgN > 0) {
                    const double scale = static_cast<double>(upSrc.size()) /
                        static_cast<double>(ecgN);
                    size_t totalSamples = 0;
                    for (const auto& p : s.ecg_bin_indexs) {
                        if (p.second < p.first) continue;
                        uint64_t a = static_cast<uint64_t>((p.first - 1) * scale) + 1;
                        uint64_t b = static_cast<uint64_t>((p.second - 1) * scale) + 1;
                        if (a > upSrc.size()) continue;
                        if (b > upSrc.size()) b = upSrc.size();
                        totalSamples += static_cast<size_t>(b - a + 1);
                    }
                    upOut.reserve(totalSamples);
                    for (const auto& p : s.ecg_bin_indexs) {
                        if (p.second < p.first) continue;
                        uint64_t a = static_cast<uint64_t>((p.first - 1) * scale) + 1;
                        uint64_t b = static_cast<uint64_t>((p.second - 1) * scale) + 1;
                        if (a > upSrc.size()) continue;
                        if (b > upSrc.size()) b = upSrc.size();
                        upOut.insert(upOut.end(),
                            upSrc.begin() + (a - 1),
                            upSrc.begin() + b);
                    }
                }
                writeVec(upOut);

                // ---- Raw (t, v) slice via time-window filter ----
                // Sentinel passthrough: missing raw channel is a single (-1, -1).
                std::vector<double> rwOut;
                const bool rawIsSentinel = (rwSrc.size() == 2 &&
                    rwSrc[0] == -1.0 && rwSrc[1] == -1.0);
                if (rawIsSentinel) {
                    rwOut = { -1.0, -1.0 };
                }
                else if (!rwSrc.empty()) {
                    size_t k = rawCursor[ch];
                    const size_t N = rwSrc.size();
                    for (const auto& p : s.ecg_bin_indexs) {
                        if (p.second < p.first) continue;
                        const double t0 = static_cast<double>(p.first - 1) / sr;
                        const double t1 = static_cast<double>(p.second - 1) / sr;
                        while (k + 1 < N && rwSrc[k] < t0) k += 2;
                        while (k + 1 < N && rwSrc[k] <= t1) {
                            rwOut.push_back(rwSrc[k]);
                            rwOut.push_back(rwSrc[k + 1]);
                            k += 2;
                        }
                    }
                    rawCursor[ch] = k;
                }
                const uint64_t nPairs = rwOut.size() / 2;
                out.write(reinterpret_cast<const char*>(&nPairs), 8);
                if (!rwOut.empty())
                    out.write(reinterpret_cast<const char*>(rwOut.data()),
                        rwOut.size() * 8);
            }
        }
    }

}   // anonymous namespace

// ============================================================================
// Public entry points
// ============================================================================

bool annealOneFile(const std::filesystem::path& binPath,
    const std::filesystem::path& noisePath,
    const std::filesystem::path& outPath,
    double binLengthMin)
{
    try {
        RawData raw;
        Extras  extras;
        read_data_bin(binPath, raw, extras);

        NoiseMarkings noise;
        if (std::filesystem::exists(noisePath))
            noise = read_noise_bin(noisePath);
        else
            std::cerr << "  no noise file at " << noisePath
            << " -- annealing with no exclusions\n";

        auto results = AnnealSegments(raw, noise, binLengthMin);
        write_output_bin(outPath, results, extras);

        std::cerr << "  -> " << results.size() << " bins -> "
            << outPath.filename() << "\n";
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "  ERROR: " << e.what() << "\n";
        return false;
    }
}

