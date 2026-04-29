/**
 * @file   3_anneal_segments.cpp
 * @brief  Takes a v2 raw data .bin (512-byte header, 41 channels) and a
 *         noise markings .bin, outputs 1-minute segments with marked noise
 *         removed. The output preserves all 41 input channels (upsampled
 *         + raw (t,v) blocks) sliced to each segment's time window.
 *
 *         Noise exclusion logic:
 *           - ECG noise is only excluded if ALL 3 ECG channels have
 *             overlapping noise marks for a given time region.
 *           - PPG noise is excluded independently.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-22
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstdint>
#include <cstring>
#include "AnnealSegments.hpp"

namespace fs = std::filesystem;
double bin_length = 1.0;

// Side-channel data carried alongside RawData. The algorithm only reads
// RawData; this struct lets the writer emit all 41 channels without the
// algorithm needing to know about them.
static constexpr int NUM_CHANNELS = 41;
struct Extras {
    // NB: parenthesis-init via assignment to invoke the count constructor.
    // Brace-init {NUM_CHANNELS} on a vector<vector<double>> hits the
    // initializer_list ctor, producing a 1-element outer vector instead of
    // a NUM_CHANNELS-element one -- which silently corrupts memory the
    // moment we index past slot 0.
    std::vector<std::vector<double>> upsampled =
        std::vector<std::vector<double>>(NUM_CHANNELS);   // per-slot upsampled samples
    std::vector<std::vector<double>> rawFlat =
        std::vector<std::vector<double>>(NUM_CHANNELS);   // per-slot interleaved (t,v,t,v,...)
    std::vector<float> nativeRates = std::vector<float>(NUM_CHANNELS, 0.0f);
    uint32_t signal_rate = 0, boolean_rate = 0, pacemaker_rate = 0, sleep_rate = 0;
};

// ============================================================================
// Binary I/O
// ============================================================================

/**
 * @brief Read a noise markings .bin produced by step 2 (noise_marking_gui).
 *
 * Format: [uint64 count] then count × [6 doubles per row]
 *   Row: startSample, endSample, startSec, endSec, labelId, typeId
 *   labelId: 0=unknown, 1=PPG, 2=ECG1, 3=ECG2, 4=ECG3
 */
static NoiseMarkings read_noise_bin(const std::string& path) {
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
        default: break;
        }
    }
    return m;
}

/**
 * @brief Read the v2 data .bin produced by file_to_bin (step 1).
 *
 * v2 layout: 512-byte header (4 uint32 rates + 41 upsampled-sizes (uint32) +
 * 41 raw-sizes (uint32) + 41 native-rates (float32) + 1 sleep-size (uint32)),
 * then per slot {upsampled doubles, raw (t,v)-pair doubles}, then sleep doubles.
 *
 * Slot order: 0=Timestamp, 1=ECG1, 2=ECG2, 3=ECG3, 4=PPG, 5..40=other.
 * The algorithm only needs ECG1..PPG + sleep + rates -- those go in `data`.
 * Everything else (all 41 upsampled + raw blocks, native rates, and the
 * file-level rate scalars) goes in `extras` so the writer can re-emit them.
 */
static void read_data_bin(const std::string& path, RawData& data, Extras& extras) {
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("Cannot open bin: " + path);
    f.seekg(0, std::ios::end);
    const uint64_t fileSize = static_cast<uint64_t>(f.tellg());
    f.seekg(0, std::ios::beg);

    constexpr size_t NHF = 4 + 3 * NUM_CHANNELS + 1;   // = 128 fields
    constexpr size_t HDR = NHF * 4;                    // = 512 bytes
    if (fileSize < HDR) throw std::runtime_error("Bin too small: " + path);

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

    // Sequential read: header is followed by 41 x {upsampled, raw} blocks,
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

/**
 * @brief Write the annealed output .bin (all 41 channels preserved).
 *
 * Layout:
 *   Header: [uint64 nSegments][double ppgSR][double ecgSR][double epochSec]
 *           [uint32 nChannels=41][41 x float32 nativeRates]
 *   Per segment: ppg_bin_indexs, ecg_bin_indexs,
 *                ppg, ecg1, ecg2, ecg3, sleep,
 *                then 41 x {upsampled_slice, raw_slice} blocks.
 *
 * The first 7 fields per segment are preserved verbatim from the original
 * schema. The trailing 41 x {upsampled, raw} blocks carry every input
 * channel sliced to this segment's time window:
 *   - Upsampled slice: indices proportional to ecg_bin_indexs. Since every
 *     channel covers the same total duration, ECG-index i in an N-sample
 *     ECG block maps to channel-X-index ceil(i * M / N) in an M-sample
 *     channel-X block. This avoids any per-channel rate bookkeeping.
 *   - Raw slice: filter (t, v) pairs whose t falls inside the segment's
 *     ECG time window. Timestamps stay in absolute seconds-from-recording-
 *     start (same time space as the input).
 */
static void write_output_bin(const std::string& path,
    const std::vector<FinalSegment>& segs,
    const Extras& extras)
{
    std::ofstream out(path, std::ios::binary);
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

    // ECG block size is the reference for proportional slicing. signal_rate
    // gives the time axis for raw filtering.
    const size_t ecgN = extras.upsampled[1].size();   // slot 1 = ECG1
    const double sr = (extras.signal_rate > 0)
        ? static_cast<double>(extras.signal_rate) : 1000.0;

    // Per-channel forward cursors into rwSrc, preserved across segments.
    // Raw pairs are time-sorted and segments are time-ordered, so we never
    // need to revisit pairs we've already passed. Without this, MESA files
    // (8h of data, dense raw blocks, ~500 segments) re-scanned each raw
    // vector from the start for every segment -- quadratic in the number
    // of segments.
    std::vector<size_t> rawCursor(NUM_CHANNELS, 0);

    for (const auto& s : segs) {
        writePairs(s.ppg_bin_indexs);
        writePairs(s.ecg_bin_indexs);
        writeVec(s.ppg);
        writeVec(s.ecg1);
        writeVec(s.ecg2);
        writeVec(s.ecg3);
        writeVec(s.sleep_stages);

        // Per-channel slices. Indices in s.ecg_bin_indexs are 1-based at sr.
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
                // Pre-size from the index pairs so we do one allocation
                // instead of n geometric growths.
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
                    // Bulk insert is one memcpy versus (b - a + 1) push_backs.
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
                // Walk forward through rwSrc using rawCursor[ch], which is
                // preserved across segments. For each window we advance to
                // the first pair with t >= t0, then scan until t > t1, and
                // leave the cursor parked at the next pair to consider for
                // the next window/segment.
                size_t k = rawCursor[ch];
                const size_t N = rwSrc.size();
                for (const auto& p : s.ecg_bin_indexs) {
                    if (p.second < p.first) continue;
                    const double t0 = static_cast<double>(p.first - 1) / sr;
                    const double t1 = static_cast<double>(p.second - 1) / sr;
                    // Advance past pairs strictly before t0.
                    while (k + 1 < N && rwSrc[k] < t0) k += 2;
                    // Collect pairs with t in [t0, t1].
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

// ============================================================================
// Config parsing
// ============================================================================

static std::vector<ProjectConfig> readConfig(const std::string& path) {
    std::vector<ProjectConfig> projects;
    std::ifstream cfg(path);
    if (!cfg.is_open()) {
        std::cerr << "Could not open " << path << std::endl;
        return projects;
    }

    std::string line;
    std::getline(cfg, line); // skip header
    while (std::getline(cfg, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string col;
        std::vector<std::string> cols;
        while (std::getline(ss, col, ',')) cols.push_back(col);

        ProjectConfig pc;
        pc.dataType = cols.size() > 0 ? cols[0] : "unknown";
        pc.binPath = cols.size() > 4 ? cols[4] : "";
        pc.noisePath = cols.size() > 5 ? cols[5] : "";
        pc.annealedPath = cols.size() > 6 ? cols[6] : "";
        projects.push_back(pc);
    }
    return projects;
}

// ============================================================================
// Main
// ============================================================================

int main() {
    auto projects = readConfig("config.csv");
    if (projects.empty()) return 1;

    std::cout << "Select Dataset:\n";
    for (size_t i = 0; i < projects.size(); ++i)
        std::cout << i + 1 << ". " << projects[i].dataType << "\n";

    int choice;
    std::cin >> choice;
    if (choice < 1 || choice >(int)projects.size()) return 1;
    const auto& sel = projects[choice - 1];

    if (!fs::exists(sel.annealedPath))
        fs::create_directories(sel.annealedPath);

    int processed = 0;
    for (const auto& entry : fs::recursive_directory_iterator(sel.binPath)) {
        if (!entry.is_regular_file()) continue;

        std::string ext = entry.path().extension().string();
        std::string name = entry.path().filename().string();
        if (name.find("_noise_markings") != std::string::npos) continue;
        if (ext != ".bin" && ext != ".BIN") continue;

        std::string id = entry.path().stem().string();
        try {
            RawData raw;
            Extras  extras;
            read_data_bin(entry.path().string(), raw, extras);

            NoiseMarkings noise;
            std::string npath = sel.noisePath + "/" + id + "_noise_markings.bin";
            if (fs::exists(npath))
                noise = read_noise_bin(npath);

            auto results = AnnealSegments(raw, noise, bin_length);
            write_output_bin(sel.annealedPath + "/" + id + ".bin", results, extras);

            ++processed;
            std::cerr << "  [" << processed << "] " << id
                << " -> " << results.size() << " bins\n";
        }
        catch (const std::exception& e) {
            std::cerr << "Error processing " << id << ": " << e.what() << "\n";
        }
    }

    std::cerr << "Done. Processed " << processed << " files.\n";
    return 0;
}