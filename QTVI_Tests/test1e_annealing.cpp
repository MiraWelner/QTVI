// ============================================================================
// test1c_anneal.cpp -- QTVI_Tests (Google Test)
//
// ACCEPTANCE TEST: the anneal step (anneal_handler.cpp) against its spec.
//   - Bin the recording at bin_size_minutes intervals.
//   - Remove samples only where R-peak noise is marked on every present
//     markable channel (intersection rule). Accelerometer marks never drive
//     removal.
//   - Fragments >= half a bin are kept. Fragments < half a bin merge into the
//     adjacent bin on the side where their mass lies. Edge fragments with no
//     neighbour are discarded.
//   - Output: annealed bins in the 36-channel slot layout, with sentinels for
//     absent channels.
//   - Parity: bin stride = bin_size + 1 samples; banker's rounding.
// Walkthrough cases (spec 3.3): clean full bin, noise-split bin, edge discard,
// channel intersection.
//
// HOW. Everything inside the anneal is private to anneal_handler.cpp, so the
// test goes through the real entry point, anneal_one_file: it writes a small
// recording and its noise marks, anneals them, and reads the result back with
// read_input_binfile. Each ECG sample's VALUE is its own 0-based sample number
// (at 200 Hz), so every annealed bin says exactly which original samples it
// kept -- no tolerance, no guessing.
//
// No fixture files: every recording is made here, in a temp folder.
//
// PROJECT SETTINGS: as test1b, plus anneal_handler.cpp in the test project
// (Add -> Existing Item, Add As Link). user_annotation_handler.cpp (added for
// test1b) writes the noise marks.
// ============================================================================
#include "pch.h"
#include "prep_for_peakfinding/anneal_handler.hpp"                  // anneal_one_file
#include "peak_finding/peakfinding_io.hpp"               // read_input_binfile
#include "noise_marking_gui/user_annotation_handler.h"   // annotation_handler (marks writer)

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {
    constexpr double kEcgHz = 200.0, kPpgHz = 100.0, kBinMin = 1.0;
    constexpr long kBinSize = static_cast<long>(kEcgHz * 60 * kBinMin);   // 12000 samples

    // ---- a recording -------------------------------------------------------
    // file_to_bin v1 layout: ECG1 and PPG present (each sample = its index),
    // the timestamp channel on the ECG grid, the other 33 slots absent.
    void writeDataBin(const fs::path& p, double totalSec) {
        constexpr int N = 36;
        const uint32_t nE = static_cast<uint32_t>(std::llround(totalSec * kEcgHz));
        const uint32_t nP = static_cast<uint32_t>(std::llround(totalSec * kPpgHz));
        uint32_t up[N], raw[N]; float nat[N], upr[N];
        for (int i = 0; i < N; ++i) { up[i] = 1; raw[i] = 1; nat[i] = 0; upr[i] = 0; }
        for (int c : { 0, 1 }) { up[c] = nE; raw[c] = nE; nat[c] = upr[c] = (float)kEcgHz; }
        up[4] = nP; raw[4] = nP; nat[4] = upr[4] = (float)kPpgHz;
        std::ofstream f(p, std::ios::binary);
        auto w32 = [&](uint32_t v) { f.write((const char*)&v, 4); };
        w32(1); w32(N); w32(30);
        f.write((const char*)up, 4 * N); f.write((const char*)raw, 4 * N);
        f.write((const char*)nat, 4 * N); f.write((const char*)upr, 4 * N);
        w32(1);
        const double t0 = 1.6e12;   // a real epoch-ms start
        for (int c = 0; c < N; ++c) {
            if (upr[c] == 0.0f) {   // absent: the sentinels file_to_bin writes
                const double m = -1, mm[2] = { -1, -1 };
                f.write((const char*)&m, 8); f.write((const char*)mm, 16);
                continue;
            }
            const double hz = upr[c];
            for (uint32_t i = 0; i < up[c]; ++i) {
                const double v = (c == 0) ? t0 + i * 1000.0 / hz : (double)i;
                f.write((const char*)&v, 8);
            }
            for (uint32_t i = 0; i < raw[c]; ++i) {
                const double t = t0 + i * 1000.0 / hz, pr[2] = { t, (c == 0) ? t : (double)i };
                f.write((const char*)pr, 16);
            }
        }
        const double s = -1; f.write((const char*)&s, 8);   // no sleep staging
    }

    // ---- noise marks, as the GUI writes them on Finish ----------------------
    struct NoiseMark { const char* channel; const char* type; double startSec, endSec, rateHz; };
    void writeNoiseBin(const fs::path& p, const std::vector<NoiseMark>& marks) {
        annotation_handler h;
        for (const NoiseMark& m : marks)
            h.addSegment((int)std::llround(m.startSec * m.rateHz), (int)std::llround(m.endSec * m.rateHz),
                m.channel, m.type, m.rateHz);
        h.export_marking_binfile(p.string());
    }
    NoiseMark rNoise(const char* ch, double a, double b) {
        return { ch, "1) R Peak Noise", a, b, std::string(ch) == "PPG" ? kPpgHz : kEcgHz };
    }

    // ---- run the real anneal; summarise what it kept -------------------------
    struct Result {
        bool ok = false;
        AnnealedData data;
        std::vector<long> binFirst, binLast, binSize;   // per bin, from the sample values
        std::map<long, int> uses;                       // sample -> how many bins hold it
        int holds(long s) const { auto it = uses.find(s); return it == uses.end() ? 0 : it->second; }
        bool keptAll(long a, long b) const { for (long s = a; s <= b; ++s) if (!holds(s)) return false; return true; }
        bool keptNone(long a, long b) const { for (long s = a; s <= b; ++s) if (holds(s)) return false; return true; }
        long duplicated() const { long n = 0; for (auto& [s, k] : uses) if (k > 1) ++n; return n; }
    };
    Result anneal(const std::string& name, double totalSec, const std::vector<NoiseMark>& marks) {
        Result r;
        const fs::path dir = fs::temp_directory_path() / "qtvi_test1c" / name;
        fs::create_directories(dir);
        const fs::path data = dir / "rec.bin", noise = dir / "rec_noise_markings.bin", out = dir / "rec_annealed.bin";
        fs::remove(noise); fs::remove(out);
        writeDataBin(data, totalSec);
        if (!marks.empty()) writeNoiseBin(noise, marks);
        if (!anneal_one_file(data, noise, out, kBinMin, 0.0, false, false, false)) return r;
        r.data = read_input_binfile(out.string());
        r.ok = true;
        for (const auto& b : r.data.bins) {
            if (b.ecg_signal_1.empty()) { r.binFirst.push_back(-1); r.binLast.push_back(-1); r.binSize.push_back(0); continue; }
            r.binFirst.push_back(std::lround(b.ecg_signal_1.front()));
            r.binLast.push_back(std::lround(b.ecg_signal_1.back()));
            r.binSize.push_back((long)b.ecg_signal_1.size());
            for (double v : b.ecg_signal_1) ++r.uses[std::lround(v)];
        }
        return r;
    }
    long S(double sec) { return std::lround(sec * kEcgHz); }   // seconds -> 0-based ECG sample
}

// ============================================================================
// WALKTHROUGH CASES (spec 3.3)
// ============================================================================

TEST(test1c_anneal, CleanFullBins) {
    // 3 min, 1-min bins, no marks: three bins, the whole record, nothing twice.
    const Result r = anneal("clean", 180, {});
    ASSERT_TRUE(r.ok);
    EXPECT_EQ(r.data.bins.size(), 3u);
    EXPECT_TRUE(r.keptAll(0, S(180) - 1)) << "a sample of a clean record is missing";
    EXPECT_EQ(r.duplicated(), 0) << "samples sit in two bins at once";
}

TEST(test1c_anneal, NoiseSplitSmallFragmentsMergeTowardTheirMass) {
    // R noise on ECG1 + PPG at 80-100 s, inside bin 2 (60-120 s): 20 s either
    // side, both under half a bin. The left one merges into bin 1, the right
    // one into bin 3; 80-100 s is gone.
    const Result r = anneal("split_small", 180, { rNoise("ECG1", 80, 100), rNoise("PPG", 80, 100) });
    ASSERT_TRUE(r.ok);
    ASSERT_EQ(r.data.bins.size(), 2u);
    EXPECT_TRUE(r.keptNone(S(80) + 1, S(100) - 1)) << "the marked noise was not removed";
    EXPECT_TRUE(r.keptAll(0, S(80))) << "bin 1 + the left fragment are not intact";
    EXPECT_TRUE(r.keptAll(S(100), S(180) - 1)) << "the right fragment + bin 3 are not intact";
    EXPECT_EQ(r.binFirst[0], 0);
    EXPECT_EQ(r.binLast[0], S(80)) << "the left fragment should end bin 1";
    EXPECT_EQ(r.binFirst[1], S(100)) << "the right fragment should start the next bin";
}

TEST(test1c_anneal, NoiseSplitBigFragmentKeptAsItsOwnBin) {
    // R noise at 60-70 s: bin 2 keeps 70-120 s (50 s >= half a bin) as a bin.
    const Result r = anneal("split_big", 180, { rNoise("ECG1", 60, 70), rNoise("PPG", 60, 70) });
    ASSERT_TRUE(r.ok);
    ASSERT_EQ(r.data.bins.size(), 3u);
    EXPECT_TRUE(r.keptNone(S(60) + 2, S(70) - 1)) << "the marked noise was not removed";
    EXPECT_EQ(r.binFirst[1], S(70)) << "the kept fragment should be a bin of its own";
    EXPECT_TRUE(r.keptAll(S(70), S(120))) << "the kept fragment is not intact";
}

TEST(test1c_anneal, EdgeFragmentAtStartDiscardedRestMergesRight) {
    // R noise at 10-45 s in bin 1: the 10 s before it has no left neighbour
    // (discarded); the 15 s after it merges right into bin 2.
    const Result r = anneal("edge_start", 180, { rNoise("ECG1", 10, 45), rNoise("PPG", 10, 45) });
    ASSERT_TRUE(r.ok);
    EXPECT_TRUE(r.keptNone(0, S(45) - 1)) << "the edge fragment before the noise was not discarded";
    EXPECT_TRUE(r.keptAll(S(45), S(120))) << "the fragment after the noise did not merge into bin 2";
    EXPECT_EQ(r.binFirst[0], S(45));
}

TEST(test1c_anneal, EdgeFragmentAtEndDiscardedRestMergesLeft) {
    // R noise at 135-170 s in the last bin: the 15 s before it (120-135 s) is
    // under half a bin with its mass on the left, so it merges into bin 2; the
    // 10 s after it has no right neighbour (discarded).
    const Result r = anneal("edge_end", 180, { rNoise("ECG1", 135, 170), rNoise("PPG", 135, 170) });
    ASSERT_TRUE(r.ok);
    EXPECT_TRUE(r.keptNone(S(170) + 1, S(180) - 1)) << "the edge fragment after the noise was not discarded";
    EXPECT_TRUE(r.keptAll(S(120), S(135))) << "the 120-135 s fragment was lost instead of merging into bin 2";
}

TEST(test1c_anneal, ChannelIntersectionRemovesOnlyTheOverlap) {
    // ECG1 noisy 80-100 s, PPG noisy 90-110 s: only 90-100 s is noisy on both.
    const Result r = anneal("intersection", 180, { rNoise("ECG1", 80, 100), rNoise("PPG", 90, 110) });
    ASSERT_TRUE(r.ok);
    EXPECT_TRUE(r.keptAll(S(80), S(90))) << "80-90 s is noisy on ECG1 only, and must stay";
    EXPECT_TRUE(r.keptNone(S(90) + 1, S(100) - 1)) << "90-100 s is noisy on both, and must go";
    EXPECT_TRUE(r.keptAll(S(100), S(110))) << "100-110 s is noisy on PPG only, and must stay";
}

TEST(test1c_anneal, UnmarkedPresentChannelVetoesRemoval) {
    // ECG1 noisy 80-100 s; PPG is IN the recording but has no noise marks. The
    // spec's "every present markable channel" includes it, so nothing is noisy
    // on every channel and nothing may be removed.
    const Result r = anneal("unmarked_present", 180, { rNoise("ECG1", 80, 100) });
    ASSERT_TRUE(r.ok);
    EXPECT_TRUE(r.keptAll(S(80), S(100))) << "removed although PPG, present and clean, was not marked";
}

TEST(test1c_anneal, MarksOnAChannelTheRecordingLacksAreIgnored) {
    // The recording has ECG1 and PPG, not ECG2. ECG2's marks neither block a
    // removal (80-100 s is marked on both present channels) nor cause one
    // (0-30 s is marked on ECG2 alone).
    const Result r = anneal("absent_channel", 180, { rNoise("ECG1", 80, 100), rNoise("PPG", 80, 100),
        rNoise("ECG2", 0, 30), rNoise("ECG2", 70, 110) });
    ASSERT_TRUE(r.ok);
    EXPECT_TRUE(r.keptNone(S(80) + 1, S(100) - 1)) << "an absent channel blocked the removal";
    EXPECT_TRUE(r.keptAll(0, S(30))) << "marks on an absent channel removed data";
}

TEST(test1c_anneal, AccelerometerMarksNeverDriveRemoval) {
    const Result r = anneal("accel", 180, { rNoise("ACCEL", 0, 180) });
    ASSERT_TRUE(r.ok);
    EXPECT_EQ(r.data.bins.size(), 3u);
    EXPECT_TRUE(r.keptAll(0, S(180) - 1));
}

TEST(test1c_anneal, OnlyRPeakNoiseDrivesRemoval) {
    const NoiseMark m1{ "ECG1", "2) Minor Noise", 80, 100, kEcgHz }, m2{ "PPG", "2) Minor Noise", 80, 100, kPpgHz };
    const Result r = anneal("minor", 180, { m1, m2 });
    ASSERT_TRUE(r.ok);
    EXPECT_TRUE(r.keptAll(0, S(180) - 1)) << "Minor Noise removed samples";
}

// ============================================================================
// THE RECORDING'S TAIL (a fragment like any other)
// ============================================================================

TEST(test1c_anneal, TailOfHalfABinOrMoreKeptOnce) {
    // 3 min 40 s: the 40 s tail is >= half a bin -> kept, once.
    const Result r = anneal("tail_long", 220, {});
    ASSERT_TRUE(r.ok);
    EXPECT_TRUE(r.keptAll(0, S(220) - 1)) << "part of the record is missing";
    EXPECT_EQ(r.duplicated(), 0) << "the tail bin repeats samples of the bin before it";
}

TEST(test1c_anneal, TailUnderHalfABinLeavesNoGap) {
    // 3 min 20 s: the 20 s tail is < half a bin -> merged left into bin 3 (or,
    // if that is the rule, discarded) -- but nothing BEFORE it may be lost.
    const Result r = anneal("tail_short", 200, {});
    ASSERT_TRUE(r.ok);
    EXPECT_TRUE(r.keptAll(0, S(180) - 1)) << "samples before the tail were dropped";
}

// ============================================================================
// IMPLEMENTATION PARITY
// ============================================================================

TEST(test1c_anneal, BinStrideIsBinSizePlusOne) {
    const Result r = anneal("stride", 180, {});
    ASSERT_TRUE(r.ok);
    ASSERT_GE(r.data.bins.size(), 2u);
    EXPECT_EQ(r.binFirst[0], 0);
    EXPECT_EQ(r.binSize[0], kBinSize + 1);
    EXPECT_EQ(r.binFirst[1], kBinSize + 1);
}

TEST(test1c_anneal, BankersRoundingAtHalfSampleBoundaries) {
    // ECG1 marks at 400 Hz: 80.0025 s and 100.0025 s are 16000.5 and 20000.5
    // samples at the anneal's 200 Hz. Half to even gives 16000 and 20000;
    // half away from zero would give 16001 and 20001.
    const NoiseMark e{ "ECG1", "1) R Peak Noise", 80.0025, 100.0025, 400.0 };
    const Result r = anneal("bankers", 180, { e, rNoise("PPG", 79, 101) });
    ASSERT_TRUE(r.ok);
    ASSERT_EQ(r.data.bins.size(), 2u);
    EXPECT_EQ(r.binLast[0], 16000);
    EXPECT_EQ(r.binFirst[1], 20000);
}

// ============================================================================
// OUTPUT LAYOUT
// ============================================================================

TEST(test1c_anneal, OutputHas36SlotsWithSentinelsForAbsentChannels) {
    const Result r = anneal("slots", 180, {});
    ASSERT_TRUE(r.ok);
    ASSERT_FALSE(r.data.bins.empty());
    for (size_t b = 0; b < r.data.bins.size(); ++b) {
        SCOPED_TRACE("bin " + std::to_string(b));
        const auto& up = r.data.bins[b].all_upsampled;
        ASSERT_EQ(up.size(), 36u);
        for (int c = 0; c < 36; ++c) {
            SCOPED_TRACE("slot " + std::to_string(c));
            const bool present = (c == 0 || c == 1 || c == 4);   // timestamp, ECG1, PPG
            const bool sentinel = (up[c].size() == 1 && up[c][0] == -1.0);
            EXPECT_EQ(sentinel, !present);
        }
    }
}