// test2f_cgm_bin_creation.cpp -- acceptance tests for the CGM pipeline.
//
// One test per statement of the spec's acceptance test:
//   "Warm-up readings are dropped, single gaps interpolated, TIR plus TAR plus
//    TBR equals 100 percent, and short/sparse bins are correctly marked
//    invalid."
//
// DATA. data/cgm/, next to this .cpp file, holds:
//   sample_cgm_file.csv                 the test input
//   TEST_sample_cgm_file_baseline.bin   the expected output, made from it
//   TEST_sample_cgm_file.bin            written by each run
// The test compares what it writes with the baseline.
//
// sample_cgm_file.csv is ECG201V3_glucose_9-22-2026.csv with these edits, so
// that every statement has something to test:
//   removed 09-11-2026 03:51 AM (192)          -> single gap at slot 9
//   removed 09-11-2026 06:51 AM, 07:06 AM      -> 2-reading gap, slots 21-22
//   changed 09-11-2026 02:51 PM 175 -> 65      -> slot 53, a reading below 70
//   removed the last three readings            -> last bin holds one reading
// Slot 0 is the first reading after warm-up (09-11-2026 01:36 AM). Bins are
// 60 min (4 slots), so slot k is in bin k / 4.
//
// data/cgm is found next to THIS .cpp file (via __FILE__), so the test works
// whatever folder the test runner starts in. To keep the data elsewhere, define
// CGM_DATA_DIR as an absolute path in the project's preprocessor definitions.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <string>

#include "cgm/cgm_io.hpp"
#include "cgm/cgm_pipeline.hpp"
#include "cgm/cgm_selection.hpp"   // needs the Qt include paths (QFileDialog)

using namespace cgm_pipeline;
namespace fs = std::filesystem;

namespace {

#ifdef CGM_DATA_DIR
    const fs::path input_output_testdata_dir = CGM_DATA_DIR;
#else
    const fs::path input_output_testdata_dir = fs::path(__FILE__).parent_path() / "data" / "cgm";
#endif
    const fs::path inputfile_path = input_output_testdata_dir / "sample_cgm_file.csv";
    const fs::path baseline_path = input_output_testdata_dir / "TEST_sample_cgm_file.bin";
    const std::string bittiumfile_samplename = "TEST";

    // Process the CSV once, write the .bin, read it back. Shared by all tests.
    // (Not named Run: inside a TEST body, testing::Test::Run would hide it.)
    struct CgmRun {
        bool ok = false;
        std::vector<CgmReading> csv;          // readings as parsed from the CSV
        cgm_selection::CgmOutcome outcome;    // pipeline result + written file
        cgm_bin_io::CgmBinFile bin;           // the .bin as read back from disk
    };

    const CgmRun& cgmRun() {
        static const CgmRun r = [] {
            CgmRun x;
            fs::create_directories(input_output_testdata_dir);
            x.csv = readLibreCsv(inputfile_path.string());
            x.outcome = cgm_selection::processCgmFile(bittiumfile_samplename, inputfile_path, input_output_testdata_dir);
            x.ok = x.outcome.written &&
                cgm_bin_io::readCgmBin(x.outcome.output.string(), x.bin);
            return x;
            }();
        return r;
    }

}  // namespace

// "Warm-up readings are dropped"
TEST(test2f_cgm_bin_creation, WarmupReadingsAreDropped) {
    const CgmRun& r = cgmRun();
    ASSERT_TRUE(r.ok) << "could not process " << fs::absolute(inputfile_path).string()
        << (fs::exists(inputfile_path) ? "" : " (file not found)")
        << (r.outcome.error.empty() ? "" : " (" + r.outcome.error + ")");
    const double firstMs = r.csv.front().tMin * 60000.0;
    const double warmupMs = 12.0 * 60.0 * 60000.0;

    // 09-10 01:21 PM through 09-11 01:21 AM inclusive: 49 readings.
    EXPECT_EQ(r.outcome.result.nWarmupDropped, 49);
    // Nothing in the file is from the first 12 hours.
    for (double t : r.bin.time) EXPECT_GT(t - firstMs, warmupMs);
    // The file starts with the first reading after warm-up: 01:36 AM, 216.
    double t0 = 0;
    ASSERT_TRUE(parseTimestampMinutes("09-11-2026 01:36 AM", t0));
    EXPECT_DOUBLE_EQ(r.bin.time.front(), t0 * 60000.0);
    EXPECT_DOUBLE_EQ(r.bin.raw.front(), 216.0);
}

// "single gaps interpolated"
TEST(test2f_cgm_bin_creation, SingleGapsAreInterpolated) {
    const CgmRun& r = cgmRun();
    ASSERT_TRUE(r.ok);
    const CgmGrid& G = r.outcome.result.grid;

    // Slot 9 (03:51 AM) was removed: filled with the average of 191 and 206.
    EXPECT_EQ(r.bin.flags[9], cgm_bin_io::kFlagInterpolated);
    EXPECT_TRUE(std::isnan(r.bin.raw[9]));
    EXPECT_DOUBLE_EQ(G.value[9], (191.0 + 206.0) / 2.0);

    // Slots 21-22 were removed together: a 2-reading gap stays missing.
    for (int k : { 21, 22 }) {
        EXPECT_EQ(r.bin.flags[k], cgm_bin_io::kFlagMissing);
        EXPECT_TRUE(std::isnan(G.value[k]));
    }
}

// "TIR plus TAR plus TBR equals 100 percent"
TEST(test2f_cgm_bin_creation, TirTarTbrSumTo100) {
    const CgmRun& r = cgmRun();
    ASSERT_TRUE(r.ok);

    // Whole recording (the .bin header).
    EXPECT_NEAR(r.bin.tir + r.bin.tar + r.bin.tbr, 100.0, 1e-9);

    // Every valid bin.
    for (const CgmBin& b : r.outcome.result.bins)
        if (b.raw.valid)
            EXPECT_NEAR(b.raw.tir + b.raw.tar + b.raw.tbr, 100.0, 1e-9) << "bin " << b.index;

    // Bin 13 has one reading in each range (193, 65, 161, 164), so all three
    // parts are non-zero there and still add up.
    const CgmBin& b13 = r.outcome.result.bins[13];
    ASSERT_TRUE(b13.raw.valid);
    EXPECT_DOUBLE_EQ(b13.raw.tar, 25.0);
    EXPECT_DOUBLE_EQ(b13.raw.tbr, 25.0);
    EXPECT_DOUBLE_EQ(b13.raw.tir, 50.0);
}

// "short ... bins are correctly marked invalid"
TEST(test2f_cgm_bin_creation, ShortBinsAreInvalid) {
    const CgmRun& r = cgmRun();
    ASSERT_TRUE(r.ok);

    // The last three readings were removed, so the last bin holds one reading:
    // 15 min, under the 30-min minimum, though nothing in it is missing.
    const CgmBin& last = r.outcome.result.bins.back();
    EXPECT_EQ(last.expected, 1);
    EXPECT_EQ(last.n, 1);
    EXPECT_DOUBLE_EQ(last.missingFrac, 0.0);
    EXPECT_FALSE(last.raw.valid);
}

// "sparse bins are correctly marked invalid"
TEST(test2f_cgm_bin_creation, SparseBinsAreInvalid) {
    const CgmRun& r = cgmRun();
    ASSERT_TRUE(r.ok);
    const auto& bins = r.outcome.result.bins;

    // Bin 5 (slots 20-23) lost slots 21-22: 50% missing, over the 20% limit.
    EXPECT_DOUBLE_EQ(bins[5].missingFrac, 0.5);
    EXPECT_FALSE(bins[5].raw.valid);

    // Bin 2 (slots 8-11) lost only slot 9, which was filled: 0% missing, valid.
    EXPECT_DOUBLE_EQ(bins[2].missingFrac, 0.0);
    EXPECT_TRUE(bins[2].raw.valid);

    // These are the only two invalid bins.
    int invalid = 0;
    for (const CgmBin& b : bins) invalid += b.raw.valid ? 0 : 1;
    EXPECT_EQ(invalid, 2);
}

// Not a spec statement: the .bin written now matches the baseline .bin that
// was made from the same CSV, so any change in output is caught. Values are
// compared to 1e-9 (NaN equals NaN) rather than byte for byte, so compiler
// rounding differences do not count as a change.
TEST(test2f_cgm_bin_creation, OutputMatchesBaseline) {
    const CgmRun& r = cgmRun();
    ASSERT_TRUE(r.ok);
    cgm_bin_io::CgmBinFile base;
    std::string err;
    ASSERT_TRUE(cgm_bin_io::readCgmBin(baseline_path.string(), base, &err))
        << baseline_path.string() << ": " << err;

    auto same = [](double a, double b) {
        return (std::isnan(a) && std::isnan(b)) || std::fabs(a - b) <= 1e-9 * std::max(1.0, std::fabs(b));
        };
    const cgm_bin_io::CgmBinFile& now = r.bin;
    EXPECT_TRUE(same(now.mean, base.mean));
    EXPECT_TRUE(same(now.sd, base.sd));
    EXPECT_TRUE(same(now.cv, base.cv));
    EXPECT_TRUE(same(now.roc, base.roc));
    EXPECT_TRUE(same(now.tar, base.tar));
    EXPECT_TRUE(same(now.tbr, base.tbr));
    EXPECT_TRUE(same(now.tir, base.tir));
    ASSERT_EQ(now.size(), base.size());
    for (std::size_t k = 0; k < now.size(); ++k) {
        EXPECT_TRUE(same(now.time[k], base.time[k])) << "time, slot " << k;
        EXPECT_TRUE(same(now.raw[k], base.raw[k])) << "raw, slot " << k;
        EXPECT_TRUE(same(now.normalized[k], base.normalized[k])) << "normalized, slot " << k;
        EXPECT_EQ(now.flags[k], base.flags[k]) << "flags, slot " << k;
    }
}