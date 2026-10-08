// test2g_accel_pipeline.cpp -- acceptance tests for the accelerometry pipeline.
//
// One test per statement of the spec's acceptance test:
//   "Known resting versus active segments classify correctly, and VM above
//    1.05 g raises the motion flag that lowers the corresponding ECG and PPG
//    SQI (Task A)."
//
// DATA. data/accel/sample_accel_file.bin, next to this .cpp file: the real
// Bittium recording 024-00-11-07-29.bin (248 s, accelerometer at 25 Hz in
// milli-g), with its accelerometer edited so every case is present:
//
//   epoch  time       edit                                   result
//   0      0-30 s     real movement, amplified 2.5x          MODERATE_VIGOROUS, motion
//   1-3    30-120 s   none                                   LIGHT, motion
//   4-5    120-180 s  real rest, 1-s moving average          SLEEP_SUPINE, no motion
//                     (removes the ~0.02 g sensor noise)
//   6-7    180-240 s  as 4-5, then Y and Z swapped           SEATED_REST, no motion
//                     (gravity on Y instead of Z)
//   8      240-248 s  none; only 8 s of samples              invalid epoch
//
// The edits are to slots 5-7 only (upsampled and raw values alike); the rest
// of the file is as recorded. Bittium stores milli-g, so unitsPerG = 1000.

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include "accel/accel_pipeline.hpp"
#include "logging/sqi_ecg.hpp"
#include "logging/sqi_ppg.hpp"

using namespace accel_pipeline;
namespace fs = std::filesystem;

namespace {

#ifdef ACCEL_DATA_DIR
    const fs::path kDataDir = ACCEL_DATA_DIR;
#else
    const fs::path kDataDir = fs::path(__FILE__).parent_path() / "data" / "accel";
#endif
    const fs::path kFile = kDataDir / "sample_accel_file.bin";

    AccelParams bittium() {
        AccelParams p;
        p.unitsPerG = 1000.0;   // milli-g
        return p;
    }

    const AccelData& data() {
        static const AccelData a = readAccelFromBin(kFile.string());
        return a;
    }
    const AccelResult& result() {
        static const AccelResult r = runAccelPipeline(data(), bittium());
        return r;
    }

    // Every test starts with this, so a missing file fails with its path
    // instead of indexing into empty data.
#define REQUIRE_FILE() \
    ASSERT_TRUE(data().present()) << "could not read " << fs::absolute(kFile).string(); \
    ASSERT_EQ(result().epochs.size(), 9u)

    // Time of the first sample flagged as motion, or -1.
    double firstMotionSec(const AccelResult& r) {
        for (std::size_t i = 0; i < r.motion.size(); ++i)
            if (r.motion[i]) return static_cast<double>(i) / r.fs;
        return -1.0;
    }

    // A clean synthetic ECG beat and its segments, so computeEcgSQI scores it
    // near 1 before any motion term.
    void cleanBeat(std::vector<double>& beat, Segments& seg, double fs) {
        const int n = static_cast<int>(0.8 * fs);
        beat.assign(n, 0.0);
        auto bump = [&](double c, double w, double a) {
            for (int i = 0; i < n; ++i) beat[i] += a * std::exp(-0.5 * std::pow((i / fs - c) / w, 2));
            };
        bump(0.10, 0.020, 0.15);
        bump(0.25, 0.008, 1.20);
        bump(0.45, 0.040, 0.30);
        seg.pLo = static_cast<int>(0.06 * fs); seg.pHi = static_cast<int>(0.14 * fs);
        seg.qrsLo = static_cast<int>(0.22 * fs); seg.qrsHi = static_cast<int>(0.28 * fs);
        seg.stLo = seg.qrsHi; seg.stHi = static_cast<int>(0.38 * fs);
        seg.tHi = static_cast<int>(0.60 * fs); seg.nextPLo = seg.tHi;
    }

}  // namespace

// "Known resting versus active segments classify correctly"
TEST(test2g_accel_pipeline, KnownSegmentsClassifyCorrectly) {
    REQUIRE_FILE();   // 8 full 30-s epochs + 8 s
    const AccelResult& r = result();

    const Activity want[8] = {
        Activity::MODERATE_VIGOROUS,
        Activity::LIGHT, Activity::LIGHT, Activity::LIGHT,
        Activity::SLEEP_SUPINE, Activity::SLEEP_SUPINE,
        Activity::SEATED_REST, Activity::SEATED_REST };
    for (int e = 0; e < 8; ++e) {
        ASSERT_TRUE(r.epochs[e].valid) << "epoch " << e;
        EXPECT_EQ(r.epochs[e].activity, want[e])
            << "epoch " << e << ": got " << activityName(r.epochs[e].activity)
            << ", meanVM " << r.epochs[e].meanVM << ", sdVM " << r.epochs[e].sdVM;
    }

    // activityAt() gives the same answer by time.
    Activity a;
    ASSERT_TRUE(r.activityAt(10.0, a));  EXPECT_EQ(a, Activity::MODERATE_VIGOROUS);
    ASSERT_TRUE(r.activityAt(150.0, a)); EXPECT_EQ(a, Activity::SLEEP_SUPINE);
    ASSERT_TRUE(r.activityAt(200.0, a)); EXPECT_EQ(a, Activity::SEATED_REST);
}

// The last 8 s are under 80% of a 30-s epoch: invalid, and no activity is
// given for that time.
TEST(test2g_accel_pipeline, PartialEpochIsInvalid) {
    REQUIRE_FILE();
    const AccelResult& r = result();
    const Epoch& last = r.epochs.back();
    EXPECT_EQ(last.n, 200);          // 8 s at 25 Hz
    EXPECT_EQ(last.expected, 750);   // 30 s at 25 Hz
    EXPECT_FALSE(last.valid);
    Activity a;
    EXPECT_FALSE(r.activityAt(245.0, a));
}

// "VM above 1.05 g raises the motion flag"
TEST(test2g_accel_pipeline, VmAbove105RaisesMotionFlag) {
    REQUIRE_FILE();
    const AccelResult& r = result();

    // VM is sqrt(X^2 + Y^2 + Z^2), in g after the milli-g conversion.
    const AccelData& a = data();
    const auto vm = vectorMagnitude(a, bittium());
    const std::size_t k = 800;   // an unedited sample
    const double x = a.x[k] / 1000.0, y = a.y[k] / 1000.0, z = a.z[k] / 1000.0;
    EXPECT_DOUBLE_EQ(vm[k], std::sqrt(x * x + y * y + z * z));

    // Exactly at the threshold is not motion; just above is.
    AccelData edge; edge.fs = 1.0;
    edge.x = { 0.0, 0.0 }; edge.y = { 0.0, 0.0 }; edge.z = { 1.05, 1.06 };
    const AccelResult e = runAccelPipeline(edge);
    EXPECT_EQ(e.motion[0], 0);
    EXPECT_EQ(e.motion[1], 1);

    // Every flagged sample has VM > 1.05 g, and every sample above is flagged.
    for (std::size_t i = 0; i < r.vm.size(); ++i)
        EXPECT_EQ(r.motion[i] != 0, r.vm[i] > 1.05) << "sample " << i;

    // Active epochs flag motion; resting epochs do not.
    for (int ep = 0; ep < 4; ++ep) EXPECT_TRUE(r.epochs[ep].motion) << "epoch " << ep;
    for (int ep = 4; ep < 8; ++ep) EXPECT_FALSE(r.epochs[ep].motion) << "epoch " << ep;

    // Per interval, in the SQI convention.
    const double t = firstMotionSec(r);
    ASSERT_GE(t, 0.0);
    EXPECT_EQ(r.sqiMotionFlag(t - 0.4, t + 0.4), 0);       // around a flagged sample
    EXPECT_EQ(r.sqiMotionFlag(150.0, 150.8), 1);           // resting
    EXPECT_EQ(r.sqiMotionFlag(500.0, 500.8), -1);          // past the recording

    // Share of an interval that is motion (the PPG SQI's input).
    EXPECT_DOUBLE_EQ(r.motionFraction(150.0, 180.0), 0.0);
    EXPECT_NEAR(r.motionFraction(0.0, 30.0 - 1.0 / r.fs), r.epochs[0].motionFrac, 1e-12);
    EXPECT_GT(r.motionFraction(0.0, 30.0), 0.0);
    EXPECT_LT(r.motionFraction(0.0, 30.0), 1.0);
    EXPECT_TRUE(std::isnan(r.motionFraction(500.0, 501.0)));
}

// "... the motion flag that lowers the corresponding ECG ... SQI"
TEST(test2g_accel_pipeline, MotionFlagLowersEcgSqi) {
    REQUIRE_FILE();
    const AccelResult& r = result();
    const double fs = 500.0;
    std::vector<double> beat;
    Segments seg;
    cleanBeat(beat, seg, fs);

    // The same perfect beat, once at a resting time and once over a flagged sample.
    const double t = firstMotionSec(r);
    ASSERT_GE(t, 0.0);
    const int restFlag = r.sqiMotionFlag(150.0, 150.8);
    const int motionFlag = r.sqiMotionFlag(t - 0.4, t + 0.4);
    ASSERT_EQ(restFlag, 1);
    ASSERT_EQ(motionFlag, 0);

    const BeatSQI rest = computeEcgSQI(beat, beat, beat, seg, restFlag, fs);
    const BeatSQI moving = computeEcgSQI(beat, beat, beat, seg, motionFlag, fs);
    EXPECT_GT(rest.composite, 0.9);
    EXPECT_EQ(rest.handling, BeatSQI::INCLUDE);
    EXPECT_LT(moving.composite, rest.composite);
    EXPECT_DOUBLE_EQ(moving.composite, 0.0);
    EXPECT_EQ(moving.handling, BeatSQI::EXCLUDE);

    // No accelerometer coverage (-1) leaves the score as it was.
    const BeatSQI unknown = computeEcgSQI(beat, beat, beat, seg, r.sqiMotionFlag(500.0, 500.8), fs);
    EXPECT_DOUBLE_EQ(unknown.composite, rest.composite);
}

// "... the motion flag that lowers the corresponding ... PPG SQI"
//
// The real computePpgSQI, fed the accelerometer's motion share of each pulse
// window (motionFraction), as writePpgSQICsv does. The same clean pulse scores
// lower over a window with flagged samples; only the motion term changes.
TEST(test2g_accel_pipeline, MotionFlagLowersPpgSqi) {
    REQUIRE_FILE();
    const AccelResult& r = result();
    const double t = firstMotionSec(r);
    ASSERT_GE(t, 0.0);

    // A clean synthetic pulse, identical to its template (125 Hz, 0.8 s).
    std::vector<double> pulse(100);
    for (int i = 0; i < 100; ++i) {
        const double s = i / 125.0;
        pulse[i] = 0.15 * std::exp(-0.5 * std::pow((s - 0.20) / 0.06, 2))
            + 0.06 * std::exp(-0.5 * std::pow((s - 0.42) / 0.05, 2));
    }
    const std::vector<double> tmpl200 = sqi_ppg::resampleTo(pulse, 200);

    const double restShare = r.motionFraction(150.0, 151.0);     // resting
    const double moveShare = r.motionFraction(t - 0.4, t + 0.4);  // around a flagged sample
    ASSERT_DOUBLE_EQ(restShare, 0.0);
    ASSERT_GT(moveShare, 0.0);

    const sqi_ppg::PulseSQI rest = sqi_ppg::computePpgSQI(pulse, tmpl200, -1e9, 1e9, restShare);
    const sqi_ppg::PulseSQI moving = sqi_ppg::computePpgSQI(pulse, tmpl200, -1e9, 1e9, moveShare);
    EXPECT_DOUBLE_EQ(rest.motion, 1.0);
    EXPECT_DOUBLE_EQ(moving.motion, 1.0 - moveShare);
    EXPECT_LT(moving.composite, rest.composite);
    EXPECT_NEAR(rest.composite - moving.composite, 0.15 * moveShare, 1e-12);   // motion weight
    EXPECT_DOUBLE_EQ(moving.templateCorr, rest.templateCorr);                   // nothing else moves
    EXPECT_DOUBLE_EQ(moving.clipping, rest.clipping);

    // No accelerometer coverage: no penalty, as writePpgSQICsv treats it.
    EXPECT_TRUE(std::isnan(r.motionFraction(500.0, 501.0)));
}

// The composite multiplier for a motion flag (sqiMotionTerm), the rule the
// ECG SQI applies: 1 for clean (1) and unavailable (-1), 0 for motion (0).
// Fed the real recording's flags at rest, around a flagged sample, and past
// the end.
TEST(test2g_accel_pipeline, SqiMotionTermFromFlag) {
    REQUIRE_FILE();
    const AccelResult& r = result();
    const double t = firstMotionSec(r);
    ASSERT_GE(t, 0.0);
    const double score = 0.85;   // any composite before motion
    const double rest = score * sqiMotionTerm(r.sqiMotionFlag(150.0, 151.0));
    const double moving = score * sqiMotionTerm(r.sqiMotionFlag(t - 0.4, t + 0.4));
    const double unknown = score * sqiMotionTerm(r.sqiMotionFlag(500.0, 501.0));
    EXPECT_DOUBLE_EQ(rest, score);
    EXPECT_LT(moving, rest);
    EXPECT_DOUBLE_EQ(moving, 0.0);
    EXPECT_DOUBLE_EQ(unknown, score);
}

// The Goal's percentile normalization: epoch meanVM and sdVM rescaled by the
// recording's own 2nd / 98th percentiles over valid epochs; invalid epochs
// get none.
TEST(test2g_accel_pipeline, PercentileNormalization) {
    REQUIRE_FILE();
    const AccelResult& r = result();
    for (const Epoch& e : r.epochs) {
        if (!e.valid) {
            EXPECT_TRUE(std::isnan(e.nMeanVM));
            EXPECT_TRUE(std::isnan(e.nSdVM));
            continue;
        }
        EXPECT_NEAR(e.nMeanVM, (e.meanVM - r.meanP2) / (r.meanP98 - r.meanP2), 1e-12);
        EXPECT_NEAR(e.nSdVM, (e.sdVM - r.sdP2) / (r.sdP98 - r.sdP2), 1e-12);
    }
    // Vigorous movement at the top of the scale, rest at the bottom.
    EXPECT_GE(r.epochs[0].nSdVM, 1.0);
    for (int ep = 4; ep < 8; ++ep) EXPECT_LT(r.epochs[ep].nSdVM, 0.05) << "epoch " << ep;
}

// Reading the real file: slots 5-7 at the header's rate, in milli-g.
TEST(test2g_accel_pipeline, ReadsAccelFromDataBin) {
    REQUIRE_FILE();
    const AccelData& a = data();
    EXPECT_DOUBLE_EQ(a.fs, 25.0);
    EXPECT_EQ(a.x.size(), 6200u);   // 248 s
    // An unedited sample, as recorded.
    EXPECT_DOUBLE_EQ(a.x[800], 460.5);
    EXPECT_DOUBLE_EQ(a.y[800], 628.5);
    EXPECT_DOUBLE_EQ(a.z[800], -584.0);
    // Resting VM is about 1 g once converted from milli-g (about 1000 if not).
    EXPECT_NEAR(result().epochs[4].meanVM, 1.0, 0.05);
    EXPECT_NEAR(runAccelPipeline(a).epochs[4].meanVM, 1000.0, 50.0);
}

// The config settings: epoch length and valid-epoch percent change the epochs.
TEST(test2g_accel_pipeline, EpochSettingsFromConfig) {
    REQUIRE_FILE();
    AccelParams p = bittium();
    p.epochSec = 60.0;                       // accel_epoch_sec
    const AccelResult r60 = runAccelPipeline(data(), p);
    ASSERT_EQ(r60.epochs.size(), 5u);        // 4 x 60 s + 8 s
    EXPECT_EQ(r60.epochs[0].expected, 1500);
    EXPECT_FALSE(r60.epochs.back().valid);

    p = bittium();
    p.minEpochCoverage = 0.20;               // accel_valid_epoch_pct = 20
    const AccelResult r20 = runAccelPipeline(data(), p);
    EXPECT_TRUE(r20.epochs.back().valid);    // 200 of 750 samples = 27%
}