// test2a_sqi.cpp -- acceptance tests for the ECG and PPG SQIs (Task A).
//
// ECG: "On a clean record, most beats score above 0.70. Inject baseline wander
// and noise into known beats; confirm their composite SQI drops into the
// substitute or exclude bands as expected, and that the subsegmental
// chi-squared rises in the affected segment only."
//
// PPG: "A clean PPG record scores above 0.70 on most pulses. Injecting a
// saturating amplitude step drops SQI_clipping specifically while leaving
// SQI_template high, which confirms the components are independent rather
// than co-varying."
//
// DATA. Synthetic, built below, so the right answer is known.

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "logging/sqi_ecg.hpp"
#include "logging/sqi_ppg.hpp"

namespace {

    constexpr double kPi = 3.14159265358979323846;
    constexpr double kEcgFs = 500.0;

    // Deterministic pseudo-random noise in [-1, 1].
    struct Lcg {
        unsigned s;
        double next() { s = s * 1103515245u + 12345u; return ((s >> 8) & 0xFFFF) / 32767.5 - 1.0; }
    };

    // ---- ECG --------------------------------------------------------------
    // One beat, 0.8 s: P at 0.10 s, R at 0.25 s, T at 0.45 s, flat TP after.
    std::vector<double> ecgTemplate() {
        const int n = static_cast<int>(0.8 * kEcgFs);
        std::vector<double> v(n, 0.0);
        auto bump = [&](double c, double w, double a) {
            for (int i = 0; i < n; ++i) v[i] += a * std::exp(-0.5 * std::pow((i / kEcgFs - c) / w, 2));
            };
        bump(0.10, 0.020, 0.15);
        bump(0.25, 0.008, 1.20);
        bump(0.45, 0.040, 0.30);
        return v;
    }
    Segments ecgSegments() {
        Segments s;
        auto at = [](double t) { return static_cast<int>(t * kEcgFs); };
        s.pLo = at(0.04);  s.pHi = at(0.16);
        s.qrsLo = at(0.21); s.qrsHi = at(0.29);
        s.stLo = s.qrsHi;  s.stHi = at(0.38);
        s.tHi = at(0.60);  s.nextPLo = at(0.78);   // TP segment 0.60-0.78 s
        return s;
    }
    std::vector<double> absOf(std::vector<double> v) { for (double& x : v) x = std::fabs(x); return v; }

    // A clean beat: the template plus 5 uV of noise.
    std::vector<double> cleanEcgBeat(Lcg& r) {
        std::vector<double> b = ecgTemplate();
        for (double& x : b) x += 0.005 * r.next();
        return b;
    }

    BeatSQI score(const std::vector<double>& beat) {
        const std::vector<double> t = ecgTemplate();
        return computeEcgSQI(beat, t, absOf(t), ecgSegments(), -1, kEcgFs);
    }

    // ---- PPG --------------------------------------------------------------
    // Zero-centred pulse, 0.8 s at 125 Hz: systolic peak plus dicrotic wave,
    // quantized to a 1/1024 step on an ADC spanning +-0.5.
    constexpr double kPpgFs = 125.0, kRail = 0.5, kStep = 1.0 / 1024.0;
    std::vector<double> ppgPulse(double gain, Lcg& r) {
        const int n = static_cast<int>(0.8 * kPpgFs);
        std::vector<double> v(n);
        for (int i = 0; i < n; ++i) {
            const double t = i / kPpgFs;
            double x = 0.15 * std::exp(-0.5 * std::pow((t - 0.20) / 0.06, 2))
                + 0.06 * std::exp(-0.5 * std::pow((t - 0.42) / 0.05, 2)) - 0.06;
            x = gain * x + 0.002 * r.next();
            x = std::max(-kRail, std::min(kRail, x));           // ADC clips at the rails
            v[i] = std::round(x / kStep) * kStep;                 // ADC quantizes
        }
        return v;
    }

}  // namespace

// ===========================================================================
//  ECG
// ===========================================================================

// "On a clean record, most beats score above 0.70."
TEST(test2a_sqi, EcgCleanRecordMostBeatsAbove070) {
    Lcg r{ 7 };
    int above = 0;
    const int nBeats = 200;
    for (int b = 0; b < nBeats; ++b) above += score(cleanEcgBeat(r)).composite > 0.70;
    EXPECT_GE(above, nBeats * 9 / 10);
}

// "Inject baseline wander ... drops into the substitute or exclude bands."
TEST(test2a_sqi, EcgBaselineWanderDropsBand) {
    Lcg r{ 11 };
    const Segments s = ecgSegments();
    auto wander = [&](double mv) {              // linear drift of mv across the beat
        std::vector<double> b = cleanEcgBeat(r);
        for (std::size_t i = 0; i < b.size(); ++i) b[i] += mv * i / (s.tHi - s.pLo);
        return score(b);
        };
    const BeatSQI clean = score(cleanEcgBeat(r));
    const BeatSQI mild = wander(0.15);          // 0.15 mV P-to-T drift
    const BeatSQI strong = wander(0.40);        // 0.4 mV
    EXPECT_EQ(clean.handling, BeatSQI::INCLUDE);
    EXPECT_EQ(mild.handling, BeatSQI::SUBSTITUTE);
    EXPECT_EQ(strong.handling, BeatSQI::EXCLUDE);
    EXPECT_LT(mild.baseline, clean.baseline);
}

// "Inject ... noise ... drops into the substitute or exclude bands."
TEST(test2a_sqi, EcgNoiseDropsBand) {
    Lcg r{ 13 };
    const Segments s = ecgSegments();
    auto noisy = [&](double sd) {               // noise added to the whole beat
        std::vector<double> b = cleanEcgBeat(r);
        for (double& x : b) x += sd * std::sqrt(3.0) * r.next();   // uniform with this SD
        return score(b);
        };
    const BeatSQI mild = noisy(0.025);
    const BeatSQI strong = noisy(0.08);
    EXPECT_EQ(mild.handling, BeatSQI::SUBSTITUTE);
    EXPECT_EQ(strong.handling, BeatSQI::EXCLUDE);
    (void)s;
}

// "... the subsegmental chi-squared rises in the affected segment only."
TEST(test2a_sqi, EcgChiSquaredRisesInAffectedSegmentOnly) {
    const std::vector<double> t = ecgTemplate();
    const Segments s = ecgSegments();
    const BeatSQI base = computeEcgSQI(t, t, absOf(t), s, -1, kEcgFs);

    // Distort one segment at a time; only that segment's chi-squared moves.
    struct Case { int lo, hi; const char* name; };
    const Case cases[3] = { { s.pLo, s.pHi, "P" }, { s.qrsLo, s.qrsHi, "QRS" }, { s.stLo, s.stHi, "ST" } };
    for (const Case& c : cases) {
        std::vector<double> b = t;
        for (int i = c.lo; i < c.hi; ++i) b[i] += 0.1;
        const BeatSQI q = computeEcgSQI(b, t, absOf(t), s, -1, kEcgFs);
        const double z[3] = { q.chiSq0_P, q.chiSq0_QRS, q.chiSq0_ST };
        const double a[3] = { q.chiSqAbs_P, q.chiSqAbs_QRS, q.chiSqAbs_ST };
        const double z0[3] = { base.chiSq0_P, base.chiSq0_QRS, base.chiSq0_ST };
        const double a0[3] = { base.chiSqAbs_P, base.chiSqAbs_QRS, base.chiSqAbs_ST };
        for (int k = 0; k < 3; ++k) {
            const bool affected = (&cases[k] == &c);
            if (affected) {
                EXPECT_GT(z[k], z0[k]) << c.name;
                EXPECT_GT(a[k], a0[k]) << c.name << " (abs)";
            }
            else {
                EXPECT_DOUBLE_EQ(z[k], z0[k]) << c.name << " moved segment " << cases[k].name;
                EXPECT_DOUBLE_EQ(a[k], a0[k]) << c.name << " moved segment " << cases[k].name << " (abs)";
            }
        }
    }
}

// Not an acceptance statement: against the absolute-value template the beat
// is compared as |beat|, so a beat identical to the template scores 0 on both
// templates, even though this beat has negative deflections.
TEST(test2a_sqi, EcgPerfectBeatScoresZeroChiSquaredOnBothTemplates) {
    std::vector<double> t = ecgTemplate();
    for (std::size_t i = 0; i < t.size(); ++i)                       // add an S wave
        t[i] -= 0.4 * std::exp(-0.5 * std::pow((i / kEcgFs - 0.27) / 0.006, 2));
    const BeatSQI q = computeEcgSQI(t, t, absOf(t), ecgSegments(), -1, kEcgFs);
    EXPECT_DOUBLE_EQ(q.chiSq0, 0.0);
    EXPECT_DOUBLE_EQ(q.chiSqAbs, 0.0);
    EXPECT_DOUBLE_EQ(q.chiSqAbs_QRS, 0.0);
}

// ===========================================================================
//  PPG
// ===========================================================================

// "A clean PPG record scores above 0.70 on most pulses."
TEST(test2a_sqi, PpgCleanRecordMostPulsesAbove070) {
    Lcg r{ 17 };
    std::vector<std::vector<double>> pulses;
    std::vector<double> record;
    for (int p = 0; p < 200; ++p) {
        pulses.push_back(ppgPulse(1.0, r));
        record.insert(record.end(), pulses.back().begin(), pulses.back().end());
    }
    const sqi_ppg::AdcRails rails = sqi_ppg::detectAdcRails(record);
    Lcg rt{ 99 };
    const std::vector<double> tmpl200 = sqi_ppg::resampleTo(ppgPulse(1.0, rt), 200);

    int above = 0;
    for (const auto& p : pulses)
        above += sqi_ppg::computePpgSQI(p, tmpl200, rails.loThreshold(), rails.hiThreshold(), 0.0).composite > 0.70;
    EXPECT_GE(above, 180);
}

// "Injecting a saturating amplitude step drops SQI_clipping specifically while
// leaving SQI_template high."
TEST(test2a_sqi, PpgSaturatingStepDropsClippingOnly) {
    Lcg r{ 23 };
    std::vector<std::vector<double>> pulses;
    std::vector<double> record;
    for (int p = 0; p < 200; ++p) {
        // Pulses 100-119: gain 8, which drives the peak past the +0.5 rail.
        pulses.push_back(ppgPulse((p >= 100 && p < 120) ? 8.0 : 1.0, r));
        record.insert(record.end(), pulses.back().begin(), pulses.back().end());
    }
    const sqi_ppg::AdcRails rails = sqi_ppg::detectAdcRails(record);
    ASSERT_NEAR(rails.hi, kRail, kStep);        // the rail is found from the record
    Lcg rt{ 99 };
    const std::vector<double> tmpl200 = sqi_ppg::resampleTo(ppgPulse(1.0, rt), 200);

    auto q = [&](int p) {
        return sqi_ppg::computePpgSQI(pulses[p], tmpl200, rails.loThreshold(), rails.hiThreshold(), 0.0);
        };
    const sqi_ppg::PulseSQI clean = q(50), clipped = q(110);
    EXPECT_DOUBLE_EQ(clean.clipping, 1.0);
    EXPECT_LT(clipped.clipping, 0.95);           // clipping drops ...
    EXPECT_GT(clipped.templateCorr, 0.90);       // ... the template term stays high
    EXPECT_GT(clean.templateCorr, 0.90);
    EXPECT_LT(clipped.clipping, clean.clipping);
    EXPECT_STREQ(sqi_ppg::exclusionReason(clean), "");
}

// Not an acceptance statement: the components feed the composite with the
// spec weights, the 0.70 cut decides inclusion, and an excluded pulse names
// the component that cost it the most.
TEST(test2a_sqi, PpgCompositeAndExclusionReason) {
    Lcg r{ 29 };
    const std::vector<double> pulse = ppgPulse(1.0, r);
    Lcg rt{ 99 };
    const std::vector<double> tmpl200 = sqi_ppg::resampleTo(ppgPulse(1.0, rt), 200);
    const sqi_ppg::PulseSQI q = sqi_ppg::computePpgSQI(pulse, tmpl200, -1e9, 1e9, 0.0);
    EXPECT_NEAR(q.composite, 0.40 * q.templateCorr + 0.25 * (1 - q.chiSq) + 0.20 * q.clipping + 0.15 * q.motion, 1e-12);
    EXPECT_TRUE(q.include);

    // A flat pulse has no shape: excluded for its template term.
    const sqi_ppg::PulseSQI flat = sqi_ppg::computePpgSQI(std::vector<double>(100, 0.0), tmpl200, -1e9, 1e9, 0.0);
    EXPECT_FALSE(flat.include);
    EXPECT_STREQ(sqi_ppg::exclusionReason(flat), "template");

    // Motion (accelerometer flag 0 = motion) costs the motion term only.
    const sqi_ppg::PulseSQI moving = sqi_ppg::computePpgSQI(pulse, tmpl200, -1e9, 1e9, sqi_ppg::ppgMotionFlag(0));
    EXPECT_DOUBLE_EQ(moving.motion, 0.0);
    EXPECT_DOUBLE_EQ(moving.templateCorr, q.templateCorr);
    EXPECT_NEAR(q.composite - moving.composite, 0.15, 1e-12);
}