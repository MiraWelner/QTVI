/**
 * @file   test1b_upsampling.cpp
 * @brief  Acceptance test for filterutils::upsample.
 *
 *   Input: sinusoid_100hz_256hz.csv -- 10 s of 1.0 * sin(2*pi*100*t + 0.5),
 *   sampled at 256 Hz, columns time_s,value. Put it next to this file, in the
 *   project folder (filetbin_tester\).
 *
 *   1. Upsampled 256 -> 1000 Hz, the output frequency, amplitude and phase
 *      match the input within the FIR passband tolerance.
 *   2. Upsampled to each modality's rate (ECG 1000 Hz, PPG 500 Hz), the
 *      output timestamps agree with the CSV's timestamps, and with each
 *      other, within one sample period of the fastest channel (1 ms).
 *
 *   FIR passband tolerance: each output sample comes from one sub-filter p of
 *   the polyphase bank; at 100 Hz that sub-filter has gain |H_p| and phase
 *   arg H_p. The tolerance is the worst of these over all p, computed from
 *   the same bank upsample() uses (buildPolyphaseBank).
 */

#include "pch.h"          // Visual Studio precompiled header: must be the first include
#include <gtest/gtest.h>

#include "../filter_utils.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

    // How sinusoid_100hz_256hz.csv was generated.
    constexpr double kFs = 256.0, kF = 100.0, kA = 1.0, kPhi = 0.5;
    constexpr double kTwoPi = 2.0 * M_PI;

    struct Csv { std::vector<double> t, v; std::string path; };

    // Where to look for the CSV: the folder this source file is in, then the
    // working directory and up to four folders above it. Visual Studio's
    // Ctrl+F5 runs from the project folder; Test Explorer runs from
    // x64\Release\, two folders below it. Either way the CSV is found.
    std::vector<std::filesystem::path> csvCandidates() {
        namespace fs = std::filesystem;
        const char* name = "sinusoid_100hz_256hz.csv";
        std::vector<fs::path> c;
        c.push_back(fs::path(__FILE__).parent_path() / name);
        std::error_code ec;
        fs::path dir = fs::current_path(ec);
        for (int up = 0; up <= 4 && !ec; ++up) {
            c.push_back(dir / name);
            if (!dir.has_parent_path() || dir.parent_path() == dir) break;
            dir = dir.parent_path();
        }
        return c;
    }

    bool readCsv(Csv& c, std::string& tried) {
        for (const auto& p : csvCandidates()) {
            tried += "\n    " + p.string();
            std::ifstream f(p);
            if (!f) continue;
            std::string line;
            std::getline(f, line);                       // header
            while (std::getline(f, line)) {
                std::stringstream ss(line);
                std::string a, b;
                if (std::getline(ss, a, ',') && std::getline(ss, b)) {
                    c.t.push_back(std::stod(a));
                    c.v.push_back(std::stod(b));
                }
            }
            c.path = p.string();
            return !c.v.empty();
        }
        return false;
    }

    // Least squares y ~ A sin(2 pi f t + phi), t = m / fs, over samples [m0, m1).
    void fitSine(const std::vector<double>& y, double fs, double f, size_t m0, size_t m1,
        double& A, double& phi)
    {
        double ss = 0, cc = 0, sc = 0, sy = 0, cy = 0;
        for (size_t m = m0; m < m1; ++m) {
            const double a = kTwoPi * f * m / fs, s = std::sin(a), c = std::cos(a);
            ss += s * s; cc += c * c; sc += s * c; sy += s * y[m]; cy += c * y[m];
        }
        const double det = ss * cc - sc * sc;
        const double ks = (sy * cc - cy * sc) / det, kc = (cy * ss - sy * sc) / det;
        A = std::hypot(ks, kc);
        phi = std::atan2(kc, ks);
    }

    double wrapPi(double x) { return std::remainder(x, kTwoPi); }

    // FIR passband tolerance at f for sourceRate -> targetRate (see file header).
    void passbandTolerance(double sourceRate, double targetRate, double f, double& dA, double& dPhi) {
        const int g = filterutils::greatest_common_divisor((int)targetRate, (int)sourceRate);
        const int P = (int)targetRate / g, Q = (int)sourceRate / g;
        const int halfLobes = std::max(16, std::max(P, Q) / 2);
        const auto bank = filterutils::buildPolyphaseBank(P, Q, halfLobes);
        const int center = halfLobes * std::max(P, Q) / P;
        const double w = kTwoPi * f / sourceRate;
        dA = dPhi = 0.0;
        for (int p = 0; p < P; ++p) {
            std::complex<double> H = 0.0;
            for (size_t k = 0; k < bank[p].size(); ++k)
                H += bank[p][k] * std::polar(1.0, w * (center - (double)k - (double)p / P));
            dA = std::max(dA, std::fabs(std::abs(H) - 1.0));
            dPhi = std::max(dPhi, std::fabs(std::arg(H)));
        }
    }

    class Upsampling : public ::testing::Test {
    protected:
        void SetUp() override {
            std::string tried;
            ASSERT_TRUE(readCsv(csv, tried)) << "sinusoid_100hz_256hz.csv not found; looked in:" << tried;
        }
        Csv csv;
    };

} // namespace

// ---- 1. Frequency, amplitude, phase -----------------------------------------

TEST_F(Upsampling, Sine100HzFrom256To1000MatchesWithinPassbandTolerance) {
    const double target = 1000.0;
    double dA = 0, dPhi = 0;
    passbandTolerance(kFs, target, kF, dA, dPhi);
    ASSERT_LT(dA, 1e-3) << "100 Hz is not in the FIR passband";
    ASSERT_LT(dPhi, 1e-3) << "100 Hz is not in the FIR passband";

    const std::vector<double> y = filterutils::upsample(csv.v, kFs, target);
    ASSERT_EQ(y.size(), (size_t)std::ceil(csv.v.size() * target / kFs));

    // Skip 1 s at each end, where the filter runs off the data.
    const size_t m0 = (size_t)target, m1 = y.size() - (size_t)target, mid = (m0 + m1) / 2;
    const double tiny = 1e-12;

    double A = 0, phi = 0;
    fitSine(y, target, kF, m0, m1, A, phi);
    EXPECT_NEAR(A, kA, kA * dA + tiny) << "amplitude";
    EXPECT_NEAR(wrapPi(phi - kPhi), 0.0, dPhi + tiny) << "phase (rad)";

    // Frequency from the phase change between the two halves: each half's
    // phase is within dA + dPhi, so the frequency is within
    // 2 (dA + dPhi) / (2 pi * time between the half centres).
    double A1, p1, A2, p2;
    fitSine(y, target, kF, m0, mid, A1, p1);
    fitSine(y, target, kF, mid, m1, A2, p2);
    const double dt = (double)(m1 - m0) / 2.0 / target;
    const double f = kF + wrapPi(p2 - p1) / (kTwoPi * dt);
    EXPECT_NEAR(f, kF, 2.0 * (dA + dPhi) / (kTwoPi * dt) + tiny) << "frequency (Hz)";

    std::printf("  tolerance dA=%.2e dPhi=%.2e rad | A=%.9f phi=%.9f f=%.9f Hz\n", dA, dPhi, A, phi, f);
}

// ---- 2. Timestamps across modalities ----------------------------------------

TEST_F(Upsampling, TimestampsAcrossModalitiesAlignWithinOneFastestSample) {
    // Modalities sampled at 256 Hz and their upsample rates (config.csv, MESA).
    struct Modality { const char* name; double rate; };
    const Modality mods[] = { { "ecg", 1000.0 }, { "ppg", 500.0 } };
    const double fastest = 1000.0, tol = 1.0 / fastest;   // one sample of the fastest channel

    // The CSV's own timestamps are the reference clock: sample n at n / 256 s.
    for (size_t n = 0; n < csv.t.size(); ++n)
        ASSERT_NEAR(csv.t[n], n / kFs, 1e-12) << "CSV timestamps are not on a 256 Hz grid";

    std::vector<double> offset;
    for (const auto& m : mods) {
        const std::vector<double> y = filterutils::upsample(csv.v, kFs, m.rate);
        // Output sample k sits at t = k / rate. Where the sine actually is on
        // that axis, relative to the CSV clock, is phase error / (2 pi f).
        double A, phi;
        fitSine(y, m.rate, kF, (size_t)m.rate, y.size() - (size_t)m.rate, A, phi);
        const double off = -wrapPi(phi - kPhi) / (kTwoPi * kF);
        offset.push_back(off);
        std::printf("  %s 256 -> %.0f Hz: offset from CSV clock %.6f ms\n", m.name, m.rate, off * 1000.0);
        EXPECT_LE(std::fabs(off), tol) << m.name << " is " << off * 1000.0 << " ms off the CSV timestamps";
    }
    EXPECT_LE(std::fabs(offset[0] - offset[1]), tol)
        << "ecg and ppg are " << (offset[0] - offset[1]) * 1000.0 << " ms apart";
}