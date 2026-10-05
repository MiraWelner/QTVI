// ============================================================================
// test1f_rpeak_consensus.cpp -- QTVI_Tests (Google Test)
//
// ACCEPTANCE TEST: run the same record with the consensus flag on and off.
//   ON:  all six detectors contribute, and the weights sum correctly at every
//        accepted peak.
//   OFF: only the single detector runs, and runtime drops.
//   Detected R peaks against a manual annotation of a clean record:
//        sensitivity must exceed 99%.
//
// WHAT RUNS. The flag is config_entry::use_consensus_rpeak (config.csv
// use_consensus_rpeak). make_beats.hpp's run_rpeak_detection is where the
// pipeline acts on it:
//   on  -> JoinedRR_full: six detectors -- rpeakdetect at thresholds 0.2, 0.1
//          and 0.4, Pan-Tompkins, ecgLms, RRsimpleSquared -- weighted 0.75,
//          0.25, 0.25, 1.25, 1.5, 0.75; detections within 8 samples merge,
//          their weights sum, and a position is accepted at >= 2.4.
//   off -> rpeakdetect at threshold 0.2, alone.
// The tests call run_rpeak_detection itself, so they test the switch the
// pipeline uses, not a copy of it. JoinedRRResult reports each accepted
// peak's votes and summed weight (peak_votes, peak_weight).
//
// FIXTURE: QTVI_Tests\data\rpeaks\3010422_20110321_4h26m00s_20s.bin
//   a data .bin (file_to_bin layout): 20 s of clean MESA 3010422, 4:26:00 -
//   4:26:20 into the recording, cut from the EDF and converted at MESA's
//   rates -- ECG 256 -> 1000 Hz, PPG 256 -> 500 Hz. The test reads ECG1's
//   upsampled block, the signal the pipeline's detectors run on.
// MANUAL ANNOTATION: kManual below -- the 26 R peaks of that ECG, marked by
//   hand from the waveform (one mark at the top of each QRS), not with any of
//   the detectors under test.
//
// OUTPUT: QTVI_Tests\test_output\test1f_rpeaks.png (Qt Gui builds only):
//   the ECG with the manual marks and both detectors' peaks, and the summed
//   weight each accepted peak was voted in on.
// ============================================================================
#include "pch.h"
#include "peak_finding/make_beats.hpp"   // run_rpeak_detection (the flag), and through it
                                         // JoinedRR_full (consensus) and rpeakdetect (single)
#include "config_file_handling/config.hpp"

#ifndef TESTS_DATA_DIR
#error "TESTS_DATA_DIR is not defined. QTVI_Tests -> Properties (All Configurations, All Platforms) -> C/C++ -> Preprocessor -> Preprocessor Definitions: add  TESTS_DATA_DIR=R\"($(ProjectDir)data)\""
#endif

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <string>
#include <vector>

#ifdef QT_GUI_LIB
#include <QColor>
#include <QFont>
#include <QGuiApplication>
#include <QImage>
#include <QPainter>
#include <QPen>
#include <QPointF>
#include <QPolygonF>
#include <QString>
#endif

namespace fs = std::filesystem;

namespace {
    constexpr double kFs = 1000.0;                // the fixture's rate = the pipeline's ECG rate
    constexpr double kMatchMs = 50.0;             // a detection within 50 ms of a mark finds it
    const std::array<double, 6> kWeights = { 0.75, 0.25, 0.25, 1.25, 1.5, 0.75 };
    constexpr double kAccept = 2.4;               // JoinedRR_full accepts at weight >= 2.399
    const char* kDetectorName[6] = { "rpeakdetect 0.2", "rpeakdetect 0.1", "rpeakdetect 0.4",
                                     "Pan-Tompkins", "ecgLms", "RRsimpleSquared" };

    const fs::path kBin = fs::path(TESTS_DATA_DIR) / "rpeaks" / "3010422_20110321_4h26m00s_20s.bin";

    // The 26 R peaks, marked by hand on kBin's ECG1 (sample numbers at 1000 Hz).
    const std::vector<std::size_t> kManual = {
          516,  1273,  2031,  2798,  3578,  4355,  5141,  5931,  6717,  7504,
         8288,  9068,  9845, 10625, 11397, 12154, 12906, 13658, 14418, 15176,
        15926, 16685, 17452, 18219, 18971, 19736 };

    // ECG1's upsampled block of a file_to_bin .bin: header (592 bytes), then
    // channel 0 (timestamp: upsampled doubles + raw pairs), then ECG1.
    std::vector<double> readEcg1(const fs::path& p, double* rateOut = nullptr) {
        std::ifstream f(p, std::ios::binary);
        std::vector<char> b((std::istreambuf_iterator<char>(f)), {});
        if (b.size() < 592) return {};
        auto u32 = [&](size_t o) { uint32_t v; std::memcpy(&v, &b[o], 4); return v; };
        float rate; std::memcpy(&rate, &b[444 + 4], 4);              // up_rates[1] = ECG1
        const size_t off = 592 + 8ull * u32(12) + 16ull * u32(156);  // past channel 0
        const size_t n = u32(12 + 4);                                // sizes_up[1]
        if (off + 8 * n > b.size()) return {};
        std::vector<double> v(n);
        std::memcpy(v.data(), &b[off], 8 * n);
        if (rateOut) *rateOut = rate;
        return v;
    }
    std::vector<double> ecg() { return readEcg1(kBin); }
    const std::vector<std::size_t>& manual() { return kManual; }

    // The pipeline's own switch: make_beats.hpp's run_rpeak_detection.
    std::vector<std::size_t> detect(const std::vector<double>& x, bool consensus) {
        config_entry cfg;
        cfg.use_consensus_rpeak = consensus;
        std::vector<std::size_t> r;
        run_rpeak_detection(x, kFs, "test1f", r, cfg, /*inverted=*/false);
        return r;
    }

    // Sensitivity (manual marks found) and PPV (detections that are real),
    // matching each mark to its nearest unused detection within kMatchMs.
    struct Score { int truth, found, detections, matched; double sens, ppv; };
    Score score(const std::vector<std::size_t>& truth, const std::vector<std::size_t>& det) {
        const long tol = static_cast<long>(kMatchMs * kFs / 1000.0);
        std::vector<bool> used(det.size(), false);
        int found = 0;
        for (std::size_t t : truth) {
            long best = -1, bestD = tol + 1;
            for (std::size_t k = 0; k < det.size(); ++k) {
                const long d = std::labs(static_cast<long>(det[k]) - static_cast<long>(t));
                if (!used[k] && d <= tol && d < bestD) { best = static_cast<long>(k); bestD = d; }
            }
            if (best >= 0) { used[best] = true; ++found; }
        }
        Score s{ (int)truth.size(), found, (int)det.size(), found, 0, 0 };
        s.sens = truth.empty() ? 0 : 100.0 * found / truth.size();
        s.ppv = det.empty() ? 0 : 100.0 * found / det.size();
        return s;
    }

    // Median wall time of `reps` calls, in ms.
    double medianMs(const std::function<void()>& fn, int reps = 15) {
        std::vector<double> t;
        for (int k = 0; k < reps; ++k) {
            const auto a = std::chrono::steady_clock::now();
            fn();
            t.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count());
        }
        std::sort(t.begin(), t.end());
        return t[t.size() / 2];
    }
}

// ============================================================================

TEST(test1f_rpeak_consensus, FixtureLoads) {
    double rate = 0;
    const std::vector<double> x = readEcg1(kBin, &rate);
    ASSERT_FALSE(x.empty()) << "cannot read ECG1 from " << kBin.string();
    EXPECT_EQ(rate, kFs);
    EXPECT_EQ(x.size(), 20000u) << "expected 20 s at 1000 Hz";
    for (std::size_t s : kManual) EXPECT_LT(s, x.size());
}

// ---- consensus ON ------------------------------------------------------------

TEST(test1f_rpeak_consensus, FlagOnRunsTheConsensus) {
    const std::vector<double> x = ecg();
    ASSERT_FALSE(x.empty());
    EXPECT_TRUE(detect(x, true) == JoinedRR_full(x, kFs, "test1f", false).peaks);
}

TEST(test1f_rpeak_consensus, ConsensusOnAllSixDetectorsContribute) {
    const std::vector<double> x = ecg();
    ASSERT_FALSE(x.empty());
    const JoinedRRResult j = JoinedRR_full(x, kFs, "test1f", false);
    ASSERT_FALSE(j.peaks.empty());
    ASSERT_EQ(j.peak_votes.size(), j.peaks.size());
    for (int d = 0; d < 6; ++d) {
        SCOPED_TRACE(kDetectorName[d]);
        EXPECT_FALSE(j.detector_peaks[d].empty()) << "detected nothing";
        int contributes = 0;
        for (const auto& v : j.peak_votes) contributes += (v[d] > 0);
        // On a clean record every detector should agree on (nearly) every beat.
        EXPECT_GE(contributes, static_cast<int>(0.95 * j.peaks.size()))
            << "voted for only " << contributes << " of " << j.peaks.size() << " accepted peaks";
    }
}

TEST(test1f_rpeak_consensus, ConsensusOnWeightsSumCorrectlyAtEveryAcceptedPeak) {
    const std::vector<double> x = ecg();
    ASSERT_FALSE(x.empty());
    const JoinedRRResult j = JoinedRR_full(x, kFs, "test1f", false);
    for (int d = 0; d < 6; ++d) EXPECT_EQ(j.detector_weights[d], kWeights[d]) << kDetectorName[d];
    ASSERT_EQ(j.peak_weight.size(), j.peaks.size());
    for (std::size_t k = 0; k < j.peaks.size(); ++k) {
        SCOPED_TRACE("accepted peak at sample " + std::to_string(j.peaks[k]));
        double sum = 0;
        for (int d = 0; d < 6; ++d) sum += j.peak_votes[k][d] * kWeights[d];
        EXPECT_NEAR(j.peak_weight[k], sum, 1e-12) << "the weight does not equal its votes";
        EXPECT_GE(j.peak_weight[k], kAccept - 1e-3) << "accepted below the 2.4 threshold";
    }
}

// ---- consensus OFF -----------------------------------------------------------

TEST(test1f_rpeak_consensus, FlagOffRunsOnlyTheSingleDetector) {
    const std::vector<double> x = ecg();
    ASSERT_FALSE(x.empty());
    EXPECT_TRUE(detect(x, false) == rpeakdetect(x, kFs, 0.2, false).r_peak_index);
}

TEST(test1f_rpeak_consensus, FlagOffIsFaster) {
    const std::vector<double> x = ecg();
    ASSERT_FALSE(x.empty());
    const double on = medianMs([&] { (void)detect(x, true); });
    const double off = medianMs([&] { (void)detect(x, false); });
    // One detector out of six: well under the consensus's time, not just under.
    EXPECT_LT(off, 0.75 * on) << "consensus " << on << " ms, single " << off << " ms";
}

// ---- against the manual annotation -------------------------------------------

TEST(test1f_rpeak_consensus, ConsensusSensitivityAbove99Percent) {
    const Score s = score(manual(), detect(ecg(), true));
    EXPECT_GT(s.sens, 99.0) << s.found << " of " << s.truth << " manual peaks found (PPV " << s.ppv << "%)";
}

TEST(test1f_rpeak_consensus, SingleDetectorSensitivityAbove99Percent) {
    const Score s = score(manual(), detect(ecg(), false));
    EXPECT_GT(s.sens, 99.0) << s.found << " of " << s.truth << " manual peaks found (PPV " << s.ppv << "%)";
}

// ============================================================================
// THE FIGURE (Qt Gui builds only) -> QTVI_Tests\test_output\test1f_rpeaks.png
// ============================================================================
#ifdef QT_GUI_LIB
namespace {
    void ensureGuiApp() {   // text needs a QGuiApplication; none is shown
        if (QCoreApplication::instance()) return;
        static int argc = 1;
        static char name[] = "qtvi_tests";
        static char* argv[] = { name, nullptr };
        static QGuiApplication app(argc, argv);
    }
    struct Frame {
        QRectF plot; double x0, x1, y0, y1;
        double px(double x) const { return plot.left() + (x - x0) / (x1 - x0) * plot.width(); }
        double py(double y) const { return plot.bottom() - (y - y0) / (y1 - y0) * plot.height(); }
    };
    Frame axes(QPainter& p, const QRectF& r, const QString& title, double x0, double x1, double y0, double y1,
        const QString& xl, const QString& yl) {
        Frame f{ r.adjusted(60, 28, -14, -34), x0, x1, y0, y1 };
        p.setPen(QPen(Qt::black, 1)); p.drawRect(f.plot);
        QFont font = p.font(); font.setPointSize(10); font.setBold(true); p.setFont(font);
        p.drawText(QRectF(r.left(), r.top() + 2, r.width(), 22), Qt::AlignHCenter, title);
        font.setBold(false); font.setPointSize(8); p.setFont(font);
        for (int k = 0; k <= 10; ++k) {
            const double xv = x0 + (x1 - x0) * k / 10;
            p.drawText(QPointF(f.px(xv) - 8, f.plot.bottom() + 14), QString::number(xv, 'f', 0));
        }
        p.drawText(QRectF(f.plot.left(), f.plot.bottom() + 16, f.plot.width(), 16), Qt::AlignHCenter, xl);
        p.drawText(QPointF(r.left() + 4, f.plot.top() + 10), QString::number(y1, 'f', 2));
        p.drawText(QPointF(r.left() + 4, f.plot.bottom()), QString::number(y0, 'f', 2));
        p.save(); p.translate(r.left() + 14, f.plot.center().y()); p.rotate(-90);
        p.drawText(QRectF(-80, -10, 160, 16), Qt::AlignHCenter, yl); p.restore();
        return f;
    }
}

TEST(test1f_rpeak_consensus, WritesFigure) {
    ensureGuiApp();
    const std::vector<double> x = ecg();
    ASSERT_FALSE(x.empty());
    const std::vector<std::size_t> m = manual(), single = detect(x, false);
    const JoinedRRResult j = JoinedRR_full(x, kFs, "test1f", false);
    const Score sc = score(m, j.peaks), ss = score(m, single);

    QImage img(1800, 760, QImage::Format_RGB32); img.fill(Qt::white);
    QPainter p(&img); p.setRenderHint(QPainter::Antialiasing);

    // ECG with the three sets of peaks
    double lo = *std::min_element(x.begin(), x.end()), hi = *std::max_element(x.begin(), x.end());
    lo -= 0.1; hi += 0.35;
    const Frame f = axes(p, QRectF(0, 0, 1800, 470),
        QString("3010422, 4:26:00 + 20 s: manual %1 | consensus found %2 (sens %3%) | single found %4 (sens %5%)")
        .arg(m.size()).arg(sc.found).arg(sc.sens, 0, 'f', 1).arg(ss.found).arg(ss.sens, 0, 'f', 1),
        0, 20, lo, hi, "s", "ECG (mV)");
    QPolygonF trace;
    for (std::size_t i = 0; i < x.size(); i += 2) trace << QPointF(f.px(i / kFs), f.py(x[i]));
    p.setClipRect(f.plot);
    p.setPen(QPen(Qt::black, 1)); p.drawPolyline(trace);
    p.setPen(QPen(QColor(214, 39, 40), 2));       // manual: red circles
    for (std::size_t s : m) p.drawEllipse(QPointF(f.px(s / kFs), f.py(x[s])), 7, 7);
    p.setPen(QPen(QColor(31, 119, 180), 2));      // consensus: blue x above
    for (std::size_t s : j.peaks) {
        const QPointF c(f.px(s / kFs), f.py(x[s]) - 18);
        p.drawLine(c + QPointF(-5, -5), c + QPointF(5, 5)); p.drawLine(c + QPointF(-5, 5), c + QPointF(5, -5));
    }
    p.setPen(QPen(QColor(44, 160, 44), 2));       // single: green + higher
    for (std::size_t s : single) {
        const QPointF c(f.px(s / kFs), f.py(x[s]) - 34);
        p.drawLine(c + QPointF(-6, 0), c + QPointF(6, 0)); p.drawLine(c + QPointF(0, -6), c + QPointF(0, 6));
    }
    p.setClipping(false);
    QFont font = p.font(); font.setPointSize(9); p.setFont(font);
    p.setPen(QColor(214, 39, 40)); p.drawText(QPointF(1460, 50), "o  manual");
    p.setPen(QColor(31, 119, 180)); p.drawText(QPointF(1460, 66), "x  consensus (flag on)");
    p.setPen(QColor(44, 160, 44)); p.drawText(QPointF(1460, 82), "+  single detector (flag off)");

    // the summed weight of each accepted peak, against the 2.4 threshold
    const Frame g = axes(p, QRectF(0, 480, 1800, 280), "consensus: summed detector weight of each accepted peak",
        0, 20, 0, 5, "s", "weight");
    p.setPen(QPen(QColor(150, 150, 150), 1, Qt::DashLine));
    p.drawLine(QPointF(g.plot.left(), g.py(kAccept)), QPointF(g.plot.right(), g.py(kAccept)));
    p.drawText(QPointF(g.plot.right() - 120, g.py(kAccept) - 4), "accept at 2.4");
    p.drawLine(QPointF(g.plot.left(), g.py(4.75)), QPointF(g.plot.right(), g.py(4.75)));
    p.drawText(QPointF(g.plot.right() - 160, g.py(4.75) - 4), "all six agree: 4.75");
    // one bar per accepted peak: its summed weight, and how many of the six voted
    font.setPointSize(8); p.setFont(font);
    for (std::size_t k = 0; k < j.peaks.size(); ++k) {
        const double xc = g.px(j.peaks[k] / kFs);
        p.fillRect(QRectF(QPointF(xc - 9, g.py(j.peak_weight[k])), QPointF(xc + 9, g.py(0))), QColor(31, 119, 180));
        int voters = 0;
        for (int d = 0; d < 6; ++d) voters += (j.peak_votes[k][d] > 0);
        p.setPen(Qt::black);
        p.drawText(QRectF(xc - 20, g.py(j.peak_weight[k]) - 16, 40, 14), Qt::AlignHCenter, QString("%1/6").arg(voters));
    }
    p.end();

    const fs::path out = fs::path(TESTS_DATA_DIR).parent_path() / "test_output";
    fs::create_directories(out);
    const fs::path png = out / "test1f_rpeaks.png";
    EXPECT_TRUE(img.save(QString::fromStdWString(png.wstring()))) << "could not write " << png.string();
}
#endif   // QT_GUI_LIB