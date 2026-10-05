// ============================================================================
// test1g_ppg_peaks_pairing.cpp -- QTVI_Tests (Google Test)
//
// ACCEPTANCE TEST: on a clean PPG record, detected peaks match a manual
// annotation with sensitivity > 95%, and pairing assigns each PPG pulse to the
// correct preceding R peak.
//
// WHAT IS TESTED -- the pipeline's own two PPG paths:
//   SegmentPPG (SegmentPPG.hpp)        segments the whole PPG into pulses;
//       its peaks (maxAmps) are the PPG peak columns of the peak-locations CSV.
//   build_pulse_template_pair_windowed (bin_pulse.hpp)  the pairing: for each
//       R pair (R_i, R_i+1) it cuts one pulse and finds its peak only in
//       [R_i, R_i+1], so the pulse belongs to R_i. Each ptt record says which
//       R pair (slice) and where its pulse peak is (peak_s) -- the records the
//       transit-time columns of the peak-locations CSV are made of.
//
// FIXTURE: QTVI_Tests\data\rpeaks\3010422_20110321_4h26m00s_20s.bin -- the
//   same 20 s as test1f (clean MESA 3010422, 4:26:00 - 4:26:20). The PPG is
//   its slot 4 at 500 Hz.
// MANUAL ANNOTATION, both marked by hand from the waveforms:
//   kManualR    the 26 R peaks (ECG samples at 1000 Hz) -- as in test1f
//   kManualPpg  the 26 PPG systolic peaks (PPG samples at 500 Hz)
//   Pulse 1 (0.16 s) belongs to an R before the recording starts, and the
//   last R's pulse falls after it ends: so 25 R pairs, paired with pulses 2-26.
//
// OUTPUT: QTVI_Tests\test_output\test1g_ppg_pairing.png (Qt Gui builds only).
//
// PROJECT SETTINGS: as the other tests. bin_pulse calls into feature_marks.cpp,
// which needs ppg_dicrotic.cpp and ppg_pipeline.cpp: add link_feature_marks.cpp,
// link_ppg_dicrotic.cpp and link_ppg_pipeline.cpp to the test project.
// ============================================================================
#include "pch.h"
#include "peak_finding/SegmentPPG.hpp"          // SegmentPPG
#include "template_generation/bin_pulse.hpp"     // build_pulse_template_pair_windowed

#ifndef TESTS_DATA_DIR
#error "TESTS_DATA_DIR is not defined. QTVI_Tests -> Properties (All Configurations, All Platforms) -> C/C++ -> Preprocessor -> Preprocessor Definitions: add  TESTS_DATA_DIR=R\"($(ProjectDir)data)\""
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>
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
    constexpr double kEcgHz = 1000.0, kPpgHz = 500.0;
    constexpr double kMatchMs = 50.0;   // a detection within 50 ms of a mark finds it
    const fs::path kBin = fs::path(TESTS_DATA_DIR) / "rpeaks" / "3010422_20110321_4h26m00s_20s.bin";

    // R peaks, marked by hand on the ECG (samples at 1000 Hz).
    const std::vector<std::size_t> kManualR = {
          516,  1273,  2031,  2798,  3578,  4355,  5141,  5931,  6717,  7504,
         8288,  9068,  9845, 10625, 11397, 12154, 12906, 13658, 14418, 15176,
        15926, 16685, 17452, 18219, 18971, 19736 };
    // PPG systolic peaks, marked by hand on the PPG (samples at 500 Hz).
    const std::vector<std::size_t> kManualPpg = {
           81,   460,   837,  1217,  1603,  1987,  2377,  2772,  3158,  3549,
         3951,  4360,  4727,  5143,  5532,  5909,  6289,  6677,  7053,  7429,
         7797,  8185,  8559,  8938,  9310,  9698 };

    // One channel's upsampled block of a file_to_bin .bin.
    std::vector<double> readSlot(const fs::path& p, int slot, double* rateOut = nullptr) {
        std::ifstream f(p, std::ios::binary);
        std::vector<char> b((std::istreambuf_iterator<char>(f)), {});
        if (b.size() < 592) return {};
        auto u32 = [&](size_t o) { uint32_t v; std::memcpy(&v, &b[o], 4); return v; };
        size_t off = 592;
        for (int c = 0; c < slot; ++c) off += 8ull * u32(12 + 4 * c) + 16ull * u32(156 + 4 * c);
        const size_t n = u32(12 + 4 * slot);
        if (off + 8 * n > b.size()) return {};
        std::vector<double> v(n);
        std::memcpy(v.data(), &b[off], 8 * n);
        if (rateOut) { float r; std::memcpy(&r, &b[444 + 4 * slot], 4); *rateOut = r; }
        return v;
    }
    std::vector<double> ppg() { return readSlot(kBin, 4); }

    // The manual R peak a pulse at PPG sample p follows: the last R before it,
    // or -1 when the recording holds none.
    int precedingR(double ppgSample) {
        const double t = ppgSample / kPpgHz;
        int idx = -1;
        for (int i = 0; i < (int)kManualR.size(); ++i)
            if (kManualR[i] / kEcgHz < t) idx = i;
        return idx;
    }
    // The manual pulse within kMatchMs of a PPG sample, or -1.
    int manualPulseAt(double ppgSample) {
        const double tol = kMatchMs * kPpgHz / 1000.0;
        int best = -1; double bd = tol + 1;
        for (int j = 0; j < (int)kManualPpg.size(); ++j) {
            const double d = std::abs(ppgSample - (double)kManualPpg[j]);
            if (d <= tol && d < bd) { best = j; bd = d; }
        }
        return best;
    }
    // Sensitivity: of `truth`, how many have a distinct detection within kMatchMs.
    struct Score { int truth = 0, found = 0, detections = 0; double sens = 0, ppv = 0; };
    Score score(const std::vector<std::size_t>& truth, const std::vector<double>& det) {
        const double tol = kMatchMs * kPpgHz / 1000.0;
        std::vector<bool> used(det.size(), false);
        Score s; s.truth = (int)truth.size(); s.detections = (int)det.size();
        for (std::size_t t : truth) {
            int best = -1; double bd = tol + 1;
            for (std::size_t k = 0; k < det.size(); ++k) {
                const double d = std::abs(det[k] - (double)t);
                if (!used[k] && d <= tol && d < bd) { best = (int)k; bd = d; }
            }
            if (best >= 0) { used[best] = true; ++s.found; }
        }
        s.sens = truth.empty() ? 0 : 100.0 * s.found / truth.size();
        s.ppv = det.empty() ? 0 : 100.0 * s.found / det.size();
        return s;
    }
    // The pipeline's pairing, run on the manual R peaks.
    std::vector<ptt_log::Beat> pairs(const std::vector<double>& x) {
        return build_pulse_template_pair_windowed(x, kPpgHz, kManualR, kEcgHz).ptt;
    }
}

// ============================================================================

TEST(test1g_ppg_peaks_pairing, FixtureLoads) {
    double rate = 0;
    const std::vector<double> x = readSlot(kBin, 4, &rate);
    ASSERT_FALSE(x.empty()) << "cannot read the PPG (slot 4) from " << kBin.string();
    EXPECT_EQ(rate, kPpgHz);
    EXPECT_EQ(x.size(), 10000u) << "expected 20 s at 500 Hz";
    for (std::size_t s : kManualPpg) EXPECT_LT(s, x.size());
}

// ---- detected peaks against the manual annotation ----------------------------

TEST(test1g_ppg_peaks_pairing, SegmentPPGSensitivityAbove95Percent) {
    const std::vector<double> x = ppg();
    ASSERT_FALSE(x.empty());
    const SegmentPPGResult seg = SegmentPPG(x, kPpgHz);
    const std::vector<double> det(seg.maxAmps.begin(), seg.maxAmps.end());
    const Score s = score(kManualPpg, det);
    EXPECT_GT(s.sens, 95.0) << s.found << " of " << s.truth << " manual peaks found (PPV " << s.ppv << "%)";
}

TEST(test1g_ppg_peaks_pairing, PairedPulsePeaksSensitivityAbove95Percent) {
    const std::vector<double> x = ppg();
    ASSERT_FALSE(x.empty());
    std::vector<double> det;
    for (const ptt_log::Beat& b : pairs(x)) det.push_back(b.peak_s * kPpgHz);
    // Only pulses with an R before them can be paired: pulse 1's R predates the recording.
    std::vector<std::size_t> pairable;
    for (std::size_t p : kManualPpg) if (precedingR((double)p) >= 0) pairable.push_back(p);
    const Score s = score(pairable, det);
    EXPECT_GT(s.sens, 95.0) << s.found << " of " << s.truth << " pairable pulses found (PPV " << s.ppv << "%)";
}

// ---- pairing ------------------------------------------------------------------

TEST(test1g_ppg_peaks_pairing, EachPulseIsPairedWithItsPrecedingR) {
    const std::vector<double> x = ppg();
    ASSERT_FALSE(x.empty());
    const std::vector<ptt_log::Beat> beats = pairs(x);
    ASSERT_FALSE(beats.empty()) << "the pairing produced no beats";
    std::set<int> claimed;
    for (const ptt_log::Beat& b : beats) {
        SCOPED_TRACE("R pair " + std::to_string(b.slice));
        ASSERT_LT(b.slice, kManualR.size());
        EXPECT_NEAR(b.r_peak_distance_from_binstart_in_s, kManualR[b.slice] / kEcgHz, 1e-9)
            << "the record names a different R than its slice";
        const int j = manualPulseAt(b.peak_s * kPpgHz);
        ASSERT_GE(j, 0) << "paired a pulse peak at " << b.peak_s << " s that is not a manual pulse";
        EXPECT_EQ(precedingR((double)kManualPpg[j]), (int)b.slice)
            << "pulse " << (j + 1) << " (" << kManualPpg[j] / kPpgHz << " s) follows R "
            << precedingR((double)kManualPpg[j]) << ", but was paired with R " << b.slice;
        EXPECT_TRUE(claimed.insert(j).second) << "pulse " << (j + 1) << " was paired twice";
    }
}

TEST(test1g_ppg_peaks_pairing, EveryPairableRGetsItsPulse) {
    const std::vector<double> x = ppg();
    ASSERT_FALSE(x.empty());
    std::set<uint32_t> paired;
    for (const ptt_log::Beat& b : pairs(x)) paired.insert(b.slice);
    int pairable = 0, got = 0;
    for (int i = 0; i < (int)kManualR.size(); ++i) {
        bool hasPulse = false;   // a manual pulse whose preceding R is this one
        for (std::size_t p : kManualPpg) hasPulse = hasPulse || (precedingR((double)p) == i);
        if (!hasPulse) continue;
        ++pairable;
        if (paired.count((uint32_t)i)) ++got; else ADD_FAILURE() << "R " << i << " got no pulse";
    }
    EXPECT_GT(100.0 * got / std::max(1, pairable), 95.0) << got << " of " << pairable << " R peaks paired";
}

TEST(test1g_ppg_peaks_pairing, PulseBeforeTheFirstRIsNotPaired) {
    const std::vector<double> x = ppg();
    ASSERT_FALSE(x.empty());
    ASSERT_EQ(precedingR((double)kManualPpg.front()), -1);   // pulse 1 predates every R here
    for (const ptt_log::Beat& b : pairs(x))
        EXPECT_NE(manualPulseAt(b.peak_s * kPpgHz), 0) << "pulse 1 was paired with R " << b.slice;
}

// ============================================================================
// THE FIGURE (Qt Gui builds only) -> QTVI_Tests\test_output\test1g_ppg_pairing.png
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
    Frame axes(QPainter& p, const QRectF& r, const QString& title, double y0, double y1, const QString& yl) {
        Frame f{ r.adjusted(60, 28, -14, -34), 0, 20, y0, y1 };
        p.setPen(QPen(Qt::black, 1)); p.drawRect(f.plot);
        QFont font = p.font(); font.setPointSize(10); font.setBold(true); p.setFont(font);
        p.drawText(QRectF(r.left(), r.top() + 2, r.width(), 22), Qt::AlignHCenter, title);
        font.setBold(false); font.setPointSize(8); p.setFont(font);
        for (int k = 0; k <= 10; ++k) p.drawText(QPointF(f.px(2.0 * k) - 6, f.plot.bottom() + 14), QString::number(2 * k));
        p.drawText(QRectF(f.plot.left(), f.plot.bottom() + 16, f.plot.width(), 16), Qt::AlignHCenter, "s");
        p.save(); p.translate(r.left() + 14, f.plot.center().y()); p.rotate(-90);
        p.drawText(QRectF(-80, -10, 160, 16), Qt::AlignHCenter, yl); p.restore();
        return f;
    }
    void trace(QPainter& p, const Frame& f, const std::vector<double>& v, double hz) {
        QPolygonF poly;
        for (std::size_t i = 0; i < v.size(); i += 2) poly << QPointF(f.px(i / hz), f.py(v[i]));
        p.setPen(QPen(Qt::black, 1)); p.drawPolyline(poly);
    }
}

TEST(test1g_ppg_peaks_pairing, WritesFigure) {
    ensureGuiApp();
    const std::vector<double> e = readSlot(kBin, 1), x = ppg();
    ASSERT_FALSE(e.empty());
    ASSERT_FALSE(x.empty());
    const SegmentPPGResult seg = SegmentPPG(x, kPpgHz);
    const std::vector<ptt_log::Beat> beats = pairs(x);
    const Score sc = score(kManualPpg, std::vector<double>(seg.maxAmps.begin(), seg.maxAmps.end()));

    QImage img(1800, 820, QImage::Format_RGB32); img.fill(Qt::white);
    QPainter p(&img); p.setRenderHint(QPainter::Antialiasing);
    const QColor red(214, 39, 40), blue(31, 119, 180), grey(130, 130, 130);

    const double eLo = *std::min_element(e.begin(), e.end()) - 0.1, eHi = *std::max_element(e.begin(), e.end()) + 0.2;
    const Frame fe = axes(p, QRectF(0, 0, 1800, 300), "ECG: R peaks marked by hand", eLo, eHi, "ECG (mV)");
    p.setClipRect(fe.plot); trace(p, fe, e, kEcgHz);
    p.setPen(QPen(red, 2));
    for (std::size_t r : kManualR) p.drawEllipse(QPointF(fe.px(r / kEcgHz), fe.py(e[r])), 6, 6);
    p.setClipping(false);

    const double xLo = *std::min_element(x.begin(), x.end()) - 0.05, xHi = *std::max_element(x.begin(), x.end()) + 0.12;
    const Frame fp = axes(p, QRectF(0, 310, 1800, 510),
        QString("PPG: manual peaks (o) %1 | SegmentPPG (x) found %2 (sens %3%) | each line joins a pulse to the R it was paired with")
        .arg(kManualPpg.size()).arg(sc.found).arg(sc.sens, 0, 'f', 1), xLo, xHi, "PPG");
    p.setClipRect(fp.plot);
    p.setPen(QPen(grey, 1, Qt::DotLine));
    for (std::size_t r : kManualR) p.drawLine(QPointF(fp.px(r / kEcgHz), fp.plot.top()), QPointF(fp.px(r / kEcgHz), fp.plot.bottom()));
    trace(p, fp, x, kPpgHz);
    p.setPen(QPen(blue, 2));   // pairing: from the R line, at the top, to the paired pulse peak
    for (const ptt_log::Beat& b : beats) {
        const double pk = b.peak_s * kPpgHz;
        const int k = std::clamp((int)std::lround(pk), 0, (int)x.size() - 1);
        p.drawLine(QPointF(fp.px(b.r_peak_distance_from_binstart_in_s), fp.py(xHi - 0.02)),
            QPointF(fp.px(b.peak_s), fp.py(x[k]) - 8));
    }
    p.setPen(QPen(red, 2));
    for (std::size_t s : kManualPpg) p.drawEllipse(QPointF(fp.px(s / kPpgHz), fp.py(x[s])), 7, 7);
    p.setPen(QPen(Qt::black, 2));
    for (std::size_t s : seg.maxAmps) {
        const QPointF c(fp.px(s / kPpgHz), fp.py(x[s]) - 16);
        p.drawLine(c + QPointF(-5, -5), c + QPointF(5, 5)); p.drawLine(c + QPointF(-5, 5), c + QPointF(5, -5));
    }
    p.setClipping(false);
    p.end();

    const fs::path out = fs::path(TESTS_DATA_DIR).parent_path() / "test_output";
    fs::create_directories(out);
    const fs::path png = out / "test1g_ppg_pairing.png";
    EXPECT_TRUE(img.save(QString::fromStdWString(png.wstring()))) << "could not write " << png.string();
}
#endif   // QT_GUI_LIB