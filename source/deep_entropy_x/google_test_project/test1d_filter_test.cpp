// ============================================================================
// test1b_noise_filters.cpp -- QTVI_Tests (Google Test)
//
// ACCEPTANCE TEST: mark several noise segments on a multi-channel record and
// verify the binary and CSV outputs match the marked regions. Toggle the high
// pass and the notch filter and confirm baseline wander / powerline noise is
// removed when on and restored when off. Writes the figures to look at, too.
//
// FIXTURE: QTVI_Tests\data\filters\
//   3010422_20110321_4h25m-4h35m_clean.bin  10 min of a MESA record (ECG 1000 Hz,
//       PPG 500 Hz), cut 4 h 25 min - 4 h 35 min into it: 3 min of clean
//       beats, then movement with several mV of baseline wander
//   3010422_20110321_4h25m-4h35m_60hz.bin   the same, plus a steady 0.10 mV
//       60 Hz sine added to the ECG only
//
// HOW THE FILTERS ARE APPLIED. In the noise-marking GUI both filters are
// display-only: when a checkbox is on, signal_renderer.cpp's filteredSpan()
// filters the visible window -- padded by max(5 s, 10 / cutoff) each side --
// high pass first, then notch, with FilterUtils' waveform_highpass and
// notch_filter. viewAsGui() below does exactly that over 30 s windows. "Off"
// means a cutoff of 0, which filteredSpan answers by drawing the original.
// Cutoffs are MESA's in config.csv: waveform_highpass_hz 0.5, notch_filter_hz 60.
//
// OUTPUT: two PNGs in QTVI_Tests\test_output\ (built with Qt Gui only), each
// OFF left and ON right, the signal on top and its spectrum below:
//   test1b_highpass.png   the original file's movement stretch, high pass
//   test1b_notch.png      the 60 Hz file, high pass on, notch
//
// PROJECT SETTINGS: as test1a, plus user_annotation_handler.cpp in the test
// project (Add -> Existing Item, from noise_marking_gui\) -- the marking tests
// call the GUI's writer, annotation_handler.
// ============================================================================
#include "pch.h"
#include "../noise_marking_gui/user_annotation_handler.h"   // annotation_handler (writer), noise_markings::loadRows (reader)
#include "../noise_marking_gui/annotation_types.hpp"        // label -> code
#include "../peak_finding/FilterUtils.hpp"                  // waveform_highpass, notch_filter

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
#include <sstream>
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
namespace nm = noise_markings;

static const double kPi = 3.14159265358979323846;

// ============================================================================
// PART 1 -- MARKED REGIONS INTO THE .bin AND THE CSV
// ============================================================================

// The regions an operator marks on the fixture record, in seconds into it.
// Several channels and every kind of span the writer stores: plain noise,
// a parameter edit (threshold / blanking), an arrhythmia mark, and one
// dragged right-to-left (end before start), which must come out ordered.
struct Mark { const char* channel; const char* type; double start, end; double thr, blk; };
static const Mark kMarks[] = {
    { "ECG1", "1) R Peak Noise",    185.000, 200.000, NAN, NAN },   // the clipped stretch
    { "ECG1", "2) Minor Noise",     240.000, 270.000, NAN, NAN },
    { "PPG",  "2) Minor Noise",     250.000, 262.500, NAN, NAN },
    { "PPG",  "4) PVC",              20.004,  20.900, NAN, NAN },
    { "ECG1", "3) Blank.+Thresh.",  300.000, 320.000, 0.6, 400.0 },
    { "ECG1", "1) R Peak Noise",    512.250, 505.000, NAN, NAN },   // dragged right-to-left
};
static double rateFor(const std::string& ch) { return ch == "PPG" ? 500.0 : 1000.0; }   // MESA upsample rates

// What the GUI's Finish does with them (main.cpp, exportMarkings): seconds to
// samples at the channel's rate with llround, then the CSV and the .bin from
// one annotation_handler.
static fs::path writeMarkings() {
    const fs::path dir = fs::temp_directory_path() / "qtvi_test1b";
    fs::create_directories(dir);
    annotation_handler h;
    for (const Mark& m : kMarks) {
        const double sr = rateFor(m.channel);
        h.addSegment(static_cast<int>(std::llround(m.start * sr)),
            static_cast<int>(std::llround(m.end * sr)), m.channel, m.type, sr, m.thr, m.blk);
    }
    const fs::path base = dir / "test1b_noise_markings";
    h.exportCSV(base.string() + ".csv");
    h.export_marking_binfile(base.string() + ".bin");
    return base;
}

// The region a mark SHOULD produce: ordered, in samples, and back in seconds.
struct Want { int channel, code; double s0, s1, t0, t1, thr, blk; };
static Want wanted(const Mark& m) {
    const double sr = rateFor(m.channel);
    const double a = std::llround(std::min(m.start, m.end) * sr), b = std::llround(std::max(m.start, m.end) * sr);
    return { nm::code_for_channel(m.channel), annotation_types::code_for_label(m.type), a, b, a / sr, b / sr, m.thr, m.blk };
}
static bool sameValue(double a, double b, double tol) { return (std::isnan(a) && std::isnan(b)) || std::abs(a - b) <= tol; }

TEST(test1b_noise_filters, MarkedRegionsWrittenToBin) {
    const fs::path base = writeMarkings();
    const nm::RowsResult r = nm::loadRows(base.string() + ".bin");
    ASSERT_TRUE(r.read) << r.error;
    ASSERT_EQ(r.rows.size(), std::size(kMarks));
    for (size_t i = 0; i < std::size(kMarks); ++i) {
        SCOPED_TRACE(std::string("mark ") + std::to_string(i + 1) + ": " + kMarks[i].channel + " " + kMarks[i].type);
        const Want w = wanted(kMarks[i]);
        const nm::Row& b = r.rows[i];
        EXPECT_EQ(int(b.channel_code), w.channel);
        EXPECT_EQ(int(b.annotation_code), w.code);
        EXPECT_EQ(b.start_sample, w.s0);
        EXPECT_EQ(b.end_sample, w.s1);
        EXPECT_NEAR(b.start_sec, w.t0, 1e-9);
        EXPECT_NEAR(b.end_sec, w.t1, 1e-9);
        // and the stored region is the marked one, to within half a sample
        EXPECT_NEAR(b.start_sec, std::min(kMarks[i].start, kMarks[i].end), 0.5 / rateFor(kMarks[i].channel));
        EXPECT_NEAR(b.end_sec, std::max(kMarks[i].start, kMarks[i].end), 0.5 / rateFor(kMarks[i].channel));
        EXPECT_TRUE(sameValue(b.threshold, w.thr, 1e-12)) << "threshold " << b.threshold;
        EXPECT_TRUE(sameValue(b.blanking_ms, w.blk, 1e-12)) << "blanking_ms " << b.blanking_ms;
    }
}

TEST(test1b_noise_filters, MarkedRegionsWrittenToCsv) {
    const fs::path base = writeMarkings();
    std::ifstream f(base.string() + ".csv");
    ASSERT_TRUE(f.good()) << "no CSV written";
    std::string line;
    std::getline(f, line);   // header
    size_t i = 0;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        ASSERT_LT(i, std::size(kMarks)) << "the CSV has more rows than marks";
        SCOPED_TRACE(std::string("row ") + std::to_string(i + 1) + ": " + kMarks[i].channel + " " + kMarks[i].type);
        std::vector<std::string> c;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) c.push_back(cell);
        while (c.size() < 8) c.push_back("");
        auto num = [](const std::string& s) { return s.empty() ? std::nan("") : std::stod(s); };
        const Want w = wanted(kMarks[i]);
        EXPECT_EQ(num(c[0]), w.s0);
        EXPECT_EQ(num(c[1]), w.s1);
        EXPECT_NEAR(num(c[2]), w.t0, 0.0005);   // the CSV prints 3 decimals
        EXPECT_NEAR(num(c[3]), w.t1, 0.0005);
        EXPECT_EQ(c[4], std::string(kMarks[i].channel));
        EXPECT_EQ(c[5], std::string(kMarks[i].type));
        EXPECT_TRUE(sameValue(num(c[6]), w.thr, 0.0005)) << "threshold '" << c[6] << "'";
        EXPECT_TRUE(sameValue(num(c[7]), w.blk, 0.0005)) << "blanking_ms '" << c[7] << "'";
        ++i;
    }
    EXPECT_EQ(i, std::size(kMarks)) << "the CSV has fewer rows than marks";
}

// ============================================================================
// PART 2 -- THE FILTERS, AS THE GUI DRAWS THEM
// ============================================================================

// ECG1's upsampled block of a file_to_bin .bin, and its rate from the header.
struct Ecg { std::vector<double> v; double fs = 0.0; };
static Ecg readEcg1(const fs::path& p) {
    Ecg out;
    std::ifstream f(p, std::ios::binary);
    if (!f) return out;
    std::vector<char> b((std::istreambuf_iterator<char>(f)), {});
    if (b.size() < 592) return out;
    auto u32 = [&](size_t o) { uint32_t v; std::memcpy(&v, &b[o], 4); return v; };
    float rate; std::memcpy(&rate, &b[444 + 4], 4);              // up_rates[1] = ECG1
    const size_t off = 592 + 8ull * u32(12) + 16ull * u32(156);  // past channel 0 (timestamp)
    const size_t n = u32(12 + 4);                                // sizes_up[1]
    if (off + 8 * n > b.size()) return out;
    out.v.resize(n);
    std::memcpy(out.v.data(), &b[off], 8 * n);
    out.fs = rate;
    return out;
}

static const fs::path kCleanBin = fs::path(TESTS_DATA_DIR) / "filters" / "3010422_20110321_4h25m-4h35m_baseline_test.bin";
static const fs::path kNoisyBin = fs::path(TESTS_DATA_DIR) / "filters" / "3010422_20110321_4h25m-4h35m_notch_test.bin";
static constexpr double kHighPassHz = 0.5, kNotchHz = 60.0, kWindowSec = 30.0;
// What the fixture HAS, measured whatever the filter is set to: the notch is
// the thing under test, so it must not also decide where we look.
static constexpr double kInjectedHz = 60.0, kInjectedMv = 0.10;

// signal_renderer.cpp's filteredSpan + the window it is called on, over the
// whole signal in kWindowSec windows. Cutoff 0 = that checkbox off.
static std::vector<double> viewAsGui(const Ecg& e, double highPassHz, double notchHz) {
    if (highPassHz <= 0.0 && notchHz <= 0.0) return e.v;      // filteredSpan: draw the original
    const int win = static_cast<int>(kWindowSec * e.fs);
    const double padSec = (highPassHz > 0.0) ? std::max(5.0, 10.0 / highPassHz) : 5.0;
    const int pad = static_cast<int>(padSec * e.fs);
    std::vector<double> out;
    out.reserve(e.v.size());
    for (int from = 0; from < (int)e.v.size(); from += win) {
        const int to = std::min((int)e.v.size(), from + win);
        const int lo = std::max(0, from - pad), hi = std::min((int)e.v.size(), to + pad);
        std::vector<double> buf(e.v.begin() + lo, e.v.begin() + hi);
        if (highPassHz > 0.0) buf = waveform_highpass(buf, highPassHz, e.fs);
        if (notchHz > 0.0) buf = notch_filter(buf, notchHz, e.fs);
        out.insert(out.end(), buf.begin() + (from - lo), buf.begin() + (from - lo) + (to - from));
    }
    return out;
}

// Amplitude of the f0 component, by least squares onto sin / cos at f0.
static double toneAmplitude(const std::vector<double>& x, double fs, double f0) {
    double ss = 0, sc = 0, cc = 0, xs = 0, xc = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        const double ph = 2 * kPi * f0 * i / fs, s = std::sin(ph), c = std::cos(ph);
        ss += s * s; sc += s * c; cc += c * c; xs += x[i] * s; xc += x[i] * c;
    }
    const double det = ss * cc - sc * sc;
    const double a = (xs * cc - xc * sc) / det, b = (xc * ss - xs * sc) / det;
    return std::hypot(a, b);
}

// Baseline wander: the signal's 1-second moving average -- what is left when
// every heartbeat-scale feature is averaged out. Its RMS is the wander.
static double wanderRms(const std::vector<double>& x, double fs, size_t i0, size_t i1) {
    const size_t w = static_cast<size_t>(fs);
    double run = 0, sum2 = 0; size_t n = 0;
    for (size_t i = i0; i < i1; ++i) {
        run += x[i];
        if (i >= i0 + w) run -= x[i - w];
        if (i + 1 >= i0 + w) { const double m = run / w; sum2 += m * m; ++n; }
    }
    double mean = 0; { double r = 0; size_t k = 0; for (size_t i = i0; i < i1; ++i) { r += x[i]; ++k; } mean = r / k; }
    return n ? std::sqrt(std::max(0.0, sum2 / n - mean * mean)) : NAN;
}

// R-peak heights above their own preceding baseline, in [i0, i1): peaks are
// found on `find` and measured on `measure`, so both views measure the same beats.
static std::vector<double> rHeights(const std::vector<double>& find, const std::vector<double>& measure,
    double fs, size_t i0, size_t i1, double thr) {
    std::vector<size_t> pk;
    const size_t refractory = static_cast<size_t>(0.3 * fs);
    for (size_t i = i0 + 1; i + 1 < i1; ++i) {
        if (find[i] <= thr || find[i] < find[i - 1] || find[i] <= find[i + 1]) continue;
        if (!pk.empty() && i - pk.back() < refractory) { if (find[i] > find[pk.back()]) pk.back() = i; continue; }
        pk.push_back(i);
    }
    std::vector<double> h;
    const size_t b0 = static_cast<size_t>(0.25 * fs), b1 = static_cast<size_t>(0.10 * fs);
    for (size_t i : pk) {
        if (i < i0 + b0) continue;
        std::vector<double> base(measure.begin() + (i - b0), measure.begin() + (i - b1));
        std::nth_element(base.begin(), base.begin() + base.size() / 2, base.end());
        h.push_back(measure[i] - base[base.size() / 2]);
    }
    return h;
}
static double median(std::vector<double> v) {
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
}

TEST(test1b_noise_filters, FixturesLoad) {
    for (const fs::path& p : { kCleanBin, kNoisyBin }) {
        SCOPED_TRACE(p.filename().string());
        const Ecg e = readEcg1(p);
        ASSERT_FALSE(e.v.empty()) << "cannot read ECG1 from " << p.string();
        EXPECT_EQ(e.fs, 1000.0);
        EXPECT_EQ(e.v.size(), 600000u);   // 10 min
    }
}

TEST(test1b_noise_filters, HighPassOffLeavesSignalUntouched) {
    const Ecg e = readEcg1(kCleanBin);
    ASSERT_FALSE(e.v.empty());
    EXPECT_TRUE(viewAsGui(e, 0.0, 0.0) == e.v);
}

TEST(test1b_noise_filters, HighPassRemovesBaselineWander) {
    const Ecg e = readEcg1(kCleanBin);
    ASSERT_FALSE(e.v.empty());
    const std::vector<double> off = viewAsGui(e, 0.0, 0.0), on = viewAsGui(e, kHighPassHz, 0.0);
    // the movement stretch, 190 s to the end
    const size_t i0 = static_cast<size_t>(190 * e.fs), i1 = e.v.size();
    const double wOff = wanderRms(off, e.fs, i0, i1), wOn = wanderRms(on, e.fs, i0, i1);
    EXPECT_GT(wOff, 0.5) << "the fixture should have real wander here; got " << wOff << " mV";
    EXPECT_LT(wOn, 0.1 * wOff) << "wander " << wOff << " -> " << wOn << " mV: less than 20 dB removed";
}

TEST(test1b_noise_filters, HighPassKeepsTheBeats) {
    const Ecg e = readEcg1(kCleanBin);
    ASSERT_FALSE(e.v.empty());
    const std::vector<double> off = viewAsGui(e, 0.0, 0.0), on = viewAsGui(e, kHighPassHz, 0.0);
    // the clean stretch, first 3 min: R heights above their own baseline
    const size_t i1 = static_cast<size_t>(180 * e.fs);
    const std::vector<double> hOff = rHeights(on, off, e.fs, 0, i1, 0.5), hOn = rHeights(on, on, e.fs, 0, i1, 0.5);
    ASSERT_GT(hOn.size(), 150u) << "too few beats found to judge";
    const double mOff = median(hOff), mOn = median(hOn);
    EXPECT_NEAR(mOn / mOff, 1.0, 0.02) << "R height " << mOff << " -> " << mOn << " mV";
}

TEST(test1b_noise_filters, NotchOffShowsInjected60Hz) {
    const Ecg e = readEcg1(kNoisyBin);
    ASSERT_FALSE(e.v.empty());
    const double a = toneAmplitude(viewAsGui(e, kHighPassHz, 0.0), e.fs, kInjectedHz);
    EXPECT_NEAR(a, kInjectedMv, 0.005) << "60 Hz with notch off: " << a << " mV (0.10 mV was injected)";
}

TEST(test1b_noise_filters, NotchOnRemoves60Hz) {
    const Ecg e = readEcg1(kNoisyBin);
    ASSERT_FALSE(e.v.empty());
    const double aOff = toneAmplitude(viewAsGui(e, kHighPassHz, 0.0), e.fs, kInjectedHz);
    const double aOn = toneAmplitude(viewAsGui(e, kHighPassHz, kNotchHz), e.fs, kInjectedHz);
    EXPECT_LT(aOn, 0.01 * aOff) << "60 Hz " << aOff << " -> " << aOn << " mV: less than 40 dB removed";
}

TEST(test1b_noise_filters, NotchOnLeavesEcgIntact) {
    const Ecg noisy = readEcg1(kNoisyBin), clean = readEcg1(kCleanBin);
    ASSERT_FALSE(noisy.v.empty());
    ASSERT_EQ(noisy.v.size(), clean.v.size());
    const std::vector<double> a = viewAsGui(noisy, kHighPassHz, kNotchHz), b = viewAsGui(clean, kHighPassHz, kNotchHz);
    double s2 = 0;
    for (size_t i = 0; i < a.size(); ++i) s2 += (a[i] - b[i]) * (a[i] - b[i]);
    const double rms = std::sqrt(s2 / a.size());
    EXPECT_LT(rms, 0.01) << "noisy and clean ECG differ by " << rms << " mV rms with the notch on";
}

TEST(test1b_noise_filters, NotchToggledOffRestoresIt) {
    const Ecg e = readEcg1(kNoisyBin);
    ASSERT_FALSE(e.v.empty());
    const std::vector<double> before = viewAsGui(e, kHighPassHz, 0.0);
    (void)viewAsGui(e, kHighPassHz, kNotchHz);                     // on ...
    const std::vector<double> after = viewAsGui(e, kHighPassHz, 0.0);   // ... and off again
    EXPECT_TRUE(after == before) << "turning the notch off did not give back what was there";
    EXPECT_NEAR(toneAmplitude(after, e.fs, kInjectedHz), kInjectedMv, 0.005);
}

// ============================================================================
// PART 3 -- THE FIGURES (Qt Gui builds only)
// ============================================================================
#ifdef QT_GUI_LIB
namespace {
    // QPainter needs a running QGuiApplication for text. The Google Test main
    // does not make one, so the first figure does; with no window shown.
    void ensureGuiApp() {
        if (QCoreApplication::instance()) return;
        static int argc = 1;
        static char name[] = "qtvi_tests";
        static char* argv[] = { name, nullptr };
        static QGuiApplication app(argc, argv);
    }

    struct Series { std::vector<double> x, y; QColor colour; };

    // One panel: box, title, min/max axis labels, the series, an optional
    // dotted vertical marker, optional log y.
    void panel(QPainter& p, const QRectF& r, const QString& title, const Series& s,
        double yLo, double yHi, bool logY, const QString& xLabel, const QString& yLabel, double markerX = NAN) {
        const QRectF plot = r.adjusted(60, 28, -10, -34);
        p.setPen(QPen(Qt::black, 1));
        p.drawRect(plot);
        QFont f = p.font(); f.setPointSize(10); f.setBold(true); p.setFont(f);
        p.drawText(QRectF(r.left(), r.top() + 2, r.width(), 22), Qt::AlignHCenter, title);
        f.setBold(false); f.setPointSize(8); p.setFont(f);
        const double x0 = s.x.front(), x1 = s.x.back();
        const double ly0 = logY ? std::log10(yLo) : yLo, ly1 = logY ? std::log10(yHi) : yHi;
        auto px = [&](double x) { return plot.left() + (x - x0) / (x1 - x0) * plot.width(); };
        auto py = [&](double y) { const double v = logY ? std::log10(std::max(y, yLo)) : y;
        return plot.bottom() - (std::clamp(v, ly0, ly1) - ly0) / (ly1 - ly0) * plot.height(); };
        p.drawText(QPointF(plot.left(), plot.bottom() + 14), QString::number(x0));
        p.drawText(QPointF(plot.right() - 40, plot.bottom() + 14), QString::number(x1));
        p.drawText(QRectF(plot.left(), plot.bottom() + 16, plot.width(), 16), Qt::AlignHCenter, xLabel);
        p.drawText(QPointF(r.left() + 2, plot.top() + 10), logY ? QString("1e%1").arg(ly1) : QString::number(yHi));
        p.drawText(QPointF(r.left() + 2, plot.bottom()), logY ? QString("1e%1").arg(ly0) : QString::number(yLo));
        p.save(); p.translate(r.left() + 12, plot.center().y()); p.rotate(-90);
        p.drawText(QRectF(-60, -10, 120, 16), Qt::AlignHCenter, yLabel); p.restore();
        if (std::isfinite(markerX)) {
            p.setPen(QPen(Qt::gray, 1, Qt::DotLine));
            p.drawLine(QPointF(px(markerX), plot.top()), QPointF(px(markerX), plot.bottom()));
        }
        QPolygonF poly;
        poly.reserve(static_cast<int>(s.x.size()));
        for (size_t i = 0; i < s.x.size(); ++i) poly << QPointF(px(s.x[i]), py(s.y[i]));
        p.setPen(QPen(s.colour, 1));
        p.setClipRect(plot);
        p.drawPolyline(poly);
        p.setClipping(false);
    }

    // Amplitude spectrum 0..fMax by direct DFT on a Hann-windowed stretch.
    Series spectrum(const std::vector<double>& x, double fs, size_t i0, size_t n, double fMax, double df, QColor c) {
        Series s; s.colour = c;
        double wsum = 0; std::vector<double> w(n);
        for (size_t k = 0; k < n; ++k) { w[k] = 0.5 - 0.5 * std::cos(2 * kPi * k / (n - 1)); wsum += w[k]; }
        double mean = 0; for (size_t k = 0; k < n; ++k) mean += x[i0 + k]; mean /= n;
        for (double f = 0; f <= fMax; f += df) {
            double re = 0, im = 0;
            for (size_t k = 0; k < n; ++k) { const double v = (x[i0 + k] - mean) * w[k], ph = 2 * kPi * f * k / fs; re += v * std::cos(ph); im -= v * std::sin(ph); }
            s.x.push_back(f); s.y.push_back(2 * std::hypot(re, im) / wsum);
        }
        return s;
    }

    fs::path outDir() {
        const fs::path d = fs::path(TESTS_DATA_DIR).parent_path() / "test_output";
        fs::create_directories(d);
        return d;
    }
}

TEST(test1b_noise_filters, WritesFigures) {
    ensureGuiApp();
    const Ecg clean = readEcg1(kCleanBin), noisy = readEcg1(kNoisyBin);
    ASSERT_FALSE(clean.v.empty());
    ASSERT_FALSE(noisy.v.empty());
    const double fs = clean.fs;

    // One stretch of signal as a series (every 2nd sample is plenty for a panel).
    auto signal = [&](const std::vector<double>& v, double t0, double t1, QColor c) {
        Series sr; sr.colour = c;
        for (size_t i = static_cast<size_t>(t0 * fs); i < static_cast<size_t>(t1 * fs) && i < v.size(); i += 2) {
            sr.x.push_back(i / fs); sr.y.push_back(v[i]);
        }
        return sr;
        };
    auto range = [](const Series& a, const Series& b, double& lo, double& hi) {
        lo = 1e9; hi = -1e9;
        for (const Series* sr : { &a, &b }) for (double y : sr->y) { lo = std::min(lo, y); hi = std::max(hi, y); }
        const double pad = 0.05 * (hi - lo); lo -= pad; hi += pad;
        };
    // Two rows, OFF left and ON right: the signal on top, its spectrum below.
    auto figure = [&](const char* file, const char* what, const char* offName, const char* onName,
        const Series& sigOff, const Series& sigOn, const Series& fftOff, const Series& fftOn, double markerHz) {
            QImage img(1600, 840, QImage::Format_RGB32); img.fill(Qt::white);
            QPainter p(&img); p.setRenderHint(QPainter::Antialiasing);
            double lo, hi; range(sigOff, sigOn, lo, hi);   // same y scale left and right
            panel(p, QRectF(0, 0, 800, 420), QString("%1: %2").arg(what, offName), sigOff, lo, hi, false, "s into the segment", "mV");
            panel(p, QRectF(800, 0, 800, 420), QString("%1: %2").arg(what, onName), sigOn, lo, hi, false, "s into the segment", "mV");
            panel(p, QRectF(0, 420, 800, 420), QString("spectrum: %1").arg(offName), fftOff, 1e-7, 1.0, true, "Hz", "amplitude (mV)", markerHz);
            panel(p, QRectF(800, 420, 800, 420), QString("spectrum: %1").arg(onName), fftOn, 1e-7, 1.0, true, "Hz", "amplitude (mV)", markerHz);
            p.end();
            const fs::path out = outDir() / file;
            EXPECT_TRUE(img.save(QString::fromStdWString(out.wstring()))) << "could not write " << out.string();
        };

    {   // HIGH PASS on the original file: the movement stretch, and 0-5 Hz where the wander lives
        const std::vector<double> off = viewAsGui(clean, 0.0, 0.0), on = viewAsGui(clean, kHighPassHz, 0.0);
        const size_t i0 = static_cast<size_t>(240 * fs), n = static_cast<size_t>(60 * fs);
        figure("test1b_highpass.png", "BASELINE WANDER, 240-270 s", "high pass OFF", "high pass ON (0.5 Hz)",
            signal(off, 240, 270, QColor(31, 119, 180)), signal(on, 240, 270, QColor(255, 127, 14)),
            spectrum(off, fs, i0, n, 5.0, 0.02, QColor(31, 119, 180)),
            spectrum(on, fs, i0, n, 5.0, 0.02, QColor(255, 127, 14)), kHighPassHz);
    }
    {   // NOTCH on the 60 Hz file, high pass on in both: a close-up, and 0-130 Hz
        const std::vector<double> off = viewAsGui(noisy, kHighPassHz, 0.0), on = viewAsGui(noisy, kHighPassHz, kNotchHz);
        const size_t i0 = static_cast<size_t>(60 * fs), n = static_cast<size_t>(20 * fs);
        figure("test1b_notch.png", "60 Hz INJECTED, 61-62.6 s, high pass ON", "notch OFF", "notch ON",
            signal(off, 61.0, 62.6, QColor(214, 39, 40)), signal(on, 61.0, 62.6, QColor(44, 160, 44)),
            spectrum(off, fs, i0, n, 130.0, 0.05, QColor(214, 39, 40)),
            spectrum(on, fs, i0, n, 130.0, 0.05, QColor(44, 160, 44)), kInjectedHz);
    }
}
#endif   // QT_GUI_LIB