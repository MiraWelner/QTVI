#pragma once
/**
 * @file   record_summary_plot.hpp
 * @brief  Task H's single-page QA output: draws a DiagnosticSummary as one
 *         SVG page, <dir>/<stem>_diagnostic_summary.svg.
 *
 *         SVG, AND NO QT. The summary is written from analysis_job, which is
 *         on the Qt-free side, and a vector page opens in any browser and
 *         prints at any size. The CSV beside it holds every number; this page
 *         is for looking.
 *
 *         LAYOUT, landscape, 3 x 2 panels plus a header and a Task J strip:
 *           chi-sq_0 vs beat     chi-sq_abs vs beat     normal/abnormal tracings
 *           RR vs PP             RR vs QT               Poincare RR(n) vs RR(n+1)
 *
 *         ROBUST AXES. Scatter ranges are the 0.5th..99.5th percentiles,
 *         padded; points outside are not drawn and are COUNTED in the panel's
 *         subtitle, so an outlier is reported rather than allowed to squash
 *         the panel. The chi-sq panels are log10, since the statistic spans
 *         orders of magnitude; non-positive values are counted the same way.
 *         The tracing panel uses its FULL range, since its extremes are the
 *         abnormal beats it exists to show.
 *         Panels with more than kMaxDrawn points are thinned evenly for
 *         drawing only -- the subtitle gives the true count.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "logging/record_summary.hpp"

namespace record_summary_plot {

    inline constexpr std::size_t kMaxDrawn = 15000;
    inline constexpr double kPageW = 1650.0, kPageH = 960.0;

    namespace detail {
        using Pts = std::vector<std::pair<double, double>>;

        inline std::string num(double v, int prec = 4) {
            char b[48];
            std::snprintf(b, sizeof b, "%.*g", prec, v);
            return b;
        }
        inline std::string esc(const std::string& s) {
            std::string o;
            for (char c : s) {
                if (c == '&') o += "&amp;"; else if (c == '<') o += "&lt;";
                else if (c == '>') o += "&gt;"; else o += c;
            }
            return o;
        }
        inline double pct(std::vector<double> v, double p) {
            if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
            const double h = (v.size() - 1) * p;
            const auto i = static_cast<std::size_t>(std::floor(h));
            std::nth_element(v.begin(), v.begin() + i, v.end());
            const double a = v[i];
            if (i + 1 >= v.size()) return a;
            const double b = *std::min_element(v.begin() + i + 1, v.end());
            return a + (h - i) * (b - a);
        }
        struct Range { double lo, hi; bool ok() const { return std::isfinite(lo) && std::isfinite(hi) && hi > lo; } };
        inline Range robust(const std::vector<double>& v, bool full = false) {
            std::vector<double> f;
            for (double x : v) if (std::isfinite(x)) f.push_back(x);
            if (f.empty()) return { NAN, NAN };
            double lo = full ? *std::min_element(f.begin(), f.end()) : pct(f, 0.005);
            double hi = full ? *std::max_element(f.begin(), f.end()) : pct(f, 0.995);
            if (!(hi > lo)) { const double d = std::max(1.0, std::abs(lo) * 0.1); lo -= d; hi += d; }
            const double pad = 0.05 * (hi - lo);
            return { lo - pad, hi + pad };
        }
        // 4-7 "nice" ticks over [lo, hi].
        inline std::vector<double> ticks(double lo, double hi) {
            std::vector<double> t;
            if (!(hi > lo)) return t;
            const double raw = (hi - lo) / 5.0;
            const double mag = std::pow(10.0, std::floor(std::log10(raw)));
            const double r = raw / mag;
            const double step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
            for (double v = std::ceil(lo / step) * step; v <= hi + 1e-9 * step; v += step) t.push_back(v);
            return t;
        }

        // One panel: frame, axes, ticks, labels; maps data to page coords.
        struct Panel {
            std::ostringstream& o;
            double x, y, w, h;   // plot area on the page
            Range xr, yr;
            bool logY = false;
            double px(double v) const { return x + (v - xr.lo) / (xr.hi - xr.lo) * w; }
            double py(double v) const {
                const double u = logY ? std::log10(v) : v;
                return y + h - (u - yr.lo) / (yr.hi - yr.lo) * h;
            }
            bool inside(double vx, double vy) const {
                if (!std::isfinite(vx) || !std::isfinite(vy)) return false;
                if (logY && !(vy > 0.0)) return false;
                const double u = logY ? std::log10(vy) : vy;
                return vx >= xr.lo && vx <= xr.hi && u >= yr.lo && u <= yr.hi;
            }
        };

        inline void frame(std::ostringstream& o, double x, double y, double w, double h,
            const std::string& title, const std::string& sub) {
            o << "<text x='" << x << "' y='" << y - 26 << "' class='t'>" << esc(title) << "</text>\n";
            o << "<text x='" << x << "' y='" << y - 10 << "' class='s'>" << esc(sub) << "</text>\n";
            o << "<rect x='" << x << "' y='" << y << "' width='" << w << "' height='" << h
                << "' class='f'/>\n";
        }
        inline void empty(std::ostringstream& o, double x, double y, double w, double h, const std::string& why) {
            o << "<text x='" << x + w / 2 << "' y='" << y + h / 2 << "' class='e'>" << esc(why) << "</text>\n";
        }
        inline void axes(const Panel& p, const std::string& xl, const std::string& yl) {
            std::ostringstream& o = p.o;
            for (double t : ticks(p.xr.lo, p.xr.hi)) {
                const double X = p.px(t);
                o << "<line x1='" << X << "' y1='" << p.y << "' x2='" << X << "' y2='" << p.y + p.h << "' class='g'/>"
                    << "<text x='" << X << "' y='" << p.y + p.h + 14 << "' class='tx'>" << num(t) << "</text>\n";
            }
            for (double t : ticks(p.yr.lo, p.yr.hi)) {
                const double Y = p.y + p.h - (t - p.yr.lo) / (p.yr.hi - p.yr.lo) * p.h;
                const std::string lbl = p.logY ? num(std::pow(10.0, t), 3) : num(t);
                o << "<line x1='" << p.x << "' y1='" << Y << "' x2='" << p.x + p.w << "' y2='" << Y << "' class='g'/>"
                    << "<text x='" << p.x - 5 << "' y='" << Y + 4 << "' class='ty'>" << lbl << "</text>\n";
            }
            o << "<text x='" << p.x + p.w / 2 << "' y='" << p.y + p.h + 32 << "' class='al'>" << esc(xl) << "</text>\n";
            o << "<text transform='translate(" << p.x - 48 << "," << p.y + p.h / 2 << ") rotate(-90)' class='al'>"
                << esc(yl) << "</text>\n";
        }
        inline std::size_t dots(const Panel& p, const Pts& pts, const char* cls, double r = 1.6) {
            std::size_t outside = 0;
            const std::size_t stride = std::max<std::size_t>(1, (pts.size() + kMaxDrawn - 1) / kMaxDrawn);
            for (std::size_t i = 0; i < pts.size(); ++i) {
                if (!p.inside(pts[i].first, pts[i].second)) { ++outside; continue; }
                if (i % stride) continue;
                p.o << "<circle cx='" << num(p.px(pts[i].first), 6) << "' cy='" << num(p.py(pts[i].second), 6)
                    << "' r='" << r << "' class='" << cls << "'/>";
            }
            p.o << "\n";
            return outside;
        }
        inline void identity(const Panel& p) {
            const double lo = std::max(p.xr.lo, p.yr.lo), hi = std::min(p.xr.hi, p.yr.hi);
            if (hi > lo)
                p.o << "<line x1='" << p.px(lo) << "' y1='" << p.py(lo) << "' x2='" << p.px(hi) << "' y2='" << p.py(hi)
                << "' class='id'/>\n";
        }

        // A scatter panel with robust ranges.
        inline void scatter(std::ostringstream& o, double x, double y, double w, double h,
            const std::string& title, const Pts& pts, const std::string& xl, const std::string& yl,
            bool idLine, const char* cls = "p", const Pts* second = nullptr, const char* cls2 = "pa") {
            std::vector<double> xs, ys;
            for (const auto& q : pts) { xs.push_back(q.first); ys.push_back(q.second); }
            if (second) for (const auto& q : *second) { xs.push_back(q.first); ys.push_back(q.second); }
            const std::size_t n = xs.size();
            const bool small = n < 50;   // few points: show all of them
            Panel p{ o, x, y, w, h, robust(xs, small), robust(ys, small) };
            if (n == 0 || !p.xr.ok() || !p.yr.ok()) {
                frame(o, x, y, w, h, title, "n = 0");
                empty(o, x, y, w, h, "no data");
                return;
            }
            std::ostringstream body;
            Panel pb{ body, x, y, w, h, p.xr, p.yr };
            axes(pb, xl, yl);
            if (idLine) identity(pb);
            std::size_t outside = dots(pb, pts, cls, second ? 4.0 : 1.6);
            if (second) outside += dots(pb, *second, cls2, 4.0);
            frame(o, x, y, w, h, title,
                "n = " + std::to_string(n) + (outside ? ", " + std::to_string(outside) + " outside range" : ""));
            o << body.str();
        }

        // chi-sq vs beat, log10 y.
        inline void series(std::ostringstream& o, double x, double y, double w, double h,
            const std::string& title, const std::vector<double>& v, const std::string& yl) {
            Pts pts;
            std::vector<double> logs;
            std::size_t nonpos = 0;
            for (std::size_t i = 0; i < v.size(); ++i) {
                if (!std::isfinite(v[i])) continue;
                if (v[i] <= 0.0) { ++nonpos; continue; }
                pts.push_back({ static_cast<double>(i), v[i] });
                logs.push_back(std::log10(v[i]));
            }
            if (pts.empty()) {
                frame(o, x, y, w, h, title, "n = 0");
                empty(o, x, y, w, h, "no scored beats");
                return;
            }
            Range xr{ 0.0, std::max(1.0, static_cast<double>(v.size())) };
            Panel p{ o, x, y, w, h, xr, robust(logs) };
            p.logY = true;
            std::ostringstream body;
            Panel pb{ body, x, y, w, h, p.xr, p.yr, true };
            axes(pb, "beat number", yl);
            const std::size_t outside = dots(pb, pts, "p", 1.2) + nonpos;
            frame(o, x, y, w, h, title, "n = " + std::to_string(pts.size() + nonpos)
                + (outside ? ", " + std::to_string(outside) + " outside range" : ""));
            o << body.str();
        }

        // Superimposed tracings: normal grey, abnormal red, drawn as polylines.
        inline void tracings(std::ostringstream& o, double x, double y, double w, double h,
            const record_summary::DiagnosticSummary& ds) {
            const std::string title = "Superimposed beats, CH" + std::to_string(ds.lead + 1);
            const std::string sub = "normal " + std::to_string(ds.normalBeats.size())
                + " (grey), abnormal " + std::to_string(ds.abnormalBeats.size()) + " (red)";
            std::vector<double> ys;
            for (const auto* set : { &ds.normalBeats, &ds.abnormalBeats })
                for (const auto& tr : *set) for (double v : tr) if (std::isfinite(v)) ys.push_back(v);
            if (ys.empty() || !std::isfinite(ds.trace_dt_ms)) {
                frame(o, x, y, w, h, title, sub);
                empty(o, x, y, w, h, "no tracings");
                return;
            }
            std::size_t len = 0;
            for (const auto* set : { &ds.normalBeats, &ds.abnormalBeats })
                for (const auto& tr : *set) len = std::max(len, tr.size());
            const Range xr{ ds.trace_t0_ms, ds.trace_t0_ms + ds.trace_dt_ms * (len ? len - 1 : 1) };
            frame(o, x, y, w, h, title, sub);
            // FULL range, not robust: the extremes here are the abnormal
            // morphologies this panel exists to show.
            Panel p{ o, x, y, w, h, xr, robust(ys, true) };
            axes(p, "ms from R", "amplitude");
            auto draw = [&](const std::vector<std::vector<double>>& set, const char* cls) {
                for (const auto& tr : set) {
                    std::string pts;
                    auto flush = [&] {
                        if (!pts.empty()) o << "<polyline points='" << pts << "' class='" << cls << "'/>";
                        pts.clear();
                        };
                    for (std::size_t i = 0; i < tr.size(); ++i) {
                        const double t = ds.trace_t0_ms + ds.trace_dt_ms * i;
                        if (!p.inside(t, tr[i])) { flush(); continue; }
                        pts += num(p.px(t), 6) + "," + num(p.py(tr[i]), 6) + " ";
                    }
                    flush();
                    o << "\n";
                }
                };
            draw(ds.normalBeats, "tn");
            draw(ds.abnormalBeats, "ta");
        }
    }

    inline std::string renderSvg(const record_summary::DiagnosticSummary& ds) {
        using namespace detail;
        std::ostringstream o;
        o << "<svg xmlns='http://www.w3.org/2000/svg' width='" << kPageW << "' height='" << kPageH
            << "' viewBox='0 0 " << kPageW << ' ' << kPageH << "' font-family='Helvetica,Arial,sans-serif'>\n"
            "<style>"
            ".t{font-size:14px;font-weight:bold}.s{font-size:11px;fill:#555}"
            ".f{fill:none;stroke:#333;stroke-width:1}.g{stroke:#e4e4e4;stroke-width:1}"
            ".tx{font-size:10px;text-anchor:middle;fill:#333}.ty{font-size:10px;text-anchor:end;fill:#333}"
            ".al{font-size:11px;text-anchor:middle}.e{font-size:13px;text-anchor:middle;fill:#999}"
            ".p{fill:#1f4e9a;fill-opacity:0.45}.pa{fill:#c0392b;fill-opacity:0.85}"
            ".id{stroke:#999;stroke-dasharray:4 3}"
            ".tn{fill:none;stroke:#777;stroke-opacity:0.18;stroke-width:1}"
            ".ta{fill:none;stroke:#c0392b;stroke-opacity:0.35;stroke-width:1}"
            "</style>\n<rect width='100%' height='100%' fill='white'/>\n";

        const std::size_t scored = std::count_if(ds.chiSq0_vs_beat.begin(), ds.chiSq0_vs_beat.end(),
            [](double v) { return std::isfinite(v); });
        o << "<text x='40' y='34' font-size='20' font-weight='bold'>Record diagnostic summary: "
            << esc(ds.recordID) << "</text>\n"
            << "<text x='40' y='54' font-size='12' fill='#555'>" << ds.n_beats << " beats; "
            << scored << " scored on CH" << ds.lead + 1 << "; " << ds.qt_templates.size() << " templates"
            << "; P waves found on " << ds.n_p_waves << " beats" << "</text>\n";

        // Grid: 3 columns x 2 rows of plot areas.
        const double left = 95, top = 110, colW = 520, rowH = 420, pw = 430, ph = 330;
        auto X = [&](int c) { return left + c * colW; };
        auto Y = [&](int r) { return top + r * rowH; };

        series(o, X(0), Y(0), pw, ph, "chi-sq_0 vs beat number", ds.chiSq0_vs_beat, "chi-sq_0 (log)");
        series(o, X(1), Y(0), pw, ph, "chi-sq_abs vs beat number", ds.chiSqAbs_vs_beat, "chi-sq_abs (log)");
        tracings(o, X(2), Y(0), pw, ph, ds);

        scatter(o, X(0), Y(1), pw, ph, "RR vs PP (P wave to P wave)", ds.rr_vs_pp, "RR (ms)", "PP (ms)", true);
        {
            Pts normal, abnormal;
            for (const auto& t : ds.qt_templates)
                if (std::isfinite(t.rr_ms) && std::isfinite(t.qt_ms))
                    (t.abnormal ? abnormal : normal).push_back({ t.rr_ms, t.qt_ms });
            scatter(o, X(1), Y(1), pw, ph, "RR vs QT, one point per template (red = ectopic)",
                normal, "median RR (ms)", "QT (ms)", false, "p", &abnormal, "pa");
        }
        scatter(o, X(2), Y(1), pw, ph, "Poincare: RR(n) vs RR(n+1)", ds.poincare, "RR(n) (ms)", "RR(n+1) (ms)", true);

        // Task J strip.
        const double sy = kPageH - 40;
        o << "<rect x='40' y='" << sy - 18 << "' width='" << kPageW - 80 << "' height='26' fill='#f4f4f4' stroke='#ccc'/>"
            << "<text x='50' y='" << sy << "' font-size='12' fill='#555'>State-stratified features: "
            << esc(ds.stratified.available ? std::string("see CSV") : ds.stratified.note) << "</text>\n";
        o << "</svg>\n";
        return o.str();
    }

    inline bool writeSummarySvg(const std::string& dir, const std::string& stem,
        const record_summary::DiagnosticSummary& ds) {
        const std::string path = dir + "/" + stem + "_diagnostic_summary.svg";
        std::ofstream f(path);
        if (!f.is_open()) {
            std::cerr << "  WARNING: could not open " << path << " for the summary page\n";
            return false;
        }
        f << renderSvg(ds);
        return static_cast<bool>(f);
    }

}  // namespace record_summary_plot