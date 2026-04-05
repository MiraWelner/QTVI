#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <functional>
#include <map>
#include <optional>
#include <stdexcept>
#include <cstdint>
#include <thread>
#include <mutex>
#include <future>
#include <queue>
#include <condition_variable>

/**
 * @file common.hpp
 * @brief Shared types, math utilities, and thread pool for the PPG feature pipeline.
 */

namespace ppg {

    // ??? Constants ???????????????????????????????????????????????????????????????

    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    constexpr double kInf = std::numeric_limits<double>::infinity();

    // ??? Scalar helpers ??????????????????????????????????????????????????????????

    /** @brief Check if value is NaN. */
    inline bool isnan(double v) { return std::isnan(v); }

    /**
     * @brief Safe array access by integer index.
     * @param arr  Source vector.
     * @param idx  Zero-based index.
     * @return Value at index, or NaN if out of bounds.
     */
    inline double safe_get(const std::vector<double>& arr, int idx) {
        if (idx < 0 || idx >= static_cast<int>(arr.size())) return kNaN;
        return arr[idx];
    }

    /**
     * @brief Safe array access by floating-point index (rounded).
     * @param arr  Source vector.
     * @param idx  Index (will be rounded). Returns NaN if NaN or out of bounds.
     */
    inline double safe_get(const std::vector<double>& arr, double idx) {
        if (std::isnan(idx)) return kNaN;
        return safe_get(arr, static_cast<int>(std::round(idx)));
    }

    // ??? Matrix (row-major) ?????????????????????????????????????????????????????

    /**
     * @brief Simple 2D row-major matrix.
     */
    struct Matrix {
        std::vector<double> data;
        int rows = 0;
        int cols = 0;

        Matrix() = default;

        /**
         * @param r    Row count.
         * @param c    Column count.
         * @param val  Fill value.
         */
        Matrix(int r, int c, double val = 0.0)
            : data(static_cast<size_t>(r)* c, val), rows(r), cols(c) {
        }

        double& operator()(int r, int c) { return data[r * cols + c]; }
        double  operator()(int r, int c) const { return data[r * cols + c]; }
    };

    // ??? Vector math ?????????????????????????????????????????????????????????????

    /** @brief Element-wise difference: out[i] = v[i+1] - v[i]. */
    inline std::vector<double> diff(const std::vector<double>& v) {
        if (v.size() < 2) return {};
        std::vector<double> d(v.size() - 1);
        for (size_t i = 0; i < d.size(); ++i)
            d[i] = v[i + 1] - v[i];
        return d;
    }

    /** @brief Euclidean (L2) distance between two equal-length vectors. */
    inline double l2_dist(const std::vector<double>& u, const std::vector<double>& v) {
        double s = 0.0;
        for (size_t i = 0; i < u.size(); ++i) {
            double d = u[i] - v[i];
            s += d * d;
        }
        return std::sqrt(s);
    }

    /**
     * @brief 1-D linear interpolation (like MATLAB interp1 'linear').
     * @param x_in   Sample x coordinates (strictly increasing).
     * @param y_in   Sample y values.
     * @param x_out  Query x coordinates.
     * @return Interpolated y values. Extrapolation is clamped to end values.
     */
    inline std::vector<double> interp1_linear(
        const std::vector<double>& x_in,
        const std::vector<double>& y_in,
        const std::vector<double>& x_out)
    {
        std::vector<double> y_out(x_out.size());
        const int n = static_cast<int>(x_in.size());
        for (size_t k = 0; k < x_out.size(); ++k) {
            double xq = x_out[k];
            if (xq <= x_in.front()) { y_out[k] = y_in.front(); continue; }
            if (xq >= x_in.back()) { y_out[k] = y_in.back();  continue; }
            // binary search for interval
            auto it = std::lower_bound(x_in.begin(), x_in.end(), xq);
            int i = static_cast<int>(it - x_in.begin());
            if (i == 0) i = 1;
            double t = (xq - x_in[i - 1]) / (x_in[i] - x_in[i - 1]);
            y_out[k] = y_in[i - 1] + t * (y_in[i] - y_in[i - 1]);
        }
        return y_out;
    }

    /**
     * @brief Cubic spline interpolation (natural spline).
     * @param x_in   Sample x coordinates (strictly increasing, ? 2 points).
     * @param y_in   Sample y values.
     * @param x_out  Query x coordinates.
     * @return Interpolated y values.
     */
    inline std::vector<double> interp1_spline(
        const std::vector<double>& x_in,
        const std::vector<double>& y_in,
        const std::vector<double>& x_out)
    {
        const int n = static_cast<int>(x_in.size());
        if (n < 2) return std::vector<double>(x_out.size(), y_in.empty() ? 0.0 : y_in[0]);

        // Compute natural cubic spline coefficients
        std::vector<double> h(n - 1), alpha(n - 1);
        for (int i = 0; i < n - 1; ++i) h[i] = x_in[i + 1] - x_in[i];
        for (int i = 1; i < n - 1; ++i)
            alpha[i] = 3.0 / h[i] * (y_in[i + 1] - y_in[i])
            - 3.0 / h[i - 1] * (y_in[i] - y_in[i - 1]);

        std::vector<double> c(n, 0.0), l(n, 1.0), mu(n, 0.0), z(n, 0.0);
        for (int i = 1; i < n - 1; ++i) {
            l[i] = 2.0 * (x_in[i + 1] - x_in[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }
        std::vector<double> b(n - 1), d(n - 1);
        for (int j = n - 2; j >= 0; --j) {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (y_in[j + 1] - y_in[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }

        // Evaluate
        std::vector<double> y_out(x_out.size());
        for (size_t k = 0; k < x_out.size(); ++k) {
            double xq = x_out[k];
            int seg = n - 2;
            if (xq <= x_in.front()) seg = 0;
            else if (xq < x_in.back()) {
                auto it = std::upper_bound(x_in.begin(), x_in.end(), xq);
                seg = static_cast<int>(it - x_in.begin()) - 1;
                if (seg < 0) seg = 0;
            }
            double dx = xq - x_in[seg];
            y_out[k] = y_in[seg] + b[seg] * dx + c[seg] * dx * dx + d[seg] * dx * dx * dx;
        }
        return y_out;
    }

    /**
     * @brief Generate linearly spaced vector.
     * @param start  First value.
     * @param stop   Last value.
     * @param n      Number of points (?1).
     */
    inline std::vector<double> linspace(double start, double stop, int n) {
        std::vector<double> v(n);
        if (n == 1) { v[0] = start; return v; }
        double step = (stop - start) / (n - 1);
        for (int i = 0; i < n; ++i) v[i] = start + i * step;
        return v;
    }

    /**
     * @brief Trapezoidal numerical integration.
     * @param x  x coordinates.
     * @param y  y values (same length as x).
     */
    inline double trapz(const std::vector<double>& x, const std::vector<double>& y) {
        double s = 0.0;
        for (size_t i = 1; i < x.size(); ++i)
            s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
        return s;
    }

    // ??? Simple fast smooth (moving-average via cumsum) ?????????????????????????

    /**
     * @brief Centered moving-average smooth (NaN-aware).
     * @param data  Input signal.
     * @param win   Window width (should be odd).
     * @return Smoothed signal, same length as input.
     */
    inline std::vector<double> fast_smooth(const std::vector<double>& data, int win) {
        const int n = static_cast<int>(data.size());
        if (n == 0 || win <= 1) return data;
        if (win % 2 == 0) ++win;
        int half = win / 2;
        std::vector<double> out(n, 0.0);
        // running sum
        double sum = 0.0;
        int cnt = 0;
        for (int i = 0; i < std::min(win, n); ++i) {
            if (!std::isnan(data[i])) { sum += data[i]; ++cnt; }
        }
        for (int i = 0; i < n; ++i) {
            int lo = i - half;
            int hi = i + half;
            // adjust window at edges
            if (i > half && hi < n) {
                if (!std::isnan(data[hi])) { sum += data[hi]; ++cnt; }
                if (lo - 1 >= 0 && !std::isnan(data[lo - 1])) { sum -= data[lo - 1]; --cnt; }
            }
            out[i] = (cnt > 0) ? sum / cnt : kNaN;
        }
        // simpler correct implementation: just brute force for clarity
        out.assign(n, 0.0);
        for (int i = 0; i < n; ++i) {
            int lo = std::max(0, i - half);
            int hi = std::min(n - 1, i + half);
            double s = 0.0; int c = 0;
            for (int j = lo; j <= hi; ++j) {
                if (!std::isnan(data[j])) { s += data[j]; ++c; }
            }
            out[i] = (c > 0) ? s / c : kNaN;
        }
        return out;
    }

    // ??? Run-length encoding ????????????????????????????????????????????????????

    struct RunLengthEntry {
        double value;
        int count;
    };

    /**
     * @brief Run-length encode a vector.
     * @param v  Input vector.
     * @return Vector of {value, count} entries.
     */
    inline std::vector<RunLengthEntry> run_length_encode(const std::vector<double>& v) {
        std::vector<RunLengthEntry> rle;
        if (v.empty()) return rle;
        double cur = v[0];
        int cnt = 1;
        auto same = [](double a, double b) {
            if (std::isnan(a) && std::isnan(b)) return true;
            return a == b;
            };
        for (size_t i = 1; i < v.size(); ++i) {
            if (same(v[i], cur)) { ++cnt; }
            else { rle.push_back({ cur, cnt }); cur = v[i]; cnt = 1; }
        }
        rle.push_back({ cur, cnt });
        return rle;
    }

    // ??? Thread pool ?????????????????????????????????????????????????????????????

    /**
     * @brief Minimal thread pool for per-bin parallel processing.
     */
    class ThreadPool {
    public:
        /**
         * @param num_threads  Number of worker threads (0 = hardware concurrency).
         */
        explicit ThreadPool(unsigned num_threads = 0) {
            if (num_threads == 0)
                num_threads = std::max(1u, std::thread::hardware_concurrency());
            for (unsigned i = 0; i < num_threads; ++i)
                workers_.emplace_back([this] { worker_loop(); });
        }

        ~ThreadPool() {
            { std::lock_guard<std::mutex> lk(mtx_); stop_ = true; }
            cv_.notify_all();
            for (auto& w : workers_) w.join();
        }

        /**
         * @brief Enqueue a callable and return a future to its result.
         * @param f  Callable (no arguments).
         * @return std::future holding the return value.
         */
        template <typename F>
        auto enqueue(F&& f) -> std::future<decltype(f())> {
            using R = decltype(f());
            auto task = std::make_shared<std::packaged_task<R()>>(std::forward<F>(f));
            std::future<R> fut = task->get_future();
            { std::lock_guard<std::mutex> lk(mtx_); tasks_.emplace([task] { (*task)(); }); }
            cv_.notify_one();
            return fut;
        }

    private:
        void worker_loop() {
            for (;;) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lk(mtx_);
                    cv_.wait(lk, [this] { return stop_ || !tasks_.empty(); });
                    if (stop_ && tasks_.empty()) return;
                    task = std::move(tasks_.front());
                    tasks_.pop();
                }
                task();
            }
        }

        std::vector<std::thread>          workers_;
        std::queue<std::function<void()>> tasks_;
        std::mutex                        mtx_;
        std::condition_variable           cv_;
        bool                              stop_ = false;
    };

} // namespace ppg