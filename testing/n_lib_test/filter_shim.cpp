// filter_shim.cpp
//
// Plain-C ABI shim around the FilterUtils header, so Python (ctypes) can call
// the same code that runs inside the C++ app. Each entry point takes raw
// double* + length, writes the result to an output buffer the caller owns.
//
// Build:
//   Linux:   g++ -std=c++17 -O2 -shared -fPIC -I<project> filter_shim.cpp -o libfilter_shim.so
//   Windows: cl /LD /std:c++17 /EHsc /O2 /I<project> filter_shim.cpp
//               /link /OUT:filter_shim.dll
//   macOS:   clang++ -std=c++17 -O2 -shared -fPIC -I<project> filter_shim.cpp
//               -o libfilter_shim.dylib
//
#include "FilterUtils.hpp"
#include <cstddef>
#include <cstring>
#include <vector>

// MSVC needs an explicit export decoration to make extern "C" symbols
// visible in the .dll import table; on GCC/Clang the empty macro is a no-op
// and default symbol visibility does the right thing.
#if defined(_WIN32) || defined(_WIN64)
#  define SHIM_API __declspec(dllexport)
#else
#  define SHIM_API
#endif

extern "C" {

// Design a Butterworth LP or HP filter, N-th order.
//   type_lp_hp: 0 = lowpass, 1 = highpass
//   b_out / a_out must have space for N+1 doubles each.
// Returns number of coefficients written (== N+1) or 0 on error.
SHIM_API int shim_butter_lp_hp(int N, double Wn, int type_lp_hp,
                      double* b_out, double* a_out) {
    if (!b_out || !a_out || N < 1) return 0;
    std::vector<double> b, a;
    const char* type = (type_lp_hp == 1) ? "high" : "low";
    butter(N, Wn, type, b, a);
    if (b.empty() || a.empty()) return 0;
    std::memcpy(b_out, b.data(), b.size() * sizeof(double));
    std::memcpy(a_out, a.data(), a.size() * sizeof(double));
    return static_cast<int>(b.size());
}

// Design a Butterworth bandpass filter, N-th order.
//   Wn_lo, Wn_hi in (0, 1) normalized to Nyquist.
//   b_out / a_out must have space for 2*N+1 doubles each.
// Returns number of coefficients written (== 2*N+1) or 0 on error.
SHIM_API int shim_butter_bp(int N, double Wn_lo, double Wn_hi,
                   double* b_out, double* a_out) {
    if (!b_out || !a_out || N < 1) return 0;
    std::vector<double> b, a;
    std::vector<double> Wn = { Wn_lo, Wn_hi };
    butter(N, Wn, b, a);
    if (b.empty() || a.empty()) return 0;
    std::memcpy(b_out, b.data(), b.size() * sizeof(double));
    std::memcpy(a_out, a.data(), a.size() * sizeof(double));
    return static_cast<int>(b.size());
}

// One-way IIR filter (Direct Form II Transposed, matches MATLAB filter()).
// y_out must have space for n doubles.
SHIM_API void shim_filter(const double* b, int nb, const double* a, int na,
                 const double* x, int n, double* y_out) {
    std::vector<double> bv(b, b + nb), av(a, a + na), xv(x, x + n);
    auto y = filter(bv, av, xv);
    std::memcpy(y_out, y.data(), y.size() * sizeof(double));
}

// Zero-phase filter (MATLAB filtfilt / SciPy filtfilt).
SHIM_API void shim_filtfilt(const double* b, int nb, const double* a, int na,
                   const double* x, int n, double* y_out) {
    std::vector<double> bv(b, b + nb), av(a, a + na), xv(x, x + n);
    auto y = filtfilt(bv, av, xv);
    std::memcpy(y_out, y.data(), y.size() * sizeof(double));
}

// Notch filter: (x - narrow_bandpass(x)) via filtfilt, NaN-aware.
SHIM_API void shim_notch_filter(const double* x, int n, double notch_hz, double fs,
                       double Q, int N, double* y_out) {
    std::vector<double> xv(x, x + n);
    auto y = notch_filter(xv, notch_hz, fs, Q, N);
    std::memcpy(y_out, y.data(), y.size() * sizeof(double));
}

// Waveform highpass: zero-phase HP via filtfilt, NaN-aware.
SHIM_API void shim_waveform_highpass(const double* x, int n, double cutoff_hz,
                            double fs, int N, double* y_out) {
    std::vector<double> xv(x, x + n);
    auto y = waveform_highpass(xv, cutoff_hz, fs, N);
    std::memcpy(y_out, y.data(), y.size() * sizeof(double));
}

}  // extern "C"
