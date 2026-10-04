#!/usr/bin/env python3
"""Put one signal through each C++ filter and its SciPy equivalent, plot both.

    python show_vs_scipy.py [folder]     # folder holds the .hpp files
"""
import ctypes, os, shutil, subprocess, sys, tempfile
import numpy as np
from scipy import signal
from scipy.ndimage import uniform_filter1d
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BIND = r'''
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include "filter_utils.hpp"
#include "nanfastsmooth.hpp"
#include "stats_utils.hpp"
#include "diff2.hpp"
#if defined(_WIN32)
#  define XP extern "C" __declspec(dllexport)
#else
#  define XP extern "C" __attribute__((visibility("default")))
#endif
static std::vector<double> V(const double* x, int n) { return std::vector<double>(x, x + n); }
static int out(const std::vector<double>& y, double* o) {
    if (!y.empty()) std::memcpy(o, y.data(), y.size() * sizeof(double));
    return (int)y.size();
}
XP int  c_upsample(const double* x, int n, double s, double t, double* o) { return out(filterutils::upsample(V(x,n), s, t), o); }
XP int  c_filtfilt(const double* x, int n, double fc, double fs, double* o) { return out(filterutils::filtfilt(filterutils::butterLP(fc, fs), V(x,n)), o); }
XP int  c_detrend (const double* x, int n, double* o) { return out(filterutils::detrend(V(x,n)), o); }
XP int  c_smooth  (const double* x, int n, double w, int ty, double* o) { return out(nanfastsmooth(V(x,n), w, ty, 0.5), o); }
XP int  c_movmean (const double* x, int n, int w, double* o) { return out(movmean(V(x,n), (size_t)w), o); }
XP int  c_diff2   (const double* x, int n, double* o) { return out(diff2(V(x,n), 1), o); }
'''


def _msvc_install():
    """Locate a Visual Studio install with the C++ toolset, via vswhere."""
    if not sys.platform.startswith("win"):
        return None
    for env in ("ProgramFiles(x86)", "ProgramFiles"):
        root = os.environ.get(env)
        if not root:
            continue
        vswhere = os.path.join(root, "Microsoft Visual Studio", "Installer",
                               "vswhere.exe")
        if not os.path.isfile(vswhere):
            continue
        p = subprocess.run([vswhere, "-latest", "-products", "*", "-requires",
                            "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
                            "-property", "installationPath"],
                           capture_output=True, text=True)
        for line in p.stdout.splitlines():
            c = line.strip()
            if c and os.path.isfile(os.path.join(c, "Common7", "Tools",
                                                 "VsDevCmd.bat")):
                return c
    return None


def build(folder):
    """Compile the bindings against the headers in `folder`. Tries MSVC first
    on Windows (via vswhere, so a plain PowerShell prompt works), then any
    g++/clang++ on PATH."""
    d = tempfile.mkdtemp(prefix="show_")
    inc = os.path.abspath(folder)
    src = os.path.join(d, "b.cpp")
    open(src, "w").write(BIND)
    lib = os.path.join(d, "b.dll" if sys.platform.startswith("win") else "b.so")

    attempts = []
    vs = _msvc_install()
    if vs:
        arch = "x64" if sys.maxsize > 2 ** 32 else "x86"
        bat = os.path.join(d, "build.bat")
        with open(bat, "w") as f:
            f.write("@echo off\r\n"
                    f'call "{os.path.join(vs, "Common7", "Tools", "VsDevCmd.bat")}"'
                    f" -arch={arch} -host_arch={arch} -no_logo || exit /b 1\r\n"
                    f'cl /nologo /LD /O2 /EHsc /std:c++17 /I"{inc}" "{src}" '
                    f'/Fe:"{lib}" /Fo:"{os.path.join(d, "b.obj")}" || exit /b 1\r\n')
        attempts.append((f"cl ({arch})", ["cmd", "/c", bat]))
    if shutil.which("cl"):
        attempts.append(("cl (on PATH)", [
            "cl", "/nologo", "/LD", "/O2", "/EHsc", "/std:c++17", f"/I{inc}",
            src, f"/Fe:{lib}", f"/Fo:{os.path.join(d, 'b.obj')}"]))
    for cxx in ("g++", "clang++"):
        if shutil.which(cxx):
            attempts.append((cxx, [cxx, "-O2", "-std=c++17", "-shared", "-fPIC",
                                   f"-I{inc}", src, "-o", lib, "-pthread"]))
    if not attempts:
        sys.exit("No C++ compiler found (looked for MSVC via vswhere, then "
                 "cl/g++/clang++ on PATH).")

    errors = []
    for name, cmd in attempts:
        p = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
        if p.returncode == 0 and os.path.isfile(lib):
            print(f"compiled via {name}")
            return ctypes.CDLL(lib), d
        errors.append(f"--- {name} ---\n{p.stdout.strip()}\n{p.stderr.strip()}")
    sys.exit("Could not build the bindings:\n\n" + "\n\n".join(errors))


def main(folder="."):
    L, tmp = build(folder)
    dp = ctypes.POINTER(ctypes.c_double)
    P = lambda a: a.ctypes.data_as(dp)

    def call(fn, x, *args, grow=1):
        x = np.ascontiguousarray(x, float)
        o = np.zeros(int(x.size * grow) + 8)
        fn.restype = ctypes.c_int
        n = fn(P(x), x.size, *args, P(o))
        return o[:n]

    # ---- one signal: trend + slow wave + fast wave + noise, 500 Hz, 10 s
    fs, n = 500.0, 5000
    t = np.arange(n) / fs
    rng = np.random.default_rng(0)
    x = (0.002 * np.arange(n) + np.sin(2 * np.pi * 1.5 * t)
         + 0.4 * np.sin(2 * np.pi * 45 * t) + 0.15 * rng.standard_normal(n))

    rows = []

    # 1. resampler 500 -> 1000 Hz
    c = call(L.c_upsample, x, ctypes.c_double(fs), ctypes.c_double(1000.0), grow=2)
    s = signal.resample_poly(x, 2, 1)
    rows.append(("upsample 500->1000 Hz", "scipy.signal.resample_poly",
                 np.arange(c.size) / 1000.0, c, s))

    # 2. butterLP + filtfilt, fc = 20 Hz
    c = call(L.c_filtfilt, x, ctypes.c_double(20.0), ctypes.c_double(fs))
    b, a = signal.butter(2, 20.0, "low", fs=fs)
    rows.append(("butterLP+filtfilt fc=20 Hz", "scipy.butter + scipy.filtfilt",
                 t, c, signal.filtfilt(b, a, x)))

    # 3. detrend
    rows.append(("detrend", "scipy.signal.detrend",
                 t, call(L.c_detrend, x), signal.detrend(x, type="linear")))

    # 4. nanfastsmooth w=31
    c = call(L.c_smooth, x, ctypes.c_double(31.0), ctypes.c_int(1))
    rows.append(("nanfastsmooth w=31", "scipy.ndimage.uniform_filter1d",
                 t, c, uniform_filter1d(x, 31, mode="nearest")))

    # 5. movmean w=31
    c = call(L.c_movmean, x, ctypes.c_int(31))
    rows.append(("movmean w=31", "scipy.ndimage.uniform_filter1d",
                 t, c, uniform_filter1d(x, 31, mode="nearest")))

    # 6. diff2  (closest SciPy analogue, NOT the same kernel)
    c = call(L.c_diff2, x)
    rows.append(("diff2", "scipy.signal.savgol_filter(5,2,deriv=1)",
                 t[:c.size], c, signal.savgol_filter(x, 5, 2, deriv=1)[:c.size]))

    fig, axes = plt.subplots(3, 2, figsize=(13.33, 7.5), dpi=150)
    fig.suptitle("same signal through each C++ filter and its SciPy equivalent",
                 fontsize=13, fontweight="bold")
    for ax, (name, ref, tt, yc, ys) in zip(axes.ravel(), rows):
        m = min(yc.size, ys.size)
        d_all = float(np.nanmax(np.abs(yc[:m] - ys[:m])))
        keep = (tt[:m] >= 2.0) & (tt[:m] <= 2.6)          # 0.6 s window
        d = float(np.nanmax(np.abs(yc[:m][keep] - ys[:m][keep])))
        ax.plot(tt[:m][keep], yc[:m][keep], color="tab:blue", lw=1.8, label="C++")
        ax.plot(tt[:m][keep], ys[:m][keep], color="tab:orange", lw=1.0,
                ls=(0, (2, 2)), label=ref)
        ax.set_title(f"{name}     max|C++ - SciPy| = {d:.2e} here, "
                     f"{d_all:.2e} incl. edges", fontsize=9)
        ax.legend(fontsize=7, loc="upper right", framealpha=.9)
        ax.grid(alpha=.3); ax.tick_params(labelsize=7)
        ax.set_xlabel("time (s)", fontsize=8)
        print(f"{name:28} vs {ref:42} max|diff| = {d:.3e}  "
              f"(whole signal incl. edges: {d_all:.3e})")
    fig.tight_layout(rect=(0, 0, 1, .96))
    out_png = os.path.join(os.path.abspath(folder), "vs_scipy.png")
    fig.savefig(out_png, dpi=170, bbox_inches="tight")
    shutil.rmtree(tmp, ignore_errors=True)
    print("wrote", out_png)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else
         os.path.dirname(os.path.abspath(__file__)))
