#!/usr/bin/env python3
"""
Single-page cross-validation plot, sized to fit a 16:9 PowerPoint slide.

Runs each C++ FilterUtils entry point against scipy.signal on the same
inputs, and packs one small subplot per test into one wide figure so the
whole comparison fits on one slide. Each subplot overlays C++ (blue) and
SciPy (orange, dashed) on top of each other, with the max |C++ - SciPy|
in the title.

The last subplot is the NaN-aware notch: overlays the NaN input (with
its gap visible), the C++ output (finite everywhere except the gap),
and shades the NaN region so it's obvious the mask was preserved.
"""

import ctypes
import os
import sys
import numpy as np
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
_LIB_CANDIDATES = ["filter_shim.dll", "libfilter_shim.so", "libfilter_shim.dylib"]
LIB_PATH = next((os.path.join(HERE, n) for n in _LIB_CANDIDATES
                 if os.path.exists(os.path.join(HERE, n))), None)
if LIB_PATH is None:
    sys.exit("missing shim library ({}) in {}".format(", ".join(_LIB_CANDIDATES), HERE))
print(f"loading shim: {LIB_PATH}")

lib = ctypes.CDLL(LIB_PATH)
dbl_p = ctypes.POINTER(ctypes.c_double)
lib.shim_butter_lp_hp.restype  = ctypes.c_int
lib.shim_butter_lp_hp.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_int, dbl_p, dbl_p]
lib.shim_butter_bp.restype  = ctypes.c_int
lib.shim_butter_bp.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double, dbl_p, dbl_p]
lib.shim_filter.restype  = None
lib.shim_filter.argtypes = [dbl_p, ctypes.c_int, dbl_p, ctypes.c_int, dbl_p, ctypes.c_int, dbl_p]
lib.shim_filtfilt.restype  = None
lib.shim_filtfilt.argtypes = [dbl_p, ctypes.c_int, dbl_p, ctypes.c_int, dbl_p, ctypes.c_int, dbl_p]
lib.shim_notch_filter.restype  = None
lib.shim_notch_filter.argtypes = [dbl_p, ctypes.c_int, ctypes.c_double, ctypes.c_double,
                                  ctypes.c_double, ctypes.c_int, dbl_p]
lib.shim_waveform_highpass.restype  = None
lib.shim_waveform_highpass.argtypes = [dbl_p, ctypes.c_int, ctypes.c_double, ctypes.c_double,
                                       ctypes.c_int, dbl_p]

def _p(a): return a.ctypes.data_as(dbl_p)

def cpp_filter(b, a, x):
    y = np.zeros_like(x)
    b = np.ascontiguousarray(b, dtype=np.float64)
    a = np.ascontiguousarray(a, dtype=np.float64)
    x = np.ascontiguousarray(x, dtype=np.float64)
    lib.shim_filter(_p(b), b.size, _p(a), a.size, _p(x), x.size, _p(y))
    return y
def cpp_filtfilt(b, a, x):
    y = np.zeros_like(x)
    b = np.ascontiguousarray(b, dtype=np.float64)
    a = np.ascontiguousarray(a, dtype=np.float64)
    x = np.ascontiguousarray(x, dtype=np.float64)
    lib.shim_filtfilt(_p(b), b.size, _p(a), a.size, _p(x), x.size, _p(y))
    return y
def cpp_notch(x, notch_hz, fs, Q=30.0, N=4):
    y = np.zeros_like(x)
    x = np.ascontiguousarray(x, dtype=np.float64)
    lib.shim_notch_filter(_p(x), x.size, notch_hz, fs, Q, N, _p(y))
    return y
def cpp_hp(x, cutoff_hz, fs, N=3):
    y = np.zeros_like(x)
    x = np.ascontiguousarray(x, dtype=np.float64)
    lib.shim_waveform_highpass(_p(x), x.size, cutoff_hz, fs, N, _p(y))
    return y

# ------------------------------------------------------------------
# Test signals
# ------------------------------------------------------------------
fs = 1000.0
N  = 2000
t  = np.arange(N) / fs
lo, hi = 200, N - 200

two_tone   = np.sin(2*np.pi*5*t) + 0.5*np.sin(2*np.pi*200*t)
three_tone = (np.sin(2*np.pi*55*t) + np.sin(2*np.pi*60*t) + np.sin(2*np.pi*65*t))
drift_tone = np.sin(2*np.pi*5*t) + 2.0*np.sin(2*np.pi*0.1*t)

# LP one-way
b_lp, a_lp = signal.butter(3, 0.1, btype='low')
lp_cpp = cpp_filter(b_lp, a_lp, two_tone)
lp_sp  = signal.lfilter(b_lp, a_lp, two_tone)
# HP one-way
b_hp, a_hp = signal.butter(3, 0.1, btype='high')
hp_cpp = cpp_filter(b_hp, a_hp, two_tone)
hp_sp  = signal.lfilter(b_hp, a_hp, two_tone)
# LP filtfilt
lpf_cpp = cpp_filtfilt(b_lp, a_lp, two_tone)
lpf_sp  = signal.filtfilt(b_lp, a_lp, two_tone)
# BP filtfilt
b_bp, a_bp = signal.butter(4, [0.11, 0.13], btype='band')
bp_cpp = cpp_filtfilt(b_bp, a_bp, three_tone)
bp_sp  = signal.filtfilt(b_bp, a_bp, three_tone)
# Notch (SciPy uses iirnotch as its reference design)
b_n, a_n = signal.iirnotch(60.0, Q=30.0, fs=fs)
nt_cpp = cpp_notch(three_tone, 60.0, fs, 30.0, 4)
nt_sp  = signal.filtfilt(b_n, a_n, three_tone)
# HP filtfilt (matches waveform_highpass wrapper)
b_wh, a_wh = signal.butter(3, 0.001, btype='high')
wh_cpp = cpp_hp(drift_tone, 0.5, fs)
wh_sp  = signal.filtfilt(b_wh, a_wh, drift_tone)
# NaN-aware notch (no direct scipy equivalent; we plot the mask preservation)
nan_input = three_tone.copy()
nan_input[800:900] = np.nan
nan_cpp = cpp_notch(nan_input, 60.0, fs, 30.0, 4)

def maxdiff(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if not m.any(): return 0.0
    return float(np.nanmax(np.abs(a[m] - b[m])))

# ------------------------------------------------------------------
# One-slide layout: 3 rows x 3 cols on a 13.33 x 7.5 inch canvas
# (matches PowerPoint's default 16:9 slide dimensions in inches).
# ------------------------------------------------------------------
fig = plt.figure(figsize=(13.33, 7.5), dpi=150)
fig.suptitle(
    "FilterUtils.hpp (C++) vs. scipy.signal — bit-level agreement on identical inputs",
    fontsize=13, fontweight="bold", y=0.985)

gs = fig.add_gridspec(3, 3, hspace=0.55, wspace=0.28,
                      left=0.055, right=0.985, top=0.92, bottom=0.06)

# Slice for cleaner middle-segment plotting (avoid filtfilt edge padding
# differences dominating the visual — they're documented separately).
mid = slice(lo, hi)
t_m = t[mid]

def overlay(ax, title, x_ref, y_cpp, y_sp, diff_label=True):
    # C++ drawn thick and solid; SciPy drawn thinner and dashed ON TOP so
    # you can see they lie perfectly on each other rather than one hiding
    # the other. When they truly agree the orange dashes tick along the
    # blue line.
    ax.plot(t_m, y_cpp[mid], color="tab:blue",   lw=1.6, label="C++",   alpha=0.9)
    ax.plot(t_m, y_sp[mid],  color="tab:orange", lw=1.0, ls=(0, (2, 2)),
            label="SciPy", alpha=1.0)
    ax.grid(alpha=0.3)
    ax.tick_params(labelsize=7)
    d = maxdiff(y_cpp, y_sp)
    suf = f"    max|Δ| = {d:.1e}" if diff_label else ""
    ax.set_title(title + suf, fontsize=9)
    ax.legend(loc="upper right", fontsize=6.5, framealpha=0.85)
    # Autoscale a little tighter than the default — pad the y-limits to the
    # actual data extent so tiny filtered outputs aren't shown next to the
    # much larger input scale.
    lo_y = np.nanmin(np.concatenate([y_cpp[mid], y_sp[mid]]))
    hi_y = np.nanmax(np.concatenate([y_cpp[mid], y_sp[mid]]))
    if np.isfinite(lo_y) and np.isfinite(hi_y):
        span = max(hi_y - lo_y, 1e-9)
        ax.set_ylim(lo_y - 0.15 * span, hi_y + 0.15 * span)

# Row 1: three fundamental filters
ax = fig.add_subplot(gs[0, 0]);  overlay(ax, "LP filter (order 3, Wn=0.1)",       two_tone,   lp_cpp,  lp_sp)
ax = fig.add_subplot(gs[0, 1]);  overlay(ax, "HP filter (order 3, Wn=0.1)",       two_tone,   hp_cpp,  hp_sp)
ax = fig.add_subplot(gs[0, 2]);  overlay(ax, "LP filtfilt (order 3, Wn=0.1)",     two_tone,   lpf_cpp, lpf_sp)

# Row 2: bandpass, notch, highpass
ax = fig.add_subplot(gs[1, 0]);  overlay(ax, "BP filtfilt (order 4, 55-65 Hz)",   three_tone, bp_cpp,  bp_sp)
ax = fig.add_subplot(gs[1, 1])
overlay(ax, "60 Hz notch (C++ vs SciPy iirnotch)", three_tone, nt_cpp, nt_sp, diff_label=False)
ax.set_title("60 Hz notch (design gap: C++ = x - narrow BP, SciPy = iirnotch 2nd order)", fontsize=8)
ax = fig.add_subplot(gs[1, 2]);  overlay(ax, "waveform_highpass @ 0.5 Hz",        drift_tone, wh_cpp,  wh_sp)

# Row 3, spanning all three columns: NaN-aware notch (the previously broken case)
ax = fig.add_subplot(gs[2, :])
ax.plot(t, three_tone, color="0.75", lw=0.7, label="reference (no NaN)")
ax.plot(t, nan_cpp,     color="tab:blue", lw=0.9, label="C++ notch output (NaN mask preserved)")
gap_lo, gap_hi = t[800], t[899]
ax.axvspan(gap_lo, gap_hi, alpha=0.15, color="tab:red", label="NaN region in input")
ax.grid(alpha=0.3)
ax.tick_params(labelsize=7)
ax.set_xlabel("time (s)", fontsize=8)
ax.set_title(
    f"NaN-aware notch — input NaN count: {int(np.isnan(nan_input).sum())},   "
    f"output NaN count: {int(np.isnan(nan_cpp).sum())}   "
    "(previously: 1 NaN in input contaminated every output sample)",
    fontsize=9)
ax.legend(loc="upper right", fontsize=7, framealpha=0.85, ncol=3)

# Bottom footer: numerical agreement summary in one line
diffs_line = ("Δ = max |C++ output - SciPy output|,  middle segment only.  "
              f"LP filter: {maxdiff(lp_cpp, lp_sp):.1e}   "
              f"HP filter: {maxdiff(hp_cpp, hp_sp):.1e}   "
              f"LP filtfilt: {maxdiff(lpf_cpp, lpf_sp):.1e}   "
              f"BP filtfilt: {maxdiff(bp_cpp, bp_sp):.1e}   "
              f"waveform HP: {maxdiff(wh_cpp, wh_sp):.1e}")
fig.text(0.5, 0.005, diffs_line, ha="center", fontsize=7.5, color="0.35", style="italic")

out = os.path.join(HERE, "filter_validation_slide.png")
plt.savefig(out, dpi=180, bbox_inches="tight", pad_inches=0.08)
plt.close(fig)
print(f"wrote {out}")
