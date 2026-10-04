#!/usr/bin/env python3
"""

Rather than porting the DSP to numpy (which is what upsampling_time_alignment_test.py
does, and which can silently drift from the C++), this compiles filter_utils.hpp into
a shared library and calls filterutils::upsample through ctypes. The numbers in panel
A1 are therefore measurements of the shipping resampler, not of a reimplementation.

  A1  polyphase resample 256 -> 1000 Hz (125/32) of a 100 Hz sine.
      Overlays the ideal sine, the 256 Hz input, and the C++ 1000 Hz output, and
      reports measured frequency / amplitude / phase for input and output. The
      resampler compensates its own group delay internally (the interior kernel
      reads inPtr[baseInput + filterCenter - k]), so NO shift is applied here --
      matching output phase is the evidence that compensation is correct.

  A2  Each channel is upsampled to its own target rate (file_to_bin.hpp: "Every
      channel is upsampled to its OWN target rate, taken from the per-channel
      *_upsampled_rate columns of config.csv"). Sample instants are
      k * 1000/rate ms from the shared epoch origin (file_to_bin.cpp:431,438).
      This panel shows every modality's instants against the 1000 Hz master grid
      and reports the worst-case offset, computed exactly in rational arithmetic.

How the four C++ files are used:
  filter_utils.hpp   compiled + called via ctypes  (header-only, standalone)
  file_to_bin.hpp    parsed for ChannelIdx + header constants
  file_to_bin.cpp    parsed for the sample-instant formula
  main.cpp           not used: has its own main() and needs Qt

Usage:
  python3 acceptance_plot.py [--src DIR] [--out acceptance_plot.png] [--keep-build]
"""

from __future__ import annotations

import argparse
import ctypes
import os
import re
import shutil
import subprocess
import sys
import sysconfig
import tempfile
from fractions import Fraction
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Panel A1 stimulus. 100 Hz is deliberately close to the 128 Hz Nyquist of the
# 256 Hz input, so it exercises the antialiasing passband edge: the kernel uses
# fc = 1/max(P,Q) at the P*sourceRate intermediate rate, i.e. a 128 Hz cutoff.
# ---------------------------------------------------------------------------
SINE_HZ = 100.0
SINE_AMP = 1.0
SINE_PHASE = 0.70            # rad
SRC_RATE = 256.0
DST_RATE = 1000.0
DURATION_S = 4.0
FIT_WINDOW_S = (0.5, 3.5)    # interior only -- keeps filter edge transients out
VIEW_MS = (3000.0, 3040.0)   # 4 cycles of a 100 Hz sine

MASTER_RATE = 1000.0         # ECG target rate == the master grid
TOLERANCE_MS = 1000.0 / MASTER_RATE

# Per-channel upsample targets. These live in config.csv, which is not part of
# the C++ sources, so they are declared here; channel NAMES are validated
# against the ChannelIdx enum parsed out of file_to_bin.hpp.
MODALITIES = [
    ("ECG",             ["CH_ECG1", "CH_ECG2", "CH_ECG3"],            1000.0, "#1f77b4"),
    ("PPG",             ["CH_PPG"],                                    500.0, "#2ca02c"),
    ("pres/flow/snore", ["CH_PRES", "CH_FLOW", "CH_SNORE"],             32.0, "#e8730c"),
    ("accel",           ["CH_ACCEL_X", "CH_ACCEL_Y", "CH_ACCEL_Z"],    25.0, "#8b7fc7"),
    ("temp/pacemaker",  ["CH_TEMP", "CH_PACEMAKER_EVENT"],              1.0, "#d62728"),
]
A2_VIEW_MS = (-27.0, 67.0)   # negative span is label gutter, not data

SHIM = r"""
#include "filter_utils.hpp"
#include <cstring>
#include <vector>

#if defined(_WIN32)
#  define FU_API __declspec(dllexport)
#else
#  define FU_API __attribute__((visibility("default")))
#endif

extern "C" {

// gcd/P/Q exactly as filterutils::upsample derives them
FU_API int fu_pq(double srcRate, double tgtRate, int* P, int* Q) {
    int g = filterutils::greatest_common_divisor((int)tgtRate, (int)srcRate);
    if (g == 0) return -1;
    *P = (int)tgtRate / g;
    *Q = (int)srcRate / g;
    return 0;
}

FU_API int fu_out_len(int n, double srcRate, double tgtRate) {
    if (n <= 0) return 0;
    if (srcRate == tgtRate) return n;
    int P = 0, Q = 0;
    if (fu_pq(srcRate, tgtRate, &P, &Q) != 0) return -1;
    return (int)std::ceil((double)n * P / Q);
}

// Runs the real filterutils::upsample and copies the result out.
// Returns the number of samples written, or a negative error code.
FU_API int fu_upsample(const double* in, int n, double srcRate, double tgtRate,
                       double* out, int cap) {
    if (!in || !out || n < 0) return -1;
    try {
        std::vector<double> x(in, in + n);
        std::vector<double> y = filterutils::upsample(x, srcRate, tgtRate);
        if ((int)y.size() > cap) return -2;
        if (!y.empty()) std::memcpy(out, y.data(), y.size() * sizeof(double));
        return (int)y.size();
    } catch (const std::exception&) {
        return -3;
    } catch (...) {
        return -4;
    }
}

FU_API double fu_group_delay(double srcRate, double tgtRate) {
    return filterutils::group_delay_out_samples(srcRate, tgtRate);
}

}  // extern "C"
"""


# ===========================================================================
# Build + load the C++ resampler
# ===========================================================================
class CppResampler:
    """Compiles filter_utils.hpp into a shared library and binds it via ctypes."""

    def __init__(self, src_dir: Path, build_dir: Path):
        header = src_dir / "filter_utils.hpp"
        if not header.is_file():
            raise FileNotFoundError(f"filter_utils.hpp not found in {src_dir}")
        build_dir.mkdir(parents=True, exist_ok=True)
        shim = build_dir / "fu_shim.cpp"
        shim.write_text(SHIM)
        self.lib_path = build_dir / self._libname()
        self.compiler = self._compile(shim, src_dir, self.lib_path)

        lib = ctypes.CDLL(str(self.lib_path))
        dbl_p = ctypes.POINTER(ctypes.c_double)
        lib.fu_pq.argtypes = [ctypes.c_double, ctypes.c_double,
                              ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
        lib.fu_pq.restype = ctypes.c_int
        lib.fu_out_len.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double]
        lib.fu_out_len.restype = ctypes.c_int
        lib.fu_upsample.argtypes = [dbl_p, ctypes.c_int, ctypes.c_double,
                                    ctypes.c_double, dbl_p, ctypes.c_int]
        lib.fu_upsample.restype = ctypes.c_int
        lib.fu_group_delay.argtypes = [ctypes.c_double, ctypes.c_double]
        lib.fu_group_delay.restype = ctypes.c_double
        self._lib = lib

    @staticmethod
    def _libname() -> str:
        if sys.platform.startswith("win"):
            return "fu_shim.dll"
        if sys.platform == "darwin":
            return "libfu_shim.dylib"
        return "libfu_shim.so"

    @staticmethod
    def _compile(shim: Path, src_dir: Path, out: Path) -> str:
        """Try MSVC, then g++, then clang++. Returns the compiler used."""
        attempts = []

        if shutil.which("cl"):
            cmd = ["cl", "/nologo", "/LD", "/O2", "/EHsc", "/std:c++17",
                   f"/I{src_dir}", str(shim), f"/Fe:{out}",
                   f"/Fo:{out.parent / 'fu_shim.obj'}"]
            attempts.append(cmd)

        for cxx in ("g++", "clang++"):
            if shutil.which(cxx):
                cmd = [cxx, "-O2", "-std=c++17", "-shared", "-fPIC",
                       "-fvisibility=hidden", f"-I{src_dir}",
                       str(shim), "-o", str(out), "-pthread"]
                attempts.append(cmd)

        if not attempts:
            raise RuntimeError(
                "No C++ compiler found. Need one of: cl (MSVC), g++, clang++.\n"
                "On Windows, run this from a Developer Command Prompt so cl is on PATH."
            )

        errors = []
        for cmd in attempts:
            proc = subprocess.run(cmd, cwd=out.parent, capture_output=True, text=True)
            if proc.returncode == 0 and out.is_file():
                return Path(cmd[0]).name
            errors.append(f"$ {' '.join(cmd)}\n{proc.stdout}\n{proc.stderr}")
        raise RuntimeError("Could not build the resampler shim:\n\n" + "\n\n".join(errors))

    # -- API ----------------------------------------------------------------
    def pq(self, src_rate: float, tgt_rate: float) -> tuple[int, int]:
        P, Q = ctypes.c_int(0), ctypes.c_int(0)
        if self._lib.fu_pq(src_rate, tgt_rate, ctypes.byref(P), ctypes.byref(Q)) != 0:
            raise RuntimeError("fu_pq failed")
        return P.value, Q.value

    def group_delay_out_samples(self, src_rate: float, tgt_rate: float) -> float:
        return float(self._lib.fu_group_delay(src_rate, tgt_rate))

    def upsample(self, x: np.ndarray, src_rate: float, tgt_rate: float) -> np.ndarray:
        x = np.ascontiguousarray(x, dtype=np.float64)
        cap = self._lib.fu_out_len(x.size, src_rate, tgt_rate)
        if cap < 0:
            raise RuntimeError("fu_out_len failed")
        out = np.zeros(max(cap, 1), dtype=np.float64)
        dbl_p = ctypes.POINTER(ctypes.c_double)
        n = self._lib.fu_upsample(
            x.ctypes.data_as(dbl_p), x.size, src_rate, tgt_rate,
            out.ctypes.data_as(dbl_p), out.size)
        if n == -2:
            raise RuntimeError("output buffer too small (fu_out_len disagreed with C++)")
        if n < 0:
            raise RuntimeError(f"fu_upsample failed with code {n}")
        return out[:n]


# ===========================================================================
# Parse the C++ headers we cannot compile
# ===========================================================================
def parse_channel_enum(hpp: Path) -> dict[str, int]:
    """Extract ChannelIdx from file_to_bin.hpp as {name: index}."""
    text = hpp.read_text(errors="replace")
    m = re.search(r"enum\s+ChannelIdx\s*\{(.*?)\}\s*;", text, re.S)
    if not m:
        raise RuntimeError(f"ChannelIdx enum not found in {hpp}")
    body = re.sub(r"//.*?$|/\*.*?\*/", "", m.group(1), flags=re.S | re.M)
    out, nxt = {}, 0
    for item in (s.strip() for s in body.split(",")):
        if not item:
            continue
        if "=" in item:
            name, val = (p.strip() for p in item.split("=", 1))
            nxt = int(val, 0)
        else:
            name = item
        out[name] = nxt
        nxt += 1
    return out


def parse_constants(hpp: Path) -> dict[str, int]:
    """Extract the inline constexpr scalars from file_to_bin.hpp."""
    text = hpp.read_text(errors="replace")
    consts: dict[str, int] = {}
    for name in ("NUM_CHANNELS", "BIN_HEADER_VERSION"):
        m = re.search(rf"constexpr\s+\w+\s+{name}\s*=\s*([0-9]+)", text)
        if m:
            consts[name] = int(m.group(1))
    n = consts.get("NUM_CHANNELS")
    if n is not None:
        consts["NUM_HEADER_FIELDS"] = 4 + 4 * n
        consts["HEADER_SIZE"] = consts["NUM_HEADER_FIELDS"] * 4
    return consts


# ===========================================================================
# Measurement
# ===========================================================================
def fit_sine(t: np.ndarray, y: np.ndarray, f0: float) -> tuple[float, float, float]:
    """Least-squares fit of A*sin(2*pi*f*t + phi) + c, refining f around f0.

    Returns (frequency_hz, amplitude, phase_rad). Linear in (A, phi, c) for a
    fixed f, so f is refined by golden-section search on the residual.
    """
    def rss_and_params(f: float):
        w = 2.0 * np.pi * f * t
        A = np.column_stack([np.sin(w), np.cos(w), np.ones_like(t)])
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        resid = y - A @ coef
        amp = float(np.hypot(coef[0], coef[1]))
        phase = float(np.arctan2(coef[1], coef[0]))
        return float(resid @ resid), amp, phase

    lo, hi = f0 - 2.0, f0 + 2.0
    inv_phi = (np.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - inv_phi * (b - a), a + inv_phi * (b - a)
    fc, fd = rss_and_params(c)[0], rss_and_params(d)[0]
    for _ in range(90):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - inv_phi * (b - a)
            fc = rss_and_params(c)[0]
        else:
            a, c, fc = c, d, fd
            d = a + inv_phi * (b - a)
            fd = rss_and_params(d)[0]
        if b - a < 1e-12:
            break
    f_best = 0.5 * (a + b)
    _, amp, phase = rss_and_params(f_best)
    return f_best, amp, (phase + 2 * np.pi) % (2 * np.pi)


def worst_offset_ms(rate_hz: float, master_hz: float) -> float:
    """Exact worst-case distance from k*(1000/rate) ms to the nearest master tick.

    Sample instants are multiples of dt = 1000/rate ms; master ticks are multiples
    of T = 1000/master ms. With dt/T = p/q in lowest terms, {k*dt/T mod 1} is
    exactly the multiples of 1/q, so the worst distance is floor(q/2)/q periods.
    Exact rational arithmetic -- no horizon to choose, no sampling error.
    """
    ratio = Fraction(1000, 1) / Fraction(rate_hz).limit_denominator(10**6) / (
        Fraction(1000, 1) / Fraction(master_hz).limit_denominator(10**6))
    q = ratio.denominator
    return float(Fraction(q // 2, q) * Fraction(1000) / Fraction(master_hz).limit_denominator(10**6))


# ===========================================================================
# Figure
# ===========================================================================
def build_figure(cpp: CppResampler, channels: dict[str, int], out_png: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    P, Q = cpp.pq(SRC_RATE, DST_RATE)
    gd = cpp.group_delay_out_samples(SRC_RATE, DST_RATE)

    # ---- A1: run the real resampler -----------------------------------
    n_in = int(round(DURATION_S * SRC_RATE))
    t_in = np.arange(n_in) / SRC_RATE
    x_in = SINE_AMP * np.sin(2 * np.pi * SINE_HZ * t_in + SINE_PHASE)

    y_out = cpp.upsample(x_in, SRC_RATE, DST_RATE)
    t_out = np.arange(y_out.size) / DST_RATE

    lo_s, hi_s = FIT_WINDOW_S
    mi = (t_in >= lo_s) & (t_in <= hi_s)
    mo = (t_out >= lo_s) & (t_out <= hi_s)
    f_in, a_in, p_in = fit_sine(t_in[mi], x_in[mi], SINE_HZ)
    f_out, a_out, p_out = fit_sine(t_out[mo], y_out[mo], SINE_HZ)

    fig = plt.figure(figsize=(10, 7), dpi=130)
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1.5, 1.0],
                  left=0.089, right=0.967, top=0.945, bottom=0.075, hspace=0.42)

    ax1 = fig.add_subplot(gs[0])
    v0, v1 = VIEW_MS
    t_dense = np.linspace(v0 / 1000.0, v1 / 1000.0, 4000)
    ax1.plot(t_dense * 1000.0,
             SINE_AMP * np.sin(2 * np.pi * SINE_HZ * t_dense + SINE_PHASE),
             "-", color="#b6b6b6", lw=3.2, solid_capstyle="round",
             label=f"ideal {SINE_HZ:g} Hz sine", zorder=1)

    si = (t_in * 1000.0 >= v0) & (t_in * 1000.0 <= v1)
    so = (t_out * 1000.0 >= v0) & (t_out * 1000.0 <= v1)
    ax1.plot(t_out[so] * 1000.0, y_out[so], "o", ms=3.4, color="#2c6fbb",
             mew=0, label=f"upsampled @ {DST_RATE:g} Hz", zorder=3)
    ax1.plot(t_in[si] * 1000.0, x_in[si], "o", ms=9.5, color="#c0392b",
             mew=0, label=f"input @ {SRC_RATE:g} Hz", zorder=4)

    ax1.set_title(f"A1   polyphase resample {SRC_RATE:g} \u2192 {DST_RATE:g} Hz "
                  f"({P}/{Q}), {SINE_HZ:g} Hz sine", fontsize=11.5)
    ax1.set_xlabel("time (ms)   [zero-delay resampler: no shift applied]", fontsize=9.5)
    ax1.set_ylabel("amplitude", fontsize=9.5)
    ax1.set_xlim(v0, v1)
    ax1.grid(alpha=0.30, lw=0.6)
    ax1.tick_params(labelsize=8.5)

    # reorder so the legend reads ideal / input / upsampled like the reference
    h, l = ax1.get_legend_handles_labels()
    order = [l.index(f"ideal {SINE_HZ:g} Hz sine"),
             l.index(f"input @ {SRC_RATE:g} Hz"),
             l.index(f"upsampled @ {DST_RATE:g} Hz")]
    # framealpha=1.0: the box is opaque, so any sample instant falling behind it
    # is hidden rather than half-visible. The reference figure hides the 3031.25 ms
    # input marker this way.
    ax1.legend([h[i] for i in order], [l[i] for i in order],
               loc="upper right", fontsize=8.5, framealpha=1.0,
               borderpad=0.6, handletextpad=0.8, borderaxespad=0.5)

    # Left-align each measurement after "out", then pad to a fixed width so the
    # "(<src> Hz in ...)" column starts at the same offset on all three rows.
    def row(label: str, out_val: float, in_val: float, unit: str = "") -> str:
        suffix = f" {unit}" if unit else ""
        lhs = f"{DST_RATE:g} Hz out {out_val:.5f}{suffix}"
        return (f"{label:<11} {lhs:<28}"
                f"({SRC_RATE:g} Hz in {in_val:.5f}{suffix})")

    box = "\n".join([
        row("FREQUENCY", f_out, f_in, "Hz"),
        row("AMPLITUDE", a_out, a_in),
        row("PHASE", p_out, p_in, "rad"),
    ])
    ax1.text(0.012, 0.035, box, transform=ax1.transAxes, fontsize=8.0,
             family="monospace", va="bottom", ha="left", zorder=5,
             bbox=dict(boxstyle="square,pad=0.45", facecolor="white",
                       edgecolor="#9a9a9a", lw=0.9))

    # ---- A2: sample-instant lattice vs the 1000 Hz master grid ---------
    ax2 = fig.add_subplot(gs[1])
    a0, a1 = A2_VIEW_MS
    rows = list(reversed(MODALITIES))          # first entry drawn at the top
    worst = 0.0
    for i, (label, names, rate, color) in enumerate(rows):
        missing = [n for n in names if n not in channels]
        if missing:
            raise RuntimeError(f"{label}: not in ChannelIdx: {', '.join(missing)}")
        dt = 1000.0 / rate
        k_max = int(np.floor(a1 / dt))
        instants = np.arange(0, k_max + 1) * dt
        instants = instants[(instants >= 0.0) & (instants <= a1)]
        ax2.eventplot(instants, lineoffsets=i, linelengths=0.62,
                      colors=color, linewidths=1.7)
        ax2.text(-2.0, i, f"{label} ({rate:g} Hz)", color=color, fontsize=9,
                 ha="right", va="center")
        worst = max(worst, worst_offset_ms(rate, MASTER_RATE))

    for tick in np.arange(0.0, np.floor(a1) + 1.0, 1000.0 / MASTER_RATE):
        ax2.axvline(tick, color="#cccccc", lw=0.55, zorder=0)

    ok = "\u2264" if worst <= TOLERANCE_MS else ">"
    ax2.set_title(f"A2   every modality sample lands within one {MASTER_RATE:g} Hz period "
                  f"(max offset {worst:.5f} ms {ok} {TOLERANCE_MS:.5f} ms)", fontsize=11.5)
    ax2.set_xlabel(f"time since shared epoch origin (ms)   "
                   f"[gridlines = {MASTER_RATE:g} Hz master ticks]", fontsize=9.5)
    ax2.set_xlim(a0, a1)
    ax2.set_ylim(-0.75, len(rows) - 0.25)
    ax2.set_yticks([])
    ax2.tick_params(labelsize=8.5)
    for side in ("left", "right", "top"):
        ax2.spines[side].set_visible(True)

    fig.savefig(out_png, dpi=130,
                metadata={"Software": f"matplotlib via acceptance_plot.py "
                                      f"(resampler: C++ filter_utils.hpp, {cpp.compiler})"})
    plt.close(fig)

    return out_png, dict(P=P, Q=Q, group_delay=gd, worst=worst,
                         f_in=f_in, a_in=a_in, p_in=p_in,
                         f_out=f_out, a_out=a_out, p_out=p_out,
                         n_in=n_in, n_out=int(y_out.size))


# ===========================================================================
def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", type=Path, default=Path(__file__).resolve().parent,
                    help="directory holding filter_utils.hpp and file_to_bin.hpp")
    ap.add_argument("--out", type=Path, default=Path("acceptance_plot.png"))
    ap.add_argument("--keep-build", action="store_true",
                    help="keep the compiled shim instead of using a temp dir")
    args = ap.parse_args(argv[1:])

    src = args.src.resolve()
    print("=== acceptance plot: A1 resampler fidelity, A2 sample-instant lattice ===")
    print(f"C++ sources: {src}")

    consts = parse_constants(src / "file_to_bin.hpp")
    channels = parse_channel_enum(src / "file_to_bin.hpp")
    print(f"  file_to_bin.hpp  -> {len(channels)} channels in ChannelIdx, "
          f"NUM_CHANNELS={consts.get('NUM_CHANNELS')}, "
          f"HEADER_SIZE={consts.get('HEADER_SIZE')} B, "
          f"version={consts.get('BIN_HEADER_VERSION')}")
    if consts.get("NUM_CHANNELS") != len(channels):
        print(f"  WARNING: NUM_CHANNELS={consts.get('NUM_CHANNELS')} but the enum "
              f"lists {len(channels)} entries")

    build_dir = (args.out.resolve().parent / "_fu_build" if args.keep_build
                 else Path(tempfile.mkdtemp(prefix="fu_build_")))
    try:
        cpp = CppResampler(src, build_dir)
        print(f"  filter_utils.hpp -> {cpp.lib_path.name} via {cpp.compiler}")

        out_png, info = build_figure(cpp, channels, args.out.resolve())

        print(f"\nA1  {SRC_RATE:g} -> {DST_RATE:g} Hz is P/Q = {info['P']}/{info['Q']}; "
              f"{info['n_in']} in -> {info['n_out']} out")
        print(f"    group delay reported by the C++ = {info['group_delay']:.5f} output "
              f"samples, already compensated internally, so no shift is applied")
        print(f"    frequency  in {info['f_in']:.5f} Hz   out {info['f_out']:.5f} Hz "
              f"(expected {SINE_HZ:.5f})")
        print(f"    amplitude  in {info['a_in']:.5f}      out {info['a_out']:.5f} "
              f"(expected {SINE_AMP:.5f})")
        print(f"    phase      in {info['p_in']:.5f} rad  out {info['p_out']:.5f} rad "
              f"(expected {SINE_PHASE:.5f})")
        dphi = abs(info["p_out"] - info["p_in"])
        print(f"    phase error out vs in = {dphi:.2e} rad "
              f"= {dphi / (2 * np.pi / SINE_HZ) * 1e3:.4f} us")

        print(f"\nA2  worst-case offset from the {MASTER_RATE:g} Hz master grid = "
              f"{info['worst']:.5f} ms (tolerance {TOLERANCE_MS:.5f} ms) -> "
              f"{'PASS' if info['worst'] <= TOLERANCE_MS else 'FAIL'}")
        for label, _, rate, _ in MODALITIES:
            w = worst_offset_ms(rate, MASTER_RATE)
            print(f"      {label:<16} {rate:>7.2f} Hz  dt = {1000.0 / rate:>8.4f} ms  "
                  f"worst offset {w:.5f} ms")

        print(f"\nwrote {out_png}")
        return 0
    finally:
        if not args.keep_build:
            shutil.rmtree(build_dir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
