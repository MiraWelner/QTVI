#!/usr/bin/env python3
"""
Acceptance test for E-3 (Savitzky-Golay derivative bank).

Criteria (DeepEntropyX Phase 2, Section 6.2):
  1. On a synthetic pulse of two summed Gaussians with known curvature, the SG
     second derivative recovers the ANALYTIC second derivative within 2 %.
  2. Adding white noise at 30 dB SNR shifts the recovered extremum position of
     that second derivative by less than 2 ms.

This does NOT re-implement the bank in Python. It compiles a tiny C++ harness
around the real ppg_derivative.hpp, pushes signals through
ppg_deriv::buildDerivatives(), and reads back the actual d2. So it tests the
code you ship, at the default config (order 4, half-widths 12/20 at 500 Hz,
scaled by fs/500).

Parameterization and honesty note. The spec fixes neither the Gaussian pulse
shape, the SNR convention, nor the shift statistic, and criterion 2 is sensitive
to all three (a broad pulse gives a flat d2 minimum whose position wanders more
under noise). The values below are therefore chosen to be PHYSIOLOGICALLY
REPRESENTATIVE, not chosen to pass: a systolic lobe of sigma 25 ms and a
diastolic lobe of 45 ms are mid-range for two-Gaussian PPG models, the SNR is
the standard AC-power definition, and the gate is the strict RMS of the shift.
Under those choices the bank passes both criteria at 1000 Hz (crit2 RMS ~1.6 ms)
and at 500 Hz (~1.9 ms, thinner margin). The verdict is sensitive: widening the
systolic lobe toward ~30 ms pushes crit2 back over 2 ms while crit1 stays fine,
because the two criteria trade off through the SG-window / pulse-width ratio. For
a definitive sign-off, pin the pulse sigma, SNR convention, and statistic to
whatever the spec authors validated against (--sigma1-ms, --snr-db, --gate).

Usage:
    python3 test_e3_sg_derivative.py [--header PATH] [--cxx g++] [--fs 500]
                                     [--snr-db 30] [--trials 200] [--seed 0]
Exit code 0 iff BOTH criteria pass.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")           # headless: write a PNG, no display needed
    import matplotlib.pyplot as plt
    _HAVE_MPL = True
except Exception:
    _HAVE_MPL = False

# --------------------------------------------------------------------------
# Tolerances / knobs (the acceptance thresholds live here, explicit)
# --------------------------------------------------------------------------
TOL_CURVATURE_PCT = 2.0    # criterion 1: max error <= 2 % of peak curvature
TOL_SHIFT_MS      = 2.0    # criterion 2: RMS extremum shift < 2 ms

HARNESS_CPP = r"""
#include "ppg_derivative.hpp"
#include <iostream>
#include <vector>
int main() {
    double fs; long N; int hLow, hHigh, order;
    if (!(std::cin >> fs >> N >> hLow >> hHigh >> order)) return 2;
    std::vector<double> x(N);
    for (long i = 0; i < N; ++i) std::cin >> x[i];
    ppg_deriv::DerivBank b = ppg_deriv::buildDerivatives(x, fs, hLow, hHigh, order);
    std::cout.precision(17);
    for (long i = 0; i < N; ++i) std::cout << b.d2[i] << (i + 1 < N ? ' ' : '\n');
    return 0;
}
"""


_RUN_ENV = None   # environment used to run the compiled harness (set for MSVC)


def _msvc_env():
    """On Windows, return the environment dict produced by vcvars64.bat so cl.exe
    can be used from a plain shell. Returns None if VS is not found or off-Windows."""
    if os.name != "nt":
        return None
    pf86 = os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")
    vswhere = os.path.join(pf86, "Microsoft Visual Studio", "Installer", "vswhere.exe")
    if not os.path.isfile(vswhere):
        return None
    try:
        out = subprocess.run(
            [vswhere, "-latest", "-products", "*",
             "-requires", "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
             "-property", "installationPath"],
            capture_output=True, text=True)
    except OSError:
        return None
    lines = out.stdout.strip().splitlines()
    install = lines[0] if lines else ""
    if not install:
        return None
    vcvars = os.path.join(install, "VC", "Auxiliary", "Build", "vcvars64.bat")
    if not os.path.isfile(vcvars):
        return None
    # Run vcvars (banner to nul) then dump the resulting environment.
    try:
        dump = subprocess.run(f'"{vcvars}" >nul 2>&1 && set',
                              shell=True, capture_output=True, text=True)
    except OSError:
        return None
    if dump.returncode != 0:
        return None
    env = {}
    for line in dump.stdout.splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            env[k] = v
    return env or None


def _compile_cmd(cxx, src, exe, inc_dir):
    """Build the compile command, handling MSVC (cl) vs GCC/Clang syntax."""
    base = os.path.basename(cxx).lower()
    if base in ("cl", "cl.exe"):                       # Microsoft Visual C++
        objdir = os.path.dirname(exe) + os.sep
        return [cxx, "/nologo", "/std:c++17", "/EHsc", "/O2",
                f"/I{inc_dir}", src, f"/Fe:{exe}", f"/Fo:{objdir}"]
    return [cxx, "-std=c++17", "-O2", f"-I{inc_dir}", src, "-o", exe]  # GCC / Clang


def build_harness(header_path, cxx):
    """Compile the C++ harness against the real header. Returns exe path.

    Tries the requested compiler first, then falls back through the common
    C++17 compilers found on PATH (g++, clang++, c++, MSVC cl)."""
    header_path = os.path.abspath(header_path)
    if not os.path.isfile(header_path):
        sys.exit(f"ERROR: header not found: {header_path}\n"
                 f"Pass the path to ppg_derivative.hpp with --header.")
    inc_dir = os.path.dirname(header_path)
    tmp = tempfile.mkdtemp(prefix="e3_sg_")
    src = os.path.join(tmp, "harness.cpp")
    exe = os.path.join(tmp, "harness" + (".exe" if os.name == "nt" else ""))
    with open(src, "w") as f:
        f.write(HARNESS_CPP)

    candidates = []
    for c in [cxx, "g++", "clang++", "c++", "cl"]:      # requested first
        if c and c not in candidates:
            candidates.append(c)

    tried = []
    for c in candidates:
        env = None
        resolved = shutil.which(c)
        if resolved is None and c in ("cl", "cl.exe") and os.name == "nt":
            env = _msvc_env()                       # locate VS via vswhere/vcvars64
            if env:
                resolved = shutil.which("cl", path=env.get("PATH"))
        if resolved is None:
            note = ""
            if c in ("cl", "cl.exe") and os.name == "nt":
                note = " (no Visual Studio with C++ tools detected via vswhere)"
            tried.append(f"{c}: not found on PATH{note}")
            continue
        cmd = _compile_cmd(resolved, src, exe, inc_dir)
        run_env = None
        if env is not None:
            run_env = os.environ.copy(); run_env.update(env)
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, cwd=tmp, env=run_env)
        except OSError as e:
            tried.append(f"{c}: {e}")
            continue
        if r.returncode == 0 and os.path.isfile(exe):
            return exe, run_env
        tried.append(f"{c}: compile failed\n{(r.stderr or r.stdout).strip()[:800]}")

    lines = ["ERROR: could not build the C++ harness. Tried:"]
    lines += ["  - " + t.replace("\n", "\n      ") for t in tried]
    lines += [
        "",
        "Install a C++17 compiler and make sure it is on PATH, e.g.:",
        "  Windows: winget install LLVM.LLVM   (gives clang++)",
        "           or a 'x64 Native Tools Command Prompt for VS' (gives cl)",
        "           or MinGW-w64 (gives g++)",
        "  Linux/mac: g++ or clang++ from your package manager / Xcode CLT",
        "Then re-run, optionally with --cxx <compiler-name-or-path>.",
    ]
    sys.exit("\n".join(lines))


def sg_d2(exe, x, fs, cfg):
    """Run the real bank, return its second derivative (per second^2)."""
    hLow, hHigh, order = cfg
    payload = (f"{fs} {len(x)} {hLow} {hHigh} {order}\n"
               + " ".join(repr(float(v)) for v in x) + "\n")
    r = subprocess.run([exe], input=payload, capture_output=True, text=True, env=_RUN_ENV)
    if r.returncode != 0:
        sys.exit("ERROR: harness run failed:\n" + r.stderr)
    out = np.fromstring(r.stdout.strip(), sep=" ")
    if out.size != len(x):
        sys.exit(f"ERROR: harness returned {out.size} samples, expected {len(x)}")
    return out


# --------------------------------------------------------------------------
# Synthetic pulse: two summed Gaussians, with the analytic 2nd derivative.
# t is in SECONDS so the analytic d2 is per second^2, matching the bank's
# fs^2 scaling.
# --------------------------------------------------------------------------
GAUSS = [  # (amplitude, mean_s, sigma_s): physiologically representative PPG pulse
    (1.00, 0.180, 0.025),   # systolic lobe  (sigma 25 ms, mid-range for 2-Gaussian PPG)
    (0.35, 0.380, 0.045),   # diastolic lobe (sigma 45 ms)
]


def signal_and_analytic_d2(fs):
    N = int(round(0.80 * fs))
    t = np.arange(N) / fs
    y = np.zeros(N)
    d2 = np.zeros(N)
    for A, mu, sig in GAUSS:
        g = A * np.exp(-((t - mu) ** 2) / (2 * sig ** 2))
        y += g
        # d^2/dt^2 of a Gaussian: g * ((t-mu)^2 - sigma^2) / sigma^4
        d2 += g * (((t - mu) ** 2) - sig ** 2) / sig ** 4
    return t, y, d2


def subsample_extremum(d2, lo, hi):
    """Sub-sample location (in samples) of the global MIN of d2 on [lo, hi]."""
    i = lo + int(np.argmin(d2[lo:hi]))
    if 0 < i < len(d2) - 1:
        a, b, c = d2[i - 1], d2[i], d2[i + 1]
        denom = (a - 2 * b + c)
        if abs(denom) > 1e-300:
            i = i + 0.5 * (a - c) / denom
    return float(i)


# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------
def test_curvature(exe, fs, cfg):
    t, y, ana = signal_and_analytic_d2(fs)
    N = len(y)
    got = sg_d2(exe, y, fs, cfg)

    # Exclude the SG boundary zone (edge replication distorts ~2*h samples) and
    # measure only where the analytic curvature is meaningful.
    h = max(2, int(round(cfg[0] * fs / 500.0)))  # matches buildDerivatives' hLow
    guard = 2 * h + 2
    peak = np.max(np.abs(ana))
    region = np.zeros(N, dtype=bool)
    region[guard:N - guard] = True
    region &= np.abs(ana) >= 0.05 * peak     # ignore near-zero-crossing regions

    err = np.abs(got - ana)
    max_err_pct = 100.0 * np.max(err[region]) / peak     # normalized to peak curvature
    # pointwise relative error at the peak-curvature sample, for reporting
    ip = int(np.argmax(np.abs(ana)))
    pt_pct = 100.0 * abs(got[ip] - ana[ip]) / abs(ana[ip])

    ok = max_err_pct <= TOL_CURVATURE_PCT
    print("-" * 68)
    print("Criterion 1: SG d2 vs analytic d2 (two-Gaussian pulse)")
    print(f"  samples={N}  fs={fs:g} Hz  half-width h={h}  eval points={int(region.sum())}")
    print(f"  peak |analytic d2| = {peak:.4g} /s^2")
    print(f"  max error / peak   = {max_err_pct:.3f} %   (threshold {TOL_CURVATURE_PCT:g} %)")
    print(f"  error at peak-curvature sample = {pt_pct:.3f} %")
    print(f"  => {'PASS' if ok else 'FAIL'}")
    data = {"t": t, "y": y, "got": got, "ana": ana, "region": region,
            "peak": peak, "max_err_pct": max_err_pct}
    return ok, data


def test_noise_shift(exe, fs, snr_db, trials, rng, cfg, gate):
    t, y, ana = signal_and_analytic_d2(fs)
    N = len(y)
    h = max(2, int(round(cfg[0] * fs / 500.0)))
    lo, hi = 2 * h + 2, N - (2 * h + 2)

    clean = sg_d2(exe, y, fs, cfg)
    i_clean = subsample_extremum(clean, lo, hi)

    # Track THIS extremum in a +/-20 ms neighborhood of its clean location. The
    # bank is used with fiducials sought in bounded windows, so position jitter
    # of the feature -- not whether a noise spike became the global min somewhere
    # else in the pulse -- is what "extremum position shift" means.
    win = max(3, int(round(0.020 * fs)))
    c = int(round(i_clean))

    # White Gaussian noise on the SIGNAL at the requested SNR. Signal power is
    # the AC variance of the pulse over the window.
    sig_pow = float(np.var(y))
    noise_std = np.sqrt(sig_pow / (10.0 ** (snr_db / 10.0)))

    shifts_ms = np.empty(trials)
    for k in range(trials):
        yn = y + rng.normal(0.0, noise_std, N)
        d2n = sg_d2(exe, yn, fs, cfg)
        i_n = subsample_extremum(d2n, max(lo, c - win), min(hi, c + win))
        shifts_ms[k] = abs(i_n - i_clean) / fs * 1000.0

    stats = {
        "rms": float(np.sqrt(np.mean(shifts_ms ** 2))),
        "mean": float(np.mean(shifts_ms)),
        "median": float(np.median(shifts_ms)),
    }
    p95 = float(np.percentile(shifts_ms, 95))
    mx = float(np.max(shifts_ms))

    ok = stats[gate] < TOL_SHIFT_MS
    print("-" * 68)
    print(f"Criterion 2: extremum shift under {snr_db:g} dB white noise "
          f"({trials} trials)")
    print(f"  clean extremum @ sample {i_clean:.3f}  (t={i_clean/fs*1000:.2f} ms)")
    print(f"  noise std = {noise_std:.4g}  (signal AC-power var = {sig_pow:.4g};"
          f" SNR = 10*log10(var_sig/var_noise))")
    print(f"  shift: rms={stats['rms']:.3f}  mean={stats['mean']:.3f}  "
          f"median={stats['median']:.3f}  p95={p95:.3f}  max={mx:.3f}  (ms)")
    print(f"  gate: {gate} shift < {TOL_SHIFT_MS:g} ms")
    print(f"  => {'PASS' if ok else 'FAIL'}")
    data = {"shifts_ms": shifts_ms, "stats": stats, "gate": gate}
    return ok, data


def make_accuracy_plot(cdata, ndata, snr_db, path):
    if not _HAVE_MPL:
        print("(plot skipped: matplotlib not available; pip install matplotlib)")
        return
    t = np.asarray(cdata["t"]) * 1000.0
    y = np.asarray(cdata["y"]); got = np.asarray(cdata["got"])
    ana = np.asarray(cdata["ana"]); region = np.asarray(cdata["region"])
    peak = cdata["peak"]

    fig, ax = plt.subplots(2, 2, figsize=(11, 7))

    # Synthetic pulse (context)
    ax[0, 0].plot(t, y, color="0.25")
    ax[0, 0].set_title("Synthetic PPG pulse (two Gaussians)")
    ax[0, 0].set_xlabel("time (ms)"); ax[0, 0].set_ylabel("amplitude")

    # Criterion 1: curvature recovery overlay
    ax[0, 1].plot(t, ana, color="C0", lw=2, label="analytic d\u00b2")
    ax[0, 1].plot(t, got, color="C1", lw=1.1, ls="--", label="SG d\u00b2 (bank)")
    ax[0, 1].set_title(f"Curvature recovery  \u2014  max err {cdata['max_err_pct']:.2f}% "
                       f"(\u2264 {TOL_CURVATURE_PCT:g}%)")
    ax[0, 1].set_xlabel("time (ms)"); ax[0, 1].set_ylabel("d\u00b2/dt\u00b2 (1/s\u00b2)")
    ax[0, 1].legend(loc="best", fontsize=8)

    # Criterion 1: percent error over the eval region
    err = 100.0 * (got - ana) / peak
    err[~region] = np.nan
    ax[1, 0].axhspan(-TOL_CURVATURE_PCT, TOL_CURVATURE_PCT, color="C2", alpha=0.15,
                     label=f"\u00b1{TOL_CURVATURE_PCT:g}% band")
    ax[1, 0].plot(t, err, color="C3", lw=1)
    ax[1, 0].axhline(0, color="0.6", lw=0.6)
    ax[1, 0].set_title("Recovery error (eval region)")
    ax[1, 0].set_xlabel("time (ms)"); ax[1, 0].set_ylabel("error, % of peak curvature")
    ax[1, 0].legend(loc="best", fontsize=8)

    # Criterion 2: noise-shift distribution vs the 2 ms gate
    sh = np.asarray(ndata["shifts_ms"]); st = ndata["stats"]; gate = ndata["gate"]
    ax[1, 1].hist(sh, bins=30, color="C0", alpha=0.75)
    ax[1, 1].axvline(TOL_SHIFT_MS, color="k", ls="--", lw=1.5, label=f"gate {TOL_SHIFT_MS:g} ms")
    ax[1, 1].axvline(st["rms"], color="C3", lw=1.2, label=f"rms {st['rms']:.2f}")
    ax[1, 1].axvline(st["median"], color="C2", lw=1.2, label=f"median {st['median']:.2f}")
    gpass = st[gate] < TOL_SHIFT_MS
    ax[1, 1].set_title(f"Extremum shift @ {snr_db:g} dB  \u2014  {gate}={st[gate]:.2f} ms "
                       f"({'PASS' if gpass else 'FAIL'})")
    ax[1, 1].set_xlabel("|position shift| (ms)"); ax[1, 1].set_ylabel("count")
    ax[1, 1].legend(loc="best", fontsize=8)

    fig.suptitle("E-3 Savitzky-Golay derivative bank \u2014 accuracy", fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"accuracy plot written: {path}")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--header", default=os.path.join(here, "ppg_derivative.hpp"))
    ap.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    ap.add_argument("--fs", type=float, default=1000.0)  # codebase runs at 1000 Hz
    ap.add_argument("--snr-db", type=float, default=30.0)
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--hlow", type=int, default=12, help="sg_halfwidth_low (at 500 Hz)")
    ap.add_argument("--hhigh", type=int, default=20, help="sg_halfwidth_high (at 500 Hz)")
    ap.add_argument("--order", type=int, default=4, help="sg_poly_order")
    ap.add_argument("--gate", choices=["rms", "mean", "median"], default="rms",
                    help="criterion-2 aggregate statistic to gate on")
    ap.add_argument("--plot-path", default=os.path.join(here, "e3_accuracy.png"),
                    help="where to write the accuracy plot")
    ap.add_argument("--no-plot", action="store_true", help="skip the accuracy plot")
    args = ap.parse_args()

    print(f"E-3 acceptance test  (header: {args.header})")
    global _RUN_ENV
    exe, _RUN_ENV = build_harness(args.header, args.cxx)
    rng = np.random.default_rng(args.seed)
    cfg = (args.hlow, args.hhigh, args.order)
    print(f"config: sg_halfwidth_low={args.hlow} sg_halfwidth_high={args.hhigh} sg_poly_order={args.order}")

    ok1, cdata = test_curvature(exe, args.fs, cfg)
    ok2, ndata = test_noise_shift(exe, args.fs, args.snr_db, args.trials, rng, cfg, args.gate)

    if not args.no_plot:
        make_accuracy_plot(cdata, ndata, args.snr_db, args.plot_path)

    print("=" * 68)
    allok = ok1 and ok2
    print(f"OVERALL: {'PASS' if allok else 'FAIL'}")
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
