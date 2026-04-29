"""
plot_bin_signals.py

Plot MATLAB (Daniel) and C++ (Mira) ECG signals for the same bin overlaid
on the same axes, so any time-shift between them is visible directly.

Outputs three PNGs in OUT_DIR:
    <id>_bin<N>_overlay_full.png    Full bin overlaid
    <id>_bin<N>_overlay_zoom.png    First 5 seconds zoomed
    <id>_bin<N>_diff.png            Sample-by-sample difference

And one CSV:
    <id>_bin<N>_signals.csv         time_s, matlab_ecg, cpp_ecg, diff

Usage:
    python plot_bin_signals.py <subject_id> <bin_number>
    e.g. python plot_bin_signals.py 3010023 10
"""

import csv
import os
import struct
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

# ============================================================================
# Config -- edit these to match your environment
# ============================================================================

BIN_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_mira"
)
MAT_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_daniel"
)
OUT_DIR = Path(r"D:\USERS\MiraWelner\QTVI\testing\4_wave_finding_tests\signal_overlay")

SR = 1000.0
ZOOM_SECONDS = 5.0   # zoom window for the second overlay plot

MAX_SANE_IDX = 50_000_000
MAX_SANE_SIG = 200_000_000
PASSTHROUGH_NUM_CHANNELS = 41


# ============================================================================
# C++ wave_markings.bin reader -- pulls bin N's ch1_raw signal + R-peaks
# ============================================================================

def read_cpp_bin(path, target_bin):
    """Walk the bin file until target_bin, then return its ch1_raw signal and
    R-peak indices. Skips all other content."""
    file_size = os.path.getsize(path)
    with open(path, "rb") as f:
        num_bins = struct.unpack("<Q", f.read(8))[0]
        if target_bin >= num_bins:
            raise ValueError(f"Bin {target_bin} out of range (have {num_bins} bins)")

        for i in range(num_bins):
            if f.tell() >= file_size:
                raise ValueError(f"File ended before bin {target_bin}")
            sig, peaks = _read_one_bin(f, want=(i == target_bin))
            if i == target_bin:
                return np.asarray(sig, dtype=np.float64), np.asarray(peaks, dtype=np.intp)

    raise ValueError(f"Bin {target_bin} not found in file")


def _read_one_bin(f, want):
    """Read one bin's contents from the file. If want=True, return
    (ch1_raw_signal, ch1_raw_peak_indices); else return ([], []) and skip."""

    def read_idx_arr():
        sz = struct.unpack("<Q", f.read(8))[0]
        if sz == 0:
            return []
        if sz > MAX_SANE_IDX:
            raise ValueError(f"idx count {sz} exceeds limit")
        data = f.read(sz * 8)
        return [v - 1 if v >= 1 else 0 for v in struct.unpack(f"<{sz}Q", data)]

    def skip_idx():
        sz = struct.unpack("<Q", f.read(8))[0]
        if sz <= MAX_SANE_IDX:
            f.seek(sz * 8, 1)

    def read_signal_arr():
        sz = struct.unpack("<Q", f.read(8))[0]
        if sz == 0:
            return []
        if sz > MAX_SANE_SIG:
            raise ValueError(f"signal size {sz} exceeds limit")
        data = f.read(sz * 8)
        return list(struct.unpack(f"<{sz}d", data))

    def skip_signal():
        sz = struct.unpack("<Q", f.read(8))[0]
        if sz <= MAX_SANE_SIG:
            f.seek(sz * 8, 1)

    def skip_raw_pairs():
        n = struct.unpack("<Q", f.read(8))[0]
        if n <= MAX_SANE_IDX:
            f.seek(n * 16, 1)

    def skip_pair_vec():
        sz = struct.unpack("<Q", f.read(8))[0]
        if sz <= MAX_SANE_IDX:
            f.seek(sz * 16, 1)

    # 9 R-peak idx arrays (3 methods x 3 channels). We only want ch1_raw_idx.
    if want:
        peaks = read_idx_arr()
    else:
        skip_idx()
        peaks = []
    for _ in range(8):
        skip_idx()

    # 2 PPG idx arrays
    skip_idx()
    skip_idx()

    # 4 raw signals: ppg, ch1_raw, ch2_raw, ch3_raw. We want ch1_raw.
    skip_signal()                     # ppg
    if want:
        sig = read_signal_arr()       # ch1_raw
    else:
        skip_signal()
        sig = []
    skip_signal()                     # ch2_raw
    skip_signal()                     # ch3_raw

    # 6 preprocessed signals
    for _ in range(6):
        skip_signal()

    # 9 noise flag bytes
    f.read(9)

    # raw pairs at end of bin
    num_pairs = struct.unpack("<Q", f.read(8))[0]
    if 0 < num_pairs <= MAX_SANE_IDX:
        f.seek(num_pairs * 16, 1)

    # 2 bin index pair vectors
    skip_pair_vec()
    skip_pair_vec()

    # pass-through region: 41 x (signal + raw pairs)
    for _ in range(PASSTHROUGH_NUM_CHANNELS):
        skip_signal()
        skip_raw_pairs()

    return sig, peaks


# ============================================================================
# MATLAB wave_data.mat reader
# ============================================================================

def read_mat_bin(path, bin_idx):
    """Return (ecg_signal, peak_indices) for one bin of MATLAB wave_data.mat."""
    mat = sio.loadmat(str(path), squeeze_me=False)
    wave_data = mat["wave_data"]
    cells = list(wave_data.flat)
    if bin_idx >= len(cells):
        raise ValueError(f"Bin {bin_idx} out of range (have {len(cells)} bins)")

    cell = cells[bin_idx]
    obj = cell
    while hasattr(obj, "shape") and obj.shape == (1, 1):
        obj = obj[0, 0]

    def get_field(name):
        if name not in obj.dtype.names:
            return np.array([], dtype=np.float64)
        val = obj[name]
        while hasattr(val, "shape") and val.ndim > 1 and 1 in val.shape:
            val = val.squeeze()
        return np.array(val, dtype=np.float64).flatten()

    # MATLAB stores 1-based indices; subtract 1 for 0-based local samples.
    peaks = (get_field("ecgRIndex") - 1.0).astype(np.intp)

    if "ecgSeg" in obj.dtype.names:
        ecg = get_field("ecgSeg")
    elif "ecg" in obj.dtype.names:
        ecg = get_field("ecg")
    else:
        ecg = get_field("ecgSignal")

    return ecg, peaks


# ============================================================================
# Plot + CSV helpers
# ============================================================================

def write_csv(out_path, mat_sig, cpp_sig, sr):
    """Write a CSV with time, matlab_ecg, cpp_ecg, diff for the overlap range."""
    n = min(len(mat_sig), len(cpp_sig))
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_s", "matlab_ecg", "cpp_ecg", "diff_cpp_minus_mat"])
        for i in range(n):
            t = i / sr
            m = float(mat_sig[i])
            c = float(cpp_sig[i])
            w.writerow([f"{t:.6f}", f"{m:.6f}", f"{c:.6f}", f"{c - m:.6f}"])


def plot_overlay(mat_sig, cpp_sig, mat_peaks, cpp_peaks, sr,
                 t_start, t_end, title, out_path):
    """Plot MATLAB and C++ signals on the same axes, with peaks marked."""
    fig, ax = plt.subplots(figsize=(14, 5))

    n_mat = len(mat_sig)
    n_cpp = len(cpp_sig)
    if n_mat == 0 and n_cpp == 0:
        ax.text(0.5, 0.5, "No signal in either source", transform=ax.transAxes,
                ha="center", va="center", fontsize=14, color="gray")
    else:
        if t_end is None:
            t_end_eff = max(n_mat, n_cpp) / sr
        else:
            t_end_eff = t_end

        # MATLAB trace
        if n_mat > 0:
            si = max(0, int(t_start * sr))
            ei = min(n_mat, int(t_end_eff * sr))
            if ei > si:
                t_mat = np.arange(si, ei) / sr
                ax.plot(t_mat, mat_sig[si:ei], color="tab:red",
                        linewidth=0.9, alpha=0.85, label=f"MATLAB ({n_mat} samples)")

        # C++ trace
        if n_cpp > 0:
            si = max(0, int(t_start * sr))
            ei = min(n_cpp, int(t_end_eff * sr))
            if ei > si:
                t_cpp = np.arange(si, ei) / sr
                ax.plot(t_cpp, cpp_sig[si:ei], color="tab:blue",
                        linewidth=0.9, alpha=0.85, label=f"C++ ({n_cpp} samples)")

        # MATLAB R-peaks (red circles)
        if len(mat_peaks) > 0 and n_mat > 0:
            visible = mat_peaks[(mat_peaks >= int(t_start * sr)) &
                                (mat_peaks < min(int(t_end_eff * sr), n_mat))]
            if len(visible) > 0:
                ax.scatter(visible / sr, mat_sig[visible],
                           marker="o", facecolors="none", edgecolors="darkred",
                           s=60, linewidths=1.5, zorder=5,
                           label=f"MATLAB R-peaks ({len(visible)})")

        # C++ R-peaks (blue squares)
        if len(cpp_peaks) > 0 and n_cpp > 0:
            visible = cpp_peaks[(cpp_peaks >= int(t_start * sr)) &
                                (cpp_peaks < min(int(t_end_eff * sr), n_cpp))]
            if len(visible) > 0:
                ax.scatter(visible / sr, cpp_sig[visible],
                           marker="s", facecolors="none", edgecolors="darkblue",
                           s=60, linewidths=1.5, zorder=5,
                           label=f"C++ R-peaks ({len(visible)})")

        ax.set_xlim(t_start, t_end_eff)

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("ECG amplitude")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_diff(mat_sig, cpp_sig, sr, title, out_path):
    """Plot sample-by-sample difference (C++ minus MATLAB) over the full bin."""
    n = min(len(mat_sig), len(cpp_sig))
    fig, ax = plt.subplots(figsize=(14, 4))

    if n == 0:
        ax.text(0.5, 0.5, "Cannot compute diff (one signal empty)",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=14, color="gray")
    else:
        t = np.arange(n) / sr
        diff = cpp_sig[:n] - mat_sig[:n]
        ax.plot(t, diff, color="purple", linewidth=0.6)
        ax.axhline(0, color="black", linewidth=0.5, alpha=0.5)
        ax.set_xlim(0, n / sr)
        rms = float(np.sqrt(np.mean(diff ** 2)))
        peak = float(np.max(np.abs(diff)))
        ax.set_title(f"{title}\nRMS diff = {rms:.6g}   peak |diff| = {peak:.6g}",
                     fontsize=12, fontweight="bold")

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("C++ − MATLAB")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# Main
# ============================================================================

def main():
    if len(sys.argv) < 3:
        print("Usage: python plot_bin_signals.py <subject_id> <bin_number>")
        print("  e.g. python plot_bin_signals.py 3010023 10")
        sys.exit(1)

    subject_id = sys.argv[1]
    bin_idx = int(sys.argv[2])

    bin_path = BIN_DIR / f"{subject_id}_wave_markings.bin"
    mat_path = MAT_DIR / f"{subject_id}_wave_data.mat"

    print(f"Subject:   {subject_id}")
    print(f"Bin index: {bin_idx}")
    print(f"C++ file:  {bin_path}  ({'exists' if bin_path.exists() else 'MISSING'})")
    print(f"MATLAB:    {mat_path}  ({'exists' if mat_path.exists() else 'MISSING'})")

    if not bin_path.exists() or not mat_path.exists():
        print("\nMissing input file(s). Exiting.")
        sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\nReading C++ bin...")
    cpp_sig, cpp_peaks = read_cpp_bin(str(bin_path), bin_idx)
    print(f"  C++ signal:  {len(cpp_sig)} samples ({len(cpp_sig) / SR:.3f} s)")
    print(f"  C++ peaks:   {len(cpp_peaks)}")

    print("Reading MATLAB mat...")
    mat_sig, mat_peaks = read_mat_bin(str(mat_path), bin_idx)
    print(f"  MATLAB signal: {len(mat_sig)} samples ({len(mat_sig) / SR:.3f} s)")
    print(f"  MATLAB peaks:  {len(mat_peaks)}")

    n_overlap = min(len(mat_sig), len(cpp_sig))
    if n_overlap > 0:
        diff = np.asarray(cpp_sig[:n_overlap]) - np.asarray(mat_sig[:n_overlap])
        rms = float(np.sqrt(np.mean(diff ** 2)))
        peak = float(np.max(np.abs(diff)))
        print(f"\nSample-by-sample diff over {n_overlap} samples:")
        print(f"  RMS  = {rms:.6g}")
        print(f"  peak = {peak:.6g}")

    base = f"{subject_id}_bin{bin_idx:03d}"

    print("\nWriting outputs...")
    csv_path = OUT_DIR / f"{base}_signals.csv"
    write_csv(csv_path, mat_sig, cpp_sig, SR)
    print(f"  {csv_path}")

    full_path = OUT_DIR / f"{base}_overlay_full.png"
    plot_overlay(mat_sig, cpp_sig, mat_peaks, cpp_peaks, SR,
                 t_start=0, t_end=None,
                 title=f"{subject_id}  -  Bin {bin_idx}  -  Full overlay",
                 out_path=full_path)
    print(f"  {full_path}")

    zoom_path = OUT_DIR / f"{base}_overlay_zoom.png"
    plot_overlay(mat_sig, cpp_sig, mat_peaks, cpp_peaks, SR,
                 t_start=0, t_end=ZOOM_SECONDS,
                 title=f"{subject_id}  -  Bin {bin_idx}  -  First {ZOOM_SECONDS:.0f}s",
                 out_path=zoom_path)
    print(f"  {zoom_path}")

    diff_path = OUT_DIR / f"{base}_diff.png"
    plot_diff(mat_sig, cpp_sig, SR,
              title=f"{subject_id}  -  Bin {bin_idx}  -  C++ minus MATLAB",
              out_path=diff_path)
    print(f"  {diff_path}")

    print(f"\nAll outputs in: {OUT_DIR}")


if __name__ == "__main__":
    main()