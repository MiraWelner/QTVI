"""
compare_single.py

Compare R-peak detection across all three sources for a single bin of one
subject:
    - Mira (C++ _wave_markings.bin)
    - Daniel (MATLAB _wave_data.mat)
    - Deep (QTVI per-subject folder of Rel*_Abs*.csv files)

Outputs (all in OUT_DIR):
    <id>_bin<N>_mira.png        ECG (Mira's ch1 raw) with Mira's R-peaks
    <id>_bin<N>_daniel.png      Daniel's ECG segment with Daniel's R-peaks
    <id>_bin<N>_deep.png        Mira's ch1 raw signal with Deep's QTVI peaks
    <id>_bin<N>_indices.csv     R-peak times (sec) per source, side by side
    <id>_bin<N>_ssd.csv         Pairwise sum-of-squared R-location differences

Usage:
    python compare_single.py <subject_id> <bin_number>
    e.g. python compare_single.py 3010023 33
"""

import csv
import os
import re
import struct
import sys
from itertools import zip_longest
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

# ============================================================================
# Config
# ============================================================================

BIN_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_mira"
)
MAT_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_daniel"
)
CSV_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_deep"
)
OUT_DIR = Path(r"D:\USERS\MiraWelner\QTVI\testing\4_wave_finding_tests\single_bin")

# Time range in seconds to display in plots. Use None for full bin.
T_START = 0
T_END = 60

SR = 1000.0
SIGNAL_SCATTER_SIZE = 0.5

MAX_SANE_IDX = 50_000_000
MAX_SANE_SIG = 200_000_000
PASSTHROUGH_NUM_CHANNELS = 41


# ============================================================================
# Mira reader: full bin file (we need the full read to skip past pass-through)
# ============================================================================


def read_wave_bin(path):
    """Read all bins from Mira's wave_markings.bin and return per-bin dicts.

    Each bin keeps the ch1 raw R-peak indices and the ch1 raw ECG signal,
    which is what the plotting code needs.
    """
    bins = []
    file_size = os.path.getsize(path)
    with open(path, "rb") as f:
        num_bins = struct.unpack("<Q", f.read(8))[0]

        for i in range(num_bins):
            if f.tell() >= file_size:
                break
            b = {"ecgFs": SR, "ppgFs": SR}
            try:
                _read_bin_contents(f, b)
            except (ValueError, struct.error) as e:
                print(f"  ERROR at bin {i}: {e}")
                break
            bins.append(b)
    return bins


def _read_bin_contents(f, b):
    """Read all fields for one bin from Mira's binary format."""

    def read_idx():
        sz = struct.unpack("<Q", f.read(8))[0]
        if sz == 0:
            return []
        if sz > MAX_SANE_IDX:
            raise ValueError(f"idx count {sz} exceeds limit")
        data = f.read(sz * 8)
        return list(struct.unpack(f"<{sz}Q", data))

    def read_signal():
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

    # 9 R-peak idx arrays (3 methods x 3 channels)
    b["ch1_raw_idx"] = read_idx()
    b["ch1_sq_idx"] = read_idx()
    b["ch1_abs_idx"] = read_idx()
    b["ch2_raw_idx"] = read_idx()
    b["ch2_sq_idx"] = read_idx()
    b["ch2_abs_idx"] = read_idx()
    b["ch3_raw_idx"] = read_idx()
    b["ch3_sq_idx"] = read_idx()
    b["ch3_abs_idx"] = read_idx()

    # PPG indices
    b["ppgMaxAmps"] = read_idx()
    b["ppgMinAmps"] = read_idx()

    # 4 raw signals
    b["ppgSignal"] = read_signal()
    b["ch1_raw_sig"] = read_signal()
    b["ch2_raw_sig"] = read_signal()
    b["ch3_raw_sig"] = read_signal()

    # 6 preprocessed signals
    for _ in range(6):
        skip_signal()

    # 9 noise flag bytes
    f.read(9)

    # pairs
    num_pairs = struct.unpack("<Q", f.read(8))[0]
    if 0 < num_pairs <= MAX_SANE_IDX:
        f.seek(num_pairs * 16, 1)

    # bin index pair vectors
    skip_pair_vec()
    skip_pair_vec()

    # pass-through region
    for _ in range(PASSTHROUGH_NUM_CHANNELS):
        skip_signal()
        skip_raw_pairs()


# ============================================================================
# Daniel reader: MATLAB .mat
# ============================================================================


def read_daniel_bin(mat_path, bin_idx):
    """Return (r_peak_indices_int, r_peak_indices_float, ecg_signal) for one
    bin of Daniel's wave_data.mat. r_peak indices are 0-based local samples.
    """
    mat = sio.loadmat(str(mat_path), squeeze_me=False)
    wave_data = mat["wave_data"]
    cells = list(wave_data.flat)
    if bin_idx >= len(cells):
        return (
            np.array([], dtype=np.intp),
            np.array([], dtype=np.float64),
            np.array([], dtype=np.float64),
        )

    cell = cells[bin_idx]
    if cell is None or (hasattr(cell, "size") and cell.size == 0):
        return (
            np.array([], dtype=np.intp),
            np.array([], dtype=np.float64),
            np.array([], dtype=np.float64),
        )

    obj = cell
    while hasattr(obj, "shape") and obj.shape == (1, 1):
        obj = obj[0, 0]
    if not (hasattr(obj, "dtype") and obj.dtype.names):
        return (
            np.array([], dtype=np.intp),
            np.array([], dtype=np.float64),
            np.array([], dtype=np.float64),
        )

    def get_field(name):
        if name not in obj.dtype.names:
            return np.array([], dtype=np.float64)
        val = obj[name]
        while hasattr(val, "shape") and val.ndim > 1 and 1 in val.shape:
            val = val.squeeze()
        return np.array(val, dtype=np.float64).flatten()

    # Daniel's MATLAB pipeline (FindWaveBounds.m) already converts ecgRIndex
    # from 1-based to 0-based before saving, so use the indices as-is.
    r_f = get_field("ecgRIndex")
    r_int = r_f.astype(np.intp)

    if "ecgSeg" in obj.dtype.names:
        ecg = get_field("ecgSeg")
    elif "ecg" in obj.dtype.names:
        ecg = get_field("ecg")
    else:
        ecg = get_field("ecgSignal")

    return r_int, r_f, ecg


# ============================================================================
# Deep reader: per-subject folder of Rel*_Abs*.csv files
# ============================================================================


def _rel_sort_key(name):
    m = re.search(r"Rel(\d+)_Abs(\d+)", name)
    return (int(m.group(1)), int(m.group(2))) if m else (9999, 9999)


def read_deep_bin(folder_path, bin_idx):
    """Return R-peak indices (0-based local samples) for one bin of Deep's
    per-subject folder.

    The folder contains Rel*_Abs*.csv files; sort them by (Rel, Abs) so
    Rel1 = bin 0, Rel2 = bin 1, etc. Each CSV's first column is time_sec
    from the bin's epoch start; local sample = round(time_sec * SR).
    """
    if folder_path is None:
        return np.array([], dtype=np.intp)
    folder = Path(folder_path)
    if not folder.exists() or not folder.is_dir():
        return np.array([], dtype=np.intp)

    csv_paths = sorted(folder.glob("*.csv"), key=lambda p: _rel_sort_key(p.name))
    if bin_idx >= len(csv_paths):
        return np.array([], dtype=np.intp)

    times = []
    with open(csv_paths[bin_idx], "r", newline="") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(",")]
            try:
                times.append(float(parts[0]))
            except (ValueError, IndexError):
                continue  # header / malformed

    return np.round(np.array(times, dtype=np.float64) * SR).astype(np.intp)


# ============================================================================
# Plot helpers
# ============================================================================


def _slice_window(signal, fs, t_start, t_end):
    """Return (start_i, end_i, t_arr, sig_sub) clipped to the signal length."""
    n = len(signal)
    if n == 0:
        return 0, 0, np.array([]), np.array([])
    if t_end is None:
        t_end_eff = n / fs
    else:
        t_end_eff = t_end
    start_i = max(0, min(int(t_start * fs), n - 1))
    end_i = max(0, min(int(t_end_eff * fs), n))
    if end_i <= start_i:
        return start_i, end_i, np.array([]), np.array([])
    t_arr = np.arange(start_i, end_i) / fs
    return start_i, end_i, t_arr, signal[start_i:end_i]


def plot_one_source(ecg, peaks, fs, t_start, t_end, title, marker_color, out_path):
    """Plot one ECG trace with R-peak markers and save as PNG."""
    fig, ax = plt.subplots(figsize=(12, 4.5))

    if len(ecg) == 0:
        ax.text(
            0.5,
            0.5,
            "No ECG signal available",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=14,
            color="gray",
        )
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        start_i, end_i, t_arr, ecg_sub = _slice_window(ecg, fs, t_start, t_end)
        ax.scatter(t_arr, ecg_sub, s=SIGNAL_SCATTER_SIZE, c="0.2", edgecolors="none")

        if len(ecg_sub) > 0:
            y_lo = float(np.percentile(ecg_sub, 0.1))
            y_hi = float(np.percentile(ecg_sub, 99.9))
        else:
            y_lo, y_hi = -1.0, 1.0

        if len(peaks) > 0:
            peaks_arr = np.asarray(peaks, dtype=np.intp)
            mask = (peaks_arr >= start_i) & (peaks_arr < end_i) & (peaks_arr < len(ecg))
            visible = peaks_arr[mask]
            ax.scatter(
                visible / fs,
                ecg[visible],
                marker="o",
                c=marker_color,
                s=40,
                zorder=5,
                edgecolors="black",
                linewidths=0.5,
                label=f"R-peaks ({len(visible)} visible / {len(peaks)} total)",
            )
            if len(visible) > 0:
                y_lo = min(y_lo, float(ecg[visible].min()))
                y_hi = max(y_hi, float(ecg[visible].max()))
            ax.legend(loc="upper right", fontsize=10)

        pad = (y_hi - y_lo) * 0.15 if y_hi > y_lo else 0.1
        ax.set_ylim(y_lo - pad, y_hi + pad)
        ax.set_xlim(t_start, t_end if t_end is not None else end_i / fs)

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("ECG")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# CSV outputs
# ============================================================================


def write_indices_csv(out_path, mira_peaks, daniel_peaks, deep_peaks, fs):
    """Three columns of R-peak times in seconds, one per source."""
    cols = [
        [idx / fs for idx in mira_peaks],
        [idx / fs for idx in daniel_peaks],
        [idx / fs for idx in deep_peaks],
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mira_sec", "daniel_sec", "deep_sec"])
        for row in zip_longest(*cols, fillvalue=""):
            formatted = [f"{v:.6f}" if isinstance(v, float) else v for v in row]
            w.writerow(formatted)


def pairwise_ssd(a, b, fs):
    """Sum of squared differences between two R-peak arrays (in s^2).

    Compares the first min(len(a), len(b)) peaks element-wise. Returns
    (ssd_seconds_squared, n_compared, n_a, n_b).
    """
    a_arr = np.asarray(a, dtype=np.float64)
    b_arr = np.asarray(b, dtype=np.float64)
    n = min(len(a_arr), len(b_arr))
    if n == 0:
        return float("nan"), 0, len(a_arr), len(b_arr)
    diffs = (a_arr[:n] - b_arr[:n]) / fs
    return float(np.sum(diffs**2)), n, len(a_arr), len(b_arr)


def write_ssd_csv(out_path, mira_peaks, daniel_peaks, deep_peaks, fs):
    """Write pairwise SSDs of R-peak locations between the three sources."""
    pairs = [
        ("Mira", "Daniel", mira_peaks, daniel_peaks),
        ("Mira", "Deep", mira_peaks, deep_peaks),
        ("Daniel", "Deep", daniel_peaks, deep_peaks),
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "source_a",
                "source_b",
                "n_a",
                "n_b",
                "n_compared",
                "ssd_s2",
                "rmse_s",
            ]
        )
        for la, lb, a, b in pairs:
            ssd, n_cmp, na, nb = pairwise_ssd(a, b, fs)
            if n_cmp > 0:
                rmse = float(np.sqrt(ssd / n_cmp))
                w.writerow([la, lb, na, nb, n_cmp, f"{ssd:.6f}", f"{rmse:.6f}"])
            else:
                w.writerow([la, lb, na, nb, 0, "", ""])


# ============================================================================
# Main
# ============================================================================


def find_subject_files(subject_id):
    """Locate Mira's bin, Daniel's mat, and Deep's folder for this subject."""
    bin_matches = list(BIN_DIR.glob(f"{subject_id}_wave_markings.bin"))
    mat_matches = list(MAT_DIR.glob(f"{subject_id}_wave_data.mat"))

    deep_folder = None
    if CSV_DIR.exists():
        for suffix in [
            "_part1_ECG_fs1000_everyRRQTinputIntoEntropy_csv",
            "_ECG_fs1000_everyRRQTinputIntoEntropy_csv",
            "_everyRRQTinputIntoEntropy_csv",
        ]:
            candidate = CSV_DIR / f"{subject_id}{suffix}"
            if candidate.is_dir():
                deep_folder = candidate
                break

    return (
        bin_matches[0] if bin_matches else None,
        mat_matches[0] if mat_matches else None,
        deep_folder,
    )


def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_single.py <subject_id> <bin_number>")
        print("  e.g. python compare_single.py 3010023 33")
        sys.exit(1)

    subject_id = sys.argv[1]
    bin_idx = int(sys.argv[2])

    bin_path, mat_path, deep_folder = find_subject_files(subject_id)
    print(f"Subject:     {subject_id}")
    print(f"Bin index:   {bin_idx}")
    print(f"Mira bin:    {bin_path}")
    print(f"Daniel mat:  {mat_path}")
    print(f"Deep folder: {deep_folder}")

    if bin_path is None and mat_path is None:
        print("\nNo data files found for this subject. Exiting.")
        sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # --- Mira ---
    mira_peaks = np.array([], dtype=np.intp)
    mira_signal = np.array([], dtype=np.float64)
    if bin_path is not None:
        mira_data = read_wave_bin(str(bin_path))
        if 0 <= bin_idx < len(mira_data):
            mira_peaks = np.asarray(mira_data[bin_idx]["ch1_raw_idx"], dtype=np.intp)
            mira_signal = np.asarray(
                mira_data[bin_idx]["ch1_raw_sig"], dtype=np.float64
            )
        else:
            print(f"  Bin {bin_idx} out of range for Mira (max {len(mira_data) - 1})")

    # --- Daniel ---
    daniel_peaks = np.array([], dtype=np.intp)
    daniel_signal = np.array([], dtype=np.float64)
    if mat_path is not None:
        daniel_peaks, _daniel_peaks_f, daniel_signal = read_daniel_bin(
            str(mat_path), bin_idx
        )

    # --- Deep ---
    deep_peaks = read_deep_bin(deep_folder, bin_idx)

    # --- Plot 3 PNGs ---
    plot_one_source(
        mira_signal,
        mira_peaks,
        SR,
        T_START,
        T_END,
        f"{subject_id}  -  Bin {bin_idx}  -  Mira (C++ ch1 raw)",
        "tab:blue",
        OUT_DIR / f"{subject_id}_bin{bin_idx:03d}_mira.png",
    )
    plot_one_source(
        daniel_signal,
        daniel_peaks,
        SR,
        T_START,
        T_END,
        f"{subject_id}  -  Bin {bin_idx}  -  Daniel (MATLAB)",
        "tab:red",
        OUT_DIR / f"{subject_id}_bin{bin_idx:03d}_daniel.png",
    )
    # Deep has no own ECG signal — overlay onto Mira's ch1 raw (same source signal).
    plot_one_source(
        mira_signal,
        deep_peaks,
        SR,
        T_START,
        T_END,
        f"{subject_id}  -  Bin {bin_idx}  -  Deep (QTVI peaks on Mira's ch1 raw)",
        "tab:green",
        OUT_DIR / f"{subject_id}_bin{bin_idx:03d}_deep.png",
    )

    # --- CSV outputs ---
    indices_path = OUT_DIR / f"{subject_id}_bin{bin_idx:03d}_indices.csv"
    write_indices_csv(indices_path, mira_peaks, daniel_peaks, deep_peaks, SR)

    ssd_path = OUT_DIR / f"{subject_id}_bin{bin_idx:03d}_ssd.csv"
    write_ssd_csv(ssd_path, mira_peaks, daniel_peaks, deep_peaks, SR)

    # --- Console summary ---
    print()
    print(
        f"R-peak counts: Mira={len(mira_peaks)}  "
        f"Daniel={len(daniel_peaks)}  Deep={len(deep_peaks)}"
    )
    print()
    print("Pairwise SSD (s^2):")
    for la, lb, a, b in [
        ("Mira", "Daniel", mira_peaks, daniel_peaks),
        ("Mira", "Deep", mira_peaks, deep_peaks),
        ("Daniel", "Deep", daniel_peaks, deep_peaks),
    ]:
        ssd, n_cmp, _, _ = pairwise_ssd(a, b, SR)
        if n_cmp > 0:
            rmse = float(np.sqrt(ssd / n_cmp))
            print(f"  {la:>6} vs {lb:<6} n={n_cmp:>4}  SSD={ssd:.6f}  RMSE={rmse:.6f}s")
        else:
            print(f"  {la:>6} vs {lb:<6} (no overlap)")

    print()
    print(f"Outputs in: {OUT_DIR}")


if __name__ == "__main__":
    main()
