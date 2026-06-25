#!/usr/bin/env python3
"""
accel_histograms.py

Walk a folder of .edf recordings, pull the accelerometer channels out of each,
and write a histogram figure (one subplot per axis) per file.

Usage:
    python accel_histograms.py /path/to/edf_folder [-o /path/to/output] [--bins 100]

Dependencies:
    pip install pyedflib numpy matplotlib scipy
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")  # no display needed; we only write PNGs
import matplotlib.pyplot as plt
import pyedflib
from scipy.stats import kurtosis, skew

# Exact channel labels per axis (these are fixed in the dataset).
ACCEL_LABELS = {
    "X": "Accelerometer_X",
    "Y": "Accelerometer_Y",
    "Z": "Accelerometer_Z",
}
EXPECTED_RATE_HZ = 25.0

# Histogram x-axis range: None = autoscale per file, or a (lo, hi) tuple to fix it.
X_RANGE = None


def read_accel(edf_path):
    """Return {axis: samples} for whichever accel axes are present in the file."""
    out = {}
    with pyedflib.EdfReader(str(edf_path)) as f:
        labels = f.getSignalLabels()
        for axis, label in ACCEL_LABELS.items():
            if label not in labels:
                continue
            idx = labels.index(label)
            rate = f.getSampleFrequency(idx)
            if abs(rate - EXPECTED_RATE_HZ) > 1e-6:
                print(
                    f"  WARNING {edf_path.name} {label}: {rate} Hz "
                    f"(expected {EXPECTED_RATE_HZ})"
                )
            out[axis] = f.readSignal(idx)
    return out


def write_stats_csv(results, out_dir):
    """One row per (file, axis): count, mean, median, min, max, std, kurtosis, skew."""
    out_path = out_dir / "accel_stats.csv"
    fields = [
        "file",
        "axis",
        "num_samples",
        "mean",
        "median",
        "min",
        "max",
        "std",
        "kurtosis",
        "skew",
    ]
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for stem, accel in results:
            for axis in ("X", "Y", "Z"):
                data = accel.get(axis)
                if data is None or not data.size:
                    w.writerow([stem, axis, 0, "", "", "", "", "", "", ""])
                    continue
                w.writerow(
                    [
                        stem,
                        axis,
                        data.size,
                        np.mean(data),
                        np.median(data),
                        np.min(data),
                        np.max(data),
                        np.std(data),  # population std (ddof=0)
                        kurtosis(data),  # Fisher (normal -> 0), bias-corrected
                        skew(data),
                    ]
                )
    print(f"Wrote stats -> {out_path}")
    return out_path


def plot_all(results, out_dir, bins):
    """One figure: a row of X/Y/Z histograms per file, all in a single PNG.

    `results` is a list of (stem, {axis: samples}) for files that had accel data.
    """
    if not results:
        print("No accelerometer data found in any file; nothing to plot.")
        return None

    axis_order = ("X", "Y", "Z")
    n_rows = len(results)
    fig, subplots = plt.subplots(n_rows, 3, figsize=(15, 3 * n_rows), squeeze=False)

    for row, (stem, accel) in enumerate(results):
        for col, axis_name in enumerate(axis_order):
            ax = subplots[row][col]
            data = accel.get(axis_name)
            if data is not None and data.size:
                ax.hist(
                    data,
                    bins=bins,
                    range=X_RANGE,
                    color="steelblue",
                    edgecolor="none",
                    zorder=2,
                )
                if X_RANGE is not None:
                    ax.set_xlim(X_RANGE)
                ax.axvline(0, color="red", linewidth=1, zorder=0)
            else:
                ax.text(
                    0.5,
                    0.5,
                    "no data",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    color="gray",
                )
            if row == 0:
                ax.set_title(f"Accel {axis_name}")
            if col == 0:
                ax.set_ylabel(f"{stem}\ncount", fontsize=8)
            if row == n_rows - 1:
                ax.set_xlabel("value")

    fig.tight_layout()
    out_path = out_dir / "accel_histograms_all.png"
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"\nWrote combined figure -> {out_path}")
    return out_path


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "folder",
        type=Path,
        nargs="?",
        default=Path(r"D:\USERS\MiraWelner\QTVI\data\bittium_raw_files"),
        help="folder containing .edf files "
        r"(default: D:\USERS\MiraWelner\QTVI\data\bittium_raw_files)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="output folder for PNGs (default: <folder>/accel_histograms)",
    )
    parser.add_argument("--bins", type=int, default=100, help="histogram bin count")
    parser.add_argument(
        "-n",
        "--limit",
        type=int,
        default=None,
        help="process only the first N files (default: all)",
    )
    parser.add_argument(
        "--recursive", action="store_true", help="also search subfolders"
    )
    args = parser.parse_args(argv)

    if not args.folder.is_dir():
        parser.error(f"not a folder: {args.folder}")

    out_dir = args.output or (args.folder / "accel_histograms")
    out_dir.mkdir(parents=True, exist_ok=True)

    pattern = "**/*.edf" if args.recursive else "*.edf"
    edf_files = sorted(p for p in args.folder.glob(pattern) if p.is_file())
    if not edf_files:
        print(f"No .edf files found in {args.folder}")
        return 0

    if args.limit is not None:
        edf_files = edf_files[: args.limit]
        print(f"Limiting to first {len(edf_files)} file(s).")

    results = []
    for edf_path in edf_files:
        print(f"=== {edf_path.name} ===")
        try:
            accel = read_accel(edf_path)
        except Exception as e:  # corrupt/unreadable EDF
            print(f"  ERROR reading {edf_path.name}: {e}")
            continue
        if any(v.size for v in accel.values()):
            results.append((edf_path.stem, accel))
        else:
            print(f"  no accelerometer channels found; skipping")

    plot_all(results, out_dir, args.bins)
    write_stats_csv(results, out_dir)
    print(f"\nDone. {len(results)} of {len(edf_files)} files had accel data.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
