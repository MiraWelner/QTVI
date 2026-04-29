"""
plot_rpeaks.py

Reads a _wave_markings.bin file produced by the C++ pipeline and generates
a .png for each requested bin showing ECG/PPG channels with R-peak markers
from the QTVI CSV zip file.

Usage:
    python plot_rpeaks.py <subject_id> <bin0> [bin1] ...
    python plot_rpeaks.py 3010023 33 34 473
"""

import csv
import os
import re
import struct
import sys
import zipfile
from itertools import zip_longest
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# USER CONFIGURATION
# =============================================================================

BIN_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_mira"
)
CSV_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\1_bin_mat_files\mesa_bin_deep\output"
    r"\qtvi_results_08-Apr-2026_EpSzSec60\RRQT_results\EveryRRQTintoEntropy"
)
OUT_DIR = Path(r"D:\USERS\MiraWelner\QTVI\testing\4_wave_finding_tests\results_csv")

# Time range in seconds to display. Use None for full signal.
T_START = 0
T_END = 60

# Which preprocessing method(s) to show as columns.
# Options: "all", "raw", "squared", "absval"
SHOW_METHOD = "raw"

# Whether to include the PPG row (with min/max amp markers) when PPG exists.
SHOW_PPG = False

# Scatter point size for the signal trace (small dots).
SIGNAL_SCATTER_SIZE = 0.5

MAX_SANE_IDX = 50_000_000
MAX_SANE_SIG = 200_000_000
PASSTHROUGH_NUM_CHANNELS = 41
SR = 1000.0

# =============================================================================
# QTVI CSV reader
# =============================================================================


def _rel_sort_key(name):
    m = re.search(r"Rel(\d+)_Abs(\d+)", name)
    return (int(m.group(1)), int(m.group(2))) if m else (9999, 9999)


def read_qtvi_bin(zip_path, bin_index):
    """
    Read one bin's worth of QTVI peaks from the zip.
    Returns array of sample indices (time_sec * SR, 0-based within bin).
    """
    if not Path(zip_path).exists():
        return np.array([], dtype=np.intp)

    with zipfile.ZipFile(zip_path, "r") as zf:
        inner_names = sorted(
            [n for n in zf.namelist() if n.lower().endswith(".csv")],
            key=_rel_sort_key,
        )
        if bin_index >= len(inner_names):
            return np.array([], dtype=np.intp)

        name = inner_names[bin_index]
        times = []
        with zf.open(name) as fh:
            for raw in fh:
                line = raw.decode("utf-8", errors="replace").strip()
                if not line:
                    continue
                parts = [p.strip() for p in line.split(",")]
                try:
                    times.append(float(parts[0]))
                except (ValueError, IndexError):
                    continue

    return np.round(np.array(times) * SR).astype(np.intp)


# =============================================================================
# Bin file reader
# =============================================================================


def read_wave_bin(path):
    bins = []
    file_size = os.path.getsize(path)
    with open(path, "rb") as f:
        num_bins = struct.unpack("<Q", f.read(8))[0]
        for i in range(num_bins):
            if f.tell() >= file_size:
                break
            b = {}
            b["ecgFs"] = SR
            b["ppgFs"] = SR
            try:
                _read_bin_contents(f, b, i)
            except (ValueError, struct.error) as e:
                print(f"  ERROR at bin {i}: {e}")
                break
            bins.append(b)
    return bins


def _read_bin_contents(f, b, bin_index):
    def read_idx():
        raw = f.read(8)
        if len(raw) < 8:
            raise ValueError("EOF reading idx count")
        sz = struct.unpack("<Q", raw)[0]
        if sz == 0:
            return []
        if sz > MAX_SANE_IDX:
            raise ValueError(f"idx count {sz} too large")
        data = f.read(sz * 8)
        return [v - 1 if v >= 1 else 0 for v in struct.unpack(f"<{sz}Q", data)]

    def read_signal():
        raw = f.read(8)
        if len(raw) < 8:
            raise ValueError("EOF reading signal count")
        sz = struct.unpack("<Q", raw)[0]
        if sz == 0:
            return []
        if sz > MAX_SANE_SIG:
            raise ValueError(f"signal size {sz} too large")
        data = f.read(sz * 8)
        return list(struct.unpack(f"<{sz}d", data))

    def skip_idx():
        sz = struct.unpack("<Q", f.read(8))[0]
        if sz <= MAX_SANE_IDX:
            f.seek(sz * 8, 1)

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

    b["ch1_raw_idx"] = read_idx()
    b["ch1_sq_idx"] = read_idx()
    b["ch1_abs_idx"] = read_idx()
    b["ch2_raw_idx"] = read_idx()
    b["ch2_sq_idx"] = read_idx()
    b["ch2_abs_idx"] = read_idx()
    b["ch3_raw_idx"] = read_idx()
    b["ch3_sq_idx"] = read_idx()
    b["ch3_abs_idx"] = read_idx()

    b["ppgMaxAmps"] = read_idx()
    b["ppgMinAmps"] = read_idx()

    b["ppgSignal"] = read_signal()
    b["ch1_raw_sig"] = read_signal()
    b["ch2_raw_sig"] = read_signal()
    b["ch3_raw_sig"] = read_signal()

    b["ch1_sq_sig"] = read_signal()
    b["ch1_abs_sig"] = read_signal()
    b["ch2_sq_sig"] = read_signal()
    b["ch2_abs_sig"] = read_signal()
    b["ch3_sq_sig"] = read_signal()
    b["ch3_abs_sig"] = read_signal()

    f.read(9)  # noise flags

    num_pairs = struct.unpack("<Q", f.read(8))[0]
    pairs = []
    if 0 < num_pairs <= MAX_SANE_IDX:
        pair_data = f.read(num_pairs * 16)
        for j in range(num_pairs):
            p0, p1 = struct.unpack_from("<qq", pair_data, j * 16)
            pairs.append((p0, p1))
    b["pairs"] = pairs

    skip_pair_vec()
    skip_pair_vec()

    for _ in range(PASSTHROUGH_NUM_CHANNELS):
        skip_signal()
        skip_raw_pairs()


# =============================================================================
# Grid / plot helpers
# =============================================================================

ALL_CELLS = [
    ("ch1_raw_sig", "Ch1", "raw", "red"),
    ("ch1_sq_sig", "Ch1", "squared", "red"),
    ("ch1_abs_sig", "Ch1", "absval", "red"),
    ("ch2_raw_sig", "Ch2", "raw", "red"),
    ("ch2_sq_sig", "Ch2", "squared", "red"),
    ("ch2_abs_sig", "Ch2", "absval", "red"),
    ("ch3_raw_sig", "Ch3", "raw", "green"),
    ("ch3_sq_sig", "Ch3", "squared", "green"),
    ("ch3_abs_sig", "Ch3", "absval", "green"),
]
METHOD_LABELS = {"raw": "Raw", "squared": "Squared", "absval": "Abs-Value"}


def build_grid(b):
    methods = (
        ["raw", "squared", "absval"]
        if SHOW_METHOD == "all"
        else [SHOW_METHOD]
        if SHOW_METHOD in METHOD_LABELS
        else ["raw", "squared", "absval"]
    )
    ecg_channels = [
        ch for ch in ["Ch1", "Ch2", "Ch3"] if len(b.get(f"ch{ch[-1]}_raw_sig", [])) > 0
    ]
    rows = [(ch, "ecg") for ch in ecg_channels]
    has_ppg = SHOW_PPG and len(b.get("ppgSignal", [])) > 0
    if has_ppg:
        rows.append(("PPG", "ppg"))
    ecg_cells = [c for c in ALL_CELLS if c[1] in ecg_channels and c[2] in methods]
    return rows, methods, ecg_cells, has_ppg


def plot_ecg_cell(ax, b, sig_key, qtvi_peaks, color, t_start, t_end, fs):
    ecg = np.array(b[sig_key])
    if len(ecg) == 0:
        ax.text(
            0.5,
            0.5,
            "No signal",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=12,
            color="gray",
        )
        ax.set_xticks([])
        ax.set_yticks([])
        return 0

    start_i = max(0, min(int(t_start * fs), len(ecg) - 1))
    end_i = max(0, min(int(t_end * fs), len(ecg)))
    t_arr = np.arange(start_i, end_i) / fs
    ecg_sub = ecg[start_i:end_i]

    ax.scatter(t_arr, ecg_sub, s=SIGNAL_SCATTER_SIZE, c="0.2", edgecolors="none")
    y_lo = float(np.percentile(ecg_sub, 0.1))
    y_hi = float(np.percentile(ecg_sub, 99.9))

    n_visible = 0
    if len(qtvi_peaks) > 0:
        mask = (qtvi_peaks >= start_i) & (qtvi_peaks < end_i) & (qtvi_peaks < len(ecg))
        vr = qtvi_peaks[mask]
        ax.scatter(
            vr / fs, ecg[vr], marker="o", c=color, s=35, zorder=5, linewidths=0.5
        )
        if len(vr) > 0:
            y_lo = min(y_lo, float(ecg[vr].min()))
            y_hi = max(y_hi, float(ecg[vr].max()))
        n_visible = len(vr)

    pad = (y_hi - y_lo) * 0.15 if y_hi > y_lo else 0.1
    ax.set_ylim(y_lo - pad, y_hi + pad)
    ax.set_xlim(t_start, t_end)
    return n_visible


def plot_ppg_row(axes_row, b, methods, t_start, t_end):
    fs = b["ppgFs"]
    ppg = np.array(b["ppgSignal"])
    mins = np.array(b["ppgMinAmps"])
    maxs = np.array(b["ppgMaxAmps"])
    start_i = max(0, min(int(t_start * fs), len(ppg) - 1))
    end_i = max(0, min(int(t_end * fs), len(ppg)))
    t_arr = np.arange(start_i, end_i) / fs
    ppg_sub = ppg[start_i:end_i]

    for col, ax in enumerate(axes_row):
        ax.scatter(t_arr, ppg_sub, s=SIGNAL_SCATTER_SIZE, c="0.4", edgecolors="none")
        y_lo = float(np.percentile(ppg_sub, 0.1))
        y_hi = float(np.percentile(ppg_sub, 99.9))

        for arr, marker, c, label in [
            (mins, "v", "purple", "Valleys"),
            (maxs, "^", "orange", "Peaks"),
        ]:
            if len(arr) > 0:
                mask = (arr >= start_i) & (arr < end_i) & (arr < len(ppg))
                v = arr[mask]
                ax.scatter(
                    v / fs,
                    ppg[v],
                    marker=marker,
                    c=c,
                    s=25,
                    zorder=5,
                    linewidths=0.5,
                    label=f"{label} ({len(arr)})",
                )
                if label == "Valleys" and len(v) > 0:
                    y_lo = min(y_lo, float(ppg[v].min()))
                if label == "Peaks" and len(v) > 0:
                    y_hi = max(y_hi, float(ppg[v].max()))

        pad = (y_hi - y_lo) * 0.15 if y_hi > y_lo else 0.1
        ax.set_ylim(y_lo - pad, y_hi + pad)
        ax.set_xlim(t_start, t_end)
        ax.set_title(f"PPG  (valleys:{len(mins)}  peaks:{len(maxs)})", fontsize=11)
        if col == 0:
            ax.set_ylabel("PPG", fontsize=12, fontweight="bold")
            ax.legend(loc="upper right", fontsize=8)


def export_indices_csv(b, qtvi_peaks, bin_idx, file_id, out_dir):
    idx_columns = [
        ("ch1_raw_idx", "ch1_raw_idx"),
        ("ch1_sq_idx", "ch1_sq_idx"),
        ("ch1_abs_idx", "ch1_abs_idx"),
        ("ch2_raw_idx", "ch2_raw_idx"),
        ("ch2_sq_idx", "ch2_sq_idx"),
        ("ch2_abs_idx", "ch2_abs_idx"),
        ("ch3_raw_idx", "ch3_raw_idx"),
        ("ch3_sq_idx", "ch3_sq_idx"),
        ("ch3_abs_idx", "ch3_abs_idx"),
        ("ppgMaxAmps", "ppg_max_amps"),
        ("ppgMinAmps", "ppg_min_amps"),
    ]
    ecg_fs = b["ecgFs"]
    ppg_fs = b["ppgFs"]
    headers = [label + "_sec" for _, label in idx_columns] + ["qtvi_sec"]
    columns = []
    for key, label in idx_columns:
        fs = ppg_fs if key.startswith("ppg") else ecg_fs
        columns.append([idx / fs for idx in b.get(key, [])])
    columns.append([idx / ecg_fs for idx in qtvi_peaks])

    csv_path = os.path.join(out_dir, f"{file_id}_bin{bin_idx:03d}_indices.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for row in zip_longest(*columns, fillvalue=""):
            writer.writerow(row)
    print(f"  Saved CSV: {csv_path}")


def plot_bin(b, qtvi_peaks, bin_idx, file_id, out_dir):
    fs = b["ecgFs"]
    max_len = max(
        len(b.get("ch1_raw_sig", [])),
        len(b.get("ch2_raw_sig", [])),
        len(b.get("ch3_raw_sig", [])),
        len(b.get("ppgSignal", [])),
        1,
    )
    t_start = T_START if T_START is not None else 0
    t_end = T_END if T_END is not None else max_len / fs

    rows, methods, ecg_cells, has_ppg = build_grid(b)
    if not rows:
        print(f"  Bin {bin_idx}: no channels to plot, skipping.")
        return

    n_rows = len(rows)
    n_cols = len(methods)
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(7 * n_cols, 4 * n_rows), squeeze=False
    )

    ecg_row_map = {ch: r for r, (ch, rt) in enumerate(rows) if rt == "ecg"}
    method_col = {m: c for c, m in enumerate(methods)}

    for sig_key, ch, method, color in ecg_cells:
        row = ecg_row_map[ch]
        col = method_col[method]
        ax = axes[row][col]
        n_peaks = plot_ecg_cell(ax, b, sig_key, qtvi_peaks, color, t_start, t_end, fs)
        ax.set_title(
            f"{ch} - {METHOD_LABELS[method]}  (QTVI: {n_peaks} peaks)", fontsize=11
        )
        if col == 0:
            ax.set_ylabel(ch, fontsize=12, fontweight="bold")
        if row == n_rows - 1:
            ax.set_xlabel("Time (s)")

    if has_ppg:
        plot_ppg_row(axes[n_rows - 1], b, methods, t_start, t_end)
        for col in range(n_cols):
            axes[n_rows - 1][col].set_xlabel("Time (s)")

    if n_cols > 1:
        for col, m in enumerate(methods):
            axes[0][col].annotate(
                METHOD_LABELS[m],
                xy=(0.5, 1.15),
                xycoords="axes fraction",
                ha="center",
                fontsize=13,
                fontweight="bold",
            )

    time_label = f"{t_start}-{t_end}s" if T_END is not None else "full"
    fig.suptitle(
        f"{file_id}  -  Bin {bin_idx}  [{time_label}]  [QTVI peaks]",
        fontsize=15,
        fontweight="bold",
        y=1.0,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    method_tag = SHOW_METHOD if SHOW_METHOD != "all" else "3x3"
    out_path = os.path.join(
        out_dir, f"{file_id}_bin{bin_idx:03d}_{method_tag}_qtvi.png"
    )
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# =============================================================================
# Main
# =============================================================================


def main():
    if len(sys.argv) < 3:
        print("Usage: python plot_rpeaks.py <subject_id> <bin0> [bin1] ...")
        print("Example: python plot_rpeaks.py 3010023 33 34 473")
        sys.exit(1)

    subject_id = sys.argv[1]
    bin_indices = [int(x) for x in sys.argv[2:]]

    matches = list(BIN_DIR.glob(f"{subject_id}_*_wave_markings.bin"))
    if not matches:
        matches = list(BIN_DIR.glob(f"{subject_id}_wave_markings.bin"))
    if not matches:
        print(f"No wave_markings.bin found for ID {subject_id} in {BIN_DIR}")
        sys.exit(1)
    bin_path = matches[0]

    zip_matches = list(CSV_DIR.glob(f"{subject_id}_*everyRRQTinputIntoEntropy_csv.zip"))
    if not zip_matches:
        print(f"WARNING: No QTVI zip found for {subject_id} in {CSV_DIR}")
        zip_path = None
    else:
        zip_path = zip_matches[0]
        print(f"QTVI zip: {zip_path.name}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    file_id = bin_path.stem.replace("_wave_markings", "")

    print(f"Reading bin: {bin_path}")
    data = read_wave_bin(str(bin_path))
    print(f"Loaded {len(data)} bins.\n")

    for b_idx in bin_indices:
        if b_idx < 0 or b_idx >= len(data):
            print(f"  Bin {b_idx} out of range (0-{len(data) - 1}), skipping.")
            continue

        b = data[b_idx]
        qtvi_peaks = (
            read_qtvi_bin(zip_path, b_idx) if zip_path else np.array([], dtype=np.intp)
        )

        print(f"Bin {b_idx}:")
        print(f"  QTVI peaks: {len(qtvi_peaks)}")
        print(f"  Ch1 raw signal length: {len(b.get('ch1_raw_sig', []))}")

        plot_bin(b, qtvi_peaks, b_idx, file_id, OUT_DIR)
        export_indices_csv(b, qtvi_peaks, b_idx, file_id, OUT_DIR)

    print("\nDone.")


if __name__ == "__main__":
    main()
