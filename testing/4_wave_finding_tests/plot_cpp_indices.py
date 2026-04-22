"""
plot_rpeaks.py

Reads a _wave_markings.bin file produced by the C++ pipeline and generates
a .png for each requested bin showing ECG/PPG channels with R-peak markers.
Also exports a CSV of indices for each requested bin.

Grid adapts to available channels and selected preprocessing method(s).

Usage:
    python plot_rpeaks.py path/to/wave_markings.bin 33 34 473
"""

import csv
import os
import struct
import sys
from itertools import zip_longest
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Time range in seconds to display.  Use None for full signal.
T_START = 0
T_END = 60

# Which preprocessing method(s) to show as columns.
# Options: "all", "raw", "squared", "absval"
SHOW_METHOD = "raw"

# Whether to include the PPG row (with min/max amp markers) when PPG exists.
SHOW_PPG = False

# Scatter point size for the signal trace (small dots).
SIGNAL_SCATTER_SIZE = 0.5

# Max sane sizes - at 2000 Hz, 8 hours = 57.6M samples
MAX_SANE_IDX = 50_000_000
MAX_SANE_SIG = 200_000_000

# =============================================================================


def read_wave_bin(path):
    bins = []
    file_size = os.path.getsize(path)
    with open(path, "rb") as f:
        num_bins = struct.unpack("<Q", f.read(8))[0]

        for i in range(num_bins):
            if f.tell() >= file_size:
                print(
                    f"  WARNING: Reached end of file at bin {i}/{num_bins}, stopping."
                )
                break

            b = {}
            b["ecgFs"] = 1000.0
            b["ppgFs"] = 1000.0

            try:
                _read_bin_contents(f, b, i)
            except (ValueError, struct.error) as e:
                print(f"  ERROR at bin {i}: {e}")
                print(f"  Successfully read {len(bins)} bins before error.")
                break

            bins.append(b)

    return bins


def _read_bin_contents(f, b, bin_index):
    """Read all fields for one bin. Raises ValueError on corrupt data."""

    def read_idx():
        raw = f.read(8)
        if len(raw) < 8:
            raise ValueError("Unexpected EOF reading idx count")
        sz = struct.unpack("<Q", raw)[0]
        if sz == 0:
            return []
        if sz > MAX_SANE_IDX:
            raise ValueError(f"idx count {sz} exceeds limit at offset {f.tell()}")
        data = f.read(sz * 8)
        if len(data) < sz * 8:
            raise ValueError(f"Truncated idx array, expected {sz * 8} got {len(data)}")
        vals = struct.unpack(f"<{sz}Q", data)
        return [v - 1 if v >= 1 else 0 for v in vals]

    def read_signal():
        raw = f.read(8)
        if len(raw) < 8:
            raise ValueError("Unexpected EOF reading signal count")
        sz = struct.unpack("<Q", raw)[0]
        if sz == 0:
            return []
        if sz > MAX_SANE_SIG:
            raise ValueError(f"signal size {sz} exceeds limit at offset {f.tell()}")
        data = f.read(sz * 8)
        if len(data) < sz * 8:
            raise ValueError(f"Truncated signal, expected {sz * 8} got {len(data)}")
        return list(struct.unpack(f"<{sz}d", data))

    def skip_idx():
        raw = f.read(8)
        if len(raw) < 8:
            raise ValueError("Unexpected EOF in skip_idx")
        sz = struct.unpack("<Q", raw)[0]
        if sz > MAX_SANE_IDX:
            raise ValueError(f"skip_idx: count {sz} exceeds limit")
        f.seek(sz * 8, 1)

    def skip_pair_vec():
        raw = f.read(8)
        if len(raw) < 8:
            raise ValueError("Unexpected EOF in skip_pair_vec")
        sz = struct.unpack("<Q", raw)[0]
        if sz > MAX_SANE_IDX:
            raise ValueError(f"skip_pair_vec: count {sz} exceeds limit")
        f.seek(sz * 16, 1)

    # -- 9 R-peak index arrays (3 methods x 3 channels) ----------
    b["ch1_raw_idx"] = read_idx()
    b["ch1_sq_idx"] = read_idx()
    b["ch1_abs_idx"] = read_idx()
    b["ch2_raw_idx"] = read_idx()
    b["ch2_sq_idx"] = read_idx()
    b["ch2_abs_idx"] = read_idx()
    b["ch3_raw_idx"] = read_idx()
    b["ch3_sq_idx"] = read_idx()
    b["ch3_abs_idx"] = read_idx()

    # -- PPG indices ----------------------------------------------
    b["ppgMaxAmps"] = read_idx()
    b["ppgMinAmps"] = read_idx()

    # -- 4 raw signals --------------------------------------------
    b["ppgSignal"] = read_signal()
    b["ch1_raw_sig"] = read_signal()
    b["ch2_raw_sig"] = read_signal()
    b["ch3_raw_sig"] = read_signal()

    # -- 6 preprocessed signals (squared + absval per channel) ----
    b["ch1_sq_sig"] = read_signal()
    b["ch1_abs_sig"] = read_signal()
    b["ch2_sq_sig"] = read_signal()
    b["ch2_abs_sig"] = read_signal()
    b["ch3_sq_sig"] = read_signal()
    b["ch3_abs_sig"] = read_signal()

    # -- 9 noise flag bytes ---------------------------------------
    flags_raw = f.read(9)
    if len(flags_raw) < 9:
        raise ValueError("Truncated noise flags")

    # -- pairs ----------------------------------------------------
    raw = f.read(8)
    if len(raw) < 8:
        raise ValueError("Unexpected EOF reading pair count")
    num_pairs = struct.unpack("<Q", raw)[0]
    pairs = []
    if num_pairs > MAX_SANE_IDX:
        raise ValueError(f"Pair count {num_pairs} exceeds limit")
    if num_pairs > 0:
        pair_data = f.read(num_pairs * 16)
        if len(pair_data) < num_pairs * 16:
            raise ValueError(
                f"Truncated pairs, expected {num_pairs * 16} got {len(pair_data)}"
            )
        for j in range(num_pairs):
            p0, p1 = struct.unpack_from("<qq", pair_data, j * 16)
            pairs.append((p0, p1))
    b["pairs"] = pairs

    # -- ppg/ecg bin index pairs ----------------------------------
    skip_pair_vec()
    skip_pair_vec()

    return bin


# Full 3x3 cell definitions:  (signal_key, idx_key, channel, method, color)
ALL_CELLS = [
    ("ch1_raw_sig", "ch1_raw_idx", "Ch1", "raw", "blue"),
    ("ch1_sq_sig", "ch1_sq_idx", "Ch1", "squared", "blue"),
    ("ch1_abs_sig", "ch1_abs_idx", "Ch1", "absval", "blue"),
    ("ch2_raw_sig", "ch2_raw_idx", "Ch2", "raw", "red"),
    ("ch2_sq_sig", "ch2_sq_idx", "Ch2", "squared", "red"),
    ("ch2_abs_sig", "ch2_abs_idx", "Ch2", "absval", "red"),
    ("ch3_raw_sig", "ch3_raw_idx", "Ch3", "raw", "green"),
    ("ch3_sq_sig", "ch3_sq_idx", "Ch3", "squared", "green"),
    ("ch3_abs_sig", "ch3_abs_idx", "Ch3", "absval", "green"),
]

METHOD_LABELS = {"raw": "Raw", "squared": "Squared", "absval": "Abs-Value"}


def build_grid(b):
    """Return (rows, col_methods) based on available data and config."""

    if SHOW_METHOD == "all":
        methods = ["raw", "squared", "absval"]
    elif SHOW_METHOD in METHOD_LABELS:
        methods = [SHOW_METHOD]
    else:
        print(f"  WARNING: Unknown SHOW_METHOD '{SHOW_METHOD}', defaulting to 'all'.")
        methods = ["raw", "squared", "absval"]

    ecg_channels = []
    for ch in ["Ch1", "Ch2", "Ch3"]:
        sig_key = f"ch{ch[-1]}_raw_sig"
        if len(b.get(sig_key, [])) > 0:
            ecg_channels.append(ch)

    rows = []
    for ch in ecg_channels:
        rows.append((ch, "ecg"))

    has_ppg = SHOW_PPG and len(b.get("ppgSignal", [])) > 0
    if has_ppg:
        rows.append(("PPG", "ppg"))

    ecg_cells = [c for c in ALL_CELLS if c[2] in ecg_channels and c[3] in methods]

    return rows, methods, ecg_cells, has_ppg


def plot_ecg_cell(ax, b, sig_key, idx_key, color, t_start, t_end, fs):
    """Plot one ECG signal + R-peak overlay."""
    ecg = np.array(b[sig_key])
    r_peaks = np.array(b[idx_key])

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

    if len(r_peaks) > 0:
        mask = (r_peaks >= start_i) & (r_peaks < end_i) & (r_peaks < len(ecg))
        vr = r_peaks[mask]
        ax.scatter(
            vr / fs, ecg[vr], marker="o", c=color, s=35, zorder=5, linewidths=0.5
        )
        if len(vr) > 0:
            y_lo = min(y_lo, float(ecg[vr].min()))
            y_hi = max(y_hi, float(ecg[vr].max()))

    pad = (y_hi - y_lo) * 0.15 if y_hi > y_lo else 0.1
    ax.set_ylim(y_lo - pad, y_hi + pad)
    ax.set_xlim(t_start, t_end)

    return len(r_peaks)


def plot_ppg_row(axes_row, b, methods, t_start, t_end):
    """Plot the PPG signal across all method columns."""
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

        if len(mins) > 0:
            mask = (mins >= start_i) & (mins < end_i) & (mins < len(ppg))
            vm = mins[mask]
            ax.scatter(
                vm / fs,
                ppg[vm],
                marker="v",
                c="purple",
                s=25,
                zorder=5,
                linewidths=0.5,
                label=f"Valleys ({len(mins)})",
            )
            if len(vm) > 0:
                y_lo = min(y_lo, float(ppg[vm].min()))

        if len(maxs) > 0:
            mask = (maxs >= start_i) & (maxs < end_i) & (maxs < len(ppg))
            vx = maxs[mask]
            ax.scatter(
                vx / fs,
                ppg[vx],
                marker="^",
                c="orange",
                s=25,
                zorder=5,
                linewidths=0.5,
                label=f"Peaks ({len(maxs)})",
            )
            if len(vx) > 0:
                y_hi = max(y_hi, float(ppg[vx].max()))

        pad = (y_hi - y_lo) * 0.15 if y_hi > y_lo else 0.1
        ax.set_ylim(y_lo - pad, y_hi + pad)
        ax.set_xlim(t_start, t_end)
        ax.set_title(f"PPG  (valleys:{len(mins)}  peaks:{len(maxs)})", fontsize=11)

        if col == 0:
            ax.set_ylabel("PPG", fontsize=12, fontweight="bold")
            ax.legend(loc="upper right", fontsize=8)


def export_indices_csv(b, bin_idx, file_id, out_dir):
    """Write a CSV with one column per index array for this bin."""
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

    headers = [label + "_sec" for _, label in idx_columns]
    columns = []
    for key, label in idx_columns:
        fs = ppg_fs if key.startswith("ppg") else ecg_fs
        columns.append([idx / fs for idx in b.get(key, [])])

    csv_path = os.path.join(out_dir, f"{file_id}_bin{bin_idx:03d}_indices.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for row in zip_longest(*columns, fillvalue=""):
            writer.writerow(row)

    print(f"  Saved CSV: {csv_path}")


def plot_bin(b, bin_idx, file_id, out_dir):
    fs = b["ecgFs"]

    max_len_samples = max(
        len(b.get("ch1_raw_sig", [])),
        len(b.get("ch2_raw_sig", [])),
        len(b.get("ch3_raw_sig", [])),
        len(b.get("ppgSignal", [])),
        1,
    )
    t_start = T_START if T_START is not None else 0
    t_end = T_END if T_END is not None else max_len_samples / fs

    rows, methods, ecg_cells, has_ppg = build_grid(b)

    if not rows:
        print(f"  Bin {bin_idx}: no channels to plot, skipping.")
        return

    n_rows = len(rows)
    n_cols = len(methods)

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(7 * n_cols, 4 * n_rows), squeeze=False
    )

    ecg_row_map = {ch: r for r, (ch, rtype) in enumerate(rows) if rtype == "ecg"}
    method_col = {m: c for c, m in enumerate(methods)}

    for sig_key, idx_key, ch, method, color in ecg_cells:
        row = ecg_row_map[ch]
        col = method_col[method]
        ax = axes[row][col]

        n_peaks = plot_ecg_cell(ax, b, sig_key, idx_key, color, t_start, t_end, fs)
        ax.set_title(f"{ch} - {METHOD_LABELS[method]}  ({n_peaks} peaks)", fontsize=11)

        if col == 0:
            ax.set_ylabel(ch, fontsize=12, fontweight="bold")
        if row == n_rows - 1:
            ax.set_xlabel("Time (s)")

    if has_ppg:
        ppg_row = n_rows - 1
        plot_ppg_row(axes[ppg_row], b, methods, t_start, t_end)
        for col in range(n_cols):
            axes[ppg_row][col].set_xlabel("Time (s)")

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
        f"{file_id}  -  Bin {bin_idx}  [{time_label}]",
        fontsize=15,
        fontweight="bold",
        y=1.0,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    method_tag = SHOW_METHOD if SHOW_METHOD != "all" else "3x3"
    out_path = os.path.join(out_dir, f"{file_id}_bin{bin_idx:03d}_{method_tag}.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def main():
    if len(sys.argv) < 3:
        print("Usage: python plot_rpeaks.py <id> <bin0> [bin1] ...")
        print("Example: python plot_rpeaks.py 3010023 33 34 473")
        sys.exit(1)

    subject_id = sys.argv[1]
    bin_indices = [int(x) for x in sys.argv[2:]]

    BIN_DIR = Path(
        r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_mira"
    )
    matches = list(BIN_DIR.glob(f"{subject_id}_wave_markings.bin"))

    if not matches:
        print(f"No wave_markings.bin found for ID {subject_id} in {BIN_DIR}")
        sys.exit(1)
    bin_path = matches[0]

    out_dir = r"D:\USERS\MiraWelner\QTVI\testing\4_wave_finding_tests\cpp_output"
    os.makedirs(out_dir, exist_ok=True)
    file_id = bin_path.stem.replace("_wave_markings", "")

    print(f"Reading: {bin_path}")
    print(
        f"Config:  SHOW_METHOD={SHOW_METHOD}  T_START={T_START}  T_END={T_END}  SHOW_PPG={SHOW_PPG}"
    )
    data = read_wave_bin(str(bin_path))
    print(f"Loaded {len(data)} bins.\n")

    for b_idx in bin_indices:
        if b_idx < 0 or b_idx >= len(data):
            print(f"  Bin {b_idx} out of range (0-{len(data) - 1}), skipping.")
            continue

        b = data[b_idx]
        print(f"Bin {b_idx}:")
        print(
            f"  Ch1 peaks - raw:{len(b['ch1_raw_idx'])}  sq:{len(b['ch1_sq_idx'])}  abs:{len(b['ch1_abs_idx'])}"
        )
        print(
            f"  Ch2 peaks - raw:{len(b['ch2_raw_idx'])}  sq:{len(b['ch2_sq_idx'])}  abs:{len(b['ch2_abs_idx'])}"
        )
        print(
            f"  Ch3 peaks - raw:{len(b['ch3_raw_idx'])}  sq:{len(b['ch3_sq_idx'])}  abs:{len(b['ch3_abs_idx'])}"
        )
        print(
            f"  PPG max:{len(b['ppgMaxAmps'])}  min:{len(b['ppgMinAmps'])}  pairs:{len(b['pairs'])}"
        )

        plot_bin(b, b_idx, file_id, out_dir)
        export_indices_csv(b, b_idx, file_id, out_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
