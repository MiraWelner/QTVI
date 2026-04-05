"""
plot_rpeaks.py

Reads a _wave_markings.bin file produced by the C++ pipeline and generates
a .png for each requested bin showing ECG/PPG channels with R-peak markers.

Grid adapts to available channels and selected preprocessing method(s).

Usage:
    python plot_rpeaks.py path/to/wave_markings.bin 33 34 473
"""

import os
import struct
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Time range in seconds to display.  Use None for full signal.
#   Examples:  T_START, T_END = 0, 60       (first minute)
#              T_START, T_END = 30, 90       (30s–90s)
#              T_START, T_END = None, None   (entire signal)
T_START = 0
T_END = 60

# Which preprocessing method(s) to show as columns.
# Options: "all", "raw", "squared", "absval"
#   "all"      -> 3 columns (Raw | Squared | Abs-Value)
#   "raw"      -> 1 column  (Raw only)
#   "squared"  -> 1 column  (Squared only)
#   "absval"   -> 1 column  (Abs-Value only)
SHOW_METHOD = "all"

# Whether to include the PPG row (with min/max amp markers) when PPG exists.
SHOW_PPG = False

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
            b["ecgFs"] = 2000.0
            b["ppgFs"] = 2000.0

            def read_idx():
                raw = f.read(8)
                if len(raw) < 8:
                    return []
                sz = struct.unpack("<Q", raw)[0]
                if sz > 10_000_000:
                    print(f"  WARNING: Suspicious idx count {sz}, skipping.")
                    return []
                vals = []
                for _ in range(sz):
                    raw = f.read(8)
                    if len(raw) < 8:
                        break
                    v = struct.unpack("<Q", raw)[0]
                    vals.append(v - 1 if v >= 1 else 0)
                return vals

            def read_signal():
                raw = f.read(8)
                if len(raw) < 8:
                    return []
                sz = struct.unpack("<Q", raw)[0]
                if sz == 0:
                    return []
                if sz > 100_000_000:
                    print(f"  WARNING: Suspicious signal size {sz}, skipping.")
                    return []
                raw = f.read(sz * 8)
                if len(raw) < sz * 8:
                    return list(
                        struct.unpack(f"<{len(raw) // 8}d", raw[: len(raw) // 8 * 8])
                    )
                return list(struct.unpack(f"<{sz}d", raw))

            def skip_signal():
                raw = f.read(8)
                if len(raw) < 8:
                    return
                sz = struct.unpack("<Q", raw)[0]
                if sz > 100_000_000:
                    return
                f.seek(sz * 8, 1)

            def skip_pair_vec():
                raw = f.read(8)
                if len(raw) < 8:
                    return
                sz = struct.unpack("<Q", raw)[0]
                if sz > 10_000_000:
                    return
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

            # -- 6 preprocessed signals -----------------------------------
            b["ch1_sq_sig"] = read_signal()
            b["ch1_abs_sig"] = read_signal()
            b["ch2_sq_sig"] = read_signal()
            b["ch2_abs_sig"] = read_signal()
            b["ch3_sq_sig"] = read_signal()
            b["ch3_abs_sig"] = read_signal()

            # -- 9 noise flag bytes ---------------------------------------
            f.read(9)

            # -- pairs ----------------------------------------------------
            raw = f.read(8)
            num_pairs = struct.unpack("<Q", raw)[0] if len(raw) == 8 else 0
            pairs = []
            if num_pairs < 1_000_000:
                for _ in range(num_pairs):
                    raw = f.read(16)
                    if len(raw) < 16:
                        break
                    p0, p1 = struct.unpack("<qq", raw)
                    pairs.append((p0, p1))
            b["pairs"] = pairs

            skip_pair_vec()
            skip_pair_vec()

            bins.append(b)

    return bins


# Full 3x3 cell definitions:  (signal_key, idx_key, channel, method, color)
ALL_CELLS = [
    # Ch1
    ("ch1_raw_sig", "ch1_raw_idx", "Ch1", "raw", "blue"),
    ("ch1_sq_sig", "ch1_sq_idx", "Ch1", "squared", "blue"),
    ("ch1_abs_sig", "ch1_abs_idx", "Ch1", "absval", "blue"),
    # Ch2
    ("ch2_raw_sig", "ch2_raw_idx", "Ch2", "raw", "red"),
    ("ch2_sq_sig", "ch2_sq_idx", "Ch2", "squared", "red"),
    ("ch2_abs_sig", "ch2_abs_idx", "Ch2", "absval", "red"),
    # Ch3
    ("ch3_raw_sig", "ch3_raw_idx", "Ch3", "raw", "green"),
    ("ch3_sq_sig", "ch3_sq_idx", "Ch3", "squared", "green"),
    ("ch3_abs_sig", "ch3_abs_idx", "Ch3", "absval", "green"),
]

METHOD_LABELS = {"raw": "Raw", "squared": "Squared", "absval": "Abs-Value"}


def build_grid(b):
    """Return (rows, col_methods) based on available data and config."""

    # Determine which methods (columns) to show
    if SHOW_METHOD == "all":
        methods = ["raw", "squared", "absval"]
    elif SHOW_METHOD in METHOD_LABELS:
        methods = [SHOW_METHOD]
    else:
        print(f"  WARNING: Unknown SHOW_METHOD '{SHOW_METHOD}', defaulting to 'all'.")
        methods = ["raw", "squared", "absval"]

    # Determine which ECG channels are present (have a raw signal)
    ecg_channels = []
    for ch in ["Ch1", "Ch2", "Ch3"]:
        sig_key = f"ch{ch[-1]}_raw_sig"
        if len(b.get(sig_key, [])) > 0:
            ecg_channels.append(ch)

    # Build row list: each ECG channel, then optionally PPG
    rows = []  # list of (row_label, row_type)
    for ch in ecg_channels:
        rows.append((ch, "ecg"))

    has_ppg = SHOW_PPG and len(b.get("ppgSignal", [])) > 0
    if has_ppg:
        rows.append(("PPG", "ppg"))

    # Build the cell list for ECG rows
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

    ax.plot(t_arr, ecg_sub, color="0.4", linewidth=0.5)

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
    """Plot the PPG signal across all method columns (same signal, same markers)."""
    fs = b["ppgFs"]
    ppg = np.array(b["ppgSignal"])
    mins = np.array(b["ppgMinAmps"])
    maxs = np.array(b["ppgMaxAmps"])

    start_i = max(0, min(int(t_start * fs), len(ppg) - 1))
    end_i = max(0, min(int(t_end * fs), len(ppg)))
    t_arr = np.arange(start_i, end_i) / fs
    ppg_sub = ppg[start_i:end_i]

    for col, ax in enumerate(axes_row):
        ax.plot(t_arr, ppg_sub, color="0.4", linewidth=0.5)

        y_lo = float(np.percentile(ppg_sub, 0.1))
        y_hi = float(np.percentile(ppg_sub, 99.9))

        # Valleys (minAmps) in purple
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

        # Peaks (maxAmps) in orange
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


def plot_bin(b, bin_idx, file_id, out_dir):
    fs = b["ecgFs"]

    # Resolve time range
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

    # --- ECG rows ---
    ecg_row_map = {ch: r for r, (ch, rtype) in enumerate(rows) if rtype == "ecg"}
    method_col = {m: c for c, m in enumerate(methods)}

    for sig_key, idx_key, ch, method, color in ecg_cells:
        row = ecg_row_map[ch]
        col = method_col[method]
        ax = axes[row][col]

        n_peaks = plot_ecg_cell(ax, b, sig_key, idx_key, color, t_start, t_end, fs)
        ax.set_title(f"{ch} — {METHOD_LABELS[method]}  ({n_peaks} peaks)", fontsize=11)

        if col == 0:
            ax.set_ylabel(ch, fontsize=12, fontweight="bold")
        if row == n_rows - 1:
            ax.set_xlabel("Time (s)")

    # --- PPG row (last row, spans all columns) ---
    if has_ppg:
        ppg_row = n_rows - 1
        plot_ppg_row(axes[ppg_row], b, methods, t_start, t_end)
        for col in range(n_cols):
            axes[ppg_row][col].set_xlabel("Time (s)")

    # Column headers
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

    time_label = f"{t_start}–{t_end}s" if T_END is not None else "full"
    fig.suptitle(
        f"{file_id}  —  Bin {bin_idx}  [{time_label}]",
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
        print("Usage: python plot_rpeaks.py <wave_markings.bin> <bin0> [bin1] ...")
        print(
            "Example: python plot_rpeaks.py results/3010112_wave_markings.bin 33 34 473"
        )
        sys.exit(1)

    bin_path = sys.argv[1]
    out_dir = "cpp_output"
    bin_indices = [int(x) for x in sys.argv[2:]]

    os.makedirs(out_dir, exist_ok=True)
    file_id = Path(bin_path).stem.replace("_wave_markings", "")

    print(f"Reading: {bin_path}")
    print(
        f"Config:  SHOW_METHOD={SHOW_METHOD}  T_START={T_START}  T_END={T_END}  SHOW_PPG={SHOW_PPG}"
    )
    data = read_wave_bin(bin_path)
    print(f"Loaded {len(data)} bins.\n")

    for b_idx in bin_indices:
        if b_idx < 0 or b_idx >= len(data):
            print(f"  Bin {b_idx} out of range (0–{len(data) - 1}), skipping.")
            continue

        b = data[b_idx]
        print(f"Bin {b_idx}:")
        print(
            f"  Ch1 peaks — raw:{len(b['ch1_raw_idx'])}  sq:{len(b['ch1_sq_idx'])}  abs:{len(b['ch1_abs_idx'])}"
        )
        print(
            f"  Ch2 peaks — raw:{len(b['ch2_raw_idx'])}  sq:{len(b['ch2_sq_idx'])}  abs:{len(b['ch2_abs_idx'])}"
        )
        print(
            f"  Ch3 peaks — raw:{len(b['ch3_raw_idx'])}  sq:{len(b['ch3_sq_idx'])}  abs:{len(b['ch3_abs_idx'])}"
        )
        print(
            f"  PPG max:{len(b['ppgMaxAmps'])}  min:{len(b['ppgMinAmps'])}  pairs:{len(b['pairs'])}"
        )

        plot_bin(b, b_idx, file_id, out_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
