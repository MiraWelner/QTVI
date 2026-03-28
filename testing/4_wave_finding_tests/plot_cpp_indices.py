"""
plot_rpeaks.py

Reads a _wave_data.bin file produced by the C++ pipeline and generates
a .png for each requested bin showing all ECG channels with R-peak markers
from 3 detection methods (raw, squared, abs) in a grid layout.

Rows = channels, Columns = methods.
Noisy channels are highlighted with a red background.

Usage:
    python plot_rpeaks.py path/to/wave_data.bin 33 34 473
"""

import os
import struct
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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

            # 9 R-peak index arrays: ch1 raw/sq/abs, ch2 raw/sq/abs, ch3 raw/sq/abs
            b["ch1_raw"] = read_idx()
            b["ch1_squared"] = read_idx()
            b["ch1_absval"] = read_idx()
            b["ch2_raw"] = read_idx()
            b["ch2_squared"] = read_idx()
            b["ch2_absval"] = read_idx()
            b["ch3_raw"] = read_idx()
            b["ch3_squared"] = read_idx()
            b["ch3_absval"] = read_idx()

            # PPG indices
            b["ppgMaxAmps"] = read_idx()
            b["ppgMinAmps"] = read_idx()

            # Signals
            b["ppgSignal"] = read_signal()
            b["ecgSignal"] = read_signal()
            b["ecgSignal2"] = read_signal()
            b["ecgSignal3"] = read_signal()

            # 9 noise flags
            raw = f.read(9)
            if len(raw) == 9:
                b["ch1_raw_noisy"] = raw[0] != 0
                b["ch1_squared_noisy"] = raw[1] != 0
                b["ch1_absval_noisy"] = raw[2] != 0
                b["ch2_raw_noisy"] = raw[3] != 0
                b["ch2_squared_noisy"] = raw[4] != 0
                b["ch2_absval_noisy"] = raw[5] != 0
                b["ch3_raw_noisy"] = raw[6] != 0
                b["ch3_squared_noisy"] = raw[7] != 0
                b["ch3_absval_noisy"] = raw[8] != 0
            else:
                for k in [
                    "ch1_raw_noisy",
                    "ch1_squared_noisy",
                    "ch1_absval_noisy",
                    "ch2_raw_noisy",
                    "ch2_squared_noisy",
                    "ch2_absval_noisy",
                    "ch3_raw_noisy",
                    "ch3_squared_noisy",
                    "ch3_absval_noisy",
                ]:
                    b[k] = False

            # Pairs
            raw = f.read(8)
            if len(raw) == 8:
                num_pairs = struct.unpack("<Q", raw)[0]
            else:
                num_pairs = 0
            pairs = []
            if num_pairs < 1_000_000:
                for _ in range(num_pairs):
                    raw = f.read(16)
                    if len(raw) < 16:
                        break
                    p0, p1 = struct.unpack("<qq", raw)
                    pairs.append((p0, p1))
            b["pairs"] = pairs

            bins.append(b)

    return bins


# (signal_key, channel_label, dot_color,
#  raw_rpeak_key, raw_noise_key,
#  sq_rpeak_key, sq_noise_key,
#  abs_rpeak_key, abs_noise_key)
CHANNELS = [
    (
        "ecgSignal",
        "Ch1",
        "blue",
        "ch1_raw",
        "ch1_raw_noisy",
        "ch1_squared",
        "ch1_squared_noisy",
        "ch1_absval",
        "ch1_absval_noisy",
    ),
    (
        "ecgSignal2",
        "Ch2",
        "red",
        "ch2_raw",
        "ch2_raw_noisy",
        "ch2_squared",
        "ch2_squared_noisy",
        "ch2_absval",
        "ch2_absval_noisy",
    ),
    (
        "ecgSignal3",
        "Ch3",
        "green",
        "ch3_raw",
        "ch3_raw_noisy",
        "ch3_squared",
        "ch3_squared_noisy",
        "ch3_absval",
        "ch3_absval_noisy",
    ),
]

METHODS = ["Raw", "Squared", "Abs"]


def plot_bin(b, bin_idx, file_id, out_dir):
    fs = b["ecgFs"]
    t_start, t_end = 0, 60

    active = [ch for ch in CHANNELS if len(b[ch[0]]) > 0]

    if not active:
        print(f"  Bin {bin_idx}: no active ECG channels, skipping plot.")
        return

    n_rows = len(active)
    n_cols = 3  # raw, squared, abs

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(8 * n_cols, 4 * n_rows), squeeze=False
    )

    for row, ch_def in enumerate(active):
        (
            sig_key,
            ch_label,
            color,
            raw_rp,
            raw_noise,
            sq_rp,
            sq_noise,
            abs_rp,
            abs_noise,
        ) = ch_def

        ecg = np.array(b[sig_key])
        start_idx = max(0, min(int(t_start * fs), len(ecg) - 1))
        end_idx = max(0, min(int(t_end * fs), len(ecg)))
        time_subset = np.arange(start_idx, end_idx) / fs
        ecg_subset = ecg[start_idx:end_idx]

        method_keys = [
            (raw_rp, raw_noise),
            (sq_rp, sq_noise),
            (abs_rp, abs_noise),
        ]

        for col, (rp_key, noise_key) in enumerate(method_keys):
            ax = axes[row, col]
            r_peaks = np.array(b[rp_key])
            is_noisy = b[noise_key]

            if is_noisy:
                ax.set_facecolor("#ffe0e0")

            ax.plot(time_subset, ecg_subset, color="0.4", linewidth=0.5)

            y_min = float(np.percentile(ecg_subset, 0.1))
            y_max = float(np.percentile(ecg_subset, 99.9))

            if len(r_peaks) > 0:
                valid_mask = (r_peaks >= start_idx) & (r_peaks < end_idx)
                valid_r = r_peaks[valid_mask]
                valid_r_times = valid_r / fs

                ax.scatter(
                    valid_r_times,
                    ecg[valid_r],
                    marker="o",
                    c=color,
                    s=30,
                    zorder=5,
                    linewidths=0.5,
                )

                if len(valid_r) > 0:
                    y_min = min(y_min, float(ecg[valid_r].min()))
                    y_max = max(y_max, float(ecg[valid_r].max()))

            padding = (y_max - y_min) * 0.15
            ax.set_ylim(y_min - padding, y_max + padding)
            ax.set_xlim(t_start, t_end)

            noise_tag = " [NOISY]" if is_noisy else ""
            title_color = "red" if is_noisy else "black"
            ax.set_title(
                f"{ch_label} {METHODS[col]} — {len(r_peaks)} Peaks{noise_tag}",
                color=title_color,
                fontsize=10,
            )

            if col == 0:
                ax.set_ylabel(ch_label)
            if row == n_rows - 1:
                ax.set_xlabel("Time (s)")

    fig.tight_layout()
    out_path = os.path.join(out_dir, f"{file_id}_bin{bin_idx:03d}_cpp.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def main():
    if len(sys.argv) < 3:
        print("Usage: python plot_rpeaks.py <wave_data.bin> <bin0> [bin1] [bin2] ...")
        print("Example: python plot_rpeaks.py results/3010112_wave_data.bin 33 34 473")
        sys.exit(1)

    bin_path = sys.argv[1]
    out_dir = "cpp_output"
    args = sys.argv[2:]
    bin_indices = [int(x) for x in args]

    os.makedirs(out_dir, exist_ok=True)
    file_id = Path(bin_path).stem.replace("_wave_data", "")

    print(f"Reading: {bin_path}")
    data = read_wave_bin(bin_path)
    print(f"Loaded {len(data)} bins.\n")

    for b_idx in bin_indices:
        if b_idx < 0 or b_idx >= len(data):
            print(f"Bin {b_idx} out of range (0–{len(data) - 1}), skipping.")
            continue

        b = data[b_idx]
        plot_bin(b, b_idx, file_id, out_dir)

        # Summary
        for ch_name, rp_raw, rp_sq, rp_abs in [
            ("Ch1", "ch1_raw", "ch1_squared", "ch1_absval"),
            ("Ch2", "ch2_raw", "ch2_squared", "ch2_absval"),
            ("Ch3", "ch3_raw", "ch3_squared", "ch3_absval"),
        ]:
            print(
                f"  {ch_name}: raw={len(b[rp_raw])}, sq={len(b[rp_sq])}, abs={len(b[rp_abs])}"
            )

        # CSV export
        fs = b["ecgFs"]
        columns = {}
        for ch_name, rp_raw, rp_sq, rp_abs in [
            ("ch1", "ch1_raw", "ch1_squared", "ch1_absval"),
            ("ch2", "ch2_raw", "ch2_squared", "ch2_absval"),
            ("ch3", "ch3_raw", "ch3_squared", "ch3_absval"),
        ]:
            columns[f"{ch_name}_raw_sec"] = [idx / fs for idx in b[rp_raw]]
            columns[f"{ch_name}_squared_sec"] = [idx / fs for idx in b[rp_sq]]
            columns[f"{ch_name}_absval_sec"] = [idx / fs for idx in b[rp_abs]]

        max_len = max(len(v) for v in columns.values()) or 1
        for k in columns:
            columns[k] += [""] * (max_len - len(columns[k]))

        pd.DataFrame(columns).to_csv(
            f"{out_dir}/{file_id}_bin_{b_idx}_indices.csv",
            index=False,
            float_format="%.6f",
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
