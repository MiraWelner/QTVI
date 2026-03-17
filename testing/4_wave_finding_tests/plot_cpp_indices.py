"""
plot_rpeaks.py

Reads a _wave_data.bin file produced by the C++ pipeline and generates
a .png for each requested bin showing all ECG channels with R-peak markers.

Usage:
    python plot_rpeaks.py path/to/wave_data.bin 33 34 473
    python plot_rpeaks.py path/to/wave_data.bin 33 34 473 --outdir plots
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
    with open(path, "rb") as f:
        num_bins = struct.unpack("<Q", f.read(8))[0]

        for _ in range(num_bins):
            b = {}
            b["ecgFs"] = struct.unpack("<d", f.read(8))[0]
            b["ppgFs"] = struct.unpack("<d", f.read(8))[0]

            def read_idx():
                sz = struct.unpack("<Q", f.read(8))[0]
                vals = []
                for _ in range(sz):
                    v = struct.unpack("<Q", f.read(8))[0]
                    vals.append(v - 1 if v >= 1 else 0)  # 1-based to 0-based
                return vals

            # --- 3 R-index arrays (new format) ---
            b["ecgRIndex"] = read_idx()
            b["ecgRIndex2"] = read_idx()
            b["ecgRIndex3"] = read_idx()

            b["ppgMaxAmps"] = read_idx()
            b["ppgMinAmps"] = read_idx()

            def read_signal():
                sz = struct.unpack("<Q", f.read(8))[0]
                if sz > 0:
                    return list(struct.unpack(f"<{sz}d", f.read(sz * 8)))
                return []

            b["ppgSignal"] = read_signal()

            b["ecgSignal"] = read_signal()
            b["ecgSignal2"] = read_signal()
            b["ecgSignal3"] = read_signal()

            num_pairs = struct.unpack("<Q", f.read(8))[0]
            pairs = []
            for _ in range(num_pairs):
                p0 = struct.unpack("<q", f.read(8))[0]
                p1 = struct.unpack("<q", f.read(8))[0]
                pairs.append((p0, p1))
            b["pairs"] = pairs

            bins.append(b)

    return bins


# Channel definitions: (signal_key, rpeak_key, label, color)
CHANNELS = [
    ("ecgSignal", "ecgRIndex", "ECG Ch1", "blue"),
    ("ecgSignal2", "ecgRIndex2", "ECG Ch2", "red"),
    ("ecgSignal3", "ecgRIndex3", "ECG Ch3", "green"),
]


def plot_bin(b, bin_idx, file_id, out_dir):
    fs = b["ecgFs"]
    t_start, t_end = 0, 60

    # Figure out how many channels actually have data
    active = [
        (sig_k, rp_k, label, color)
        for sig_k, rp_k, label, color in CHANNELS
        if len(b[sig_k]) > 0
    ]

    n_ch = len(active)
    fig, axes = plt.subplots(n_ch, 1, figsize=(14, 4 * n_ch), squeeze=False)

    for row, (sig_key, rp_key, label, color) in enumerate(active):
        ax = axes[row, 0]
        ecg = np.array(b[sig_key])
        r_peaks = np.array(b[rp_key])

        start_idx = max(0, min(int(t_start * fs), len(ecg) - 1))
        end_idx = max(0, min(int(t_end * fs), len(ecg)))

        time_subset = np.arange(start_idx, end_idx) / fs
        ecg_subset = ecg[start_idx:end_idx]

        ax.plot(time_subset, ecg_subset, color="0.4", linewidth=0.5, label=label)

        if len(r_peaks) > 0:
            valid_mask = (r_peaks >= start_idx) & (r_peaks < end_idx)
            valid_r = r_peaks[valid_mask]
            valid_r_times = valid_r / fs

            ax.scatter(
                valid_r_times,
                ecg[valid_r],
                marker="o",
                c=color,
                s=50,
                zorder=5,
                linewidths=0.5,
            )

        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Amplitude")
        ax.set_title(f"{file_id} — Bin {bin_idx} | {label} | {len(r_peaks)} Peaks")
        ax.set_xlim(t_start, t_end)

    fig.tight_layout()
    out_path = os.path.join(out_dir, f"{file_id}_bin{bin_idx:03d}_cpp.png")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"Saved: {out_path}")


def main():
    if len(sys.argv) < 3:
        print(
            "Usage: python plot_rpeaks.py <wave_data.bin> <bin0> [bin1] [bin2] ... [--outdir DIR]"
        )
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

    for b in bin_indices:
        if b < 0 or b >= len(data):
            print(f"Bin {b} out of range (0–{len(data) - 1}), skipping.")
            continue
        plot_bin(data[b], b, file_id, out_dir)
        print(f"Bin {b}: ecgRIndex has {len(data[b]['ecgRIndex'])} values")
        pd.DataFrame({"r_peak_index": data[b]["ecgRIndex"]}).to_csv(
            f"{out_dir}/{file_id}_bin_{b}_indices.csv", index=False
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
