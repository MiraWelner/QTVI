"""
plot_rpeaks.py

Reads a _wave_data.bin file produced by the C++ pipeline and generates
a .png for each requested bin showing the ECG with R-peak markers.

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

            b["ecgRIndex"] = read_idx()
            b["ppgMaxAmps"] = read_idx()
            b["ppgMinAmps"] = read_idx()

            def read_signal():
                sz = struct.unpack("<Q", f.read(8))[0]
                if sz > 0:
                    return list(struct.unpack(f"<{sz}d", f.read(sz * 8)))
                return []

            b["ppgSignal"] = read_signal()
            b["ecgSignal"] = read_signal()

            num_pairs = struct.unpack("<Q", f.read(8))[0]
            pairs = []
            for _ in range(num_pairs):
                p0 = struct.unpack("<q", f.read(8))[0]
                p1 = struct.unpack("<q", f.read(8))[0]
                pairs.append((p0, p1))
            b["pairs"] = pairs

            bins.append(b)

    return bins


def plot_bin(b, bin_idx, file_id, out_dir):
    # 1. Use the actual sampling rate from the file
    fs = b["ecgFs"]

    # Define time range in seconds
    t_start, t_end = 0.0, 60.0

    # 2. Convert time to sample indices
    start_idx = int(t_start * fs)
    end_idx = int(t_end * fs)

    ecg = np.array(b["ecgSignal"])
    r_peaks = np.array(b["ecgRIndex"])

    if len(ecg) == 0:
        print(f"Bin {bin_idx}: empty ECG, skipping.")
        return

    # Ensure indices are within bounds of the signal
    start_idx = max(0, min(start_idx, len(ecg) - 1))
    end_idx = max(0, min(end_idx, len(ecg)))

    # 3. Create time vector in SECONDS for the x-axis
    # This matches the Matlab: time = (idx_start:idx_end-1) / fs
    time_subset = np.arange(start_idx, end_idx) / fs
    ecg_subset = ecg[start_idx:end_idx]

    # 4. Filter R-peaks that fall within this index range
    valid_mask = (r_peaks >= start_idx) & (r_peaks < end_idx)
    valid_r_indices = r_peaks[valid_mask]

    # Convert filtered R-peak indices to seconds for plotting
    valid_r_times = valid_r_indices / fs

    fig, ax = plt.subplots(figsize=(14, 4))  # Adjusted height to match Matlab 1400x400

    # Plot ECG against time in seconds
    ax.plot(time_subset, ecg_subset, color="0.4", linewidth=0.5, label="ECG")

    # Plot R-peaks against time in seconds
    if len(valid_r_times) > 0:
        ax.scatter(
            valid_r_times,
            ecg[valid_r_indices],
            marker="v",
            c="red",
            s=50,
            zorder=5,
            edgecolors="black",  # Match Matlab 'k'
            linewidths=0.5,
            label="R-peaks",
        )

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude")
    ax.set_title(
        f"{file_id} — Bin {bin_idx} | {len(valid_r_indices)} R-peaks | Fs = {fs:.1f} Hz"
    )

    # Set x-limits to the requested seconds
    ax.set_xlim(t_start, t_end)
    # ax.set_ylim(-0.4, 1.4)

    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.15)
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

    # Parse optional --outdir
    out_dir = "rpeaks_plots"
    args = sys.argv[2:]
    if "--outdir" in args:
        idx = args.index("--outdir")
        out_dir = args[idx + 1]
        args = args[:idx] + args[idx + 2 :]

    bin_indices = [int(x) for x in args]

    os.makedirs(out_dir, exist_ok=True)

    # Extract file ID
    file_id = Path(bin_path).stem.replace("_wave_data", "")

    print(f"Reading: {bin_path}")
    data = read_wave_bin(bin_path)
    print(f"Loaded {len(data)} bins.\n")

    for b in bin_indices:
        if b < 0 or b >= len(data):
            print(f"Bin {b} out of range (0–{len(data) - 1}), skipping.")
            continue
        plot_bin(data[b], b, file_id, out_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
