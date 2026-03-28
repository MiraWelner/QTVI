"""
Visualize the output of writeAnnealedBin (annealed segment .bin files).

Binary format:
  Header:
    uint64  — number of segments
    float64 — PPG sample rate
    float64 — ECG sample rate
    float64 — scoring epoch size (seconds)

  Per segment:
    uint64 nPpgPairs, then nPpgPairs × (uint64 start, uint64 end)
    uint64 nEcgPairs, then nEcgPairs × (uint64 start, uint64 end)
    uint64 nPpg,   then float64[nPpg]
    uint64 nEcg1,  then float64[nEcg1]
    uint64 nEcg2,  then float64[nEcg2]
    uint64 nEcg3,  then float64[nEcg3]
    uint64 nSleep, then float64[nSleep]
"""

import struct
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ── Reader ───────────────────────────────────────────────────────────────────


def read_uint64(f):
    return struct.unpack("<Q", f.read(8))[0]


def read_float64(f):
    return struct.unpack("<d", f.read(8))[0]


def read_float64_array(f, n):
    return np.frombuffer(f.read(n * 8), dtype=np.float64).copy()


def read_index_pairs(f):
    n = read_uint64(f)
    pairs = []
    for _ in range(n):
        s = read_uint64(f)
        e = read_uint64(f)
        pairs.append((s, e))
    return pairs


def read_annealed_bin(path: str):
    """Return (header_dict, list_of_segment_dicts)."""
    segments = []
    with open(path, "rb") as f:
        n_segments = read_uint64(f)
        if n_segments == 0:
            return None, []

        header = {
            "n_segments": n_segments,
            "ppg_sr": read_float64(f),
            "ecg_sr": read_float64(f),
            "epoch_sec": read_float64(f),
        }

        for _ in range(n_segments):
            seg = {}
            seg["ppg_idx"] = read_index_pairs(f)
            seg["ecg_idx"] = read_index_pairs(f)

            n_ppg = read_uint64(f)
            seg["ppg"] = read_float64_array(f, n_ppg)

            n_ecg1 = read_uint64(f)
            seg["ecg1"] = read_float64_array(f, n_ecg1)

            n_ecg2 = read_uint64(f)
            seg["ecg2"] = read_float64_array(f, n_ecg2)

            n_ecg3 = read_uint64(f)
            seg["ecg3"] = read_float64_array(f, n_ecg3)

            n_sleep = read_uint64(f)
            seg["sleep"] = read_float64_array(f, n_sleep)

            segments.append(seg)

    return header, segments


# ── Plotting ─────────────────────────────────────────────────────────────────


def time_axis(signal, sr):
    """Create a time axis in seconds."""
    return np.arange(len(signal)) / sr


def plot_segment(header, seg, seg_index, show=True):
    """Plot all available channels for a single segment."""
    ppg_sr = header["ppg_sr"]
    ecg_sr = header["ecg_sr"]

    channels = []
    if len(seg["ppg"]) > 0:
        channels.append(("PPG", seg["ppg"], ppg_sr))
    if len(seg["ecg1"]) > 0:
        channels.append(("ECG Ch1", seg["ecg1"], ecg_sr))
    if len(seg["ecg2"]) > 0:
        channels.append(("ECG Ch2", seg["ecg2"], ecg_sr))
    if len(seg["ecg3"]) > 0:
        channels.append(("ECG Ch3", seg["ecg3"], ecg_sr))

    n_plots = len(channels) + (1 if len(seg["sleep"]) > 0 else 0)
    if n_plots == 0:
        print(f"Segment {seg_index}: no data to plot.")
        return

    fig, axes = plt.subplots(n_plots, 1, figsize=(14, 3 * n_plots), sharex=False)
    if n_plots == 1:
        axes = [axes]

    fig.suptitle(
        f"Segment {seg_index}  "
        f"(PPG: {len(seg['ppg'])} samples, "
        f"ECG: {len(seg['ecg1'])} samples)",
        fontsize=13,
    )

    for i, (label, data, sr) in enumerate(channels):
        t = time_axis(data, sr)
        axes[i].plot(t, data, linewidth=0.4)
        axes[i].set_ylabel(label)
        axes[i].set_xlabel("Time (s)")
        axes[i].grid(True, alpha=0.3)

    # Sleep stages as a step plot
    if len(seg["sleep"]) > 0:
        ax = axes[-1]
        epoch_sec = header["epoch_sec"]
        t_sleep = np.arange(len(seg["sleep"])) * epoch_sec
        ax.step(t_sleep, seg["sleep"], where="post", linewidth=1.2, color="purple")
        ax.set_ylabel("Sleep Stage")
        ax.set_xlabel("Time (s)")
        ax.invert_yaxis()  # convention: deeper sleep = higher number
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    if show:
        plt.show()
    return fig


def plot_overview(header, segments):
    """Show a timeline of all segments with their index ranges."""
    fig, ax = plt.subplots(figsize=(14, max(3, len(segments) * 0.5)))

    for i, seg in enumerate(segments):
        # Use PPG index pairs if available, else ECG
        pairs = seg["ppg_idx"] if seg["ppg_idx"] else seg["ecg_idx"]
        sr = header["ppg_sr"] if seg["ppg_idx"] else header["ecg_sr"]
        for start, end in pairs:
            t0 = (start - 1) / sr / 60.0  # minutes
            t1 = (end - 1) / sr / 60.0
            ax.barh(i, t1 - t0, left=t0, height=0.6, color=f"C{i % 10}", alpha=0.7)

    ax.set_yticks(range(len(segments)))
    ax.set_yticklabels([f"Seg {i}" for i in range(len(segments))])
    ax.set_xlabel("Time (minutes)")
    ax.set_title(f"Annealed Segments Overview — {len(segments)} segments")
    ax.grid(True, axis="x", alpha=0.3)
    plt.tight_layout()
    plt.show()
    return fig


# ── Main ─────────────────────────────────────────────────────────────────────


def main():
    if len(sys.argv) < 2:
        print(
            "Usage: python visualize_annealed.py <path_to_annealed.bin> [segment_index]"
        )
        print("  No index  → overview + all segments")
        print("  index N   → plot only segment N")
        sys.exit(1)

    path = sys.argv[1]
    header, segments = read_annealed_bin(path)

    if not segments:
        print("No segments found in file.")
        sys.exit(0)

    print(f"File: {path}")
    print(f"  Segments:    {header['n_segments']}")
    print(f"  PPG SR:      {header['ppg_sr']} Hz")
    print(f"  ECG SR:      {header['ecg_sr']} Hz")
    print(f"  Epoch size:  {header['epoch_sec']} s")

    for i, seg in enumerate(segments):
        dur_ppg = len(seg["ppg"]) / header["ppg_sr"] if len(seg["ppg"]) else 0
        dur_ecg = len(seg["ecg1"]) / header["ecg_sr"] if len(seg["ecg1"]) else 0
        print(
            f"  Seg {i}: PPG {dur_ppg:.1f}s ({len(seg['ppg'])} samples), "
            f"ECG {dur_ecg:.1f}s ({len(seg['ecg1'])} samples), "
            f"Sleep epochs: {len(seg['sleep'])}"
        )

    if len(sys.argv) >= 3:
        idx = int(sys.argv[2])
        if 0 <= idx < len(segments):
            plot_segment(header, segments[idx], idx)
        else:
            print(f"Segment index {idx} out of range [0, {len(segments) - 1}]")
    else:
        plot_overview(header, segments)
        for i, seg in enumerate(segments):
            plot_segment(header, segments[i], i)


if __name__ == "__main__":
    main()
