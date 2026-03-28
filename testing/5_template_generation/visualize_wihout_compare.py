"""
visualize_templates_multichannel.py

Visualizes C++ template outputs for Bittium (3 ECG, no PPG) and
CHAOS (3 ECG + PPG). No MATLAB comparison / overlay.

Usage:
    python visualize_templates_multichannel.py <dataset> <subject> <bin>
    dataset: "bittium" or "chaos"
"""

import argparse
import os
import struct
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ============================================================================
# Config — fill in paths
# ============================================================================
BITTIUM_DIR = (
    "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\5_generate_template_files\\bittium"
)

CHAOS_DIR = (
    "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\5_generate_template_files\\chaos"
)

RESULTS_DIR = (
    "D:\\USERS\\MiraWelner\\QTVI\\testing\\5_template_generation\\results\\plots"
)

CPP_SR = 2000.0
METHODS = ["raw", "squared", "absval"]
METHOD_COLORS = {"raw": "red", "squared": "orange", "absval": "purple"}

TITLE_SIZE = 20
LABEL_SIZE = 18
TICK_SIZE = 15
LEGEND_SIZE = 15
SUPTITLE_SIZE = 24


# ============================================================================
# Reader (same binary format as MESA)
# ============================================================================
def read_cpp_template_bin(path):
    templates = []
    with open(path, "rb") as f:

        def read_u64():
            d = f.read(8)
            return struct.unpack("<Q", d)[0] if len(d) == 8 else None

        def read_f64():
            d = f.read(8)
            return struct.unpack("<d", d)[0] if len(d) == 8 else None

        def read_u8():
            d = f.read(1)
            return struct.unpack("<B", d)[0] if len(d) == 1 else None

        def read_double_vec():
            sz = read_u64()
            if sz is None or sz == 0:
                return np.array([], dtype=np.float64)
            if sz > 50000000:
                raise ValueError(f"bad double_vec size: {sz} at pos {f.tell()}")
            return np.frombuffer(f.read(sz * 8), dtype=np.float64).copy()

        def read_pair_vec():
            sz = read_u64()
            if sz is None or sz == 0:
                return np.array([], dtype=np.uint64).reshape(0, 2)
            if sz > 50000000:
                raise ValueError(f"bad pair_vec size: {sz} at pos {f.tell()}")
            return np.frombuffer(f.read(sz * 16), dtype=np.uint64).reshape(sz, 2).copy()

        num_bins = read_u64()
        if num_bins is None:
            return templates
        for _ in range(num_bins):
            info = {}
            info["index"] = read_u64()
            info["ppg_bin_indexs"] = read_pair_vec()
            info["ecg_bin_indexs"] = read_pair_vec()
            info["bad_segment"] = bool(read_u8())
            for ch in ["ch1", "ch2", "ch3"]:
                info[f"{ch}_ecgTemplate_raw"] = read_double_vec()
                info[f"{ch}_ecgTemplate_squared"] = read_double_vec()
                info[f"{ch}_ecgTemplate_absval"] = read_double_vec()
                info[f"{ch}_alignment_point_raw"] = read_f64()
                info[f"{ch}_alignment_point_squared"] = read_f64()
                info[f"{ch}_alignment_point_absval"] = read_f64()
                info[f"{ch}_avg_r_expand_raw"] = read_f64()
                info[f"{ch}_avg_r_expand_squared"] = read_f64()
                info[f"{ch}_avg_r_expand_absval"] = read_f64()
            info["ppgTemplate"] = read_double_vec()
            templates.append(info)
    return templates


# ============================================================================
# Plotting
# ============================================================================
def plot_bin_multichannel(subject_id, dataset, bin_idx, cpp, save_path, has_ppg):
    channels = ["ch1", "ch2", "ch3"]
    n_cols = len(METHODS)
    if has_ppg:
        n_cols += 1

    n_rows = len(channels)

    fig = plt.figure(figsize=(6 * n_cols, 5 * n_rows), constrained_layout=True)
    fig.suptitle(
        f"{dataset.upper()}  —  {subject_id}  —  Bin {bin_idx}",
        fontsize=SUPTITLE_SIZE,
        fontweight="bold",
    )
    gs = GridSpec(n_rows, n_cols, figure=fig)

    for row_idx, ch in enumerate(channels):
        for col_idx, method in enumerate(METHODS):
            ax = fig.add_subplot(gs[row_idx, col_idx])
            key = f"{ch}_ecgTemplate_{method}"
            vec = cpp[key]

            if vec.size > 0:
                t = np.arange(vec.size) / CPP_SR * 1000
                ax.plot(t, vec, color=METHOD_COLORS[method], linewidth=1)

            if row_idx == 0:
                ax.set_title(f"ECG ({method})", fontsize=TITLE_SIZE)
            if col_idx == 0:
                ax.set_ylabel(f"Ch {row_idx + 1}\nAmplitude", fontsize=LABEL_SIZE)
            if row_idx == n_rows - 1:
                ax.set_xlabel("Time (ms)", fontsize=LABEL_SIZE)
            ax.tick_params(labelsize=TICK_SIZE)

        # PPG column (only for CHAOS)
        if has_ppg:
            ax = fig.add_subplot(gs[row_idx, len(METHODS)])
            if row_idx == 0:
                vec = cpp["ppgTemplate"]
                if vec.size > 0:
                    t = np.arange(vec.size) / CPP_SR * 1000
                    ax.plot(t, vec, "b-", linewidth=1)
                ax.set_title("PPG", fontsize=TITLE_SIZE)
                ax.set_ylabel("Amplitude", fontsize=LABEL_SIZE)
                ax.tick_params(labelsize=TICK_SIZE)
            else:
                ax.axis("off")

    fig.savefig(save_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {save_path}")


def plot_summary_multichannel(subject_id, dataset, cpp_list, save_path, has_ppg):
    n = len(cpp_list)
    channels = ["ch1", "ch2", "ch3"]

    n_plot_rows = len(channels) * len(METHODS)
    if has_ppg:
        n_plot_rows += 1

    fig, axes = plt.subplots(n_plot_rows, 1, figsize=(14, 3 * n_plot_rows), sharex=True)
    fig.suptitle(
        f"{dataset.upper()}  —  {subject_id}  —  Template Lengths ({n} bins)",
        fontsize=SUPTITLE_SIZE,
    )

    x = np.arange(n)
    ax_idx = 0

    for ch in channels:
        for method in METHODS:
            key = f"{ch}_ecgTemplate_{method}"
            lengths = [cpp[key].size if cpp is not None else 0 for cpp in cpp_list]
            ax = axes[ax_idx]
            colors = ["green" if l > 0 else "red" for l in lengths]
            ax.bar(x, lengths, color=colors, width=1.0, edgecolor="none")
            n_valid = sum(1 for l in lengths if l > 0)
            ax.set_title(
                f"{ch} {method}  —  valid: {n_valid}/{n} ({100 * n_valid / n:.1f}%)",
                fontsize=TITLE_SIZE - 4,
            )
            ax.set_ylabel("Length", fontsize=LABEL_SIZE - 4)
            ax.tick_params(labelsize=TICK_SIZE - 2)
            ax_idx += 1

    if has_ppg:
        lengths = [
            cpp["ppgTemplate"].size if cpp is not None else 0 for cpp in cpp_list
        ]
        ax = axes[ax_idx]
        colors = ["green" if l > 0 else "red" for l in lengths]
        ax.bar(x, lengths, color=colors, width=1.0, edgecolor="none")
        n_valid = sum(1 for l in lengths if l > 0)
        ax.set_title(
            f"PPG  —  valid: {n_valid}/{n} ({100 * n_valid / n:.1f}%)",
            fontsize=TITLE_SIZE - 4,
        )
        ax.set_ylabel("Length", fontsize=LABEL_SIZE - 4)
        ax.tick_params(labelsize=TICK_SIZE - 2)

    axes[-1].set_xlabel("Bin Index", fontsize=LABEL_SIZE)
    plt.tight_layout()
    fig.savefig(save_path, dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {save_path}")


# ============================================================================
# Main
# ============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Visualize Bittium/CHAOS C++ template outputs (3 channels)"
    )
    parser.add_argument(
        "dataset",
        type=str,
        choices=["bittium", "chaos"],
        help="Dataset type",
    )
    parser.add_argument("subject", type=str, help="Subject ID")
    parser.add_argument("bin", type=int, help="Bin index (-1 for summary only)")
    args = parser.parse_args()

    has_ppg = args.dataset == "chaos"
    data_dir = BITTIUM_DIR if args.dataset == "bittium" else CHAOS_DIR

    if not data_dir:
        print(
            f"ERROR: Set {'BITTIUM_DIR' if args.dataset == 'bittium' else 'CHAOS_DIR'} path first"
        )
        return 1

    # Find the file
    bin_file = None
    for f in Path(data_dir).glob("*_template_info.bin"):
        idx = f.stem.find("_template_info")
        if idx >= 0 and f.stem[:idx] == args.subject:
            bin_file = f
            break

    if bin_file is None:
        available = []
        for f in Path(data_dir).glob("*_template_info.bin"):
            idx = f.stem.find("_template_info")
            if idx >= 0:
                available.append(f.stem[:idx])
        print(f"Subject '{args.subject}' not found in {data_dir}")
        if available:
            print(f"Available: {', '.join(sorted(available))}")
        return 1

    cpp_list = read_cpp_template_bin(str(bin_file))
    n = len(cpp_list)
    print(f"Loaded {n} bins from {bin_file.name}")

    os.makedirs(RESULTS_DIR or "results", exist_ok=True)
    results_dir = RESULTS_DIR or "results"

    # Summary plot
    summary_path = os.path.join(
        results_dir, f"{args.dataset}_{args.subject}_summary.svg"
    )
    plot_summary_multichannel(
        args.subject, args.dataset, cpp_list, summary_path, has_ppg
    )

    # Single bin plot
    if args.bin >= 0:
        if args.bin >= n:
            print(f"Bin {args.bin} out of range (0-{n - 1})")
            return 1

        cpp = cpp_list[args.bin]
        if cpp is None:
            print(f"Bin {args.bin} is None")
            return 1

        save_path = os.path.join(
            results_dir,
            f"{args.dataset}_{args.subject}_bin_{args.bin:04d}.svg",
        )
        plot_bin_multichannel(
            args.subject, args.dataset, args.bin, cpp, save_path, has_ppg
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
