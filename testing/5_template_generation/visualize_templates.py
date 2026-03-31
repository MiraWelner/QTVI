"""
visualize_templates.py

Generates side-by-side visualizations of MATLAB vs C++ template outputs.
Saves all plots to a results folder. No interactive display.

Compares C++ ECG raw method against MATLAB (ch1).
Shows squared and absval C++ templates without MATLAB overlay.
Resamples C++ (2000 Hz) to MATLAB (256 Hz) for overlay comparison.
"""

import os
import struct
import sys
from pathlib import Path

import numpy as np
from scipy import signal as scipy_signal

try:
    import scipy.io as sio
except ImportError:
    print("ERROR: scipy required: pip install scipy")
    sys.exit(1)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

MATLAB_DIR = (
    "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\5_generate_template_files\\matlab"
)
CPP_DIR = (
    "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\5_generate_template_files\\mesa"
)
RESULTS_DIR = (
    "D:\\USERS\\MiraWelner\\QTVI\\testing\\5_template_generation\\results\\plots"
)
MATLAB_SR = 256.0
CPP_SR = 2000.0
SR_RATIO = CPP_SR / MATLAB_SR
METHODS = ["raw", "squared", "absval"]
METHOD_COLORS = {"raw": "red", "squared": "orange", "absval": "purple"}


# ============================================================================
# Readers
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
            return np.frombuffer(f.read(sz * 8), dtype=np.float64).copy()

        def read_pair_vec():
            sz = read_u64()
            if sz is None or sz == 0:
                return np.array([], dtype=np.uint64).reshape(0, 2)
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


def extract_scalar(val):
    while isinstance(val, np.ndarray):
        if val.size == 0:
            return None
        if val.size == 1:
            val = val.flat[0]
        else:
            return val
    return val


def extract_vector(val):
    while isinstance(val, np.ndarray) and val.dtype == object:
        if val.size == 1:
            val = val.flat[0]
        else:
            break
    if isinstance(val, np.ndarray):
        return np.array(val, dtype=np.float64).flatten()
    return np.array([], dtype=np.float64)


def read_matlab_template_mat(path):
    mat = sio.loadmat(path, squeeze_me=False)
    if "template_info" not in mat:
        raise ValueError(f"No 'template_info' in {path}")
    raw = mat["template_info"]
    templates = []
    cells = raw.flatten()
    for i in range(len(cells)):
        cell = cells[i]
        while isinstance(cell, np.ndarray) and cell.dtype == object and cell.ndim == 0:
            cell = cell.flat[0]
        if cell is None or (isinstance(cell, np.ndarray) and cell.size == 0):
            templates.append(None)
            continue
        info = {}
        try:
            v = extract_scalar(cell["index"])
            info["index"] = int(v) if v is not None else i
        except:
            info["index"] = i
        try:
            info["bad_segment"] = bool(extract_scalar(cell["bad_segment"]))
        except:
            info["bad_segment"] = False
        for field in ["ecgTemplate", "ppgTemplate"]:
            try:
                info[field] = extract_vector(cell[field])
            except Exception as e:
                info[field] = np.array([], dtype=np.float64)
        try:
            v = extract_scalar(cell["alignment_point"])
            info["alignment_point"] = float(v) if v is not None else 0.0
        except:
            info["alignment_point"] = 0.0
        try:
            v = extract_scalar(cell["avg_r_expand"])
            info["avg_r_expand"] = float(v) if v is not None else 0.0
        except:
            info["avg_r_expand"] = 0.0
        templates.append(info)
    return templates


# ============================================================================
# Scoring
# ============================================================================
def resample_and_corr(mat_vec, cpp_vec):
    if mat_vec.size == 0 or cpp_vec.size == 0:
        return 0.0, np.array([]), np.array([]), np.array([])
    cpp_resampled = scipy_signal.resample(cpp_vec, mat_vec.size)
    mat_r = np.ptp(mat_vec)
    cpp_r = np.ptp(cpp_resampled)

    mat_n = (
        (mat_vec - np.min(mat_vec)) / mat_r if mat_r > 1e-10 else np.zeros_like(mat_vec)
    )
    cpp_n = (
        (cpp_resampled - np.min(cpp_resampled)) / cpp_r
        if cpp_r > 1e-10
        else np.zeros_like(cpp_resampled)
    )

    if np.std(mat_n) > 1e-10 and np.std(cpp_n) > 1e-10:
        corr = np.corrcoef(mat_n, cpp_n)[0, 1]
    else:
        corr = 0.0
    return corr, cpp_resampled, mat_n, cpp_n


def score_bin(mat, cpp):
    scores = {}
    corr, resampled, mn, cn = resample_and_corr(mat["ppgTemplate"], cpp["ppgTemplate"])
    scores["ppgTemplate"] = {
        "corr": corr,
        "mat_len": mat["ppgTemplate"].size,
        "cpp_len": cpp["ppgTemplate"].size,
        "resampled": resampled,
        "mat_norm": mn,
        "cpp_norm": cn,
    }
    for method in METHODS:
        mat_ecg = mat["ecgTemplate"]
        cpp_ecg = cpp[f"ch1_ecgTemplate_{method}"]
        corr, resampled, mn, cn = resample_and_corr(mat_ecg, cpp_ecg)
        scores[f"ecgTemplate_{method}"] = {
            "corr": corr,
            "mat_len": mat_ecg.size,
            "cpp_len": cpp_ecg.size,
            "resampled": resampled,
            "mat_norm": mn,
            "cpp_norm": cn,
        }
    for method in METHODS:
        ap_cpp = cpp[f"ch1_alignment_point_{method}"]
        ar_cpp = cpp[f"ch1_avg_r_expand_{method}"]
        scores[f"alignment_point_{method}"] = {
            "mat": mat["alignment_point"],
            "cpp": ap_cpp,
            "cpp_scaled": ap_cpp / SR_RATIO if SR_RATIO != 0 else 0,
        }
        scores[f"avg_r_expand_{method}"] = {
            "mat": mat["avg_r_expand"],
            "cpp": ar_cpp,
            "cpp_scaled": ar_cpp / SR_RATIO if SR_RATIO != 0 else 0,
        }
    return scores


# ============================================================================
# Plotting
# ============================================================================
def plot_bin(subject_id, bin_idx, mat, cpp, scores, save_path):
    fig = plt.figure(figsize=(20, 14), constrained_layout=True)
    fig.suptitle(f"{subject_id}  —  Bin {bin_idx}", fontsize=18, fontweight="bold")
    gs = GridSpec(3, 4, figure=fig)

    col_configs = []
    for ci, method in enumerate(METHODS):
        col_configs.append(
            (ci, f"ecgTemplate_{method}", f"ECG ({method})", METHOD_COLORS[method])
        )
    col_configs.append((3, "ppgTemplate", "PPG", "blue"))

    TITLE_SIZE = 20
    LABEL_SIZE = 20
    TICK_SIZE = 11
    LEGEND_SIZE = 11

    for col, key, label, color in col_configs:
        sc = scores[key]
        corr = sc["corr"]
        if key == "ppgTemplate":
            mat_vec = mat["ppgTemplate"]
            cpp_vec = cpp["ppgTemplate"]
        else:
            mat_vec = mat["ecgTemplate"]
            cpp_vec = cpp[f"ch1_{key}"]
        mat_norm = sc["mat_norm"]
        cpp_norm = sc["cpp_norm"]

        # Row 0: MATLAB — only for raw ECG (col 0) and PPG (col 3)
        ax0 = fig.add_subplot(gs[0, col])
        if col == 0:
            if mat_vec.size > 0:
                t = np.arange(mat_vec.size) / MATLAB_SR * 1000
                ax0.plot(t, mat_vec, "b-", linewidth=1)
            ax0.set_title(f"MATLAB ECG", fontsize=TITLE_SIZE)
            ax0.set_ylabel("Daniel's Code", fontsize=LABEL_SIZE)
            ax0.tick_params(labelsize=TICK_SIZE)
        elif col == 3:
            if mat_vec.size > 0:
                t = np.arange(mat_vec.size) / MATLAB_SR * 1000
                ax0.plot(t, mat_vec, "b-", linewidth=1)
            ax0.set_title(f"MATLAB PPG", fontsize=TITLE_SIZE)
            ax0.tick_params(labelsize=TICK_SIZE)
        else:
            ax0.axis("off")

        # Row 1: C++ templates
        ax1 = fig.add_subplot(gs[1, col])
        if cpp_vec.size > 0:
            t = np.arange(cpp_vec.size) / CPP_SR * 1000
            ax1.plot(t, cpp_vec, color=color, linewidth=1)
        ax1.set_title(f"C++ {label}", fontsize=TITLE_SIZE)
        if col == 0:
            ax1.set_ylabel("My Code", fontsize=LABEL_SIZE)
        ax1.tick_params(labelsize=TICK_SIZE)

        # Row 2: Normalized overlay — only for raw ECG (col 0) and PPG (col 3)
        ax2 = fig.add_subplot(gs[2, col])
        if col == 0 or col == 3:
            if mat_norm.size > 0 and cpp_norm.size > 0:
                t = np.arange(mat_norm.size) / MATLAB_SR * 1000
                ax2.plot(t, mat_norm, "b-", linewidth=1.2, label="MATLAB", alpha=0.8)
                ax2.plot(
                    t,
                    cpp_norm,
                    "--",
                    color=color,
                    linewidth=1.2,
                    label=f"C++ ({key.split('_')[-1] if '_' in key else 'ppg'})",
                    alpha=0.8,
                )
                ax2.legend(fontsize=LEGEND_SIZE)
            ax2.set_title(
                f"corr={corr:.4f}",
                fontsize=TITLE_SIZE,
                color="green" if corr >= 0.9 else "orange" if corr >= 0.7 else "red",
            )
            ax2.set_xlabel("Time (ms)", fontsize=LABEL_SIZE)
            if col == 0:
                ax2.set_ylabel("Overlayed", fontsize=LABEL_SIZE)
            ax2.tick_params(labelsize=TICK_SIZE)
        else:
            ax2.axis("off")

    fig.savefig(save_path, dpi=300)
    plt.close(fig)


def plot_summary(subject_id, mat_list, cpp_list, save_path):
    n = min(len(mat_list), len(cpp_list))
    corrs = {m: [] for m in METHODS}
    ppg_corrs = []

    for i in range(n):
        mat = mat_list[i]
        cpp = cpp_list[i]
        if mat is None or cpp is None:
            for m in METHODS:
                corrs[m].append(0)
            ppg_corrs.append(0)
            continue
        for method in METHODS:
            c, _, _, _ = resample_and_corr(
                mat["ecgTemplate"], cpp[f"ch1_ecgTemplate_{method}"]
            )
            corrs[method].append(c)
        pc, _, _, _ = resample_and_corr(mat["ppgTemplate"], cpp["ppgTemplate"])
        ppg_corrs.append(pc)

    fig, axes = plt.subplots(4, 1, figsize=(4, 2), sharex=True)
    fig.suptitle(f"{subject_id}  —  Template Correlations ({n} bins)", fontsize=13)
    x = np.arange(n)

    for ax_idx, (method, ax) in enumerate(zip(METHODS, axes[:3])):
        arr = np.array(corrs[method])
        colors = ["green" if c >= 0.9 else "orange" if c >= 0.7 else "red" for c in arr]
        ax.bar(x, arr, color=colors, width=1.0, edgecolor="none")
        ax.axhline(0.9, color="green", linestyle="--", alpha=0.5)
        ax.set_ylabel("Corr")
        above90 = np.sum(arr >= 0.9)
        above80 = np.sum(arr >= 0.8)
        ax.set_title(
            f"ECG {method}  —  ≥0.9: {above90}/{n} ({100 * above90 / n:.1f}%)  "
            f"≥0.8: {above80}/{n} ({100 * above80 / n:.1f}%)  "
            f"mean={np.mean(arr):.3f}",
            fontsize=9,
        )
        ax.set_ylim(-0.2, 1.05)

    ppg_arr = np.array(ppg_corrs)
    colors = ["green" if c >= 0.9 else "orange" if c >= 0.7 else "red" for c in ppg_arr]
    axes[3].bar(x, ppg_arr, color=colors, width=1.0, edgecolor="none")
    axes[3].axhline(0.9, color="green", linestyle="--", alpha=0.5)
    axes[3].set_ylabel("Corr", fontsize=9)
    axes[3].set_xlabel("Bin Index")
    above90 = np.sum(ppg_arr >= 0.9)
    above80 = np.sum(ppg_arr >= 0.8)
    axes[3].set_title(
        f"PPG  —  ≥0.9: {above90}/{n} ({100 * above90 / n:.1f}%)  "
        f"≥0.8: {above80}/{n} ({100 * above80 / n:.1f}%)  "
        f"mean={np.mean(ppg_arr):.3f}",
        fontsize=9,
    )
    axes[3].set_ylim(-0.2, 1.05)

    plt.tight_layout()
    fig.savefig(save_path, dpi=100, bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# File matching
# ============================================================================
def find_matching_files():
    mat_files = {}
    for f in Path(MATLAB_DIR).glob("*_template_info.mat"):
        idx = f.stem.find("_template_info")
        if idx >= 0:
            mat_files[f.stem[:idx]] = f

    cpp_files = {}
    for f in Path(CPP_DIR).glob("*_template_info.bin"):
        idx = f.stem.find("_template_info")
        if idx >= 0:
            cpp_files[f.stem[:idx]] = f

    common = sorted(set(mat_files.keys()) & set(cpp_files.keys()))
    return common, mat_files, cpp_files


# ============================================================================
# Main
# ============================================================================
def main():
    import argparse

    parser = argparse.ArgumentParser(description="Visualize MATLAB vs C++ templates")
    parser.add_argument("subject", type=str, help="Subject ID (e.g. 3010023_20110817)")
    parser.add_argument("bin", type=int, help="Bin index")
    args = parser.parse_args()

    common, mat_files, cpp_files = find_matching_files()

    if args.subject not in common:
        print(f"Subject {args.subject} not found. Available: {', '.join(common)}")
        return 1

    mat_list = read_matlab_template_mat(str(mat_files[args.subject]))
    cpp_list = read_cpp_template_bin(str(cpp_files[args.subject]))
    n = min(len(mat_list), len(cpp_list))

    if args.bin >= n or args.bin < 0:
        print(f"Bin {args.bin} out of range (0-{n - 1})")
        return 1

    mat = mat_list[args.bin]
    cpp = cpp_list[args.bin]
    if mat is None or cpp is None:
        print(f"Bin {args.bin} is None")
        return 1

    os.makedirs(RESULTS_DIR, exist_ok=True)
    scores = score_bin(mat, cpp)
    save_path = os.path.join(RESULTS_DIR, f"{args.subject}_bin_{args.bin:04d}.svg")
    plot_bin(args.subject, args.bin, mat, cpp, scores, save_path)
    print(f"Saved: {save_path}")
    return 0


if __name__ == "__main__":
    main()
