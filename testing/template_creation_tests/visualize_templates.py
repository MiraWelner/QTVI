"""
visualize_templates.py

Compare per-bin ECG templates from three sources side-by-side:
  - Daniel:  original MATLAB pipeline (_template_info.mat, ecgTemplate field)
  - Mira:    C++ pipeline (_templates.bin, ch1_raw_template)
  - Deep:    new MATLAB pipeline (_fullresults.mat, ALLavgstdcycAsStruct[i].avgcyc)

Outputs N plots per subject, each containing up to 12 per-bin overlays.
Each source is plotted on its own honest time axis (no display-time
stretching), so a 281-sample Daniel template ends at 281 ms while a
757-sample Deep template ends at 757 ms. All panels in a chunk share
the same x-range so morphology shifts across the recording are visible.

If a Mira `_beats.bin` file is present alongside the `_templates.bin`,
the individual ch1-raw kept beats that went into each bin's mean are
overlaid in translucent gray behind the Mira mean template.
"""

import math
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

try:
    import h5py

    HAVE_H5PY = True
except ImportError:
    HAVE_H5PY = False

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

daniel_path = r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\5_generate_template_files\mesa_templates_daniel"
mira_path = r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\5_generate_template_files\mesa_templates_mira"
deep_path = r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\5_generate_template_files\mesa_templates_deep"

RESULTS_DIR = r"D:\USERS\MiraWelner\QTVI\testing\template_creation_tests"

DANIEL_SR = 1000.0
MIRA_SR = 1000.0
DEEP_SR_DEFAULT = 1000.0

DANIEL_COLOR = "red"
MIRA_COLOR = "black"
DEEP_COLOR = "blue"
BEAT_COLOR = "0.6"  # light gray for individual beats


# ============================================================================
# Mira reader (C++ _templates.bin)
# ============================================================================
def read_mira_template_bin(path):
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

        def skip_double_vec():
            sz = read_u64()
            if sz and sz > 0:
                f.seek(sz * 8, 1)

        def skip_method_block():
            skip_double_vec()
            f.seek(16, 1)  # alignment_point + avg_r_expand

        num_bins = read_u64()
        if num_bins is None:
            return templates

        for i in range(num_bins):
            info = {"index": i}
            info["ch1_raw_template"] = read_double_vec()
            info["ch1_raw_alignment_point"] = read_f64()
            info["ch1_raw_avg_r_expand"] = read_f64()
            for _ in range(11):  # ch1_squared, ch1_absval, ch1_unfiltered, ch2x4, ch3x4
                skip_method_block()
            skip_double_vec()  # PPG
            info["bad_segment"] = bool(read_u8())
            templates.append(info)
    return templates


# ============================================================================
# Mira beats reader (C++ _beats.bin)
#
# Format (matches write_beats_binfile in template_io.cpp):
#   [u64 n_bins]
#   per bin:
#     [u8  bad_segment]
#     [u64 n_beats]
#     per beat:
#       [u64 sz][sz x f64]
#
# Returns a list (one entry per bin) of lists of np.ndarray (one per beat).
# Bad segments come back as empty lists.
# ============================================================================
def read_mira_beats_bin(path):
    per_bin_beats = []
    with open(path, "rb") as f:

        def read_u64():
            d = f.read(8)
            return struct.unpack("<Q", d)[0] if len(d) == 8 else None

        def read_u8():
            d = f.read(1)
            return struct.unpack("<B", d)[0] if len(d) == 1 else None

        n_bins = read_u64()
        if n_bins is None:
            return per_bin_beats

        for _ in range(n_bins):
            bad = read_u8()
            n_beats = read_u64()
            beats = []
            if n_beats is None:
                per_bin_beats.append(beats)
                continue
            for _ in range(n_beats):
                sz = read_u64()
                if sz is None or sz == 0:
                    beats.append(np.array([], dtype=np.float64))
                    continue
                arr = np.frombuffer(f.read(sz * 8), dtype=np.float64).copy()
                beats.append(arr)
            per_bin_beats.append(beats)
    return per_bin_beats


def find_mira_beats_path(templates_path):
    """Given .../<stem>_templates.bin, return .../<stem>_beats.bin if present."""
    p = Path(templates_path)
    stem = p.stem
    if stem.endswith("_templates"):
        stem = stem[: -len("_templates")]
    candidate = p.with_name(f"{stem}_beats.bin")
    return candidate if candidate.exists() else None


# ============================================================================
# Daniel reader (original MATLAB _template_info.mat)
# ============================================================================
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


def read_daniel_template_mat(path):
    mat = sio.loadmat(path, squeeze_me=False)
    if "template_info" not in mat:
        raise ValueError(f"No 'template_info' in {path}")
    cells = mat["template_info"].flatten()
    templates = []
    for i in range(len(cells)):
        cell = cells[i]
        while isinstance(cell, np.ndarray) and cell.dtype == object and cell.ndim == 0:
            cell = cell.flat[0]
        if cell is None or (isinstance(cell, np.ndarray) and cell.size == 0):
            templates.append(None)
            continue
        info = {"index": i}
        try:
            info["bad_segment"] = bool(extract_scalar(cell["bad_segment"]))
        except Exception:
            info["bad_segment"] = False
        try:
            info["ecgTemplate"] = extract_vector(cell["ecgTemplate"])
        except Exception:
            info["ecgTemplate"] = np.array([], dtype=np.float64)
        templates.append(info)
    return templates


# ============================================================================
# Deep reader (new MATLAB _fullresults.mat)
#
# Per-epoch averaged ECG lives in ALLavgstdcycAsStruct[i].avgcyc.
# v7.3 .mat files (HDF5) need h5py; older .mat files use scipy.io.
# ============================================================================
def _is_hdf5_mat(path):
    """v7.3 .mat files start with a 128-byte MATLAB description header
    that includes the substring 'MATLAB 7.3'. The actual HDF5 magic only
    appears starting at byte 128, so a naive 4-byte check against
    b'\\x89HDF' incorrectly returns False for v7.3 .mat files."""
    try:
        with open(path, "rb") as f:
            header = f.read(128)
        return b"MATLAB 7.3" in header
    except OSError:
        return False


def _read_deep_h5py(path):
    """v7.3 .mat: read ALLavgstdcycAsStruct[i].avgcyc lazily via HDF5 refs."""
    out = {"sample_rate": DEEP_SR_DEFAULT, "templates": []}
    with h5py.File(path, "r") as h5:
        if "ALLavgstdcycSampRate" in h5:
            sr = np.asarray(h5["ALLavgstdcycSampRate"]).flatten()
            if sr.size > 0:
                out["sample_rate"] = float(sr[0])

        if "ALLavgstdcycAsStruct" not in h5:
            return out
        struct_grp = h5["ALLavgstdcycAsStruct"]

        # In v7.3 MATLAB stores struct arrays as a group whose fields are
        # datasets of HDF5 references (one per array element).
        if "avgcyc" not in struct_grp:
            return out
        refs = np.asarray(struct_grp["avgcyc"]).ravel()

        for ref in refs:
            try:
                arr = np.asarray(h5[ref]).ravel().astype(np.float64)
                out["templates"].append({"avgcyc": arr})
            except Exception:
                out["templates"].append({"avgcyc": np.array([], dtype=np.float64)})
    return out


def _read_deep_scipy(path):
    """Older .mat: load the whole file, unpack the struct array."""
    out = {"sample_rate": DEEP_SR_DEFAULT, "templates": []}
    mat = sio.loadmat(path, squeeze_me=False)

    if "ALLavgstdcycSampRate" in mat:
        sr = np.asarray(mat["ALLavgstdcycSampRate"]).flatten()
        if sr.size > 0:
            out["sample_rate"] = float(sr[0])

    if "ALLavgstdcycAsStruct" not in mat:
        return out
    struct_arr = mat["ALLavgstdcycAsStruct"].flatten()

    for cell in struct_arr:
        while isinstance(cell, np.ndarray) and cell.dtype == object and cell.ndim == 0:
            cell = cell.flat[0]
        if cell is None or (isinstance(cell, np.ndarray) and cell.size == 0):
            out["templates"].append({"avgcyc": np.array([], dtype=np.float64)})
            continue
        try:
            out["templates"].append({"avgcyc": extract_vector(cell["avgcyc"])})
        except Exception:
            out["templates"].append({"avgcyc": np.array([], dtype=np.float64)})
    return out


def read_deep_template_mat(path):
    if _is_hdf5_mat(path):
        if not HAVE_H5PY:
            raise RuntimeError(
                f"{path} is MATLAB v7.3 (HDF5). Install h5py: pip install h5py"
            )
        return _read_deep_h5py(path)
    return _read_deep_scipy(path)


# ============================================================================
# Scoring
# ============================================================================
def normalize(vec):
    if vec.size == 0:
        return vec
    r = np.ptp(vec)
    return (vec - np.min(vec)) / r if r > 1e-10 else np.zeros_like(vec)


def normalize_to_ref(vec, ref_min, ref_range):
    """Scale `vec` using the same min and range used to normalize a
    reference template. This keeps individual beats on the same vertical
    scale as the mean they were averaged into, so the gray traces sit
    cleanly around the red Mira line instead of drifting up or down."""
    if vec.size == 0 or ref_range < 1e-10:
        return vec
    return (vec - ref_min) / ref_range


def resample_to_length(vec, target_len):
    if vec.size == 0 or target_len <= 0:
        return np.array([], dtype=np.float64)
    if vec.size == target_len:
        return vec
    return scipy_signal.resample(vec, target_len)


def safe_corr(a, b):
    if a.size == 0 or b.size == 0 or a.size != b.size:
        return np.nan
    if np.std(a) < 1e-10 or np.std(b) < 1e-10:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def corr_pair_aligned(a, b, sr_a=1000.0, sr_b=1000.0, max_post_r_ms=None):
    """Correlation between two normalized templates, aligned at the R-peak.

    If `max_post_r_ms` is given, the comparison window is capped at that
    many ms after the R-peak. This matches a visually cropped plot so the
    PCC reflects only the morphology actually shown.
    """
    if a.size == 0 or b.size == 0:
        return np.nan

    r_a = int(np.argmax(a))
    r_b = int(np.argmax(b))

    a_left_ms = -r_a / sr_a * 1000
    a_right_ms = (a.size - 1 - r_a) / sr_a * 1000
    b_left_ms = -r_b / sr_b * 1000
    b_right_ms = (b.size - 1 - r_b) / sr_b * 1000

    overlap_left = max(a_left_ms, b_left_ms)
    overlap_right = min(a_right_ms, b_right_ms)

    # Cap the post-R window so the score reflects only the cropped region.
    if max_post_r_ms is not None:
        overlap_right = min(overlap_right, max_post_r_ms)

    if overlap_right <= overlap_left:
        return np.nan

    a_start = int(round(r_a + overlap_left * sr_a / 1000))
    a_end = int(round(r_a + overlap_right * sr_a / 1000)) + 1
    b_start = int(round(r_b + overlap_left * sr_b / 1000))
    b_end = int(round(r_b + overlap_right * sr_b / 1000)) + 1

    a_crop = a[max(0, a_start) : min(a.size, a_end)]
    b_crop = b[max(0, b_start) : min(b.size, b_end)]

    if a_crop.size < 2 or b_crop.size < 2:
        return np.nan

    n = min(a_crop.size, b_crop.size)
    return safe_corr(
        resample_to_length(a_crop, n),
        resample_to_length(b_crop, n),
    )


def find_r_peak_idx(vec):
    """Return the index of the R-peak in a normalized template.
    Heuristic: the global maximum, since QRS is the largest deflection.
    NaN-safe so it works on raw beats with NaN tails."""
    if vec.size == 0:
        return 0
    if np.all(np.isnan(vec)):
        return 0
    return int(np.nanargmax(vec))


# ============================================================================
# Plot: chunks of 9 bins per file
# ============================================================================
def plot_chunk(
    subject_id,
    daniel_chunk,
    mira_chunk,
    deep_chunk,
    beats_chunk,
    start_idx,
    save_path,
    deep_sr,
    cols_per_row=3,
):
    n = max(len(daniel_chunk), len(mira_chunk), len(deep_chunk))
    rows = math.ceil(n / cols_per_row)

    panel_w, panel_h = 8, 5
    fig, axes = plt.subplots(
        rows,
        cols_per_row,
        figsize=(panel_w * cols_per_row, panel_h * rows),
        squeeze=False,
    )

    # Track the longest x-extent across all panels so we can give the chunk
    # a single shared x-range at the end. Bins of similar morphology should
    # be visually comparable side by side; auto-scaled axes would let panels
    # differ in width and obscure that.
    chunk_max_x = 0.0
    chunk_min_x = 0.0

    legend_added = False

    for i in range(rows * cols_per_row):
        ax = axes[i // cols_per_row, i % cols_per_row]
        ax.set_xticks([])
        ax.set_yticks([])

        if i >= n:
            ax.axis("off")
            continue

        bin_idx = start_idx + i

        d = daniel_chunk[i] if i < len(daniel_chunk) else None
        m = mira_chunk[i] if i < len(mira_chunk) else None
        p = deep_chunk[i] if i < len(deep_chunk) else None
        beats = beats_chunk[i] if i < len(beats_chunk) else []

        d_vec = d["ecgTemplate"] if d is not None else np.array([])
        m_vec = m["ch1_raw_template"] if m is not None else np.array([])
        p_vec = p["avgcyc"] if p is not None else np.array([])

        m_bad = m is not None and m.get("bad_segment", False)

        any_data = d_vec.size > 0 or m_vec.size > 0 or p_vec.size > 0

        if not any_data or m_bad:
            ax.set_facecolor("#f0f0f0")
            label = "no data" if not any_data else "(mira bad)"
            ax.text(
                0.5,
                0.5,
                f"bin {bin_idx}\n{label}",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=10,
                color="gray",
            )
            continue

        # Each source plotted on its OWN time axis -- no display resampling.
        # Daniel/Mira/Deep templates can cover different physical windows
        # (e.g. 281 ms vs 757 ms around the R-peak), and stretching them
        # to a common sample count would visually distort the shorter ones.
        d_n = normalize(d_vec)
        m_n = normalize(m_vec)
        p_n = normalize(p_vec)

        # ---- Individual beats first, so the colored means draw on top ----
        # Each beat is normalized using the SAME min/range as the Mira mean
        # template, not its own min/range. Per-beat normalization would
        # squash every beat into [0, 1] and hide amplitude variation, which
        # is exactly the thing the gray cloud is supposed to show.
        #
        # Alignment: the beats are ALREADY aligned to the same column inside
        # EnsembleTemplate (that's how the median collapses cleanly). The
        # Mira template is a column-wise median of the same matrix, so its
        # argmax IS the alignment column. We use it as t=0 for every beat.
        # Calling find_r_peak_idx per beat would run argmax independently
        # on each one, and for beats where a flanking R-peak happens to be
        # slightly taller than the central one, argmax returns the wrong
        # column, shifting that beat by ~one R-R interval and producing
        # ghost QRS smudges in the gray cloud.
        if m_vec.size > 0 and beats:
            m_min = float(np.nanmin(m_vec))
            m_range = float(np.nanmax(m_vec) - m_min)
            if m_range > 1e-10:
                m_normed_for_align = (m_vec - m_min) / m_range
                r_align = find_r_peak_idx(m_normed_for_align)
                for beat in beats:
                    if beat.size == 0:
                        continue
                    beat_scaled = normalize_to_ref(beat, m_min, m_range)
                    t_b = (np.arange(beat.size) - r_align) / MIRA_SR * 1000
                    ax.plot(
                        t_b,
                        beat_scaled,
                        linewidth=0.5,
                        color=BEAT_COLOR,
                        alpha=0.25,
                        zorder=1,
                    )
                    if t_b.size > 0:
                        chunk_max_x = max(chunk_max_x, float(t_b[-1]))
                        chunk_min_x = min(chunk_min_x, float(t_b[0]))

        if d_n.size > 0:
            r_d = find_r_peak_idx(d_n)
            t_d = (np.arange(d_n.size) - r_d) / DANIEL_SR * 1000
            ax.plot(
                t_d,
                d_n,
                linewidth=1.0,
                color=DANIEL_COLOR,
                alpha=0.85,
                label="Daniel",
                zorder=3,
            )
            chunk_max_x = max(chunk_max_x, float(t_d[-1]))
            chunk_min_x = min(chunk_min_x, float(t_d[0]))

        if m_n.size > 0:
            r_m = find_r_peak_idx(m_n)
            t_m = (np.arange(m_n.size) - r_m) / MIRA_SR * 1000
            ax.plot(
                t_m,
                m_n,
                linewidth=1.0,
                color=MIRA_COLOR,
                alpha=0.85,
                label="Mira",
                zorder=4,
            )
            chunk_max_x = max(chunk_max_x, float(t_m[-1]))
            chunk_min_x = min(chunk_min_x, float(t_m[0]))

        if p_n.size > 0:
            r_p = find_r_peak_idx(p_n)
            t_p = (np.arange(p_n.size) - r_p) / deep_sr * 1000
            ax.plot(
                t_p,
                p_n,
                linewidth=1.0,
                color=DEEP_COLOR,
                alpha=0.85,
                label="Deep",
                zorder=5,
            )
            chunk_max_x = max(chunk_max_x, float(t_p[-1]))
            chunk_min_x = min(chunk_min_x, float(t_p[0]))

        # Pairwise correlations: still need equal length, but truncate
        # rather than stretch so we don't compare invented samples.
        c_dm = corr_pair_aligned(d_n, m_n)
        c_dp = corr_pair_aligned(d_n, p_n)
        c_mp = corr_pair_aligned(m_n, p_n)

        if np.isnan(c_dm):
            tint = "#f0f0f0"
        elif c_dm >= 1.0 - 1e-6:
            tint = "#e8f5e9"  # green: perfect match to 6 decimals
        elif c_dm >= 0.9:
            tint = "#fffde7"  # yellow: close
        else:
            tint = "#ffebee"  # red: diverging
        ax.set_facecolor(tint)

        def fmt(label, c):
            if np.isnan(c):
                return f"{label} -"
            return f"{label} {c:.6f}"

        n_beats = sum(1 for b in beats if b.size > 0) if beats else 0

        title = f"Bin {bin_idx}: {fmt('Daniel-Mira', c_dm)} {fmt('Daniel-Deep', c_dp)} {fmt('Mira-Deep', c_mp)}"
        ax.set_title(title, fontsize=12)
        ax.set_ylim(-0.3, 1.3)  # widened so individual beats aren't clipped

        if not legend_added:
            ax.legend(fontsize=8, loc="upper right")
            legend_added = True

    last_idx = start_idx + n - 1
    fig.suptitle(
        f"{subject_id}  -  bins {start_idx}-{last_idx}",
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.subplots_adjust(wspace=0.30, hspace=0.40)
    fig.savefig(save_path, format="svg", bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)


def plot_in_chunks(
    subject_id,
    daniel_list,
    mira_list,
    deep_list,
    beats_list,
    deep_sr,
    results_dir,
    bins_per_plot=9,
    cols_per_row=3,
):
    n = max(len(daniel_list), len(mira_list), len(deep_list))
    if n == 0:
        print("Nothing to plot.")
        return []

    n_chunks = math.ceil(n / bins_per_plot)
    paths = []
    for chunk_idx in range(n_chunks):
        start = chunk_idx * bins_per_plot
        end = min(start + bins_per_plot, n)
        save_path = os.path.join(
            results_dir,
            f"{subject_id}_bins_{start:04d}-{end - 1:04d}.svg",
        )
        plot_chunk(
            subject_id,
            daniel_list[start:end],
            mira_list[start:end],
            deep_list[start:end],
            beats_list[start:end] if beats_list else [],
            start,
            save_path,
            deep_sr=deep_sr,
            cols_per_row=cols_per_row,
        )
        paths.append(save_path)
    return paths


# ============================================================================
# File matching
#
# Filename conventions per source:
#   Daniel: <subject>_template_info.mat                  -> strip "_template_info"
#   Mira:   <subject>_<rate>_<binlen>_templates.bin      -> strip suffix + last two tokens
#   Deep:   <subject>_part1_ECG_fs<rate>_fullresults.mat -> strip "_part1_ECG_..."
# ============================================================================
def _daniel_subject_id(stem):
    idx = stem.find("_template_info")
    return stem[:idx] if idx >= 0 else None


def _mira_subject_id(stem):
    if not stem.endswith("_templates"):
        return None
    s = stem[: -len("_templates")]
    parts = s.rsplit("_", 2)
    if len(parts) == 3:
        return parts[0]
    return s


def _deep_subject_id(stem):
    if stem.endswith("_fullresults"):
        s = stem[: -len("_fullresults")]
    else:
        s = stem
    idx = s.find("_part1")
    return s[:idx] if idx >= 0 else s


def find_matching_files():
    daniel_files = {}
    for f in Path(daniel_path).glob("*_template_info.mat"):
        sid = _daniel_subject_id(f.stem)
        if sid:
            daniel_files[sid] = f

    mira_files = {}
    for f in Path(mira_path).glob("*_templates.bin"):
        sid = _mira_subject_id(f.stem)
        if sid:
            mira_files[sid] = f

    deep_files = {}
    for f in Path(deep_path).glob("*_fullresults.mat"):
        sid = _deep_subject_id(f.stem)
        if sid:
            deep_files[sid] = f

    all_subjects = sorted(set(daniel_files) | set(mira_files) | set(deep_files))
    return all_subjects, daniel_files, mira_files, deep_files


# ============================================================================
# Main
# ============================================================================
def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Plot Daniel/Mira/Deep template overlays in 12-bin chunks"
    )
    parser.add_argument("subject", type=str, help="Subject ID (e.g. 3010023_20110817)")
    parser.add_argument(
        "--bins-per-plot",
        type=int,
        default=4,
        help="Bins per output file (default: 9, i.e. 3x3)",
    )
    parser.add_argument(
        "--cols",
        type=int,
        default=2,
        help="Panels per row within a plot (default: 3)",
    )
    args = parser.parse_args()

    all_subjects, daniel_files, mira_files, deep_files = find_matching_files()
    if args.subject not in all_subjects:
        print(f"Subject {args.subject} not found.")
        print(f"Available: {', '.join(all_subjects)}")
        return 1

    daniel_list = (
        read_daniel_template_mat(str(daniel_files[args.subject]))
        if args.subject in daniel_files
        else []
    )

    mira_list = []
    beats_list = []
    if args.subject in mira_files:
        mira_templates_path = mira_files[args.subject]
        mira_list = read_mira_template_bin(str(mira_templates_path))

        beats_path = find_mira_beats_path(mira_templates_path)
        if beats_path is not None:
            try:
                beats_list = read_mira_beats_bin(str(beats_path))
            except Exception as e:
                print(f"  WARNING: failed to read beats file {beats_path}: {e}")
                beats_list = []

    if args.subject in deep_files:
        deep_data = read_deep_template_mat(str(deep_files[args.subject]))
        deep_list = deep_data["templates"]
        deep_sr = deep_data.get("sample_rate", DEEP_SR_DEFAULT)
    else:
        deep_list = []
        deep_sr = DEEP_SR_DEFAULT

    print(
        f"  Daniel: {len(daniel_list)} bins  "
        f"({'present' if daniel_list else 'MISSING'})"
    )
    print(f"  Mira:   {len(mira_list)} bins  ({'present' if mira_list else 'MISSING'})")
    if beats_list:
        total_beats = sum(len(b) for b in beats_list)
        print(f"  Beats:  {len(beats_list)} bins, {total_beats} total beats")
    else:
        print("  Beats:  (none — _beats.bin not found or empty)")
    print(
        f"  Deep:   {len(deep_list)} bins  "
        f"({'present' if deep_list else 'MISSING'})  "
        f"sr={deep_sr}"
    )

    os.makedirs(RESULTS_DIR, exist_ok=True)
    paths = plot_in_chunks(
        args.subject,
        daniel_list,
        mira_list,
        deep_list,
        beats_list,
        deep_sr,
        RESULTS_DIR,
        bins_per_plot=args.bins_per_plot,
        cols_per_row=args.cols,
    )
    print(f"\nSaved {len(paths)} plot(s):")
    for p in paths:
        print(f"  {p}")
    return 0


if __name__ == "__main__":
    main()
