"""
compare_single.py

Compare R-peak detection across all sources for a single bin of one
subject:
    - Mira at 500/1000/2000 Hz (C++ {id}_{rate}_wave_markings.bin)
    - Daniel (MATLAB _wave_data.mat)
    - Deep (QTVI per-subject whole-case _wholecaseRRiQTi.csv)

Outputs (all in OUT_DIR):
    <id>_bin<N>_overlay.png        14x7 wide overlay, full window (T_START..T_END)
    <id>_bin<N>_overlay_zoom.png   7x7 square overlay, zoomed window
                                   (ZOOM_T_START..ZOOM_T_END)
    <id>_bin<N>_results.csv        R-peak times (sec) per source side by side,
                                   with all-pairs SSD columns appended on first row

Usage:
    python compare_single.py <subject_id> <bin_number> [--raw]
    e.g. python compare_single.py 3010023 33
         python compare_single.py 3010023 33 --raw
"""

import csv
import os
import struct
import sys
from itertools import combinations, zip_longest
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

# ============================================================================
# Config
# ============================================================================

BIN_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\output_mira\mesa\r_peak_finding_output"
)
MAT_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\output_daniel\r_peak_finding_output"
)
CSV_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\output_deep\mesa\r_peak_finding_output"
)
DEEP_FILE_SUFFIX = "_wholecaseRRiQTi.csv"
OUT_DIR = Path(r"D:\USERS\MiraWelner\QTVI\testing\r_peak_finding_tests\single_bin")

T_START = 0
T_END = 300
ZOOM_T_START = 20
ZOOM_T_END = 30

# Mira files come at three sample rates. Each rate has its own .bin file
# named {id}_{rate}_wave_markings.bin and its own native fs for slicing,
# peak-time conversion, and SSD scoring.
MIRA_RATES = (500, 1000, 2000)

BIN_LENGTH = "005"

# Visual style per Mira rate. Daniel/Deep stay constant (red square / green
# triangle) to match the previous look.
MIRA_COLORS = {500: "tab:blue", 1000: "black", 2000: "tab:cyan"}
MIRA_MARKERS = {500: "o", 1000: "p", 2000: "D"}  # circle, pentagon, diamond

# Deep timestamps are wall-clock seconds and need an fs to convert into
# sample indices. The conversion has to happen in the same fs space as
# whichever Mira file we use to derive bin boundaries; we anchor on 1000 Hz
# when available, falling back to whatever rate is present.
DEEP_ANCHOR_RATE_PREFERENCE = (1000, 2000, 500)

SIGNAL_SCATTER_SIZE = 3
MISS_PENALTY_S = 1.0
MATCH_TOLERANCE_S = 0.15
RAW_ECG_CH1_SLOT = 1
PASSTHROUGH_NUM_CHANNELS = 40
RAW_NATIVE_SR = 256.0
MAX_SANE_IDX = 50_000_000
MAX_SANE_SIG = 200_000_000


# ============================================================================
# Mira reader
# ============================================================================


def read_wave_bin(path, max_bins=None):
    """Read up to `max_bins` bins from a wave_markings.bin.

    Optimization: only the LAST bin we read fully populates its signal
    arrays. Earlier bins record `len(ch1_raw_sig)` (needed by the Deep
    boundary derivation) but skip past every other heavy field. This
    turns "load bin N" from O(N) full parses into O(N) cheap header
    walks — the difference between minutes and seconds at large N.
    """
    bins = []
    file_size = os.path.getsize(path)
    with open(path, "rb") as f:
        num_bins = struct.unpack("<Q", f.read(8))[0]
        if max_bins is not None:
            num_bins = min(num_bins, max_bins)
        for i in range(num_bins):
            if f.tell() >= file_size:
                break
            b = {}
            keep = i == num_bins - 1
            try:
                _read_bin_contents(f, b, keep=keep)
            except (ValueError, struct.error) as e:
                print(f"  ERROR at bin {i}: {e}")
                break
            bins.append(b)
    return bins


def _read_bin_contents(f, b, keep=True):
    """Parse one bin's contents starting at the current file offset.

    If `keep` is False, we skip past every field we don't need for the
    Deep-boundary derivation (which only needs `len(ch1_raw_sig)` per
    bin). This is the difference between allocating ~1 GB of Python
    floats per skipped bin and just doing some seeks.
    """

    def read_idx():
        sz = struct.unpack("<Q", f.read(8))[0]
        if sz == 0:
            return []
        if sz > MAX_SANE_IDX:
            raise ValueError(f"idx count {sz} exceeds limit")
        if keep:
            return list(struct.unpack(f"<{sz}Q", f.read(sz * 8)))
        f.seek(sz * 8, 1)
        return []

    def read_signal_np():
        """Read a signal block. Returns a numpy array if keep=True;
        otherwise seeks past it and returns the sample count."""
        sz = struct.unpack("<Q", f.read(8))[0]
        if sz > MAX_SANE_SIG:
            raise ValueError(f"signal size {sz} exceeds limit")
        if keep:
            if sz == 0:
                return np.empty(0, dtype=np.float64)
            return np.frombuffer(f.read(sz * 8), dtype=np.float64).copy()
        f.seek(sz * 8, 1)
        return sz  # return count so callers can record length without parsing

    def skip_signal():
        sz = struct.unpack("<Q", f.read(8))[0]
        if sz > MAX_SANE_SIG:
            raise ValueError(
                f"signal size {sz} exceeds limit (file desync — wrong channel count?)"
            )
        f.seek(sz * 8, 1)

    def skip_raw_pairs():
        n = struct.unpack("<Q", f.read(8))[0]
        if n > MAX_SANE_IDX:
            raise ValueError(
                f"raw pairs count {n} exceeds limit (file desync — wrong channel count?)"
            )
        f.seek(n * 16, 1)

    def read_raw_pairs():
        n = struct.unpack("<Q", f.read(8))[0]
        if n == 0:
            return np.empty((0, 2), dtype=np.float64)
        if n > MAX_SANE_IDX:
            raise ValueError(f"raw pairs count {n} exceeds limit")
        if keep:
            return np.frombuffer(f.read(n * 16), dtype=np.float64).reshape(n, 2).copy()
        f.seek(n * 16, 1)
        return np.empty((0, 2), dtype=np.float64)

    def skip_pair_vec():
        sz = struct.unpack("<Q", f.read(8))[0]
        if sz > MAX_SANE_IDX:
            raise ValueError(
                f"pair vec size {sz} exceeds limit (file desync — wrong channel count?)"
            )
        f.seek(sz * 16, 1)

    # 11 idx arrays (R-peak indices). Only kept bins materialize them.
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

    # Raw signals. ppg + 3 ECG channels. Earlier bins only need to know
    # ch1_raw_sig's length (Deep cumulates lengths to find bin offsets).
    if keep:
        b["ppgSignal"] = read_signal_np()
        b["ch1_raw_sig"] = read_signal_np()
        b["ch2_raw_sig"] = read_signal_np()
        b["ch3_raw_sig"] = read_signal_np()
    else:
        skip_signal()  # ppgSignal
        ch1_len = read_signal_np()  # ch1_raw_sig: count only
        b["ch1_raw_sig"] = _LenProxy(ch1_len)  # supports len() cheaply
        skip_signal()  # ch2_raw_sig
        skip_signal()  # ch3_raw_sig

    # 6 preprocessed signals: always skipped (not used by this script).
    for _ in range(6):
        skip_signal()
    f.read(9)  # 9 noise flag bytes
    num_pairs = struct.unpack("<Q", f.read(8))[0]
    if 0 < num_pairs <= MAX_SANE_IDX:
        f.seek(num_pairs * 16, 1)

    # ppg_bin_indexs and ecg_bin_indexs: each is a uint64 count followed by
    # `count` (start, end) pairs at the file's signal_rate. We only need
    # ecg_bin_indexs[0].first for cross-rate alignment (gives us the bin's
    # window start in signal-rate sample units; caller divides by signal_rate
    # to get absolute seconds-from-recording-start). PPG indexs are skipped.
    skip_pair_vec()  # ppg_bin_indexs
    if keep:
        ecg_pair_count = struct.unpack("<Q", f.read(8))[0]
        if ecg_pair_count > MAX_SANE_IDX:
            raise ValueError(f"ecg_bin_indexs count {ecg_pair_count} exceeds limit")
        if ecg_pair_count > 0:
            # First pair is (start_sample, end_sample) as uint64.
            first_start = struct.unpack("<Q", f.read(8))[0]
            f.seek(8, 1)  # skip the end half
            # Skip the rest of the pairs.
            if ecg_pair_count > 1:
                f.seek((ecg_pair_count - 1) * 16, 1)
            b["ecg_bin_first_sample"] = first_start
        else:
            b["ecg_bin_first_sample"] = None
    else:
        skip_pair_vec()  # ecg_bin_indexs

    # Per-channel pass-through: { upsampled, raw (t,v) pairs } x N channels.
    # Only the kept bin's slot-1 raw pairs are loaded; all upsampled blocks
    # and other slots' raw pairs are skipped.
    b["ch1_raw_pairs"] = np.empty((0, 2), dtype=np.float64)
    for ch in range(PASSTHROUGH_NUM_CHANNELS):
        skip_signal()
        if keep and ch == RAW_ECG_CH1_SLOT:
            b["ch1_raw_pairs"] = read_raw_pairs()
        else:
            skip_raw_pairs()


class _LenProxy:
    """Tiny stand-in for a sample array that only knows its length.
    Used for non-target bins so Deep's bin-offset accumulation can call
    len() without us having parsed the full signal."""

    __slots__ = ("_n",)

    def __init__(self, n):
        self._n = int(n)

    def __len__(self):
        return self._n


# ============================================================================
# Daniel reader
# ============================================================================


def _is_hdf5_mat(path):
    """v7.3 .mat files are HDF5 underneath, identified by the magic bytes."""
    try:
        with open(path, "rb") as f:
            return f.read(8)[:4] == b"\x89HDF"
    except OSError:
        return False


def _read_daniel_bin_h5py(mat_path, bin_idx):
    """Fast path for v7.3 .mat: lazy-load only the requested cell.

    v7.3 stores cell arrays as HDF5 references; we resolve the bin_idx-th
    reference and read only that group's fields. This avoids loading the
    full overnight signal into memory just to get one minute of it.
    """
    import h5py

    empty = (
        np.array([], dtype=np.intp),
        np.array([], dtype=np.float64),
        np.array([], dtype=np.float64),
    )
    with h5py.File(str(mat_path), "r") as h5:
        if "wave_data" not in h5:
            return empty
        wd = h5["wave_data"]
        # wave_data is a dataset of HDF5 object references; flatten to 1-D.
        refs = np.asarray(wd).ravel()
        if bin_idx >= len(refs):
            return empty
        cell_grp = h5[refs[bin_idx]]

        # cell_grp is the struct for this bin. Read only the fields we need.
        def get_field(name):
            if name not in cell_grp:
                return np.array([], dtype=np.float64)
            arr = np.asarray(cell_grp[name]).ravel().astype(np.float64)
            return arr

        r_f = get_field("ecgRIndex")
        if "ecgSeg" in cell_grp:
            ecg = get_field("ecgSeg")
        elif "ecg" in cell_grp:
            ecg = get_field("ecg")
        else:
            ecg = get_field("ecgSignal")
    return r_f.astype(np.intp), r_f, ecg


def read_daniel_bin(mat_path, bin_idx):
    """Read one bin from Daniel's wave_data.mat.

    Uses h5py for v7.3 files (lazy, fast) and falls back to scipy.io for
    older formats (loads the entire file into memory)."""
    if _is_hdf5_mat(mat_path):
        try:
            return _read_daniel_bin_h5py(mat_path, bin_idx)
        except Exception as e:
            print(f"  WARNING: h5py read failed ({e}); falling back to loadmat")

    mat = sio.loadmat(str(mat_path), squeeze_me=False)
    wave_data = mat["wave_data"]
    cells = list(wave_data.flat)
    empty = (
        np.array([], dtype=np.intp),
        np.array([], dtype=np.float64),
        np.array([], dtype=np.float64),
    )
    if bin_idx >= len(cells):
        return empty
    cell = cells[bin_idx]
    if cell is None or (hasattr(cell, "size") and cell.size == 0):
        return empty
    obj = cell
    while hasattr(obj, "shape") and obj.shape == (1, 1):
        obj = obj[0, 0]
    if not (hasattr(obj, "dtype") and obj.dtype.names):
        return empty

    def get_field(name):
        if name not in obj.dtype.names:
            return np.array([], dtype=np.float64)
        val = obj[name]
        while hasattr(val, "shape") and val.ndim > 1 and 1 in val.shape:
            val = val.squeeze()
        return np.array(val, dtype=np.float64).flatten()

    r_f = get_field("ecgRIndex")
    if "ecgSeg" in obj.dtype.names:
        ecg = get_field("ecgSeg")
    elif "ecg" in obj.dtype.names:
        ecg = get_field("ecg")
    else:
        ecg = get_field("ecgSignal")
    return r_f.astype(np.intp), r_f, ecg


# ============================================================================
# Deep reader
# ============================================================================


def _parse_deep_timestamps_s(csv_path):
    """Read column 1 (0-based) as absolute timestamps in seconds.
    File is comma-separated; header row skipped because it fails float-parse."""
    timestamps = []
    with open(csv_path, "r", newline="", encoding="utf-8-sig") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) < 2:
                continue
            try:
                timestamps.append(float(parts[1].strip()))
            except (ValueError, IndexError):
                continue
    return np.array(timestamps, dtype=np.float64)


def _derive_deep_global_peaks(csv_path, fs):
    """Timestamps are seconds-from-midnight, same reference as Mira.
    Multiply by fs to get global sample indices in that fs's space."""
    if csv_path is None:
        return np.array([], dtype=np.intp)
    csv_path = Path(csv_path)
    if not csv_path.exists() or not csv_path.is_file():
        return np.array([], dtype=np.intp)
    timestamps_s = _parse_deep_timestamps_s(csv_path)
    if len(timestamps_s) == 0:
        return np.array([], dtype=np.intp)
    return np.round(timestamps_s * fs).astype(np.intp)


def read_deep_bin(csv_path, bin_idx, mira_bins, fs):
    """Convert Deep's wall-clock peaks to bin-local sample indices in
    the same fs space as the supplied mira_bins. fs must match the rate
    that mira_bins was produced at."""
    if csv_path is None or not mira_bins:
        return np.array([], dtype=np.intp)
    if bin_idx < 0 or bin_idx >= len(mira_bins):
        return np.array([], dtype=np.intp)
    global_peaks = _derive_deep_global_peaks(csv_path, fs)
    if len(global_peaks) == 0:
        return np.array([], dtype=np.intp)
    cum = 0
    bin_offsets = []
    for b in mira_bins:
        bin_offsets.append(cum)
        cum += len(b.get("ch1_raw_sig", []))
    bin_start = bin_offsets[bin_idx]
    bin_length = len(mira_bins[bin_idx].get("ch1_raw_sig", []))
    if bin_length <= 0:
        return np.array([], dtype=np.intp)
    mask = (global_peaks >= bin_start) & (global_peaks < bin_start + bin_length)
    return (global_peaks[mask] - bin_start).astype(np.intp)


# ============================================================================
# Plot helpers
# ============================================================================


def _slice_window(signal, fs, t_start, t_end):
    n = len(signal)
    if n == 0:
        return 0, 0, np.array([]), np.array([])
    t_end_eff = n / fs if t_end is None else t_end
    start_i = max(0, min(int(t_start * fs), n - 1))
    end_i = max(0, min(int(t_end_eff * fs), n))
    if end_i <= start_i:
        return start_i, end_i, np.array([]), np.array([])
    return start_i, end_i, np.arange(start_i, end_i) / fs, signal[start_i:end_i]


def plot_overlay(
    signals,
    peak_sets,
    t_start,
    t_end,
    title,
    out_path,
    figsize=(7, 7),
    raw_pairs=None,
):
    """signals: list of (sig, color, label, fs, bin_offset). Each signal is
    sliced and drawn at its own fs; bin_offset (seconds) shifts the x axis
    so the bin starts at t=0 in plot coordinates.
    peak_sets: list of (peaks, color, marker, label, sig_idx). Peaks are
    sample indices in the fs of signals[sig_idx]; they pick up the same
    bin_offset as that signal.
    raw_pairs: optional (pairs_array, t_origin) tuple. The pairs array is
    Nx2 absolute (t, v); t_origin is subtracted to put it in plot coords."""
    fig, ax = plt.subplots(figsize=figsize)

    if not signals or all(len(s[0]) == 0 for s in signals):
        ax.text(
            0.5,
            0.5,
            "No ECG signal available",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=14,
            color="gray",
        )
    else:
        n_sigs = len(signals)
        y_lo, y_hi = np.inf, -np.inf
        t_end_max = 0.0

        for idx, sig_tuple in enumerate(signals):
            sig, color, label, fs = sig_tuple[:4]
            bin_offset = sig_tuple[4] if len(sig_tuple) > 4 else 0.0
            sig = np.asarray(sig, dtype=np.float64)
            if len(sig) == 0:
                continue
            # Slice in bin-local time (sig is indexed from its own bin start),
            # then add bin_offset to put it on the shared plot x axis.
            sig_t_start = t_start - bin_offset
            sig_t_end = (t_end - bin_offset) if t_end is not None else None
            si, ei, ts, ss = _slice_window(sig, fs, sig_t_start, sig_t_end)
            if ei > 0:
                t_end_max = max(t_end_max, ei / fs + bin_offset)
            # Dense signal trace: a thin line is visually identical to a
            # 1-pixel scatter at this density, and renders 50-100x faster.
            ax.plot(
                ts + bin_offset,
                ss,
                color=color,
                linewidth=0.5,
                label="_nolegend_",
                rasterized=True,
                zorder=2 + (n_sigs - 1 - idx),
            )
            if len(ss) > 0:
                y_lo = min(y_lo, float(np.percentile(ss, 0.1)))
                y_hi = max(y_hi, float(np.percentile(ss, 99.9)))

        if raw_pairs is not None:
            # raw_pairs may be just an Nx2 array (back-compat) or a
            # (pairs, t_origin) tuple. The tuple form lets callers anchor
            # the raw overlay to the same shared t=0 as the signals.
            if isinstance(raw_pairs, tuple):
                rp_data, rp_origin = raw_pairs
            else:
                rp_data, rp_origin = raw_pairs, None
            rp = np.asarray(rp_data, dtype=np.float64)
            if len(rp) > 0:
                t_vals, v_vals = rp[:, 0], rp[:, 1]
                # If we were given an explicit origin, anchor to it; otherwise
                # fall back to the legacy "subtract first sample" behavior.
                origin = (
                    rp_origin
                    if rp_origin is not None
                    else (t_vals[0] if len(t_vals) > 0 else 0.0)
                )
                t_sec = t_vals - origin
                mask = (t_sec >= t_start) & (
                    t_sec <= (t_end if t_end is not None else t_sec.max() + 1)
                )
                if mask.any():
                    ax.scatter(
                        t_sec[mask],
                        v_vals[mask],
                        s=SIGNAL_SCATTER_SIZE,
                        c="purple",
                        edgecolors="none",
                        label=f"raw ({mask.sum()})",
                        rasterized=True,
                        zorder=2 + n_sigs,
                    )
                    y_lo = min(y_lo, float(np.percentile(v_vals[mask], 0.1)))
                    y_hi = max(y_hi, float(np.percentile(v_vals[mask], 99.9)))

        if not np.isfinite(y_lo):
            y_lo, y_hi = -1.0, 1.0

        for i, (peaks, color, marker, label, sig_idx) in enumerate(peak_sets):
            if len(peaks) == 0 or sig_idx >= len(signals):
                continue
            ref_tuple = signals[sig_idx]
            ref_signal = np.asarray(ref_tuple[0], dtype=np.float64)
            ref_fs = ref_tuple[3]
            ref_offset = ref_tuple[4] if len(ref_tuple) > 4 else 0.0
            if len(ref_signal) == 0:
                continue
            # Peaks are indexed into ref_signal, which is bin-local. Convert
            # the (plot) t_start/t_end window back into the signal's local
            # time range before sample-index masking.
            sig_t_start = t_start - ref_offset
            sig_t_end = (t_end - ref_offset) if t_end is not None else None
            si, ei, _, _ = _slice_window(ref_signal, ref_fs, sig_t_start, sig_t_end)
            peaks_arr = np.asarray(peaks, dtype=np.intp)
            mask = (peaks_arr >= si) & (peaks_arr < ei) & (peaks_arr < len(ref_signal))
            visible = peaks_arr[mask]
            if len(visible) == 0:
                continue
            ax.scatter(
                visible / ref_fs + ref_offset,
                ref_signal[visible],
                marker=marker,
                facecolors="none",
                edgecolors=color,
                # When the same beat is detected by 500/1000/2000Hz they
                # land on top of each other; vary size by enumeration order
                # so each marker's outline peeks out from the next.
                s=260 - i * 30,
                linewidths=1.8,
                zorder=20 + i,
                label=f"{label} ({len(visible)} / {len(peaks)})",
                rasterized=True,
            )
            y_lo = min(y_lo, float(ref_signal[visible].min()))
            y_hi = max(y_hi, float(ref_signal[visible].max()))

        pad = (y_hi - y_lo) * 0.15 if y_hi > y_lo else 0.1
        ax.set_ylim(y_lo - pad, y_hi + pad)
        ax.set_xlim(t_start, t_end if t_end is not None else t_end_max)
        ax.legend(loc="lower right", fontsize=10)

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("ECG")
    ax.set_title(title, fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, format="png", bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# SSD
# ============================================================================


def pairwise_ssd(
    a,
    b,
    fs_a,
    fs_b,
    offset_a=0.0,
    offset_b=0.0,
    tolerance_s=MATCH_TOLERANCE_S,
    miss_penalty_s=MISS_PENALTY_S,
):
    """SSD between two peak sets, each in their own fs's sample-index space.
    Each side is converted to absolute seconds (offset + index/fs) before
    matching, so a 500 Hz Mira and a 1000 Hz Mira can be compared
    meaningfully across the per-rate bin-start quantization gap."""
    a_arr = np.sort(np.asarray(a, dtype=np.float64) / fs_a + offset_a)
    b_arr = np.sort(np.asarray(b, dtype=np.float64) / fs_b + offset_b)
    na, nb = len(a_arr), len(b_arr)
    if na == 0 and nb == 0:
        return float("nan")
    if na == 0:
        return float(nb * miss_penalty_s**2)
    if nb == 0:
        return float(na * miss_penalty_s**2)
    pos = np.searchsorted(b_arr, a_arr)
    if nb > 1:
        pos = np.clip(pos, 1, nb - 1)
        left, right = b_arr[pos - 1], b_arr[pos]
        nearest_idx = np.where((a_arr - left) <= (right - a_arr), pos - 1, pos)
    else:
        nearest_idx = np.zeros_like(pos)
    diff = a_arr - b_arr[nearest_idx]
    within = np.abs(diff) <= tolerance_s
    candidates = np.where(within)[0]
    order = candidates[np.argsort(np.abs(diff[candidates]))]
    b_taken = np.zeros(nb, dtype=bool)
    matched_diffs = []
    for ai in order:
        bi = int(nearest_idx[ai])
        if not b_taken[bi]:
            b_taken[bi] = True
            matched_diffs.append(diff[ai])
    matched_diffs = np.asarray(matched_diffs, dtype=np.float64)
    n_matched = len(matched_diffs)
    ssd = float(np.sum(matched_diffs**2))
    ssd += (na - n_matched + nb - n_matched) * miss_penalty_s**2
    return ssd


# ============================================================================
# Combined results CSV
# ============================================================================


def write_results_csv(out_path, sources):
    """Write peak-times and all-pairs SSD to a single CSV.

    sources: ordered dict of {name: (peaks_array, fs, offset)}.
    Each source's absolute peak time is offset + index/fs, used for both
    the printed _sec column and the cross-source SSDs.

    Columns: <name>_sec for each source, then ssd_<a>_<b> for every
    unordered pair. SSD values appear only on the first data row.
    """
    names = list(sources.keys())
    cols = [
        [float(idx) / fs + offset for idx in peaks]
        for (peaks, fs, offset) in sources.values()
    ]

    pair_names = []
    pair_values = []
    for a, b in combinations(names, 2):
        pa, fa, oa = sources[a]
        pb, fb, ob = sources[b]
        pair_names.append(f"ssd_{a}_{b}")
        pair_values.append(pairwise_ssd(pa, pb, fa, fb, oa, ob))

    def fmt_ssd(v):
        return f"{v:.6f}" if not (isinstance(v, float) and np.isnan(v)) else ""

    header = [f"{n}_sec" for n in names] + pair_names

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        first = True
        for row in zip_longest(*cols, fillvalue=""):
            peak_cells = [f"{v:.6f}" if isinstance(v, float) else v for v in row]
            ssd_cells = (
                [fmt_ssd(v) for v in pair_values]
                if first
                else ["" for _ in pair_values]
            )
            first = False
            w.writerow(peak_cells + ssd_cells)


# ============================================================================
# File discovery
# ============================================================================


def find_subject_files(subject_id):
    """Return ({rate: path} for Mira, daniel_mat_path, deep_csv_path)."""
    mira_paths = {}
    for rate in MIRA_RATES:
        m = list(BIN_DIR.glob(f"{subject_id}_{rate}_wave_markings.bin"))
        if m:
            mira_paths[rate] = m[0]
    mat_matches = list(MAT_DIR.glob(f"{subject_id}_wave_data.mat"))
    deep_csv = None
    if CSV_DIR.exists():
        candidates = sorted(CSV_DIR.glob(f"{subject_id}*{DEEP_FILE_SUFFIX}"))
        if candidates:
            deep_csv = candidates[0]
    return (
        mira_paths,
        mat_matches[0] if mat_matches else None,
        deep_csv,
    )


# ============================================================================
# Main
# ============================================================================


def main():
    args = sys.argv[1:]
    show_raw = False
    skip_daniel = False
    if "--raw" in args:
        show_raw = True
        args.remove("--raw")
    if "--skip-daniel" in args:
        skip_daniel = True
        args.remove("--skip-daniel")

    if len(args) < 2:
        print(
            "Usage: python compare_single.py <subject_id> <bin_number> "
            "[--raw] [--skip-daniel]"
        )
        print("  e.g. python compare_single.py 3010023 33")
        sys.exit(1)

    subject_id = args[0]
    bin_idx = int(args[1])

    mira_paths, mat_path, deep_csv = find_subject_files(subject_id)
    print(f"Subject:     {subject_id}")
    print(f"Bin index:   {bin_idx}")
    for rate in MIRA_RATES:
        print(f"Mira {rate}Hz:  {mira_paths.get(rate)}")
    print(f"Daniel mat:  {mat_path}")
    print(f"Deep CSV:    {deep_csv}")
    print(f"Show raw:    {show_raw}")

    if not mira_paths and mat_path is None:
        print("\nNo data files found for this subject. Exiting.")
        sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # --- Mira (per rate) ---
    # mira_per_rate[rate] = {peaks, signal, raw_pairs}
    # mira_data_by_rate[rate] is the full list of bins for that rate, used
    # to derive Deep's bin boundaries (we anchor on whichever rate we have).
    # Only bins 0..bin_idx are read, and only the target bin parses signal
    # arrays; earlier bins skip past them.
    mira_per_rate = {}
    mira_data_by_rate = {}
    for rate, path in mira_paths.items():
        print(f"  reading Mira {rate}Hz...", flush=True)
        data = read_wave_bin(str(path), max_bins=bin_idx + 1)
        print(f"    done ({len(data)} bins read)", flush=True)
        mira_data_by_rate[rate] = data
        if not (0 <= bin_idx < len(data)):
            print(
                f"  Bin {bin_idx} out of range for Mira {rate}Hz (max {len(data) - 1})"
            )
            continue
        bin_obj = data[bin_idx]
        print(f"  Mira {rate}Hz bin {bin_idx} peak counts:")
        print(
            f"    ch1: raw={len(bin_obj['ch1_raw_idx'])} "
            f"sq={len(bin_obj['ch1_sq_idx'])} "
            f"abs={len(bin_obj['ch1_abs_idx'])}"
        )
        print(
            f"    ch2: raw={len(bin_obj['ch2_raw_idx'])} "
            f"sq={len(bin_obj['ch2_sq_idx'])} "
            f"abs={len(bin_obj['ch2_abs_idx'])}"
        )
        print(
            f"    ch3: raw={len(bin_obj['ch3_raw_idx'])} "
            f"sq={len(bin_obj['ch3_sq_idx'])} "
            f"abs={len(bin_obj['ch3_abs_idx'])}"
        )
        print(
            f"    ppg: max={len(bin_obj['ppgMaxAmps'])} "
            f"min={len(bin_obj['ppgMinAmps'])}"
        )
        print(f"    signal len={len(bin_obj['ch1_raw_sig'])}")

        raw_pairs_full = bin_obj.get("ch1_raw_pairs")
        # Bin window start (absolute seconds from recording start). The
        # ANNEALER stamps each segment with its (start_sample, end_sample)
        # in signal_rate units -- that's the authoritative bin-window
        # boundary. We divide by the file's rate (encoded in its filename)
        # to get the absolute time. NOTE: indices in ecg_bin_indexs are
        # 1-based for MATLAB compatibility, so subtract 1 first.
        ecg_first_sample = bin_obj.get("ecg_bin_first_sample")
        if ecg_first_sample is not None:
            bin_start_t = (ecg_first_sample - 1) / float(rate)
        elif raw_pairs_full is not None and len(raw_pairs_full) > 0:
            # Fallback: the first raw (t, v) pair lands at most one
            # native-sample period inside the window, so this is approximate
            # but better than nothing.
            print(
                f"  WARNING: Mira {rate}Hz bin {bin_idx} missing ecg_bin_indexs; "
                f"falling back to first raw timestamp (alignment may drift)."
            )
            bin_start_t = float(raw_pairs_full[0, 0])
        else:
            print(
                f"  WARNING: Mira {rate}Hz bin {bin_idx} has no anchor; "
                f"can't align x-axis. Treating as 0."
            )
            bin_start_t = 0.0
        raw_pairs = raw_pairs_full if show_raw else None
        if show_raw and (raw_pairs is None or len(raw_pairs) == 0):
            print(
                f"  WARNING: --raw requested but Mira {rate}Hz bin {bin_idx} has no raw pairs"
            )
        mira_per_rate[rate] = {
            "peaks": np.asarray(bin_obj["ch1_raw_idx"], dtype=np.intp),
            "signal": np.asarray(bin_obj["ch1_raw_sig"], dtype=np.float64),
            "raw_pairs": raw_pairs,
            "bin_start_t": bin_start_t,
        }

    # --- Daniel ---
    print("  reading Daniel...", flush=True)
    daniel_peaks = np.array([], dtype=np.intp)
    daniel_signal = np.array([], dtype=np.float64)
    if skip_daniel:
        print("    skipped (--skip-daniel)", flush=True)
    elif mat_path is not None:
        daniel_peaks, _, daniel_signal = read_daniel_bin(str(mat_path), bin_idx)
        print(
            f"    done ({len(daniel_peaks)} peaks, signal len {len(daniel_signal)})",
            flush=True,
        )
    else:
        print("    no .mat file found", flush=True)

    # --- Deep ---
    # Deep needs an anchor rate to map wall-clock timestamps -> sample indices
    # and to find which bin a peak belongs to. Prefer 1000Hz, then 2000, then 500.
    print("  reading Deep...", flush=True)
    deep_peaks = np.array([], dtype=np.intp)
    deep_anchor_rate = None
    if deep_csv is not None:
        for r in DEEP_ANCHOR_RATE_PREFERENCE:
            if r in mira_data_by_rate and mira_data_by_rate[r]:
                deep_anchor_rate = r
                break
        if deep_anchor_rate is None:
            print(
                "  WARNING: Deep CSV present but no Mira data available; "
                "cannot determine bin boundaries."
            )
        else:
            deep_peaks = read_deep_bin(
                deep_csv,
                bin_idx,
                mira_data_by_rate[deep_anchor_rate],
                deep_anchor_rate,
            )
    print(f"    done ({len(deep_peaks)} peaks)", flush=True)

    # --- Build plot inputs ---
    # signals: (sig, color, label, fs, bin_offset) per drawn trace.
    # overlay_sets: (peaks, color, marker, label, sig_idx).
    panel_signals = []
    overlay_sets = []
    raw_pairs_for_plot = None

    # Cross-rate alignment (Fix B): each Mira rate's bin starts at a slightly
    # different absolute time because the annealer's bin boundaries are
    # quantized to 1/final_sampling_rate. We pick the EARLIEST bin start
    # across all loaded rates as plot t=0 and shift each rate's signal by
    # its (own_start - earliest) offset. After this shift the three Mira
    # waveforms overlay precisely.
    if mira_per_rate:
        t_origin = min(e["bin_start_t"] for e in mira_per_rate.values())
    else:
        t_origin = 0.0
    print(f"  alignment t_origin = {t_origin:.6f} s; per-rate offsets:", flush=True)
    for rate in MIRA_RATES:
        if rate in mira_per_rate:
            print(
                f"    Mira {rate}Hz: +{mira_per_rate[rate]['bin_start_t'] - t_origin:.6f} s",
                flush=True,
            )

    rate_to_sig_idx = {}
    for rate in MIRA_RATES:
        if rate not in mira_per_rate:
            continue
        sig_idx = len(panel_signals)
        rate_to_sig_idx[rate] = sig_idx
        entry = mira_per_rate[rate]
        bin_offset = entry["bin_start_t"] - t_origin
        panel_signals.append(
            (
                entry["signal"],
                MIRA_COLORS[rate],
                f"Mira {rate}Hz",
                float(rate),
                bin_offset,
            )
        )
        overlay_sets.append(
            (
                entry["peaks"],
                MIRA_COLORS[rate],
                MIRA_MARKERS[rate],
                f"Mira {rate}",
                sig_idx,
            )
        )
        # Use raw pairs from the first Mira rate that has them. They're in
        # absolute seconds-from-recording-start; pass t_origin so they land
        # on the same shared x axis as the signals.
        if (
            raw_pairs_for_plot is None
            and entry["raw_pairs"] is not None
            and len(entry["raw_pairs"]) > 0
        ):
            raw_pairs_for_plot = (entry["raw_pairs"], t_origin)

    # Daniel's wave_data bin doesn't carry an absolute start timestamp, so
    # we pin Daniel's bin start by aligning its FIRST R-peak with the first
    # R-peak of the Mira anchor rate (preferring 1000Hz). Both algorithms
    # almost always agree on the first peak in a 60s segment, so this gives
    # a sub-millisecond anchor. Any drift after the first peak reflects real
    # algorithmic disagreement, which is what we actually want SSD to capture.
    #
    # NOTE: there are TWO offset spaces in this script and they're easy to
    # confuse:
    #   - Absolute time (s from recording start): used by sources dict and
    #     pairwise SSD. CSV columns print absolute peak times.
    #   - Plot-relative time (s from t_origin): used by panel_signals so
    #     all traces sit in the [0, 60]-ish plot window. Plot offsets are
    #     just `absolute_offset - t_origin`.
    daniel_sig_idx = None
    daniel_abs_offset = t_origin  # default: anchor to plot t=0 if no Mira peaks
    if len(daniel_signal) > 0:
        anchor_for_daniel = None
        for r in DEEP_ANCHOR_RATE_PREFERENCE:  # same preference order as Deep
            if r in mira_per_rate and len(mira_per_rate[r]["peaks"]) > 0:
                anchor_for_daniel = r
                break
        if anchor_for_daniel is not None and len(daniel_peaks) > 0:
            mira_first_peak_abs = (
                int(mira_per_rate[anchor_for_daniel]["peaks"][0])
                / float(anchor_for_daniel)
                + mira_per_rate[anchor_for_daniel]["bin_start_t"]
            )
            daniel_first_peak_local = int(daniel_peaks[0]) / 1000.0
            daniel_abs_offset = mira_first_peak_abs - daniel_first_peak_local
            print(
                f"  Daniel anchor: first peak from Mira {anchor_for_daniel}Hz at "
                f"{mira_first_peak_abs:.6f} s, Daniel local at "
                f"{daniel_first_peak_local:.6f} s -> abs offset "
                f"{daniel_abs_offset:.6f} s (plot offset "
                f"{daniel_abs_offset - t_origin:.6f} s)",
                flush=True,
            )
        else:
            print(
                "  Daniel anchor: no Mira peaks available; "
                "Daniel will use t_origin (may drift)."
            )

        daniel_sig_idx = len(panel_signals)
        panel_signals.append(
            (
                daniel_signal,
                "tab:red",
                "Daniel signal",
                1000.0,
                daniel_abs_offset - t_origin,
            )
        )
        overlay_sets.append((daniel_peaks, "tab:red", "s", "Daniel", daniel_sig_idx))

    # Deep has no signal of its own; overlay it on whichever Mira rate we
    # used as the anchor for boundary derivation. Falls back to Daniel if
    # Mira is unavailable.
    if len(deep_peaks) > 0:
        anchor_idx = (
            rate_to_sig_idx.get(deep_anchor_rate)
            if deep_anchor_rate is not None
            else daniel_sig_idx
        )
        if anchor_idx is not None:
            overlay_sets.append((deep_peaks, "tab:green", "^", "Deep", anchor_idx))

    print("  plotting full window...", flush=True)
    plot_overlay(
        panel_signals,
        overlay_sets,
        T_START,
        T_END,
        title=f"{subject_id}  -  Bin {bin_idx}  ({T_START}-{T_END}s)",
        out_path=OUT_DIR / f"{subject_id}_bin{bin_idx:03d}_overlay.png",
        figsize=(14, 7),
        raw_pairs=raw_pairs_for_plot,
    )

    print("  plotting zoom...", flush=True)
    plot_overlay(
        panel_signals,
        overlay_sets,
        ZOOM_T_START,
        ZOOM_T_END,
        title=f"{subject_id}  -  Bin {bin_idx}  ({ZOOM_T_START}-{ZOOM_T_END}s)",
        out_path=OUT_DIR / f"{subject_id}_bin{bin_idx:03d}_overlay_zoom.png",
        figsize=(7, 7),
        raw_pairs=raw_pairs_for_plot,
    )

    # --- Results CSV: peak times for each Mira rate + Daniel + Deep,
    #     plus all-pairs SSD on the first row. Each source carries an
    #     absolute-time offset so peaks compare on a common time axis.
    print("  writing CSV...", flush=True)
    sources = {}
    for rate in MIRA_RATES:
        if rate in mira_per_rate:
            entry = mira_per_rate[rate]
            # Offset = bin_start_t (absolute s from recording start). Peaks
            # are bin-local indices into the upsampled signal.
            sources[f"mira{rate}"] = (entry["peaks"], float(rate), entry["bin_start_t"])
    if len(daniel_peaks) > 0:
        # Daniel's offset (absolute s from recording start) was computed
        # above by first-beat alignment. CSV columns and SSDs work in
        # absolute time, so we pass it as-is.
        sources["daniel"] = (daniel_peaks, 1000.0, daniel_abs_offset)
    if len(deep_peaks) > 0 and deep_anchor_rate is not None:
        # Deep's bin-local peaks were derived against the anchor rate's
        # mira_data, so they share that rate's bin_start_t.
        deep_offset = (
            mira_per_rate[deep_anchor_rate]["bin_start_t"]
            if deep_anchor_rate in mira_per_rate
            else t_origin
        )
        sources["deep"] = (deep_peaks, float(deep_anchor_rate), deep_offset)

    results_path = OUT_DIR / f"{subject_id}_bin{bin_idx:03d}_results.csv"
    write_results_csv(results_path, sources)

    print()
    counts = "  ".join(f"{name}={len(p)}" for name, (p, _, _) in sources.items())
    print(f"R-peak counts: {counts}")
    print(f"Outputs in: {OUT_DIR}")


if __name__ == "__main__":
    main()
