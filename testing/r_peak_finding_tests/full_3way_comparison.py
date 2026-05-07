"""
compare_wave_bounds.py
Compare wave bound outputs between C++ (_wave_markings.bin),
MATLAB (_wave_data.mat), and Deep's whole-case _wholecaseRRiQTi.csv.
"""

import csv
import math as _math
import struct
from pathlib import Path

import numpy as np
import scipy.io as sio

# ============================================================================
# Config
# ============================================================================

BIN_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_mira"
)
MAT_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_daniel"
)
CSV_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_deep"
)
OUTPUT_DIR = Path(r"D:\USERS\MiraWelner\QTVI\testing\4_wave_finding_tests\results")

SR = 1000.0
MAX_SANE = 50_000_000
PASSTHROUGH_NUM_CHANNELS = 41

# Deep's data is stored as one whole-case tab-separated CSV per subject,
# e.g. 3010023_20110817_part1_ECG_fs1000_wholecaseRRiQTi.csv. Columns include:
#     timestamp (i) ; seconds from 1st Midnight
# We read the absolute timestamp column, subtract the first timestamp so the
# series starts at zero, then convert to sample indices (multiply by SR).
# The peak series is then sliced per bin using Mira's bin signal lengths.
DEEP_FILE_SUFFIX = "_wholecaseRRiQTi.csv"

# Per-position penalty (in seconds) added in quadrature for every entry
# in the longer array that has no counterpart in the shorter one.
MISS_PENALTY_S = 1.0

# Maximum offset (in seconds) at which two R-peaks from different
# detectors are considered the same beat (matched-NN R-peak SSD only).
MATCH_TOLERANCE_S = 0.15


# ============================================================================
# Read C++ _wave_markings.bin
# ============================================================================


def read_bin_file(path):
    results = []
    with open(path, "rb") as f:

        def read_u64():
            return struct.unpack("<Q", f.read(8))[0]

        def read_idx():
            sz = read_u64()
            if sz > MAX_SANE or sz == 0:
                return np.array([], dtype=np.intp)
            return (np.frombuffer(f.read(sz * 8), dtype=np.uint64).copy() - 1).astype(
                np.intp
            )

        def skip_idx():
            sz = read_u64()
            if sz <= MAX_SANE:
                f.seek(sz * 8, 1)

        def read_signal():
            sz = read_u64()
            if sz > MAX_SANE or sz == 0:
                return np.array([], dtype=np.float64)
            return np.frombuffer(f.read(sz * 8), dtype=np.float64).copy()

        def skip_signal():
            sz = read_u64()
            if sz <= MAX_SANE:
                f.seek(sz * 8, 1)

        def skip_raw_pairs():
            n_pairs = read_u64()
            if n_pairs <= MAX_SANE:
                f.seek(n_pairs * 16, 1)

        def read_pair_vec():
            sz = read_u64()
            if sz > MAX_SANE or sz == 0:
                return np.empty((0, 2), dtype=np.uint64)
            return np.frombuffer(f.read(sz * 16), dtype=np.uint64).reshape(sz, 2).copy()

        num_bins = read_u64()
        if num_bins > 100_000:
            return results

        for _ in range(num_bins):
            b = {}
            b["ecgRIndex"] = read_idx()
            for _ in range(9):
                skip_idx()
            b["ppgMinAmps"] = read_idx()
            b["ppgSignal"] = read_signal()
            b["ecgSignal"] = read_signal()
            skip_signal()
            skip_signal()
            for _ in range(6):
                skip_signal()
            f.read(9)

            num_pairs = read_u64()
            if 0 < num_pairs < MAX_SANE:
                raw = (
                    np.frombuffer(f.read(num_pairs * 16), dtype=np.int64)
                    .reshape(num_pairs, 2)
                    .copy()
                )
                b["pairs"] = np.where(raw == -1, -1.0, raw - 1).astype(np.float64)
            else:
                b["pairs"] = np.empty((0, 2), dtype=np.float64)

            read_pair_vec()
            read_pair_vec()

            for _ in range(PASSTHROUGH_NUM_CHANNELS):
                skip_signal()
                skip_raw_pairs()

            results.append(b)

    return results


# ============================================================================
# Read MATLAB _wave_data.mat
# ============================================================================


def read_mat_file(path):
    results = []
    mat = sio.loadmat(str(path), squeeze_me=False)
    wave_data = mat["wave_data"]

    for cell in wave_data.flat:
        r = {
            "ecgRIndex": np.array([]),
            "ppgMinAmps": np.array([]),
            "ecgSignal": np.array([]),
            "ppgSignal": np.array([]),
        }

        if cell is None or (hasattr(cell, "size") and cell.size == 0):
            results.append(r)
            continue

        obj = cell
        while hasattr(obj, "shape") and obj.shape == (1, 1):
            obj = obj[0, 0]
        if not (hasattr(obj, "dtype") and obj.dtype.names):
            results.append(r)
            continue

        def get_field(name):
            if name not in obj.dtype.names:
                return np.array([], dtype=np.float64)
            val = obj[name]
            while hasattr(val, "shape") and val.ndim > 1 and 1 in val.shape:
                val = val.squeeze()
            return np.array(val, dtype=np.float64).flatten()

        r["ecgRIndex"] = get_field("ecgRIndex").astype(np.intp) - 1
        r["ecgRIndex_f"] = get_field("ecgRIndex") - 1.0
        r["ppgMinAmps"] = get_field("ppgMinAmps").astype(np.intp) - 1
        results.append(r)

    return results


# ============================================================================
# Read Deep whole-case CSV (_wholecaseRRiQTi.csv)
# ============================================================================
#
# Deep's pipeline emits one tab-separated file per subject covering the whole
# recording. We read the absolute timestamp column ("seconds from 1st
# Midnight"), subtract the first timestamp so the series starts at zero, then
# convert to sample indices by multiplying by SR. The resulting global peak
# array is sliced per bin using Mira's bin signal lengths.


def _parse_deep_timestamps_s(csv_path, debug=False):
    """Parse Deep's whole-case CSV and return absolute timestamps in seconds.

    The file is tab-separated. Column 1 (0-based) is
    "timestamp (i) ; seconds from 1st Midnight". The header row is
    skipped automatically because its cell fails to float-parse.

    Returns a numpy array of timestamps in seconds (not yet zeroed).
    """
    timestamps = []
    n_skipped = 0

    with open(csv_path, "r", newline="", encoding="utf-8-sig") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) < 2:
                n_skipped += 1
                continue
            try:
                timestamps.append(float(parts[1].strip()))
            except (ValueError, IndexError):
                # Header row or malformed cell — skip.
                n_skipped += 1
                continue

    if debug:
        print(
            f"\n  [deep parse] {csv_path}\n"
            f"    timestamps parsed: {len(timestamps)}, skipped: {n_skipped}"
        )

    return np.array(timestamps, dtype=np.float64)


def _derive_deep_global_peaks(csv_path, fs):
    """Return Deep's R-peak indices in global sample coordinates.

    Reads absolute timestamps (seconds from first midnight), converts to sample indices by multiplying by fs.
    Both Deep and Mira use seconds-from-midnight as their time reference,
    so no zeroing offset is needed.
    """
    if csv_path is None:
        return np.array([], dtype=np.intp)
    csv_path = Path(csv_path)
    if not csv_path.exists() or not csv_path.is_file():
        return np.array([], dtype=np.intp)

    timestamps_s = _parse_deep_timestamps_s(csv_path)
    if len(timestamps_s) == 0:
        return np.array([], dtype=np.intp)

    return np.round(timestamps_s * fs).astype(np.intp)


def read_deep_file(csv_path, mira_bins, n_bins_expected):
    """Return a per-bin list of Deep dicts, sliced from the whole-case CSV.

    csv_path:        Path to the *_wholecaseRRiQTi.csv (or None / missing).
    mira_bins:       Mira's bin list, used to compute each bin's start and
                     length in the global sample timeline.
    n_bins_expected: Number of bins to return. Trailing bins beyond Mira's
                     range are returned empty.

    Each output dict mirrors the schema from read_bin_file / read_mat_file.
    Only ecgRIndex / ecgRIndex_f are populated; signals are empty (Deep's
    file carries no waveform).
    """
    empty = {
        "ecgRIndex": np.array([], dtype=np.intp),
        "ecgRIndex_f": np.array([], dtype=np.float64),
        "ppgMinAmps": np.array([], dtype=np.intp),
        "ecgSignal": np.array([], dtype=np.float64),
        "ppgSignal": np.array([], dtype=np.float64),
    }

    if csv_path is None or not Path(csv_path).is_file() or not mira_bins:
        return [dict(empty) for _ in range(n_bins_expected)]

    global_peaks = _derive_deep_global_peaks(csv_path, SR)
    if len(global_peaks) == 0:
        return [dict(empty) for _ in range(n_bins_expected)]

    # Build bin start offsets from Mira's signal lengths.
    bin_offsets = []
    cum = 0
    for b in mira_bins:
        bin_offsets.append(cum)
        cum += len(b.get("ecgSignal", []))

    results = []
    for bin_idx in range(n_bins_expected):
        if bin_idx >= len(mira_bins):
            results.append(dict(empty))
            continue
        bin_start = bin_offsets[bin_idx]
        bin_length = len(mira_bins[bin_idx].get("ecgSignal", []))
        if bin_length <= 0:
            results.append(dict(empty))
            continue
        mask = (global_peaks >= bin_start) & (global_peaks < bin_start + bin_length)
        local = (global_peaks[mask] - bin_start).astype(np.intp)
        results.append(
            {
                "ecgRIndex": local,
                "ecgRIndex_f": local.astype(np.float64),
                "ppgMinAmps": np.array([], dtype=np.intp),
                "ecgSignal": np.array([], dtype=np.float64),
                "ppgSignal": np.array([], dtype=np.float64),
            }
        )

    return results


# ============================================================================
# Per-bin stats
# ============================================================================


def rr_stats(indices):
    idx = np.array(indices, dtype=np.float64)
    if len(idx) < 2:
        return dict(
            mean=np.nan,
            median=np.nan,
            q1=np.nan,
            q3=np.nan,
            rate_mean=np.nan,
            rate_median=np.nan,
            rate_q1=np.nan,
            rate_q3=np.nan,
        )
    rr_ms = np.diff(idx) / SR * 1000.0
    rate = 60000.0 / rr_ms
    return dict(
        mean=float(np.mean(rr_ms)),
        median=float(np.median(rr_ms)),
        q1=float(np.percentile(rr_ms, 25)),
        q3=float(np.percentile(rr_ms, 75)),
        rate_mean=float(np.mean(rate)),
        rate_median=float(np.median(rate)),
        rate_q1=float(np.percentile(rate, 25)),
        rate_q3=float(np.percentile(rate, 75)),
    )


def amp_stats(indices, signal):
    idx = np.array(indices, dtype=np.intp)
    idx = idx[(idx >= 0) & (idx < len(signal))]
    if len(idx) == 0 or len(signal) == 0:
        return dict(mean=np.nan, median=np.nan, q1=np.nan, q3=np.nan)
    amps = signal[idx]
    return dict(
        mean=float(np.mean(amps)),
        median=float(np.median(amps)),
        q1=float(np.percentile(amps, 25)),
        q3=float(np.percentile(amps, 75)),
    )


def _rr_array_ms(indices):
    idx = np.asarray(indices, dtype=np.float64)
    if len(idx) < 2:
        return np.array([], dtype=np.float64)
    return np.diff(idx) / SR * 1000.0


def _peak_amplitudes(indices, signal):
    sig = np.asarray(signal, dtype=np.float64)
    if len(sig) == 0 or len(indices) == 0:
        return np.array([], dtype=np.float64)
    ii = np.asarray(indices, dtype=np.intp)
    ii = ii[(ii >= 0) & (ii < len(sig))]
    if len(ii) == 0:
        return np.array([], dtype=np.float64)
    return sig[ii]


def rpeak_ssd(
    a, b, fs=SR, tolerance_s=MATCH_TOLERANCE_S, miss_penalty_s=MISS_PENALTY_S
):
    """R-peak SSD (s^2) using nearest-neighbor matching within tolerance."""
    a_arr = np.sort(np.asarray(a, dtype=np.float64)) / fs
    b_arr = np.sort(np.asarray(b, dtype=np.float64)) / fs
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
        left = b_arr[pos - 1]
        right = b_arr[pos]
        pick_left = (a_arr - left) <= (right - a_arr)
        nearest_idx = np.where(pick_left, pos - 1, pos)
    else:
        nearest_idx = np.zeros_like(pos)

    nearest_val = b_arr[nearest_idx]
    diff = a_arr - nearest_val
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
    n_a_unmatched = na - n_matched
    n_b_unmatched = nb - n_matched

    ssd = float(np.sum(matched_diffs**2))
    ssd += (n_a_unmatched + n_b_unmatched) * (miss_penalty_s**2)
    return ssd


def rr_ssd(a, b, fs=SR):
    a_arr = np.asarray(a, dtype=np.float64)
    b_arr = np.asarray(b, dtype=np.float64)

    rr_a = np.diff(a_arr) / fs if len(a_arr) >= 2 else np.array([], dtype=np.float64)
    rr_b = np.diff(b_arr) / fs if len(b_arr) >= 2 else np.array([], dtype=np.float64)

    n = min(len(rr_a), len(rr_b))
    if n == 0:
        return float("nan")

    diffs = rr_a[:n] - rr_b[:n]
    return float(np.sum(diffs**2))


def compare_file(file_id, bin_data, mat_data):
    total = min(len(bin_data), len(mat_data))
    bins = []
    for i in range(total):
        b = bin_data[i]
        m = mat_data[i]

        b_rr = rr_stats(b["ecgRIndex"])
        m_rr = rr_stats(m["ecgRIndex"])
        b_ppg = rr_stats(b["ppgMinAmps"])
        m_ppg = rr_stats(m["ppgMinAmps"])

        b_ecg_amp = amp_stats(b["ecgRIndex"], b["ecgSignal"])
        m_ecg_amp = amp_stats(m["ecgRIndex"], m["ecgSignal"])
        b_ppg_amp = amp_stats(b["ppgMinAmps"], b["ppgSignal"])
        m_ppg_amp = amp_stats(m["ppgMinAmps"], m["ppgSignal"])

        b_rr_arr_ms = _rr_array_ms(b["ecgRIndex"])
        m_rr_arr_ms = _rr_array_ms(m["ecgRIndex"])
        b_amp_arr = _peak_amplitudes(b["ecgRIndex"], b.get("ecgSignal", []))
        m_amp_arr = _peak_amplitudes(m["ecgRIndex"], m.get("ecgSignal", []))

        b_rr_raw = np.diff(np.array(b["ecgRIndex"], dtype=np.float64))
        m_rr_raw = np.diff(np.array(m["ecgRIndex"], dtype=np.float64))
        if len(b_rr_raw) == len(m_rr_raw) and len(b_rr_raw) > 0:
            rr_identical = np.array_equal(b_rr_raw, m_rr_raw)
        else:
            rr_identical = False

        b_rloc = np.array(b["ecgRIndex"], dtype=np.float64)
        m_rloc_f = np.array(m.get("ecgRIndex_f", m["ecgRIndex"]), dtype=np.float64)
        rloc_ssd_val = rpeak_ssd(b_rloc, m_rloc_f, SR)
        bin_rr_ssd_val = rr_ssd(b["ecgRIndex"], m["ecgRIndex"], SR)

        if len(b_rloc) == len(m_rloc_f) and len(b_rloc) > 0:
            signed_f = b_rloc - m_rloc_f
            if np.all(np.abs(signed_f - signed_f[0]) < 1e-6):
                bin_offset_s = float(signed_f[0]) / SR
            else:
                bin_offset_s = -0.1
        else:
            bin_offset_s = float("nan")

        bins.append(
            {
                "bin": i,
                "n_ecg_cpp": len(b["ecgRIndex"]),
                "n_ecg_mat": len(m["ecgRIndex"]),
                "n_ppg_cpp": len(b["ppgMinAmps"]),
                "n_ppg_mat": len(m["ppgMinAmps"]),
                "rr_cpp": b_rr,
                "rr_mat": m_rr,
                "ppg_rr_cpp": b_ppg,
                "ppg_rr_mat": m_ppg,
                "ecg_amp_cpp": b_ecg_amp,
                "ecg_amp_mat": m_ecg_amp,
                "ppg_amp_cpp": b_ppg_amp,
                "ppg_amp_mat": m_ppg_amp,
                "rr_cpp_arr_ms": b_rr_arr_ms,
                "rr_mat_arr_ms": m_rr_arr_ms,
                "amp_cpp_arr": b_amp_arr,
                "amp_mat_arr": m_amp_arr,
                "rr_ssd": bin_rr_ssd_val,
                "rr_identical": rr_identical,
                "rloc_ssd": rloc_ssd_val,
                "bin_offset_s": bin_offset_s,
            }
        )
    return {"id": file_id, "total_bins": total, "bins": bins}


# ============================================================================
# Summary CSV
# ============================================================================


def nanmean(vals):
    v = [x for x in vals if not np.isnan(x)]
    return float(np.mean(v)) if v else np.nan


def nanmedian(vals):
    v = [x for x in vals if not np.isnan(x)]
    return float(np.median(v)) if v else np.nan


def nanpercentile(vals, p):
    v = [x for x in vals if not np.isnan(x)]
    return float(np.percentile(v, p)) if v else np.nan


def fmt(v, d=6):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    if isinstance(v, str):
        return v
    return round(v, d)


def iqr_str(q1, q3, d=6):
    if np.isnan(q1) or np.isnan(q3):
        return ""
    return f"{round(q1, d)},{round(q3, d)}"


def _pool_mean(arr):
    return float(np.mean(arr)) if len(arr) else float("nan")


def _pool_median(arr):
    return float(np.median(arr)) if len(arr) else float("nan")


def _pool_pct(arr, p):
    return float(np.percentile(arr, p)) if len(arr) else float("nan")


def write_summary_csv(path, results, label_a="A", label_b="B"):
    all_bins = [b for fr in results for b in fr["bins"]]
    n_subjects = len(results)

    rr_cpp_pool = (
        np.concatenate(
            [b.get("rr_cpp_arr_ms", np.array([], dtype=np.float64)) for b in all_bins]
        )
        if all_bins
        else np.array([], dtype=np.float64)
    )
    rr_mat_pool = (
        np.concatenate(
            [b.get("rr_mat_arr_ms", np.array([], dtype=np.float64)) for b in all_bins]
        )
        if all_bins
        else np.array([], dtype=np.float64)
    )
    rate_cpp_pool = (
        (60000.0 / rr_cpp_pool[rr_cpp_pool > 0]) if len(rr_cpp_pool) else np.array([])
    )
    rate_mat_pool = (
        (60000.0 / rr_mat_pool[rr_mat_pool > 0]) if len(rr_mat_pool) else np.array([])
    )

    amp_cpp_pool = (
        np.concatenate(
            [b.get("amp_cpp_arr", np.array([], dtype=np.float64)) for b in all_bins]
        )
        if all_bins
        else np.array([], dtype=np.float64)
    )
    amp_mat_pool = (
        np.concatenate(
            [b.get("amp_mat_arr", np.array([], dtype=np.float64)) for b in all_bins]
        )
        if all_bins
        else np.array([], dtype=np.float64)
    )

    extra_a = sum(max(b["n_ecg_cpp"] - b["n_ecg_mat"], 0) for b in all_bins)
    extra_b = sum(max(b["n_ecg_mat"] - b["n_ecg_cpp"], 0) for b in all_bins)

    cols = [label_a, label_b]

    def make_row(label, vals):
        return [label] + [fmt(v) for v in vals]

    def make_iqr_row(label, vals_q1, vals_q3):
        return [label] + [iqr_str(q1, q3) for q1, q3 in zip(vals_q1, vals_q3)]

    rows = [
        ["Metric"] + cols,
        make_row("N Subjects", [n_subjects, n_subjects]),
        make_row(
            "N Beats Total",
            [
                sum(b["n_ecg_cpp"] for b in all_bins),
                sum(b["n_ecg_mat"] for b in all_bins),
            ],
        ),
        make_row(f"Beats in {label_a} not in {label_b}", [extra_a, ""]),
        make_row(f"Beats in {label_b} not in {label_a}", ["", extra_b]),
        make_row(
            "Mean N Beats Per Patient",
            [
                nanmean([sum(b["n_ecg_cpp"] for b in fr["bins"]) for fr in results]),
                nanmean([sum(b["n_ecg_mat"] for b in fr["bins"]) for fr in results]),
            ],
        ),
        make_row(
            "Mean Rate (bpm)", [_pool_mean(rate_cpp_pool), _pool_mean(rate_mat_pool)]
        ),
        make_row(
            "Median Rate (bpm)",
            [_pool_median(rate_cpp_pool), _pool_median(rate_mat_pool)],
        ),
        make_iqr_row(
            "Rate IQR (bpm)",
            [_pool_pct(rate_cpp_pool, 25), _pool_pct(rate_mat_pool, 25)],
            [_pool_pct(rate_cpp_pool, 75), _pool_pct(rate_mat_pool, 75)],
        ),
        make_row(
            "Mean Interval (ms)", [_pool_mean(rr_cpp_pool), _pool_mean(rr_mat_pool)]
        ),
        make_row(
            "Median Interval (ms)",
            [_pool_median(rr_cpp_pool), _pool_median(rr_mat_pool)],
        ),
        make_iqr_row(
            "Interval IQR (ms)",
            [_pool_pct(rr_cpp_pool, 25), _pool_pct(rr_mat_pool, 25)],
            [_pool_pct(rr_cpp_pool, 75), _pool_pct(rr_mat_pool, 75)],
        ),
        make_row(
            "Mean Amplitude", [_pool_mean(amp_cpp_pool), _pool_mean(amp_mat_pool)]
        ),
        make_row(
            "Median Amplitude", [_pool_median(amp_cpp_pool), _pool_median(amp_mat_pool)]
        ),
        make_iqr_row(
            "Amplitude IQR",
            [_pool_pct(amp_cpp_pool, 25), _pool_pct(amp_mat_pool, 25)],
            [_pool_pct(amp_cpp_pool, 75), _pool_pct(amp_mat_pool, 75)],
        ),
    ]

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        for row in rows:
            w.writerow(row)


# ============================================================================
# Three-way summary
# ============================================================================


def summarize_paired_source(pair_file_results, side):
    rr_key = f"rr_{side}_arr_ms"
    amp_key = f"amp_{side}_arr"
    n_key = f"n_ecg_{side}"

    n_subjects = len(pair_file_results)
    n_beats_tot = 0
    per_pat = []
    rr_chunks = []
    amp_chunks = []

    for fr in pair_file_results:
        pat_total = 0
        for b in fr["bins"]:
            pat_total += b[n_key]
            rr_chunks.append(b.get(rr_key, np.array([], dtype=np.float64)))
            amp_chunks.append(b.get(amp_key, np.array([], dtype=np.float64)))
        n_beats_tot += pat_total
        per_pat.append(pat_total)

    rr_pool = np.concatenate(rr_chunks) if rr_chunks else np.array([], dtype=np.float64)
    amp_pool = (
        np.concatenate(amp_chunks) if amp_chunks else np.array([], dtype=np.float64)
    )
    rate_pool = (60000.0 / rr_pool[rr_pool > 0]) if len(rr_pool) else rr_pool

    return {
        "n_subjects": n_subjects,
        "n_beats_total": n_beats_tot,
        "mean_n_per_pat": nanmean(per_pat),
        "rate_mean": _pool_mean(rate_pool),
        "rate_median": _pool_median(rate_pool),
        "rate_iqr_q1": _pool_pct(rate_pool, 25),
        "rate_iqr_q3": _pool_pct(rate_pool, 75),
        "rr_mean": _pool_mean(rr_pool),
        "rr_median": _pool_median(rr_pool),
        "rr_iqr_q1": _pool_pct(rr_pool, 25),
        "rr_iqr_q3": _pool_pct(rr_pool, 75),
        "amp_mean": _pool_mean(amp_pool),
        "amp_median": _pool_median(amp_pool),
        "amp_iqr_q1": _pool_pct(amp_pool, 25),
        "amp_iqr_q3": _pool_pct(amp_pool, 75),
    }


def write_three_way_summary_csv(path, summaries, pair_diffs=None):
    labels = list(summaries.keys())
    pair_diffs = pair_diffs or {}

    def row(metric, key):
        return [metric] + [fmt(summaries[lab][key]) for lab in labels]

    def iqr_row(metric, q1, q3):
        return [metric] + [
            iqr_str(summaries[lab][q1], summaries[lab][q3]) for lab in labels
        ]

    def diff_row(src_label):
        cells = []
        for col in labels:
            if col == src_label:
                cells.append("")
            else:
                v = pair_diffs.get((src_label, col), "")
                cells.append(fmt(v) if v != "" else "")
        return [f"Beats Diff From {src_label}"] + cells

    rows = [
        ["Metric"] + labels,
        row("N Subjects", "n_subjects"),
        row("N Beats Total", "n_beats_total"),
        row("Mean N Beats Per Patient", "mean_n_per_pat"),
    ]
    for lab in labels:
        rows.append(diff_row(lab))
    rows.extend(
        [
            row("Mean Rate (bpm)", "rate_mean"),
            row("Median Rate (bpm)", "rate_median"),
            iqr_row("Rate IQR (bpm)", "rate_iqr_q1", "rate_iqr_q3"),
            row("Mean Interval (ms)", "rr_mean"),
            row("Median Interval (ms)", "rr_median"),
            iqr_row("Interval IQR (ms)", "rr_iqr_q1", "rr_iqr_q3"),
            row("Mean Amplitude", "amp_mean"),
            row("Median Amplitude", "amp_median"),
            iqr_row("Amplitude IQR", "amp_iqr_q1", "amp_iqr_q3"),
        ]
    )

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        for r in rows:
            w.writerow(r)


# ============================================================================
# File summary CSV
# ============================================================================


def write_file_summary_csv(path, results):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "MESA ID",
                "MATLAB R Peak N",
                "C++ R Peak N",
                "Extra ECG MATLAB",
                "Extra ECG C++",
                "MATLAB PPG Peak N",
                "C++ PPG Peak N",
                "Extra PPG MATLAB",
                "Extra PPG C++",
            ]
        )
        for fr in results:
            mat_ecg = sum(b["n_ecg_mat"] for b in fr["bins"])
            cpp_ecg = sum(b["n_ecg_cpp"] for b in fr["bins"])
            mat_ppg = sum(b["n_ppg_mat"] for b in fr["bins"])
            cpp_ppg = sum(b["n_ppg_cpp"] for b in fr["bins"])
            extra_ecg_mat = sum(
                max(b["n_ecg_mat"] - b["n_ecg_cpp"], 0) for b in fr["bins"]
            )
            extra_ecg_cpp = sum(
                max(b["n_ecg_cpp"] - b["n_ecg_mat"], 0) for b in fr["bins"]
            )
            extra_ppg_mat = sum(
                max(b["n_ppg_mat"] - b["n_ppg_cpp"], 0) for b in fr["bins"]
            )
            extra_ppg_cpp = sum(
                max(b["n_ppg_cpp"] - b["n_ppg_mat"], 0) for b in fr["bins"]
            )
            w.writerow(
                [
                    fr["id"],
                    mat_ecg,
                    cpp_ecg,
                    extra_ecg_mat,
                    extra_ecg_cpp,
                    mat_ppg,
                    cpp_ppg,
                    extra_ppg_mat,
                    extra_ppg_cpp,
                ]
            )


# ============================================================================
# Per-file CSV
# ============================================================================


def write_per_file_csv(path, fr, label_a="C++", label_b="MATLAB"):
    write_summary_csv(path, [fr], label_a, label_b)

    bins = fr["bins"]
    summary_path = Path(path)
    per_bin_path = summary_path.with_name(
        summary_path.stem + "_per_bin" + summary_path.suffix
    )

    def col_mean(key_fn):
        vals = [key_fn(b) for b in bins if not np.isnan(key_fn(b))]
        return nanmean(vals)

    def col_sum(key_fn):
        vals = [key_fn(b) for b in bins if not np.isnan(key_fn(b))]
        return float(np.sum(vals)) if vals else np.nan

    sum_row = {
        "rpeaks_cpp": col_sum(lambda b: b["n_ecg_cpp"]),
        "rpeaks_mat": col_sum(lambda b: b["n_ecg_mat"]),
        "mean_rr_cpp": col_mean(lambda b: b["rr_cpp"]["mean"]),
        "mean_rr_mat": col_mean(lambda b: b["rr_mat"]["mean"]),
        "med_rr_cpp": col_mean(lambda b: b["rr_cpp"]["median"]),
        "med_rr_mat": col_mean(lambda b: b["rr_mat"]["median"]),
        "q1_cpp": nanpercentile([b["rr_cpp"]["q1"] for b in bins], 25),
        "q3_cpp": nanpercentile([b["rr_cpp"]["q3"] for b in bins], 75),
        "q1_mat": nanpercentile([b["rr_mat"]["q1"] for b in bins], 25),
        "q3_mat": nanpercentile([b["rr_mat"]["q3"] for b in bins], 75),
        "rr_ssd": col_sum(lambda b: b["rr_ssd"]),
        "rloc_ssd": col_sum(lambda b: b.get("rloc_ssd", float("nan"))),
    }
    mean_row = {
        "rpeaks_cpp": col_mean(lambda b: b["n_ecg_cpp"]),
        "rpeaks_mat": col_mean(lambda b: b["n_ecg_mat"]),
        "mean_rr_cpp": col_mean(lambda b: b["rr_cpp"]["mean"]),
        "mean_rr_mat": col_mean(lambda b: b["rr_mat"]["mean"]),
        "med_rr_cpp": col_mean(lambda b: b["rr_cpp"]["median"]),
        "med_rr_mat": col_mean(lambda b: b["rr_mat"]["median"]),
        "q1_cpp": sum_row["q1_cpp"],
        "q3_cpp": sum_row["q3_cpp"],
        "q1_mat": sum_row["q1_mat"],
        "q3_mat": sum_row["q3_mat"],
        "rr_ssd": col_mean(lambda b: b["rr_ssd"]),
        "rloc_ssd": col_mean(lambda b: b.get("rloc_ssd", float("nan"))),
    }

    with open(per_bin_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "bin",
                f"R Peaks {label_a}",
                f"R Peaks {label_b}",
                f"Mean RR {label_a} (ms)",
                f"Mean RR {label_b} (ms)",
                f"Median RR {label_a} (ms)",
                f"Median RR {label_b} (ms)",
                f"IQR RR {label_a} (q1,q3)",
                f"IQR RR {label_b} (q1,q3)",
                f"RR SSD {label_a} vs {label_b} (s^2)",
                f"R-Loc SSD {label_a} vs {label_b} (s^2)",
            ]
        )

        for label, d in [("SUM", sum_row), ("MEAN", mean_row)]:
            w.writerow(
                [
                    label,
                    fmt(d["rpeaks_cpp"]),
                    fmt(d["rpeaks_mat"]),
                    fmt(d["mean_rr_cpp"]),
                    fmt(d["mean_rr_mat"]),
                    fmt(d["med_rr_cpp"]),
                    fmt(d["med_rr_mat"]),
                    iqr_str(d["q1_cpp"], d["q3_cpp"]),
                    iqr_str(d["q1_mat"], d["q3_mat"]),
                    fmt(d["rr_ssd"]),
                    fmt(d["rloc_ssd"]),
                ]
            )

        for b in bins:
            w.writerow(
                [
                    b["bin"],
                    b["n_ecg_cpp"],
                    b["n_ecg_mat"],
                    fmt(b["rr_cpp"]["mean"]),
                    fmt(b["rr_mat"]["mean"]),
                    fmt(b["rr_cpp"]["median"]),
                    fmt(b["rr_mat"]["median"]),
                    iqr_str(b["rr_cpp"]["q1"], b["rr_cpp"]["q3"]),
                    iqr_str(b["rr_mat"]["q1"], b["rr_mat"]["q3"]),
                    fmt(b["rr_ssd"]),
                    fmt(b.get("rloc_ssd", float("nan"))),
                ]
            )


# ============================================================================
# SVG helpers
# ============================================================================


def _make_log_y_axis(lines, cx, cy, chart_h, max_count):
    log_max = _math.log10(max_count + 1) * 1.3 if max_count > 0 else 1.0

    def bar_h(count):
        if count <= 0:
            return 0.0
        return chart_h * _math.log10(count + 1) / log_max

    tick_vals = [0]
    v = 1
    while v <= max_count + 1:
        tick_vals.append(v)
        v *= 10
    if max_count not in tick_vals:
        tick_vals.append(max_count)

    for tv in tick_vals:
        yp = cy + chart_h - bar_h(tv)
        if yp < cy - 2:
            continue
        is_max = tv == max_count and max_count not in [0, 1, 10, 100, 1000, 10000]
        lines.append(
            f'<line x1="{cx - 3}" y1="{yp:.1f}" x2="{cx + 1}" y2="{yp:.1f}" '
            f'stroke="#333" stroke-width="1"/>'
        )
        lines.append(
            f'<text x="{cx - 5}" y="{yp + 4:.1f}" '
            f'text-anchor="end" class="tick"'
            f"{' fill="#e74c3c" font-weight="bold"' if is_max else ''}>"
            f"{tv}</text>"
        )
        if tv > 0:
            lines.append(
                f'<line x1="{cx + 1}" y1="{yp:.1f}" x2="{cx + 9999}" y2="{yp:.1f}" '
                f'stroke="#eee" stroke-width="0.5" stroke-dasharray="2,2"/>'
            )
    return bar_h


def _append_zero_ref_panel(
    lines,
    ox,
    oy,
    cell_w,
    cell_h,
    pad_l,
    pad_r,
    pad_t,
    pad_b,
    n_total,
    ref_max_val,
    x_label,
    num_buckets=50,
    fill_color="#2ecc71",
):
    chart_w = cell_w - pad_l - pad_r
    chart_h = cell_h - pad_t - pad_b
    cx = ox + pad_l
    cy = oy + pad_t

    lines.append(
        f'<text x="{ox + cell_w // 2}" y="{oy + 14}" '
        f'text-anchor="middle" class="subtitle" fill="#27ae60">'
        f"ZERO-DIFF REFERENCE (n={n_total})</text>"
    )
    lines.append(
        f'<text x="{ox + cell_w // 2}" y="{oy + 26}" '
        f'text-anchor="middle" class="stats">'
        f"all diffs=0 | x-range=[0, {ref_max_val:.6f}]</text>"
    )

    lines.append(
        f'<line x1="{cx}" y1="{cy}" x2="{cx}" y2="{cy + chart_h}" '
        f'stroke="#333" stroke-width="1"/>'
    )
    lines.append(
        f'<line x1="{cx}" y1="{cy + chart_h}" '
        f'x2="{cx + chart_w}" y2="{cy + chart_h}" '
        f'stroke="#333" stroke-width="1"/>'
    )
    lines.append(
        f'<text x="{ox + 18}" y="{cy + chart_h // 2}" '
        f'text-anchor="middle" class="y-axis-label" '
        f'transform="rotate(-90 {ox + 10},{cy + chart_h // 2})">'
        f"# bins</text>"
    )

    y_axis_max = n_total if n_total > 0 else 1
    for t in range(4):
        y_val = int(round(y_axis_max * t / 3.0))
        y_pos = cy + chart_h - (y_val / y_axis_max) * chart_h
        lines.append(
            f'<text x="{cx - 3}" y="{y_pos + 3}" '
            f'text-anchor="end" class="tick">{y_val}</text>'
        )
        if t > 0:
            lines.append(
                f'<line x1="{cx + 1}" y1="{y_pos}" '
                f'x2="{cx + chart_w}" y2="{y_pos}" '
                f'stroke="#eee" stroke-width="0.5"/>'
            )

    bar_w = chart_w / num_buckets
    bx = cx
    bar_h = chart_h
    by = cy
    lines.append(
        f'<rect x="{bx + 0.5:.1f}" y="{by:.1f}" '
        f'width="{bar_w - 1:.1f}" height="{bar_h:.1f}" '
        f'fill="{fill_color}" fill-opacity="0.8">'
        f"<title>[0, first_bucket): {n_total}</title></rect>"
    )
    lines.append(
        f'<text x="{bx + bar_w / 2:.1f}" y="{cy - 2}" '
        f'text-anchor="middle" class="clipped" fill="{fill_color}">{n_total}</text>'
    )

    for t in range(3):
        val = ref_max_val * t / 2.0
        x_pos = cx + chart_w * t / 2.0
        lines.append(
            f'<text x="{x_pos:.1f}" y="{cy + chart_h + 16}" '
            f'text-anchor="middle" class="x-tick">{val:.6f}</text>'
        )

    lines.append(
        f'<text x="{cx + chart_w // 2}" y="{cy + chart_h + 30}" '
        f'text-anchor="middle" class="axis-label">{x_label}</text>'
    )


# ============================================================================
# Combined single-chart histogram
# ============================================================================


def write_combined_histogram(path, results, field, chart_title, num_buckets=80):
    """One big histogram pooling `field` across every bin of every subject."""
    all_vals = []
    for fr in results:
        for b in fr["bins"]:
            v = b.get(field, float("nan"))
            if not np.isnan(v):
                all_vals.append(round(v, 6))

    if not all_vals:
        return

    arr = np.array(all_vals)
    mn = float(np.mean(arr))
    md = float(np.median(arr))
    q1 = float(np.percentile(arr, 25))
    q3 = float(np.percentile(arr, 75))
    zero_count = sum(1 for v in all_vals if v == 0.0)
    zero_pct = zero_count / len(all_vals) * 100.0 if all_vals else 0.0

    svg_w, svg_h = 900, 500
    pad_l, pad_r, pad_t, pad_b = 70, 40, 70, 60
    chart_w = svg_w - pad_l - pad_r
    chart_h = svg_h - pad_t - pad_b

    nonzero = arr[arr > 0]
    max_val = float(nonzero.max()) if len(nonzero) > 0 else 1.0
    min_val = 0.0
    n_nonzero_bins = num_buckets - 1
    bucket_width = max_val / n_nonzero_bins if n_nonzero_bins > 0 else 1.0

    counts = np.zeros(num_buckets, dtype=int)
    for v in all_vals:
        if v == 0.0:
            counts[0] += 1
        else:
            bi = 1 + int((v - 1e-300) / bucket_width)
            if bi >= num_buckets:
                bi = num_buckets - 1
            counts[bi] += 1

    max_count = int(counts.max()) if counts.max() > 0 else 1
    bar_w = chart_w / num_buckets
    cx, cy = pad_l, pad_t

    lines = []
    lines.append(
        f'<svg width="{svg_w}" height="{svg_h}" xmlns="http://www.w3.org/2000/svg">'
    )
    lines.append(
        "<style>\n"
        "  text { font-family: Consolas, 'Courier New', monospace; }\n"
        "  .title { font-size: 16px; font-weight: bold; }\n"
        "  .stats { font-size: 12px; fill: #555; }\n"
        "  .tick { font-size: 11px; }\n"
        "  .axis-label { font-size: 13px; }\n"
        "  .y-axis-label { font-size: 12px; }\n"
        "</style>"
    )
    lines.append('<rect width="100%" height="100%" fill="#fcfcfc"/>')

    lines.append(
        f'<text x="{svg_w // 2}" y="24" text-anchor="middle" class="title">'
        f"{chart_title} - All {len(results)} Subjects "
        f"({len(all_vals)} bins, {zero_pct:.1f}% 0.000000 err)</text>"
    )
    lines.append(
        f'<text x="{svg_w // 2}" y="44" text-anchor="middle" class="stats">'
        f"mean={mn:.6f}  median={md:.6f}  IQR=[{q1:.6f}, {q3:.6f}]</text>"
    )

    lines.append(
        f'<line x1="{cx}" y1="{cy}" x2="{cx}" y2="{cy + chart_h}" '
        f'stroke="#333" stroke-width="1"/>'
    )
    lines.append(
        f'<line x1="{cx}" y1="{cy + chart_h}" '
        f'x2="{cx + chart_w}" y2="{cy + chart_h}" '
        f'stroke="#333" stroke-width="1"/>'
    )
    lines.append(
        f'<text x="18" y="{cy + chart_h // 2}" '
        f'text-anchor="middle" class="y-axis-label" '
        f'transform="rotate(-90 18,{cy + chart_h // 2})">'
        f"# bins (log)</text>"
    )

    bar_h_fn = _make_log_y_axis(lines, cx, cy, chart_h, max_count)

    for b_idx in range(num_buckets):
        c = counts[b_idx]
        if c == 0:
            continue
        bh = bar_h_fn(c)
        bx = cx + b_idx * bar_w
        by = cy + chart_h - bh
        rs = min_val + b_idx * bucket_width
        re = min_val + (b_idx + 1) * bucket_width
        tip = f"[{rs:.6f}, {re:.6f}): {c}"
        lines.append(
            f'<rect x="{bx + 0.5:.1f}" y="{by:.1f}" '
            f'width="{bar_w - 1:.1f}" height="{bh:.1f}" '
            f'fill="#3498db" fill-opacity="0.85">'
            f"<title>{tip}</title></rect>"
        )

        n_x_ticks = 5
        for t in range(n_x_ticks + 1):
            val = max_val * t / n_x_ticks
            x_pos = cx + chart_w * t / n_x_ticks
        lines.append(
            f'<text x="{x_pos:.1f}" y="{cy + chart_h + 16}" '
            f'text-anchor="middle" class="tick">{val:.6f}</text>'
        )
    lines.append(
        f'<text x="{cx + chart_w // 2}" y="{cy + chart_h + 38}" '
        f'text-anchor="middle" class="axis-label">SSD (s^2)</text>'
    )
    lines.append("</svg>")

    with open(path, "w", encoding="utf-8") as fout:
        fout.write("\n".join(lines))


# ============================================================================
# Per-subject histograms (unified for both RR and R-loc)
# ============================================================================


def write_per_subject_histograms(
    path,
    results,
    field,
    chart_title,
    bar_color="#1a5fa8",
    ref_color="#2ecc71",
    num_buckets=50,
):
    """Per-subject grid (5 columns) of `field` (e.g. 'rr_ssd' or 'rloc_ssd'),
    plus a zero-difference reference panel in the next free slot.
    Bin 0 collects only exact-zero values; remaining bins divide (0, max]."""
    if not results:
        return

    num_files = len(results)
    n_cols = 5
    n_rows = (num_files + 1 + n_cols - 1) // n_cols

    cell_w, cell_h = 260, 275
    pad_l, pad_r, pad_t, pad_b = 50, 25, 48, 42
    chart_w = cell_w - pad_l - pad_r
    chart_h = cell_h - pad_t - pad_b

    svg_w = n_cols * cell_w
    svg_h = n_rows * cell_h + 55

    lines = []
    lines.append(
        f'<svg width="{svg_w}" height="{svg_h}" xmlns="http://www.w3.org/2000/svg">'
    )
    lines.append(
        "<style>\n"
        "  text { font-family: Consolas, 'Courier New', monospace; }\n"
        "  .title { font-size: 16px; font-weight: bold; }\n"
        "  .subtitle { font-size: 12px; font-weight: bold; }\n"
        "  .stats { font-size: 12px; fill: #555; }\n"
        "  .tick { font-size: 12px; }\n"
        "  .axis-label { font-size: 12px; }\n"
        "  .y-axis-label { font-size: 12px; }\n"
        "  .x-tick { font-size: 12px; }\n"
        "  .clipped { font-size: 12px; fill: #e74c3c; font-weight: bold; }\n"
        "  .empty { font-size: 12px; fill: #999; }\n"
        "</style>"
    )
    lines.append('<rect width="100%" height="100%" fill="#fcfcfc"/>')
    lines.append(
        f'<text x="{svg_w // 2}" y="30" text-anchor="middle" class="title">'
        f"{chart_title} - all {num_files} files</text>"
    )

    all_max_vals = []
    all_n = []

    for f_idx, fr in enumerate(results):
        col = f_idx % n_cols
        row = f_idx // n_cols
        ox = col * cell_w
        oy = row * cell_h + 50

        file_id = fr["id"]
        all_bins = fr["bins"]

        vals = [
            round(b.get(field, float("nan")), 6)
            for b in all_bins
            if not np.isnan(b.get(field, float("nan")))
        ]

        n_total = len(vals)
        zero_count = sum(1 for v in vals if v == 0.0)
        zero_pct = (zero_count / len(all_bins) * 100.0) if all_bins else 0.0

        lines.append(
            f'<text x="{ox + cell_w // 2}" y="{oy + 14}" '
            f'text-anchor="middle" class="subtitle">'
            f"{file_id} ({zero_pct:.1f}% 0.000000 err)</text>"
        )

        if not vals:
            lines.append(
                f'<text x="{ox + cell_w // 2}" y="{oy + cell_h // 2}" '
                f'text-anchor="middle" class="empty">(no bins)</text>'
            )
            continue

        arr = np.array(vals)
        mn = float(np.mean(arr))
        md = float(np.median(arr))

        nonzero = arr[arr > 0]
        max_val = float(nonzero.max()) if len(nonzero) > 0 else 1.0

        all_max_vals.append(max_val)
        all_n.append(n_total)

        lines.append(
            f'<text x="{ox + cell_w // 2}" y="{oy + 26}" '
            f'text-anchor="middle" class="stats">'
            f"n={n_total}  mean={mn:.6f}  med={md:.6f}</text>"
        )

        n_nonzero_bins = num_buckets - 1
        bucket_width = max_val / n_nonzero_bins if n_nonzero_bins > 0 else 1.0

        counts = np.zeros(num_buckets, dtype=int)
        for v in vals:
            if v == 0.0:
                counts[0] += 1
            else:
                bi = 1 + int((v - 1e-300) / bucket_width)
                if bi >= num_buckets:
                    bi = num_buckets - 1
                counts[bi] += 1

        max_count = int(counts.max()) if counts.max() > 0 else 1

        cx = ox + pad_l
        cy = oy + pad_t
        bar_w = chart_w / num_buckets

        lines.append(
            f'<line x1="{cx}" y1="{cy}" x2="{cx}" y2="{cy + chart_h}" '
            f'stroke="#333" stroke-width="1"/>'
        )
        lines.append(
            f'<line x1="{cx}" y1="{cy + chart_h}" '
            f'x2="{cx + chart_w}" y2="{cy + chart_h}" '
            f'stroke="#333" stroke-width="1"/>'
        )
        lines.append(
            f'<text x="{ox + 18}" y="{cy + chart_h // 2}" '
            f'text-anchor="middle" class="y-axis-label" '
            f'transform="rotate(-90 {ox + 10},{cy + chart_h // 2})">'
            f"# bins (log)</text>"
        )

        bar_h_fn = _make_log_y_axis(lines, cx, cy, chart_h, max_count)

        for b_idx in range(num_buckets):
            c = counts[b_idx]
            if c == 0:
                continue
            bh = bar_h_fn(c)
            bx = cx + b_idx * bar_w
            by = cy + chart_h - bh
            if b_idx == 0:
                tip = f"[0, 0]: {c}"
            else:
                rs = (b_idx - 1) * bucket_width
                re = b_idx * bucket_width
                tip = f"({rs:.6f}, {re:.6f}]: {c}"
            lines.append(
                f'<rect x="{bx + 0.5:.1f}" y="{by:.1f}" '
                f'width="{bar_w - 1:.1f}" height="{bh:.1f}" '
                f'fill="{bar_color}" fill-opacity="0.92">'
                f"<title>{tip}</title></rect>"
            )

        for t in range(3):
            val = max_val * t / 2.0
            x_pos = cx + chart_w * t / 2.0
            lines.append(
                f'<text x="{x_pos:.1f}" y="{cy + chart_h + 16}" '
                f'text-anchor="middle" class="x-tick">{val:.6f}</text>'
            )

        lines.append(
            f'<text x="{cx + chart_w // 2}" y="{cy + chart_h + 30}" '
            f'text-anchor="middle" class="axis-label">SSD (s^2)</text>'
        )

    if all_max_vals:
        global_max = float(max(all_max_vals))
        mean_n = int(round(sum(all_n) / len(all_n)))
        ref_idx = num_files
        ref_col = ref_idx % n_cols
        ref_row = ref_idx // n_cols
        ox_ref = ref_col * cell_w
        oy_ref = ref_row * cell_h + 50
        _append_zero_ref_panel(
            lines,
            ox=ox_ref,
            oy=oy_ref,
            cell_w=cell_w,
            cell_h=cell_h,
            pad_l=pad_l,
            pad_r=pad_r,
            pad_t=pad_t,
            pad_b=pad_b,
            n_total=mean_n,
            ref_max_val=global_max,
            x_label="SSD (s^2)",
            num_buckets=num_buckets,
            fill_color=ref_color,
        )

    lines.append("</svg>")

    with open(path, "w", encoding="utf-8") as fout:
        fout.write("\n".join(lines))


# ============================================================================
# Per-bin offset scatter
# ============================================================================


def write_allfiles_offset_scatter(path, results):
    if not results:
        return

    num_files = len(results)
    n_cols = 5
    n_rows = (num_files + n_cols - 1) // n_cols

    cell_w, cell_h = 260, 275
    pad_l, pad_r, pad_t, pad_b = 55, 25, 48, 50
    chart_w = cell_w - pad_l - pad_r
    chart_h = cell_h - pad_t - pad_b

    svg_w = n_cols * cell_w
    svg_h = n_rows * cell_h + 55

    lines = []
    lines.append(
        f'<svg width="{svg_w}" height="{svg_h}" xmlns="http://www.w3.org/2000/svg">'
    )
    lines.append(
        "<style>\n"
        "  text { font-family: Consolas, 'Courier New', monospace; }\n"
        "  .title { font-size: 16px; font-weight: bold; }\n"
        "  .subtitle { font-size: 12px; font-weight: bold; }\n"
        "  .stats { font-size: 12px; fill: #555; }\n"
        "  .tick { font-size: 11px; }\n"
        "  .axis-label { font-size: 12px; }\n"
        "  .y-axis-label { font-size: 12px; }\n"
        "  .empty { font-size: 12px; fill: #999; }\n"
        "</style>"
    )
    lines.append('<rect width="100%" height="100%" fill="#fcfcfc"/>')
    lines.append(
        f'<text x="{svg_w // 2}" y="30" text-anchor="middle" class="title">'
        f"Per-Bin R-Location Offset (C++ minus MATLAB) for all {num_files} files</text>"
    )

    for f_idx, fr in enumerate(results):
        col = f_idx % n_cols
        row = f_idx // n_cols
        ox = col * cell_w
        oy = row * cell_h + 50

        file_id = fr["id"]
        all_bins = fr["bins"]

        points = []
        for b in all_bins:
            off = b.get("bin_offset_s", float("nan"))
            if not (off != off):
                points.append((b["bin"], off))

        lines.append(
            f'<text x="{ox + cell_w // 2}" y="{oy + 14}" '
            f'text-anchor="middle" class="subtitle">{file_id}</text>'
        )

        if not points:
            lines.append(
                f'<text x="{ox + cell_w // 2}" y="{oy + cell_h // 2}" '
                f'text-anchor="middle" class="empty">(no data)</text>'
            )
            continue

        n_bins_total = len(all_bins)
        offsets = [p[1] for p in points]
        sentinel = -0.1
        display_offsets = [0.0 if v == sentinel else v for v in offsets]
        real_offsets = [v for v in offsets if v != sentinel]

        y_min = min(display_offsets)
        y_max = max(display_offsets)
        if y_max - y_min < 1e-9:
            y_min -= 0.005
            y_max += 0.005
        y_lo = y_min - (y_max - y_min) * 0.15
        y_hi = y_max + (y_max - y_min) * 0.15

        def to_y(val):
            return oy + pad_t + chart_h - (val - y_lo) / (y_hi - y_lo) * chart_h

        def to_x(bin_idx):
            return ox + pad_l + bin_idx / max(n_bins_total - 1, 1) * chart_w

        cx_ax = ox + pad_l
        cy_ax = oy + pad_t

        lines.append(
            f'<line x1="{cx_ax}" y1="{cy_ax}" x2="{cx_ax}" y2="{cy_ax + chart_h}" '
            f'stroke="#333" stroke-width="1"/>'
        )
        lines.append(
            f'<line x1="{cx_ax}" y1="{cy_ax + chart_h}" '
            f'x2="{cx_ax + chart_w}" y2="{cy_ax + chart_h}" '
            f'stroke="#333" stroke-width="1"/>'
        )

        lines.append(
            f'<text x="{ox + 14}" y="{cy_ax + chart_h // 2}" '
            f'text-anchor="middle" class="y-axis-label" '
            f'transform="rotate(-90 {ox + 14},{cy_ax + chart_h // 2})">'
            f"offset (s)</text>"
        )

        n_y_ticks = 4
        y_range = y_max - y_min
        for t in range(n_y_ticks + 1):
            val = y_min + y_range * t / n_y_ticks
            yp = to_y(val)
            lines.append(
                f'<line x1="{cx_ax - 3}" y1="{yp:.1f}" x2="{cx_ax + chart_w}" y2="{yp:.1f}" '
                f'stroke="#eee" stroke-width="0.5"/>'
            )
            lines.append(
                f'<text x="{cx_ax - 5}" y="{yp + 3:.1f}" '
                f'text-anchor="end" class="tick">{val:.6f}</text>'
            )

        REF_LINE = -0.0015
        y_ref_raw = to_y(REF_LINE)
        y_ref_px = max(oy + pad_t, min(oy + pad_t + chart_h, y_ref_raw))
        lines.append(
            f'<line x1="{cx_ax}" y1="{y_ref_px:.1f}" x2="{cx_ax + chart_w}" y2="{y_ref_px:.1f}" '
            f'stroke="#e67e22" stroke-width="1.2" stroke-dasharray="4,3"/>'
        )
        lines.append(
            f'<text x="{cx_ax + chart_w - 2}" y="{y_ref_px - 3:.1f}" '
            f'text-anchor="end" class="tick" fill="#e67e22">-0.0015s</text>'
        )

        n_x_ticks = min(5, n_bins_total)
        for t in range(n_x_ticks + 1):
            bin_i = int(round((n_bins_total - 1) * t / n_x_ticks))
            xp = to_x(bin_i)
            lines.append(
                f'<text x="{xp:.1f}" y="{cy_ax + chart_h + 14}" '
                f'text-anchor="middle" class="tick">{bin_i}</text>'
            )

        lines.append(
            f'<text x="{cx_ax + chart_w // 2}" y="{cy_ax + chart_h + 28}" '
            f'text-anchor="middle" class="axis-label">bin</text>'
        )

        n_inconsistent = sum(1 for v in offsets if v == sentinel)
        n_consistent = len(offsets) - n_inconsistent
        med_off = float(np.median(real_offsets)) if real_offsets else float("nan")
        lines.append(
            f'<text x="{ox + cell_w // 2}" y="{oy + 26}" '
            f'text-anchor="middle" class="stats">'
            f"ok={n_consistent}  incons={n_inconsistent}  med={med_off:.6f}s</text>"
        )

        if y_lo < 0 < y_hi:
            y_zero = to_y(0.0)
            lines.append(
                f'<line x1="{cx_ax}" y1="{y_zero:.1f}" x2="{cx_ax + chart_w}" y2="{y_zero:.1f}" '
                f'stroke="#aaa" stroke-width="0.8" stroke-dasharray="2,2"/>'
            )

        for (bin_i, off), disp in zip(points, display_offsets):
            xp = to_x(bin_i)
            yp = to_y(disp)
            color = "#e74c3c" if off == sentinel else "#1a5fa8"
            label = (
                f"bin {bin_i}: inconsistent (plotted at 0)"
                if off == sentinel
                else f"bin {bin_i}: {off:.6f}s"
            )
            lines.append(
                f'<circle cx="{xp:.1f}" cy="{yp:.1f}" r="3" '
                f'fill="{color}" fill-opacity="0.85">'
                f"<title>{label}</title></circle>"
            )

    lines.append("</svg>")

    with open(path, "w", encoding="utf-8") as fout:
        fout.write("\n".join(lines))


# ============================================================================
# Main
# ============================================================================


def run_pair(pair_label, ids, load_a, load_b, out_dir, label_a, label_b):
    out_dir.mkdir(parents=True, exist_ok=True)
    print("=" * 70)
    print(f"PAIR: {pair_label}")
    print("=" * 70)
    print(f"Output dir: {out_dir}")
    print(f"Subjects:   {len(ids)}\n")

    results = []
    for file_id in ids:
        print(f"  {file_id}...", end="", flush=True)
        try:
            data_a = load_a(file_id)
            data_b = load_b(file_id, len(data_a) if data_a else 0)
        except Exception as e:
            print(f" FAILED ({e})")
            continue
        if not data_a or not data_b:
            print(" SKIPPED (empty data)")
            continue
        fr = compare_file(file_id, data_a, data_b)
        results.append(fr)
        write_per_file_csv(
            out_dir / f"{file_id}_wave_comparison.csv", fr, label_a, label_b
        )
        print(f" {fr['total_bins']} bins")

    if not results:
        print("  No subjects produced results.\n")
        return results

    write_summary_csv(out_dir / "summary.csv", results, label_a, label_b)
    write_file_summary_csv(out_dir / "file_summary.csv", results)

    write_per_subject_histograms(
        out_dir / "rr_interval_histograms.svg",
        results,
        field="rr_ssd",
        chart_title=f"Comparison of R Peak Intervals between {label_a} and {label_b}'s Peak Finding Code",
        bar_color="#1a5fa8",
        ref_color="#2ecc71",
    )
    write_per_subject_histograms(
        out_dir / "rloc_histograms.svg",
        results,
        field="rloc_ssd",
        chart_title=f"Comparison of R Peak Locations between {label_a} and {label_b}'s Peak Finding Code",
        bar_color="#e67e22",
        ref_color="#9b59b6",
    )

    write_combined_histogram(
        out_dir / "combined_rr_histogram.svg",
        results,
        field="rr_ssd",
        chart_title=f"{label_a} vs {label_b}: RR Interval SSD (pooled)",
    )
    write_combined_histogram(
        out_dir / "combined_rloc_histogram.svg",
        results,
        field="rloc_ssd",
        chart_title=f"{label_a} vs {label_b}: R-Location SSD (pooled)",
    )

    write_allfiles_offset_scatter(out_dir / "rloc_offset_scatter.svg", results)

    print(f"  Wrote outputs to {out_dir}\n")
    return results


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("BIN_DIR exists:", BIN_DIR.exists())
    print("MAT_DIR exists:", MAT_DIR.exists())
    print("CSV_DIR exists:", CSV_DIR.exists())

    mira = {}
    for p in BIN_DIR.glob("*_wave_markings.bin"):
        mira[p.stem.replace("_wave_markings", "")] = p
    daniel = {}
    if MAT_DIR.exists():
        for p in MAT_DIR.glob("*_wave_data.mat"):
            daniel[p.stem.replace("_wave_data", "")] = p
    deep = {}
    deep_full_suffixes = [
        "_part1_ECG_fs1000_wholecaseRRiQTi.csv",
        "_ECG_fs1000_wholecaseRRiQTi.csv",
        "_wholecaseRRiQTi.csv",
    ]
    if CSV_DIR.exists():
        for p in CSV_DIR.glob(f"*{DEEP_FILE_SUFFIX}"):
            if not p.is_file():
                continue
            name = p.name
            subject_id = None
            for suf in deep_full_suffixes:
                if name.endswith(suf):
                    subject_id = name[: -len(suf)]
                    break
            if subject_id is None:
                subject_id = name[: -len(DEEP_FILE_SUFFIX)]
            deep.setdefault(subject_id, p)

    print(f"\nMira subjects:   {len(mira)}")
    print(f"Daniel subjects: {len(daniel)}")
    print(f"Deep subjects:   {len(deep)}\n")

    mira_bins_cache = {}
    daniel_bins_cache = {}
    deep_bins_cache = {}

    def get_mira_bins(fid):
        if fid not in mira_bins_cache:
            if fid in mira:
                mira_bins_cache[fid] = read_bin_file(mira[fid])
            else:
                mira_bins_cache[fid] = []
        return mira_bins_cache[fid]

    def get_daniel_bins(fid):
        if fid not in daniel_bins_cache:
            if fid in daniel:
                daniel_bins_cache[fid] = read_mat_file(daniel[fid])
            else:
                daniel_bins_cache[fid] = []
        return daniel_bins_cache[fid]

    def get_deep_bins(fid):
        if fid not in deep_bins_cache:
            if fid in deep:
                # Still pass mira_bins so we can compute per-bin offsets for
                # slicing; Deep's timestamps are now self-contained (no anchor).
                bins = get_mira_bins(fid)
                deep_bins_cache[fid] = read_deep_file(deep[fid], bins, len(bins))
            else:
                deep_bins_cache[fid] = []
        return deep_bins_cache[fid]

    def min_bins_across_sources(fid):
        counts = []
        if fid in mira:
            counts.append(len(get_mira_bins(fid)))
        if fid in daniel:
            counts.append(len(get_daniel_bins(fid)))
        if fid in deep:
            counts.append(len(get_deep_bins(fid)))
        return min(counts) if counts else 0

    def load_mira(fid, *_):
        return get_mira_bins(fid)[: min_bins_across_sources(fid)]

    def load_daniel(fid, *_):
        return get_daniel_bins(fid)[: min_bins_across_sources(fid)]

    def load_deep(fid, *_):
        return get_deep_bins(fid)[: min_bins_across_sources(fid)]

    pair_results = {}
    ids = sorted(set(mira) & set(daniel))
    if ids:
        r = run_pair(
            "Mira vs Daniel",
            ids,
            load_mira,
            load_daniel,
            OUTPUT_DIR / "mira_vs_daniel",
            "Mira",
            "Daniel",
        )
        if r:
            pair_results[("Mira", "Daniel")] = r

    ids = sorted(set(mira) & set(deep))
    if ids:
        r = run_pair(
            "Mira vs Deep",
            ids,
            load_mira,
            load_deep,
            OUTPUT_DIR / "mira_vs_deep",
            "Mira",
            "Deep",
        )
        if r:
            pair_results[("Mira", "Deep")] = r

    ids = sorted(set(daniel) & set(deep))
    if ids:
        r = run_pair(
            "Daniel vs Deep",
            ids,
            load_daniel,
            load_deep,
            OUTPUT_DIR / "daniel_vs_deep",
            "Daniel",
            "Deep",
        )
        if r:
            pair_results[("Daniel", "Deep")] = r

    summaries = {}
    if ("Mira", "Daniel") in pair_results:
        summaries["Mira"] = summarize_paired_source(
            pair_results[("Mira", "Daniel")], "cpp"
        )
    elif ("Mira", "Deep") in pair_results:
        summaries["Mira"] = summarize_paired_source(
            pair_results[("Mira", "Deep")], "cpp"
        )

    if ("Mira", "Daniel") in pair_results:
        summaries["Daniel"] = summarize_paired_source(
            pair_results[("Mira", "Daniel")], "mat"
        )
    elif ("Daniel", "Deep") in pair_results:
        summaries["Daniel"] = summarize_paired_source(
            pair_results[("Daniel", "Deep")], "cpp"
        )

    if ("Mira", "Deep") in pair_results:
        summaries["Deep"] = summarize_paired_source(
            pair_results[("Mira", "Deep")], "mat"
        )
    elif ("Daniel", "Deep") in pair_results:
        summaries["Deep"] = summarize_paired_source(
            pair_results[("Daniel", "Deep")], "mat"
        )

    if summaries:
        pair_diffs = {}
        for (src_a, src_b), pr in pair_results.items():
            extra_a = sum(
                max(b["n_ecg_cpp"] - b["n_ecg_mat"], 0) for fr in pr for b in fr["bins"]
            )
            extra_b = sum(
                max(b["n_ecg_mat"] - b["n_ecg_cpp"], 0) for fr in pr for b in fr["bins"]
            )
            pair_diffs[(src_a, src_b)] = extra_a
            pair_diffs[(src_b, src_a)] = extra_b

        write_three_way_summary_csv(
            OUTPUT_DIR / "three_way_summary.csv", summaries, pair_diffs
        )
        print(f"Wrote three_way_summary.csv with sources: {list(summaries.keys())}")

    print("=" * 70)
    print("ALL PAIRS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
