"""
compare_wave_bounds.py
Compare wave bound outputs between C++ (_wave_markings.bin)
and MATLAB (_wave_data.mat).
"""

import csv
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
OUTPUT_DIR = Path(r"D:\USERS\MiraWelner\QTVI\testing\4_wave_finding_tests\results")

SR = 1000.0
MAX_SANE = 50_000_000
PASSTHROUGH_NUM_CHANNELS = 41

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
            skip_idx()
            skip_idx()
            skip_idx()
            skip_idx()
            skip_idx()
            skip_idx()
            skip_idx()
            skip_idx()
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

        r["ecgRIndex"] = get_field("ecgRIndex").astype(np.intp)
        r["ecgRIndex_f"] = get_field(
            "ecgRIndex"
        )  # float64, preserves sub-sample precision
        r["ppgMinAmps"] = get_field("ppgMinAmps").astype(np.intp)
        r["ecgSignal"] = (
            get_field("ecg")
            if "ecg" in (obj.dtype.names or [])
            else get_field("ecgSignal")
        )
        r["ppgSignal"] = (
            get_field("po")
            if "po" in (obj.dtype.names or [])
            else get_field("ppgSignal")
        )
        results.append(r)

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


def rr_ssd(bin_idx, mat_idx):
    rr_b = np.diff(np.array(bin_idx, dtype=np.float64))
    rr_m = np.diff(np.array(mat_idx, dtype=np.float64))
    if len(rr_b) == 0 or len(rr_m) == 0:
        return np.nan
    n = min(len(rr_b), len(rr_m))
    return float(np.sum((rr_b[:n] - rr_m[:n]) ** 2)) / (SR * SR)


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

        b_rr_raw = np.diff(np.array(b["ecgRIndex"], dtype=np.float64))
        m_rr_raw = np.diff(np.array(m["ecgRIndex"], dtype=np.float64))
        if len(b_rr_raw) == len(m_rr_raw) and len(b_rr_raw) > 0:
            rr_identical = np.array_equal(b_rr_raw, m_rr_raw)
        else:
            rr_identical = False

        # Per-bin SSD of R-location differences (s^2)
        # sum((r_cpp[i] - r_mat[i])^2) / SR^2 across matched peaks
        b_rloc = np.array(b["ecgRIndex"], dtype=np.float64)
        m_rloc_f = np.array(m.get("ecgRIndex_f", m["ecgRIndex"]), dtype=np.float64)
        n_rloc = min(len(b_rloc), len(m_rloc_f))
        if n_rloc > 0:
            rloc_ssd = float(np.sum((b_rloc[:n_rloc] - m_rloc_f[:n_rloc]) ** 2)) / (
                SR * SR
            )
        else:
            rloc_ssd = float("nan")

        # Per-bin offset (for scatter): constant R-location offset if consistent
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
                "rr_ssd": rr_ssd(b["ecgRIndex"], m["ecgRIndex"]),
                "rr_identical": rr_identical,
                "rloc_ssd": rloc_ssd,
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


def fmt(v, d=3):
    return "" if (v is None or (isinstance(v, float) and np.isnan(v))) else round(v, d)


def iqr_str(q1, q3, d=4):
    if np.isnan(q1) or np.isnan(q3):
        return ""
    return f"{round(q1, d)},{round(q3, d)}"


def write_summary_csv(path, results):
    all_bins = [b for fr in results for b in fr["bins"]]

    def collect(key_fn):
        return [key_fn(b) for b in all_bins if not np.isnan(key_fn(b))]

    cols = ["C++ ECG", "MATLAB ECG", "C++ PPG", "MATLAB PPG"]

    def make_row(label, vals):
        return [label] + [fmt(v) for v in vals]

    def make_iqr_row(label, vals_q1, vals_q3):
        return [label] + [iqr_str(q1, q3) for q1, q3 in zip(vals_q1, vals_q3)]

    n_subjects = len(results)
    n_beats = [
        sum(b["n_ecg_cpp"] for b in all_bins),
        sum(b["n_ecg_mat"] for b in all_bins),
        sum(b["n_ppg_cpp"] for b in all_bins),
        sum(b["n_ppg_mat"] for b in all_bins),
    ]

    rows = [
        ["Metric"] + cols,
        make_row("N Subjects", [n_subjects] * 4),
        make_row("N Beats Total", n_beats),
        make_row(
            "Extra Beats MATLAB",
            [
                sum(max(b["n_ecg_mat"] - b["n_ecg_cpp"], 0) for b in all_bins),
                sum(max(b["n_ecg_mat"] - b["n_ecg_cpp"], 0) for b in all_bins),
                sum(max(b["n_ppg_mat"] - b["n_ppg_cpp"], 0) for b in all_bins),
                sum(max(b["n_ppg_mat"] - b["n_ppg_cpp"], 0) for b in all_bins),
            ],
        ),
        make_row(
            "Extra Beats C++",
            [
                sum(max(b["n_ecg_cpp"] - b["n_ecg_mat"], 0) for b in all_bins),
                sum(max(b["n_ecg_cpp"] - b["n_ecg_mat"], 0) for b in all_bins),
                sum(max(b["n_ppg_cpp"] - b["n_ppg_mat"], 0) for b in all_bins),
                sum(max(b["n_ppg_cpp"] - b["n_ppg_mat"], 0) for b in all_bins),
            ],
        ),
        make_row(
            "Mean N Beats Per Bin",
            [
                nanmean(collect(lambda b: b["n_ecg_cpp"])),
                nanmean(collect(lambda b: b["n_ecg_mat"])),
                nanmean(collect(lambda b: b["n_ppg_cpp"])),
                nanmean(collect(lambda b: b["n_ppg_mat"])),
            ],
        ),
        make_row(
            "Mean N Beats Per Patient",
            [
                nanmean([sum(b["n_ecg_cpp"] for b in fr["bins"]) for fr in results]),
                nanmean([sum(b["n_ecg_mat"] for b in fr["bins"]) for fr in results]),
                nanmean([sum(b["n_ppg_cpp"] for b in fr["bins"]) for fr in results]),
                nanmean([sum(b["n_ppg_mat"] for b in fr["bins"]) for fr in results]),
            ],
        ),
        make_row(
            "Mean Rate (bpm)",
            [
                nanmean(collect(lambda b: b["rr_cpp"]["rate_mean"])),
                nanmean(collect(lambda b: b["rr_mat"]["rate_mean"])),
                nanmean(collect(lambda b: b["ppg_rr_cpp"]["rate_mean"])),
                nanmean(collect(lambda b: b["ppg_rr_mat"]["rate_mean"])),
            ],
        ),
        make_row(
            "Median Rate (bpm)",
            [
                nanmedian(collect(lambda b: b["rr_cpp"]["rate_median"])),
                nanmedian(collect(lambda b: b["rr_mat"]["rate_median"])),
                nanmedian(collect(lambda b: b["ppg_rr_cpp"]["rate_median"])),
                nanmedian(collect(lambda b: b["ppg_rr_mat"]["rate_median"])),
            ],
        ),
        make_iqr_row(
            "Rate IQR (bpm)",
            [
                nanpercentile(collect(lambda b: b["rr_cpp"]["rate_q1"]), 25),
                nanpercentile(collect(lambda b: b["rr_mat"]["rate_q1"]), 25),
                nanpercentile(collect(lambda b: b["ppg_rr_cpp"]["rate_q1"]), 25),
                nanpercentile(collect(lambda b: b["ppg_rr_mat"]["rate_q1"]), 25),
            ],
            [
                nanpercentile(collect(lambda b: b["rr_cpp"]["rate_q3"]), 75),
                nanpercentile(collect(lambda b: b["rr_mat"]["rate_q3"]), 75),
                nanpercentile(collect(lambda b: b["ppg_rr_cpp"]["rate_q3"]), 75),
                nanpercentile(collect(lambda b: b["ppg_rr_mat"]["rate_q3"]), 75),
            ],
        ),
        make_row(
            "Mean Interval (ms)",
            [
                nanmean(collect(lambda b: b["rr_cpp"]["mean"])),
                nanmean(collect(lambda b: b["rr_mat"]["mean"])),
                nanmean(collect(lambda b: b["ppg_rr_cpp"]["mean"])),
                nanmean(collect(lambda b: b["ppg_rr_mat"]["mean"])),
            ],
        ),
        make_row(
            "Median Interval (ms)",
            [
                nanmedian(collect(lambda b: b["rr_cpp"]["median"])),
                nanmedian(collect(lambda b: b["rr_mat"]["median"])),
                nanmedian(collect(lambda b: b["ppg_rr_cpp"]["median"])),
                nanmedian(collect(lambda b: b["ppg_rr_mat"]["median"])),
            ],
        ),
        make_iqr_row(
            "Interval IQR (ms)",
            [
                nanpercentile(collect(lambda b: b["rr_cpp"]["q1"]), 25),
                nanpercentile(collect(lambda b: b["rr_mat"]["q1"]), 25),
                nanpercentile(collect(lambda b: b["ppg_rr_cpp"]["q1"]), 25),
                nanpercentile(collect(lambda b: b["ppg_rr_mat"]["q1"]), 25),
            ],
            [
                nanpercentile(collect(lambda b: b["rr_cpp"]["q3"]), 75),
                nanpercentile(collect(lambda b: b["rr_mat"]["q3"]), 75),
                nanpercentile(collect(lambda b: b["ppg_rr_cpp"]["q3"]), 75),
                nanpercentile(collect(lambda b: b["ppg_rr_mat"]["q3"]), 75),
            ],
        ),
    ]

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        for row in rows:
            w.writerow(row)


# ============================================================================
# File summary CSV (per-subject R-peak and PPG counts)
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


def write_per_file_csv(path, fr):
    bins = fr["bins"]

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

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "bin",
                "R Peaks C++",
                "R Peaks MATLAB",
                "Mean RR C++ (ms)",
                "Mean RR MATLAB (ms)",
                "Median RR C++ (ms)",
                "Median RR MATLAB (ms)",
                "IQR RR C++ (q1,q3)",
                "IQR RR MATLAB (q1,q3)",
                "RR SSD C++ vs MATLAB (s²)",
                "R-Loc SSD C++ vs MATLAB (s²)",
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
# Shared log-y-axis helper for histograms
# ============================================================================

import math as _math


def _make_log_y_axis(lines, cx, cy, chart_h, max_count):
    """
    Render log-scale y-axis ticks at 0, 1, 10, 100, ... up to max_count.
    Returns bar_h(count) -> SVG bar height in pixels.
    The gridlines are clipped to the chart area (cx .. cx+chart_w is not
    known here, so gridlines are drawn long and clipped by the SVG viewport).
    """
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

    # Add the actual max_count as an explicit tick if it doesn't fall on a power of 10
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


# ============================================================================
# Combined SSD Histogram (all files pooled into one large chart)
# ============================================================================


def write_combined_ssd_histogram(path, results, num_buckets=80):
    all_ssd = []
    for fr in results:
        for b in fr["bins"]:
            if not np.isnan(b["rr_ssd"]):
                all_ssd.append(round(b["rr_ssd"], 3))

    if not all_ssd:
        return

    arr = np.array(all_ssd)
    mn = float(np.mean(arr))
    md = float(np.median(arr))
    q1 = float(np.percentile(arr, 25))
    q3 = float(np.percentile(arr, 75))
    zero_count = sum(1 for v in all_ssd if v == 0.0)
    total_bins_count = len(all_ssd)
    zero_pct = zero_count / total_bins_count * 100.0 if total_bins_count > 0 else 0.0

    svg_w, svg_h = 900, 500
    pad_l, pad_r, pad_t, pad_b = 70, 40, 70, 60
    chart_w = svg_w - pad_l - pad_r
    chart_h = svg_h - pad_t - pad_b

    min_val = float(arr.min())
    max_val = float(arr.max())
    if max_val - min_val < 1e-15:
        max_val = min_val + 1.0
    bucket_width = (max_val - min_val) / num_buckets

    counts = np.zeros(num_buckets, dtype=int)
    for v in all_ssd:
        bi = int((v - min_val) / bucket_width)
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
        f"Combined ECG SSD Histogram - All {len(results)} Subjects "
        f"({len(all_ssd)} bins, {zero_pct:.1f}% zero err)</text>"
    )
    lines.append(
        f'<text x="{svg_w // 2}" y="44" text-anchor="middle" class="stats">'
        f"mean={mn:.6f}  median={md:.6f}  "
        f"IQR=[{q1:.6f}, {q3:.6f}]</text>"
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
        range_start = min_val + b_idx * bucket_width
        range_end = min_val + (b_idx + 1) * bucket_width
        tip = f"[{range_start:.4f}, {range_end:.4f}): {c}"
        lines.append(
            f'<rect x="{bx + 0.5:.1f}" y="{by:.1f}" '
            f'width="{bar_w - 1:.1f}" height="{bh:.1f}" '
            f'fill="#3498db" fill-opacity="0.85">'
            f"<title>{tip}</title></rect>"
        )

    n_x_ticks = 5
    for t in range(n_x_ticks + 1):
        val = min_val + (max_val - min_val) * t / n_x_ticks
        x_pos = cx + chart_w * t / n_x_ticks
        lines.append(
            f'<text x="{x_pos:.1f}" y="{cy + chart_h + 16}" '
            f'text-anchor="middle" class="tick">{val:.4f}</text>'
        )
    lines.append(
        f'<text x="{cx + chart_w // 2}" y="{cy + chart_h + 38}" '
        f'text-anchor="middle" class="axis-label">SSD (s^2)</text>'
    )
    lines.append("</svg>")

    with open(path, "w", encoding="utf-8") as fout:
        fout.write("\n".join(lines))


# ============================================================================
# Helper: render a single "zero-difference" reference histogram panel into
# an existing lines list.  All N values are placed in bin 0; the x-axis
# spans [0, ref_max_val] so it matches the neighbouring real histogram.
# ============================================================================


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

    # Title row
    lines.append(
        f'<text x="{ox + cell_w // 2}" y="{oy + 14}" '
        f'text-anchor="middle" class="subtitle" fill="#27ae60">'
        f"ZERO-DIFF REFERENCE (n={n_total})</text>"
    )
    lines.append(
        f'<text x="{ox + cell_w // 2}" y="{oy + 26}" '
        f'text-anchor="middle" class="stats">'
        f"all diffs=0 | x-range=[0, {ref_max_val:.4f}]</text>"
    )

    # Axes
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

    # Y-axis ticks
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

    # Single bar in bin 0
    bar_w = chart_w / num_buckets
    bx = cx  # bin 0
    bar_h = chart_h  # full height (all values here)
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

    # X-axis ticks (same scale as real chart: 0 to ref_max_val)
    for t in range(3):
        val = ref_max_val * t / 2.0
        x_pos = cx + chart_w * t / 2.0
        lines.append(
            f'<text x="{x_pos:.1f}" y="{cy + chart_h + 16}" '
            f'text-anchor="middle" class="x-tick">{val:.4f}</text>'
        )

    lines.append(
        f'<text x="{cx + chart_w // 2}" y="{cy + chart_h + 30}" '
        f'text-anchor="middle" class="axis-label">{x_label}</text>'
    )


# ============================================================================
# Per-file SSD Histograms with zero-reference panels (RR interval SSD)
# ============================================================================


def write_allfiles_ssd_histograms(path, results, num_buckets=50):
    """
    One chart per subject (5 columns), plus a single zero-difference reference
    panel placed in the next available grid cell after the last subject.
    Bin 0 collects only exact-zero SSD values; remaining bins divide (0, max].
    """
    if not results:
        return

    num_files = len(results)
    n_cols = 5
    # +1 for the zero-ref panel placed in the next grid slot
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
        f"MESA Summary Statistics histograms for all {num_files} files</text>"
    )

    all_max_vals = []
    all_n_ssd = []

    for f_idx, fr in enumerate(results):
        col = f_idx % n_cols
        row = f_idx // n_cols
        ox = col * cell_w
        oy = row * cell_h + 50

        file_id = fr["id"]
        all_bins = fr["bins"]

        ssd_vals = [
            round(b["rr_ssd"], 3) for b in all_bins if not np.isnan(b["rr_ssd"])
        ]
        zero_count = sum(1 for v in ssd_vals if v == 0.0)
        zero_pct = (zero_count / len(all_bins) * 100.0) if all_bins else 0.0

        lines.append(
            f'<text x="{ox + cell_w // 2}" y="{oy + 14}" '
            f'text-anchor="middle" class="subtitle">'
            f"{file_id} ({zero_pct:.1f}% zero err)</text>"
        )

        if not ssd_vals:
            lines.append(
                f'<text x="{ox + cell_w // 2}" y="{oy + cell_h // 2}" '
                f'text-anchor="middle" class="empty">(no bins)</text>'
            )
            continue

        arr = np.array(ssd_vals)
        mn = float(np.mean(arr))
        md = float(np.median(arr))

        nonzero = arr[arr > 0]
        max_val = float(nonzero.max()) if len(nonzero) > 0 else 1.0

        all_max_vals.append(max_val)
        all_n_ssd.append(len(ssd_vals))

        lines.append(
            f'<text x="{ox + cell_w // 2}" y="{oy + 26}" '
            f'text-anchor="middle" class="stats">'
            f"n={len(ssd_vals)}  mean={mn:.4f}  med={md:.4f}</text>"
        )

        # Bin 0 = exact zeros; bins 1..num_buckets-1 divide (0, max_val]
        n_nonzero_bins = num_buckets - 1
        bucket_width = max_val / n_nonzero_bins if n_nonzero_bins > 0 else 1.0

        counts = np.zeros(num_buckets, dtype=int)
        for v in ssd_vals:
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
                tip = f"({rs:.4f}, {re:.4f}]: {c}"
            lines.append(
                f'<rect x="{bx + 0.5:.1f}" y="{by:.1f}" '
                f'width="{bar_w - 1:.1f}" height="{bh:.1f}" '
                f'fill="#1a5fa8" fill-opacity="0.92">'
                f"<title>{tip}</title></rect>"
            )

        for t in range(3):
            val = max_val * t / 2.0
            x_pos = cx + chart_w * t / 2.0
            lines.append(
                f'<text x="{x_pos:.1f}" y="{cy + chart_h + 16}" '
                f'text-anchor="middle" class="x-tick">{val:.4f}</text>'
            )

        lines.append(
            f'<text x="{cx + chart_w // 2}" y="{cy + chart_h + 30}" '
            f'text-anchor="middle" class="axis-label">SSD (s^2)</text>'
        )

    # ---- Zero-reference panel in next grid slot ----
    if all_max_vals:
        global_max = float(max(all_max_vals))
        mean_n = int(round(sum(all_n_ssd) / len(all_n_ssd)))
        ref_idx = num_files  # next slot after all subjects
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
            fill_color="#2ecc71",
        )

    lines.append("</svg>")

    with open(path, "w", encoding="utf-8") as fout:
        fout.write("\n".join(lines))


# ============================================================================
# Per-file R-location difference Histograms (raw sample offsets, not intervals)
# ============================================================================


def write_allfiles_rloc_histograms(path, results, num_buckets=50):
    """
    One chart per subject (5 columns) of |r_cpp - r_mat| (raw R-location diffs),
    plus a single zero-difference reference panel in the next grid slot.
    Bin 0 collects only exact-zero differences; remaining bins divide (0, max].
    """
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
        f"MESA Per-Bin R-Location SSD Histograms for all {num_files} files</text>"
    )

    all_max_vals = []
    all_n_beats = []

    for f_idx, fr in enumerate(results):
        col = f_idx % n_cols
        row = f_idx // n_cols
        ox = col * cell_w
        oy = row * cell_h + 50

        file_id = fr["id"]
        all_bins = fr["bins"]

        rloc_vals = [
            round(b["rloc_ssd"], 3)
            for b in all_bins
            if not np.isnan(b.get("rloc_ssd", float("nan")))
        ]

        n_total_beats = len(rloc_vals)  # actually n bins here
        zero_count = sum(1 for v in rloc_vals if v == 0.0)
        zero_pct = (zero_count / n_total_beats * 100.0) if n_total_beats > 0 else 0.0

        lines.append(
            f'<text x="{ox + cell_w // 2}" y="{oy + 14}" '
            f'text-anchor="middle" class="subtitle">'
            f"{file_id} ({zero_pct:.1f}% zero err)</text>"
        )

        if not rloc_vals:
            lines.append(
                f'<text x="{ox + cell_w // 2}" y="{oy + cell_h // 2}" '
                f'text-anchor="middle" class="empty">(no data)</text>'
            )
            continue

        arr = np.array(rloc_vals)
        mn = float(np.mean(arr))
        md = float(np.median(arr))

        nonzero_rloc = arr[arr > 0]
        max_val = float(nonzero_rloc.max()) if len(nonzero_rloc) > 0 else 1.0

        all_max_vals.append(max_val)
        all_n_beats.append(n_total_beats)

        lines.append(
            f'<text x="{ox + cell_w // 2}" y="{oy + 26}" '
            f'text-anchor="middle" class="stats">'
            f"n={n_total_beats}  mean={mn:.5f}  med={md:.5f}</text>"
        )

        # Bin 0 = exact zeros; bins 1..num_buckets-1 divide (0, max_val]
        n_nonzero_bins = num_buckets - 1
        bucket_width = max_val / n_nonzero_bins if n_nonzero_bins > 0 else 1.0

        counts = np.zeros(num_buckets, dtype=int)
        for v in rloc_vals:
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
                tip = f"({rs:.5f}, {re:.5f}]: {c}"
            lines.append(
                f'<rect x="{bx + 0.5:.1f}" y="{by:.1f}" '
                f'width="{bar_w - 1:.1f}" height="{bh:.1f}" '
                f'fill="#e67e22" fill-opacity="0.8">'
                f"<title>{tip}</title></rect>"
            )

        for t in range(3):
            val = max_val * t / 2.0
            x_pos = cx + chart_w * t / 2.0
            lines.append(
                f'<text x="{x_pos:.1f}" y="{cy + chart_h + 16}" '
                f'text-anchor="middle" class="x-tick">{val:.4f}</text>'
            )

        lines.append(
            f'<text x="{cx + chart_w // 2}" y="{cy + chart_h + 30}" '
            f'text-anchor="middle" class="axis-label">SSD (s^2)</text>'
        )

    # ---- Zero-reference panel in next grid slot ----
    if all_max_vals:
        global_max = float(max(all_max_vals))
        mean_n = int(round(sum(all_n_beats) / len(all_n_beats)))
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
            fill_color="#9b59b6",
        )

    lines.append("</svg>")

    with open(path, "w", encoding="utf-8") as fout:
        fout.write("\n".join(lines))


# ============================================================================
# Per-file R-location offset scatterplots (one small chart per subject)
# ============================================================================


def write_allfiles_offset_scatter(path, results):
    """
    One scatterplot per subject (5 columns).
    X axis = bin index, Y axis = signed R-location offset in seconds.
    Bins with inconsistent offsets are plotted at y = -0.1 in red.
    Bins with no data are skipped.
    """
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
        f"MESA Per-Bin R-Location Offset (C++ minus MATLAB) for all {num_files} files</text>"
    )

    for f_idx, fr in enumerate(results):
        col = f_idx % n_cols
        row = f_idx // n_cols
        ox = col * cell_w
        oy = row * cell_h + 50

        file_id = fr["id"]
        all_bins = fr["bins"]

        # Collect (bin_index, offset_s) pairs, skip nan
        points = []
        for b in all_bins:
            off = b.get("bin_offset_s", float("nan"))
            if not (off != off):  # skip nan
                points.append((b["bin"], off))

        lines.append(
            f'<text x="{ox + cell_w // 2}" y="{oy + 14}" '
            f'text-anchor="middle" class="subtitle">'
            f"{file_id}</text>"
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
        # Inconsistent bins display at 0
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

        # Axes
        lines.append(
            f'<line x1="{cx_ax}" y1="{cy_ax}" x2="{cx_ax}" y2="{cy_ax + chart_h}" '
            f'stroke="#333" stroke-width="1"/>'
        )
        lines.append(
            f'<line x1="{cx_ax}" y1="{cy_ax + chart_h}" '
            f'x2="{cx_ax + chart_w}" y2="{cy_ax + chart_h}" '
            f'stroke="#333" stroke-width="1"/>'
        )

        # Y axis label
        lines.append(
            f'<text x="{ox + 14}" y="{cy_ax + chart_h // 2}" '
            f'text-anchor="middle" class="y-axis-label" '
            f'transform="rotate(-90 {ox + 14},{cy_ax + chart_h // 2})">'
            f"offset (s)</text>"
        )

        # Y ticks: show a few real values
        n_y_ticks = 4
        y_range = y_max - y_min
        # Choose precision based on range magnitude
        if y_range < 0.001:
            tick_fmt = ".5f"
        elif y_range < 0.01:
            tick_fmt = ".4f"
        else:
            tick_fmt = ".3f"
        for t in range(n_y_ticks + 1):
            val = y_min + y_range * t / n_y_ticks
            yp = to_y(val)
            lines.append(
                f'<line x1="{cx_ax - 3}" y1="{yp:.1f}" x2="{cx_ax + chart_w}" y2="{yp:.1f}" '
                f'stroke="#eee" stroke-width="0.5"/>'
            )
            lines.append(
                f'<text x="{cx_ax - 5}" y="{yp + 3:.1f}" '
                f'text-anchor="end" class="tick">{val:{tick_fmt}}</text>'
            )

        # Dashed reference line at -0.0015 s (always drawn; clamped to chart if outside range)
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

        # X ticks
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

        # Stats line
        n_inconsistent = sum(1 for v in offsets if v == sentinel)
        n_consistent = len(offsets) - n_inconsistent
        med_off = float(np.median(real_offsets)) if real_offsets else float("nan")
        lines.append(
            f'<text x="{ox + cell_w // 2}" y="{oy + 26}" '
            f'text-anchor="middle" class="stats">'
            f"ok={n_consistent}  incons={n_inconsistent}  med={med_off:.4f}s</text>"
        )

        # Draw zero line if in range
        if y_lo < 0 < y_hi:
            y_zero = to_y(0.0)
            lines.append(
                f'<line x1="{cx_ax}" y1="{y_zero:.1f}" x2="{cx_ax + chart_w}" y2="{y_zero:.1f}" '
                f'stroke="#aaa" stroke-width="0.8" stroke-dasharray="2,2"/>'
            )

        # Points (inconsistent bins plot at y=0, shown in orange)
        for (bin_i, off), disp in zip(points, display_offsets):
            xp = to_x(bin_i)
            yp = to_y(disp)
            color = "#e74c3c" if off == sentinel else "#1a5fa8"
            label = (
                f"bin {bin_i}: inconsistent (plotted at 0)"
                if off == sentinel
                else f"bin {bin_i}: {off:.5f}s"
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


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("BIN_DIR exists:", BIN_DIR.exists())
    print("MAT_DIR exists:", MAT_DIR.exists())

    ids = []
    for bin_file in BIN_DIR.glob("*_wave_markings.bin"):
        file_id = bin_file.stem.replace("_wave_markings", "")
        if (MAT_DIR / f"{file_id}_wave_data.mat").exists():
            ids.append(file_id)
    ids.sort()

    print("=" * 70)
    print("WAVE BOUNDS COMPARISON: C++ vs MATLAB")
    print("=" * 70)
    print(f"C++ dir:    {BIN_DIR}")
    print(f"MATLAB dir: {MAT_DIR}")
    print(f"Output dir: {OUTPUT_DIR}")
    print(f"Matched:    {len(ids)} files\n")

    results = []
    for file_id in ids:
        print(f"  {file_id}...", end="", flush=True)
        bin_data = read_bin_file(BIN_DIR / f"{file_id}_wave_markings.bin")
        mat_data = read_mat_file(MAT_DIR / f"{file_id}_wave_data.mat")

        if not bin_data or not mat_data:
            print(" SKIPPED (empty data)")
            continue

        fr = compare_file(file_id, bin_data, mat_data)
        results.append(fr)
        write_per_file_csv(OUTPUT_DIR / f"{file_id}_wave_comparison.csv", fr)
        print(f" {fr['total_bins']} bins")

    write_summary_csv(OUTPUT_DIR / "summary.csv", results)
    write_file_summary_csv(OUTPUT_DIR / "file_summary.csv", results)
    write_allfiles_ssd_histograms(OUTPUT_DIR / "allfiles_ssd_histograms.svg", results)
    write_allfiles_rloc_histograms(
        OUTPUT_DIR / "allfiles_ssd_histograms_rr.svg", results
    )
    write_allfiles_offset_scatter(OUTPUT_DIR / "allfiles_offset_scatter.svg", results)
    write_combined_ssd_histogram(OUTPUT_DIR / "combined_ssd_histogram.svg", results)

    print(f"\nSummary CSV:       {OUTPUT_DIR / 'summary.csv'}")
    print(f"File Summary CSV:  {OUTPUT_DIR / 'file_summary.csv'}")
    print(f"SSD Per-File:      {OUTPUT_DIR / 'allfiles_ssd_histograms.svg'}")
    print(f"RLoc Per-File:     {OUTPUT_DIR / 'allfiles_ssd_histograms_rr.svg'}")
    print(f"Offset Scatter:    {OUTPUT_DIR / 'allfiles_offset_scatter.svg'}")
    print(f"SSD Combined:      {OUTPUT_DIR / 'combined_ssd_histogram.svg'}")
    print(f"Per-file CSVs:     {OUTPUT_DIR}/<subject>_wave_comparison.csv")


if __name__ == "__main__":
    main()
