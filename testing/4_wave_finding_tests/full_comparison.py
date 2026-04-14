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

BIN_SR = 1000.0
MAT_SR = 1000.0
MAX_SANE = 50_000_000

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
            b["ecgRIndex"] = read_idx()  # ch1 raw
            skip_idx()
            skip_idx()  # ch1 sq, abs
            skip_idx()
            skip_idx()
            skip_idx()  # ch2
            skip_idx()
            skip_idx()
            skip_idx()  # ch3
            skip_idx()  # ppgMaxAmps
            b["ppgMinAmps"] = read_idx()  # ppgMinAmps
            b["ppgSignal"] = read_signal()
            b["ecgSignal"] = read_signal()
            skip_signal()
            skip_signal()
            for _ in range(6):
                skip_signal()
            f.read(9)  # noise flags

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

            read_pair_vec()  # ppg_bin_indexs
            read_pair_vec()  # ecg_bin_indexs
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
            "ecgSamplingRate": MAT_SR,
            "ppgSamplingRate": MAT_SR,
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

        eRate = get_field("ecgSamplingRate")
        r["ecgSamplingRate"] = float(eRate[0]) if len(eRate) > 0 else MAT_SR
        pRate = get_field("ppgSamplingRate")
        r["ppgSamplingRate"] = float(pRate[0]) if len(pRate) > 0 else MAT_SR
        r["ecgRIndex"] = get_field("ecgRIndex").astype(np.intp)
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


def rr_stats(indices, sr):
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
    rr_ms = np.diff(idx) / sr * 1000.0
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


def rr_ssd(bin_idx, bin_sr, mat_idx, mat_sr):
    rr_s = np.diff(np.array(bin_idx, dtype=np.float64)) / bin_sr
    m_s = np.diff(np.array(mat_idx, dtype=np.float64)) / mat_sr
    if len(rr_s) == 0 or len(m_s) == 0:
        return np.nan
    n = min(len(rr_s), len(m_s))
    return float(np.sum((rr_s[:n] - m_s[:n]) ** 2))


def compare_file(file_id, bin_data, mat_data):
    total = min(len(bin_data), len(mat_data))
    bins = []
    for i in range(total):
        b = bin_data[i]
        m = mat_data[i]
        mat_ecg_sr = m["ecgSamplingRate"] if m["ecgSamplingRate"] > 0 else MAT_SR
        mat_ppg_sr = m["ppgSamplingRate"] if m["ppgSamplingRate"] > 0 else MAT_SR

        b_rr = rr_stats(b["ecgRIndex"], BIN_SR)
        m_rr = rr_stats(m["ecgRIndex"], mat_ecg_sr)
        b_ppg = rr_stats(b["ppgMinAmps"], BIN_SR)
        m_ppg = rr_stats(m["ppgMinAmps"], mat_ppg_sr)

        b_ecg_amp = amp_stats(b["ecgRIndex"], b["ecgSignal"])
        m_ecg_amp = amp_stats(m["ecgRIndex"], m["ecgSignal"])
        b_ppg_amp = amp_stats(b["ppgMinAmps"], b["ppgSignal"])
        m_ppg_amp = amp_stats(m["ppgMinAmps"], m["ppgSignal"])

        # Check if RR intervals are identical within quantization limits.
        # C++ is at 1000 Hz, MATLAB at 256 Hz — intervals can differ by up to
        # half a sample at the lower rate.  Use 1 sample at the lower rate as
        # the tolerance (≈3.9 ms for 256 Hz).
        b_idx = np.array(b["ecgRIndex"], dtype=np.float64)
        m_idx = np.array(m["ecgRIndex"], dtype=np.float64)
        b_rr_sec = np.diff(b_idx) / BIN_SR
        m_rr_sec = np.diff(m_idx) / mat_ecg_sr
        if len(b_rr_sec) == len(m_rr_sec) and len(b_rr_sec) > 0:
            tol = 1.0 / min(BIN_SR, mat_ecg_sr)  # 1 sample at lower rate
            rr_identical = np.allclose(b_rr_sec, m_rr_sec, atol=tol, rtol=0)
        else:
            rr_identical = False

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
                "rr_ssd": rr_ssd(b["ecgRIndex"], BIN_SR, m["ecgRIndex"], mat_ecg_sr),
                "rr_identical": rr_identical,
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


def fmt(v, d=4):
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

    abs_diff_ecg = sum(abs(b["n_ecg_cpp"] - b["n_ecg_mat"]) for b in all_bins)
    abs_diff_ppg = sum(abs(b["n_ppg_cpp"] - b["n_ppg_mat"]) for b in all_bins)

    rows = [
        ["Metric"] + cols,
        make_row("N Subjects", [n_subjects] * 4),
        make_row("N Beats Total", n_beats),
        make_row(
            "Abs R Peak Diff", [abs_diff_ecg, abs_diff_ecg, abs_diff_ppg, abs_diff_ppg]
        ),
        make_row(
            "Mean N Beats Per Patient",
            [
                nanmean(collect(lambda b: b["n_ecg_cpp"])),
                nanmean(collect(lambda b: b["n_ecg_mat"])),
                nanmean(collect(lambda b: b["n_ppg_cpp"])),
                nanmean(collect(lambda b: b["n_ppg_mat"])),
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
                "R Peak Diff",
                "MATLAB PPG Peak N",
                "C++ PPG Peak N",
                "PPG Diff",
            ]
        )
        for fr in results:
            mat_ecg = sum(b["n_ecg_mat"] for b in fr["bins"])
            cpp_ecg = sum(b["n_ecg_cpp"] for b in fr["bins"])
            mat_ppg = sum(b["n_ppg_mat"] for b in fr["bins"])
            cpp_ppg = sum(b["n_ppg_cpp"] for b in fr["bins"])
            w.writerow(
                [
                    fr["id"],
                    mat_ecg,
                    cpp_ecg,
                    cpp_ecg - mat_ecg,
                    mat_ppg,
                    cpp_ppg,
                    cpp_ppg - mat_ppg,
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
                ]
            )


# ============================================================================
# Per-file SSD Histograms (one small chart per subject, 4 columns for PPT)
# ============================================================================


def write_allfiles_ssd_histograms(path, results, num_buckets=50):
    """
    Generate a single SVG with one small SSD histogram per subject,
    laid out 4 across to fit on a PowerPoint slide.
    """
    if not results:
        return

    num_files = len(results)
    n_cols = 5
    n_rows = (num_files + n_cols - 1) // n_cols

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

    for f_idx, fr in enumerate(results):
        col = f_idx % n_cols
        row = f_idx // n_cols
        ox = col * cell_w
        oy = row * cell_h + 50

        file_id = fr["id"]
        all_bins = fr["bins"]

        ssd_vals = [b["rr_ssd"] for b in all_bins if not np.isnan(b["rr_ssd"])]

        # "Zero error" = bins where all RR intervals are identical
        zero_count = sum(1 for b in all_bins if b["rr_identical"])
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

        lines.append(
            f'<text x="{ox + cell_w // 2}" y="{oy + 26}" '
            f'text-anchor="middle" class="stats">'
            f"n={len(ssd_vals)}  mean={mn:.4f}  med={md:.4f}</text>"
        )

        min_val = float(arr.min())
        max_val = float(arr.max())
        if max_val - min_val < 1e-15:
            max_val = min_val + 1.0
        bucket_width = (max_val - min_val) / num_buckets

        counts = np.zeros(num_buckets, dtype=int)
        for v in ssd_vals:
            bi = int((v - min_val) / bucket_width)
            if bi >= num_buckets:
                bi = num_buckets - 1
            counts[bi] += 1

        max_count = int(counts.max()) if counts.max() > 0 else 1

        sorted_counts = np.sort(counts)
        y_axis_max = max_count
        for k in range(len(sorted_counts) - 2, -1, -1):
            if sorted_counts[k] < max_count and sorted_counts[k] > 0:
                y_axis_max = int(np.ceil(sorted_counts[k] * 1.4))
                break
        if y_axis_max <= 0:
            y_axis_max = max_count

        cx = ox + pad_l
        cy = oy + pad_t
        bar_w = chart_w / num_buckets

        # Y axis
        lines.append(
            f'<line x1="{cx}" y1="{cy}" x2="{cx}" y2="{cy + chart_h}" '
            f'stroke="#333" stroke-width="1"/>'
        )
        # X axis
        lines.append(
            f'<line x1="{cx}" y1="{cy + chart_h}" '
            f'x2="{cx + chart_w}" y2="{cy + chart_h}" '
            f'stroke="#333" stroke-width="1"/>'
        )

        # Y-axis label
        lines.append(
            f'<text x="{ox + 18}" y="{cy + chart_h // 2}" '
            f'text-anchor="middle" class="y-axis-label" '
            f'transform="rotate(-90 {ox + 10},{cy + chart_h // 2})">'
            f"# bins</text>"
        )

        # Y ticks
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

        # Bars
        for b_idx in range(num_buckets):
            c = counts[b_idx]
            if c == 0:
                continue

            clipped = c > y_axis_max
            display_h = min(c, y_axis_max)
            bar_h = (display_h / y_axis_max) * chart_h
            bx = cx + b_idx * bar_w
            by = cy + chart_h - bar_h

            range_start = min_val + b_idx * bucket_width
            range_end = min_val + (b_idx + 1) * bucket_width
            tip = f"[{range_start:.4f}, {range_end:.4f}): {c}"

            lines.append(
                f'<rect x="{bx + 0.5:.1f}" y="{by:.1f}" '
                f'width="{bar_w - 1:.1f}" height="{bar_h:.1f}" '
                f'fill="#3498db" fill-opacity="0.8">'
                f"<title>{tip}</title></rect>"
            )

            if clipped:
                lines.append(
                    f'<text x="{bx + bar_w / 2:.1f}" y="{cy - 2}" '
                    f'text-anchor="middle" class="clipped">{c}</text>'
                )

        # X ticks
        for t in range(3):
            val = min_val + (max_val - min_val) * t / 2.0
            x_pos = cx + chart_w * t / 2.0
            lines.append(
                f'<text x="{x_pos:.1f}" y="{cy + chart_h + 10}" '
                f'text-anchor="middle" class="x-tick">{val:.4f}</text>'
            )

        # X-axis label
        lines.append(
            f'<text x="{cx + chart_w // 2}" y="{cy + chart_h + 22}" '
            f'text-anchor="middle" class="axis-label">ECG SSD (s)</text>'
        )

    lines.append("</svg>")

    with open(path, "w") as fout:
        fout.write("\n".join(lines))


# ============================================================================
# Combined SSD Histogram (all files pooled into one large chart)
# ============================================================================


def write_combined_ssd_histogram(path, results, num_buckets=80):
    """
    Generate a single large histogram SVG pooling all per-bin SSD values
    across every subject.
    """
    all_ssd = []
    for fr in results:
        for b in fr["bins"]:
            if not np.isnan(b["rr_ssd"]):
                all_ssd.append(b["rr_ssd"])

    if not all_ssd:
        return

    arr = np.array(all_ssd)
    mn = float(np.mean(arr))
    md = float(np.median(arr))
    q1 = float(np.percentile(arr, 25))
    q3 = float(np.percentile(arr, 75))
    zero_count = 0
    total_bins_count = 0
    for fr in results:
        for b in fr["bins"]:
            total_bins_count += 1
            if b["rr_identical"]:
                zero_count += 1
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

    # Clip y-axis like per-file version
    sorted_counts = np.sort(counts)
    y_axis_max = max_count
    for k in range(len(sorted_counts) - 2, -1, -1):
        if sorted_counts[k] < max_count and sorted_counts[k] > 0:
            y_axis_max = int(np.ceil(sorted_counts[k] * 1.2))
            break
    if y_axis_max <= 0:
        y_axis_max = max_count

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
        "  .clipped { font-size: 11px; fill: #e74c3c; font-weight: bold; }\n"
        "</style>"
    )
    lines.append('<rect width="100%" height="100%" fill="#fcfcfc"/>')

    # Title
    lines.append(
        f'<text x="{svg_w // 2}" y="24" text-anchor="middle" class="title">'
        f"Combined ECG SSD Histogram - All {len(results)} Subjects "
        f"({len(all_ssd)} bins, {zero_pct:.1f}% zero err)</text>"
    )

    # Stats line
    lines.append(
        f'<text x="{svg_w // 2}" y="44" text-anchor="middle" class="stats">'
        f"mean={mn:.6f}  median={md:.6f}  "
        f"IQR=[{q1:.6f}, {q3:.6f}]</text>"
    )

    # Y axis
    lines.append(
        f'<line x1="{cx}" y1="{cy}" x2="{cx}" y2="{cy + chart_h}" '
        f'stroke="#333" stroke-width="1"/>'
    )
    # X axis
    lines.append(
        f'<line x1="{cx}" y1="{cy + chart_h}" '
        f'x2="{cx + chart_w}" y2="{cy + chart_h}" '
        f'stroke="#333" stroke-width="1"/>'
    )

    # Y-axis label
    lines.append(
        f'<text x="18" y="{cy + chart_h // 2}" '
        f'text-anchor="middle" class="y-axis-label" '
        f'transform="rotate(-90 18,{cy + chart_h // 2})">'
        f"# bins</text>"
    )

    # Y ticks
    n_y_ticks = 5
    for t in range(n_y_ticks + 1):
        y_val = int(round(y_axis_max * t / n_y_ticks))
        y_pos = cy + chart_h - (y_val / y_axis_max) * chart_h
        lines.append(
            f'<text x="{cx - 5}" y="{y_pos + 4}" '
            f'text-anchor="end" class="tick">{y_val}</text>'
        )
        if t > 0:
            lines.append(
                f'<line x1="{cx + 1}" y1="{y_pos}" '
                f'x2="{cx + chart_w}" y2="{y_pos}" '
                f'stroke="#eee" stroke-width="0.5"/>'
            )

    # Bars
    for b_idx in range(num_buckets):
        c = counts[b_idx]
        if c == 0:
            continue

        clipped = c > y_axis_max
        display_h = min(c, y_axis_max)
        bar_h = (display_h / y_axis_max) * chart_h
        bx = cx + b_idx * bar_w
        by = cy + chart_h - bar_h

        range_start = min_val + b_idx * bucket_width
        range_end = min_val + (b_idx + 1) * bucket_width
        tip = f"[{range_start:.4f}, {range_end:.4f}): {c}"

        lines.append(
            f'<rect x="{bx + 0.5:.1f}" y="{by:.1f}" '
            f'width="{bar_w - 1:.1f}" height="{bar_h:.1f}" '
            f'fill="#3498db" fill-opacity="0.85">'
            f"<title>{tip}</title></rect>"
        )

        if clipped:
            lines.append(
                f'<text x="{bx + bar_w / 2:.1f}" y="{cy - 4}" '
                f'text-anchor="middle" class="clipped">{c}</text>'
            )

    # X ticks
    n_x_ticks = 5
    for t in range(n_x_ticks + 1):
        val = min_val + (max_val - min_val) * t / n_x_ticks
        x_pos = cx + chart_w * t / n_x_ticks
        lines.append(
            f'<text x="{x_pos:.1f}" y="{cy + chart_h + 16}" '
            f'text-anchor="middle" class="tick">{val:.4f}</text>'
        )

    # X-axis label
    lines.append(
        f'<text x="{cx + chart_w // 2}" y="{cy + chart_h + 38}" '
        f'text-anchor="middle" class="axis-label">ECG SSD (s)</text>'
    )

    lines.append("</svg>")

    with open(path, "w") as fout:
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
    write_combined_ssd_histogram(OUTPUT_DIR / "combined_ssd_histogram.svg", results)

    print(f"\nSummary CSV:       {OUTPUT_DIR / 'summary.csv'}")
    print(f"File Summary CSV:  {OUTPUT_DIR / 'file_summary.csv'}")
    print(f"SSD Per-File:      {OUTPUT_DIR / 'allfiles_ssd_histograms.svg'}")
    print(f"SSD Combined:      {OUTPUT_DIR / 'combined_ssd_histogram.svg'}")
    print(f"Per-file CSVs:     {OUTPUT_DIR}/<subject>_wave_comparison.csv")


if __name__ == "__main__":
    main()
