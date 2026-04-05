"""
compare_templates.py

Compares MATLAB (256 Hz) vs C++ (2000 Hz) template outputs.
Only compares raw ECG method and PPG.
Outputs per-subject CSVs and a summary CSV with mean, median, and IQR.
Also outputs filtered versions that exclude bins where MATLAB and C++
have different R-peak counts.
"""

import csv
import os
import struct
import sys
import traceback
from pathlib import Path

import numpy as np
import scipy.io as sio
from scipy import signal as scipy_signal

MATLAB_DIR = (
    "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\5_generate_template_files\\matlab"
)
CPP_DIR = (
    "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\5_generate_template_files\\mesa"
)
MATLAB_WAVE_DIR = (
    "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\4_wave_bound_files\\matlab"
)
CPP_WAVE_DIR = (
    "D:\\USERS\\MiraWelner\\QTVI\\QTVI-data-files\\4_wave_bound_files\\mesa_files"
)
OUTPUT_DIR = "D:\\USERS\\MiraWelner\\QTVI\\testing\\5_template_generation\\results"
MATLAB_SR = 256.0
CPP_SR = 2000.0
SR_RATIO = CPP_SR / MATLAB_SR


# =============================================================================
# Read C++ template_info .bin
# =============================================================================
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
        for bi in range(num_bins):
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


# =============================================================================
# Read C++ wave_markings .bin  (only ch1.raw R-peak count per bin)
# =============================================================================
def read_cpp_wave_rpeak_counts(path):
    counts = []
    MAX_SANE = 50000000

    with open(path, "rb") as f:

        def ru64():
            d = f.read(8)
            if len(d) < 8:
                return None
            return struct.unpack("<Q", d)[0]

        def skip_idx_array():
            sz = ru64()
            if sz is None:
                return 0
            if sz > 0 and sz <= MAX_SANE:
                f.seek(sz * 8, 1)
            return sz

        def skip_signal_array():
            sz = ru64()
            if sz is None:
                return
            if sz > 0 and sz <= MAX_SANE:
                f.seek(sz * 8, 1)

        def skip_pair_vec():
            sz = ru64()
            if sz is None:
                return
            if sz > 0 and sz <= MAX_SANE:
                f.seek(sz * 16, 1)

        num_bins = ru64()
        if num_bins is None:
            return counts

        for i in range(num_bins):
            ch1_raw_count = skip_idx_array()  # ch1.raw
            skip_idx_array()  # ch1.squared
            skip_idx_array()  # ch1.absval
            skip_idx_array()  # ch2.raw
            skip_idx_array()  # ch2.squared
            skip_idx_array()  # ch2.absval
            skip_idx_array()  # ch3.raw
            skip_idx_array()  # ch3.squared
            skip_idx_array()  # ch3.absval
            skip_idx_array()  # ppgMaxAmps
            skip_idx_array()  # ppgMinAmps
            skip_signal_array()  # ppgSignal
            skip_signal_array()  # ecgSignal
            skip_signal_array()  # ecgSignal2
            skip_signal_array()  # ecgSignal3
            skip_signal_array()  # ch1.squared_signal
            skip_signal_array()  # ch1.absval_signal
            skip_signal_array()  # ch2.squared_signal
            skip_signal_array()  # ch2.absval_signal
            skip_signal_array()  # ch3.squared_signal
            skip_signal_array()  # ch3.absval_signal
            f.read(9)  # 9 noise flags
            num_pairs = ru64()
            if num_pairs is not None and num_pairs > 0 and num_pairs <= MAX_SANE:
                f.seek(num_pairs * 16, 1)
            skip_pair_vec()  # ppg_bin_indexs
            skip_pair_vec()  # ecg_bin_indexs
            counts.append(ch1_raw_count)

    return counts


# =============================================================================
# Read MATLAB wave_data .mat  (ecgRIndex count per bin)
# =============================================================================
def read_matlab_wave_rpeak_counts(path):
    mat = sio.loadmat(path, squeeze_me=False)
    if "wave_data" not in mat:
        raise ValueError(f"No 'wave_data' in {path}")
    raw = mat["wave_data"].flatten()
    counts = []
    for i in range(len(raw)):
        cell = raw[i]
        while isinstance(cell, np.ndarray) and cell.dtype == object and cell.ndim == 0:
            cell = cell.flat[0]
        if cell is None or (isinstance(cell, np.ndarray) and cell.size == 0):
            counts.append(0)
            continue
        try:
            r_idx = cell["ecgRIndex"]
            while isinstance(r_idx, np.ndarray) and r_idx.dtype == object:
                if r_idx.size == 1:
                    r_idx = r_idx.flat[0]
                else:
                    break
            if isinstance(r_idx, np.ndarray):
                flat = r_idx.flatten()
                count = int(np.sum(~np.isnan(flat.astype(np.float64))))
                counts.append(count)
            else:
                counts.append(0)
        except (ValueError, KeyError):
            counts.append(0)
    return counts


# =============================================================================
# Read MATLAB template_info .mat
# =============================================================================
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
            except:
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


# =============================================================================
# Correlation helpers
# =============================================================================
def resample_corr(mat_vec, cpp_vec):
    if mat_vec.size == 0 or cpp_vec.size == 0:
        return np.nan
    cpp_resampled = scipy_signal.resample(cpp_vec, mat_vec.size)
    mat_r = np.ptp(mat_vec)
    cpp_r = np.ptp(cpp_resampled)
    if mat_r < 1e-10 or cpp_r < 1e-10:
        return np.nan
    mat_n = (mat_vec - np.min(mat_vec)) / mat_r
    cpp_n = (cpp_resampled - np.min(cpp_resampled)) / cpp_r
    if np.std(mat_n) < 1e-10 or np.std(cpp_n) < 1e-10:
        return np.nan
    return np.corrcoef(mat_n, cpp_n)[0, 1]


def compute_bin_correlations(mat_list, cpp_list):
    n = min(len(mat_list), len(cpp_list))
    rows = []
    for i in range(n):
        mat = mat_list[i]
        cpp = cpp_list[i]
        row = {"bin": i}
        if mat is None or cpp is None:
            row["ecg_raw"] = np.nan
            row["ppg"] = np.nan
            rows.append(row)
            continue
        row["ecg_raw"] = resample_corr(mat["ecgTemplate"], cpp["ch1_ecgTemplate_raw"])
        row["ppg"] = resample_corr(mat["ppgTemplate"], cpp["ppgTemplate"])
        rows.append(row)
    return rows


# =============================================================================
# Stats
# =============================================================================
def compute_stats(vals):
    clean = [v for v in vals if not np.isnan(v)]
    if not clean:
        return np.nan, np.nan, np.nan, np.nan
    arr = np.array(clean)
    q1 = np.percentile(arr, 25)
    q3 = np.percentile(arr, 75)
    return np.mean(arr), np.median(arr), q1, q3


def compute_summary(rows):
    stats = {}
    for f in ["ecg_raw", "ppg"]:
        vals = [r[f] for r in rows]
        mean_val, median_val, q1_val, q3_val = compute_stats(vals)
        stats[f"{f}_mean"] = mean_val
        stats[f"{f}_median"] = median_val
        stats[f"{f}_q1"] = q1_val
        stats[f"{f}_q3"] = q3_val
    return stats


def fmt(val):
    return f"{val:.3f}" if not np.isnan(val) else ""


def fmt_iqr(q1, q3):
    if np.isnan(q1) or np.isnan(q3):
        return ""
    return f"({q1:.3f}, {q3:.3f})"


# =============================================================================
# CSV writers
# =============================================================================
def write_subject_csv(path, rows, mat_rpeak_counts=None, cpp_rpeak_counts=None):
    has_beats = mat_rpeak_counts is not None and cpp_rpeak_counts is not None
    fields = ["bin"]
    if has_beats:
        fields += ["matlab_beats", "cpp_beats"]
    fields += ["ecg_raw", "ppg"]
    value_fields = ["ecg_raw", "ppg"]
    stats = compute_summary(rows)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        mean_row = ["MEAN"]
        if has_beats:
            mean_row += ["", ""]
        mean_row += [fmt(stats["ecg_raw_mean"]), fmt(stats["ppg_mean"])]
        writer.writerow(mean_row)
        median_row = ["MEDIAN"]
        if has_beats:
            median_row += ["", ""]
        median_row += [fmt(stats["ecg_raw_median"]), fmt(stats["ppg_median"])]
        writer.writerow(median_row)
        iqr_row = ["IQR"]
        if has_beats:
            iqr_row += ["", ""]
        iqr_row += [
            fmt_iqr(stats["ecg_raw_q1"], stats["ecg_raw_q3"]),
            fmt_iqr(stats["ppg_q1"], stats["ppg_q3"]),
        ]
        writer.writerow(iqr_row)
        for row in rows:
            i = row["bin"]
            csv_row = [i]
            if has_beats:
                mat_b = mat_rpeak_counts[i] if i < len(mat_rpeak_counts) else ""
                cpp_b = cpp_rpeak_counts[i] if i < len(cpp_rpeak_counts) else ""
                csv_row += [mat_b, cpp_b]
            for field in value_fields:
                csv_row.append(fmt(row[field]))
            writer.writerow(csv_row)


def write_summary_csv(path, summary_rows, include_beats=False):
    fields = ["Subject", "N Bins"]
    if include_beats:
        fields += ["MATLAB Beats", "C++ Beats"]
    fields += [
        "Mean ECG PCC",
        "Median ECG PCC",
        "IQR ECG PCC",
        "Mean PPG PCC",
        "Median PPG PCC",
        "IQR PPG PCC",
    ]

    grand = {}
    for f in [
        "ecg_raw_mean",
        "ecg_raw_median",
        "ecg_raw_q1",
        "ecg_raw_q3",
        "ppg_mean",
        "ppg_median",
        "ppg_q1",
        "ppg_q3",
    ]:
        vals = [r[f] for r in summary_rows if not np.isnan(r[f])]
        grand[f] = np.mean(vals) if vals else np.nan

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        avg_row = ["AVERAGE", ""]
        if include_beats:
            total_mat = sum(r.get("mat_beats", 0) for r in summary_rows)
            total_cpp = sum(r.get("cpp_beats", 0) for r in summary_rows)
            avg_row += [total_mat, total_cpp]
        avg_row += [
            fmt(grand["ecg_raw_mean"]),
            fmt(grand["ecg_raw_median"]),
            fmt_iqr(grand["ecg_raw_q1"], grand["ecg_raw_q3"]),
            fmt(grand["ppg_mean"]),
            fmt(grand["ppg_median"]),
            fmt_iqr(grand["ppg_q1"], grand["ppg_q3"]),
        ]
        writer.writerow(avg_row)
        for row in summary_rows:
            csv_row = [row["subject"], row["n_bins"]]
            if include_beats:
                csv_row += [row.get("mat_beats", 0), row.get("cpp_beats", 0)]
            csv_row += [
                fmt(row["ecg_raw_mean"]),
                fmt(row["ecg_raw_median"]),
                fmt_iqr(row["ecg_raw_q1"], row["ecg_raw_q3"]),
                fmt(row["ppg_mean"]),
                fmt(row["ppg_median"]),
                fmt_iqr(row["ppg_q1"], row["ppg_q3"]),
            ]
            writer.writerow(csv_row)


# =============================================================================
# File matching
# =============================================================================
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
    mat_only = sorted(set(mat_files.keys()) - set(cpp_files.keys()))
    cpp_only = sorted(set(cpp_files.keys()) - set(mat_files.keys()))
    return common, mat_files, cpp_files, mat_only, cpp_only


def find_wave_files(mat_dir, cpp_dir, subject_id):
    mat_file = None
    for f in Path(mat_dir).glob(f"*{subject_id}*_wave_data.mat"):
        mat_file = f
        break
    if mat_file is None:
        for f in Path(mat_dir).glob(f"*{subject_id}*.mat"):
            mat_file = f
            break

    cpp_file = None
    for f in Path(cpp_dir).glob(f"*{subject_id}*_wave_markings.bin"):
        cpp_file = f
        break
    if cpp_file is None:
        for f in Path(cpp_dir).glob(f"*{subject_id}*.bin"):
            cpp_file = f
            break

    return mat_file, cpp_file


# =============================================================================
# Main
# =============================================================================
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    common, mat_files, cpp_files, mat_only, cpp_only = find_matching_files()

    print("=" * 70)
    print(f"TEMPLATE COMPARISON: MATLAB ({MATLAB_SR:.0f} Hz) vs C++ ({CPP_SR:.0f} Hz)")
    print("=" * 70)
    print(f"MATLAB template dir:  {MATLAB_DIR}")
    print(f"C++    template dir:  {CPP_DIR}")
    print(f"MATLAB wave dir:      {MATLAB_WAVE_DIR}")
    print(f"C++    wave dir:      {CPP_WAVE_DIR}")
    print(f"Output dir:           {OUTPUT_DIR}")
    print(f"Matched:              {len(common)}")
    if mat_only:
        print(f"MATLAB only: {', '.join(mat_only)}")
    if cpp_only:
        print(f"C++ only:    {', '.join(cpp_only)}")
    print()

    summary_rows_all = []
    summary_rows_filtered = []

    for sid in common:
        print(f"  {sid}...", end="", flush=True)

        # ---- Step 1: Read templates ----
        try:
            mat_templates = read_matlab_template_mat(str(mat_files[sid]))
        except Exception:
            print(f" MATLAB TEMPLATE READ ERROR:")
            traceback.print_exc()
            continue

        try:
            cpp_templates = read_cpp_template_bin(str(cpp_files[sid]))
        except Exception:
            print(f" CPP TEMPLATE READ ERROR:")
            traceback.print_exc()
            continue

        # ---- Step 2: Compute all-bin correlations ----
        try:
            rows_all = compute_bin_correlations(mat_templates, cpp_templates)
        except Exception:
            print(f" COMPARISON ERROR:")
            traceback.print_exc()
            continue

        # ---- Step 3: Write unfiltered per-subject CSV (written later after rpeak counts are loaded) ----
        stats_all = compute_summary(rows_all)

        # ---- Step 4: Read R-peak counts ----
        mat_wave_file, cpp_wave_file = find_wave_files(
            MATLAB_WAVE_DIR, CPP_WAVE_DIR, sid
        )

        mat_rpeak_counts = None
        cpp_rpeak_counts = None
        rows_filtered = None
        n_excluded = 0

        if mat_wave_file and cpp_wave_file:
            try:
                mat_rpeak_counts = read_matlab_wave_rpeak_counts(str(mat_wave_file))
                cpp_rpeak_counts = read_cpp_wave_rpeak_counts(str(cpp_wave_file))

                n_bins = min(
                    len(rows_all), len(mat_rpeak_counts), len(cpp_rpeak_counts)
                )

                rows_filtered = []
                n_excluded = 0
                for idx in range(n_bins):
                    if mat_rpeak_counts[idx] == cpp_rpeak_counts[idx]:
                        rows_filtered.append(rows_all[idx])
                    else:
                        n_excluded += 1

            except Exception:
                print(f" RPEAK COUNT ERROR:")
                traceback.print_exc()
                rows_filtered = None
        else:
            missing = []
            if not mat_wave_file:
                missing.append("MATLAB")
            if not cpp_wave_file:
                missing.append("C++")
            print(f" (missing wave files: {', '.join(missing)})", end="")

        # ---- Step 5: Compute beat totals, write unfiltered CSV, store summary ----
        total_mat_beats = 0
        total_cpp_beats = 0
        if mat_rpeak_counts is not None:
            total_mat_beats = sum(mat_rpeak_counts)
        if cpp_rpeak_counts is not None:
            total_cpp_beats = sum(cpp_rpeak_counts)

        subject_csv_all = os.path.join(OUTPUT_DIR, f"{sid}_correlations.csv")
        write_subject_csv(subject_csv_all, rows_all, mat_rpeak_counts, cpp_rpeak_counts)

        summary_rows_all.append(
            {
                "subject": sid,
                "n_bins": len(rows_all),
                "mat_beats": total_mat_beats,
                "cpp_beats": total_cpp_beats,
                **stats_all,
            }
        )

        # ---- Step 6: Write filtered per-subject CSV ----
        if rows_filtered is not None and len(rows_filtered) > 0:
            subject_csv_filt = os.path.join(
                OUTPUT_DIR, f"{sid}_correlations_matched_peaks.csv"
            )
            write_subject_csv(
                subject_csv_filt, rows_filtered, mat_rpeak_counts, cpp_rpeak_counts
            )

            stats_filt = compute_summary(rows_filtered)
            summary_rows_filtered.append(
                {"subject": sid, "n_bins": len(rows_filtered), **stats_filt}
            )

            print(
                f" {len(rows_all)} bins "
                f"(matched: {len(rows_filtered)}, excluded: {n_excluded})"
            )
            print(
                f"    ALL:  ECG[mean={stats_all['ecg_raw_mean']:.3f} "
                f"med={stats_all['ecg_raw_median']:.3f} "
                f"IQR={fmt_iqr(stats_all['ecg_raw_q1'], stats_all['ecg_raw_q3'])}] "
                f"PPG[mean={stats_all['ppg_mean']:.3f} "
                f"med={stats_all['ppg_median']:.3f} "
                f"IQR={fmt_iqr(stats_all['ppg_q1'], stats_all['ppg_q3'])}]"
            )
            print(
                f"    FILT: ECG[mean={stats_filt['ecg_raw_mean']:.3f} "
                f"med={stats_filt['ecg_raw_median']:.3f} "
                f"IQR={fmt_iqr(stats_filt['ecg_raw_q1'], stats_filt['ecg_raw_q3'])}] "
                f"PPG[mean={stats_filt['ppg_mean']:.3f} "
                f"med={stats_filt['ppg_median']:.3f} "
                f"IQR={fmt_iqr(stats_filt['ppg_q1'], stats_filt['ppg_q3'])}]"
            )
        else:
            if rows_filtered is not None:
                print(f" {len(rows_all)} bins (all excluded by peak filter)")
            else:
                print(f" {len(rows_all)} bins (no wave files for peak filtering)")

    # --- Write summaries ---
    summary_csv_all = os.path.join(OUTPUT_DIR, "summary.csv")
    write_summary_csv(summary_csv_all, summary_rows_all, include_beats=True)

    summary_csv_filt = os.path.join(OUTPUT_DIR, "summary_matched_peaks.csv")
    if summary_rows_filtered:
        write_summary_csv(summary_csv_filt, summary_rows_filtered, include_beats=False)

    # --- Print summary tables ---
    for label, s_rows, csv_path, show_beats in [
        ("ALL BINS", summary_rows_all, summary_csv_all, True),
        ("MATCHED R-PEAK BINS ONLY", summary_rows_filtered, summary_csv_filt, False),
    ]:
        if not s_rows:
            continue
        print(f"\n{'=' * 70}")
        print(f"SUMMARY: {label}")
        print(f"{'=' * 70}")
        hdr = f"{'Subject':<25s} {'Bins':>5s} "
        if show_beats:
            hdr += f"{'MAT beats':>10s} {'CPP beats':>10s} "
        hdr += (
            f"{'ECG mean':>9s} {'ECG med':>9s} {'ECG IQR':>16s} "
            f"{'PPG mean':>9s} {'PPG med':>9s} {'PPG IQR':>16s}"
        )
        print(hdr)
        sep_len = 126 if show_beats else 106
        print("-" * sep_len)

        for row in s_rows:
            line = f"{row['subject']:<25s} {row['n_bins']:>5d} "
            if show_beats:
                line += (
                    f"{row.get('mat_beats', 0):>10d} {row.get('cpp_beats', 0):>10d} "
                )
            line += (
                f"{row['ecg_raw_mean']:>9.3f} {row['ecg_raw_median']:>9.3f} "
                f"{fmt_iqr(row['ecg_raw_q1'], row['ecg_raw_q3']):>16s} "
                f"{row['ppg_mean']:>9.3f} {row['ppg_median']:>9.3f} "
                f"{fmt_iqr(row['ppg_q1'], row['ppg_q3']):>16s}"
            )
            print(line)

        print("-" * sep_len)
        grand = {}
        for f in [
            "ecg_raw_mean",
            "ecg_raw_median",
            "ecg_raw_q1",
            "ecg_raw_q3",
            "ppg_mean",
            "ppg_median",
            "ppg_q1",
            "ppg_q3",
        ]:
            vals = [r[f] for r in s_rows if not np.isnan(r[f])]
            grand[f] = np.mean(vals) if vals else np.nan
        avg_line = f"{'AVERAGE':<25s} {'':>5s} "
        if show_beats:
            total_mat = sum(r.get("mat_beats", 0) for r in s_rows)
            total_cpp = sum(r.get("cpp_beats", 0) for r in s_rows)
            avg_line += f"{total_mat:>10d} {total_cpp:>10d} "
        avg_line += (
            f"{grand['ecg_raw_mean']:>9.3f} {grand['ecg_raw_median']:>9.3f} "
            f"{fmt_iqr(grand['ecg_raw_q1'], grand['ecg_raw_q3']):>16s} "
            f"{grand['ppg_mean']:>9.3f} {grand['ppg_median']:>9.3f} "
            f"{fmt_iqr(grand['ppg_q1'], grand['ppg_q3']):>16s}"
        )
        print(avg_line)
        print(f"CSV: {csv_path}")

    print(f"\nPer-subject CSVs:          {OUTPUT_DIR}/<subject>_correlations.csv")
    print(
        f"Per-subject filtered CSVs: {OUTPUT_DIR}/<subject>_correlations_matched_peaks.csv"
    )
    print(f"Summary CSV:               {summary_csv_all}")
    if summary_rows_filtered:
        print(f"Filtered summary CSV:      {summary_csv_filt}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
