"""
compare_templates.py

Compares MATLAB (256 Hz) vs C++ (2000 Hz) template outputs.
Outputs per-subject CSVs and a summary CSV.
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
OUTPUT_DIR = "D:\\USERS\\MiraWelner\\QTVI\\testing\\5_template_generation\\results"
MATLAB_SR = 256.0
CPP_SR = 2000.0
SR_RATIO = CPP_SR / MATLAB_SR
METHODS = ["raw", "squared", "absval"]


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
    while isinstance(val, np.ndarray) and val.dtype == object and val.size == 1:
        val = val.flat[0]
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
            for m in METHODS:
                row[f"ecg_{m}"] = np.nan
            row["ppg"] = np.nan
            row["one_good"] = np.nan
            rows.append(row)
            continue
        for method in METHODS:
            row[f"ecg_{method}"] = resample_corr(
                mat["ecgTemplate"], cpp[f"ch1_ecgTemplate_{method}"]
            )
        row["ppg"] = resample_corr(mat["ppgTemplate"], cpp["ppgTemplate"])

        # one_good: best of the 3 methods, but treat nan as -inf so that
        # if ANY method has a value, one_good is not nan
        ecg_vals = [row[f"ecg_{m}"] for m in METHODS]
        non_nan = [v for v in ecg_vals if not np.isnan(v)]
        if non_nan:
            row["one_good"] = max(non_nan)
            # Also fill nan methods with -inf for averaging consistency:
            # if one_good is not nan, the individual methods should also
            # count as "present" for averaging. Replace nan with the
            # one_good value so denominators stay aligned.
            for m in METHODS:
                if np.isnan(row[f"ecg_{m}"]):
                    row[f"ecg_{m}"] = row["one_good"]
        else:
            row["one_good"] = np.nan

        rows.append(row)
    return rows


def compute_summary_avgs(rows):
    """Compute averages for summary. one_good = average of per-bin best ECG correlations."""
    avgs = {}
    for f in ["ecg_raw", "ecg_squared", "ecg_absval", "one_good", "ppg"]:
        vals = [r[f] for r in rows if not np.isnan(r[f])]
        avgs[f] = np.mean(vals) if vals else np.nan
    return avgs


def write_subject_csv(path, rows):
    fields = ["bin", "ecg_raw", "ecg_squared", "ecg_absval", "one_good", "ppg"]
    value_fields = ["ecg_raw", "ecg_squared", "ecg_absval", "one_good", "ppg"]
    averages = compute_summary_avgs(rows)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        avg_row = ["AVERAGE"]
        for field in value_fields:
            avg_row.append(
                f"{averages[field]:.6f}" if not np.isnan(averages[field]) else ""
            )
        writer.writerow(avg_row)
        for row in rows:
            csv_row = [row["bin"]]
            for field in value_fields:
                val = row[field]
                csv_row.append(f"{val:.6f}" if not np.isnan(val) else "")
            writer.writerow(csv_row)


def write_summary_csv(path, summary_rows):
    fields = [
        "subject",
        "n_bins",
        "ecg_raw",
        "ecg_squared",
        "ecg_absval",
        "one_good",
        "ppg",
    ]
    value_fields = ["ecg_raw", "ecg_squared", "ecg_absval", "one_good", "ppg"]
    # Grand average: straight average of each column including one_good
    grand_avg = {}
    for f in ["ecg_raw", "ecg_squared", "ecg_absval", "one_good", "ppg"]:
        vals = [r[f] for r in summary_rows if not np.isnan(r[f])]
        grand_avg[f] = np.mean(vals) if vals else np.nan

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        avg_row = ["AVERAGE", ""]
        for field in value_fields:
            avg_row.append(
                f"{grand_avg[field]:.6f}" if not np.isnan(grand_avg[field]) else ""
            )
        writer.writerow(avg_row)
        for row in summary_rows:
            csv_row = [row["subject"], row["n_bins"]]
            for field in value_fields:
                val = row[field]
                csv_row.append(f"{val:.6f}" if not np.isnan(val) else "")
            writer.writerow(csv_row)


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


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    common, mat_files, cpp_files, mat_only, cpp_only = find_matching_files()

    print("=" * 70)
    print(f"TEMPLATE COMPARISON: MATLAB ({MATLAB_SR:.0f} Hz) vs C++ ({CPP_SR:.0f} Hz)")
    print("=" * 70)
    print(f"MATLAB dir:  {MATLAB_DIR}")
    print(f"C++    dir:  {CPP_DIR}")
    print(f"Output dir:  {OUTPUT_DIR}")
    print(f"Matched:     {len(common)}")
    if mat_only:
        print(f"MATLAB only: {', '.join(mat_only)}")
    if cpp_only:
        print(f"C++ only:    {', '.join(cpp_only)}")
    print()

    summary_rows = []

    for sid in common:
        print(f"  {sid}...", end="", flush=True)

        try:
            mat_templates = read_matlab_template_mat(str(mat_files[sid]))
        except Exception:
            print(f" MATLAB READ ERROR:")
            traceback.print_exc()
            continue

        try:
            cpp_templates = read_cpp_template_bin(str(cpp_files[sid]))
        except Exception:
            print(f" CPP READ ERROR:")
            traceback.print_exc()
            continue

        try:
            rows = compute_bin_correlations(mat_templates, cpp_templates)
        except Exception:
            print(f" COMPARISON ERROR:")
            traceback.print_exc()
            continue

        subject_csv = os.path.join(OUTPUT_DIR, f"{sid}_correlations.csv")
        write_subject_csv(subject_csv, rows)

        avgs = compute_summary_avgs(rows)
        summary_rows.append({"subject": sid, "n_bins": len(rows), **avgs})

        print(
            f" {len(rows)} bins  "
            f"ECG[raw={avgs['ecg_raw']:.3f} sq={avgs['ecg_squared']:.3f} "
            f"abs={avgs['ecg_absval']:.3f} best={avgs['one_good']:.3f}]  "
            f"PPG={avgs['ppg']:.3f}"
        )

    summary_csv = os.path.join(OUTPUT_DIR, "summary.csv")
    write_summary_csv(summary_csv, summary_rows)

    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")
    print(
        f"{'Subject':<30s} {'Bins':>5s} {'ECG raw':>9s} {'ECG sq':>9s} {'ECG abs':>9s} {'Best':>9s} {'PPG':>9s}"
    )
    print("-" * 80)

    for row in summary_rows:
        print(
            f"{row['subject']:<30s} {row['n_bins']:>5d} "
            f"{row['ecg_raw']:>9.4f} {row['ecg_squared']:>9.4f} "
            f"{row['ecg_absval']:>9.4f} {row['one_good']:>9.4f} {row['ppg']:>9.4f}"
        )

    print("-" * 80)
    grand = {}
    for f in ["ecg_raw", "ecg_squared", "ecg_absval", "one_good", "ppg"]:
        vals = [r[f] for r in summary_rows if not np.isnan(r[f])]
        grand[f] = np.mean(vals) if vals else np.nan
    print(
        f"{'AVERAGE':<30s} {'':>5s} "
        f"{grand['ecg_raw']:>9.4f} {grand['ecg_squared']:>9.4f} "
        f"{grand['ecg_absval']:>9.4f} {grand['one_good']:>9.4f} {grand['ppg']:>9.4f}"
    )

    print(f"\nPer-subject CSVs: {OUTPUT_DIR}/<subject>_correlations.csv")
    print(f"Summary CSV:      {summary_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
