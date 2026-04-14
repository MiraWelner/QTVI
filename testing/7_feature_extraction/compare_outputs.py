#!/usr/bin/env python3
"""
compare_outputs.py — Compare MATLAB .mat feature outputs vs C++ .bin feature outputs.

Usage:
    python compare_outputs.py

Outputs:
    comparison_report.csv  — Per-field, per-subject comparison results.
    comparison_summary.csv — Per-field summary across all subjects.
"""

import os
import sys
import struct
import glob
import re
import csv
import numpy as np
from pathlib import Path

# ──────────────────────────────────────────────────────────────────────────────
# Configuration — edit these paths
# ──────────────────────────────────────────────────────────────────────────────
MATLAB_DIR = r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\7_generate_features\mesa_features_daniel"
CPP_DIR    = r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\7_generate_features\mesa_features_mira"

REPORT_CSV  = "comparison_report.csv"
SUMMARY_CSV = "comparison_summary.csv"

# Tolerance for floating-point comparison
RTOL = 1e-6   # relative tolerance
ATOL = 1e-9   # absolute tolerance


# ──────────────────────────────────────────────────────────────────────────────
# .mat loader (requires scipy)
# ──────────────────────────────────────────────────────────────────────────────
def load_mat_features(path):
    """Load a MATLAB _feature_output.mat file. Returns dict of {field_name: np.array}."""
    try:
        import scipy.io as sio
        mat = sio.loadmat(path, squeeze_me=True)
    except NotImplementedError:
        # v7.3 .mat files require h5py
        import h5py
        mat = {}
        with h5py.File(path, 'r') as f:
            if 'beats_flattened' in f:
                group = f['beats_flattened']
                for key in group.keys():
                    val = group[key][()]
                    if isinstance(val, np.ndarray):
                        mat[key] = val.flatten().astype(np.float64)
                    else:
                        mat[key] = np.array([val], dtype=np.float64)
            else:
                for key in f.keys():
                    if key.startswith('#'):
                        continue
                    val = f[key][()]
                    if isinstance(val, np.ndarray):
                        mat[key] = val.flatten().astype(np.float64)
        return mat

    # scipy.io.loadmat: beats_flattened is typically a struct
    fields = {}
    if 'beats_flattened' in mat:
        bf = mat['beats_flattened']
        if hasattr(bf, 'dtype') and bf.dtype.names:
            # Structured array
            for name in bf.dtype.names:
                val = bf[name]
                if isinstance(val, np.ndarray):
                    fields[name] = val.flatten().astype(np.float64)
                else:
                    fields[name] = np.array(val).flatten().astype(np.float64)
        elif isinstance(bf, np.ndarray) and bf.shape == ():
            # 0-d array wrapping a void/struct
            item = bf.item()
            if hasattr(item, '_fieldnames'):
                for name in item._fieldnames:
                    val = getattr(item, name)
                    fields[name] = np.array(val).flatten().astype(np.float64)
            elif hasattr(bf.dtype, 'names') and bf.dtype.names:
                for name in bf.dtype.names:
                    val = bf[name]
                    if isinstance(val, np.ndarray):
                        fields[name] = val.flatten().astype(np.float64)
                    else:
                        fields[name] = np.array(val).flatten().astype(np.float64)
        else:
            # Try direct field access
            try:
                for name in dir(bf):
                    if name.startswith('_'):
                        continue
                    val = getattr(bf, name)
                    if isinstance(val, (np.ndarray, list)):
                        fields[name] = np.array(val).flatten().astype(np.float64)
            except Exception:
                pass
    else:
        # Fields at top level
        for key, val in mat.items():
            if key.startswith('__'):
                continue
            if isinstance(val, np.ndarray):
                fields[key] = val.flatten().astype(np.float64)

    return fields


# ──────────────────────────────────────────────────────────────────────────────
# .bin loader
# ──────────────────────────────────────────────────────────────────────────────
def load_bin_features(path):
    """Load a C++ _feature_output.bin file. Returns dict of {field_name: np.array}."""
    fields = {}
    with open(path, 'rb') as f:
        data = f.read()

    pos = 0

    def read_u64():
        nonlocal pos
        val = struct.unpack_from('<Q', data, pos)[0]
        pos += 8
        return val

    def read_doubles(n):
        nonlocal pos
        arr = np.frombuffer(data, dtype='<f8', count=n, offset=pos).copy()
        pos += n * 8
        return arr

    def read_vec():
        n = read_u64()
        return read_doubles(n)

    def read_field():
        name_len = read_u64()
        name = data[pos:pos + name_len].decode('ascii')
        nonlocal pos
        pos += name_len
        vec = read_vec()
        return name, vec

    n_beats = read_u64()
    n_fields = read_u64()

    for _ in range(n_fields):
        name, vec = read_field()
        fields[name] = vec

    # Trailing ppg_wout_noise
    if pos < len(data):
        try:
            ppg = read_vec()
            fields['ppg_wout_noise'] = ppg
        except Exception:
            pass

    fields['__nBeats__'] = np.array([n_beats])
    return fields


# ──────────────────────────────────────────────────────────────────────────────
# Comparison
# ──────────────────────────────────────────────────────────────────────────────
def compare_fields(mat_fields, bin_fields):
    """
    Compare two dicts of {name: np.array}.
    Returns list of dicts with comparison results per field.
    """
    all_keys = sorted(set(list(mat_fields.keys()) + list(bin_fields.keys())))
    results = []

    for key in all_keys:
        if key.startswith('__'):
            continue

        row = {'field': key}

        if key not in mat_fields:
            row['status'] = 'CPP_ONLY'
            row['mat_len'] = 0
            row['cpp_len'] = len(bin_fields[key])
            row['max_abs_diff'] = ''
            row['max_rel_diff'] = ''
            row['mean_abs_diff'] = ''
            row['pct_match'] = ''
            row['nan_count_mat'] = ''
            row['nan_count_cpp'] = ''
            results.append(row)
            continue

        if key not in bin_fields:
            row['status'] = 'MAT_ONLY'
            row['mat_len'] = len(mat_fields[key])
            row['cpp_len'] = 0
            row['max_abs_diff'] = ''
            row['max_rel_diff'] = ''
            row['mean_abs_diff'] = ''
            row['pct_match'] = ''
            row['nan_count_mat'] = ''
            row['nan_count_cpp'] = ''
            results.append(row)
            continue

        a = mat_fields[key]
        b = bin_fields[key]

        row['mat_len'] = len(a)
        row['cpp_len'] = len(b)
        row['nan_count_mat'] = int(np.sum(np.isnan(a)))
        row['nan_count_cpp'] = int(np.sum(np.isnan(b)))

        if len(a) != len(b):
            row['status'] = 'LENGTH_MISMATCH'
            min_len = min(len(a), len(b))
            if min_len > 0:
                aa, bb = a[:min_len], b[:min_len]
                both_nan = np.isnan(aa) & np.isnan(bb)
                either_nan = np.isnan(aa) | np.isnan(bb)
                valid = ~either_nan
                if np.sum(valid) > 0:
                    diff = np.abs(aa[valid] - bb[valid])
                    row['max_abs_diff'] = f"{np.max(diff):.10e}"
                    denom = np.maximum(np.abs(aa[valid]), np.abs(bb[valid]))
                    denom[denom == 0] = 1
                    row['max_rel_diff'] = f"{np.max(diff / denom):.10e}"
                    row['mean_abs_diff'] = f"{np.mean(diff):.10e}"
                    matching = np.sum(both_nan) + np.sum(np.isclose(aa[valid], bb[valid], rtol=RTOL, atol=ATOL))
                    row['pct_match'] = f"{100.0 * matching / min_len:.4f}"
                else:
                    row['max_abs_diff'] = 'N/A (all NaN)'
                    row['max_rel_diff'] = 'N/A'
                    row['mean_abs_diff'] = 'N/A'
                    row['pct_match'] = f"{100.0 * np.sum(both_nan) / min_len:.4f}"
            results.append(row)
            continue

        # Same length — element-wise comparison
        both_nan = np.isnan(a) & np.isnan(b)
        either_nan = np.isnan(a) | np.isnan(b)
        only_mat_nan = np.isnan(a) & ~np.isnan(b)
        only_cpp_nan = ~np.isnan(a) & np.isnan(b)
        valid = ~either_nan
        n = len(a)

        if np.sum(valid) > 0:
            diff = np.abs(a[valid] - b[valid])
            max_abs = np.max(diff)
            denom = np.maximum(np.abs(a[valid]), np.abs(b[valid]))
            denom[denom == 0] = 1
            max_rel = np.max(diff / denom)
            mean_abs = np.mean(diff)

            close = np.isclose(a[valid], b[valid], rtol=RTOL, atol=ATOL)
            matching = int(np.sum(both_nan)) + int(np.sum(close))
            pct = 100.0 * matching / n

            row['max_abs_diff'] = f"{max_abs:.10e}"
            row['max_rel_diff'] = f"{max_rel:.10e}"
            row['mean_abs_diff'] = f"{mean_abs:.10e}"
            row['pct_match'] = f"{pct:.4f}"

            if pct == 100.0:
                row['status'] = 'IDENTICAL'
            elif pct >= 99.9:
                row['status'] = 'NEAR_MATCH'
            elif pct >= 95.0:
                row['status'] = 'CLOSE'
            else:
                row['status'] = 'DIFFERENT'
        else:
            # All NaN in at least one
            nan_match = np.sum(both_nan)
            nan_mismatch = np.sum(only_mat_nan) + np.sum(only_cpp_nan)
            row['max_abs_diff'] = 'N/A (all NaN)'
            row['max_rel_diff'] = 'N/A'
            row['mean_abs_diff'] = 'N/A'
            row['pct_match'] = f"{100.0 * nan_match / n:.4f}" if n > 0 else '0'
            row['status'] = 'IDENTICAL' if nan_mismatch == 0 else 'NAN_MISMATCH'

        results.append(row)

    return results


# ──────────────────────────────────────────────────────────────────────────────
# File discovery
# ──────────────────────────────────────────────────────────────────────────────
def find_uuids(directory, ext):
    """Find all *_feature_output.{ext} files and extract UUIDs."""
    pattern = os.path.join(directory, f"*_feature_output.{ext}")
    files = glob.glob(pattern)
    # Also search subdirectories
    pattern2 = os.path.join(directory, "**", f"*_feature_output.{ext}")
    files += glob.glob(pattern2, recursive=True)
    files = list(set(files))

    uuid_map = {}
    for fpath in files:
        basename = os.path.basename(fpath)
        # Extract UUID: everything before _feature_output
        m = re.match(r'(.+?)_feature_output\.' + ext, basename)
        if m:
            uuid_map[m.group(1)] = fpath
    return uuid_map


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────
def main():
    print(f"MATLAB dir: {MATLAB_DIR}")
    print(f"C++ dir:    {CPP_DIR}")
    print()

    # Find files
    mat_files = find_uuids(MATLAB_DIR, "mat")
    bin_files = find_uuids(CPP_DIR, "bin")

    print(f"MATLAB .mat files found: {len(mat_files)}")
    print(f"C++ .bin files found:    {len(bin_files)}")

    common = sorted(set(mat_files.keys()) & set(bin_files.keys()))
    mat_only = sorted(set(mat_files.keys()) - set(bin_files.keys()))
    cpp_only = sorted(set(bin_files.keys()) - set(mat_files.keys()))

    print(f"Common UUIDs: {len(common)}")
    if mat_only:
        print(f"MATLAB-only UUIDs ({len(mat_only)}): {mat_only[:5]}{'...' if len(mat_only)>5 else ''}")
    if cpp_only:
        print(f"C++-only UUIDs ({len(cpp_only)}): {cpp_only[:5]}{'...' if len(cpp_only)>5 else ''}")
    print()

    if not common:
        print("ERROR: No matching UUIDs found. Check paths and file naming.")
        sys.exit(1)

    # Per-file comparison
    all_rows = []
    field_stats = {}  # field_name -> list of pct_match values

    report_cols = ['uuid', 'field', 'status', 'mat_len', 'cpp_len',
                   'max_abs_diff', 'max_rel_diff', 'mean_abs_diff',
                   'pct_match', 'nan_count_mat', 'nan_count_cpp']

    for idx, uuid in enumerate(common):
        print(f"[{idx+1}/{len(common)}] Comparing {uuid}...")

        try:
            mat_data = load_mat_features(mat_files[uuid])
        except Exception as e:
            print(f"  ERROR loading MATLAB file: {e}")
            continue

        try:
            bin_data = load_bin_features(bin_files[uuid])
        except Exception as e:
            print(f"  ERROR loading C++ file: {e}")
            continue

        results = compare_fields(mat_data, bin_data)

        for row in results:
            row['uuid'] = uuid
            all_rows.append(row)

            fname = row['field']
            if fname not in field_stats:
                field_stats[fname] = {'match_pcts': [], 'statuses': []}
            field_stats[fname]['statuses'].append(row['status'])
            try:
                field_stats[fname]['match_pcts'].append(float(row['pct_match']))
            except (ValueError, TypeError):
                pass

    # Write detailed report
    print(f"\nWriting {REPORT_CSV}...")
    with open(REPORT_CSV, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=report_cols, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(all_rows)

    # Write summary
    print(f"Writing {SUMMARY_CSV}...")
    summary_cols = ['field', 'n_subjects', 'n_identical', 'n_near_match', 'n_close',
                    'n_different', 'n_length_mismatch', 'n_mat_only', 'n_cpp_only',
                    'n_nan_mismatch', 'avg_pct_match', 'min_pct_match', 'verdict']

    summary_rows = []
    for fname in sorted(field_stats.keys()):
        stats = field_stats[fname]
        statuses = stats['statuses']
        n = len(statuses)
        n_id = statuses.count('IDENTICAL')
        n_nm = statuses.count('NEAR_MATCH')
        n_cl = statuses.count('CLOSE')
        n_df = statuses.count('DIFFERENT')
        n_lm = statuses.count('LENGTH_MISMATCH')
        n_mo = statuses.count('MAT_ONLY')
        n_co = statuses.count('CPP_ONLY')
        n_nan = statuses.count('NAN_MISMATCH')

        pcts = stats['match_pcts']
        avg_pct = np.mean(pcts) if pcts else 0
        min_pct = np.min(pcts) if pcts else 0

        if n_id == n:
            verdict = 'PERFECT'
        elif n_id + n_nm == n:
            verdict = 'EXCELLENT'
        elif n_id + n_nm + n_cl == n:
            verdict = 'GOOD'
        elif n_df > 0 or n_lm > 0:
            verdict = 'NEEDS_REVIEW'
        else:
            verdict = 'CHECK'

        summary_rows.append({
            'field': fname,
            'n_subjects': n,
            'n_identical': n_id,
            'n_near_match': n_nm,
            'n_close': n_cl,
            'n_different': n_df,
            'n_length_mismatch': n_lm,
            'n_mat_only': n_mo,
            'n_cpp_only': n_co,
            'n_nan_mismatch': n_nan,
            'avg_pct_match': f"{avg_pct:.4f}",
            'min_pct_match': f"{min_pct:.4f}",
            'verdict': verdict,
        })

    with open(SUMMARY_CSV, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=summary_cols)
        writer.writeheader()
        writer.writerows(summary_rows)

    # Print quick summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    perfect = sum(1 for r in summary_rows if r['verdict'] == 'PERFECT')
    excellent = sum(1 for r in summary_rows if r['verdict'] == 'EXCELLENT')
    good = sum(1 for r in summary_rows if r['verdict'] == 'GOOD')
    review = sum(1 for r in summary_rows if r['verdict'] == 'NEEDS_REVIEW')
    check = sum(1 for r in summary_rows if r['verdict'] == 'CHECK')
    total_f = len(summary_rows)

    print(f"  Total fields compared:  {total_f}")
    print(f"  PERFECT  (100%% match):  {perfect}")
    print(f"  EXCELLENT (>=99.9%%):    {excellent}")
    print(f"  GOOD     (>=95%%):       {good}")
    print(f"  NEEDS_REVIEW (<95%%):    {review}")
    print(f"  CHECK    (other):       {check}")
    print()

    if review > 0:
        print("Fields needing review:")
        for r in summary_rows:
            if r['verdict'] == 'NEEDS_REVIEW':
                print(f"  {r['field']:50s}  avg={r['avg_pct_match']}%  min={r['min_pct_match']}%")

    print(f"\nDetailed report: {REPORT_CSV}")
    print(f"Summary report:  {SUMMARY_CSV}")


if __name__ == '__main__':
    main()
