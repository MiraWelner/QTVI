"""
compare_summary_stats.py

Compute summary statistics for R-R intervals from:
  - Folder N (Deep's output): per-bin CSV files
  - Folder M (Mira's output): wave_markings.bin files

Outputs a side-by-side summary table for all 13 subjects combined.
"""

import csv
import math
import os
import re
import statistics
import struct
from pathlib import Path

# ============================================================================
# CONSTANTS — SET THESE TO YOUR PATHS
# ============================================================================
FOLDER_N = r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_deep"  # CSV folder root
FOLDER_M = r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_mira"  # wave_markings.bin folder

ECG_SAMPLE_RATE = 1000.0  # Both systems now at 1 kHz


# ============================================================================
# Binary reader for wave_markings.bin
# ============================================================================


def read_uint64(f):
    data = f.read(8)
    if len(data) < 8:
        return None
    return struct.unpack("<Q", data)[0]


def read_idx_array(f):
    """Read a uint64 count followed by that many uint64 values."""
    sz = read_uint64(f)
    if sz is None:
        return None
    if sz == 0:
        return []
    data = f.read(sz * 8)
    return list(struct.unpack(f"<{sz}Q", data))


def skip_idx_array(f):
    sz = read_uint64(f)
    if sz is None:
        return False
    if sz > 0:
        f.seek(sz * 8, 1)
    return True


def skip_signal_array(f):
    sz = read_uint64(f)
    if sz is None:
        return False
    if sz > 0:
        f.seek(sz * 8, 1)
    return True


def skip_pair_vec(f):
    sz = read_uint64(f)
    if sz is None:
        return False
    if sz > 0:
        f.seek(sz * 16, 1)
    return True


def read_wave_markings(path):
    """
    Read wave_markings.bin, extract ch1 raw R-peak indices per bin.
    Returns list of lists: ch1_raw R-peak indices per bin.
    """
    bins = []
    with open(path, "rb") as f:
        num_bins = read_uint64(f)
        if num_bins is None:
            return bins

        for _ in range(num_bins):
            ch1_raw = read_idx_array(f)
            if ch1_raw is None:
                break

            for _ in range(8):
                if not skip_idx_array(f):
                    return bins

            if not skip_idx_array(f):
                break
            if not skip_idx_array(f):
                break

            for _ in range(4):
                if not skip_signal_array(f):
                    return bins

            for _ in range(6):
                if not skip_signal_array(f):
                    return bins

            flags = f.read(9)
            if len(flags) < 9:
                break

            num_pairs_data = f.read(8)
            if len(num_pairs_data) < 8:
                break
            num_pairs = struct.unpack("<Q", num_pairs_data)[0]
            if num_pairs > 0:
                f.seek(num_pairs * 16, 1)

            if not skip_pair_vec(f):
                break
            if not skip_pair_vec(f):
                break

            bins.append(ch1_raw)

    return bins


# ============================================================================
# CSV helpers
# ============================================================================


def extract_rel_number(filename):
    m = re.search(r"_Rel(\d+)_", filename)
    return int(m.group(1)) if m else -1


def extract_study_id_from_folder(folder_name):
    pos = folder_name.find("_ECG_fs1000")
    if pos != -1:
        folder_name = folder_name[:pos]
    return re.sub(r"_part\d+$", "", folder_name)


def extract_study_id_from_bin(filename):
    pos = filename.find("_wave_markings")
    return filename[:pos] if pos != -1 else filename


def read_csv_file(path):
    """Read tab-delimited CSV, return (headers, columns)."""
    headers = []
    columns = []
    with open(path, "r", newline="") as f:
        first_line = f.readline()
        f.seek(0)

        if "\t" in first_line:
            reader = csv.reader(f, delimiter="\t")
        else:
            reader = csv.reader(f)

        try:
            headers = [h.strip() for h in next(reader)]
        except StopIteration:
            return headers, columns
        columns = [[] for _ in headers]
        for row in reader:
            for i, val in enumerate(row):
                if i >= len(columns):
                    break
                val = val.strip()
                if not val:
                    continue
                try:
                    columns[i].append(float(val))
                except ValueError:
                    columns[i].append(float("nan"))
    return headers, columns


# ============================================================================
# Quartile helpers
# ============================================================================


def quartiles(vals):
    """Return (Q1, median, Q3) using statistics.quantiles (Python 3.8+)."""
    if len(vals) < 2:
        return (float("nan"), float("nan"), float("nan"))
    q1, med, q3 = statistics.quantiles(vals, n=4)
    return (q1, med, q3)


# ============================================================================
# Main
# ============================================================================


def main():
    folder_n = Path(FOLDER_N)
    folder_m = Path(FOLDER_M)

    if not folder_n.exists():
        print(f"Folder N does not exist: {folder_n}")
        return
    if not folder_m.exists():
        print(f"Folder M does not exist: {folder_m}")
        return

    # Index CSV study folders in N (Deep's output)
    csv_study_folders = {}
    for entry in folder_n.iterdir():
        if entry.is_dir():
            study_id = extract_study_id_from_folder(entry.name)
            csv_study_folders[study_id] = entry

    # Collect all RR intervals per subject and globally
    # Structure: {study_id: [list of RR intervals in seconds]}
    mira_rr_by_subject = {}
    deep_rr_by_subject = {}

    # Also track beat counts (R-peaks, not intervals) per subject
    mira_beats_by_subject = {}
    deep_beats_by_subject = {}

    for bin_path in sorted(folder_m.glob("*.bin")):
        study_id = extract_study_id_from_bin(bin_path.stem)
        print(f"Processing study: {study_id}")

        # ---- Mira's data (bin files) ----
        bin_data = read_wave_markings(str(bin_path))
        if not bin_data:
            print(f"  No bins read from {bin_path}")
            continue

        mira_rr_subject = []
        mira_beats_subject = 0
        for ch1_raw in bin_data:
            mira_beats_subject += len(ch1_raw)
            if len(ch1_raw) >= 2:
                for i in range(1, len(ch1_raw)):
                    rr_sec = (ch1_raw[i] - ch1_raw[i - 1]) / ECG_SAMPLE_RATE
                    mira_rr_subject.append(rr_sec)

        mira_rr_by_subject[study_id] = mira_rr_subject
        mira_beats_by_subject[study_id] = mira_beats_subject

        # ---- Deep's data (CSV files) ----
        if study_id not in csv_study_folders:
            print(f"  No matching CSV folder for study: {study_id}")
            continue

        csv_folder = csv_study_folders[study_id]
        csv_by_rel = {}
        for csv_path in csv_folder.glob("*.csv"):
            rel = extract_rel_number(csv_path.name)
            if rel >= 1:
                csv_by_rel[rel] = csv_path

        deep_rr_subject = []
        deep_beats_subject = 0

        for bin_idx in range(len(bin_data)):
            rel_number = bin_idx + 1
            if rel_number not in csv_by_rel:
                continue

            headers, columns = read_csv_file(str(csv_by_rel[rel_number]))

            # Auto-detect RR or R-peak column
            csv_rr = []
            rr_col = -1
            rpeak_col = -1

            for c, h in enumerate(headers):
                hl = h.lower()
                if "rr" in hl and rr_col == -1:
                    rr_col = c
                if (
                    any(k in hl for k in ("r_peak", "rpeak", "r_index"))
                    and rpeak_col == -1
                ):
                    rpeak_col = c

            if rr_col >= 0 and rr_col < len(columns):
                csv_rr = [v for v in columns[rr_col] if not math.isnan(v)]
                # Each row is one beat with its RR interval
                deep_beats_subject += len(csv_rr)
            elif rpeak_col >= 0 and rpeak_col < len(columns):
                peaks = [v for v in columns[rpeak_col] if not math.isnan(v)]
                deep_beats_subject += len(peaks)
                csv_rr = [peaks[i] - peaks[i - 1] for i in range(1, len(peaks))]
            elif columns:
                col0 = [v for v in columns[0] if not math.isnan(v)]
                if col0 and all(v > 0 for v in col0):
                    csv_rr = col0
                    deep_beats_subject += len(csv_rr)

            deep_rr_subject.extend(csv_rr)

        deep_rr_by_subject[study_id] = deep_rr_subject
        deep_beats_by_subject[study_id] = deep_beats_subject

    # ========================================================================
    # Compute summary statistics
    # ========================================================================

    def compute_stats(rr_by_subject, beats_by_subject, label):
        """Compute all requested stats. RR intervals are in seconds."""
        all_rr = []
        for sid in rr_by_subject:
            all_rr.extend(rr_by_subject[sid])

        total_beats = sum(beats_by_subject.get(sid, 0) for sid in beats_by_subject)
        n_subjects = len(rr_by_subject)
        mean_beats_per_patient = total_beats / n_subjects if n_subjects > 0 else 0

        # Convert RR (seconds) to ms and bpm
        rr_ms = [rr * 1000.0 for rr in all_rr]
        rates_bpm = [60.0 / rr for rr in all_rr if rr > 0]

        # Interval stats (ms)
        mean_interval = statistics.mean(rr_ms) if rr_ms else float("nan")
        q1_int, med_int, q3_int = quartiles(rr_ms)

        # Rate stats (bpm)
        mean_rate = statistics.mean(rates_bpm) if rates_bpm else float("nan")
        q1_rate, med_rate, q3_rate = quartiles(rates_bpm)

        # Delta diff: sum of absolute differences isn't applicable for a single
        # source — we'll compute this separately between the two sources.

        return {
            "label": label,
            "n_beats_total": total_beats,
            "n_subjects": n_subjects,
            "mean_beats_per_patient": mean_beats_per_patient,
            "mean_rate_bpm": mean_rate,
            "median_rate_bpm": med_rate,
            "rate_iqr_q1": q1_rate,
            "rate_iqr_q3": q3_rate,
            "mean_interval_ms": mean_interval,
            "median_interval_ms": med_int,
            "interval_iqr_q1": q1_int,
            "interval_iqr_q3": q3_int,
        }

    mira_stats = compute_stats(mira_rr_by_subject, mira_beats_by_subject, "Mira")
    deep_stats = compute_stats(deep_rr_by_subject, deep_beats_by_subject, "Deep")

    # Delta diff: sum of |mira_beat_count - deep_beat_count| across all subjects
    common_subjects = sorted(
        set(mira_beats_by_subject.keys()) & set(deep_beats_by_subject.keys())
    )
    delta_diff = sum(
        abs(mira_beats_by_subject[sid] - deep_beats_by_subject[sid])
        for sid in common_subjects
    )

    # ========================================================================
    # Print results
    # ========================================================================

    def fmt(val, decimals=2):
        if isinstance(val, float) and math.isnan(val):
            return "N/A"
        if isinstance(val, float):
            return f"{val:.{decimals}f}"
        return str(val)

    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS — ALL 13 SUBJECTS")
    print("=" * 70)

    header = f"{'Metric':<35} {'Mira':>15} {'Deep':>15}"
    print(header)
    print("-" * 70)

    rows = [
        (
            "N Beats Total",
            fmt(mira_stats["n_beats_total"]),
            fmt(deep_stats["n_beats_total"]),
        ),
        ("Delta Diff (|beat count diff| sum)", fmt(delta_diff), ""),
        (
            "Mean N Beats Per Patient",
            fmt(mira_stats["mean_beats_per_patient"]),
            fmt(deep_stats["mean_beats_per_patient"]),
        ),
        (
            "Mean Rate (beat/min)",
            fmt(mira_stats["mean_rate_bpm"]),
            fmt(deep_stats["mean_rate_bpm"]),
        ),
        (
            "Median Rate (beat/min)",
            fmt(mira_stats["median_rate_bpm"]),
            fmt(deep_stats["median_rate_bpm"]),
        ),
        (
            "Rate IQR (beat/min)",
            f"({fmt(mira_stats['rate_iqr_q1'])}, {fmt(mira_stats['rate_iqr_q3'])})",
            f"({fmt(deep_stats['rate_iqr_q1'])}, {fmt(deep_stats['rate_iqr_q3'])})",
        ),
        (
            "Mean Interval (ms)",
            fmt(mira_stats["mean_interval_ms"]),
            fmt(deep_stats["mean_interval_ms"]),
        ),
        (
            "Median Interval (ms)",
            fmt(mira_stats["median_interval_ms"]),
            fmt(deep_stats["median_interval_ms"]),
        ),
        (
            "Interval IQR (ms)",
            f"({fmt(mira_stats['interval_iqr_q1'])}, {fmt(mira_stats['interval_iqr_q3'])})",
            f"({fmt(deep_stats['interval_iqr_q1'])}, {fmt(deep_stats['interval_iqr_q3'])})",
        ),
    ]

    for metric, mira_val, deep_val in rows:
        print(f"{metric:<35} {mira_val:>15} {deep_val:>15}")

    # ========================================================================
    # Also write to CSV
    # ========================================================================
    out_path = "rr_summary_stats.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Metric", "Mira", "Deep"])
        for metric, mira_val, deep_val in rows:
            writer.writerow([metric, mira_val, deep_val])

    print(f"\nResults also written to {out_path}")


if __name__ == "__main__":
    main()
