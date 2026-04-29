"""
compare_summary_stats.py

Compute summary statistics for R-R intervals from:
  - Folder N (Deep's output): per-bin CSV files
  - Folder M (Mira's output): wave_markings.bin files

Outputs a side-by-side summary table for all 13 subjects combined.

"N Beats" = number of R-peak detections (not RR intervals).
"Delta Diff" = sum of |mira_peaks - deep_peaks| per bin across all subjects.
"""

import csv
import math
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
PRECISION = 4  # Decimal places for output


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
# Formatting
# ============================================================================

D = PRECISION


def fmt(val):
    if isinstance(val, float) and math.isnan(val):
        return "N/A"
    if isinstance(val, float):
        return f"{val:.{D}f}"
    return str(val)


def iqr_str(q1, q3):
    if isinstance(q1, float) and math.isnan(q1):
        return "N/A"
    return f"({fmt(q1)}, {fmt(q3)})"


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

    # -----------------------------------------------------------------------
    # Per-bin data collection
    # -----------------------------------------------------------------------

    all_mira_rr = []  # all Mira RR intervals (seconds) pooled
    all_deep_rr = []  # all Deep RR intervals (seconds) pooled
    total_mira_beats = 0  # total R-peak count
    total_deep_beats = 0
    delta_diff = 0  # sum of |mira_peaks - deep_peaks| per bin
    n_subjects = 0
    per_subject_mira_beats = []
    per_subject_deep_beats = []

    # Debug: show what IDs were found
    bin_ids = {}
    for bin_path in sorted(folder_m.glob("*.bin")):
        study_id = extract_study_id_from_bin(bin_path.stem)
        bin_ids[study_id] = bin_path

    print(f"Bin files found: {len(bin_ids)}")
    print(f"CSV folders found: {len(csv_study_folders)}")
    matched = sorted(set(bin_ids.keys()) & set(csv_study_folders.keys()))
    unmatched_bin = sorted(set(bin_ids.keys()) - set(csv_study_folders.keys()))
    unmatched_csv = sorted(set(csv_study_folders.keys()) - set(bin_ids.keys()))
    print(f"Matched subjects: {len(matched)}")
    if unmatched_bin:
        print(f"  Bin-only (no CSV match): {unmatched_bin}")
    if unmatched_csv:
        print(f"  CSV-only (no bin match): {unmatched_csv}")
    print()

    for study_id in matched:
        bin_path = bin_ids[study_id]
        print(f"Processing study: {study_id}")

        # ---- Mira's data (bin files) ----
        bin_data = read_wave_markings(str(bin_path))
        if not bin_data:
            print(f"  No bins read from {bin_path}")
            continue

        # ---- Deep's data (CSV files) ----
        csv_folder = csv_study_folders[study_id]
        csv_by_rel = {}
        for csv_path in csv_folder.glob("*.csv"):
            rel = extract_rel_number(csv_path.name)
            if rel >= 1:
                csv_by_rel[rel] = csv_path

        n_subjects += 1
        subject_mira_beats = 0
        subject_deep_beats = 0

        for bin_idx, ch1_raw in enumerate(bin_data):
            rel_number = bin_idx + 1  # 1-based

            # Mira: R-peak count and RR intervals for this bin
            mira_n_peaks = len(ch1_raw)
            subject_mira_beats += mira_n_peaks
            mira_rr_bin = []
            if mira_n_peaks >= 2:
                mira_rr_bin = [
                    (ch1_raw[i] - ch1_raw[i - 1]) / ECG_SAMPLE_RATE
                    for i in range(1, mira_n_peaks)
                ]
            all_mira_rr.extend(mira_rr_bin)

            # Deep: read CSV for this bin (0 if missing)
            deep_n_peaks = 0
            deep_rr_bin = []

            if rel_number in csv_by_rel:
                headers, columns = read_csv_file(str(csv_by_rel[rel_number]))

                rr_col = -1
                for c, h in enumerate(headers):
                    hl = h.lower()
                    if "rr" in hl and rr_col == -1:
                        rr_col = c

                if rr_col >= 0 and rr_col < len(columns):
                    deep_rr_bin = [v for v in columns[rr_col] if not math.isnan(v)]
                    # Each row = one beat with its RR interval to the previous beat.
                    # Number of R-peaks = number of RR intervals + 1
                    # (the first peak has no preceding interval).
                    deep_n_peaks = len(deep_rr_bin) + 1

            subject_deep_beats += deep_n_peaks
            all_deep_rr.extend(deep_rr_bin)

            # Delta diff: per-bin absolute difference in peak counts
            delta_diff += abs(mira_n_peaks - deep_n_peaks)

        total_mira_beats += subject_mira_beats
        total_deep_beats += subject_deep_beats
        per_subject_mira_beats.append(subject_mira_beats)
        per_subject_deep_beats.append(subject_deep_beats)

        print(f"  Mira beats: {subject_mira_beats}  Deep beats: {subject_deep_beats}")

    # ========================================================================
    # Compute summary statistics
    # ========================================================================

    def compute_pooled_stats(all_rr_sec):
        """Compute rate and interval stats from pooled RR intervals (seconds)."""
        rr_ms = [rr * 1000.0 for rr in all_rr_sec]
        rates_bpm = [60.0 / rr for rr in all_rr_sec if rr > 0]

        mean_interval = statistics.mean(rr_ms) if rr_ms else float("nan")
        q1_int, med_int, q3_int = quartiles(rr_ms)

        mean_rate = statistics.mean(rates_bpm) if rates_bpm else float("nan")
        q1_rate, med_rate, q3_rate = quartiles(rates_bpm)

        return {
            "mean_rate_bpm": mean_rate,
            "median_rate_bpm": med_rate,
            "rate_iqr_q1": q1_rate,
            "rate_iqr_q3": q3_rate,
            "mean_interval_ms": mean_interval,
            "median_interval_ms": med_int,
            "interval_iqr_q1": q1_int,
            "interval_iqr_q3": q3_int,
        }

    mira = compute_pooled_stats(all_mira_rr)
    deep = compute_pooled_stats(all_deep_rr)

    mean_mira_bpp = total_mira_beats / n_subjects if n_subjects > 0 else 0
    mean_deep_bpp = total_deep_beats / n_subjects if n_subjects > 0 else 0

    # ========================================================================
    # Print results
    # ========================================================================

    print("\n" + "=" * 75)
    print(f"SUMMARY STATISTICS — ALL {n_subjects} SUBJECTS")
    print("=" * 75)

    header = f"{'Metric':<35} {'Mira':>18} {'Deep':>18}"
    print(header)
    print("-" * 75)

    rows = [
        ("N Beats Total", fmt(float(total_mira_beats)), fmt(float(total_deep_beats))),
        ("Delta Diff (sum |peak diff|/bin)", fmt(float(delta_diff)), ""),
        ("Mean N Beats Per Patient", fmt(mean_mira_bpp), fmt(mean_deep_bpp)),
        (
            "Mean Rate (beat/min)",
            fmt(mira["mean_rate_bpm"]),
            fmt(deep["mean_rate_bpm"]),
        ),
        (
            "Median Rate (beat/min)",
            fmt(mira["median_rate_bpm"]),
            fmt(deep["median_rate_bpm"]),
        ),
        (
            "Rate IQR (beat/min)",
            iqr_str(mira["rate_iqr_q1"], mira["rate_iqr_q3"]),
            iqr_str(deep["rate_iqr_q1"], deep["rate_iqr_q3"]),
        ),
        (
            "Mean Interval (ms)",
            fmt(mira["mean_interval_ms"]),
            fmt(deep["mean_interval_ms"]),
        ),
        (
            "Median Interval (ms)",
            fmt(mira["median_interval_ms"]),
            fmt(deep["median_interval_ms"]),
        ),
        (
            "Interval IQR (ms)",
            iqr_str(mira["interval_iqr_q1"], mira["interval_iqr_q3"]),
            iqr_str(deep["interval_iqr_q1"], deep["interval_iqr_q3"]),
        ),
    ]

    for metric, mira_val, deep_val in rows:
        print(f"{metric:<35} {mira_val:>18} {deep_val:>18}")

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
