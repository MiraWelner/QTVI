"""
compare_rr_distances.py

Compare R-R interval distances between wave_markings.bin files (folder M)
and per-bin CSV files (folder N).

Folder N structure:
  N/<studyID>_ECG_fs1000_everyRRQTinputIntoEntropy_csv/
    <studyID>_ECG_fs1000_everyRRQTinputIntoEntropy_Rel<bin>_Abs<x>.csv

Folder M structure:
  M/<studyID>_wave_markings.bin

The Rel number in the CSV filename corresponds to the bin index (1-based).
"""

import csv
import math
import os
import re
import statistics
import struct
from dataclasses import dataclass, field
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

    Binary layout per bin (from write_output_binfile):
      9 index arrays (ch1 raw/sq/abs, ch2 raw/sq/abs, ch3 raw/sq/abs)
      2 PPG index arrays (maxAmps, minAmps)
      4 raw signal arrays (ppg, ecg1, ecg2, ecg3)
      6 preprocessed signal arrays
      9 noise flags (1 byte each)
      pairs array
      ppg_bin_indexs
      ecg_bin_indexs

    Returns list of lists: ch1_raw R-peak indices (1-based) per bin.
    """
    bins = []
    with open(path, "rb") as f:
        num_bins = read_uint64(f)
        if num_bins is None:
            return bins

        for _ in range(num_bins):
            # ch1 raw — the one we want
            ch1_raw = read_idx_array(f)
            if ch1_raw is None:
                break

            # skip remaining 8 index arrays
            for _ in range(8):
                if not skip_idx_array(f):
                    return bins

            # 2 PPG index arrays
            if not skip_idx_array(f):
                break
            if not skip_idx_array(f):
                break

            # 4 raw signal arrays
            for _ in range(4):
                if not skip_signal_array(f):
                    return bins

            # 6 preprocessed signal arrays
            for _ in range(6):
                if not skip_signal_array(f):
                    return bins

            # 9 noise flags
            flags = f.read(9)
            if len(flags) < 9:
                break

            # pairs array
            num_pairs_data = f.read(8)
            if len(num_pairs_data) < 8:
                break
            num_pairs = struct.unpack("<Q", num_pairs_data)[0]
            if num_pairs > 0:
                f.seek(num_pairs * 16, 1)

            # ppg_bin_indexs, ecg_bin_indexs
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
    """Extract the Rel number from a CSV filename."""
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
        # Auto-detect delimiter (tab-separated with possible whitespace)
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
# Stats
# ============================================================================


def safe_mean(v):
    return statistics.mean(v) if v else 0.0


def safe_stdev(v):
    return statistics.stdev(v) if len(v) > 1 else 0.0


# ============================================================================
# Main
# ============================================================================


@dataclass
class ComparisonResult:
    study_id: str = ""
    bin_index: int = 0
    bin_rr_count: int = 0
    csv_rr_count: int = 0
    mean_bin_rr: float = 0.0
    mean_csv_rr: float = 0.0
    std_bin_rr: float = 0.0
    std_csv_rr: float = 0.0
    mean_abs_diff: float = float("nan")
    max_abs_diff: float = float("nan")
    counts_match: bool = False
    has_csv: bool = False


def main():
    folder_n = Path(FOLDER_N)
    folder_m = Path(FOLDER_M)

    if not folder_n.exists():
        print(f"Folder N does not exist: {folder_n}")
        return
    if not folder_m.exists():
        print(f"Folder M does not exist: {folder_m}")
        return

    # 1. Index CSV study folders in N
    csv_study_folders = {}
    for entry in folder_n.iterdir():
        if entry.is_dir():
            study_id = extract_study_id_from_folder(entry.name)
            csv_study_folders[study_id] = entry

    # 2. Process each .bin in M
    all_results = []

    for bin_path in sorted(folder_m.glob("*.bin")):
        study_id = extract_study_id_from_bin(bin_path.stem)
        print(f"Processing study: {study_id}")

        bin_data = read_wave_markings(str(bin_path))
        if not bin_data:
            print(f"  No bins read from {bin_path}")
            continue

        # Find matching CSV folder
        if study_id not in csv_study_folders:
            print(f"  No matching CSV folder for study: {study_id}")
            continue

        # Index CSVs by Rel number
        csv_folder = csv_study_folders[study_id]
        csv_by_rel = {}
        for csv_path in csv_folder.glob("*.csv"):
            rel = extract_rel_number(csv_path.name)
            if rel >= 1:
                csv_by_rel[rel] = csv_path

        # 3. Compare per bin
        for bin_idx, ch1_raw in enumerate(bin_data):
            rel_number = bin_idx + 1  # 1-based

            cr = ComparisonResult(study_id=study_id, bin_index=rel_number)

            # Bin R-R intervals: convert from samples to seconds
            bin_rr = []
            if len(ch1_raw) >= 2:
                bin_rr = [
                    (ch1_raw[i] - ch1_raw[i - 1]) / ECG_SAMPLE_RATE
                    for i in range(1, len(ch1_raw))
                ]
            cr.bin_rr_count = len(bin_rr)
            cr.mean_bin_rr = safe_mean(bin_rr)
            cr.std_bin_rr = safe_stdev(bin_rr)

            # Find CSV
            if rel_number not in csv_by_rel:
                cr.has_csv = False
                all_results.append(cr)
                continue

            cr.has_csv = True
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
            elif rpeak_col >= 0 and rpeak_col < len(columns):
                peaks = [v for v in columns[rpeak_col] if not math.isnan(v)]
                csv_rr = [peaks[i] - peaks[i - 1] for i in range(1, len(peaks))]
            elif columns:
                # Fallback: first column
                col0 = [v for v in columns[0] if not math.isnan(v)]
                if col0 and all(v > 0 for v in col0):
                    csv_rr = col0
                else:
                    print(
                        f"  Bin {rel_number}: Can't identify RR column. Headers: {headers}"
                    )

            cr.csv_rr_count = len(csv_rr)
            cr.mean_csv_rr = safe_mean(csv_rr)
            cr.std_csv_rr = safe_stdev(csv_rr)
            cr.counts_match = len(bin_rr) == len(csv_rr)

            # Element-wise comparison
            compare_count = min(len(bin_rr), len(csv_rr))
            if compare_count > 0:
                abs_diffs = [abs(bin_rr[k] - csv_rr[k]) for k in range(compare_count)]
                cr.mean_abs_diff = statistics.mean(abs_diffs)
                cr.max_abs_diff = max(abs_diffs)

            all_results.append(cr)

    # ========================================================================
    # Output to CSV
    # ========================================================================
    out_path = "rr_comparison.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "Study",
                "Bin",
                "BinRR#",
                "CsvRR#",
                "CountsMatch",
                "MeanBinRR_sec",
                "MeanCsvRR_sec",
                "StdBinRR_sec",
                "StdCsvRR_sec",
                "MeanAbsDiff_sec",
                "MaxAbsDiff_sec",
            ]
        )
        for cr in all_results:
            writer.writerow(
                [
                    cr.study_id,
                    cr.bin_index,
                    cr.bin_rr_count,
                    cr.csv_rr_count,
                    "YES"
                    if cr.counts_match
                    else ("NO_CSV" if not cr.has_csv else "NO"),
                    f"{cr.mean_bin_rr:.6f}",
                    f"{cr.mean_csv_rr:.6f}",
                    f"{cr.std_bin_rr:.6f}",
                    f"{cr.std_csv_rr:.6f}",
                    f"{cr.mean_abs_diff:.6f}"
                    if not math.isnan(cr.mean_abs_diff)
                    else "",
                    f"{cr.max_abs_diff:.6f}" if not math.isnan(cr.max_abs_diff) else "",
                ]
            )

    print(f"Results written to {out_path}")


if __name__ == "__main__":
    main()
