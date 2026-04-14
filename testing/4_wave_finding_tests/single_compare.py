"""
compare_single_bin.py

Given a study ID and bin number, output a CSV with side-by-side
R-peak locations and R-R intervals from the wave_markings.bin and
the reference CSV.

Usage:
    python compare_single_bin.py <study_id> <bin_number>
    e.g. python compare_single_bin.py 3010023_20110817 5
"""

import csv
import math
import os
import re
import struct
import sys
from pathlib import Path

# ============================================================================
# CONSTANTS — SET THESE TO YOUR PATHS
# ============================================================================
FOLDER_N = r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_deep"  # CSV folder root
FOLDER_M = r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_mira"  # wave_markings.bin folder

ECG_SAMPLE_RATE = 1000.0


# ============================================================================
# Binary reader — read only the target bin
# ============================================================================


def read_uint64(f):
    data = f.read(8)
    if len(data) < 8:
        return None
    return struct.unpack("<Q", data)[0]


def read_idx_array(f):
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


def skip_one_bin(f):
    """Skip all fields for one bin. Returns False on read failure."""
    for _ in range(9):
        if not skip_idx_array(f):
            return False
    if not skip_idx_array(f):
        return False
    if not skip_idx_array(f):
        return False
    for _ in range(4):
        if not skip_signal_array(f):
            return False
    for _ in range(6):
        if not skip_signal_array(f):
            return False
    if f.read(9) == b"":
        return False
    num_pairs_data = f.read(8)
    if len(num_pairs_data) < 8:
        return False
    num_pairs = struct.unpack("<Q", num_pairs_data)[0]
    if num_pairs > 0:
        f.seek(num_pairs * 16, 1)
    if not skip_pair_vec(f):
        return False
    if not skip_pair_vec(f):
        return False
    return True


def read_bin_rpeaks(bin_path, target_bin):
    """Read ch1 raw R-peaks for a single bin (0-based index)."""
    with open(bin_path, "rb") as f:
        num_bins = read_uint64(f)
        if num_bins is None or target_bin >= num_bins:
            return None

        # Skip preceding bins
        for _ in range(target_bin):
            if not skip_one_bin(f):
                return None

        # Read ch1 raw from target bin
        return read_idx_array(f)


# ============================================================================
# CSV reader
# ============================================================================


def find_csv_file(folder_n, study_id, rel_number):
    """Find the CSV file matching study_id and Rel number."""
    for entry in Path(folder_n).iterdir():
        if not entry.is_dir():
            continue
        name = entry.name
        pos = name.find("_ECG_fs1000")
        if pos == -1:
            continue
        folder_id = re.sub(r"_part\d+$", "", name[:pos])
        if folder_id != study_id:
            continue

        for csv_path in entry.glob("*.csv"):
            m = re.search(r"_Rel(\d+)_", csv_path.name)
            if m and int(m.group(1)) == rel_number:
                return csv_path
    return None


def read_csv_rr(csv_path):
    """Read R-R intervals (sec) and times (sec) from the tab-delimited CSV."""
    times = []
    rr_intervals = []
    with open(csv_path, "r", newline="") as f:
        first_line = f.readline()
        f.seek(0)
        delim = "\t" if "\t" in first_line else ","
        reader = csv.reader(f, delimiter=delim)
        headers = [h.strip().lower() for h in next(reader)]

        # Find columns
        time_col = -1
        rr_col = -1
        for i, h in enumerate(headers):
            if "time" in h and time_col == -1:
                time_col = i
            if "rr" in h and rr_col == -1:
                rr_col = i

        for row in reader:
            if rr_col >= 0 and rr_col < len(row):
                val = row[rr_col].strip()
                if val:
                    try:
                        rr_intervals.append(float(val))
                    except ValueError:
                        rr_intervals.append(float("nan"))
            if time_col >= 0 and time_col < len(row):
                val = row[time_col].strip()
                if val:
                    try:
                        times.append(float(val))
                    except ValueError:
                        times.append(float("nan"))

    return times, rr_intervals


# ============================================================================
# Main
# ============================================================================


def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_single_bin.py <study_id> <bin_number>")
        print("  e.g. python compare_single_bin.py 3010023_20110817 5")
        sys.exit(1)

    study_id = sys.argv[1]
    bin_number = int(sys.argv[2])  # 1-based

    # Read bin file
    bin_path = Path(FOLDER_M) / f"{study_id}_wave_markings.bin"
    if not bin_path.exists():
        print(f"Bin file not found: {bin_path}")
        sys.exit(1)

    ch1_raw = read_bin_rpeaks(str(bin_path), bin_number - 1)
    if ch1_raw is None:
        print(f"Could not read bin {bin_number} from {bin_path}")
        sys.exit(1)

    # Bin R-peak locations in seconds, R-R in seconds
    # Indices are stored 1-based; convert to 0-based then to seconds
    bin_r_sec = [(idx - 1) / ECG_SAMPLE_RATE for idx in ch1_raw]
    bin_rr_sec = [bin_r_sec[i] - bin_r_sec[i - 1] for i in range(1, len(bin_r_sec))]

    # Find and read CSV
    csv_path = find_csv_file(FOLDER_N, study_id, bin_number)
    if csv_path is None:
        print(f"No CSV found for study {study_id}, Rel {bin_number}")
        sys.exit(1)

    csv_times, csv_rr = read_csv_rr(csv_path)

    # Build output: align by beat index
    max_rows = max(len(bin_r_sec), len(csv_times), len(bin_rr_sec) + 1, len(csv_rr) + 1)

    out_path = os.path.join(FOLDER_M, f"{study_id}_bin{bin_number}_comparison.csv")
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "Beat#",
                "Bin_RPeak_sec",
                "CSV_Time_sec",
                "RPeak_Diff_sec",
                "Bin_RR_sec",
                "CSV_RR_sec",
                "RR_Diff_sec",
            ]
        )

        for i in range(max_rows):
            beat = i + 1
            bin_r = f"{bin_r_sec[i]:.6f}" if i < len(bin_r_sec) else ""
            csv_t = f"{csv_times[i]:.6f}" if i < len(csv_times) else ""

            # R-peak location difference
            r_diff = ""
            if i < len(bin_r_sec) and i < len(csv_times):
                d = bin_r_sec[i] - csv_times[i]
                r_diff = f"{d:.6f}"

            # R-R intervals (shifted by 1: first beat has no RR)
            bin_rr_val = ""
            csv_rr_val = ""
            rr_diff = ""
            if i > 0:
                if i - 1 < len(bin_rr_sec):
                    bin_rr_val = f"{bin_rr_sec[i - 1]:.6f}"
                if i - 1 < len(csv_rr):
                    csv_rr_val = f"{csv_rr[i - 1]:.6f}"
                if i - 1 < len(bin_rr_sec) and i - 1 < len(csv_rr):
                    rr_diff = f"{bin_rr_sec[i - 1] - csv_rr[i - 1]:.6f}"

            writer.writerow(
                [beat, bin_r, csv_t, r_diff, bin_rr_val, csv_rr_val, rr_diff]
            )

    print(f"Written to {out_path}")
    print(f"  Bin R-peaks: {len(ch1_raw)}, CSV beats: {len(csv_times)}")


if __name__ == "__main__":
    main()
