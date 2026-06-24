#!/usr/bin/env python3
"""
bin_channel_lengths.py

Scan a folder of Bittium .bin files (produced by file_to_bin.cpp) and write a
CSV summarizing the recorded duration of each channel, per file.

Output CSV layout (files are COLUMNS, fields are ROWS):

    filename               a.bin          b.bin       ...
    ecg1_len               DD:HH:MM:SS    DD:HH:MM:SS ...
    ecg2_len               ...
    ...
    num_file               <count>
    files_missing_accel    <count>
    files_missing_temp     <count>
    files_missing_mark     <count>

Duration per channel = (raw sample count) / (native Hz).  Raw sample counts are
read from the .bin header; the Hz used for the conversion are the constants
below.  A channel that is absent in a file (native rate stored as 0.0 in the
header) is reported as 0 length -> 00:00:00.

Usage:
    python bin_channel_lengths.py <folder_with_bins> [output.csv]
"""

import csv
import glob
import os
import struct
import sys

# ---------------------------------------------------------------------------
# Channel sampling rates (Hz) -- EDIT THESE to match your Bittium config.csv.
# These drive the sample-count -> duration conversion.
# ---------------------------------------------------------------------------
ECG_HZ = 500.0  # ecg_1 / ecg_2 / ecg_3
ACCEL_HZ = 25.0  # accel_x / accel_y / accel_z
TEMP_HZ = 1.0  # temp (DEV_Temperature)
MARKER_HZ = 1.0  # marker

# ---------------------------------------------------------------------------
# .bin header layout (see file_to_bin.hpp).  500-byte header:
#   offset   0 : 4  x uint32   scalars (signal_rate, boolean_rate, ...)
#   offset  16 : 40 x uint32   upsampled sizes
#   offset 176 : 40 x uint32   raw-pair sizes  (count of native samples)
#   offset 336 : 40 x float32  native sampling rates (0.0 = channel absent)
#   offset 496 : 1  x uint32   sleep size
# ---------------------------------------------------------------------------
NUM_CHANNELS = 40
HEADER_SIZE = 500
RAW_SIZES_OFFSET = 176
NATIVE_RATE_OFFSET = 336

# Channel slot indices (must match the ChannelIdx enum in file_to_bin.hpp).
CH_ECG1, CH_ECG2, CH_ECG3 = 1, 2, 3
CH_ACCEL_X, CH_ACCEL_Y, CH_ACCEL_Z = 5, 6, 7
CH_MARKER, CH_TEMP = 8, 9

# (row label, channel index, Hz) in the order requested.
FIELDS = [
    ("ecg1_len", CH_ECG1, ECG_HZ),
    ("ecg2_len", CH_ECG2, ECG_HZ),
    ("ecg3_len", CH_ECG3, ECG_HZ),
    ("accel_x_len", CH_ACCEL_X, ACCEL_HZ),
    ("accel_y_len", CH_ACCEL_Y, ACCEL_HZ),
    ("accel_z_len", CH_ACCEL_Z, ACCEL_HZ),
    ("temp_len", CH_TEMP, TEMP_HZ),
    ("mark_len", CH_MARKER, MARKER_HZ),
]


def fmt_dhms(seconds: float) -> str:
    """Seconds -> DD:HH:MM:SS."""
    total = int(round(seconds))
    d, rem = divmod(total, 86400)
    h, rem = divmod(rem, 3600)
    m, s = divmod(rem, 60)
    return f"{d:02d}:{h:02d}:{m:02d}:{s:02d}"


def is_absent(native_rates, ch_idx):
    """A slot stamped with native rate 0.0 by write_missing() is absent."""
    return native_rates[ch_idx] <= 0.0


def read_header(path):
    """Return (raw_sizes, native_rates) lists of length NUM_CHANNELS."""
    with open(path, "rb") as f:
        header = f.read(HEADER_SIZE)
    if len(header) < HEADER_SIZE:
        raise ValueError(f"file shorter than {HEADER_SIZE}-byte header")

    raw_sizes = struct.unpack_from(f"<{NUM_CHANNELS}I", header, RAW_SIZES_OFFSET)
    native_rates = struct.unpack_from(f"<{NUM_CHANNELS}f", header, NATIVE_RATE_OFFSET)
    return raw_sizes, native_rates


def channel_length(raw_sizes, native_rates, ch_idx, hz):
    """Duration string for one channel; absent channel (native rate 0) -> 0."""
    if is_absent(native_rates, ch_idx):  # write_missing() stamped this slot
        return fmt_dhms(0)
    return fmt_dhms(raw_sizes[ch_idx] / hz)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    folder = sys.argv[1]
    out_csv = sys.argv[2] if len(sys.argv) > 2 else "bin_channel_lengths.csv"

    bin_files = sorted(glob.glob(os.path.join(folder, "*.bin")))
    if not bin_files:
        print(f"No .bin files found in {folder!r}")
        sys.exit(1)

    # Each entry: {"filename": ..., "ecg1_len": ..., ...}
    columns = []
    missing_accel = 0
    missing_temp = 0
    missing_mark = 0
    for path in bin_files:
        name = os.path.basename(path)
        col = {"filename": name}
        try:
            raw_sizes, native_rates = read_header(path)
            for label, ch_idx, hz in FIELDS:
                col[label] = channel_length(raw_sizes, native_rates, ch_idx, hz)
            # A file "misses" accel only if all three axes are absent.
            if all(
                is_absent(native_rates, c) for c in (CH_ACCEL_X, CH_ACCEL_Y, CH_ACCEL_Z)
            ):
                missing_accel += 1
            if is_absent(native_rates, CH_TEMP):
                missing_temp += 1
            if is_absent(native_rates, CH_MARKER):
                missing_mark += 1
        except Exception as e:
            print(f"WARNING: skipping {name}: {e}")
            for label, _, _ in FIELDS:
                col[label] = ""
        columns.append(col)

    # Write transposed: files are columns, fields are rows.
    row_labels = ["filename"] + [label for label, _, _ in FIELDS]
    blanks = [""] * (len(columns) - 1)  # summary scalars sit in the first column
    summary_rows = [
        ["num_file", len(columns)] + blanks,
        ["files_missing_accel", missing_accel] + blanks,
        ["files_missing_temp", missing_temp] + blanks,
        ["files_missing_mark", missing_mark] + blanks,
    ]
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        for label in row_labels:
            writer.writerow([label] + [c[label] for c in columns])
        for row in summary_rows:
            writer.writerow(row)

    print(
        f"Wrote {out_csv} "
        f"({len(columns)} files, {len(row_labels) + len(summary_rows)} rows)"
    )


if __name__ == "__main__":
    main()
