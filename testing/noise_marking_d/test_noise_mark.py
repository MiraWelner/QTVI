#!/usr/bin/env python3
"""
print_noise_bin.py -- dump the contents of a *_noise_markings.bin file.

Format (see user_annotation_handler.cpp exportBinary):
  uint64  count
  count x [6 x float64]:  start_sample, end_sample, start_sec, end_sec,
                           label_code, marking_type_code

Usage:
  python print_noise_bin.py path/to/file_noise_markings.bin
"""
import sys, struct

# annotation_types.hpp: label codes used by exportBinary's labelMap
LABEL_CODES = {
    1.0: "PPG", 2.0: "ECG1", 3.0: "ECG2", 4.0: "ECG3",
    5.0: "ABP", 6.0: "ACCEL", 7.0: "ART", 8.0: "ART_PULM",
}
# annotation_types.hpp: noise_types[].code -> label
MARK_TYPE_CODES = {
    1: "1) R Peak Noise", 2: "2) Minor Noise", 3: "3) Blank.+Thresh.",
    4: "4) PVC", 5: "5) PAC", 6: "6) Cond. Delay", 7: "7) AF",
    8: "8) SVT", 9: "9) VT", 10: "Benign Arr.", 11: "Sig. Arr.",
    12: "Other", 13: "Invert/Noninvert",
}

def main(path):
    with open(path, "rb") as f:
        (count,) = struct.unpack("<Q", f.read(8))
        print(f"{path}: {count} segment(s)\n")
        header = f"{'#':>4}  {'label':<9} {'marking_type':<18} {'start_sample':>12} {'end_sample':>12} {'start_sec':>10} {'end_sec':>10}"
        print(header)
        print("-" * len(header))
        for i in range(count):
            row = f.read(48)
            if len(row) < 48:
                print(f"  [truncated: expected {count} rows, file ended at row {i}]")
                break
            start_s, end_s, start_sec, end_sec, label_code, mark_code = struct.unpack("<6d", row)
            label = LABEL_CODES.get(label_code, f"?({label_code:g})")
            mark = MARK_TYPE_CODES.get(int(mark_code), f"?({mark_code:g})")
            print(f"{i:>4}  {label:<9} {mark:<18} {start_s:>12.0f} {end_s:>12.0f} "
                  f"{start_sec:>10.3f} {end_sec:>10.3f}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python print_noise_bin.py <file_noise_markings.bin>")
        sys.exit(1)
    main(sys.argv[1])
