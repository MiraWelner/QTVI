#!/usr/bin/env python3
"""
given a *noise_markings.bin file, dump out the contents into the terminal. This is neccecary to visually inspect the contents to verify they match
with what you marked, and with what is listed on the .csv
"""
import os, sys, struct
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
def ask_path(prompt, exts={".bin"}):
    """Prompt until the answer names an existing file with an accepted
    extension. Surrounding quotes are stripped because Windows drag-and-drop
    and 'Copy as path' both add them."""
    while True:
        try:
            answer = input(prompt).strip()
        except (EOFError, KeyboardInterrupt):
            print("\naborted.")
            sys.exit(1)
        if not answer:
            print("        (required)")
            continue
        answer = answer.strip('"').strip("'").strip()
        path = os.path.expanduser(os.path.expandvars(answer))
        ext = os.path.splitext(path)[1].lower()
        if ext not in exts:
            print(f"        need {' or '.join(sorted(exts))}, got "
                  f"'{ext or '(no extension)'}'")
            continue
        if not os.path.isfile(path):
            print(f"        not found: {path}")
            continue
        return path
def main(path):
    with open(path, "rb") as f:
        head = f.read(8)
        if len(head) < 8:
            print(f"{path}: only {len(head)} byte(s) -- too short to hold a "
                  f"segment count; is this a noise_markings .bin?")
            return 1
        (count,) = struct.unpack("<Q", head)
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
    return 0
if __name__ == "__main__":
    sys.exit(main(ask_path("  *_noise_markings.bin : ")))
