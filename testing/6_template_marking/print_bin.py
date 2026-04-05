"""
compare_markings.py  (v3)

Loads a C++ _template_markings.bin and a MATLAB _template_markings.mat,
writes each to its own CSV, then produces a third CSV comparing
onset, peak, dicrotic, and end values — converted to seconds so that
different sampling rates are handled correctly.

Usage:
    python compare_markings.py [options] <markings.bin> <markings.mat>

Options:
    --cpp-rate HZ     PPG sampling rate used by C++ side   (default: 2000)
    --mat-rate HZ     PPG sampling rate used by MATLAB side (default: 256)
    --mat-offset N    Value to subtract from MATLAB bin index to align
                      with C++ 0-based indices (default: 1)
    --debug           Print raw parsing details

Outputs (written to the current working directory):
    cpp_markings.csv       (raw sample indices)
    matlab_markings.csv    (raw sample indices)
    comparison.csv         (matched bins, values in seconds, diffs in ms)
"""

import csv
import math
import os
import struct
import sys
from pathlib import Path

try:
    import scipy.io as sio
except ImportError:
    sys.exit("scipy is required:  pip install scipy")

try:
    import numpy as np
except ImportError:
    sys.exit("numpy is required:  pip install numpy")

DEBUG = False
COMPARE_FIELDS = ["onset", "peak", "dicrotic", "end"]


# ── helpers ────────────────────────────────────────────────────────────────


def _nan_if_neg1(v: int) -> float:
    return float("nan") if v == -1 else float(v)


def _nan_passthrough(v) -> float:
    if v is None:
        return float("nan")
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def _deep_squeeze(v):
    if not isinstance(v, np.ndarray):
        return v
    while isinstance(v, np.ndarray):
        if v.dtype == object:
            if v.size == 1:
                v = v.flat[0]
            else:
                break
        elif v.ndim == 0:
            return v.item()
        elif v.size == 1:
            v = v.flat[0]
        else:
            break
    return v


def _get_field(struct_val, name):
    if isinstance(struct_val, np.void):
        try:
            return _deep_squeeze(struct_val[name])
        except (KeyError, IndexError, ValueError):
            pass
    if isinstance(struct_val, np.ndarray) and struct_val.dtype.names:
        if name in struct_val.dtype.names:
            return _deep_squeeze(struct_val[name])
    if isinstance(struct_val, dict) and name in struct_val:
        return _deep_squeeze(struct_val[name])
    return None


# ── read C++ _template_markings.bin ────────────────────────────────────────

_BIN_RECORD_SIZE = 28


def read_bin(path: str) -> list[dict]:
    file_size = os.path.getsize(path)
    rows = []

    with open(path, "rb") as f:
        (n,) = struct.unpack("<Q", f.read(8))
        expected = 8 + n * _BIN_RECORD_SIZE

        if file_size != expected:
            print(
                f"  *** WARNING: file size {file_size} != expected {expected} "
                f"for {n} bins ***"
            )
            if file_size > expected * 5:
                print(
                    f"      Likely the wrong file (template_info instead of "
                    f"template_markings)."
                )

        for i in range(n):
            raw = f.read(_BIN_RECORD_SIZE)
            if len(raw) < _BIN_RECORD_SIZE:
                print(f"  WARNING: EOF at bin {i}")
                break

            (index,) = struct.unpack_from("<Q", raw, 0)
            br1, br2, br3, ppg_iss = struct.unpack_from("<4B", raw, 8)
            dic, ons, pk, end_ = struct.unpack_from("<4i", raw, 12)

            if DEBUG and i < 3:
                print(
                    f"  [debug] bin[{i}] hex={raw.hex()}  "
                    f"idx={index} dic={dic} ons={ons} pk={pk} end={end_}"
                )

            rows.append(
                {
                    "index": int(index),
                    "bad_r_ch1": bool(br1),
                    "bad_r_ch2": bool(br2),
                    "bad_r_ch3": bool(br3),
                    "ppg_issue": int(ppg_iss),
                    "onset": _nan_if_neg1(ons),
                    "peak": _nan_if_neg1(pk),
                    "dicrotic": _nan_if_neg1(dic),
                    "end": _nan_if_neg1(end_),
                }
            )
    return rows


# ── read MATLAB _template_markings.mat ─────────────────────────────────────


def read_mat(path: str) -> list[dict]:
    for squeeze in (True, False):
        mat = sio.loadmat(path, squeeze_me=squeeze)

        if "template_info" in mat:
            data = mat["template_info"]
        else:
            candidates = [k for k in mat if not k.startswith("__")]
            if len(candidates) == 1:
                data = mat[candidates[0]]
            else:
                sys.exit(f"Cannot find 'template_info'. Keys: {candidates}")

        data = np.asarray(data).flatten()

        if DEBUG:
            print(
                f"  [debug] mat: squeeze={squeeze}, shape={data.shape}, "
                f"dtype={data.dtype}"
            )
            if len(data) > 0:
                e = _deep_squeeze(data[0])
                if hasattr(e, "dtype") and e.dtype.names:
                    print(f"  [debug] mat: fields={e.dtype.names}")

        if len(data) > 0:
            test = _deep_squeeze(data[0])
            if _get_field(test, "Dicrotic") is not None:
                break
    else:
        pass

    rows = []
    for i, cell in enumerate(data):
        s = _deep_squeeze(cell)

        def _field(name, default=float("nan")):
            v = _get_field(s, name)
            if v is None:
                v = _get_field(s, name.lower())
            if v is None:
                if DEBUG and i == 0:
                    print(f"  [debug] mat: '{name}' not found in bin {i}")
                return default
            result = _nan_passthrough(v)
            if DEBUG and i == 0:
                print(f"  [debug] mat: bin[0].{name} = {v} → {result}")
            return result

        idx = _field("index", default=float(i))
        bad_ppg = _field("bad_ppg_templates", default=0.0)
        tmpl_bad = _field("TemplateBad", default=0.0)

        if not math.isnan(bad_ppg) and bad_ppg:
            ppg_issue = 1
        elif not math.isnan(tmpl_bad) and tmpl_bad:
            ppg_issue = 2
        else:
            ppg_issue = 0

        rows.append(
            {
                "index": int(idx) if not math.isnan(idx) else i,
                "bad_r_templates": bool(bad_ppg if not math.isnan(bad_ppg) else 0),
                "bad_ppg_templates": bool(bad_ppg if not math.isnan(bad_ppg) else 0),
                "ppg_issue": ppg_issue,
                "onset": _field("Onset"),
                "peak": _field("Peak"),
                "dicrotic": _field("Dicrotic"),
                "end": _field("End"),
            }
        )

    # Dump fields if everything is NaN
    if rows and all(math.isnan(r["dicrotic"]) for r in rows):
        print(f"\n  *** WARNING: all dicrotic NaN — dumping struct fields:")
        e = _deep_squeeze(data[0])
        if hasattr(e, "dtype") and e.dtype.names:
            for fn in e.dtype.names:
                print(f"    .{fn} = {repr(_deep_squeeze(e[fn]))[:100]}")
        print()

    return rows


# ── write CSV ──────────────────────────────────────────────────────────────


def write_csv(path: str, rows: list[dict]) -> None:
    if not rows:
        print(f"  (no rows for {path})")
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  Wrote {len(rows)} rows → {path}")


# ── comparison (matched by bin index, values in seconds) ───────────────────


def samp_to_sec(val: float, rate: float) -> float:
    """Convert a sample index to seconds.  NaN stays NaN."""
    if math.isnan(val):
        return float("nan")
    return val / rate


def compare(
    cpp_rows: list[dict],
    mat_rows: list[dict],
    cpp_rate: float,
    mat_rate: float,
    mat_offset: int,
) -> list[dict]:
    """
    Match bins by index.  mat_offset is subtracted from the MATLAB index
    before matching (handles MATLAB 1-based → C++ 0-based).
    Values are converted to seconds; diffs reported in milliseconds.
    """
    cpp_by_idx = {r["index"]: r for r in cpp_rows}
    mat_by_idx = {r["index"] - mat_offset: r for r in mat_rows}

    all_indices = sorted(set(cpp_by_idx.keys()) | set(mat_by_idx.keys()))

    matched = 0
    cpp_only = 0
    mat_only = 0
    out = []

    for idx in all_indices:
        cr = cpp_by_idx.get(idx)
        mr = mat_by_idx.get(idx)

        row: dict = {"bin_index": idx}

        if cr and mr:
            matched += 1
            row["status"] = "both"
        elif cr:
            cpp_only += 1
            row["status"] = "cpp_only"
        else:
            mat_only += 1
            row["status"] = "mat_only"

        for field in COMPARE_FIELDS:
            cv_samp = cr[field] if cr else float("nan")
            mv_samp = mr[field] if mr else float("nan")

            cv_sec = samp_to_sec(cv_samp, cpp_rate)
            mv_sec = samp_to_sec(mv_samp, mat_rate)

            row[f"cpp_{field}_samp"] = cv_samp
            row[f"cpp_{field}_sec"] = (
                round(cv_sec, 6) if not math.isnan(cv_sec) else cv_sec
            )
            row[f"mat_{field}_samp"] = mv_samp
            row[f"mat_{field}_sec"] = (
                round(mv_sec, 6) if not math.isnan(mv_sec) else mv_sec
            )

            if math.isnan(cv_sec) or math.isnan(mv_sec):
                row[f"diff_{field}_ms"] = float("nan")
                row[f"match_{field}"] = math.isnan(cv_sec) and math.isnan(mv_sec)
            else:
                diff_ms = (cv_sec - mv_sec) * 1000.0
                row[f"diff_{field}_ms"] = round(diff_ms, 3)
                row[f"match_{field}"] = abs(diff_ms) < 1.0  # within 1 ms

        out.append(row)

    print(
        f"\n  Bin matching:  {matched} matched,  "
        f"{cpp_only} C++ only,  {mat_only} MATLAB only"
    )

    return out


# ── main ───────────────────────────────────────────────────────────────────


def main() -> None:
    global DEBUG

    args = sys.argv[1:]
    cpp_rate = 2000.0
    mat_rate = 256.0
    mat_offset = 1  # subtract from MATLAB index to get 0-based

    # Parse options
    positional = []
    i = 0
    while i < len(args):
        if args[i] == "--debug":
            DEBUG = True
        elif args[i] == "--cpp-rate" and i + 1 < len(args):
            i += 1
            cpp_rate = float(args[i])
        elif args[i] == "--mat-rate" and i + 1 < len(args):
            i += 1
            mat_rate = float(args[i])
        elif args[i] == "--mat-offset" and i + 1 < len(args):
            i += 1
            mat_offset = int(args[i])
        elif args[i].startswith("--"):
            sys.exit(f"Unknown option: {args[i]}")
        else:
            positional.append(args[i])
        i += 1

    if len(positional) != 2:
        print(__doc__)
        sys.exit(1)

    bin_path, mat_path = positional

    if not Path(bin_path).is_file():
        sys.exit(f"File not found: {bin_path}")
    if not Path(mat_path).is_file():
        sys.exit(f"File not found: {mat_path}")

    print(f"C++ PPG rate : {cpp_rate} Hz")
    print(f"MATLAB PPG rate: {mat_rate} Hz")
    print(f"MATLAB index offset: {mat_offset}  (subtracted to align with C++)")
    print(f"Rate ratio (C++/MAT): {cpp_rate / mat_rate:.4f}\n")

    print(f"Reading C++ binary : {bin_path}  ({os.path.getsize(bin_path)} bytes)")
    cpp_rows = read_bin(bin_path)
    print(f"  → {len(cpp_rows)} bins")

    print(f"Reading MATLAB .mat: {mat_path}")
    mat_rows = read_mat(mat_path)
    print(f"  → {len(mat_rows)} bins")

    print("\nWriting CSVs...")
    write_csv("cpp_markings.csv", cpp_rows)
    write_csv("matlab_markings.csv", mat_rows)

    print("\nComparing (values converted to seconds)...")
    comp = compare(cpp_rows, mat_rows, cpp_rate, mat_rate, mat_offset)
    write_csv("comparison.csv", comp)

    # ── summary stats ──
    matched_rows = [r for r in comp if r["status"] == "both"]
    if matched_rows:
        print(f"\n{'─' * 60}")
        print(f"  Summary for {len(matched_rows)} matched bins:")
        print(f"{'─' * 60}")
        print(
            f"  {'Field':>10s}  {'<1ms':>6s}  {'<5ms':>6s}  {'<10ms':>6s}  "
            f"{'Mean±Std (ms)':>18s}  {'Max (ms)':>10s}"
        )

        for field in COMPARE_FIELDS:
            diffs = [
                r[f"diff_{field}_ms"]
                for r in matched_rows
                if not math.isnan(r[f"diff_{field}_ms"])
            ]
            if not diffs:
                print(f"  {field:>10s}  (all NaN)")
                continue

            abs_diffs = [abs(d) for d in diffs]
            n = len(diffs)
            lt1 = sum(1 for d in abs_diffs if d < 1.0)
            lt5 = sum(1 for d in abs_diffs if d < 5.0)
            lt10 = sum(1 for d in abs_diffs if d < 10.0)
            mean_d = sum(diffs) / n
            std_d = (sum((d - mean_d) ** 2 for d in diffs) / n) ** 0.5
            max_d = max(abs_diffs)

            print(
                f"  {field:>10s}  {lt1:>4d}/{n:<1d}  {lt5:>4d}/{n:<1d}  "
                f"{lt10:>4d}/{n:<1d}  {mean_d:>+8.2f} ± {std_d:<6.2f}  "
                f"{max_d:>8.2f}"
            )

        print(f"{'─' * 60}")


if __name__ == "__main__":
    main()
