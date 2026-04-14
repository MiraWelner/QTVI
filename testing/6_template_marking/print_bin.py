"""
compare_markings.py  (v6)

Compares C++ and MATLAB template markings for a given subject ID.

Usage:
    python compare_markings.py [--debug] <ID> [BIN]

Outputs:
    {ID}_comparison.csv
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

CPP_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\6_template_marking_files\mesa_templatemarking_mira"
)
MAT_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\6_template_marking_files\mesa_templatemarking_daniel"
)


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


def _get_ppg_template(struct_val):
    v = _get_field(struct_val, "ppgTemplate")
    if v is None:
        v = _get_field(struct_val, "ppgtemplate")
    if v is None:
        return None
    if isinstance(v, np.ndarray):
        return v.flatten().astype(float)
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
                return default
            return _nan_passthrough(v)

        idx = _field("index", default=float(i))

        ppg_tmpl = _get_ppg_template(s)
        if ppg_tmpl is not None and len(ppg_tmpl) > 0:
            ppg_amplitude = float(np.nanmax(ppg_tmpl) - np.nanmin(ppg_tmpl))
        else:
            ppg_amplitude = float("nan")

        rows.append(
            {
                "index": int(idx) if not math.isnan(idx) else i,
                "onset": _field("Onset"),
                "peak": _field("Peak"),
                "dicrotic": _field("Dicrotic"),
                "end": _field("End"),
                "ppg_amplitude": ppg_amplitude,
            }
        )

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


# ── comparison ─────────────────────────────────────────────────────────────


SAMPLING_RATE = 1000.0


def _to_sec(v: float) -> float:
    if math.isnan(v):
        return float("nan")
    return v / SAMPLING_RATE


def _diff(a: float, b: float) -> float:
    if math.isnan(a) or math.isnan(b):
        return float("nan")
    return a - b


def compare(cpp_rows: list[dict], mat_rows: list[dict]) -> list[dict]:
    cpp_by_idx = {r["index"]: r for r in cpp_rows}
    mat_by_idx = {r["index"]: r for r in mat_rows}

    all_indices = sorted(set(cpp_by_idx.keys()) | set(mat_by_idx.keys()))
    nan = float("nan")

    out = []
    for idx in all_indices:
        cr = cpp_by_idx.get(idx)
        mr = mat_by_idx.get(idx)

        co = _to_sec(cr["onset"] if cr else nan)
        mo = _to_sec(mr["onset"] if mr else nan)
        cp = _to_sec(cr["peak"] if cr else nan)
        mp = _to_sec(mr["peak"] if mr else nan)
        ce = _to_sec(cr["end"] if cr else nan)
        me = _to_sec(mr["end"] if mr else nan)

        row = {
            "bin_index": idx,
            "cpp_onset_sec": co,
            "mat_onset_sec": mo,
            "diff_onset_sec": _diff(co, mo),
            "cpp_peak_sec": cp,
            "mat_peak_sec": mp,
            "diff_peak_sec": _diff(cp, mp),
            "cpp_end_sec": ce,
            "mat_end_sec": me,
            "diff_end_sec": _diff(ce, me),
            "mat_ppg_amplitude": mr["ppg_amplitude"] if mr else nan,
        }
        out.append(row)

    matched = sum(1 for idx in all_indices if idx in cpp_by_idx and idx in mat_by_idx)
    cpp_only = sum(
        1 for idx in all_indices if idx in cpp_by_idx and idx not in mat_by_idx
    )
    mat_only = sum(
        1 for idx in all_indices if idx not in cpp_by_idx and idx in mat_by_idx
    )
    print(
        f"\n  Bin matching:  {matched} matched,  "
        f"{cpp_only} C++ only,  {mat_only} MATLAB only"
    )

    return out


# ── main ───────────────────────────────────────────────────────────────────


def main() -> None:
    global DEBUG

    args = sys.argv[1:]

    positional = []
    i = 0
    while i < len(args):
        if args[i] == "--debug":
            DEBUG = True
        elif args[i].startswith("--"):
            sys.exit(f"Unknown option: {args[i]}")
        else:
            positional.append(args[i])
        i += 1

    if len(positional) < 1 or len(positional) > 2:
        print(__doc__)
        sys.exit(1)

    subject_id = positional[0]
    bin_filter = int(positional[1]) if len(positional) == 2 else None

    bin_path = CPP_DIR / f"{subject_id}_template_markings.bin"
    mat_path = MAT_DIR / f"{subject_id}_template_markings.mat"

    if not bin_path.is_file():
        sys.exit(f"File not found: {bin_path}")
    if not mat_path.is_file():
        sys.exit(f"File not found: {mat_path}")

    print(f"Subject ID: {subject_id}")
    if bin_filter is not None:
        print(f"Bin filter: {bin_filter}")

    print(f"Reading C++ binary...")
    cpp_rows = read_bin(str(bin_path))
    print(f"  → {len(cpp_rows)} bins")

    print(f"Reading MATLAB .mat...")
    mat_rows = read_mat(str(mat_path))
    print(f"  → {len(mat_rows)} bins")

    # Normalize MATLAB 1-based indices to 0-based
    for r in mat_rows:
        r["index"] = r["index"] - 1

    if bin_filter is not None:
        cpp_rows = [r for r in cpp_rows if r["index"] == bin_filter]
        mat_rows = [r for r in mat_rows if r["index"] == bin_filter]
        print(
            f"\n  Filtered to bin {bin_filter}: {len(cpp_rows)} C++, {len(mat_rows)} MATLAB"
        )

    print("\nComparing...")
    comp = compare(cpp_rows, mat_rows)
    write_csv(f"{subject_id}_comparison.csv", comp)

    # ── summary stats ──
    matched_rows = [
        r
        for r in comp
        if not math.isnan(r["cpp_onset_sec"]) and not math.isnan(r["mat_onset_sec"])
    ]
    if matched_rows:
        fields = ["onset", "peak", "end"]
        print(f"\n{'─' * 60}")
        print(f"  Summary for {len(matched_rows)} matched bins:")
        print(f"{'─' * 60}")
        print(
            f"  {'Field':>10s}  {'Exact':>8s}  {'<5ms':>8s}  "
            f"{'Mean±Std (sec)':>20s}  {'Max (sec)':>12s}"
        )

        for field in fields:
            diffs = [
                r[f"diff_{field}_sec"]
                for r in matched_rows
                if not math.isnan(r[f"diff_{field}_sec"])
            ]
            if not diffs:
                print(f"  {field:>10s}  (all NaN)")
                continue

            abs_diffs = [abs(d) for d in diffs]
            n = len(diffs)
            exact = sum(1 for d in abs_diffs if d == 0)
            lt5ms = sum(1 for d in abs_diffs if d < 0.005)
            mean_d = sum(diffs) / n
            std_d = (sum((d - mean_d) ** 2 for d in diffs) / n) ** 0.5
            max_d = max(abs_diffs)

            print(
                f"  {field:>10s}  {exact:>4d}/{n:<3d}  {lt5ms:>4d}/{n:<3d}  "
                f"{mean_d:>+9.6f} ± {std_d:<9.6f}  {max_d:>10.6f}"
            )

        print(f"{'─' * 60}")


if __name__ == "__main__":
    main()
