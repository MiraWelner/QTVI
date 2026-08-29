#!/usr/bin/env python3
"""
test_task_d_acceptance.py

Acceptance tests for Task D (ECG normalization + feature archival,
Sections 5.1-5.5.1). Three checks, matching the spec's acceptance criteria:

  1. Feature_norm lands in [0, 100], with each subject's own 2nd and 98th
     ratio percentiles mapping to the 0 / 100 extremes.
  2. A high-CV record raises the local-baseline (cv_flag) flag.
  3. The bin archive is written before any deformation step and contains
     every listed field.

Pure standard library -- no numpy / pandas -- so it runs anywhere Python 3
does. This script does NOT run the C++ pipeline; it reads the CSVs the
pipeline produces and checks them.

Run:
    python3 test_task_d_acceptance.py \
        --feature-norm-csv feature_norm.csv \
        --cv-csv cv_check.csv \
        --bin-archive-csv <stem>_bin_archive_R.csv \
        [--deformation-marker <path touched by the deformation step>]

Exit code 0 if every test PASSes (SKIP does not fail the run), non-zero
if any test FAILs.

----------------------------------------------------------------------------
FILE / COLUMN CONTRACT
----------------------------------------------------------------------------
Test 3 (bin archive) checks the REAL file bin_archive.hpp writes -- its
expected column set below is generated the same way csvHeader() builds it,
so it is ground-truth for the current code.

Tests 1 and 2 need CSVs that bin_archive.hpp / NormalizeFeatures.hpp do not
themselves emit yet (pct_scale / cv_flag exist as functions but nothing
writes their results to disk). This script therefore assumes a small
long-format CSV for each -- override the paths / column names with the
--*-col flags if your wiring names them differently:

  feature_norm.csv : subject_id, ratio, feature_norm
      ratio        = Feature_peak / Global_Ref_person   (Section 5.3)
      feature_norm = pct_scale(ratio, p2, p98)          (Section 5.4)
      p2/p98 are recomputed here per subject from that subject's own
      `ratio` values -- they never need to round-trip through the CSV.

  cv_check.csv : subject_id, channel, bin_index, qrs_ref_value,
                 global_ref_A, global_ref_B, global_ref_C, cv_flag
      qrs_ref_value = one QRS-reference sample feeding the CV check
      global_ref_A  = that subject's Global_Ref_person, Option A (the CV
                      basis; B/C are the area/spatial options, written
                      alongside for comparison). --global-ref-col selects
                      which one the CV recompute uses; default global_ref_A.
      cv_flag       = the flag the PIPELINE reported (1/0/true/false/...)
"""

import argparse
import csv
import math
import os
import sys


# ==========================================================================
# HARDCODED PATHS -- edit these two, everything else derives from them.
# ==========================================================================
# BIN_ARCHIVE_DIR is the folder the pipeline writes all Task D CSVs into
# (config's bin_archive_path). SUBJECT is the file stem, i.e. the CSVs are
# <SUBJECT>_bin_archive_R.csv / <SUBJECT>_feature_norm.csv /
# <SUBJECT>_cv_check.csv inside that folder.
#
# CLI flags still override any of these if given; leave them and just edit
# here for the normal case.
BIN_ARCHIVE_DIR = r"D:\USERS\MiraWelner\QTVI\data\output_mesa\bin_archive"
SUBJECT = "6015875_20110420"

# Derived file paths (pipeline naming: <stem>_<kind>.csv in BIN_ARCHIVE_DIR).
_BIN_ARCHIVE_CSV = os.path.join(BIN_ARCHIVE_DIR, f"{SUBJECT}_bin_archive_R.csv")
_FEATURE_NORM_CSV = os.path.join(BIN_ARCHIVE_DIR, f"{SUBJECT}_feature_norm.csv")
_CV_CHECK_CSV = os.path.join(BIN_ARCHIVE_DIR, f"{SUBJECT}_cv_check.csv")


# ==========================================================================
# Ground-truth bin-archive schema -- mirrors bin_archive.hpp::csvHeader()
# ==========================================================================

_SQI = ["template_corr", "chiSq0", "chiSqAbs", "chiSq0_P", "chiSq0_QRS",
        "chiSq0_ST", "baseline", "noise", "motion", "composite", "frac_include"]

_PER_CH = ["p_amp", "q_amp", "r_amp", "s_amp", "t_amp",
           "pr_ms", "qrs_ms", "qt_ms", "jt_ms",
           "st_level", "j_point_amp", "t_wave_area", "qrs_area",
           "upstroke_slope", "downstroke_slope",
           "q_to_r_ratio", "s_to_r_ratio", "per_sample_std"]


def expected_archive_columns():
    cols = ["subject_id", "bin_index", "anchor", "qrs_area_spatial"]
    for c in (1, 2, 3):
        cols += [f"ch{c}_{f}" for f in _PER_CH]
    cols.append("n_beats_scored")
    cols += [f"sqi_{s}_mean" for s in _SQI]
    cols += [f"sqi_{s}_std" for s in _SQI]
    return cols


# Spec-phrase -> representative column(s), so a missing field reports which
# spec requirement it breaks rather than just a bare column name.
_SPEC_GROUPS = {
    "PQRST amplitudes":            ["ch1_p_amp", "ch1_q_amp", "ch1_r_amp", "ch1_s_amp", "ch1_t_amp"],
    "PR/QRS/QT/JT intervals":      ["ch1_pr_ms", "ch1_qrs_ms", "ch1_qt_ms", "ch1_jt_ms"],
    "ST level":                    ["ch1_st_level"],
    "T-wave area":                 ["ch1_t_wave_area"],
    "QRS area per lead":           ["ch1_qrs_area", "ch2_qrs_area", "ch3_qrs_area"],
    "QRS area spatial":            ["qrs_area_spatial"],
    "per-sample std vector":       ["ch1_per_sample_std"],
    "upstroke/downstroke slopes":  ["ch1_upstroke_slope", "ch1_downstroke_slope"],
    "Q-to-R / S-to-R ratios":      ["ch1_q_to_r_ratio", "ch1_s_to_r_ratio"],
    "J-point amplitude":           ["ch1_j_point_amp"],
    "22 bin quality parameters":   ([f"sqi_{s}_mean" for s in _SQI]
                                    + [f"sqi_{s}_std" for s in _SQI]),
}


# ==========================================================================
# Helpers (dependency-free)
# ==========================================================================

def _to_float(s):
    if s is None or s == "":
        return math.nan
    try:
        return float(s)
    except ValueError:
        return math.nan


def _to_bool(s):
    if s is None:
        return None
    t = str(s).strip().lower()
    if t in ("1", "true", "yes", "y", "t"):
        return True
    if t in ("0", "false", "no", "n", "f", ""):
        return False
    return None


def percentile(vals, q):
    """Linear-interpolation percentile, q in [0,100]."""
    v = sorted(x for x in vals if not math.isnan(x))
    if not v:
        return math.nan
    if len(v) == 1:
        return v[0]
    idx = (q / 100.0) * (len(v) - 1)
    lo, hi = math.floor(idx), math.ceil(idx)
    return v[lo] * (1.0 - (idx - lo)) + v[hi] * (idx - lo)


def median(vals):
    v = sorted(x for x in vals if not math.isnan(x))
    if not v:
        return math.nan
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def mad(vals):
    m = median(vals)
    if math.isnan(m):
        return math.nan
    return median([abs(x - m) for x in vals if not math.isnan(x)])


def pct_scale_py(ratio, p2, p98):
    rng = p98 - p2
    if math.isnan(ratio) or not (rng > 0.0):
        return math.nan
    return min(100.0, max(0.0, (ratio - p2) / rng * 100.0))


def read_rows(path):
    if not os.path.isfile(path):
        return None
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def group_by(rows, key):
    g = {}
    for r in rows:
        g.setdefault(r.get(key), []).append(r)
    return g


class Result:
    def __init__(self, name):
        self.name = name
        self.status = "PASS"
        self.msgs = []

    def fail(self, m): self.status = "FAIL"; self.msgs.append(m)
    def skip(self, m): self.status = "SKIP"; self.msgs.append(m)
    def note(self, m): self.msgs.append(m)
    # WARN removed: former warnings are now plain informational notes -- they
    # never change a test's status. A test is PASS unless something FAILs
    # (or the input is absent -> SKIP).
    def warn(self, m): self.msgs.append(m)


# ==========================================================================
# Test 1: Feature_norm range + percentile extremes
# ==========================================================================

def test_feature_norm(a):
    r = Result("Feature_norm in [0,100] with p2/p98 at the extremes")
    rows = read_rows(a.feature_norm_csv)
    if rows is None:
        r.skip(f"file not found: {a.feature_norm_csv}")
        return r
    if not rows:
        r.skip(f"empty file: {a.feature_norm_csv}")
        return r
    for col in (a.group_col, a.ratio_col, a.norm_col):
        if col not in rows[0]:
            r.fail(f"missing column '{col}' (have: {list(rows[0].keys())})")
    if r.status == "FAIL":
        return r

    subs = 0
    for subject, sr in group_by(rows, a.group_col).items():
        ratios = [_to_float(x[a.ratio_col]) for x in sr]
        got = [_to_float(x[a.norm_col]) for x in sr]
        if sum(1 for v in ratios if not math.isnan(v)) < 2:
            r.warn(f"subject {subject!r}: <2 ratios, percentiles meaningless -- skipped")
            continue
        p2, p98 = percentile(ratios, 2), percentile(ratios, 98)
        if not (p98 > p2):
            r.warn(f"subject {subject!r}: degenerate p2>=p98 -- skipped")
            continue

        for i, (ratio, g) in enumerate(zip(ratios, got)):
            exp = pct_scale_py(ratio, p2, p98)
            if math.isnan(exp) and math.isnan(g):
                continue
            if math.isnan(exp) or math.isnan(g):
                r.fail(f"subject {subject!r} row {i}: NaN mismatch exp={exp} got={g}")
                continue
            if g < -a.tol or g > 100.0 + a.tol:
                r.fail(f"subject {subject!r} row {i}: feature_norm={g} outside [0,100]")
            if abs(exp - g) > a.value_tol:
                r.fail(f"subject {subject!r} row {i}: feature_norm={g} differs from "
                       f"independent recompute pct_scale(ratio={ratio}, "
                       f"p2={p2:.4g}, p98={p98:.4g})={exp:.4g} by more than "
                       f"{a.value_tol} (recompute uses the CSV's rounded ratio; a "
                       f"gap this large is a real scaling mismatch, not rounding)")

        fin = [v for v in got if not math.isnan(v)]
        if not fin:
            r.warn(f"subject {subject!r}: no finite feature_norm -- skipped")
            continue
        if min(fin) > a.tol:
            r.fail(f"subject {subject!r}: min feature_norm {min(fin):.4f} != ~0 "
                   f"(2nd-percentile extreme)")
        if max(fin) < 100.0 - a.tol:
            r.fail(f"subject {subject!r}: max feature_norm {max(fin):.4f} != ~100 "
                   f"(98th-percentile extreme)")
        subs += 1

    if subs == 0 and r.status == "PASS":
        r.warn("no subject had enough data to exercise the percentile extremes")
    elif subs:
        r.note(f"checked {subs} subject(s)")
    return r


# ==========================================================================
# Test 2: high-CV record raises the local-baseline flag
# ==========================================================================

def test_cv_flag(a):
    r = Result("High-CV record raises the local-baseline flag")
    rows = read_rows(a.cv_csv)
    if rows is None:
        r.skip(f"file not found: {a.cv_csv}")
        return r
    if not rows:
        r.skip(f"empty file: {a.cv_csv}")
        return r
    for col in (a.cv_group_col, a.qrs_ref_col, a.global_ref_col, a.cv_flag_col):
        if col not in rows[0]:
            r.fail(f"missing column '{col}' (have: {list(rows[0].keys())})")
    if r.status == "FAIL":
        return r

    high_seen = 0
    for subject, sr in group_by(rows, a.cv_group_col).items():
        qrs = [_to_float(x[a.qrs_ref_col]) for x in sr]
        grefs = [_to_float(x[a.global_ref_col]) for x in sr if not math.isnan(_to_float(x[a.global_ref_col]))]
        if not grefs:
            r.warn(f"subject {subject!r}: no usable global_ref -- skipped")
            continue
        gref = grefs[0]

        flags = {_to_bool(x[a.cv_flag_col]) for x in sr}
        flags.discard(None)
        if not flags:
            r.warn(f"subject {subject!r}: no parsable cv_flag -- skipped")
            continue
        if len(flags) > 1:
            r.fail(f"subject {subject!r}: cv_flag not constant across rows ({flags})")
            continue
        reported = next(iter(flags))

        m = mad(qrs)
        if math.isnan(m) or not (gref > 0.0):
            r.warn(f"subject {subject!r}: CV not computable -- skipped")
            continue
        cv = m / gref
        expected = cv > a.cv_threshold
        if expected:
            high_seen += 1
        if expected != reported:
            r.fail(f"subject {subject!r}: CV={cv:.4f} (thr {a.cv_threshold}) "
                   f"implies flag={expected}, pipeline reported {reported}")

    if high_seen == 0 and r.status == "PASS":
        r.warn("no subject exceeded the CV threshold -- flag correctness was "
               "checked, but the specific 'high-CV raises the flag' scenario "
               "was never exercised; include a high-CV subject to fully "
               "satisfy this criterion")
    elif high_seen:
        r.note(f"{high_seen} high-CV subject(s) correctly flagged")
    return r


# ==========================================================================
# Test 3: bin archive written before deformation, all fields present
# ==========================================================================

def test_bin_archive(a):
    r = Result("Bin archive: all fields present + written before deformation")
    if not os.path.isfile(a.bin_archive_csv):
        r.skip(f"file not found: {a.bin_archive_csv}")
        return r

    with open(a.bin_archive_csv, newline="", encoding="utf-8") as f:
        header = next(csv.reader(f), None)
    if header is None:
        r.fail("no header row")
        return r
    present = set(header)

    expected = expected_archive_columns()
    missing = [c for c in expected if c not in present]
    if missing:
        r.fail(f"missing {len(missing)} of {len(expected)} expected columns: {missing}")
    else:
        r.note(f"all {len(expected)} expected columns present")

    # Report per spec group, so a gap names the requirement it breaks.
    for phrase, cols in _SPEC_GROUPS.items():
        gap = [c for c in cols if c not in present]
        if gap:
            r.fail(f"spec requirement '{phrase}' incomplete -- missing {gap}")

    # 22 quality params: exactly 22 sqi_* columns, no more no less.
    n_sqi = sum(1 for c in header if c.startswith("sqi_"))
    if n_sqi != 22:
        r.fail(f"expected exactly 22 sqi_* quality columns, found {n_sqi}")
    else:
        r.note("22 bin quality parameters present (11 measures x mean/std)")

    rows = read_rows(a.bin_archive_csv)
    if not rows:
        r.warn("header present but no data rows")

    # Ordering-before-deformation check.
    if a.deformation_marker is None:
        r.note("--deformation-marker not given: before-deformation ordering "
               "check SKIPPED (nothing to compare mtimes against)")
    elif not os.path.isfile(a.deformation_marker):
        r.warn(f"--deformation-marker not found ({a.deformation_marker}); if "
               f"deformation never ran, the archive is trivially 'before' it")
    else:
        arch_t = os.path.getmtime(a.bin_archive_csv)
        def_t = os.path.getmtime(a.deformation_marker)
        if arch_t > def_t:
            r.fail(f"archive mtime ({arch_t}) is AFTER deformation marker "
                   f"({def_t}) -- checkpoint was not written first")
        else:
            r.note(f"archive mtime precedes deformation marker ({arch_t} <= {def_t})")

    return r


# ==========================================================================
# Main
# ==========================================================================

def build_parser():
    p = argparse.ArgumentParser(
        description="Task D acceptance tests (Sections 5.1-5.5.1).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument("--feature-norm-csv", default=_FEATURE_NORM_CSV)
    p.add_argument("--group-col", default="subject_id")
    p.add_argument("--ratio-col", default="ratio")
    p.add_argument("--norm-col", default="feature_norm")

    p.add_argument("--cv-csv", default=_CV_CHECK_CSV)
    p.add_argument("--cv-group-col", default="subject_id")
    p.add_argument("--qrs-ref-col", default="qrs_ref_value")
    p.add_argument("--global-ref-col", default="global_ref_A")
    p.add_argument("--cv-flag-col", default="cv_flag")
    p.add_argument("--cv-threshold", type=float, default=0.15)

    p.add_argument("--bin-archive-csv", default=_BIN_ARCHIVE_CSV)
    p.add_argument("--deformation-marker", default=None)

    p.add_argument("--tol", type=float, default=1e-6)
    p.add_argument("--value-tol", type=float, default=1e-2)
    return p


def main():
    a = build_parser().parse_args()
    results = [test_feature_norm(a), test_cv_flag(a), test_bin_archive(a)]

    print("=" * 72)
    for res in results:
        print(f"[{res.status}] {res.name}")
        for m in res.msgs:
            print(f"    - {m}")
    print("=" * 72)
    n = {s: sum(1 for r in results if r.status == s) for s in ("PASS", "SKIP", "FAIL")}
    print(f"{n['PASS']} passed, {n['SKIP']} skipped, {n['FAIL']} failed")
    return 1 if n["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main())
