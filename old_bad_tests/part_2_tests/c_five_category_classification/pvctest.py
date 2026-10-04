#!/usr/bin/env python3
"""
pvctest.py -- the first Task C acceptance test, read off the C++ outputs.

  "On a record with known PVCs, the prematurity filter plus voting flags them,
   they are excluded from the reference template but retained with flags, and
   per-bin category percentages are produced. The substituted beat is a smooth
   blend, not a copy."

Nothing here processes a signal. The C++ driver reads the record and writes
CSVs; this file reads those CSVs and the record's ground truth, and checks the
four clauses.

  python3 make_test_edfs.py                              # writes records/
  <your pipeline> records/pvc_isolated.edf                # writes the CSVs
  python3 pvctest.py                                      # this file reads them

EXPECTED OUTPUT FILES, in records/, named {stem}_*.csv:

  pvc_isolated_beats.csv        one row per beat. Columns read here:
                                beat, time_sec, epoch, premature,
                                vote_confirmed, category, handling, retained
  pvc_isolated_bins.csv         one row per 30 s epoch: n_beats,
                                pct_cat1..pct_cat5, pvc_burden_pct
  pvc_isolated_reference.csv    the reference average: column, reference,
                                n_contributing
  pvc_isolated_substituted.csv  every substituted beat, one row per column:
                                beat, column, original, reference, substituted

Categories are 1..5 per Section 4.5; handling is 0 INCLUDE, 1 SUBSTITUTE,
2 EXCLUDE.
"""

import bisect
import csv
import os
import sys

# ---------------------------------------------------------------------------
# EDIT OUTPUT_FOLDER AND TRUTH_DIR. Nothing else should need changing.
#
# OUTPUT_FOLDER is the config.csv "output_folder" column -- for the mesa row,
# D:\USERS\MiraWelner\QTVI\data\output_bigiamy. The subfolder names below are
# exactly the ones deriveSubpaths() in config_loader.cpp builds from it.
#
# The files do NOT all land in one place. post_process.hpp writes per-beat and
# peak data under r_peak_data_path, the template bin under template_path, and
# the quality products under quality_metric, so Task C's outputs follow the
# same split.
# ---------------------------------------------------------------------------
OUTPUT_FOLDER = r"D:\USERS\MiraWelner\QTVI\data\output_bigiamy"

BEATS_DIR       = os.path.join(OUTPUT_FOLDER, "r_peak_finding_output")   # {stem}_beats.csv
TEMPLATES_DIR   = os.path.join(OUTPUT_FOLDER, "template_outputs")       # {stem}_templates.csv
BINS_DIR        = os.path.join(OUTPUT_FOLDER, "quality_metric")         # {stem}_bins.csv
NSVT_DIR        = os.path.join(OUTPUT_FOLDER, "quality_metric")         # {stem}_nsvt.csv
REFERENCE_DIR   = os.path.join(OUTPUT_FOLDER, "quality_metric")         # {stem}_reference.csv
SUBSTITUTED_DIR = os.path.join(OUTPUT_FOLDER, "quality_metric")         # {stem}_substituted.csv

# The records and their ground truth: wherever the EDFs were put.
TRUTH_DIR = r"D:\USERS\MiraWelner\QTVI\data\mesa_raw_files\synth"

DIR_OF = {
    "beats": BEATS_DIR,
    "templates": TEMPLATES_DIR,
    "bins": BINS_DIR,
    "nsvt": NSVT_DIR,
    "reference": REFERENCE_DIR,
    "substituted": SUBSTITUTED_DIR,
    "truth": TRUTH_DIR,
}
STEM = "pvc_isolated"

CAT_BONAFIDE, CAT_ABNORMAL = 1, 2
H_INCLUDE, H_SUBSTITUTE, H_EXCLUDE = 0, 1, 2


def load(kind):
    path = os.path.join(DIR_OF[kind], f"{STEM}_{kind}.csv")
    if not os.path.exists(path):
        d = DIR_OF[kind]
        have = sorted(f for f in os.listdir(d) if f.endswith(".csv")) \
            if os.path.isdir(d) else "(directory does not exist)"
        sys.exit(f"missing: {path}\n"
                 f"expected {STEM}_{kind}.csv in the {kind} folder.\n"
                 f"CSVs there now: {have}")
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def truth():
    path = os.path.join(TRUTH_DIR, f"{STEM}_truth.csv")
    if not os.path.exists(path):
        sys.exit(f"{path} missing; regenerate the record.")
    with open(path, newline="") as f:
        return [(float(r["time_sec"]), r["kind"]) for r in csv.DictReader(f)]


def match_beats(beats, tr, tol=0.15):
    """Pair each C++ beat row with the nearest ground-truth beat."""
    times = [t for t, _ in tr]
    out = []
    for r in beats:
        t = float(r["time_sec"])
        i = bisect.bisect_left(times, t)
        cands = range(max(0, i - 2), min(len(tr), i + 3))
        j = min(cands, key=lambda k: abs(times[k] - t), default=None)
        if j is not None and abs(times[j] - t) < tol:
            out.append((r, tr[j][1]))
    return out


# ---------------------------------------------------------------------------

def clause_prematurity_flags_them(pairs):
    pvc = [r for r, k in pairs if k == "pvc_a"]
    sinus = [r for r, k in pairs if k == "sinus"]
    assert pvc, "no known PVCs matched in the output"

    flagged = sum(int(r["premature"]) for r in pvc)
    assert flagged == len(pvc), \
        f"prematurity flagged {flagged} of {len(pvc)} known PVCs"
    false_pos = sum(int(r["premature"]) for r in sinus)
    assert false_pos <= 0.02 * len(sinus), \
        f"prematurity flagged {false_pos} of {len(sinus)} sinus beats"

    # Voting cannot confirm isolated ectopy: the literal rule needs five flags
    # inside an eight-slot window and an isolated PVC has one flag in its
    # neighbourhood. Asserted here so the limit is recorded; the CONFIRMING
    # half of the clause is exercised on nsvt.edf, where consecutive premature
    # intervals do reach the threshold (see bigeminy.py).
    confirmed = sum(int(r["vote_confirmed"]) for r in pvc)
    assert confirmed == 0, (
        f"5-of-8 voting confirmed {confirmed} of {len(pvc)} isolated PVCs; "
        "the literal rule cannot confirm isolated ectopy")


def clause_excluded_from_reference_but_retained(pairs, reference):
    pvc = [r for r, k in pairs if k == "pvc_a"]
    cats = [int(r["category"]) for r in pvc]
    assert all(c == CAT_ABNORMAL for c in cats), \
        f"known PVCs classified {sorted(set(cats))}, expected all {CAT_ABNORMAL}"
    assert all(int(r["handling"]) == H_EXCLUDE for r in pvc), \
        "excluded: a category-2 beat must not reach the reference"
    assert all(int(r["retained"]) == 1 for r in pvc), \
        "retained with flags: arrhythmia burden needs them"

    n_ref = int(reference[0]["n_contributing"])
    n_include = sum(1 for r, _ in pairs if int(r["handling"]) == H_INCLUDE)
    assert n_ref == n_include, \
        f"reference built from {n_ref} beats, {n_include} were admitted"
    assert n_ref > 0, "the reference template is empty"


def clause_per_bin_category_percentages(bins):
    assert bins, "no per-bin rows written"
    for b in bins:
        pcts = [float(b[f"pct_cat{k}"]) for k in range(1, 6)]
        assert abs(sum(pcts) - 100.0) < 1e-3, \
            f"epoch {b['epoch']}: percentages sum to {sum(pcts):.3f}"
        assert float(b["pvc_burden_pct"]) >= 0.0
    total = sum(int(b["n_beats"]) for b in bins)
    assert total > 0
    with_pvc = [b for b in bins if float(b["pvc_burden_pct"]) > 0.0]
    assert with_pvc, "no epoch reports any PVC burden on a record with PVCs"


def clause_substituted_beat_is_a_blend(subs):
    assert subs, ("no substituted beats written -- the substitution path was "
                  "never exercised on this record, so the clause is untested")
    by_beat = {}
    for r in subs:
        by_beat.setdefault(int(r["beat"]), []).append(r)
    for beat, rows in by_beat.items():
        o = [float(r["original"]) for r in rows]
        ref = [float(r["reference"]) for r in rows]
        out = [float(r["substituted"]) for r in rows]
        assert out != o, f"beat {beat}: substituted row is a copy of the beat"
        assert out != ref, f"beat {beat}: substituted row is a copy of the reference"
        for a, b, c in zip(o, ref, out):
            lo, hi = min(a, b), max(a, b)
            assert lo - 1e-9 <= c <= hi + 1e-9, \
                f"beat {beat}: a sample lies outside the two inputs"
        # alpha = 1/8 on the reference
        for a, b, c in zip(o, ref, out):
            assert abs(c - (0.875 * b + 0.125 * a)) < 1e-6, \
                f"beat {beat}: blend weight is not 1/8"


# ---------------------------------------------------------------------------

def main():
    for kind, d in DIR_OF.items():
        if not os.path.isdir(d):
            sys.exit(f"the {kind} folder does not exist:\n  {d}\n"
                     f"edit the constants at the top of this file.")
    beats = load("beats")
    bins = load("bins")
    reference = load("reference")
    subs = load("substituted")
    pairs = match_beats(beats, truth())

    print(f"\n{STEM}: {len(beats)} beats in the C++ output, "
          f"{len(pairs)} matched to ground truth")
    kinds = {}
    for _, k in pairs:
        kinds[k] = kinds.get(k, 0) + 1
    print("  truth kinds matched:", dict(sorted(kinds.items())))

    checks = [
        ('"the prematurity filter ... flags them" (isolated PVCs)',
         lambda: clause_prematurity_flags_them(pairs)),
        ('"excluded from the reference template but retained with flags"',
         lambda: clause_excluded_from_reference_but_retained(pairs, reference)),
        ('"per-bin category percentages are produced"',
         lambda: clause_per_bin_category_percentages(bins)),
        ('"the substituted beat is a smooth blend, not a copy"',
         lambda: clause_substituted_beat_is_a_blend(subs)),
    ]
    passed = failed = 0
    print()
    for name, fn in checks:
        try:
            fn()
        except AssertionError as exc:
            print(f"  [FAIL] {name}\n         {exc}")
            failed += 1
        except Exception as exc:                            # noqa: BLE001
            print(f"  [ERROR] {name}\n          {type(exc).__name__}: {exc}")
            failed += 1
        else:
            print(f"  [PASS] {name}")
            passed += 1
    print(f"\n{'ALL PASS' if not failed else 'FAILURES'}: "
          f"{passed}/{passed + failed} clauses")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
