#!/usr/bin/env python3
"""
bigeminy.py -- the second Task C acceptance test, read off the pipeline outputs.

  "On a record with known bigeminy, the bank converges to exactly two templates
   and the PVC share approaches 50 percent. On a record with a documented NSVT
   run, detectNsvt recovers the run with the correct onset beat and length, and
   a record with isolated unifocal PVCs produces no runs."

Three clauses, three records -- which is how the clause is written ("on a
record ... on a record ... and a record"):

  bigeminy.edf         every other beat a PVC   -> two templates, 50 percent
  nsvt.edf             sinus with 3 documented NSVT runs (6, 4, 9 beats)
  pvc_isolated.edf     a few isolated PVCs      -> no runs

Nothing here processes a signal. Run the pipeline over each record, then run
this file; it reads the CSVs and each record's ground truth.

EXPECTED OUTPUT FILES, in records/, named {stem}_*.csv:

  {stem}_beats.csv      one row per beat. Columns read here: beat, time_sec,
                        epoch, template_id
  {stem}_templates.csv  one row per template per bin: epoch, n_templates,
                        template_id, label, count, share_pct
  {stem}_nsvt.csv       one row per detected run: epoch, start_beat, length,
                        template_id, mean_cycle_ms, rate_bpm, sustained

`label` is the template class: sinus, pvc_a, pvc_b, pac, vt, artifact.

THE VOTE IS TESTED HERE TOO. "The prematurity filter plus voting flags them"
cannot be satisfied by isolated PVCs or by bigeminy -- five flags have to fall
inside an eight-slot window, and neither rhythm ever puts them there. The NSVT
runs do: consecutive premature intervals. Measured on these records:

  pvc_isolated   29 of 29 ectopic beats flagged, 0 confirmed
  bigeminy      359 of 364 flagged,               0 confirmed
  nsvt           14 of 19 flagged,                9 confirmed

So the CONFIRMING half of that clause is checked on nsvt.edf, and the 4-beat
run is deliberately below the threshold: long enough for detectNsvt (minimum
three) and too short for the vote.
"""

import csv
import os
import sys

REC = "records"
BIG = "bigeminy"          # every other beat a PVC
RUN = "nsvt"              # sinus with documented NSVT runs
ISO = "pvc_isolated"      # a few isolated unifocal PVCs
EPOCH = 30.0


def load(stem, name):
    path = os.path.join(REC, f"{stem}_{name}.csv")
    if not os.path.exists(path):
        sys.exit(f"{path} missing. Run the pipeline over "
                 f"{REC}/{stem}.edf first.")
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def truth(stem):
    path = os.path.join(REC, f"{stem}_truth.csv")
    if not os.path.exists(path):
        sys.exit(f"{path} missing; regenerate the record.")
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def epochs_of(tr, regime):
    """Epochs whose beats are entirely inside one regime."""
    per = {}
    for r in tr:
        per.setdefault(int(r["epoch"]), set()).add(r["regime"])
    return sorted(e for e, regs in per.items() if regs == {regime})


# ---------------------------------------------------------------------------

def clause_bank_converges_to_two_templates(templates, big_epochs):
    counts = {}
    for r in templates:
        ep = int(r["epoch"])
        if ep in big_epochs:
            counts[ep] = int(r["n_templates"])
    assert counts, "no template rows for the bigeminy epochs"
    two = [e for e, n in counts.items() if n == 2]
    assert len(two) >= 0.8 * len(counts), (
        f"{len(two)} of {len(counts)} bigeminy bins converged to exactly two "
        f"templates; counts seen: {sorted(set(counts.values()))}")


def clause_pvc_share_approaches_fifty(templates, big_epochs):
    """The ventricular-labelled template should hold about half the bin."""
    shares = []
    for r in templates:
        ep = int(r["epoch"])
        if ep in big_epochs and r["label"].startswith("pvc"):
            shares.append(float(r["share_pct"]))
    assert shares, "no template in any bigeminy bin was labelled ventricular"
    near = [s for s in shares if abs(s - 50.0) < 10.0]
    assert len(near) >= 0.8 * len(shares), (
        f"{len(near)} of {len(shares)} ventricular templates are near 50 "
        f"percent; shares: {[round(s) for s in shares[:8]]}")


def clause_voting_confirms_consecutive_ectopy(beats, tr):
    """The confirming half of test 1's first clause, on the record that can
    satisfy it: runs of consecutive premature intervals."""
    times = [float(r["time_sec"]) for r in tr]
    run_idx = {i for i, r in enumerate(tr) if r["regime"] == "nsvt_run"}
    assert run_idx, "the record documents no runs"

    confirmed = [r for r in beats if int(r["vote_confirmed"]) == 1]
    assert confirmed, ("5-of-8 voting confirmed nothing on a record containing "
                       "runs of consecutive premature beats")
    # every confirmation should sit on or beside a documented run
    for r in confirmed:
        t = float(r["time_sec"])
        near = min(abs(t - times[i]) for i in run_idx)
        assert near < 2.0, \
            f"a confirmation at {t:.1f} s is not near any documented run"


def clause_run_recovered_with_onset_and_length(nsvt, beats, tr):
    """Each documented run must come back with the right onset AND the right
    length. Three runs of different length, because reporting the onset right
    and the length wrong is the easy failure mode and one run cannot catch
    it."""
    # documented runs, from ground truth: contiguous blocks of nsvt_run beats
    docs, cur = [], []
    for r in tr:
        if r["regime"] == "nsvt_run":
            cur.append(float(r["time_sec"]))
        elif cur:
            docs.append(cur); cur = []
    if cur:
        docs.append(cur)
    assert docs, "the record documents no runs"

    assert nsvt, f"no runs detected; {len(docs)} are documented"
    found = []
    for r in nsvt:
        row = next((b for b in beats if int(b["beat"]) == int(r["start_beat"])), None)
        if row is not None:
            found.append((float(row["time_sec"]), int(r["length"]),
                          float(r["rate_bpm"]), int(r["sustained"])))
    found.sort()

    for doc in docs:
        t0, n = doc[0], len(doc)
        hit = [f for f in found if abs(f[0] - t0) < 1.5]
        assert hit, (f"the {n}-beat run at {t0:.1f} s was not detected; "
                     f"detected onsets: {[round(f[0], 1) for f in found]}")
        t, length, rate, sustained = hit[0]
        assert length == n, \
            f"run at {t0:.1f} s: detected length {length}, documented {n}"
        assert rate > 100.0, \
            f"run at {t0:.1f} s: rate {rate:.0f} bpm must exceed 100"
        assert sustained == 0, \
            f"run at {t0:.1f} s is {n} beats, well under 30 s; not sustained"


def clause_isolated_pvcs_produce_no_runs(nsvt, beats, iso_epochs):
    """A few isolated unifocal PVCs: one morphology, never two in a row, so no
    run of three consecutive ectopic beats exists to find."""
    if not nsvt:
        return
    reported = []
    for r in nsvt:
        row = next((b for b in beats if int(b["beat"]) == int(r["start_beat"])), None)
        if row is None or int(row["epoch"]) in iso_epochs:
            reported.append(r)
    assert not reported, (
        f"{len(reported)} run(s) reported on a record of isolated unifocal "
        "PVCs; there is no run of three consecutive ectopic beats in it")


# ---------------------------------------------------------------------------

def main():
    big_t, big_b, big_tpl = truth(BIG), load(BIG, "beats"), load(BIG, "templates")
    run_t, run_b, run_ns = truth(RUN), load(RUN, "beats"), load(RUN, "nsvt")
    iso_t, iso_b, iso_ns = truth(ISO), load(ISO, "beats"), load(ISO, "nsvt")

    big_epochs = set(epochs_of(big_t, "bigeminy"))
    iso_epochs = set(int(r["epoch"]) for r in iso_t)

    print(f"\nrecords read:")
    print(f"  {BIG:18s} {len(big_b):5d} beats, {len(big_epochs)} bigeminy bins")
    print(f"  {RUN:18s} {len(run_b):5d} beats, {len(run_ns)} run(s) reported")
    print(f"  {ISO:18s} {len(iso_b):5d} beats, {len(iso_ns)} run(s) reported")

    checks = [
        ('"the bank converges to exactly two templates"',
         lambda: clause_bank_converges_to_two_templates(big_tpl, big_epochs)),
        ('"the PVC share approaches 50 percent"',
         lambda: clause_pvc_share_approaches_fifty(big_tpl, big_epochs)),
        ('"the prematurity filter plus voting flags them" (consecutive ectopy)',
         lambda: clause_voting_confirms_consecutive_ectopy(run_b, run_t)),
        ('"detectNsvt recovers the run with the correct onset beat and length"',
         lambda: clause_run_recovered_with_onset_and_length(run_ns, run_b, run_t)),
        ('"isolated unifocal PVCs produces no runs"',
         lambda: clause_isolated_pvcs_produce_no_runs(iso_ns, iso_b, iso_epochs)),
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
