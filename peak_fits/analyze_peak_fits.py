#!/usr/bin/env python3
"""
Analyse P AND R PEAK POSITIONS across the peak fit models.

Three numbers per template:

    p_peak_abs_ms    P peak position along the template
    r_peak_abs_ms    R peak position along the template
    rp_distance_ms   r_peak_abs_ms - p_peak_abs_ms, so POSITIVE, P before R

The two absolute positions are columns along a template and every template has
R at its own column, so their spread is dominated by that offset rather than by
anything measured -- read them as provenance, and read rp_distance_ms as the
measurement. Clinically the P-to-R distance runs about 100-200 ms; it is the PR
interval measured peak-to-peak rather than onset-to-onset, so it sits below a
conventional PR.

PER-TEMPLATE SCHEMA. One row per (bin_index, channel, template) -- one row per
column you can see and mark in the grid. `template` is PQRST_A / PQRST_B / ...,
the name the grid shows minus its channel prefix, and it is only unique WITHIN
a (bin, channel). The row key is the TRIPLE, not bin_index: two templates in
one bin are different waveforms with their own marks.

Writes three CSVs:
    peak_pos_per_row.csv    one row per (model, bin, channel, template)
    peak_pos_summary.csv    mean / sd per (model, metric), across all rows
    peak_pos_per_row_stats.csv
                            mean / sd per (bin, channel, template, metric),
                            across the MODELS (auto, 5_point, quadratic,
                            cubic) -- how far the choice of peak model moves
                            that one template's landmark

Columns read (R-aligned auto):
    p_peak_x_ms_auto_R, r_peak_x_ms_auto_R

NOTE ON FILENAMES
    The peak radio exports as 5_point / quadratic / cubic. 'cubic' is also an
    on/offset model with the same filename, so keep the peak exports in their
    own directory or the cubic column will be whichever file was written last.

Standard library only. Usage:
    python analyze_peak_fits.py [input_dir] [output_dir]
"""

import csv
import os
import sys
from statistics import mean, stdev

# Model name -> filename prefix. Add to this if more exports appear.
#
# 'auto' is the radio's default -- the BIC contest -- so it is not a fourth
# independent model: per template it reports whichever of the other three won.
# That means it DUPLICATES one of them in the sensitivity sd below, pulling the
# sd down a little. Read the sd as "how much does forcing a model move this
# landmark", and use auto's own row to see which model the detector picked.
MODELS = ["auto", "5_point", "quadratic", "cubic"]

P_ABS = "p_peak_abs_ms"
R_ABS = "r_peak_abs_ms"
RP_DIST = "rp_distance_ms"

# Metric-major, models alphabetical inside each metric -- same layout as the
# on/offset summary so the two read the same way. The distance is last because
# it is the one to look at; the two absolutes above it say where it came from.
METRICS = (P_ABS, R_ABS, RP_DIST)

KEY = ("bin_index", "channel", "template")

# ---- ROWS TO LEAVE OUT ---------------------------------------------------
#
# bad_ecg is the grid's X mark, and it is the reason to skip a row: the export
# writes every template and records the flag as a column rather than dropping
# it, so the filtering is the reader's job. An earlier version of this script
# looked for 'bad_r', which does not exist in this schema -- so the flag was
# silently ignored and a template marked bad went into the statistics anyway.
#
# bad_ppg is NOT consulted. Every metric here is an ECG landmark, so a bad
# pulse says nothing about whether the P and R peaks are trustworthy.
BAD_FLAG = "bad_ecg"

# Ad-hoc exclusions by bin_index, for a bin you know is wrong but have not
# marked. Empty is the normal state -- prefer the X mark, which travels with
# the data and is visible to everything downstream.
EXCLUDE_BINS = {"30"}


def find_file(indir, prefix):
    """Locate the export for one model.

    Require the character after 'prefix_' to be a digit (the subject id), so a
    model name that is a prefix of another cannot match the wrong file.
    """
    hits = []
    for name in sorted(os.listdir(indir)):
        if not name.endswith(".csv"):
            continue
        if not name.startswith(prefix + "_"):
            continue
        rest = name[len(prefix) + 1:]
        if rest and rest[0].isdigit():
            hits.append(name)
    if len(hits) != 1:
        raise SystemExit(
            f"expected exactly one file for '{prefix}' in {indir}, found {hits}")
    return os.path.join(indir, hits[0])


def num(row, col):
    """Float from a cell, or None when absent/blank/unparseable.

    Blank cells are normal: a template with no average for an alignment gets a
    blanked block rather than a fallback to the R-aligned one.
    """
    v = row.get(col)
    if v is None:
        return None
    v = v.strip()
    if not v:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def read_model(path, model):
    """One record per (bin, channel, template) with an R peak.

    R is REQUIRED: it is the anchor, and without it there is nothing to measure
    the P peak against. P is optional -- a template with no detectable P wave
    reports its R position and a blank distance, which is the honest answer and
    keeps the row available for the R statistics.
    """
    rows = list(csv.DictReader(open(path, newline="")))
    if not rows:
        raise SystemExit(f"{os.path.basename(path)}: empty")
    missing = [k for k in KEY if k not in rows[0]]
    if missing:
        raise SystemExit(
            f"{os.path.basename(path)}: no {missing} column -- this looks like "
            "a pre-per-template export. Use the older script for it.")

    out = []
    n_bad = n_excl = 0
    for row in rows:
        if row.get("bin_index", "").strip() in EXCLUDE_BINS:
            n_excl += 1
            continue
        if row.get(BAD_FLAG, "").strip() == "1":
            n_bad += 1
            continue

        members = num(row, "n_members")
        r_abs = num(row, "r_peak_x_ms_auto_R")
        if r_abs is None:
            continue          # no anchor, nothing to measure against
        p_abs = num(row, "p_peak_x_ms_auto_R")

        rec = {"fit_model": model}
        rec.update({k: row.get(k, "").strip() for k in KEY})
        rec["n_members"] = int(members) if members is not None else ""
        rec[P_ABS] = p_abs
        rec[R_ABS] = r_abs
        rec[RP_DIST] = None if p_abs is None else r_abs - p_abs
        out.append(rec)
    return out, len(rows), n_bad, n_excl


def summarize_by_model(rows):
    """mean / sd per (metric, model), across every row. Sample sd (n-1)."""
    groups = {}
    for r in rows:
        groups.setdefault(r["fit_model"], []).append(r)
    out = []
    for metric in METRICS:
        for model in sorted(groups):
            vals = [r[metric] for r in groups[model] if r[metric] is not None]
            out.append({
                "metric": metric,
                "Model": model,
                "n_rows": len(vals),
                "Mean ms": round(mean(vals), 4) if vals else "",
                "Std ms": round(stdev(vals), 4) if len(vals) > 1 else "",
                "Min ms": round(min(vals), 4) if vals else "",
                "Max ms": round(max(vals), 4) if vals else "",
            })
    return out


def summarize_by_row(rows):
    """mean / sd per (bin, channel, template, metric), ACROSS THE MODELS.

    The model-sensitivity view, per template rather than per bin: a large sd_ms
    means the choice of peak model moves THAT template's landmark. Grouping by
    bin alone would mix a bin's several morphologies together.
    """
    groups = {}
    for r in rows:
        groups.setdefault(tuple(r[k] for k in KEY), []).append(r)

    def order(k):
        try:
            return (0, int(k[0]), k[1], k[2])
        except ValueError:
            return (1, 0, k[1], k[2])

    out = []
    for k in sorted(groups, key=order):
        recs = groups[k]
        for metric in METRICS:
            vals = [r[metric] for r in recs if r[metric] is not None]
            rec = dict(zip(KEY, k))
            rec.update({
                "metric": metric,
                "n_models": len(vals),
                "n_members": recs[0]["n_members"],
                "mean_ms": round(mean(vals), 4) if vals else "",
                "sd_ms": round(stdev(vals), 4) if len(vals) > 1 else "",
                "min_ms": round(min(vals), 4) if vals else "",
                "max_ms": round(max(vals), 4) if vals else "",
                "range_ms": round(max(vals) - min(vals), 4) if vals else "",
            })
            out.append(rec)
    return out


def write_csv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path}  ({len(rows)} rows)")


def main():
    indir = sys.argv[1] if len(sys.argv) > 1 else "."
    outdir = sys.argv[2] if len(sys.argv) > 2 else "."
    os.makedirs(outdir, exist_ok=True)

    rows = []
    for model in MODELS:
        path = find_file(indir, model)
        got, total, n_bad, n_excl = read_model(path, model)
        rows.extend(got)
        skipped = ""
        if n_bad or n_excl:
            bits = []
            if n_bad:
                bits.append(f"{n_bad} {BAD_FLAG}")
            if n_excl:
                bits.append(f"{n_excl} in EXCLUDE_BINS")
            skipped = "  (skipped " + ", ".join(bits) + ")"
        print(f"{model:13s} {os.path.basename(path):55s} "
              f"{len(got):4d} of {total} rows usable{skipped}")

    if not rows:
        raise SystemExit("no usable rows: no template had an r_peak")

    write_csv(os.path.join(outdir, "peak_pos_per_row.csv"), rows,
              ["fit_model", *KEY, "n_members", *METRICS])
    write_csv(os.path.join(outdir, "peak_pos_summary.csv"),
              summarize_by_model(rows),
              ["metric", "Model", "n_rows",
               "Mean ms", "Std ms", "Min ms", "Max ms"])
    byrow = summarize_by_row(rows)
    write_csv(os.path.join(outdir, "peak_pos_per_row_stats.csv"), byrow,
              [*KEY, "metric", "n_models", "n_members",
               "mean_ms", "sd_ms", "min_ms", "max_ms", "range_ms"])

    print()
    print("SUMMARY, across all rows")
    hdr = (f"{'metric':18s} {'Model':10s} {'n':>5s} {'Mean ms':>10s} "
           f"{'Std ms':>9s} {'Min ms':>10s} {'Max ms':>10s}")
    print(hdr)
    print("-" * len(hdr))
    for r in summarize_by_model(rows):
        sd = f"{r['Std ms']:9.4f}" if r["Std ms"] != "" else f"{'-':>9s}"
        print(f"{r['metric']:18s} {r['Model']:10s} {r['n_rows']:5d} "
              f"{r['Mean ms']:10.4f} {sd} {r['Min ms']:10.4f} {r['Max ms']:10.4f}")

    print()
    print(f"R-P DISTANCE PER TEMPLATE, across the {len(MODELS)} "
          "models -- sd is model sensitivity")
    hdr = (f"{'bin':>4s} {'ch':4s} {'template':9s} {'n':>3s} "
           f"{'beats':>6s} {'mean':>10s} {'sd':>8s} {'range':>8s}")
    print(hdr)
    print("-" * len(hdr))
    for r in byrow:
        if r["metric"] != RP_DIST:
            continue          # the absolutes are in the CSV
        sd = f"{r['sd_ms']:8.3f}" if r["sd_ms"] != "" else f"{'-':>8s}"
        mean_s = f"{r['mean_ms']:10.3f}" if r["mean_ms"] != "" else f"{'-':>10s}"
        rng = f"{r['range_ms']:8.3f}" if r["range_ms"] != "" else f"{'-':>8s}"
        print(f"{r['bin_index']:>4s} {r['channel']:4s} {r['template']:9s} "
              f"{r['n_models']:3d} {str(r['n_members']):>6s} "
              f"{mean_s} {sd} {rng}")

    # A landmark whose per-row sd is 0 everywhere is not "stable", it is
    # unplumbed -- the radio never reached the code that places it. Worth
    # stating outright rather than leaving a column of zeros to be read as
    # agreement between the models. R is expected to show exactly that: it is
    # the alignment anchor, re-derived from r_col rather than fitted, so the
    # peak radio cannot move it. P moving while R does not is the correct
    # result; NEITHER moving means the radio is not reaching the peak finder.
    print()
    print("model sensitivity, over all templates")
    for metric in METRICS:
        per = [r for r in byrow if r["metric"] == metric]
        sds = [r["sd_ms"] for r in per if r["sd_ms"] not in ("", 0.0)]
        tail = (f"  median sd {sorted(sds)[len(sds) // 2]:.3f} ms"
                if sds else "   <-- identical in every template")
        print(f"  {metric:18s} moved in {len(sds):4d}/{len(per)} templates{tail}")


if __name__ == "__main__":
    main()
