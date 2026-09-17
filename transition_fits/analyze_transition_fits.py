#!/usr/bin/env python3
"""
Compare auto-detected ECG landmark timings across on/offset fit models.

Reads the *_template_markings.csv exports (one per fit model), pulls the
R-aligned AUTO landmark positions, and derives the two intervals:

    QRS = s_end - q_onset
    QT  = t_end - q_onset

UPDATED FOR THE PER-TEMPLATE SCHEMA. The export used to be one row per bin
with the channel as a column suffix (q_onset_ch1_x_ms_auto_R). It is now one
row per (bin_index, channel, template) -- one row per column you can see and
mark in the grid -- and the suffix is gone (q_onset_x_ms_auto_R). So:

  * the row key is the TRIPLE, not bin_index. Two templates in one bin are
    different waveforms with their own marks; averaging them would be the same
    mistake the old export made by reporting slot 0 only.
  * channel is read from its column instead of being sniffed from which
    _chN columns happen to be populated.
  * `template` is PQRST_A / PQRST_B / ... -- the name the grid shows, minus its
    channel prefix. It is only unique WITHIN a (bin, channel).

Writes three CSVs:
    fit_model_per_row.csv    one row per (model, bin, channel, template)
    fit_model_summary.csv    mean / sd per (model, metric), across all rows
    fit_model_per_row_stats.csv
                             mean / sd per (bin, channel, template, metric),
                             across the MODELS -- how far the choice of
                             on/offset model moves that one template's landmark.
                             Auto is one of the models, so its distance from
                             each forced model is readable directly.

Columns read (R-aligned auto):
    q_onset_x_ms_auto_R, s_end_x_ms_auto_R, t_end_x_ms_auto_R

Standard library only. Usage:
    python analyze_transition_fits.py [input_dir] [output_dir]
"""

import csv
import os
import sys
from statistics import mean, stdev

# Model name -> filename prefix. Add to this if more exports appear.
# "auto" is the BIC selector's own choice, exported the same way as a forced
# model. It belongs in the comparison rather than outside it: the question the
# five forced exports exist to answer is how far each model moves a landmark
# AWAY FROM what Auto picked, and that needs Auto in the same table. Models are
# sorted alphabetically inside each metric, so it lands first.
MODELS = ["auto", "linear", "sigmoid", "cubic", "cubic_spline", "frac_poly"]

METRICS = ("q_onset_ms", "s_end_ms", "qrs_ms", "t_end_ms", "qt_ms")

# Row order in the summary CSV: metric-major, positions before intervals, with
# the models alphabetical inside each metric. Grouping by metric is what makes
# the models directly comparable down a column.
SUMMARY_METRIC_ORDER = ("q_onset_ms", "s_end_ms", "t_end_ms", "qrs_ms", "qt_ms")

KEY = ("bin_index", "channel", "template")


def find_file(indir, prefix):
    """Locate the export for one model.

    'cubic' and 'cubic_spline' share a prefix, so an exact-prefix match is not
    enough -- 'cubic_*' would also match the spline file. Require the character
    after the prefix to be '_' followed by a digit (the subject id), which
    'cubic_spline_...' fails.
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
    """One record per (bin, channel, template) that has all three landmarks."""
    rows = list(csv.DictReader(open(path, newline="")))
    if not rows:
        raise SystemExit(f"{os.path.basename(path)}: empty")
    missing = [k for k in KEY if k not in rows[0]]
    if missing:
        raise SystemExit(
            f"{os.path.basename(path)}: no {missing} column -- this looks like "
            "a pre-per-template export. Use the older script for it.")

    out = []
    for row in rows:
        members = num(row, "n_members")
        q = num(row, "q_onset_x_ms_auto_R")
        s = num(row, "s_end_x_ms_auto_R")
        t = num(row, "t_end_x_ms_auto_R")
        # All three required: a partial row would give an interval measured
        # against a missing endpoint.
        if q is None or s is None or t is None:
            continue
        rec = {"fit_model": model}
        rec.update({k: row.get(k, "").strip() for k in KEY})
        rec.update({
            "n_members": int(members) if members is not None else "",
            "bad_r": row.get("bad_r", "").strip(),
            "q_onset_ms": q,
            "s_end_ms": s,
            "qrs_ms": s - q,
            "t_end_ms": t,
            "qt_ms": t - q,
        })
        out.append(rec)
    return out, len(rows)


def summarize_by_model(rows):
    """mean / sd per (metric, model), across every row. Sample sd (n-1)."""
    groups = {}
    for r in rows:
        groups.setdefault(r["fit_model"], []).append(r)
    out = []
    for metric in SUMMARY_METRIC_ORDER:
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

    The model-sensitivity view, now per template rather than per bin: a large
    sd_ms means the choice of on/offset model moves THAT template's landmark.
    Grouping by bin alone would mix a bin's several morphologies together.
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
        got, total = read_model(path, model)
        rows.extend(got)
        print(f"{model:13s} {os.path.basename(path):55s} "
              f"{len(got):4d} of {total} rows usable")

    if not rows:
        raise SystemExit("no usable rows: no template had all of "
                         "q_onset / s_end / t_end populated")

    write_csv(os.path.join(outdir, "fit_model_per_row.csv"), rows,
              ["fit_model", *KEY, "n_members", "bad_r",
               "q_onset_ms", "s_end_ms", "qrs_ms", "t_end_ms", "qt_ms"])
    write_csv(os.path.join(outdir, "fit_model_summary.csv"),
              summarize_by_model(rows),
              ["metric", "Model", "n_rows",
               "Mean ms", "Std ms", "Min ms", "Max ms"])
    byrow = summarize_by_row(rows)
    write_csv(os.path.join(outdir, "fit_model_per_row_stats.csv"), byrow,
              [*KEY, "metric", "n_models", "n_members",
               "mean_ms", "sd_ms", "min_ms", "max_ms", "range_ms"])

    print()
    print("SUMMARY, across all rows")
    hdr = (f"{'metric':11s} {'Model':13s} {'n':>5s} {'Mean ms':>10s} "
           f"{'Std ms':>9s} {'Min ms':>10s} {'Max ms':>10s}")
    print(hdr)
    print("-" * len(hdr))
    for r in summarize_by_model(rows):
        sd = f"{r['Std ms']:9.4f}" if r["Std ms"] != "" else f"{'-':>9s}"
        print(f"{r['metric']:11s} {r['Model']:13s} {r['n_rows']:5d} "
              f"{r['Mean ms']:10.4f} {sd} {r['Min ms']:10.4f} {r['Max ms']:10.4f}")

    print()
    print(f"QRS / QT PER TEMPLATE, across the {len(MODELS)} models "
          "-- sd is model sensitivity")
    hdr = (f"{'bin':>4s} {'ch':4s} {'template':9s} {'metric':8s} {'n':>3s} "
           f"{'beats':>6s} {'mean':>10s} {'sd':>8s} {'range':>8s}")
    print(hdr)
    print("-" * len(hdr))
    for r in byrow:
        if r["metric"] not in ("qrs_ms", "qt_ms"):
            continue          # the two intervals; full table is in the CSV
        sd = f"{r['sd_ms']:8.3f}" if r["sd_ms"] != "" else f"{'-':>8s}"
        print(f"{r['bin_index']:>4s} {r['channel']:4s} {r['template']:9s} "
              f"{r['metric']:8s} {r['n_models']:3d} {str(r['n_members']):>6s} "
              f"{r['mean_ms']:10.3f} {sd} {r['range_ms']:8.3f}")


if __name__ == "__main__":
    main()
