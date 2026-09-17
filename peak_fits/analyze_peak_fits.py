#!/usr/bin/env python3
"""
Analyse PEAK LANDMARK POSITIONS across the peak fit models.

Each landmark is reported as a DISTANCE FROM THE R PEAK in ms:

    p_peak_from_r, q_peak_from_r, s_peak_from_r, t_peak_from_r

Negative means before R. The raw x_ms_auto_R columns are absolute positions
along a template, and every template has R at its own column, so their spread
is dominated by that offset rather than by anything measured. Subtracting R
removes it and makes the numbers comparable to clinical ranges (P peak about
R-200 to R-100, Q peak R-45 to R-20, S peak R+20 to R+50, T peak R+180 to
R+320). r_peak_abs_ms rides along as the raw anchor, since R's distance from
itself is zero.

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
    peak_pos_per_row.csv    one row per (model, bin, channel, template)
    peak_pos_summary.csv    mean / sd per (model, metric), across all rows
    peak_pos_per_row_stats.csv
                            mean / sd per (bin, channel, template, metric),
                            across the MODELS -- how far the choice of peak
                            model moves that one template's landmark

Columns read (R-aligned auto):
    p_peak_x_ms_auto_R, q_peak_x_ms_auto_R, s_peak_x_ms_auto_R,
    t_peak_x_ms_auto_R, r_peak_x_ms_auto_R

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
MODELS = ["5_point", "quadratic", "cubic"]

# metric -> CSV column stem, for the landmarks reported as an offset from R.
OFFSET_LANDMARKS = [
    ("p_peak_from_r_ms", "p_peak"),
    ("q_peak_from_r_ms", "q_peak"),
    ("s_peak_from_r_ms", "s_peak"),
    ("t_peak_from_r_ms", "t_peak"),
]
R_ABS = "r_peak_abs_ms"
METRICS = tuple(m for m, _ in OFFSET_LANDMARKS) + (R_ABS,)

# Metric-major, models alphabetical inside each metric -- same layout as the
# on/offset summary so the two read the same way.
SUMMARY_METRIC_ORDER = METRICS

KEY = ("bin_index", "channel", "template")


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

    R is REQUIRED -- every other value is measured relative to it. Beyond that
    a missing landmark stays None and drops out of that metric's statistics
    only; the five are independent positions, not interval endpoints.
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
    for row in rows:
        members = num(row, "n_members")
        r_abs = num(row, "r_peak_x_ms_auto_R")
        if r_abs is None:
            continue          # no anchor, no offsets
        rec = {"fit_model": model}
        rec.update({k: row.get(k, "").strip() for k in KEY})
        rec["n_members"] = int(members) if members is not None else ""
        rec["bad_r"] = row.get("bad_r", "").strip()
        for metric, stem in OFFSET_LANDMARKS:
            v = num(row, f"{stem}_x_ms_auto_R")
            rec[metric] = None if v is None else v - r_abs
        rec[R_ABS] = r_abs
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
        raise SystemExit("no usable rows: no template had an r_peak")

    write_csv(os.path.join(outdir, "peak_pos_per_row.csv"), rows,
              ["fit_model", *KEY, "n_members", "bad_r", *METRICS])
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
    print(f"P PEAK PER TEMPLATE (offset from R), across the {len(MODELS)} "
          "models -- sd is model sensitivity")
    hdr = (f"{'bin':>4s} {'ch':4s} {'template':9s} {'metric':17s} {'n':>3s} "
           f"{'beats':>6s} {'mean':>10s} {'sd':>8s} {'range':>8s}")
    print(hdr)
    print("-" * len(hdr))
    for r in byrow:
        if r["metric"] != "p_peak_from_r_ms":
            continue          # full table is in the CSV
        sd = f"{r['sd_ms']:8.3f}" if r["sd_ms"] != "" else f"{'-':>8s}"
        print(f"{r['bin_index']:>4s} {r['channel']:4s} {r['template']:9s} "
              f"{r['metric']:17s} {r['n_models']:3d} {str(r['n_members']):>6s} "
              f"{r['mean_ms']:10.3f} {sd} {r['range_ms']:8.3f}")

    # A landmark whose per-row sd is 0 everywhere is not "stable", it is
    # unplumbed -- the radio never reached the code that places it. Worth
    # stating outright rather than leaving a column of zeros to be read as
    # agreement between the models.
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
