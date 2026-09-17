#!/usr/bin/env python3
"""
Compare auto-detected ECG landmark timings across on/offset fit models.

Reads the *_template_markings.csv exports (one per fit model), pulls the
R-aligned AUTO landmark positions, and derives the two intervals:

    QRS = s_end - q_onset
    QT  = t_end - q_onset

Writes three CSVs:
    fit_model_per_bin.csv   one row per (model, bin)
    fit_model_summary.csv   mean / sd per (model, metric), across bins
    fit_model_per_bin_stats.csv
                            mean / sd per (bin, metric), across the 5 MODELS --
                            i.e. how far the choice of fit model moves that bin

Columns read, for channel N in 1..3:
    q_onset_chN_x_ms_auto_R
    s_end_chN_x_ms_auto_R
    t_end_chN_x_ms_auto_R

Standard library only. Usage:
    python compare_fit_models.py [input_dir] [output_dir]
"""

import csv
import os
import sys
from statistics import mean, stdev

# Model name -> filename prefix. Add to this if more exports appear.
MODELS = ["linear", "sigmoid", "cubic", "cubic_spline", "frac_poly"]

CHANNELS = (1, 2, 3)   # scanned to find the populated one; not an output column
METRICS = ("q_onset_ms", "s_end_ms", "qrs_ms", "t_end_ms", "qt_ms")

# Row order in the summary CSV: metric-major, positions before intervals, with
# the models alphabetical inside each metric. Grouping by metric is what makes
# the five models directly comparable down a column.
SUMMARY_METRIC_ORDER = ("q_onset_ms", "s_end_ms", "t_end_ms", "qrs_ms", "qt_ms")


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

    Blank cells are normal here: an export only populates the channels that
    were marked, so ch2/ch3 are often empty for every bin.
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
    """One row per bin that has all three landmarks present.

    No channel column in the output, so exactly one channel may contribute --
    otherwise two channels' rows would collide on the same bin_index and the
    per-bin statistics would silently average across leads. The populated
    channel is detected and reported; more than one is an error rather than a
    guess.
    """
    rows = list(csv.DictReader(open(path, newline="")))

    populated = []
    for ch in CHANNELS:
        if any(num(r, f"q_onset_ch{ch}_x_ms_auto_R") is not None for r in rows):
            populated.append(ch)
    if not populated:
        raise SystemExit(f"{os.path.basename(path)}: no channel has q_onset data")
    if len(populated) > 1:
        raise SystemExit(
            f"{os.path.basename(path)}: channels {populated} both populated. "
            "The output has no channel column, so pick one by editing CHANNELS.")
    ch = populated[0]

    out = []
    for row in rows:
        q = num(row, f"q_onset_ch{ch}_x_ms_auto_R")
        s = num(row, f"s_end_ch{ch}_x_ms_auto_R")
        t = num(row, f"t_end_ch{ch}_x_ms_auto_R")
        # All three required: a partial row would give an interval computed
        # against a missing endpoint.
        if q is None or s is None or t is None:
            continue
        out.append({
            "fit_model": model,
            "bin_index": row.get("bin_index", "").strip(),
            "q_onset_ms": q,
            "s_end_ms": s,
            "qrs_ms": s - q,
            "t_end_ms": t,
            "qt_ms": t - q,
        })
    return out, ch


def summarize_by_model(rows):
    """mean / sd per (metric, model), taken ACROSS BINS. Sample sd (n-1).

    Metric-major: all five models for q_onset, then all five for s_end, and so
    on, so one metric's models read straight down the sheet.
    """
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
                "Mean ms": round(mean(vals), 4) if vals else "",
                "Std ms": round(stdev(vals), 4) if len(vals) > 1 else "",
                "Min ms": round(min(vals), 4) if vals else "",
                "Max ms": round(max(vals), 4) if vals else "",
            })
    return out


def summarize_by_bin(rows):
    """mean / sd per (bin, metric), taken ACROSS THE FIT MODELS.

    This is the model-sensitivity view: a large sd_ms here means the choice of
    on/offset model moves that bin's landmark a lot, which is the quantity the
    five exports exist to compare. Bins are ordered numerically where the index
    parses as an integer.
    """
    groups = {}
    for r in rows:
        groups.setdefault(r["bin_index"], []).append(r)

    def key(b):
        try:
            return (0, int(b))
        except ValueError:
            return (1, 0)

    out = []
    for b in sorted(groups, key=key):
        recs = groups[b]
        for metric in METRICS:
            vals = [r[metric] for r in recs if r[metric] is not None]
            out.append({
                "bin_index": b,
                "metric": metric,
                "n_models": len(vals),
                "mean_ms": round(mean(vals), 4) if vals else "",
                "sd_ms": round(stdev(vals), 4) if len(vals) > 1 else "",
                "min_ms": round(min(vals), 4) if vals else "",
                "max_ms": round(max(vals), 4) if vals else "",
                "range_ms": round(max(vals) - min(vals), 4) if vals else "",
            })
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
        got, ch = read_model(path, model)
        rows.extend(got)
        print(f"{model:13s} {os.path.basename(path):55s} "
              f"{len(got):4d} bins, ch{ch}")

    if not rows:
        raise SystemExit("no usable rows: no bin had all of "
                         "q_onset / s_end / t_end populated")

    write_csv(os.path.join(outdir, "fit_model_per_bin.csv"), rows,
              ["fit_model", "bin_index",
               "q_onset_ms", "s_end_ms", "qrs_ms", "t_end_ms", "qt_ms"])
    write_csv(os.path.join(outdir, "fit_model_summary.csv"),
              summarize_by_model(rows),
              ["metric", "Model", "Mean ms", "Std ms", "Min ms", "Max ms"])
    write_csv(os.path.join(outdir, "fit_model_per_bin_stats.csv"),
              summarize_by_bin(rows),
              ["bin_index", "metric", "n_models",
               "mean_ms", "sd_ms", "min_ms", "max_ms", "range_ms"])

    # Console view, so a run is readable without opening a file.
    print()
    print("SUMMARY, across bins")
    hdr = (f"{'metric':11s} {'Model':13s} {'Mean ms':>10s} {'Std ms':>9s} "
           f"{'Min ms':>10s} {'Max ms':>10s}")
    print(hdr); print("-" * len(hdr))
    for r in summarize_by_model(rows):
        sd = f"{r['Std ms']:9.4f}" if r["Std ms"] != "" else f"{'-':>9s}"
        print(f"{r['metric']:11s} {r['Model']:13s} {r['Mean ms']:10.4f} {sd} "
              f"{r['Min ms']:10.4f} {r['Max ms']:10.4f}")

    print()
    print("PER BIN, across the 5 models -- sd is model sensitivity")
    bybin = summarize_by_bin(rows)
    hdr = (f"{'bin':>4s} {'metric':11s} {'n':>3s} {'mean':>10s} "
           f"{'sd':>8s} {'range':>8s}")
    print(hdr); print("-" * len(hdr))
    for r in bybin:
        if r["metric"] not in ("qrs_ms", "qt_ms"):
            continue          # the two intervals; full table is in the CSV
        sd = f"{r['sd_ms']:8.3f}" if r["sd_ms"] != "" else f"{'-':>8s}"
        print(f"{r['bin_index']:>4s} {r['metric']:11s} {r['n_models']:3d} "
              f"{r['mean_ms']:10.3f} {sd} {r['range_ms']:8.3f}")


if __name__ == "__main__":
    main()
