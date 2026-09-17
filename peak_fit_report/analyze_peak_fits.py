#!/usr/bin/env python3
"""
Analyse PEAK LANDMARK POSITIONS across the peak fit models.

Third script in the set:
    compare_fit_models.py   on/offset models, the q_onset/s_end/t_end bars
                            and the QRS / QT intervals
    compare_peak_fits.py    peak models, p_peak plus those same intervals
    compare_peak_positions.py   <- this one

This one reports the landmarks the Fit-Peaks radio actually places, each as a
DISTANCE FROM THE R PEAK in ms:

    p_peak_from_r, q_peak_from_r, s_peak_from_r, t_peak_from_r

Negative means before R. The raw x_ms_auto_R columns are absolute positions
along each template, and each bin's R sits at a different column, so their sd
is dominated by that offset rather than by anything measured -- subtracting R
removes it and makes the numbers comparable to clinical ranges (P peak about
R-190 to R-100, Q peak R-40 to R-20, S peak R+20 to R+45, T peak R+180 to
R+320).

r_peak_abs_ms is carried alongside as the raw anchor, since R's own distance
from R is identically zero. The on/offset bars are deliberately not
here -- they belong to the other radio, and compare_peak_fits.py already covers
them for this model set.

Writes three CSVs:
    peak_pos_per_bin.csv        one row per (model, bin)
    peak_pos_summary.csv        mean / sd per (model, metric), across bins
    peak_pos_per_bin_stats.csv  mean / sd per (bin, metric), across the MODELS
                                -- how far the choice of peak model moves that
                                bin's landmark

NOTE ON FILENAMES
    The peak radio exports as 5_point / quadratic / cubic. 'cubic' is also an
    on/offset model with the same filename, so keep the peak exports in their
    own directory or the cubic column will be whichever file was written last.

Standard library only. Usage:
    python compare_peak_positions.py [input_dir] [output_dir]
"""

import csv
import os
import sys
from statistics import mean, stdev

MODELS = ["5_point", "quadratic", "cubic"]

CHANNELS = (1, 2, 3)   # scanned to find the populated one; not an output column

# metric name -> CSV column stem, for the landmarks reported as an offset from
# R. Order is the output column order.
OFFSET_LANDMARKS = [
    ("p_peak_from_r_ms", "p_peak"),
    ("q_peak_from_r_ms", "q_peak"),
    ("s_peak_from_r_ms", "s_peak"),
    ("t_peak_from_r_ms", "t_peak"),
]

# R's absolute position, kept as the anchor: its offset from itself is zero.
R_ABS = "r_peak_abs_ms"

METRICS = [m for m, _ in OFFSET_LANDMARKS] + [R_ABS]

# Row order in the summary CSV: metric-major, models alphabetical inside each
# metric -- same layout as the on/offset summary so the two read the same way.
SUMMARY_METRIC_ORDER = METRICS


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
    """Float from a cell, or None when absent/blank/unparseable."""
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
    """One row per bin. Returns (rows, channel_used).

    No channel column in the output, so exactly one channel may contribute --
    two leads would collide on the same bin_index and the per-bin statistics
    would silently average across them.

    R peak is REQUIRED -- every other value is measured relative to it, so a
    bin without R has no usable offsets and is skipped. Beyond that, a missing
    landmark is left as None and dropped from that metric's statistics only.
    """
    rows = list(csv.DictReader(open(path, newline="")))

    populated = [ch for ch in CHANNELS
                 if any(num(r, f"r_peak_ch{ch}_x_ms_auto_R") is not None
                        for r in rows)]
    if not populated:
        raise SystemExit(f"{os.path.basename(path)}: no channel has r_peak data")
    if len(populated) > 1:
        raise SystemExit(
            f"{os.path.basename(path)}: channels {populated} both populated. "
            "The output has no channel column, so narrow CHANNELS to one.")
    ch = populated[0]

    out = []
    for row in rows:
        r_abs = num(row, f"r_peak_ch{ch}_x_ms_auto_R")
        if r_abs is None:
            continue          # no anchor, no offsets
        rec = {
            "fit_model": model,
            "bin_index": row.get("bin_index", "").strip(),
        }
        for metric, stem in OFFSET_LANDMARKS:
            v = num(row, f"{stem}_ch{ch}_x_ms_auto_R")
            rec[metric] = None if v is None else v - r_abs
        rec[R_ABS] = r_abs
        out.append(rec)
    return out, ch


def stats(vals):
    """mean / sd / min / max / range for one metric's values. Sample sd (n-1)."""
    vals = [v for v in vals if v is not None]
    if not vals:
        return {"n": 0, "mean_ms": "", "sd_ms": "",
                "min_ms": "", "max_ms": "", "range_ms": ""}
    return {
        "n": len(vals),
        "mean_ms": round(mean(vals), 4),
        "sd_ms": round(stdev(vals), 4) if len(vals) > 1 else "",
        "min_ms": round(min(vals), 4),
        "max_ms": round(max(vals), 4),
        "range_ms": round(max(vals) - min(vals), 4),
    }


def summarize_by_model(rows):
    """mean / sd per (metric, model), taken ACROSS BINS.

    Metric-major: all three models for p_peak, then all three for q_peak, and
    so on, so one landmark's models read straight down the sheet.
    """
    groups = {}
    for r in rows:
        groups.setdefault(r["fit_model"], []).append(r)
    out = []
    for metric in SUMMARY_METRIC_ORDER:
        for model in sorted(groups):
            st = stats([r[metric] for r in groups[model]])
            out.append({
                "metric": metric,
                "Model": model,
                "Mean ms": st["mean_ms"],
                "Std ms": st["sd_ms"],
                "Min ms": st["min_ms"],
                "Max ms": st["max_ms"],
            })
    return out


def summarize_by_bin(rows):
    """mean / sd per (bin, metric), taken ACROSS THE PEAK MODELS.

    sd_ms is the model sensitivity of that landmark in that bin. An sd of
    exactly 0 in every bin means the radio never reached the landmark -- a
    plumbing result, not agreement between models.
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
            st = stats([r[metric] for r in recs])
            out.append({"bin_index": b, "metric": metric,
                        "n_models": st.pop("n"), **st})
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
        print(f"{model:10s} {os.path.basename(path):55s} "
              f"{len(got):4d} bins, ch{ch}")

    if not rows:
        raise SystemExit("no usable rows: no bin had any peak landmark populated")

    write_csv(os.path.join(outdir, "peak_pos_per_bin.csv"), rows,
              ["fit_model", "bin_index"] + METRICS)
    write_csv(os.path.join(outdir, "peak_pos_summary.csv"),
              summarize_by_model(rows),
              ["metric", "Model", "Mean ms", "Std ms", "Min ms", "Max ms"])
    bybin = summarize_by_bin(rows)
    write_csv(os.path.join(outdir, "peak_pos_per_bin_stats.csv"), bybin,
              ["bin_index", "metric", "n_models",
               "mean_ms", "sd_ms", "min_ms", "max_ms", "range_ms"])

    # ---- console view ------------------------------------------------------
    print()
    print("SUMMARY, across bins")
    hdr = (f"{'metric':18s} {'Model':10s} {'Mean ms':>10s} {'Std ms':>9s} "
           f"{'Min ms':>10s} {'Max ms':>10s}")
    print(hdr)
    print("-" * len(hdr))
    for r in summarize_by_model(rows):
        sd = f"{r['Std ms']:9.4f}" if r["Std ms"] != "" else f"{'-':>9s}"
        print(f"{r['metric']:18s} {r['Model']:10s} {r['Mean ms']:10.4f} {sd} "
              f"{r['Min ms']:10.4f} {r['Max ms']:10.4f}")

    print()
    print(f"P PEAK PER BIN (offset from R), across the {len(MODELS)} models")
    hdr = (f"{'bin':>4s} {'n':>3s} {'mean':>11s} {'sd':>8s} {'range':>8s}")
    print(hdr)
    print("-" * len(hdr))
    for r in bybin:
        if r["metric"] != "p_peak_ms":
            continue
        sd = f"{r['sd_ms']:8.3f}" if r["sd_ms"] != "" else f"{'-':>8s}"
        print(f"{r['bin_index']:>4s} {r['n_models']:3d} "
              f"{r['mean_ms']:11.3f} {sd} {r['range_ms']:8.3f}")

    # ---- did each landmark actually respond to the radio? ------------------
    # A landmark whose per-bin sd is 0 everywhere is not "stable", it is
    # unplumbed. Worth stating outright rather than leaving a column of zeros.
    print()
    print("model sensitivity, over all bins")
    for metric in METRICS:
        per = [r for r in bybin if r["metric"] == metric]
        moved = sum(1 for r in per if r["sd_ms"] not in ("", 0.0))
        sds = [r["sd_ms"] for r in per if r["sd_ms"] not in ("", 0.0)]
        tail = (f"  median sd {sorted(sds)[len(sds) // 2]:.3f} ms"
                if sds else "   <-- identical in every bin")
        print(f"  {metric:18s} moved in {moved:3d}/{len(per)} bins{tail}")


if __name__ == "__main__":
    main()
