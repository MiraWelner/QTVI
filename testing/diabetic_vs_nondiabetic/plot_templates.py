#!/usr/bin/env python3
"""
align_templates_r.py

Applies alignment.hpp's horizontal + vertical alignment (R-aligned pass) to the
bin templates in each folder, for ECG (ch1/ch2/ch3) and PPG, then plots the
aligned overlay + median per channel per folder.

WHICH PART OF alignment.hpp THIS IS
-----------------------------------
alignment.hpp's beat-level slice + Tukey rejection already ran inside the C++
(each bin IS one median beat). The two steps that remain meaningful on the
bin templates are the two ALIGNMENT steps, ported here exactly:

  * Pass 1 -- horizontal R-align (extract_beats_and_align, lines ~212-237):
      anchor every beat's R column at R_anchor = max R column over survivors,
      then integer-shift each beat right by (R_anchor - r_col), NaN-padding.
      No resampling.

  * Pass 3 -- vertical PR-baseline DC align (lines ~239-275):
      reference beat = first beat whose length == median length. Its PR-segment
      mean (a small window just before R) is the target; every beat is shifted
      vertically so its own PR mean matches. Window, per the C++:
          pr_w   = max(3, median_length // 20)
          pr_gap = max(1, median_length // 50)
          pr_lo  = max(0, R_anchor - pr_gap - pr_w)
          pr_hi  = max(pr_lo, R_anchor - pr_gap)
      baseline = mean of non-NaN samples in [pr_lo, pr_hi); needs >= 3 samples.

Here each "beat" is one bin template and "length" is that bin's non-NaN sample
count. R column comes from the r_peak_ch<c>_location_autodetect_r one-hot flag
(PPG anchors on ppg_onset ... _r, i.e. the foot, matching the PPG path's
foot/upslope alignment intent).

Edit the constants below, then run:  python align_templates_r.py
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---- Folders to process. Edit these. ----
FOLDERS = [
    Path("templates_healthy"),
    Path("templates_untreated"),
]

# R-aligned pass columns. Only ch1 and ppg are aligned/processed.
ECG_CHANNELS = ["ch1"]
PULSE_CHANNELS = ["ppg"]
VALUE_SUFFIX = "_Normalized_r"          # what we align/plot
BIN_COL = "bin_num"

# Anchor flag columns (R-pass, autodetect).
def ecg_anchor_col(chan):   # chan = "ch1" -> "1"
    return f"r_peak_ch{chan[-1]}_location_autodetect_r"

def pulse_anchor_col(chan):
    return f"{chan}_onset_location_autodetect_r"   # foot


def split_on_bin_change(df):
    changed = df[BIN_COL].ne(df[BIN_COL].shift())
    return [g for _, g in df.groupby(changed.cumsum(), sort=False)]


def _ppg_first_pulse(v, rate=1000.0):
    """Find the FIRST pulse's foot, systolic peak, and 50%-upslope crossing.

    A _template.csv bin spans several pulses, so a global argmax lands on a
    random pulse and makes up50/foot come from different beats. Instead we
    locate the first systolic peak (max within the first ~1.2 s), its preceding
    foot (min before it), then the half-height crossing between them.
    Returns (foot, peak, up50) indices into v, or (-1,-1,-1) if none."""
    n = len(v)
    if np.isfinite(v).sum() < 3:
        return -1, -1, -1
    win = min(n, int(round(1.2 * rate)))
    seg = np.where(np.isnan(v[:win]), -np.inf, v[:win])
    peak = int(np.argmax(seg))
    pre = np.where(np.isnan(v[:peak + 1]), np.inf, v[:peak + 1])
    foot = int(np.argmin(pre)) if peak > 0 else 0
    fv, pv = v[foot], v[peak]
    if not (np.isfinite(fv) and np.isfinite(pv)) or pv <= fv:
        return foot, peak, peak
    half = fv + 0.5 * (pv - fv)
    up50 = peak
    for k in range(foot + 1, peak + 1):
        if np.isnan(v[k]) or np.isnan(v[k - 1]):
            continue
        if v[k] >= half > v[k - 1]:
            up50 = k
            break
    return foot, peak, up50


def anchor_index(sub, is_ecg, value):
    """R (ECG, argmax|y-baseline|) or first-pulse 50%-upslope crossing (PPG)."""
    v = value
    finite = ~np.isnan(v)
    if not finite.any():
        return 0
    if is_ecg:
        baseline = np.nanmedian(v)
        dev = np.where(finite, np.abs(v - baseline), -np.inf)
        return int(np.argmax(dev))
    _, _, up50 = _ppg_first_pulse(v)
    return up50 if up50 >= 0 else 0


def collect_beats(df, chan, is_ecg):
    """Return list of dicts {y, r_col, length} for each GOOD bin on this
    channel. A bin is rejected if its detected anchor is not a genuine sharp
    peak (ratio of the anchor deviation to the median deviation is small),
    which is how the corrupt/flat bins get excluded before alignment."""
    vcol = f"{chan}{VALUE_SUFFIX}"
    if vcol not in df.columns:
        return []
    beats = []
    for sub in split_on_bin_change(df):
        sub = sub.reset_index(drop=True)
        y = pd.to_numeric(sub[vcol], errors="coerce").to_numpy(dtype=float)
        if np.all(np.isnan(y)):
            continue
        r_col = anchor_index(sub, is_ecg, y)

        # Reject bins whose "R" is really a noise spike / flat trace.
        base = np.nanmedian(y)
        dev = np.abs(y - base)
        med_dev = np.nanmedian(dev)
        ratio = (dev[r_col] / med_dev) if med_dev > 0 else 0.0
        if is_ecg and ratio < 8.0:
            continue                        # no sharp R -> corrupt/flat bin

        # Flip inverted-R beats upright so upright and inverted beats don't
        # cancel in the pooled median. alignment.hpp's r_peak() reports the
        # deflection sign and the pipeline orients the trace upright; here the
        # sign of (y[R] - baseline) tells us the polarity.
        if is_ecg and (y[r_col] - base) < 0:
            y = base - (y - base)          # reflect about baseline -> upright R
            dev = np.abs(y - base)         # (unchanged in magnitude, but keep consistent)

        length = int(np.sum(~np.isnan(y)))
        r_amp = float(dev[r_col])           # R height above baseline (for Tukey)
        beats.append({"y": y, "r_col": r_col, "length": length, "r_amp": r_amp})
    return beats


def _tukey_keep(beats, key, k):
    """Keep beats whose beats[i][key] is within [Q1-k*IQR, Q3+k*IQR].
    Faithful port of alignment.hpp keep_within_tukey."""
    vals = np.array([b[key] for b in beats], dtype=float)
    finite = vals[~np.isnan(vals)]
    if finite.size < 4:
        return beats
    q1, q3 = np.quantile(finite, 0.25), np.quantile(finite, 0.75)
    iqr = q3 - q1
    if iqr <= 0.0:
        return beats
    lo_b, hi_b = q1 - k * iqr, q3 + k * iqr
    return [b for b, v in zip(beats, vals)
            if not np.isnan(v) and lo_b <= v <= hi_b]


def align_r(beats):
    """alignment.hpp Pass 1 (horizontal anchor) + Pass 3 (vertical DC baseline).
    Horizontal: anchor every beat's R (ECG) / foot (PPG) column at the max
    anchor column, integer-shift + NaN-pad. Vertical: match each beat's
    pre-anchor baseline mean (PR segment for ECG, pre-foot for PPG) to the
    reference (median-length) beat's, via a constant vertical shift. Keeps the
    Tukey rejection + upright-flip from collect_beats.
    Returns (matrix [n_beats x shared_w], anchor_col)."""
    if not beats:
        return None, None

    beats = _tukey_keep(beats, key="length", k=1.5)
    # PPG rejects on systolic amplitude (peak-foot) and up50 position, as
    # alignment.hpp does; ECG rejects on R amplitude.
    is_ppg = any(b.get("foot_col", -1) >= 0 for b in beats)
    if is_ppg:
        for b in beats:
            b["_up50"] = float(b["r_col"])          # up50 column = anchor for PPG
        beats = _tukey_keep(beats, key="sys_amp", k=1.5)
        beats = _tukey_keep(beats, key="_up50", k=3.0)
    else:
        beats = _tukey_keep(beats, key="r_amp", k=1.5)
    if not beats:
        return None, None

    # ---- Pass 1: horizontal align on a shared axis ----
    anchor = max(b["r_col"] for b in beats)            # max anchor column
    max_tail = max(len(b["y"]) - b["r_col"] for b in beats)
    shared_w = anchor + max_tail

    mat = np.full((len(beats), shared_w), np.nan)
    foot_cols_shifted = []
    for i, b in enumerate(beats):
        prepend = anchor - b["r_col"]                  # >= 0 by construction
        y = b["y"]
        n = min(prepend + len(y), shared_w) - prepend
        if n > 0:
            mat[i, prepend:prepend + n] = y[:n]
        fc = b.get("foot_col", -1)
        foot_cols_shifted.append(prepend + fc if fc >= 0 else -1)

    # ---- reference beat = first beat whose length == median length ----
    lens = sorted(b["length"] for b in beats)
    median_length = lens[len(lens) // 2]
    ref_idx = next((i for i, b in enumerate(beats)
                    if b["length"] == median_length), -1)

    is_ppg = any(fc >= 0 for fc in foot_cols_shifted)

    # ---- vertical DC baseline alignment ----
    if ref_idx >= 0:
        if is_ppg:
            # alignment.hpp Pass 2: match each beat's mean over a +/- fb_w
            # window around its OWN foot column to the reference beat's.
            fb_w = max(1, median_length // 50)

            def baseline(i):
                fc = foot_cols_shifted[i]
                if fc < 0:
                    return np.nan
                lo = max(0, fc - fb_w)
                hi = min(shared_w, fc + fb_w + 1)
                seg = mat[i, lo:hi]
                seg = seg[~np.isnan(seg)]
                return np.mean(seg) if seg.size >= 1 else np.nan
        else:
            # alignment.hpp Pass 3 (ECG): PR-segment window just before R.
            pr_w = max(3, median_length // 20)
            pr_gap = max(1, median_length // 50)
            pr_lo = max(0, anchor - pr_gap - pr_w)
            pr_hi = max(pr_lo, anchor - pr_gap)

            def baseline(i):
                seg = mat[i, pr_lo:min(pr_hi, shared_w)]
                seg = seg[~np.isnan(seg)]
                return np.mean(seg) if seg.size >= 3 else np.nan

        target = baseline(ref_idx)
        if not np.isnan(target):
            for i in range(mat.shape[0]):
                b_base = baseline(i)
                if np.isnan(b_base):
                    continue
                d = b_base - target
                if d != 0.0:
                    mat[i] = mat[i] - d        # constant shift; NaNs stay NaN

    return mat, anchor


def plot_group(chan, mat, R_anchor, folder_name, out_dir, is_ecg):
    if mat is None:
        return None
    x = np.arange(mat.shape[1]) - R_anchor   # samples relative to anchor
    m = (x >= -100) & (x <= 500)
    x, mat = x[m], mat[:, m]
    fig, ax = plt.subplots(figsize=(10, 4))
    # Higher alpha + thin lines so individual sharp beats stay visible rather
    # than washing into a smooth-looking band at low alpha.
    n = mat.shape[0]
    a = max(0.15, min(0.6, 30.0 / max(1, n)))
    for i in range(n):
        ax.plot(x, mat[i], color="tab:blue", alpha=a, linewidth=0.5)
    with np.errstate(all="ignore"):
        med = np.nanmedian(mat, axis=0)
    ax.plot(x, med, color="tab:red", linewidth=2.0,
            label=f"median (n={mat.shape[0]})")
    ax.axvline(0, color="k", linewidth=0.7, linestyle="--", alpha=0.5)
    anchor = "R" if is_ecg else "foot"
    ax.set_title(f"{chan}: {folder_name} -- aligned (H+V), anchor={anchor}")
    ax.set_xlabel(f"samples relative to {anchor}")
    ax.set_ylabel("normalized")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    out_path = out_dir / f"aligned_r_{chan}_{folder_name}.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path


def beats_for_file(df, chan, is_ecg, value_suffix):
    """All bins from one file except the final 5, flipped upright. Returns the
    beat list ready for align_r. For PPG, also records the foot column so the
    vertical DC step can match foot baselines the way alignment.hpp does."""
    vcol = f"{chan}{value_suffix}"
    if vcol not in df.columns:
        return []
    bins = split_on_bin_change(df)
    bins = bins[:-5] if len(bins) > 5 else []
    beats = []
    for sub in bins:
        sub = sub.reset_index(drop=True)
        y = pd.to_numeric(sub[vcol], errors="coerce").to_numpy(dtype=float)
        if np.all(np.isnan(y)):
            continue
        r_col = anchor_index(sub, is_ecg, y)
        base = np.nanmedian(y)
        if is_ecg and (y[r_col] - base) < 0:
            y = base - (y - base)
        # PPG foot + systolic amplitude from the SAME first pulse as the up50
        # anchor (so vertical foot-baseline and up50 horizontal use one pulse).
        foot_col = -1
        sys_amp = np.nan
        if not is_ecg:
            foot_col, peak_col, _ = _ppg_first_pulse(y)
            if foot_col >= 0 and peak_col >= 0:
                sys_amp = float(y[peak_col] - y[foot_col])
        beats.append({"y": y, "r_col": r_col, "foot_col": foot_col,
                      "length": int(np.sum(~np.isnan(y))),
                      "r_amp": float(abs(y[r_col] - base)),
                      "sys_amp": sys_amp})
    return beats


def run_suffix(value_suffix, tag):
    """Build per-patient plots for one value column suffix (e.g. '_Normalized_r'
    or '_raw_mv_r'). Each patient gets ONE subplot with ch1 (ECG, left axis,
    aligned on R) and ppg (right axis, aligned on foot) overlaid. Writes
    perpatient_ecg_ppg_<tag>.png."""
    out_root = Path(".")
    group_color = {}
    palette = ["tab:green", "tab:red", "tab:blue", "tab:purple"]

    # patients = list of (group, label, ecg_mat, ecg_anchor, ppg_mat, ppg_anchor)
    patients = []

    for gi, folder in enumerate(FOLDERS):
        if not folder.is_dir():
            print(f"[skip] {folder} is not a directory", file=sys.stderr)
            continue
        group_color[folder.name] = palette[gi % len(palette)]
        files = sorted(p for p in folder.glob("*_template.csv")
                       if "_template_mark" not in p.name and "_template_snips" not in p.name)
        if not files:
            print(f"[warn] no *_template.csv in {folder}")
            continue
        print(f"\n=== [{tag}] {folder}  ({len(files)} file(s)) ===")
        for path in files:
            try:
                df = pd.read_csv(path)
            except Exception as e:
                print(f"  [error] {path.name}: {e}")
                continue

            ecg_beats = beats_for_file(df, "ch1", True, value_suffix)
            ppg_beats = beats_for_file(df, "ppg", False, value_suffix)
            ecg_mat, ecg_anchor = align_r(ecg_beats) if ecg_beats else (None, None)
            ppg_mat, ppg_anchor = align_r(ppg_beats) if ppg_beats else (None, None)
            if ecg_mat is None and ppg_mat is None:
                continue
            patients.append((folder.name, path.stem.replace("_template", ""),
                             ecg_mat, ecg_anchor, ppg_mat, ppg_anchor))

    if not patients:
        print(f"[skip] [{tag}]: no patients")
        return

    patients.sort(key=lambda e: (list(group_color).index(e[0]), e[1]))

    def windowed(mat, anchor):
        if mat is None:
            return None, None
        x = np.arange(mat.shape[1]) - anchor
        msk = (x >= -100) & (x <= 500)
        return x[msk], mat[:, msk]

    n = len(patients)
    ncols = min(5, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(3.4 * ncols, 2.4 * nrows),
                             squeeze=False, sharex=True)
    for idx in range(nrows * ncols):
        ax = axes[idx // ncols][idx % ncols]
        if idx >= n:
            ax.axis("off")
            continue
        group, label, ecg_mat, ecg_anchor, ppg_mat, ppg_anchor = patients[idx]
        gcolor = group_color[group]

        # ECG on the left axis (group color).
        xe, me = windowed(ecg_mat, ecg_anchor)
        if me is not None:
            ae = max(0.12, min(0.5, 25.0 / max(1, me.shape[0])))
            for row in me:
                ax.plot(xe, row, color=gcolor, alpha=ae, linewidth=0.4)
            with np.errstate(all="ignore"):
                ax.plot(xe, np.nanmedian(me, axis=0), color=gcolor, linewidth=1.6)
        ax.set_ylabel("ecg", color=gcolor, fontsize=6)
        ax.tick_params(axis="y", labelcolor=gcolor, labelsize=6)

        # PPG on a twin right axis, same group color as ECG (distinguished by
        # being on the right axis + dashed median).
        xp, mp = windowed(ppg_mat, ppg_anchor)
        if mp is not None:
            ax2 = ax.twinx()
            ap = max(0.12, min(0.5, 25.0 / max(1, mp.shape[0])))
            for row in mp:
                ax2.plot(xp, row, color=gcolor, alpha=ap, linewidth=0.4)
            with np.errstate(all="ignore"):
                ax2.plot(xp, np.nanmedian(mp, axis=0), color=gcolor,
                         linewidth=1.4, linestyle="--")
            ax2.set_ylabel("ppg", color=gcolor, fontsize=6)
            ax2.tick_params(axis="y", labelcolor=gcolor, labelsize=6)

        ax.axvline(0, color="k", linewidth=0.5, linestyle="--", alpha=0.4)
        ax.set_title(f"{group}\n{label}", fontsize=6)
        ax.tick_params(axis="x", labelsize=6)

    unit = "normalized" if "Normalized" in value_suffix else "raw mV"
    fig.suptitle(f"per-patient ECG (solid, anchor=R) + PPG (dashed, anchor=50%-upslope) "
                 f"[{tag}, {unit}] -- {', '.join(group_color)}", fontsize=11)
    fig.supxlabel("samples relative to anchor", fontsize=9)
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))

    out_path = out_root / f"perpatient_ecg_ppg_{tag}.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[ok] [{tag}]: {n} patients -> {out_path.name}")


def run_group_overlay(folder, value_suffix, tag):
    """Overlay ALL bins from ALL files in one folder into a single figure
    (ECG solid on left axis, PPG dashed on right axis, both aligned). Writes
    group_overlay_<folder>_<tag>.png."""
    if not folder.is_dir():
        print(f"[skip] {folder} is not a directory", file=sys.stderr)
        return
    files = sorted(p for p in folder.glob("*_template.csv")
                   if "_template_mark" not in p.name and "_template_snips" not in p.name)
    if not files:
        print(f"[warn] no *_template.csv in {folder}")
        return

    # Pool every bin across every file in this folder.
    ecg_beats, ppg_beats = [], []
    for path in files:
        try:
            df = pd.read_csv(path)
        except Exception as e:
            print(f"  [error] {path.name}: {e}")
            continue
        ecg_beats.extend(beats_for_file(df, "ch1", True, value_suffix))
        ppg_beats.extend(beats_for_file(df, "ppg", False, value_suffix))

    ecg_mat, ecg_anchor = align_r(ecg_beats) if ecg_beats else (None, None)
    ppg_mat, ppg_anchor = align_r(ppg_beats) if ppg_beats else (None, None)
    if ecg_mat is None and ppg_mat is None:
        print(f"[skip] {folder.name} [{tag}]: no data")
        return

    def windowed(mat, anchor):
        if mat is None:
            return None, None
        x = np.arange(mat.shape[1]) - anchor
        msk = (x >= -100) & (x <= 500)
        return x[msk], mat[:, msk]

    color = "tab:green" if "healthy" in folder.name.lower() else "tab:red"
    fig, ax = plt.subplots(figsize=(10, 5))

    xe, me = windowed(ecg_mat, ecg_anchor)
    if me is not None:
        ae = max(0.05, min(0.4, 15.0 / max(1, me.shape[0])))
        for row in me:
            ax.plot(xe, row, color=color, alpha=ae, linewidth=0.3)
        with np.errstate(all="ignore"):
            ax.plot(xe, np.nanmedian(me, axis=0), color=color, linewidth=2.0,
                    label=f"ECG median (n={me.shape[0]})")
    ax.set_ylabel("ecg", color=color)
    ax.tick_params(axis="y", labelcolor=color)

    xp, mp = windowed(ppg_mat, ppg_anchor)
    if mp is not None:
        ax2 = ax.twinx()
        ap = max(0.05, min(0.4, 15.0 / max(1, mp.shape[0])))
        for row in mp:
            ax2.plot(xp, row, color=color, alpha=ap, linewidth=0.3)
        with np.errstate(all="ignore"):
            ax2.plot(xp, np.nanmedian(mp, axis=0), color=color, linewidth=2.0,
                     linestyle="--", label=f"PPG median (n={mp.shape[0]})")
        ax2.set_ylabel("ppg", color=color)
        ax2.tick_params(axis="y", labelcolor=color)

    ax.axvline(0, color="k", linewidth=0.6, linestyle="--", alpha=0.4)
    unit = "normalized" if "Normalized" in value_suffix else "raw mV"
    ax.set_title(f"{folder.name}: all bins overlaid -- ECG (solid, anchor=R) + "
                 f"PPG (dashed, anchor=50%-upslope) [{unit}]")
    ax.set_xlabel("samples relative to anchor")
    fig.tight_layout()

    out_path = Path(".") / f"group_overlay_{folder.name}_{tag}.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[ok] group overlay {folder.name} [{tag}] -> {out_path.name}")


def main():
    run_suffix("_Normalized_r", "normalized")
    run_suffix("_raw_mv_r", "raw")
    # Four group-overlay PNGs: {each folder} x {normalized, raw}.
    for folder in FOLDERS:
        run_group_overlay(folder, "_Normalized_r", "normalized")
        run_group_overlay(folder, "_raw_mv_r", "raw")


if __name__ == "__main__":
    main()
