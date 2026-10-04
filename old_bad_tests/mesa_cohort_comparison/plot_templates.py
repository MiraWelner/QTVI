#!/usr/bin/env python3
"""
plot_templates_by_stage.py

Same load + R-aligned overlay look as plot_templates.py, but EACH PATIENT gets
FIVE subplots -- Wake / N1 / N2 / N3 / REM -- instead of one. A patient's bins are
split by sleep stage (read from the source .bin), and each stage's bins are
R-aligned and overlaid (ECG solid + PPG dashed) exactly as before.

Stage codes come from the .bin written by file_to_bin.cpp (with its 5->4
remap, so 4 == REM in this data):
    0 -> Wake, 1 -> N1, 2 -> N2, 3 -> N3, 4 -> REM   (-1/other -> skipped)

Pairing: a patient's .bin lives IN THE SAME FOLDER as its template CSV. For
4013204_20110418_template.csv the bin is 4013204_20110418.bin in that folder
(label = CSV stem minus "_template").

The alignment (align_r: horizontal R/upslope anchor + vertical PR/foot DC
baseline) and windowing are IDENTICAL to plot_templates.py -- only the bin
grouping (by stage) and the subplot layout (3 per patient) are new.

Run:  python plot_templates_by_stage.py
Output: perpatient_by_stage_<tag>.png  (rows = patients, cols = Wake/N1/N2/N3/REM)
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import struct
import warnings

# Benign: an all-NaN column in an overlay matrix makes np.nanmedian warn. The
# result (NaN) is fine to plot as a gap; silence just this warning.
warnings.filterwarnings("ignore", message="All-NaN slice encountered")
warnings.filterwarnings("ignore", message="Mean of empty slice")

# ---- Folders with *_template.csv. Edit these. ----
FOLDERS = [
    Path("templates_healthy"),
    Path("templates_untreated"),
]

# ---- Where the source .bin files live: IN the same FOLDERS as the templates.
# A template 4013204_20110418_template.csv pairs with 4013204_20110418.bin in
# the same folder (label = CSV stem minus "_template"). ----
BIN_SUFFIX = ".bin"          # <label> + BIN_SUFFIX  ->  the .bin filename

# cfg.bin_size_minutes used when the templates were built (bin duration).
BIN_SIZE_MINUTES = 5.0

VALUE_SUFFIX = "_Normalized_r"     # column to align/plot
BIN_COL = "bin_num"

# ---- Absolute RMS rejection thresholds (no Tukey). A bin is dropped if its
# metric exceeds the cutoff in EITHER channel. Units are in the normalized
# signal's amplitude. THESE DEFAULTS ARE PLACEHOLDERS -- set PRINT_RMS_STATS
# below to see your data's actual min/median/max per stage, then put each
# threshold just above the clean-bin cluster and below the noisy ones.
NOISE_RMS_MAX = 0.15    # HF noise: RMS of the waveform's 2nd difference
SHAPE_RMS_MAX = 0.30    # shape:    RMS distance of a bin to the group median

# Set True to print, per stage, the distribution of both RMS metrics across
# bins (min / median / max) so you can choose the two thresholds above. Does
# not change plotting -- just prints. Turn off once thresholds are set.
PRINT_RMS_STATS = False

# Normalized y-axis window shared by every plot.
Y_LIM = (-0.7, 2.0)
# X-axis window in samples relative to the anchor.
XLIM = (-300, 500)

# .bin layout constants (from file_to_bin.hpp; the "=142/=568" comments there
# are stale -- the formula gives 146 fields / 584 bytes).
NUM_CHANNELS = 36
NUM_HEADER_FIELDS = 1 + 4 * NUM_CHANNELS + 1     # 146
HEADER_SIZE = NUM_HEADER_FIELDS * 4              # 584

_X_PERIOD_S = None
_X_UNITS = "samples"

STAGE_ORDER = ["wake", "n1", "n2", "n3", "rem"]
STAGE_LABEL = {"wake": "Wake", "n1": "N1", "n2": "N2", "n3": "N3", "rem": "REM"}


def read_bin_stages(bin_path):
    """Read (epoch_seconds, stages[]) from a .bin. The stages array is the last
    thing in the file; its offset is validated two independent ways (forward
    sum of channel sizes vs. from-EOF) so a wrong layout errors instead of
    silently returning garbage. Returns (None, None) if the file has no stages."""
    data = bin_path.read_bytes()
    if len(data) < HEADER_SIZE:
        raise ValueError(f"{bin_path.name}: too small for a {HEADER_SIZE}-byte header")
    ints = struct.unpack_from(f"<{NUM_HEADER_FIELDS}I", data, 0)
    epoch_seconds = float(ints[0])
    sizes_up = ints[1:1 + NUM_CHANNELS]
    sizes_raw = ints[1 + NUM_CHANNELS:1 + 2 * NUM_CHANNELS]
    sleep_size = ints[1 + 4 * NUM_CHANNELS]
    if epoch_seconds <= 0 or sleep_size == 0:
        return None, None

    stages_bytes = sleep_size * 8
    eof_offset = len(data) - stages_bytes
    body = HEADER_SIZE
    for ch in range(NUM_CHANNELS):
        body += sizes_up[ch] * 8
        body += sizes_raw[ch] * 2 * 8
    if body != eof_offset:
        raise ValueError(
            f"{bin_path.name}: stages offset mismatch "
            f"(forward={body}, fromEOF={eof_offset}); refusing to read.")

    stages = np.frombuffer(data, dtype="<f8", count=sleep_size, offset=eof_offset)
    return epoch_seconds, np.asarray(stages, dtype=float)


def _stage_group(code):
    """Map a numeric stage code to one of the five sleep states, or None to
    skip. Codes (as written by file_to_bin.cpp, 5->4 remap so 4 == REM):
        0 -> Wake, 1 -> N1, 2 -> N2, 3 -> N3, 4 -> REM."""
    return {0: "wake", 1: "n1", 2: "n2", 3: "n3", 4: "rem"}.get(code)


def bin_stage_group(bin_index, bin_seconds, epoch_seconds, stages):
    """Stage group for a bin ONLY if every sleep epoch it overlaps maps to the
    same group -- i.e. the bin is pure. Mixed bins (spanning a stage
    transition, or containing any unknown/-1 epoch) return None and are
    dropped, so a stage's template is never blended with beats from another
    stage. A 5-min bin spans ~10 epochs of 30 s, so this drops transition bins."""
    kind, g = _classify_bin(bin_index, bin_seconds, epoch_seconds, stages)
    return g if kind == "stage" else None


def _classify_bin(bin_index, bin_seconds, epoch_seconds, stages):
    """Classify a bin by the sleep epochs it overlaps:
        ("stage", name)  -- every overlapped epoch is the SAME known state
        ("mixed", None)  -- overlaps >1 known state (a transition bin)
        ("unknown", None)-- overlaps an unknown/-1 epoch, or none at all
    'mixed' is reported separately so it can be counted per patient."""
    t0 = bin_index * bin_seconds
    t1 = (bin_index + 1) * bin_seconds
    e0 = max(0, int(t0 // epoch_seconds))
    e1 = min(len(stages), int(np.ceil(t1 / epoch_seconds)))
    if e1 <= e0:
        return ("unknown", None)
    groups = set()
    for e in range(e0, e1):
        g = _stage_group(stages[e])
        if g is None:
            return ("unknown", None)     # -1 / unrecognized epoch
        groups.add(g)
    if len(groups) == 1:
        return ("stage", next(iter(groups)))
    return ("mixed", None)               # spans more than one known state


def _detect_period_s(df):
    """Set the module x-axis scale from a dataframe's x_ms column (once)."""
    global _X_PERIOD_S, _X_UNITS
    if _X_PERIOD_S is not None:
        return
    if "x_ms" in df.columns:
        xm = pd.to_numeric(df["x_ms"], errors="coerce").to_numpy(dtype=float)
        d = np.diff(xm)
        d = d[np.isfinite(d) & (d > 0)]       # positive = within-bin step (skip resets)
        if d.size:
            _X_PERIOD_S = float(np.median(d)) / 1000.0
            _X_UNITS = "s"
            return
    _X_PERIOD_S = 1.0                          # fallback: x stays in samples
    _X_UNITS = "samples"


def _xlabel():
    return ("seconds relative to anchor" if _X_UNITS == "s"
            else "samples relative to anchor")


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

def align_r(beats, return_keep=False):
    """alignment.hpp Pass 1 (horizontal anchor) + Pass 3 (vertical DC baseline).
    Horizontal: anchor every beat's R (ECG) / 50%-upslope (PPG) column at the
    max anchor column, integer-shift + NaN-pad. Vertical: match each beat's
    pre-anchor baseline mean (PR segment for ECG, pre-foot for PPG) to the
    reference (median-length) beat's, via a constant vertical shift. Aligns
    every input beat -- no scalar Tukey; outlier rejection is the RMS shape pass.
    Returns (matrix [n_beats x shared_w], anchor_col). If return_keep=True, also
    returns keep_idx: the indices into the INPUT `beats` in matrix-row order (all
    of them now, since no beat is rejected here) -- so a caller can map each
    aligned row back to its source bin for the shape-rejection pass."""
    orig = beats                                   # keep the input for index mapping
    if not beats:
        return (None, None, []) if return_keep else (None, None)

    # No scalar Tukey here: outlier rejection is done downstream by the RMS
    # shape-distance pass (_shape_outlier_bins), which is the only Python-side
    # rejection. This function just aligns every input beat. is_ppg is still
    # needed to pick the vertical-baseline strategy below.
    is_ppg = any(b.get("foot_col", -1) >= 0 for b in beats)
    if is_ppg:
        for b in beats:
            b["_up50"] = float(b["r_col"])          # up50 column = anchor for PPG

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

    if return_keep:
        survivors = {id(b) for b in beats}     # kept dicts (identity-stable)
        keep_idx = [i for i, b in enumerate(orig) if id(b) in survivors]
        return mat, anchor, keep_idx
    return mat, anchor

def beats_for_file(df, chan, is_ecg, value_suffix):
    """Same as plot_templates.py's beats_for_file, but each beat dict also
    carries its 'bin_num' so beats can be split by sleep stage before align_r.
    All bins except the final 5; ECG flipped upright; PPG records foot col."""
    vcol = f"{chan}{value_suffix}"
    if vcol not in df.columns:
        return []
    _detect_period_s(df)
    subs = split_on_bin_change(df)
    subs = subs[:-5] if len(subs) > 5 else []
    beats = []
    for sub in subs:
        sub = sub.reset_index(drop=True)
        bnum = int(sub[BIN_COL].iloc[0])
        y = pd.to_numeric(sub[vcol], errors="coerce").to_numpy(dtype=float)
        if np.all(np.isnan(y)):
            continue
        r_col = anchor_index(sub, is_ecg, y)
        base = np.nanmedian(y)
        if is_ecg and (y[r_col] - base) < 0:
            y = base - (y - base)
        foot_col = -1
        sys_amp = np.nan
        if not is_ecg:
            foot_col, peak_col, _ = _ppg_first_pulse(y)
            if foot_col >= 0 and peak_col >= 0:
                sys_amp = float(y[peak_col] - y[foot_col])
        beats.append({"y": y, "r_col": r_col, "foot_col": foot_col,
                      "length": int(np.sum(~np.isnan(y))),
                      "r_amp": float(abs(y[r_col] - base)),
                      "sys_amp": sys_amp, "bin_num": bnum})
    return beats


def _windowed(mat, anchor):
    if mat is None:
        return None, None
    xs = np.arange(mat.shape[1]) - anchor
    msk = (xs >= XLIM[0]) & (xs <= XLIM[1])
    return xs[msk] / 1000.0, mat[:, msk]


def _find_bin_for(label, folder):
    """Locate the .bin for a patient, which lives in the SAME folder as the
    template CSV and is named <label>.bin (label = CSV stem minus '_template').
    E.g. 4013204_20110418_template.csv  ->  4013204_20110418.bin. Falls back to
    a <label>*.bin glob in that folder for minor naming variations."""
    direct = folder / (label + BIN_SUFFIX)
    if direct.is_file():
        return direct
    hits = sorted(folder.glob(f"{label}*{BIN_SUFFIX}"))
    return hits[0] if hits else None


def _split_beats_by_stage(beats, bin2stage):
    """Partition a beat list into the five sleep states by each beat's bin_num."""
    out = {g: [] for g in STAGE_ORDER}
    for bt in beats:
        g = bin2stage.get(bt["bin_num"])
        if g in out:
            out[g].append(bt)
    return out


def _row_rms_to_median(mat):
    """Per-row RMS distance to the column-wise median of the matrix, over
    columns where both the row and the median are finite. Returns an array of
    length n_rows (NaN if a row shares no finite columns with the median)."""
    if mat is None or mat.shape[0] == 0:
        return np.array([])
    with np.errstate(all="ignore"):
        med = np.nanmedian(mat, axis=0)
    out = np.full(mat.shape[0], np.nan)
    for i in range(mat.shape[0]):
        row = mat[i]
        m = ~np.isnan(row) & ~np.isnan(med)
        if m.sum() >= 3:
            d = row[m] - med[m]
            out[i] = float(np.sqrt(np.mean(d * d)))
    return out


def _shape_outlier_bins(mat, anchor, keep_idx, beats, thresh):
    """Return bin_num values whose shape RMS (distance to the group median)
    exceeds an ABSOLUTE threshold. No Tukey/IQR -- a plain RMS cutoff, so it
    works even when most bins are bad (unlike a relative outlier test)."""
    if mat is None or not keep_idx:
        return set()
    rms = _row_rms_to_median(mat)
    bad = set()
    for row_i, src_i in enumerate(keep_idx):
        v = rms[row_i]
        if not np.isnan(v) and v > thresh:
            bad.add(beats[src_i]["bin_num"])
    return bad


def _drop_bins_from_aligned(mat, keep_idx, beats, drop_bins):
    """Remove rows whose source bin_num is in drop_bins. Returns a new matrix."""
    if mat is None:
        return None
    rows = [row_i for row_i, src_i in enumerate(keep_idx)
            if beats[src_i]["bin_num"] not in drop_bins]
    if not rows:
        return None
    return mat[rows, :]


def _row_hf_noise(mat):
    """Per-row high-frequency noise metric: RMS of the discrete 2nd difference
    (a simple high-pass -- smooth signal -> near 0, HF noise -> large), computed
    within finite samples only. Returns an array of length n_rows (NaN if a row
    has too few finite samples)."""
    if mat is None or mat.shape[0] == 0:
        return np.array([])
    out = np.full(mat.shape[0], np.nan)
    for i in range(mat.shape[0]):
        row = mat[i]
        finite = row[~np.isnan(row)]
        if finite.size >= 5:
            d2 = np.diff(finite, n=2)          # 2nd difference: HF emphasis
            out[i] = float(np.sqrt(np.mean(d2 * d2)))
    return out


def _noise_outlier_bins(mat, keep_idx, beats, thresh):
    """Return bin_num values whose HF-noise RMS (RMS of the 2nd difference)
    exceeds an ABSOLUTE threshold. No Tukey -- a plain cutoff, so noisy bins
    are dropped even when noise is the norm (e.g. Wake)."""
    if mat is None or not keep_idx:
        return set()
    hf = _row_hf_noise(mat)
    bad = set()
    for row_i, src_i in enumerate(keep_idx):
        v = hf[row_i]
        if not np.isnan(v) and v > thresh:
            bad.add(beats[src_i]["bin_num"])
    return bad


def run_by_stage(value_suffix, tag):
    """One ROW per patient, five COLUMNS (Wake/N1/N2/N3/REM). Each cell is the
    R-aligned ECG (solid) + PPG (dashed) overlay for that patient's bins in
    that stage -- identical alignment/window to plot_templates.py."""
    palette_group = {}
    palette = ["tab:green", "tab:red", "tab:blue", "tab:purple"]

    # patients = list of (group, label, {stage: (ecgmat,ecganc,ppgmat,ppganc)})
    patients = []

    for gi, folder in enumerate(FOLDERS):
        if not folder.is_dir():
            print(f"[skip] {folder} not a directory", file=sys.stderr)
            continue
        palette_group[folder.name] = palette[gi % len(palette)]
        files = sorted(p for p in folder.glob("*_template.csv")
                       if "_template_mark" not in p.name and "_template_snips" not in p.name)
        if not files:
            print(f"[warn] no *_template.csv in {folder}")
            continue
        print(f"\n=== [{tag}] {folder}  ({len(files)} file(s)) ===")
        for path in files:
            label = path.stem.replace("_template", "")
            binp = _find_bin_for(label, folder)
            if binp is None:
                print(f"  [skip] {label}: no {label}{BIN_SUFFIX} in {folder}")
                continue
            try:
                epoch_s, stages = read_bin_stages(binp)
            except Exception as e:
                print(f"  [error] {binp.name}: {e}")
                continue
            if stages is None:
                print(f"  [skip] {label}: .bin has no sleep stages")
                continue
            try:
                df = pd.read_csv(path)
            except Exception as e:
                print(f"  [error] {path.name}: {e}")
                continue

            ecg_beats = beats_for_file(df, "ch1", True, value_suffix)
            ppg_beats = beats_for_file(df, "ppg", False, value_suffix)

            # Eligible bins = those that survived beats_for_file's last-5 drop
            # (union across ECG/PPG). Classify each: pure stage / mixed / unknown.
            bin_seconds = BIN_SIZE_MINUTES * 60.0
            eligible_bins = sorted({b["bin_num"] for b in ecg_beats}
                                   | {b["bin_num"] for b in ppg_beats})
            bin2stage = {}
            n_mixed = 0
            for bnum in eligible_bins:
                kind, g = _classify_bin(bnum, bin_seconds, epoch_s, stages)
                if kind == "stage":
                    bin2stage[bnum] = g
                elif kind == "mixed":
                    n_mixed += 1
                # "unknown" (-1 epochs) is neither assigned nor counted as mixed

            ecg_by = _split_beats_by_stage(ecg_beats, bin2stage)
            ppg_by = _split_beats_by_stage(ppg_beats, bin2stage)

            per_stage = {}
            rms_dropped = {}          # stage -> bins dropped by RMS shape
            noise_dropped = {}        # stage -> bins dropped for HF noise
            kept = {}                 # stage -> bins kept in the overlay
            any_stage = False
            for g in STAGE_ORDER:
                eb, pb = ecg_by[g], ppg_by[g]
                # Align each channel (anchor only; no scalar Tukey), keeping the
                # map from matrix rows back to source bins for rejection.
                em, ea, e_keep = align_r(eb, return_keep=True) if eb else (None, None, [])
                pm, pa, p_keep = align_r(pb, return_keep=True) if pb else (None, None, [])

                if PRINT_RMS_STATS:
                    def _stats(arr):
                        f = arr[~np.isnan(arr)] if arr.size else arr
                        if f.size == 0:
                            return "n/a"
                        return f"min={f.min():.4f} med={np.median(f):.4f} max={f.max():.4f} (n={f.size})"
                    if em is not None:
                        print(f"    [{g}] ECG noise {_stats(_row_hf_noise(em))} | "
                              f"shape {_stats(_row_rms_to_median(em))}")
                    if pm is not None:
                        print(f"    [{g}] PPG noise {_stats(_row_hf_noise(pm))} | "
                              f"shape {_stats(_row_rms_to_median(pm))}")

                # (1) HF-noise rejection: drop bins whose high-frequency content
                # (RMS of 2nd difference) exceeds an ABSOLUTE cutoff in EITHER
                # channel -- the messy/noisy bins, removed entirely.
                noisy = set()
                noisy |= _noise_outlier_bins(em, e_keep, eb, NOISE_RMS_MAX)
                noisy |= _noise_outlier_bins(pm, p_keep, pb, NOISE_RMS_MAX)
                noise_dropped[g] = len(noisy)

                # (2) Shape rejection: drop bins whose RMS distance to the group
                # median exceeds an ABSOLUTE cutoff in EITHER channel. Union with
                # the noisy set; drop all from BOTH overlays.
                bad = set(noisy)
                bad |= _shape_outlier_bins(em, ea, e_keep, eb, SHAPE_RMS_MAX)
                bad |= _shape_outlier_bins(pm, pa, p_keep, pb, SHAPE_RMS_MAX)
                rms_dropped[g] = len(bad) - len(noisy)   # shape-only, excl. noise

                if bad:
                    em = _drop_bins_from_aligned(em, e_keep, eb, bad)
                    pm = _drop_bins_from_aligned(pm, p_keep, pb, bad)

                # Bins kept in this stage = all bins assigned to it minus the
                # dropped (noise + shape) set.
                stage_bins = {b["bin_num"] for b in eb} | {b["bin_num"] for b in pb}
                kept[g] = len(stage_bins - bad)

                per_stage[g] = (em, ea, pm, pa)
                if em is not None or pm is not None:
                    any_stage = True
            if any_stage:
                patients.append((folder.name, label, per_stage, n_mixed,
                                 rms_dropped, noise_dropped, kept))
            n_by = {g: len(ecg_by[g]) for g in STAGE_ORDER}
            print(f"  {label}: bins/stage ECG {n_by}  | mixed dropped={n_mixed}"
                  f"  | noise dropped={ {g: noise_dropped[g] for g in STAGE_ORDER if noise_dropped[g]} }"
                  f"  | rms dropped={ {g: rms_dropped[g] for g in STAGE_ORDER if rms_dropped[g]} }"
                  f"  | kept={ {g: kept[g] for g in STAGE_ORDER if kept[g]} }")

    if not patients:
        print(f"[skip] [{tag}]: no patients")
        return

    patients.sort(key=lambda e: (list(palette_group).index(e[0]), e[1]))

    nrows = len(patients)
    ncols = len(STAGE_ORDER)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(3.6 * ncols, 2.4 * nrows),
                             squeeze=False, sharex=True)

    for ri, (group, label, per_stage, n_mixed, rms_dropped, noise_dropped, kept) in enumerate(patients):
        gcolor = palette_group[group]
        for ci, stg in enumerate(STAGE_ORDER):
            ax = axes[ri][ci]
            em, ea, pm, pa = per_stage[stg]

            xe, me = _windowed(em, ea)
            if me is not None:
                al = max(0.12, min(0.5, 25.0 / max(1, me.shape[0])))
                for row in me:
                    ax.plot(xe, row, color=gcolor, alpha=al, linewidth=0.4)
                with np.errstate(all="ignore"):
                    ax.plot(xe, np.nanmedian(me, axis=0), color=gcolor, linewidth=1.6)
            ax.set_ylim(*Y_LIM)
            ax.tick_params(axis="y", labelsize=6, labelcolor=gcolor)

            xp, mp = _windowed(pm, pa)
            if mp is not None:
                ax2 = ax.twinx()
                ap = max(0.12, min(0.5, 25.0 / max(1, mp.shape[0])))
                for row in mp:
                    ax2.plot(xp, row, color=gcolor, alpha=ap, linewidth=0.4)
                with np.errstate(all="ignore"):
                    ax2.plot(xp, np.nanmedian(mp, axis=0), color=gcolor,
                             linewidth=1.4, linestyle="--")
                ax2.set_ylim(*Y_LIM)
                ax2.tick_params(axis="y", labelsize=6, labelcolor=gcolor)

            ax.axvline(0, color="k", linewidth=0.5, linestyle="--", alpha=0.4)
            ax.tick_params(axis="x", labelbottom=True, labelsize=6)
            # Each subplot title: stage name + bins kept, and bins dropped this
            # stage for HF noise and for RMS shape.
            title = (f"{STAGE_LABEL.get(stg, stg)}  "
                     f"(kept: {kept.get(stg, 0)}, noise: {noise_dropped.get(stg, 0)}, "
                     f"rms: {rms_dropped.get(stg, 0)})")
            ax.set_title(title, fontsize=8)
            if ci == 0:
                # Row label: patient + how many bins were dropped as mixed
                # (spanning >1 sleep stage) for this patient.
                ax.set_ylabel(f"{group}\n{label}\nmixed drop: {n_mixed}", fontsize=6)

    fig.suptitle("per-patient R-aligned ECG (solid) + PPG (dashed) by sleep stage",
                 fontsize=11)
    fig.supxlabel(_xlabel(), fontsize=9)
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    out_path = Path(f"perpatient_by_stage_{tag}.png")
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[ok] [{tag}]: {nrows} patients x {ncols} stages -> {out_path.name}")


def main():
    run_by_stage(VALUE_SUFFIX, "normalized")


if __name__ == "__main__":
    main()
