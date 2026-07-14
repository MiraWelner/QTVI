"""
Templates-only combined overlay (no snips).

For each cohort folder, produce ONE PNG: a single axes with all channels
(CH1, CH2, CH3, PPG, ABP, ART, ART_PULM) on one shared millisecond timeline,
one normalized MEAN curve per subject (averaged across that subject's bins).

Data source: <subject>_templates.csv (written by the viewer). Each row is one
sample of one bin's template. Columns:
    x_ms                     shared master timeline (0 at ch1 R anchor)
    <c>_mv                   the channel's template value
    <c>_x_peak_ms            that channel's own-anchor axis (0 at its peak)
for c in ch1, ch2, ch3, ppg, abp, art, art_pulm.

Method
------
1. For each subject, for each channel, resample every bin's template onto a
   common OWN-ANCHOR ms grid (0 at that channel's peak, via <c>_x_peak_ms),
   then average the bins -> one mean curve per channel per subject.
2. Normalize (per subject), matching the Global Baseline spec:
     ECG : Feature / median_bins( |R_peak| + |S_peak| )
     PPG : Feature_local_ratio / median_bins( PI ),
           PI = 100*(sys - dia)/|dia|,  local_ratio = 100*(y - dia)/|dia|
   PPG-like channels are DC-centered, so we add +1 before PI/ratio so the
   diastolic trough is positive (keeps PI stable); this only rescales PI.
3. Place each channel's mean curve on the shared timeline using the per-bin
   ECG->channel offset from x_ms (IQR-trimmed median, rejecting beat-mispaired
   bins). x=0 is re-zeroed to the left edge of the visible window.

Output PNGs are written next to this script, named after this script.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from pathlib import Path
from collections import defaultdict

# ======================= EDIT THESE =======================
# Common own-anchor grid: step (ms) and half-width (ms) each side of the
# anchor. Bins are resampled onto [-GRID_HALF_MS, +GRID_HALF_MS] at GRID_STEP_MS
# before averaging. (The templates are already in ms via *_x_peak_ms.)
GRID_STEP_MS = 2.0
GRID_HALF_MS = 800.0

# Visible x-window in MILLISECONDS on the shared timeline (left edge -> 0).
X_WINDOW_MS = (-200.0, 1200.0)

# Physiological band (ms) for the R -> PPG systolic-peak offset. The PPG
# template spans ~2 cycles, so we pick the systolic peak whose offset from R
# lands in this window (rejecting the wrong-cycle peak that sits ~one RR off).
# A real R->peak is ~250-450 ms; the band is padded for HR variation.
PULSE_MIN_MS = 100.0
PULSE_MAX_MS = 600.0
# A re-picked peak within this many ms of a band edge is treated as "no real
# peak in band" (the true peak is outside the window and the sample just
# clamped to the boundary) -> the bin is skipped as invalid.
PULSE_EDGE_TOL_MS = 30.0
# Minimum peak prominence, as a fraction of the pulse's own amplitude above
# baseline, for a candidate to count as a real systolic peak.
PULSE_MIN_PROMINENCE = 0.30

# Half-width (ms) of the single cardiac cycle kept around the correct PPG
# peak when averaging. The raw PPG template spans ~2 cycles; we keep only
# +-PULSE_HALF_MS around the re-picked systolic peak so the average is one
# clean pulse. ~450 ms each side covers a full cycle down to ~65 bpm.
PULSE_HALF_MS = 450.0

# Empirical scale applied to the final pulse shift before placement. 0.5
# halves it. Set to 1.0 to disable.
SHIFT_SCALE = 0.5

# Where to find each subject's _templates.csv. None => auto-search the cohort
# folder and nearby csv_for_analysis folders.
TEMPLATES_DIR = None

# Cohort folders. The plot label for each is its folder name. All cohorts
# are drawn on ONE shared plot, colored by cohort (see COHORT_COLORS); ECG and
# PPG share the cohort's color.
COHORTS = [
    Path(r"D:\USERS\MiraWelner\QTVI\testing\diabetic_vs_nondiabetic\BL_DMU"),
    Path(r"D:\USERS\MiraWelner\QTVI\testing\diabetic_vs_nondiabetic\bl_healthy"),
]
# Colors assigned to cohorts in order (first folder red, second blue, ...).
COHORT_COLORS = ["#c0202a", "#1f4fb0", "#0a7a3c", "#d07a15"]
# ==========================================================

CHAN_ORDER = ["CH1", "CH2", "CH3", "PPG", "ABP", "ART", "ART_PULM"]
ECG_CHANS  = {"CH1", "CH2", "CH3"}
PPG_LIKE   = {"PPG", "ABP", "ART", "ART_PULM"}

TMPL_COL = {
    "CH1": "ch1", "CH2": "ch2", "CH3": "ch3",
    "PPG": "ppg", "ABP": "abp", "ART": "art", "ART_PULM": "art_pulm",
}
CHAN_COLOR = {
    "CH1": "#1f4fb0", "CH2": "#2f74d0", "CH3": "#4aa3e0",
    "PPG": "#b01020", "ABP": "#0a7a3c", "ART": "#7a3fb0", "ART_PULM": "#d07a15",
}

SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR  = SCRIPT_PATH.parent
SCRIPT_STEM = SCRIPT_PATH.stem

GRID = np.arange(-GRID_HALF_MS, GRID_HALF_MS + GRID_STEP_MS, GRID_STEP_MS)


def find_templates_csv(subj_id, cohort_folder):
    name = f"{subj_id}_templates.csv"
    cands = []
    if TEMPLATES_DIR is not None:
        cands.append(Path(TEMPLATES_DIR) / name)
    cands += [
        cohort_folder / name,
        cohort_folder / "csv_for_analysis" / name,
        cohort_folder.parent / "csv_for_analysis" / name,
    ]
    for c in cands:
        if c.exists():
            return c
    return None


def pulse_anchor_correction(g, col, ecg_anchor):
    """The PPG-like template spans ~two cardiac cycles, so the CSV's anchor
    (argmax) sometimes lands on the SECOND cycle's peak -- inflating the
    offset by ~one RR. Re-pick the correct systolic peak as the HIGHEST-
    amplitude sample whose shared-axis offset from the ECG R falls in the
    physiological band [PULSE_MIN_MS, PULSE_MAX_MS]. (The signal is noisy, so
    we take the max within the band, not the first local bump.)
    Returns the correction (ms) to ADD to that channel's own-anchor axis so 0
    sits on the correct-cycle peak (~0 when the argmax was already correct),
    or None if the band is empty for this bin."""
    xpk = g[f"{col}_x_peak_ms"].to_numpy(dtype=float)
    xms = g["x_ms"].to_numpy(dtype=float)
    y = g[f"{col}_mv"].to_numpy(dtype=float)
    ok = np.isfinite(xpk) & np.isfinite(xms) & np.isfinite(y)
    if ok.sum() < 5:
        return None
    xpk, xms, y = xpk[ok], xms[ok], y[ok]

    off = xms - ecg_anchor                     # R -> each sample (ms)
    band = (off >= PULSE_MIN_MS) & (off <= PULSE_MAX_MS)
    if not band.any():
        return None
    idx_band = np.where(band)[0]
    best = idx_band[np.argmax(y[idx_band])]

    # Reject a non-peak: the winner must be a genuine INTERIOR local maximum
    # (not the signal still rising/falling through the band edge) with
    # positive prominence above the bin's own baseline. Bins whose only
    # in-band maximum sits at the band edge -- i.e. the real systolic peak is
    # outside [PULSE_MIN_MS, PULSE_MAX_MS], so the sample just clamps to the
    # boundary -- return None, and the first-valid-bin search skips them.
    off_best = off[best]
    if off_best <= PULSE_MIN_MS + PULSE_EDGE_TOL_MS or \
       off_best >= PULSE_MAX_MS - PULSE_EDGE_TOL_MS:
        return None                             # peak pinned to band edge
    # Local-maximum check on the full (x_ms-sorted) trace around `best`.
    order = np.argsort(xms)
    xs, ys = xms[order], y[order]
    pos = int(np.searchsorted(xs, xms[best]))
    pos = min(max(pos, 1), len(ys) - 2)
    if not (ys[pos] >= ys[pos - 1] and ys[pos] >= ys[pos + 1]):
        return None                             # not a local max (on a slope)
    # Prominence: the peak must rise above the bin's baseline (median) by a
    # fraction of the pulse's own amplitude, else it's baseline noise.
    baseline = float(np.nanmedian(y))
    amp = float(np.nanmax(y) - baseline)
    if amp <= 0 or (y[best] - baseline) < PULSE_MIN_PROMINENCE * amp:
        return None
    return float(xpk[best])                     # its position on the x_peak axis


def bin_curve_on_grid(g, col, x_correction=0.0, window_half_ms=None):
    """Resample one bin's template (own-anchor axis) onto the common GRID.
    x_correction (ms) is subtracted from the own-anchor axis first, so a
    re-picked correct-cycle peak becomes 0 on the GRID (keeps the averaged
    shape centered on the same physiological peak the shift uses).
    window_half_ms (if set) blanks the grid beyond +-window_half_ms of the
    peak, so only ONE cardiac cycle around the correct peak survives -- this
    removes the second cycle the raw PPG template carries, giving a clean
    single-pulse average.
    Returns a GRID-length array (NaN where out of range), or None if empty."""
    x = g[f"{col}_x_peak_ms"].to_numpy(dtype=float) - x_correction
    y = g[f"{col}_mv"].to_numpy(dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return None
    x, y = x[ok], y[ok]
    order = np.argsort(x)
    x, y = x[order], y[order]
    # linear interpolation onto the grid; outside the bin's own span -> NaN
    out = np.interp(GRID, x, y, left=np.nan, right=np.nan)
    if window_half_ms is not None:
        out[np.abs(GRID) > window_half_ms] = np.nan   # keep one cycle only
    return out


def ecg_ref_from_bins(bin_curves):
    # median over bins of (|R| + |S|) = (|max| + |min|) of each bin curve
    vals = []
    for c in bin_curves:
        m = np.isfinite(c)
        if m.any():
            vals.append(abs(np.nanmax(c)) + abs(np.nanmin(c)))
    return float(np.median(vals)) if vals else np.nan


def ppg_ref_from_bins(bin_curves):
    # median over bins of PI = 100*(sys-dia)/|dia|, on +1-shifted curve
    pis = []
    for c in bin_curves:
        m = np.isfinite(c)
        if not m.any():
            continue
        cc = c[m] + 1.0
        sys_pk, dia_tr = np.max(cc), np.min(cc)
        if dia_tr != 0:
            pis.append(100.0 * (sys_pk - dia_tr) / abs(dia_tr))
    return float(np.median(pis)) if pis else np.nan


def subject_channel_mean(tmpl_df, chan, ecg_col=None):
    """FIRST-VALID-BIN curve for one channel (no averaging). Iterates bins in
    order and returns the first one that yields a usable template for this
    channel, normalized by that template's own values. PPG-like channels get
    the correct-cycle peak re-pick + single-cycle window, same as before.
    Returns (curve_on_grid, None, 1) -- se is None (single template, no band)
    -- or (None, None, 0) if no valid bin."""
    col = TMPL_COL.get(chan)
    if col is None:
        return None, None, 0
    if f"{col}_mv" not in tmpl_df.columns or f"{col}_x_peak_ms" not in tmpl_df.columns:
        return None, None, 0

    is_pulse = chan in PPG_LIKE
    epk = f"{ecg_col}_x_peak_ms" if ecg_col else None

    for _, g in tmpl_df.groupby("bin_num"):
        corr = 0.0
        if is_pulse:
            # A pulse bin is only valid if it has a real in-band systolic peak.
            if not (epk and epk in tmpl_df.columns):
                continue
            ec = (g["x_ms"] - g[epk]).dropna()
            if not len(ec):
                continue
            c2 = pulse_anchor_correction(g, col, float(ec.iloc[0]))
            if c2 is None:
                continue   # no valid peak in this bin -> skip to next bin
            corr = c2
        c = bin_curve_on_grid(g, col, x_correction=corr,
                              window_half_ms=(PULSE_HALF_MS if is_pulse else None))
        if c is None or not np.isfinite(c).any():
            continue   # empty / bad bin for this channel -> try next

        # Normalize by THIS template's own values (no cross-bin median).
        if chan in ECG_CHANS:
            ref = abs(np.nanmax(c)) + abs(np.nanmin(c))     # |R| + |S|
            if not np.isfinite(ref) or ref == 0:
                continue
            curve = c / ref
        else:
            cc = c + 1.0                                    # DC-centered shift
            dia = np.nanmin(cc)
            sys_pk = np.nanmax(cc)
            if not np.isfinite(dia) or dia == 0:
                continue
            pi = 100.0 * (sys_pk - dia) / abs(dia)          # this pulse's PI
            if not np.isfinite(pi) or pi == 0:
                continue
            local_ratio = 100.0 * (cc - dia) / abs(dia)
            curve = local_ratio / pi

        return curve, None, 1

    return None, None, 0


def pick_ecg_col(tmpl_df):
    for chan in ("CH1", "CH2", "CH3"):
        col = TMPL_COL[chan]
        xpk = f"{col}_x_peak_ms"
        if xpk in tmpl_df.columns and np.isfinite(
                (tmpl_df["x_ms"] - tmpl_df[xpk]).to_numpy(dtype=float)).any():
            return col
    return None


def channel_shift_vs_ecg(tmpl_df, chan, ecg_col):
    """Shift (ms) placing this channel's anchor on the shared timeline
    relative to the ECG anchor, using the RE-PICKED correct-cycle systolic
    peak per bin (median across bins).

    The PPG-like template spans ~2 cardiac cycles, so the CSV's argmax anchor
    sometimes sits on the 2nd cycle's peak, inflating the offset by ~one RR.
    Per bin we re-pick the systolic peak whose R-offset is physiological
    (PULSE_MIN_MS..PULSE_MAX_MS) and use that peak. This both fixes the
    'doubled' bins and rescues clamp-floored bins (they still contain a
    correct-cycle peak), so no bins are dropped."""
    col = TMPL_COL.get(chan)
    if col is None or ecg_col is None:
        return np.nan
    xpk, epk = f"{col}_x_peak_ms", f"{ecg_col}_x_peak_ms"
    mv = f"{col}_mv"
    for need in ("x_ms", "bin_num", xpk, epk, mv):
        if need not in tmpl_df.columns:
            return np.nan

    is_pulse = chan in PPG_LIKE
    for _, g in tmpl_df.groupby("bin_num"):
        ch = (g["x_ms"] - g[xpk]).dropna()   # argmax-peak shared position
        ec = (g["x_ms"] - g[epk]).dropna()   # ECG-R shared position
        if not (len(ch) and len(ec)):
            continue
        ecg_anchor = float(ec.iloc[0])
        argmax_off = float(ch.iloc[0] - ec.iloc[0])   # R -> argmax peak

        if not is_pulse:
            return argmax_off   # first valid bin

        corr = pulse_anchor_correction(g, col, ecg_anchor)
        if corr is None:
            continue            # no physiological peak in this bin -> try next
        # R -> correct peak = (R -> argmax peak) + corr.
        shift = argmax_off + corr
        return shift * SHIFT_SCALE

    return np.nan


def collect_cohort_curves(folder):
    """Gather per-subject, per-channel mean curves for one cohort folder.
    Returns (label, curves) where label is the folder name and curves[chan]
    is a list of (subj_id, x_shared_ms, mean, se)."""
    label = folder.name
    if not folder.exists():
        print(f"folder does not exist: {folder}")
        return label, {}

    tmpl_files = sorted(folder.glob("*_templates.csv"))
    if not tmpl_files and TEMPLATES_DIR is not None:
        tmpl_files = sorted(Path(TEMPLATES_DIR).glob("*_templates.csv"))
    print(f"\n==== cohort '{label}' : {len(tmpl_files)} subject(s) ====")

    curves = defaultdict(list)
    for tmpl_path in tmpl_files:
        subj_id = tmpl_path.name[: -len("_templates.csv")]
        try:
            tmpl_df = pd.read_csv(tmpl_path)
        except Exception as e:
            print(f"  [{subj_id}] read failed: {e}")
            continue

        ecg_col = pick_ecg_col(tmpl_df)
        if ecg_col is None:
            print(f"  [{subj_id}] SKIP: no ECG anchor")
            continue

        for chan in CHAN_ORDER:
            mean, se, nb = subject_channel_mean(tmpl_df, chan, ecg_col)
            if mean is None:
                continue
            shift = channel_shift_vs_ecg(tmpl_df, chan, ecg_col)
            if not np.isfinite(shift):
                continue
            curves[chan].append((subj_id, GRID + shift, mean, se))
            print(f"  [{subj_id}] {chan}: bins={nb}, shift={shift:+.1f} ms")

    return label, curves


def main():
    x_lo, x_hi = X_WINDOW_MS
    x0 = x_lo                          # left edge of window -> 0

    # Collect every cohort first.
    cohorts = []   # list of (label, color, curves)
    for i, folder in enumerate(COHORTS):
        label, curves = collect_cohort_curves(folder)
        if curves:
            color = COHORT_COLORS[i % len(COHORT_COLORS)]
            cohorts.append((label, color, curves))

    if not cohorts:
        print("nothing plottable")
        return

    fig, ax = plt.subplots(figsize=(11, 6))
    vis_lo, vis_hi = np.inf, -np.inf

    for label, color, curves in cohorts:
        present = [ch for ch in CHAN_ORDER if ch in curves]
        labelled = False   # one legend entry per cohort
        for chan in present:
            for subj_id, x_shared, mean, se in curves[chan]:
                xp = x_shared - x0
                has_se = se is not None
                if has_se:
                    band_ok = np.isfinite(mean) & np.isfinite(se)
                    if np.any(band_ok):
                        ax.fill_between(xp[band_ok],
                                        (mean - se)[band_ok],
                                        (mean + se)[band_ok],
                                        color=color, alpha=0.15, linewidth=0)
                ax.plot(xp, mean, color=color, lw=1.3, alpha=0.85,
                        label=(label if not labelled else None))
                labelled = True
                m = (xp >= 0.0) & (xp <= (x_hi - x_lo)) & np.isfinite(mean)
                if np.any(m):
                    se_m = (np.where(np.isfinite(se[m]), se[m], 0.0)
                            if has_se else 0.0)
                    vis_lo = min(vis_lo, float(np.min(mean[m] - se_m)))
                    vis_hi = max(vis_hi, float(np.max(mean[m] + se_m)))

    ax.axvline(-x0, color="k", ls="--", lw=1, alpha=0.6)   # ECG anchor
    ax.set_xlim(0.0, x_hi - x_lo)
    if np.isfinite(vis_lo) and np.isfinite(vis_hi) and vis_hi > vis_lo:
        pad = 0.05 * (vis_hi - vis_lo)
        ax.set_ylim(vis_lo - pad, vis_hi + pad)
    ax.set_xlabel("time (ms)  \u2014  0 = left edge of window")
    ax.set_ylabel("Feature / Global_Ref")
    ax.set_title("templates, first valid bin per subject, shared anchor timeline",
                 fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.legend(fontsize=9, loc="upper right", title="cohort")

    out_png = SCRIPT_DIR / f"{SCRIPT_STEM}.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"\nwrote {out_png}")


if __name__ == "__main__":
    main()
