"""
plot_pair_marks.py
==================
Plot ECG and PPG signals for a given bin, with vertical lines showing
where the PPG-ECG pair marks fall.

Directories are hard-coded at the top of this file. Pass the wave_markings
filename stem (without extension) and a bin number on the command line.
The matching annealed .bin is found automatically by subject ID.

Usage
-----
    python plot_pair_marks.py <wave_markings_stem> <bin_number> [options]

    wave_markings_stem   Filename stem, e.g.  3010724_20110811_1000_wave_markings
                         (with or without a .bin extension — both are fine)
    bin_number           0-based bin index to plot

Options
-------
    --ecg-rate FLOAT    ECG sample rate Hz  (default: from annealed header)
    --ppg-rate FLOAT    PPG sample rate Hz  (default: from annealed header)
    --save PATH         Save figure to this path instead of displaying it
    --no-rpeaks         Hide R-peak / PPG-valley scatter markers
    --unpaired          Also draw unpaired entries as grey dashed lines

Examples
--------
    python plot_pair_marks.py 3010724_20110811_1000_wave_markings 0
    python plot_pair_marks.py 3010724_20110811_1000_wave_markings 3 --save bin3.png
"""

import argparse
import glob
import os
import struct
import sys
from pathlib import Path

import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Hard-coded directory paths  ← edit these if folders move
# ---------------------------------------------------------------------------
WAVE_MARKINGS_DIR = (
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\output_mira\chaos\r_peak_finding_output"
)
ANNEALED_DIR = (
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\output_mira\chaos\annealed_output"
)


# ---------------------------------------------------------------------------
# File resolution
# ---------------------------------------------------------------------------


def resolve_wave_path(stem: str) -> Path:
    """
    Return the full path for a wave_markings file given a stem.
    Tries (in order):
      1. stem.bin / stem exactly as given
      2. stem + _wave_markings.bin  (if stem doesn't already end that way)
      3. glob for stem*_wave_markings.bin  (partial stem, e.g. subject+date only)
    """
    stem = stem.removesuffix(".bin")
    candidates = [stem + ".bin", stem]
    if not stem.endswith("_wave_markings"):
        candidates += [stem + "_wave_markings.bin", stem + "_wave_markings"]
    for candidate in candidates:
        p = Path(WAVE_MARKINGS_DIR) / candidate
        if p.exists():
            return p
    # Last resort: glob for anything starting with stem and ending _wave_markings.bin
    pattern = str(Path(WAVE_MARKINGS_DIR) / f"{stem}*_wave_markings.bin")
    matches = glob.glob(pattern)
    if len(matches) == 1:
        return Path(matches[0])
    if len(matches) > 1:
        print(f"[warn] Multiple wave_markings files match '{stem}*'; using first:")
        for m in matches:
            print(f"       {m}")
        return Path(matches[0])
    raise FileNotFoundError(
        f"Cannot find wave_markings file for stem '{stem}' in:\n  {WAVE_MARKINGS_DIR}"
    )


def resolve_annealed_path(stem: str) -> Path:
    """
    Find the annealed .bin whose filename starts with the same subject ID
    as the wave_markings stem.

    Wave markings stem:    3010724_20110811_1000_wave_markings
    Subject ID extracted:  3010724  (everything before the first '_')
    Annealed filename format: Z1011347_20210629T104500_20210629T184500_1000.bin
    """
    subject_id = stem.split("_")[0]
    # Files are named  <subjectID>_<datetime>_<epoch>.bin
    pattern = str(Path(ANNEALED_DIR) / f"{subject_id}_*.bin")
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(
            f"No original .bin found for subject ID '{subject_id}' in:\n  {ANNEALED_DIR}\n"
            f"  (searched: {pattern})"
        )
    if len(matches) > 1:
        # Narrow by epoch length (the numeric token just before '_wave_markings')
        epoch = None
        for part in reversed(stem.split("_")):
            if part.isdigit():
                epoch = part
                break
        if epoch:
            epoch_matches = [m for m in matches if f"_{epoch}_" in Path(m).name]
            if len(epoch_matches) == 1:
                return Path(epoch_matches[0])
        print(
            f"[warn] Multiple original .bin files found for subject '{subject_id}'; using first:"
        )
        for m in matches:
            print(f"       {m}")
    return Path(matches[0])


# ---------------------------------------------------------------------------
# Binary readers
# ---------------------------------------------------------------------------


def _read_idx_array(f):
    """uint64 count + count uint64 values; convert from 1-based to 0-based."""
    (sz,) = struct.unpack("<Q", f.read(8))
    if sz == 0:
        return np.array([], dtype=np.int64)
    raw = np.frombuffer(f.read(sz * 8), dtype=np.uint64).astype(np.int64)
    return raw - 1


def _read_signal(f):
    """uint64 count + count float64 samples."""
    (sz,) = struct.unpack("<Q", f.read(8))
    if sz == 0:
        return np.array([], dtype=np.float64)
    return np.frombuffer(f.read(sz * 8), dtype=np.float64).copy()


def read_wave_markings_bin(path):
    """
    Parse a wave_markings .bin (write_output_binfile layout).

    Per-bin keys returned:
        ch1_raw / ch1_squared / ch1_absval   — R-peak indices (0-based)
        ch2_* / ch3_*                        — same for channels 2 & 3
        ppg_max_amps / ppg_min_amps          — PPG event indices (0-based)
        ch*_sq_sig / ch*_abs_sig             — preprocessed signal arrays
        noise_flags                          — list[9] of bool
        pairs                                — Nx2 int64  [[ppg_idx, ecg_idx], ...]
                                               -1 means unpaired on that side
    """
    bins = []
    with open(path, "rb") as f:
        (num_bins,) = struct.unpack("<Q", f.read(8))
        for _ in range(num_bins):
            b = {}
            # 9 R-peak index arrays (3 channels x 3 methods)
            b["ch1_raw"] = _read_idx_array(f)
            b["ch1_squared"] = _read_idx_array(f)
            b["ch1_absval"] = _read_idx_array(f)
            b["ch2_raw"] = _read_idx_array(f)
            b["ch2_squared"] = _read_idx_array(f)
            b["ch2_absval"] = _read_idx_array(f)
            b["ch3_raw"] = _read_idx_array(f)
            b["ch3_squared"] = _read_idx_array(f)
            b["ch3_absval"] = _read_idx_array(f)
            # PPG event indices
            b["ppg_max_amps"] = _read_idx_array(f)
            b["ppg_min_amps"] = _read_idx_array(f)
            # 6 preprocessed signals
            b["ch1_sq_sig"] = _read_signal(f)
            b["ch1_abs_sig"] = _read_signal(f)
            b["ch2_sq_sig"] = _read_signal(f)
            b["ch2_abs_sig"] = _read_signal(f)
            b["ch3_sq_sig"] = _read_signal(f)
            b["ch3_abs_sig"] = _read_signal(f)
            # 9 noise flags
            b["noise_flags"] = list(struct.unpack("9B", f.read(9)))
            # Pairs: int64 interleaved [ppg, ecg, ppg, ecg, ...]  1-based; -1 = unpaired
            (num_pairs,) = struct.unpack("<Q", f.read(8))
            if num_pairs > 0:
                raw = np.frombuffer(f.read(num_pairs * 16), dtype=np.int64).copy()
                p = raw.reshape(num_pairs, 2)
                ppg_col = np.where(p[:, 0] == -1, -1, p[:, 0] - 1)
                ecg_col = np.where(p[:, 1] == -1, -1, p[:, 1] - 1)
                b["pairs"] = np.stack([ppg_col, ecg_col], axis=1)
            else:
                b["pairs"] = np.empty((0, 2), dtype=np.int64)
            bins.append(b)
    return bins


def read_annealed_bin(path, n_bins=None):
    """
    Parse an annealed-segments .bin produced by annealOneFile().

    Header (little-endian):
        uint64   numBins
        double   filePpgSR
        double   fileEcgSR
        double   scoringEpoch
        uint32   nChannels
        uint32   nativeSR[nChannels]   (skipped)

    Per bin:
        uint64   nPpgPairs
        (uint64,uint64) x nPpgPairs
        uint64   nEcgPairs
        (uint64,uint64) x nEcgPairs
        uint64+double[]  ppg_signal
        uint64+double[]  ecg_signal_1
        uint64+double[]  ecg_signal_2
        uint64+double[]  ecg_signal_3
        uint64+double[]  sleep_state
        nChannels x (uint64+double[] upsampled, uint64+double[2n] raw_tv)

    n_bins: if provided, the signal is split into exactly this many equal
            chunks (matched to wave_markings bin count).  If None, the
            file's own numBins is used directly.

    Returns list of dicts with keys:
        ppg_signal, ecg_signal_1/2/3  — float64 arrays at file rates
        ecg_sr, ppg_sr                — from header
        ppg_bin_indexs, ecg_bin_indexs — list of (uint64,uint64) pairs
    """
    with open(path, "rb") as f:
        num_bins_file = struct.unpack("<Q", f.read(8))[0]
        ppg_sr, ecg_sr, _epoch = struct.unpack("<ddd", f.read(24))
        n_ch = struct.unpack("<I", f.read(4))[0]
        if n_ch > 0:
            f.seek(n_ch * 4, 1)  # skip native rates

        print(
            f"  annealed header: numBins={num_bins_file}  "
            f"ppgSR={ppg_sr:.1f}  ecgSR={ecg_sr:.1f}  nChannels={n_ch}"
        )

        bins_raw = []
        for _ in range(num_bins_file):
            b = {"ppg_sr": ppg_sr, "ecg_sr": ecg_sr}

            n = struct.unpack("<Q", f.read(8))[0]
            b["ppg_bin_indexs"] = [struct.unpack("<QQ", f.read(16)) for _ in range(n)]

            n = struct.unpack("<Q", f.read(8))[0]
            b["ecg_bin_indexs"] = [struct.unpack("<QQ", f.read(16)) for _ in range(n)]

            b["ppg_signal"] = _read_signal(f)
            b["ecg_signal_1"] = _read_signal(f)
            b["ecg_signal_2"] = _read_signal(f)
            b["ecg_signal_3"] = _read_signal(f)
            b["sleep_state"] = _read_signal(f)

            for _ in range(n_ch):
                _read_signal(f)  # upsampled
                np_ = struct.unpack("<Q", f.read(8))[0]
                f.seek(np_ * 2 * 8, 1)  # raw_tv interleaved pairs

            bins_raw.append(b)

    if n_bins is None or n_bins == num_bins_file:
        return bins_raw

    # If the wave_markings has a different bin count, concatenate all signals
    # and re-slice evenly.
    print(
        f"  [warn] annealed numBins={num_bins_file} != wave numBins={n_bins}; re-slicing"
    )
    ecg1 = (
        np.concatenate([b["ecg_signal_1"] for b in bins_raw])
        if bins_raw
        else np.array([])
    )
    ecg2 = (
        np.concatenate([b["ecg_signal_2"] for b in bins_raw])
        if bins_raw
        else np.array([])
    )
    ecg3 = (
        np.concatenate([b["ecg_signal_3"] for b in bins_raw])
        if bins_raw
        else np.array([])
    )
    ppg = (
        np.concatenate([b["ppg_signal"] for b in bins_raw])
        if bins_raw
        else np.array([])
    )

    n_ecg = len(ecg1)
    n_ppg = len(ppg)
    ecg_bin = max(1, int(np.ceil(n_ecg / n_bins)))
    ppg_bin = max(1, int(np.ceil(n_ppg / n_bins)))

    def sl(arr, s, e):
        return (
            arr[s : min(e, len(arr))]
            if s < len(arr)
            else np.array([], dtype=np.float64)
        )

    return [
        {
            "ecg_signal_1": sl(ecg1, i * ecg_bin, (i + 1) * ecg_bin),
            "ecg_signal_2": sl(ecg2, i * ecg_bin, (i + 1) * ecg_bin),
            "ecg_signal_3": sl(ecg3, i * ecg_bin, (i + 1) * ecg_bin),
            "ppg_signal": sl(ppg, i * ppg_bin, (i + 1) * ppg_bin),
            "ecg_sr": ecg_sr,
            "ppg_sr": ppg_sr,
            "ppg_bin_indexs": [],
            "ecg_bin_indexs": [],
        }
        for i in range(n_bins)
    ]


def _time_axis(n, rate):
    return np.arange(n) / rate


def plot_bin(
    wave_bin,
    ann_bin,
    bin_idx,
    ecg_rate,
    ppg_rate,
    show_rpeaks=True,
    show_unpaired=False,
    save_path=None,
):
    ecg = ann_bin["ecg_signal_1"]
    ppg = ann_bin["ppg_signal"]
    pairs = wave_bin["pairs"]  # Nx2: [ppg_idx, ecg_idx]
    r_raw = wave_bin["ch1_raw"]
    ppg_v = wave_bin["ppg_min_amps"]

    # Trim to first 20 seconds for readability
    preview_samples_ecg = int(20 * ecg_rate)
    preview_samples_ppg = int(20 * ppg_rate)
    ecg = ecg[:preview_samples_ecg]
    ppg = ppg[:preview_samples_ppg]

    ecg_t = _time_axis(len(ecg), ecg_rate)
    ppg_t = _time_axis(len(ppg), ppg_rate)

    fig = plt.figure(figsize=(17, 7))
    fig.suptitle(
        f"ECG / PPG Pair Marks  —  Bin {bin_idx}", fontsize=14, fontweight="bold"
    )
    gs = gridspec.GridSpec(2, 1, hspace=0.45)
    ax_ecg = fig.add_subplot(gs[0])
    ax_ppg = fig.add_subplot(gs[1])

    cmap = matplotlib.colormaps.get_cmap("tab20")

    def normalise(sig):
        """Remove DC offset and scale to [-1, 1] for display."""
        if len(sig) == 0:
            return sig
        mn, mx = sig.min(), sig.max()
        rng = mx - mn
        if rng < 1e-12:
            return sig - mn
        return (sig - mn) / rng * 2 - 1

    ecg_disp = normalise(ecg)
    ppg_disp = normalise(ppg)

    if len(ecg):
        ax_ecg.plot(ecg_t, ecg_disp, color="#1a6eb5", lw=0.6, zorder=2)
    ax_ecg.set_ylabel("Normalised amplitude")
    ax_ecg.set_xlabel("Time (s)")
    ax_ecg.set_title("ECG (ch1) — coloured lines = paired R-peaks")
    ax_ecg.margins(x=0.005)

    if len(ppg):
        ax_ppg.plot(ppg_t, ppg_disp, color="#c0392b", lw=0.6, zorder=2)
    ax_ppg.set_ylabel("Normalised amplitude")
    ax_ppg.set_xlabel("Time (s)")
    ax_ppg.set_title("PPG — coloured lines = paired valleys  (same colour = same pair)")
    ax_ppg.margins(x=0.005)

    paired_count = 0
    for i, (ppg_idx, ecg_idx) in enumerate(pairs):
        fully_paired = (ppg_idx >= 0) and (ecg_idx >= 0)

        if fully_paired:
            colour = cmap(i % 20)
            alpha, ls, lw = 0.70, "-", 1.1
            paired_count += 1
        else:
            if not show_unpaired:
                continue
            colour = "#888888"
            alpha, ls, lw = 0.40, "--", 0.7

        if ecg_idx >= 0 and ecg_idx < len(ecg):
            ax_ecg.axvline(
                ecg_t[ecg_idx], color=colour, alpha=alpha, ls=ls, lw=lw, zorder=3
            )
        if ppg_idx >= 0 and ppg_idx < len(ppg):
            ax_ppg.axvline(
                ppg_t[ppg_idx], color=colour, alpha=alpha, ls=ls, lw=lw, zorder=3
            )

    if show_rpeaks:
        vr = r_raw[(r_raw >= 0) & (r_raw < len(ecg))]
        if len(vr):
            ax_ecg.scatter(
                ecg_t[vr],
                ecg_disp[vr],
                s=16,
                color="orange",
                zorder=5,
                label=f"R-peaks ({len(vr)})",
            )
            ax_ecg.legend(loc="upper right", fontsize=8)

        vv = ppg_v[(ppg_v >= 0) & (ppg_v < len(ppg))]
        if len(vv):
            ax_ppg.scatter(
                ppg_t[vv],
                ppg_disp[vv],
                s=16,
                color="#8e44ad",
                zorder=5,
                label=f"PPG valleys ({len(vv)})",
            )
            ax_ppg.legend(loc="upper right", fontsize=8)

    n_total = len(pairs)
    fig.text(
        0.5,
        0.005,
        f"{n_total} pairs total  |  {paired_count} fully paired  |  "
        f"{n_total - paired_count} unpaired",
        ha="center",
        fontsize=9,
        color="#555555",
    )

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved → {save_path}")
    else:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Plot ECG + PPG with pair-mark lines for one bin.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "wave_stem",
        help="Wave markings filename stem, e.g. 3010724_20110811_1000_wave_markings",
    )
    parser.add_argument("bin_number", type=int, help="0-based bin index to plot")
    parser.add_argument("--ecg-rate", type=float, default=None)
    parser.add_argument("--ppg-rate", type=float, default=None)
    parser.add_argument(
        "--save",
        default=None,
        help="Save figure here (PNG / PDF / SVG) instead of displaying",
    )
    parser.add_argument(
        "--no-rpeaks",
        action="store_true",
        help="Hide R-peak and PPG-valley scatter markers",
    )
    parser.add_argument(
        "--unpaired",
        action="store_true",
        help="Draw unpaired entries as grey dashed lines",
    )
    args = parser.parse_args()

    stem = Path(args.wave_stem).stem  # tolerate trailing .bin

    print(f"Resolving wave_markings file …")
    wave_path = resolve_wave_path(stem)
    print(f"  {wave_path}")

    print(f"Resolving annealed file …")
    ann_path = resolve_annealed_path(stem)
    print(f"  {ann_path}")

    print(f"\nReading wave_markings …")
    wave_bins = read_wave_markings_bin(wave_path)

    print(f"Reading annealed bin …")
    ann_bins = read_annealed_bin(ann_path, n_bins=len(wave_bins))

    if len(wave_bins) != len(ann_bins):
        sys.exit(
            f"ERROR: bin count mismatch — wave={len(wave_bins)}, annealed={len(ann_bins)}"
        )

    idx = args.bin_number
    if not (0 <= idx < len(wave_bins)):
        sys.exit(
            f"ERROR: bin {idx} out of range — file has {len(wave_bins)} bins (0-based)"
        )

    wb = wave_bins[idx]
    ab = ann_bins[idx]

    ecg_rate = args.ecg_rate or ab["ecg_sr"] or 256.0
    ppg_rate = args.ppg_rate or ab["ppg_sr"] or 64.0

    n_fully_paired = (
        int(np.sum((wb["pairs"][:, 0] >= 0) & (wb["pairs"][:, 1] >= 0)))
        if len(wb["pairs"])
        else 0
    )

    print(f"\nBin {idx}  ({len(wave_bins)} total in file)")
    print(f"  ECG samples : {len(ab['ecg_signal_1'])}  @ {ecg_rate:.1f} Hz")
    print(f"  PPG samples : {len(ab['ppg_signal'])}  @ {ppg_rate:.1f} Hz")
    print(f"  R-peaks     : {len(wb['ch1_raw'])}  (ch1 raw)")
    print(f"  PPG valleys : {len(wb['ppg_min_amps'])}")
    print(f"  Pairs       : {len(wb['pairs'])}  ({n_fully_paired} fully paired)\n")

    plot_bin(
        wb,
        ab,
        idx,
        ecg_rate=ecg_rate,
        ppg_rate=ppg_rate,
        show_rpeaks=not args.no_rpeaks,
        show_unpaired=args.unpaired,
        save_path=args.save,
    )


if __name__ == "__main__":
    main()
