"""
diagnose_offset.py

Standalone diagnostic for the "Mira signals are offset across rates" bug.

The .bin files at 500/1000/2000 Hz all derive from the same source ECG.
If the offset comes from the upsampler, the upsampled traces will be shifted
across rates, but the RAW (native-rate) (t, v) pairs should be identical
in all three files -- the raw block is just the source samples with their
real timestamps, untouched by the upsampler.

This script:
  1. Loads ch1_raw_pairs from all three Mira files for one bin.
  2. Reports their counts and the first few timestamps from each.
  3. Plots the raw pairs from all three files in one panel, plus a zoom.

If the raw points overlap perfectly across all three files: the upsampler
is introducing the offset. Patch resample.hpp.

If the raw points themselves are offset across files: the offset is being
introduced earlier in the pipeline (file_to_bin.cpp's timestamp logic,
the writer's raw block, etc.). The upsampler is innocent.

Usage:
    python diagnose_offset.py <subject_id> <bin_idx> [--zoom-start S] [--zoom-end S]
"""

import os
import struct
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Reuse paths from the comparison script.
BIN_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\4_wave_bound_files\mesa_rloc_mira"
)
OUT_DIR = Path(
    r"D:\USERS\MiraWelner\QTVI\testing\4_wave_finding_tests\single_bin"
)

MIRA_RATES = (500, 1000, 2000)
MIRA_COLORS = {500: "tab:blue", 1000: "black", 2000: "tab:cyan"}

PASSTHROUGH_NUM_CHANNELS = 40
RAW_ECG_CH1_SLOT = 1
MAX_SANE_IDX = 50_000_000
MAX_SANE_SIG = 200_000_000


def _read_only_target_bin_raw_pairs(path, target_idx):
    """Walk a wave_markings.bin to bin `target_idx` and return ONLY the
    ch1 raw (t, v) pairs for that bin. Skips everything else."""
    file_size = os.path.getsize(path)
    with open(path, "rb") as f:
        num_bins = struct.unpack("<Q", f.read(8))[0]
        if target_idx >= num_bins:
            print(f"  bin {target_idx} out of range (file has {num_bins})")
            return np.empty((0, 2), dtype=np.float64)

        for bin_i in range(target_idx + 1):
            keep = (bin_i == target_idx)

            # 11 idx arrays
            for _ in range(11):
                sz = struct.unpack("<Q", f.read(8))[0]
                f.seek(sz * 8, 1)

            # ppgSignal + 3 ECG signals = 4 signal arrays
            for _ in range(4):
                sz = struct.unpack("<Q", f.read(8))[0]
                f.seek(sz * 8, 1)

            # 6 preprocessed signals
            for _ in range(6):
                sz = struct.unpack("<Q", f.read(8))[0]
                f.seek(sz * 8, 1)

            # 9 noise flag bytes
            f.read(9)

            # pairs (numPairs prefix + 16 bytes each)
            n_pairs = struct.unpack("<Q", f.read(8))[0]
            f.seek(n_pairs * 16, 1)

            # 2 pair-vec skips
            for _ in range(2):
                sz = struct.unpack("<Q", f.read(8))[0]
                f.seek(sz * 16, 1)

            # 40 channel pass-through. Slot 1 (RAW_ECG_CH1_SLOT) is the
            # one we want raw pairs from -- ONLY for the target bin.
            target_pairs = np.empty((0, 2), dtype=np.float64)
            for ch in range(PASSTHROUGH_NUM_CHANNELS):
                # Upsampled block (always skipped)
                sz = struct.unpack("<Q", f.read(8))[0]
                if sz > MAX_SANE_SIG:
                    raise ValueError(f"signal size {sz} (file desync)")
                f.seek(sz * 8, 1)

                # Raw (t, v) pairs
                n = struct.unpack("<Q", f.read(8))[0]
                if n > MAX_SANE_IDX:
                    raise ValueError(f"raw pairs count {n} (file desync)")
                if keep and ch == RAW_ECG_CH1_SLOT:
                    if n > 0:
                        target_pairs = (
                            np.frombuffer(f.read(n * 16), dtype=np.float64)
                            .reshape(n, 2)
                            .copy()
                        )
                    # If n == 0 we leave target_pairs as the empty default.
                else:
                    f.seek(n * 16, 1)

            if keep:
                return target_pairs

    return np.empty((0, 2), dtype=np.float64)


def main():
    if len(sys.argv) < 3:
        print("Usage: python diagnose_offset.py <subject_id> <bin_idx> "
              "[--zoom-start S] [--zoom-end S]")
        sys.exit(1)

    subject_id = sys.argv[1]
    bin_idx = int(sys.argv[2])

    zoom_start = 20.0
    zoom_end = 22.0
    args = sys.argv[3:]
    while args:
        a = args.pop(0)
        if a == "--zoom-start":
            zoom_start = float(args.pop(0))
        elif a == "--zoom-end":
            zoom_end = float(args.pop(0))

    print(f"Subject:   {subject_id}")
    print(f"Bin idx:   {bin_idx}")
    print(f"Zoom:      {zoom_start}..{zoom_end} s")

    raw_by_rate = {}
    for rate in MIRA_RATES:
        path = BIN_DIR / f"{subject_id}_{rate}_wave_markings.bin"
        if not path.exists():
            print(f"  [{rate} Hz] missing: {path}")
            continue
        print(f"  [{rate} Hz] reading {path.name}...")
        pairs = _read_only_target_bin_raw_pairs(str(path), bin_idx)
        raw_by_rate[rate] = pairs
        if len(pairs) == 0:
            print(f"    no raw pairs in bin {bin_idx}")
            continue
        ts = pairs[:, 0]
        print(f"    {len(pairs)} pairs, "
              f"t range {ts.min():.6f}..{ts.max():.6f} s, "
              f"first 3 t's = {ts[:3]}")

    # Pairwise comparison: are the timestamps from different rates
    # aligned to within a sample?
    print()
    print("Pairwise raw-pair t-axis comparison:")
    rates_present = sorted(raw_by_rate)
    for i in range(len(rates_present)):
        for j in range(i + 1, len(rates_present)):
            ra, rb = rates_present[i], rates_present[j]
            pa, pb = raw_by_rate[ra], raw_by_rate[rb]
            if len(pa) == 0 or len(pb) == 0:
                continue
            n = min(len(pa), len(pb))
            t_diff = pa[:n, 0] - pb[:n, 0]
            v_diff = pa[:n, 1] - pb[:n, 1]
            print(f"  {ra}Hz vs {rb}Hz (first {n} pairs):")
            print(f"    t diff: max abs = {np.abs(t_diff).max():.6e} s, "
                  f"mean = {t_diff.mean():.6e} s")
            print(f"    v diff: max abs = {np.abs(v_diff).max():.6e}, "
                  f"identical: {np.array_equal(pa[:n, 1], pb[:n, 1])}")

    # Plot.
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, (ax_full, ax_zoom) = plt.subplots(2, 1, figsize=(14, 8))

    # Trim to bin window: raw timestamps are absolute (seconds from
    # recording start), so subtract the minimum t across rates to get
    # bin-local time.
    t0 = min(
        (p[:, 0].min() for p in raw_by_rate.values() if len(p) > 0),
        default=0.0,
    )
    for rate in MIRA_RATES:
        if rate not in raw_by_rate or len(raw_by_rate[rate]) == 0:
            continue
        pairs = raw_by_rate[rate]
        t_local = pairs[:, 0] - t0
        v = pairs[:, 1]
        marker_size = {500: 8, 1000: 5, 2000: 3}[rate]
        ax_full.scatter(
            t_local, v, s=marker_size, c=MIRA_COLORS[rate],
            label=f"{rate} Hz raw ({len(pairs)} pts)", alpha=0.7,
        )
        zoom_mask = (t_local >= zoom_start) & (t_local <= zoom_end)
        if zoom_mask.any():
            ax_zoom.scatter(
                t_local[zoom_mask], v[zoom_mask],
                s=marker_size + 4, c=MIRA_COLORS[rate],
                label=f"{rate} Hz raw ({zoom_mask.sum()} pts)", alpha=0.7,
            )

    ax_full.set_title(f"{subject_id} bin {bin_idx} — RAW (t, v) pairs from "
                      f"all three Mira files")
    ax_full.set_xlabel("Time (s, bin-local)")
    ax_full.set_ylabel("ECG")
    ax_full.legend()

    ax_zoom.set_title(f"Zoom: {zoom_start}..{zoom_end} s — "
                      f"if these dots don't all line up, the offset is "
                      f"in file_to_bin, not the upsampler")
    ax_zoom.set_xlabel("Time (s, bin-local)")
    ax_zoom.set_ylabel("ECG")
    ax_zoom.set_xlim(zoom_start, zoom_end)
    ax_zoom.legend()

    fig.tight_layout()
    out_path = OUT_DIR / f"{subject_id}_bin{bin_idx:03d}_raw_diagnostic.png"
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()