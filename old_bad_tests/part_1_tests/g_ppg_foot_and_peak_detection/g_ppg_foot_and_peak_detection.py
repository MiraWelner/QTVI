#!/usr/bin/env python3
"""
overlay_ppg_peaks.py -- overlay detected PPG peaks (ppgMaxAmps/ppgMinAmps in
the wave_markings .bin) on the PPG signal they were detected from (annealed
.bin), the same pairing read_output_binfile(wavePath, annealedPath) does in
peakfinding_io.hpp. Sibling of overlay_r_peaks.py -- same readers, PPG focus.

Run it with no arguments and answer the prompts. Both inputs are .bin files
produced downstream of file_to_bin -- the original .edf/.dat is NOT involved
(that pairing belongs to the file_to_bin conversion test, not this one):

  python overlay_ppg_peaks.py
    PPG peaks (*_peak_locations_all_beats.bin) : D:\\path\\to\\rec_peak_locations_all_beats.bin
    annealed  (*_annealed.bin)                 : D:\\path\\to\\rec_annealed.bin
    bin index [0-95]                           : 0
    kind [max/min/both]                        : both
    window seconds                             : 20

Formats (little-endian, matching peakfinding_io.hpp / run_find_r_peaks.hpp
exactly):

  wave_markings .bin (R-peaks-only):
    u64 numBins
    per bin: 9x [u64 count, count x u64 idx(1-based)]   (ch1/2/3 x raw/sq/abs)
             ppgMaxAmps, ppgMinAmps (same layout)
             6x [u64 count, count x f64]                (ch1/2/3 x sq/abs signal)
             9x u8 noise flags
             u64 numPairs, numPairs x 2 x i64 (interleaved ppg,ecg; -1 sentinel)

  annealed .bin:
    header: u64 numBins, f64 ppgSR, f64 ecgSR, f64 scoringEpoch,
            u32 nChannels, nChannels x u32 nativeRate(skip),
            u8 inv1, u8 inv2, u8 inv3
    per bin: u64 nPpgPairs, nPpgPairs x (u64,u64)
             u64 nEcgPairs, nEcgPairs x (u64,u64)
             5x [u64 count, count x f64]  (ppg_signal, ecg_signal_1/2/3, sleep_state)
             nChannels x [u64 nUp, nUp x f64 upsampled;
                          u64 nPairs, 2*nPairs x f64 raw t,v interleaved]

Writes the overlay to ppg_peak_overlay.png. Requires numpy and matplotlib.
"""
import sys, struct, os
import numpy as np

PLOT_NAME = "ppg_peak_overlay.png"   # written to the current working directory
KINDS = ("max", "min", "both")
MAX_BINS = 100_000   # a night at 5 min/bin is ~96; anything near this is a misread

# ---------------------------------------------------------------------------
# prompts
# ---------------------------------------------------------------------------
def _ask(prompt):
    try:
        return input(prompt).strip()
    except (EOFError, KeyboardInterrupt):
        print("\naborted.")
        sys.exit(1)

def ask_path(prompt, exts={".bin"}):
    """Prompt until the answer names an existing file with an accepted
    extension. Surrounding quotes are stripped because Windows drag-and-drop
    and 'Copy as path' both add them."""
    while True:
        answer = _ask(prompt)
        if not answer:
            print("        (required)")
            continue
        answer = answer.strip('"').strip("'").strip()
        path = os.path.expanduser(os.path.expandvars(answer))
        ext = os.path.splitext(path)[1].lower()
        if ext not in exts:
            print(f"        need {' or '.join(sorted(exts))}, got "
                  f"'{ext or '(no extension)'}'")
            continue
        if not os.path.isfile(path):
            print(f"        not found: {path}")
            continue
        return path

def ask_int(prompt, default, lo, hi):
    while True:
        answer = _ask(prompt)
        if not answer:
            return default
        try:
            v = int(answer)
        except ValueError:
            print(f"        not a whole number: '{answer}'")
            continue
        if not (lo <= v <= hi):
            print(f"        out of range, need {lo}-{hi}")
            continue
        return v

def ask_choice(prompt, options, default):
    while True:
        answer = _ask(prompt).lower()
        if not answer:
            return default
        if answer not in options:
            print(f"        need one of {'/'.join(options)}")
            continue
        return answer

def ask_float(prompt, default=None, positive=True):
    while True:
        answer = _ask(prompt)
        if not answer:
            return default
        try:
            v = float(answer)
        except ValueError:
            print(f"        not a number: '{answer}'")
            continue
        if positive and v <= 0:
            print("        must be greater than 0")
            continue
        return v

# ---------------------------------------------------------------------------
# wave_markings (.bin) reader
# ---------------------------------------------------------------------------
def _read_idx_array(f):
    (n,) = struct.unpack("<Q", f.read(8))
    if n == 0:
        return np.array([], dtype=np.int64)
    raw = np.frombuffer(f.read(n * 8), "<u8").astype(np.int64)
    return raw - 1   # writer added 1 for MATLAB compat; back to 0-based

def _read_sig_array(f):
    (n,) = struct.unpack("<Q", f.read(8))
    if n == 0:
        return np.array([], dtype=np.float64)
    return np.frombuffer(f.read(n * 8), "<f8").copy()

def read_wave_markings_bin(path):
    bins = []
    with open(path, "rb") as f:
        (numBins,) = struct.unpack("<Q", f.read(8))
        if numBins > MAX_BINS:
            raise ValueError(f"declares {numBins} bins -- not a wave_markings .bin")
        for _ in range(numBins):
            b = {}
            for ch in ("ch1", "ch2", "ch3"):
                b[ch] = {
                    "raw":     _read_idx_array(f),
                    "squared": _read_idx_array(f),
                    "absval":  _read_idx_array(f),
                }
            b["ppgMaxAmps"] = _read_idx_array(f)
            b["ppgMinAmps"] = _read_idx_array(f)
            for ch in ("ch1", "ch2", "ch3"):
                b[ch]["squared_signal"] = _read_sig_array(f)
                b[ch]["absval_signal"]  = _read_sig_array(f)
            flags = np.frombuffer(f.read(9), "u1")
            for i, ch in enumerate(("ch1", "ch2", "ch3")):
                b[ch]["raw_noisy"]     = bool(flags[i*3 + 0])
                b[ch]["squared_noisy"] = bool(flags[i*3 + 1])
                b[ch]["absval_noisy"]  = bool(flags[i*3 + 2])
            (numPairs,) = struct.unpack("<Q", f.read(8))
            if numPairs:
                tmp = np.frombuffer(f.read(numPairs * 16), "<i8").reshape(-1, 2).astype(np.float64)
                tmp[tmp == -1] = -1.0
                mask = tmp != -1.0
                tmp[mask] -= 1.0
                b["pairs"] = tmp
            else:
                b["pairs"] = np.empty((0, 2))
            bins.append(b)
    return bins

# ---------------------------------------------------------------------------
# annealed (.bin) reader
# ---------------------------------------------------------------------------
def _read_double_array(f):
    (n,) = struct.unpack("<Q", f.read(8))
    if n == 0:
        return np.array([], dtype=np.float64)
    return np.frombuffer(f.read(n * 8), "<f8").copy()

def read_annealed_bin(path, skip_passthrough=True):
    bins = []
    with open(path, "rb") as f:
        (numBins,) = struct.unpack("<Q", f.read(8))
        if numBins > MAX_BINS:
            raise ValueError(f"declares {numBins} bins -- not an annealed .bin")
        ppgSR, ecgSR, scoringEpoch = struct.unpack("<ddd", f.read(24))
        (nChannels,) = struct.unpack("<I", f.read(4))
        if nChannels:
            f.read(nChannels * 4)   # native rates, skipped
        inv1, inv2, inv3 = struct.unpack("<BBB", f.read(3))

        for _ in range(numBins):
            b = {}
            (nPpg,) = struct.unpack("<Q", f.read(8))
            b["ppg_bin_indexs"] = np.frombuffer(f.read(nPpg * 16), "<u8").reshape(-1, 2) if nPpg else np.empty((0, 2))
            (nEcg,) = struct.unpack("<Q", f.read(8))
            b["ecg_bin_indexs"] = np.frombuffer(f.read(nEcg * 16), "<u8").reshape(-1, 2) if nEcg else np.empty((0, 2))
            b["ppg_signal"] = _read_double_array(f)
            b["ecg_signal_1"] = _read_double_array(f)
            b["ecg_signal_2"] = _read_double_array(f)
            b["ecg_signal_3"] = _read_double_array(f)
            b["sleep_state_signal"] = _read_double_array(f)
            for ch in range(nChannels):
                (nUp,) = struct.unpack("<Q", f.read(8))
                up = f.read(nUp * 8)
                (nPairs,) = struct.unpack("<Q", f.read(8))
                rp = f.read(nPairs * 16)
                if not skip_passthrough:
                    b.setdefault("all_upsampled", []).append(np.frombuffer(up, "<f8").copy())
                    b.setdefault("all_raw_pairs_flat", []).append(np.frombuffer(rp, "<f8").copy())
            bins.append(b)
    return {"numBins": numBins, "ppgSR": ppgSR, "ecgSR": ecgSR,
            "scoringEpoch": scoringEpoch, "nChannels": nChannels,
            "ecg1_inverted": bool(inv1), "ecg2_inverted": bool(inv2),
            "ecg3_inverted": bool(inv3), "bins": bins}

# ---------------------------------------------------------------------------
# overlay
# ---------------------------------------------------------------------------
def overlay(waves, ann, bin_index=0, kind="both",
            rate_hz=None, window_sec=None, out_png=None):
    """waves / ann are the already-parsed structures from the two readers.
    kind: 'max' (systolic peaks), 'min' (valleys/feet), or 'both'."""
    if bin_index >= len(waves) or bin_index >= len(ann["bins"]):
        raise IndexError(f"bin {bin_index} out of range "
                          f"(wave has {len(waves)}, annealed has {len(ann['bins'])})")
    wb = waves[bin_index]
    ab = ann["bins"][bin_index]
    sig = ab["ppg_signal"]

    max_idx = wb["ppgMaxAmps"]; max_idx = max_idx[(max_idx >= 0) & (max_idx < len(sig))]
    min_idx = wb["ppgMinAmps"]; min_idx = min_idx[(min_idx >= 0) & (min_idx < len(sig))]

    fs = rate_hz if rate_hz else (ann["ppgSR"] if ann["ppgSR"] > 0 else 1.0)
    t = np.arange(len(sig)) / fs
    lo, hi = 0, len(sig)
    if window_sec:
        hi = min(len(sig), int(window_sec * fs))
        max_idx = max_idx[max_idx < hi]; min_idx = min_idx[min_idx < hi]

    print(f"bin {bin_index}: ppg_signal = {len(sig)} samples @ {fs:g} Hz "
          f"({len(sig)/fs:.1f} s), {len(max_idx)} max-amp peaks, {len(min_idx)} min-amp valleys")

    if out_png:
        try:
            import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        except Exception as e:
            print(f"  (plot skipped: {e})"); return sig, max_idx, min_idx, fs
        fig, ax = plt.subplots(figsize=(11, 3.5))
        ax.plot(t[lo:hi], sig[lo:hi], "-", color="#2c3e50", lw=0.8, label="ppg_signal")
        if kind in ("max", "both") and len(max_idx):
            ax.plot(t[max_idx], sig[max_idx], "o", color="#c0392b", ms=5,
                    label=f"max-amp peaks ({len(max_idx)})")
        if kind in ("min", "both") and len(min_idx):
            ax.plot(t[min_idx], sig[min_idx], "o", color="#2980b9", ms=5,
                    label=f"min-amp valleys ({len(min_idx)})")
        ax.set_xlabel("time (s)"); ax.set_ylabel("amplitude")
        ax.set_title(f"bin {bin_index}  PPG  {kind}")
        ax.legend(loc="upper right", fontsize=8)
        # Show exactly the range that was plotted. (This used to be a hardcoded
        # set_xlim(0, 20), which silently overrode the window setting.)
        if hi > lo:
            ax.set_xlim(t[lo], t[hi - 1])
        fig.tight_layout(); fig.savefig(out_png, dpi=120); plt.close(fig)
        print(f"  plot written: {out_png}")
    return sig, max_idx, min_idx, fs

# ---------------------------------------------------------------------------
def main():
    print("=== overlay_ppg_peaks: PPG peaks over the signal they came from ===")

    # Both inputs are .bin, so the extension cannot tell them apart and getting
    # the order wrong is easy. Parse inside the loop and re-ask on failure
    # rather than dying with a struct/overflow traceback.
    while True:
        wave_path = ask_path("  PPG peaks (*_peak_locations_all_beats.bin) : ")
        ann_path  = ask_path("  annealed  (*_annealed.bin)                 : ")
        try:
            waves = read_wave_markings_bin(wave_path)
            ann   = read_annealed_bin(ann_path)
            break
        except (struct.error, ValueError, OverflowError,
                MemoryError, IndexError) as e:
            print(f"\n  could not parse these as (wave_markings, annealed): {e}")
            print("  if you entered them in the other order, swap and retry.\n")

    nbins = min(len(waves), len(ann["bins"]))
    if nbins == 0:
        print("\nno bins in common between the two files -- nothing to overlay.")
        return 1
    if len(waves) != len(ann["bins"]):
        print(f"\n  NOTE: bin counts differ -- wave has {len(waves)}, "
              f"annealed has {len(ann['bins'])}; using the first {nbins}.")

    # Report what is actually populated so the next answers are informed.
    print(f"\n  {nbins} bin(s), ppgSR = {ann['ppgSR']:g} Hz, "
          f"ppg_signal in bin 0 = {ann['bins'][0]['ppg_signal'].size} samples")
    print(f"  bin 0 peak counts: max={waves[0]['ppgMaxAmps'].size}, "
          f"min={waves[0]['ppgMinAmps'].size}\n")

    bin_index = ask_int(f"  bin index [0-{nbins - 1}]            : ", 0, 0, nbins - 1)
    kind      = ask_choice(f"  kind [{'/'.join(KINDS)}]           : ", KINDS, "both")
    window    = ask_float("  window seconds (Enter = all)  : ", None)

    # The annealed header normally carries ppgSR; only ask when it does not.
    rate = None
    if ann["ppgSR"] <= 0:
        print("  (annealed header has no ppgSR)")
        rate = ask_float("  sample rate Hz                : ", None)
    print()

    if ann["bins"][bin_index]["ppg_signal"].size == 0:
        print(f"  NOTE: ppg_signal is empty in bin {bin_index} -- "
              f"the plot will have no trace.")

    overlay(waves, ann, bin_index=bin_index, kind=kind,
            rate_hz=rate, window_sec=window, out_png=PLOT_NAME)
    return 0

if __name__ == "__main__":
    sys.exit(main())
