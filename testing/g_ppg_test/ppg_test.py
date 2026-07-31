#!/usr/bin/env python3
"""
overlay_ppg_peaks.py -- overlay detected PPG peaks (ppgMaxAmps/ppgMinAmps in
the wave_markings .bin) on the PPG signal they were detected from (annealed
.bin), the same pairing read_output_binfile(wavePath, annealedPath) does in
peakfinding_io.hpp. Sibling of overlay_r_peaks.py -- same readers, PPG focus.

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

Usage:
  python overlay_r_peaks.py <wave_markings.bin> <annealed.bin>
                             [--bin N] [--channel 1|2|3] [--method raw|squared|absval]
                             [--plot out.png] [--window SECONDS] [--rate HZ]
  python overlay_r_peaks.py                      # self-test on synthetic files
"""
import sys, struct, os
import numpy as np

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
def overlay(wave_path, annealed_path, bin_index=0, kind="both",
            rate_hz=None, window_sec=None, out_png=None):
    """kind: 'max' (systolic peaks), 'min' (valleys/feet), or 'both'."""
    waves = read_wave_markings_bin(wave_path)
    ann = read_annealed_bin(annealed_path)
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
        fig.tight_layout()
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
        ax.set_xlim(0,20)
        fig.tight_layout(); fig.savefig(out_png, dpi=120); plt.close(fig)
        print(f"  plot written: {out_png}")
    return sig, max_idx, min_idx, fs

# ---------------------------------------------------------------------------
# self-test: build synthetic wave_markings.bin + annealed.bin, verify round-trip
# ---------------------------------------------------------------------------
def _write_annealed(path, ecg1, ppg, rate=250.0):
    with open(path, "wb") as f:
        f.write(struct.pack("<Q", 1))                       # numBins
        f.write(struct.pack("<ddd", 125.0, rate, 30.0))      # ppgSR, ecgSR, scoringEpoch
        f.write(struct.pack("<I", 0))                        # nChannels = 0 (no pass-through)
        f.write(struct.pack("<BBB", 0, 0, 0))                # inversion flags
        f.write(struct.pack("<Q", 0))                        # nPpgPairs
        f.write(struct.pack("<Q", 0))                        # nEcgPairs
        for sig in (ppg, ecg1, np.zeros(0), np.zeros(0), np.zeros(0)):  # ppg,ecg1,ecg2,ecg3,sleep
            f.write(struct.pack("<Q", len(sig)))
            f.write(np.asarray(sig, "<f8").tobytes())

def _write_wave_markings(path, r_peaks_ch1=(), ppg_max=(), ppg_min=()):
    with open(path, "wb") as f:
        f.write(struct.pack("<Q", 1))   # numBins
        def widx(v):
            v1 = (np.asarray(v, dtype=np.int64) + 1)  # 0-based -> 1-based
            f.write(struct.pack("<Q", len(v1)))
            f.write(v1.astype("<u8").tobytes())
        def wsig(v):
            f.write(struct.pack("<Q", len(v)))
            f.write(np.asarray(v, "<f8").tobytes())
        widx(r_peaks_ch1); widx([]); widx([])                # ch1 raw/sq/abs
        widx([]); widx([]); widx([])                          # ch2
        widx([]); widx([]); widx([])                          # ch3
        widx(ppg_max); widx(ppg_min)                          # ppgMaxAmps, ppgMinAmps
        wsig([]); wsig([])                                    # ch1 squared/absval signal
        wsig([]); wsig([])                                    # ch2
        wsig([]); wsig([])                                    # ch3
        f.write(bytes(9))                                     # noise flags
        f.write(struct.pack("<Q", 0))                         # numPairs

def selftest():
    import tempfile
    d = tempfile.mkdtemp()
    wave_p = os.path.join(d, "t_peak_locations_all_beats.bin")
    ann_p  = os.path.join(d, "t_annealed.bin")
    rate = 125.0
    n = 1250  # 10 s
    true_max = np.array([40, 150, 260, 370, 480, 590, 700, 810, 920, 1030])   # systolic peaks
    true_min = np.array([10, 120, 230, 340, 450, 560, 670, 780, 890, 1000])   # feet/valleys
    ppg = -0.3 * np.ones(n)
    ppg[true_max] = 1.0
    ppg[true_min] = -1.0
    ecg1 = np.sin(2*np.pi*1.0*np.arange(n)/rate)   # unused by this script, just needs to exist

    _write_annealed(ann_p, ecg1, ppg, rate=rate)
    _write_wave_markings(wave_p, ppg_max=true_max, ppg_min=true_min)

    plot_out = os.path.join(d, "selftest.png")
    sig, max_idx, min_idx, fs = overlay(wave_p, ann_p, bin_index=0, kind="both",
                                          out_png=plot_out)
    ok_rate = abs(fs - rate) < 1e-9
    ok_max = np.array_equal(np.sort(max_idx), np.sort(true_max))
    ok_min = np.array_equal(np.sort(min_idx), np.sort(true_min))
    ok_vals = np.allclose(sig[max_idx], 1.0) and np.allclose(sig[min_idx], -1.0)
    print(f"[{'PASS' if ok_rate else 'FAIL'}] recovered ppgSR == {rate}")
    print(f"[{'PASS' if ok_max else 'FAIL'}] recovered max-amp indices match ground truth")
    print(f"[{'PASS' if ok_min else 'FAIL'}] recovered min-amp indices match ground truth")
    print(f"[{'PASS' if ok_vals else 'FAIL'}] signal values at recovered peaks/valleys correct")
    return 0 if (ok_rate and ok_max and ok_min and ok_vals) else 1

def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    def opt(name, default=None, cast=str):
        for a in argv:
            if a.startswith(f"--{name}="):
                return cast(a.split("=", 1)[1])
        return default
    if len(args) < 2:
        print("(no wave/annealed paths given -> self-test)")
        return selftest()
    wave_path, annealed_path = args[0], args[1]
    overlay(wave_path, annealed_path,
            bin_index=opt("bin", 0, int),
            kind=opt("kind", "both", str),
            rate_hz=opt("rate", None, float),
            window_sec=opt("window", None, float),
            out_png=opt("plot", "ppg_peak_overlay.png", str))
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
