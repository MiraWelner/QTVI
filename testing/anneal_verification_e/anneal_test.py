#!/usr/bin/env python3
"""
list_annealed_bin_sizes.py -- list the size of every bin in an annealed .bin
(produced by annealOneFile, read by read_input_binfile in run_find_r_peaks.hpp).

Format (little-endian, matching run_find_r_peaks.hpp exactly):
  Header: u64 numBins, f64 ppgSR, f64 ecgSR, f64 scoringEpoch,
          u32 nChannels, nChannels x u32 nativeRate(skip),
          u8 inv1, u8 inv2, u8 inv3
  Per bin: u64 nPpgPairs, nPpgPairs x (u64,u64)
           u64 nEcgPairs, nEcgPairs x (u64,u64)
           5x [u64 count, count x f64]  (ppg_signal, ecg_signal_1/2/3, sleep_state)
           nChannels x [u64 nUp, nUp x f64 upsampled;
                        u64 nPairs, 2*nPairs x f64 raw t,v interleaved]

Usage:
  python list_annealed_bin_sizes.py <annealed.bin>
  python list_annealed_bin_sizes.py           # self-test on a synthetic file
"""
import sys, struct, os
import numpy as np


def _skip_double_array(f):
    (n,) = struct.unpack("<Q", f.read(8))
    if n:
        f.seek(n * 8, os.SEEK_CUR)
    return n

def list_bin_sizes(path):
    rows = []
    with open(path, "rb") as f:
        (numBins,) = struct.unpack("<Q", f.read(8))
        ppgSR, ecgSR, scoringEpoch = struct.unpack("<ddd", f.read(24))
        (nChannels,) = struct.unpack("<I", f.read(4))
        if nChannels:
            f.read(nChannels * 4)   # native rates, skipped
        inv1, inv2, inv3 = struct.unpack("<BBB", f.read(3))

        print(f"numBins={numBins}  ppgSR={ppgSR:g}  ecgSR={ecgSR:g}  "
              f"scoringEpoch={scoringEpoch:g}  nChannels={nChannels}  "
              f"inv=({inv1},{inv2},{inv3})")
        header = f"{'bin':>4}  {'minutes':>8}  {'ecg1_n':>8}  {'ppg_n':>8}  " \
                 f"{'ecg2_n':>8}  {'ecg3_n':>8}  {'sleep_n':>7}"
        print(header); print("-" * len(header))

        for i in range(numBins):
            (nPpg,) = struct.unpack("<Q", f.read(8))
            f.seek(nPpg * 16, os.SEEK_CUR)
            (nEcg,) = struct.unpack("<Q", f.read(8))
            f.seek(nEcg * 16, os.SEEK_CUR)

            n_ppg_sig = _skip_double_array(f)
            n_ecg1 = _skip_double_array(f)
            n_ecg2 = _skip_double_array(f)
            n_ecg3 = _skip_double_array(f)
            n_sleep = _skip_double_array(f)

            for _ in range(nChannels):
                _skip_double_array(f)                 # upsampled
                (nPairs,) = struct.unpack("<Q", f.read(8))
                f.seek(nPairs * 16, os.SEEK_CUR)       # raw (t,v) interleaved

            # Duration is driven by ecg_signal_1 at the file's own ecgSR --
            # the same channel/rate pairing create_ecg_ppg_pairs_raw detects
            # R-peaks on, so "bin length" here matches what the peak-finding
            # step actually sees per bin.
            minutes = (n_ecg1 / ecgSR / 60.0) if ecgSR > 0 else float("nan")
            print(f"{i:>4}  {minutes:>8.2f}  {n_ecg1:>8}  {n_ppg_sig:>8}  "
                  f"{n_ecg2:>8}  {n_ecg3:>8}  {n_sleep:>7}")
            rows.append({"bin": i, "minutes": minutes, "nPpgPairs": nPpg, "nEcgPairs": nEcg,
                         "ppg_signal": n_ppg_sig, "ecg_signal_1": n_ecg1,
                         "ecg_signal_2": n_ecg2, "ecg_signal_3": n_ecg3,
                         "sleep_state_signal": n_sleep})

        end = f.tell()
        total = os.path.getsize(path)
        print(f"\nparsed through byte {end} of {total} "
              f"({'OK, matches file size' if end == total else 'MISMATCH!'})")
    return rows


def _write_synthetic(path):
    rng = np.random.default_rng(0)
    numBins = 4
    nChannels = 0
    with open(path, "wb") as f:
        f.write(struct.pack("<Q", numBins))
        f.write(struct.pack("<ddd", 125.0, 250.0, 30.0))
        f.write(struct.pack("<I", nChannels))
        f.write(struct.pack("<BBB", 0, 0, 0))
        for i in range(numBins):
            f.write(struct.pack("<Q", 0))          # nPpgPairs
            f.write(struct.pack("<Q", 0))          # nEcgPairs
            n_ecg = 250 * (i + 1) * 60             # (i+1) minutes @ 250 Hz
            n_ppg = 125 * (i + 1) * 60
            ppg = rng.standard_normal(n_ppg)
            ecg1 = rng.standard_normal(n_ecg)
            ecg2 = rng.standard_normal(n_ecg // 2)  # deliberately different, to prove per-array sizes
            ecg3 = np.array([])                     # deliberately empty, to prove 0 handled
            sleep = rng.standard_normal(3)
            for sig in (ppg, ecg1, ecg2, ecg3, sleep):
                f.write(struct.pack("<Q", len(sig)))
                f.write(np.asarray(sig, "<f8").tobytes())
    return numBins, [125*(i+1)*60 for i in range(numBins)], \
                     [250*(i+1)*60 for i in range(numBins)], \
                     [250*(i+1)*60//2 for i in range(numBins)]


def selftest():
    import tempfile
    d = tempfile.mkdtemp()
    p = os.path.join(d, "synthetic_annealed.bin")
    numBins, exp_ppg, exp_ecg1, exp_ecg2 = _write_synthetic(p)
    print("(self-test: synthetic 4-bin annealed file, deliberately uneven sizes)\n")
    rows = list_bin_sizes(p)
    ok = True
    if len(rows) != numBins:
        ok = False
    for i, r in enumerate(rows):
        if r["ppg_signal"] != exp_ppg[i] or r["ecg_signal_1"] != exp_ecg1[i] \
           or r["ecg_signal_2"] != exp_ecg2[i] or r["ecg_signal_3"] != 0:
            ok = False
        if abs(r["minutes"] - (i + 1)) > 1e-9:   # writer built bin i as exactly (i+1) minutes
            ok = False
    print(f"\n[{'PASS' if ok else 'FAIL'}] recovered sizes AND minutes match what was written "
          f"(bin i should read back as exactly {{i+1}} minutes; including a deliberately-empty ecg_signal_3)")
    return 0 if ok else 1


def main(argv):
    if len(argv) < 2:
        return selftest()
    list_bin_sizes(argv[1])
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
