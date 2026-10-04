#!/usr/bin/env python3
"""
so this doesn't actually test anything regarding the annealing it just prints out the lengths in seconds of all the bins because that way
you can compare it to the rules and the markings you made and figure out if the annealing is correct or not
"""
import sys, struct, os

def ask_path(prompt, exts={".bin"}):
    """Prompt until the answer names an existing file with an accepted
    extension. Surrounding quotes are stripped because Windows drag-and-drop
    and 'Copy as path' both add them."""
    while True:
        try:
            answer = input(prompt).strip()
        except (EOFError, KeyboardInterrupt):
            print("\naborted.")
            sys.exit(1)
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


def _skip_double_array(f):
    (n,) = struct.unpack("<Q", f.read(8))
    if n:
        f.seek(n * 8, os.SEEK_CUR)
    return n

def _read_header(f, path):
    """Parse and sanity-check the header BEFORE anything is printed, so a wrong
    file can be rejected without leaving half a table on screen."""
    head = f.read(8)
    if len(head) < 8:
        raise ValueError(f"only {len(head)} byte(s) -- too short to hold a bin count")
    (numBins,) = struct.unpack("<Q", head)
    rest = f.read(24 + 4)
    if len(rest) < 28:
        raise ValueError("header truncated before the channel count")
    ppgSR, ecgSR, scoringEpoch = struct.unpack("<ddd", rest[:24])
    (nChannels,) = struct.unpack("<I", rest[24:])
    if nChannels > 4096:
        raise ValueError(f"declares {nChannels} channels -- not an annealed .bin")
    if nChannels:
        f.read(nChannels * 4)   # native rates, skipped
    flags = f.read(3)
    if len(flags) < 3:
        raise ValueError("header truncated before the inversion flags")
    inv1, inv2, inv3 = struct.unpack("<BBB", flags)
    return numBins, ppgSR, ecgSR, scoringEpoch, nChannels, (inv1, inv2, inv3)

def list_bin_sizes(path):
    rows = []
    with open(path, "rb") as f:
        numBins, ppgSR, ecgSR, scoringEpoch, nChannels, inv = _read_header(f, path)
        inv1, inv2, inv3 = inv

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


def main():
    print("=== list_annealed_bin_sizes: per-bin sizes in an annealed .bin ===")
    while True:
        path = ask_path("  annealed .bin : ")
        try:
            list_bin_sizes(path)
            return 0
        except (struct.error, ValueError, OverflowError, MemoryError) as e:
            print(f"\n  could not parse {os.path.basename(path)} as an annealed .bin: {e}")
            print("  (an annealed file is named like *_annealed.bin)\n")

if __name__ == "__main__":
    sys.exit(main())
