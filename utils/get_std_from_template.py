#!/usr/bin/env python3
"""
templates_sd.py -- pull every template's per-sample standard deviation out of
a <stem>_templates.bin (written by template_io::writeTemplatesBin).

    python templates_sd.py 3014720_20111011_templates.bin
    python templates_sd.py rec_templates.bin --out sd_out/ --ecg-rate 1000 --ppg-rate 500

WRITES (next to the input, or in --out):
    <stem>_template_sd_summary.csv   one row per template: channel, bin,
        template name, members, members averaged, then the SD summarised --
        mean, median and max over the template's samples, and the SD at the
        anchor column -- plus the template's own peak-to-peak amplitude, so an SD
        can be read relative to the waveform it belongs to.
    <stem>_template_sd.csv           the full curves, long format: one row per
        (template, sample) with sample index, samples and ms from the anchor,
        the template value and its SD.

THE ANCHOR is the column the writer records per bin (ChannelBlock::r_col):
the R peak for the ECG leads, but the SYSTOLIC PEAK for PPG. So on PPG rows
"from_anchor" is measured from the pulse's own peak, not from R.

WHAT THE SD IS: the per-sample sample standard deviation (ddof = 1) across the
beats the template was averaged over (members_clean). In the code the field is
called tmpl_iqr, and the viewer draws it as the "IQR band", but it is an SD --
bank_structs.hpp says so. Units are the template's own units (raw mV for ECG,
raw pulse units for PPG), not the viewer's normalised scale.

LAYOUT READ (template_io.hpp; checked, not assumed -- a mismatch stops with an
error rather than reading garbage):
    char[8] "MRPHTMPL", uint32 version (0), uint32 nBlocks
    per block: uint32 nameLen, char[nameLen] channel, uint32 width, uint64 nCols
      per column:
        TemplateRecord (40 bytes): uint32 bin, uint8 category, premature, tukey,
          confirmed, label_code, letter, landmark_marked, 1 pad, int32
          template_id, int32 r_col, uint32 n_members, uint32 n_excluded,
          uint8 too_few_beats, 3 pad, double beat_share
        double[width]                 waveform (NaN past the template's end)
        uint32 n, uint32[n]           members
        uint32 n, uint32[n]           members_clean
        TemplateTrailer (64 bytes)    mean_rr_ms is its last double
        uint32 n, double[n]           the per-sample SD
    then optionally the channel-lag trailer: [payload][uint32 bytes]["CHLAGRC1"]
"""
from __future__ import annotations

import argparse
import csv
import math
import struct
import sys
from pathlib import Path

MAGIC = b"MRPHTMPL"
VERSION = 0
LAG_MAGIC = b"CHLAGRC1"
REC = struct.Struct("<I7BxiiIIB3xd")      # TemplateRecord, 40 bytes
TRAILER_SIZE = 64
assert REC.size == 40

# label_code -> class, as the viewer and the markings CSV name templates.
CLASS = {0: "PQRST", 2: "NOISE", 4: "PVC", 5: "PAC", 9: "VT"}


class Reader:
    def __init__(self, data: bytes):
        self.b, self.p = data, 0

    def take(self, n: int) -> bytes:
        if self.p + n > len(self.b):
            raise ValueError(f"file ends early (need {n} bytes at offset {self.p}, "
                             f"file is {len(self.b)})")
        v = self.b[self.p:self.p + n]
        self.p += n
        return v

    def u32(self) -> int: return struct.unpack("<I", self.take(4))[0]
    def u64(self) -> int: return struct.unpack("<Q", self.take(8))[0]

    def doubles(self, n: int) -> list[float]:
        return list(struct.unpack(f"<{n}d", self.take(8 * n))) if n else []

    def vec_u32(self) -> list[int]:
        n = self.u32()
        return list(struct.unpack(f"<{n}I", self.take(4 * n))) if n else []

    def vec_d(self) -> list[float]:
        return self.doubles(self.u32())


def read_templates_bin(path: Path):
    r = Reader(path.read_bytes())
    if r.take(8) != MAGIC:
        raise ValueError("not a templates.bin (bad magic) -- is this the _beats.bin?")
    ver = r.u32()
    if ver != VERSION:
        raise ValueError(f"version {ver}; this script reads version {VERSION}")
    n_blocks = r.u32()
    out = []
    for _ in range(n_blocks):
        name = r.take(r.u32()).decode("ascii", "replace")
        width = r.u32()
        n_cols = r.u64()
        if width > 1_000_000 or n_cols > 1_000_000:
            raise ValueError(f"implausible block header for '{name}' "
                             f"(width {width}, columns {n_cols})")
        for _ in range(n_cols):
            (bin_, category, premature, tukey, confirmed, label, letter, lmk,
             tid, r_col, n_mem, n_exc, too_few, share) = REC.unpack(r.take(REC.size))
            wave = r.doubles(width)
            members = r.vec_u32()
            clean = r.vec_u32()
            trailer = r.take(TRAILER_SIZE)
            mean_rr = struct.unpack_from("<d", trailer, 56)[0]
            sd = r.vec_d()
            out.append(dict(channel=name, bin=bin_, template_id=tid,
                            name=f"{CLASS.get(label, f'CODE{label}')}_{chr(ord('A') + letter % 26)}",
                            r_col=r_col, n_members=n_mem,
                            n_averaged=len(clean) if clean else n_mem,
                            too_few_beats=too_few, mean_rr_ms=mean_rr,
                            wave=wave, sd=sd))
    rest = len(r.b) - r.p
    if rest and not (rest >= 12 and r.b.endswith(LAG_MAGIC)):
        raise ValueError(f"{rest} unexpected bytes after the last block")
    return out


def median(v: list[float]) -> float:
    v = sorted(v)
    m = len(v) // 2
    return v[m] if len(v) % 2 else 0.5 * (v[m - 1] + v[m])


def finite(v: list[float]) -> list[float]:
    return [x for x in v if x is not None and math.isfinite(x)]


def main() -> int:
    ap = argparse.ArgumentParser(description="Per-template SD from a templates.bin")
    ap.add_argument("input", type=Path)
    ap.add_argument("--out", type=Path, help="output folder (default: next to the input)")
    ap.add_argument("--ecg-rate", type=float, default=1000.0,
                    help="ECG template rate, Hz, for the ms column (default 1000)")
    ap.add_argument("--ppg-rate", type=float, default=500.0,
                    help="PPG template rate, Hz, for the ms column (default 500)")
    a = ap.parse_args()

    try:
        rows = read_templates_bin(a.input)
    except ValueError as e:
        print(f"error: {a.input}: {e}", file=sys.stderr)
        return 1

    out_dir = a.out or a.input.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = a.input.stem[:-len("_templates")] if a.input.stem.endswith("_templates") else a.input.stem
    summary_path = out_dir / f"{stem}_template_sd_summary.csv"
    curves_path = out_dir / f"{stem}_template_sd.csv"

    def rate_for(ch: str) -> float:
        return a.ppg_rate if ch.upper().startswith("PPG") else a.ecg_rate

    with open(summary_path, "w", newline="") as fs, open(curves_path, "w", newline="") as fc:
        ws, wc = csv.writer(fs), csv.writer(fc)
        ws.writerow(["channel", "bin", "template", "n_members", "n_averaged", "too_few_beats",
                     "mean_rr_ms", "sd_mean", "sd_median", "sd_max", "sd_at_anchor",
                     "anchor", "template_peak_to_peak"])
        wc.writerow(["channel", "bin", "template", "sample", "samples_from_anchor",
                     "ms_from_anchor", "template_value", "sd"])
        for t in rows:
            sd, wave = t["sd"], t["wave"]
            fsd = finite(sd)
            fw = finite(wave)
            sd_at_r = sd[t["r_col"]] if 0 <= t["r_col"] < len(sd) and math.isfinite(sd[t["r_col"]]) else ""
            ws.writerow([t["channel"], t["bin"], t["name"], t["n_members"], t["n_averaged"],
                         t["too_few_beats"], f"{t['mean_rr_ms']:.6g}",
                         f"{sum(fsd) / len(fsd):.6g}" if fsd else "",
                         f"{median(fsd):.6g}" if fsd else "",
                         f"{max(fsd):.6g}" if fsd else "",
                         f"{sd_at_r:.6g}" if sd_at_r != "" else "",
                         "systolic_peak" if t["channel"].upper().startswith("PPG") else "R",
                         f"{max(fw) - min(fw):.6g}" if fw else ""])
            ms_per = 1000.0 / rate_for(t["channel"])
            for s, v in enumerate(sd):
                w = wave[s] if s < len(wave) else float("nan")
                k = s - t["r_col"] if t["r_col"] >= 0 else ""
                wc.writerow([t["channel"], t["bin"], t["name"], s, k,
                             f"{k * ms_per:.6g}" if k != "" else "",
                             f"{w:.6g}" if math.isfinite(w) else "",
                             f"{v:.6g}" if math.isfinite(v) else ""])

    chans = sorted({t["channel"] for t in rows})
    print(f"{a.input.name}: {len(rows)} templates across {len(chans)} channel block(s): {', '.join(chans)}")
    print(f"  wrote {summary_path}")
    print(f"  wrote {curves_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
