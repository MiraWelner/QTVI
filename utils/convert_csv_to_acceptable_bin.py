#!/usr/bin/env python3
"""
table_to_bin.py -- convert a CSV or XLSX of (time, ECG, PPG) samples into the
36-channel .bin that file_to_bin.cpp writes, so it can go straight into the
noise-marking GUI and the analysis pipeline.

INPUT
    A table with a time column in seconds and one column per signal, e.g.

        time_s,ecg,ppg
        0.000,0.012,0.503
        0.002,0.015,0.501
        ...

    Column names are configurable (--time-col, --ecg-col, --ppg-col). The
    native rate of each signal is the table's row rate, measured from the
    time column. A signal you don't name (or that isn't in the table) is
    written as a missing-channel placeholder, exactly as file_to_bin does.

OUTPUT  -- the same on-disk layout as file_to_bin.cpp / file_to_bin.hpp
    592-byte header (little-endian):
        uint32 version (=1), uint32 n_channels (=36), uint32 sleep_state_len,
        uint32 sizes_up[36], uint32 sizes_raw[36],
        float32 native_rates[36], float32 up_rates[36], uint32 sleep_size
    then, for each of the 36 channels in ChannelIdx order:
        float64 upsampled[sizes_up]                 -- at that channel's target rate
        float64 (t_epoch_ms, value)[sizes_raw]      -- the native samples
    then float64 sleep_stages[sleep_size].

    Channel 0 is the synthetic timestamp, built as make_binfile_edf_window
    builds it: absolute Unix-epoch ms on the ECG target grid (upsampled
    block), and (t, t) pairs on the ECG native grid (raw block, one extra
    sample: floor(duration * rate) + 1). Missing channels are a single -1.0 and
    a single (-1, -1) pair with both rates 0. No sleep staging is a single
    -1.0, the placeholder sleep_data_present() recognises.

THE START TIME MUST BE A REAL EPOCH. bin_chunk_loader.cpp treats the raw
x-values as epoch-ms only when channel 0's first sample is > 1e9; below that
it reads them as SECONDS, which would misplace every raw sample. So --start
defaults to a fixed real date rather than 0.

RESAMPLING is a line-for-line port of filterutils::upsample /
polyphase_resample (filter_utils.hpp) -- the same Blackman-windowed sinc
bank, per-phase normalisation, boundary handling and output length -- so a
signal converted here matches one converted by file_to_bin. Same limits:
integer rates, P and Q each <= 1000 after reducing by their GCD.

RATES
    Target rates come from --ecg-up / --ppg-up, or from your config.csv
    (--config config.csv --dataset MESA reads ecg_upsampled_rate,
    ppg_upsampled_rate and sleepstate_length from that dataset's row, the
    same columns config_loader.cpp reads). With neither, the ECG goes to
    1000 Hz and the PPG to 500 Hz.

A PPG WHOSE BASELINE IS EXACTLY 0 (synthetic data)
    The viewer scales a pulse by dividing by its foot value, (y - foot) / |foot|.
    Real PPG never sits on exactly 0 -- MESA's pleth is centred around 0 with
    its feet near -0.13 -- so this always works on real recordings. A
    synthetic PPG built on a flat 0.0 baseline has its foot at EXACTLY 0, the
    division is impossible, and every pulse panel is hidden.
    So when the PPG's minimum is at 0, the converter makes it like MESA's
    pleth: it subtracts the signal's mean, centring it on 0. The foot then
    sits below 0, as in MESA; pulse shape and timing are unchanged. It says so
    when it does this. --ppg-offset X instead adds a constant of your choosing
    (and turns the automatic centring off); --ppg-offset 0 keeps the PPG as is.

USAGE
    python table_to_bin.py synthetic_ecg_ppg_20min.xlsx --out bins/   # ECG 1000 Hz, PPG 500 Hz
    python table_to_bin.py rec.csv --config config.csv --dataset MESA --out bins/
    python table_to_bin.py rec.csv --ecg-up 1000 --ppg-up 500 --start 2026-01-01T22:00:00
    python table_to_bin.py synthetic.xlsx --out bins/          # centred automatically if needed
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import struct
import sys
from pathlib import Path

import numpy as np

# ---- layout constants (file_to_bin.hpp) --------------------------------------
NUM_CHANNELS = 36
BIN_HEADER_VERSION = 1
HEADER_SIZE = (4 + 4 * NUM_CHANNELS) * 4          # 592
CH_TIMESTAMP, CH_ECG1, CH_PPG = 0, 1, 4           # the only slots this tool fills

# Target rates when neither --ecg-up/--ppg-up nor --config gives one.
DEFAULT_ECG_UP = 1000.0
DEFAULT_PPG_UP = 500.0


# ---- filterutils port (filter_utils.hpp) -------------------------------------
def _gcd(a: int, b: int) -> int:
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def _polyphase_bank(P: int, Q: int, half_lobes: int) -> np.ndarray:
    """buildPolyphaseBank: Blackman-windowed sinc, split into P phases, each
    phase normalised to unit DC gain."""
    max_pq = max(P, Q)
    num_taps = 2 * half_lobes * max_pq + 1
    fc = 1.0 / max_pq
    M = num_taps - 1
    n = np.arange(num_taps, dtype=np.float64)
    x = n - M / 2.0
    with np.errstate(invalid="ignore", divide="ignore"):
        sinc = np.where(np.abs(x) < 1e-12, 1.0, np.sin(np.pi * fc * x) / (np.pi * x))
    w = 0.42 - 0.5 * np.cos(2.0 * np.pi * n / M) + 0.08 * np.cos(4.0 * np.pi * n / M)
    h = sinc * w
    sub_len = (num_taps + P - 1) // P
    bank = np.zeros((P, sub_len))
    for i in range(num_taps):                      # bank[i % P][i / P] = h[i]
        bank[i % P, i // P] = h[i]
    s = bank.sum(axis=1)
    ok = np.abs(s) > 1e-15
    bank[ok] /= s[ok, None]
    return bank


def polyphase_resample(x: np.ndarray, P: int, Q: int) -> np.ndarray:
    """polyphase_resample: out[m] = sum_k bank[(mQ)%P][k] * x[(mQ)/P - k + center],
    with samples outside the input contributing zero -- the boundary and
    interior paths of the C++ compute exactly this, the interior just skips
    the bounds checks."""
    if x.size == 0:
        return x.copy()
    if P == 1 and Q == 1:
        return x.copy()
    half_lobes = max(16, max(P, Q) // 2)
    bank = _polyphase_bank(P, Q, half_lobes)
    sub_len = bank.shape[1]
    in_len = x.size
    out_len = int(math.ceil(in_len * P / Q))
    center = half_lobes * max(P, Q) // P
    m = np.arange(out_len, dtype=np.int64)
    up = m * Q
    phase = up % P
    base = up // P
    # Zero-pad so every index base - k + center is in range; padding = zeros =
    # the C++ boundary path's "inIdx outside [0, inLen) contributes nothing".
    lo_pad = sub_len
    hi_pad = center + 1
    xp = np.concatenate([np.zeros(lo_pad), x.astype(np.float64), np.zeros(hi_pad)])
    out = np.zeros(out_len)
    for k in range(sub_len):
        out += bank[phase, k] * xp[base - k + center + lo_pad]
    return out


def upsample(x: np.ndarray, source_rate: float, target_rate: float) -> np.ndarray:
    """filterutils::upsample, including its integer-rate truncation and the
    P, Q <= 1000 refusal."""
    if x.size == 0:
        return x.copy()
    if source_rate == target_rate:
        return x.copy()
    g = _gcd(int(target_rate), int(source_rate))
    P, Q = int(target_rate) // g, int(source_rate) // g
    if P > 1000 or Q > 1000:
        raise ValueError(f"Resampling ratio {P}/{Q} too large - source rate: {source_rate}")
    return polyphase_resample(x, P, Q)


# ---- input -------------------------------------------------------------------
def read_table(path: Path, cols: list[str]) -> dict[str, np.ndarray]:
    """Read the named columns of a CSV or XLSX as float64 arrays."""
    suffix = path.suffix.lower()
    if suffix in (".xlsx", ".xlsm", ".xls"):
        import openpyxl
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        ws = wb.worksheets[0]
        rows = ws.iter_rows(values_only=True)
        header = [str(h).strip() if h is not None else "" for h in next(rows)]
        idx = {c: header.index(c) for c in cols if c in header}
        data = {c: [] for c in idx}
        for r in rows:
            if r is None or all(v is None for v in r):
                continue
            for c, i in idx.items():
                v = r[i] if i < len(r) else None
                data[c].append(np.nan if v is None or v == "" else float(v))
        wb.close()
    else:
        with open(path, newline="") as f:
            rd = csv.reader(f)
            header = [h.strip() for h in next(rd)]
            idx = {c: header.index(c) for c in cols if c in header}
            data = {c: [] for c in idx}
            for r in rd:
                if not r:
                    continue
                for c, i in idx.items():
                    v = r[i].strip() if i < len(r) else ""
                    data[c].append(np.nan if v == "" else float(v))
    return {c: np.asarray(v, dtype=np.float64) for c, v in data.items()}


def native_rate(t: np.ndarray) -> float:
    """The table's row rate from its time column. Must be regular: the .bin
    has no per-sample time for the upsampled block, so an irregular grid
    would be stretched without warning."""
    d = np.diff(t[np.isfinite(t)])
    if d.size == 0 or np.median(d) <= 0:
        raise ValueError("time column is empty, constant or decreasing")
    med = float(np.median(d))
    rate = 1.0 / med
    if abs(rate - round(rate)) < 1e-6 * rate:
        rate = float(round(rate))
    jitter = float(np.max(np.abs(d - med))) / med
    if jitter > 0.01:
        raise ValueError(f"time column is not uniformly sampled (worst step is "
                         f"{100 * jitter:.1f}% off the median {med * 1000:.4f} ms); "
                         f"resample it to a fixed rate first")
    return rate


def fill_nans(x: np.ndarray, name: str) -> np.ndarray:
    """EDF and .dat sources never carry NaN, so the C++ never sees one; the
    resampler would spread a NaN across a whole filter length. Linear
    interpolation, with a warning, keeps the file usable."""
    bad = ~np.isfinite(x)
    if not bad.any():
        return x
    if bad.all():
        raise ValueError(f"column '{name}' has no numeric values")
    i = np.arange(x.size)
    y = x.copy()
    y[bad] = np.interp(i[bad], i[~bad], x[~bad])
    print(f"  warning: {bad.sum()} blank/non-numeric value(s) in '{name}' "
          f"filled by linear interpolation", file=sys.stderr)
    return y


def rates_from_config(path: Path, dataset: str) -> dict[str, float]:
    """The ecg/ppg upsampled rates and sleepstate_length from one dataset row
    of config.csv -- the same columns, case-insensitive header, as
    config_loader.cpp."""
    with open(path, newline="") as f:
        rd = csv.reader(f)
        header = [h.strip().lower() for h in next(rd)]
        col = {h: i for i, h in enumerate(header)}
        for r in rd:
            r = [v.strip() for v in r]
            if "data_type" in col and r[col["data_type"]].upper() == dataset.upper():
                def num(name):
                    i = col.get(name)
                    try:
                        return float(r[i]) if i is not None and i < len(r) and r[i] else 0.0
                    except ValueError:
                        return 0.0
                return {"ecg_up": num("ecg_upsampled_rate"),
                        "ppg_up": num("ppg_upsampled_rate"),
                        "sleepstate_length": num("sleepstate_length")}
    raise ValueError(f"dataset '{dataset}' not found in {path}")


# ---- output ------------------------------------------------------------------
def write_bin(out_path: Path, ecg: np.ndarray | None, ppg: np.ndarray | None,
              rate: float, ecg_up: float, ppg_up: float,
              start_epoch_ms: float, sleep_state_len: int) -> None:
    sizes_up = [0] * NUM_CHANNELS
    sizes_raw = [0] * NUM_CHANNELS
    native_rates = [0.0] * NUM_CHANNELS
    up_rates = [0.0] * NUM_CHANNELS
    dt_raw_ms = 1000.0 / rate

    with open(out_path, "wb") as f:
        f.write(b"\0" * HEADER_SIZE)              # header written last, as in C++

        def write_missing(ch):
            f.write(struct.pack("<d", -1.0))
            f.write(struct.pack("<dd", -1.0, -1.0))
            sizes_up[ch], sizes_raw[ch] = 1, 1
            native_rates[ch], up_rates[ch] = 0.0, 0.0

        def write_signal(ch, x, up_rate):
            if x is None or x.size == 0 or up_rate <= 0.0:
                write_missing(ch)
                return
            up = upsample(x, rate, up_rate)
            f.write(up.astype("<f8").tobytes())
            t = start_epoch_ms + np.arange(x.size, dtype=np.float64) * dt_raw_ms
            f.write(np.column_stack([t, x]).astype("<f8").tobytes())
            sizes_up[ch], sizes_raw[ch] = up.size, x.size
            native_rates[ch], up_rates[ch] = rate, up_rate

        # Channel 0: synthetic timestamp on the ECG's grids (write_synthetic_timestamp).
        if ecg is not None and ecg.size and ecg_up > 0:
            dur = ecg.size / rate
            up_len = int(math.ceil(dur * ecg_up))
            f.write((start_epoch_ms + np.arange(up_len) * (1000.0 / ecg_up)).astype("<f8").tobytes())
            raw_len = int(math.floor(dur * rate)) + 1
            t = start_epoch_ms + np.arange(raw_len) * dt_raw_ms
            f.write(np.column_stack([t, t]).astype("<f8").tobytes())
            sizes_up[0], sizes_raw[0] = up_len, raw_len
            native_rates[0], up_rates[0] = rate, ecg_up
        else:
            write_missing(CH_TIMESTAMP)

        for ch in range(1, NUM_CHANNELS):
            if ch == CH_ECG1:
                write_signal(ch, ecg, ecg_up)
            elif ch == CH_PPG:
                write_signal(ch, ppg, ppg_up)
            else:
                write_missing(ch)

        # No sleep staging: the single -1.0 placeholder.
        f.write(struct.pack("<d", -1.0))
        sleep_size = 1

        f.seek(0)
        f.write(struct.pack("<III", BIN_HEADER_VERSION, NUM_CHANNELS, int(sleep_state_len)))
        f.write(struct.pack(f"<{NUM_CHANNELS}I", *sizes_up))
        f.write(struct.pack(f"<{NUM_CHANNELS}I", *sizes_raw))
        f.write(struct.pack(f"<{NUM_CHANNELS}f", *native_rates))
        f.write(struct.pack(f"<{NUM_CHANNELS}f", *up_rates))
        f.write(struct.pack("<I", sleep_size))


def verify_bin(path: Path) -> None:
    """Read the file back the way the consumers do (header, then channel
    blocks at offsets summed in channel order) and check it ends exactly."""
    b = path.read_bytes()
    ver, nch, epoch = struct.unpack_from("<III", b, 0)
    o = 12
    s_up = struct.unpack_from(f"<{NUM_CHANNELS}I", b, o); o += 4 * NUM_CHANNELS
    s_raw = struct.unpack_from(f"<{NUM_CHANNELS}I", b, o); o += 4 * NUM_CHANNELS
    o += 8 * NUM_CHANNELS                                   # both rate arrays
    (sleep,) = struct.unpack_from("<I", b, o); o += 4
    assert o == HEADER_SIZE and ver == BIN_HEADER_VERSION and nch == NUM_CHANNELS
    body = sum(8 * u + 16 * r for u, r in zip(s_up, s_raw))
    assert HEADER_SIZE + body + 8 * sleep == len(b), "size does not match the header"
    first_ts = struct.unpack_from("<d", b, HEADER_SIZE)[0]
    assert first_ts > 1e9, "channel 0 does not start with an epoch-ms value"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path, help="CSV or XLSX file")
    ap.add_argument("--out", type=Path, default=Path("."), help="output folder (default: current)")
    ap.add_argument("--time-col", default="time_s")
    ap.add_argument("--ecg-col", default="ecg")
    ap.add_argument("--ppg-col", default="ppg")
    ap.add_argument("--ecg-up", type=float, help=f"ECG target rate, Hz (default {DEFAULT_ECG_UP:g})")
    ap.add_argument("--ppg-up", type=float, help=f"PPG target rate, Hz (default {DEFAULT_PPG_UP:g})")
    ap.add_argument("--config", type=Path, help="config.csv to take target rates from")
    ap.add_argument("--dataset", help="dataset row in --config (MESA, BITTIUM, CHAOS, SHHS1, SHHS2)")
    ap.add_argument("--sleep-epoch", type=int, help="sleep epoch length, s (default: config, else 30)")
    ap.add_argument("--ppg-offset", type=float, default=None,
                    help="constant added to every PPG sample; overrides the automatic "
                         "centring of a PPG whose baseline is exactly 0 (0 = leave as is)")
    ap.add_argument("--start", default="2026-01-01T00:00:00+00:00",
                    help="recording start, ISO 8601 (default 2026-01-01T00:00:00Z)")
    a = ap.parse_args()

    cfg = {}
    if a.config:
        if not a.dataset:
            ap.error("--config needs --dataset")
        cfg = rates_from_config(a.config, a.dataset)

    cols = read_table(a.input, [a.time_col, a.ecg_col, a.ppg_col])
    if a.time_col not in cols:
        ap.error(f"no '{a.time_col}' column in {a.input}")
    rate = native_rate(cols[a.time_col])
    ecg = fill_nans(cols[a.ecg_col], a.ecg_col) if a.ecg_col in cols else None
    ppg = fill_nans(cols[a.ppg_col], a.ppg_col) if a.ppg_col in cols else None
    ppg_note = ""
    if ppg is not None:
        if a.ppg_offset is not None:
            if a.ppg_offset:
                ppg = ppg + a.ppg_offset
                ppg_note = f", PPG offset {a.ppg_offset:+g}"
        elif abs(float(np.min(ppg))) <= 1e-9 * max(float(np.max(np.abs(ppg))), 1.0):
            mean = float(np.mean(ppg))
            ppg = ppg - mean
            ppg_note = f", PPG centred on 0 (minus its mean {mean:g})"
            print(f"  note: the PPG's baseline is exactly 0, which hides every pulse "
                  f"panel (the viewer divides by the foot value). Centred it on 0 like "
                  f"MESA's pleth: foot now {float(np.min(ppg)):g}, peak {float(np.max(ppg)):g}.",
                  file=sys.stderr)
    if ecg is None:
        print("  warning: no ECG column -- the timestamp channel is also written "
              "as missing, as file_to_bin does", file=sys.stderr)

    ecg_up = a.ecg_up or cfg.get("ecg_up") or DEFAULT_ECG_UP
    ppg_up = a.ppg_up or cfg.get("ppg_up") or DEFAULT_PPG_UP
    sleep_len = a.sleep_epoch or int(cfg.get("sleepstate_length") or 30)

    start = dt.datetime.fromisoformat(a.start.replace("Z", "+00:00"))
    if start.tzinfo is None:
        start = start.replace(tzinfo=dt.timezone.utc)
    start_ms = start.timestamp() * 1000.0

    a.out.mkdir(parents=True, exist_ok=True)
    out_path = a.out / (a.input.stem + ".bin")      # <stem>.bin, as file_to_bin names it
    write_bin(out_path, ecg, ppg, rate, ecg_up, ppg_up, start_ms, sleep_len)
    verify_bin(out_path)

    n = len(cols[a.time_col])
    print(f"{a.input.name}: {n} rows at {rate:g} Hz ({n / rate / 60:.1f} min)")
    print(f"  ECG {'-> ' + format(ecg_up, 'g') + ' Hz' if ecg is not None else 'missing'}, "
          f"PPG {'-> ' + format(ppg_up, 'g') + ' Hz' if ppg is not None else 'missing'}, "
          f"sleep epoch {sleep_len} s, start {start.isoformat()}"
          + ppg_note)
    print(f"  wrote {out_path} ({out_path.stat().st_size:,} bytes), layout verified")
    return 0


if __name__ == "__main__":
    sys.exit(main())
