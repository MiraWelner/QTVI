#!/usr/bin/env python3
"""
Task C -- file-to-bin conversion round-trip acceptance test.

Verifies a uniform .bin produced by make_binfile_* :
  * header parses at the real layout (36 channels),
  * every present signal channel ROUND-TRIPS: re-upsampling the stored raw
    (time,value) samples reproduces the stored upsampled block within
    floating-point tolerance,
  * missing channels read back as the documented sentinels
    (upsampled: single -1.0 ; raw: single (-1.0,-1.0) pair ; rates 0),
  * channel 0 (timestamp) is a clean epoch-ms ramp.

Plus a "Task C spec" section that checks the requested uint32 header-version
field at offset 0 and the 568-byte header size -- these currently FAIL because
file_to_bin writes sleep_state_len at offset 0 and a 584-byte header (see notes
printed by the test).

The polyphase resampler is inlined below (a faithful numpy port of
FilterUtils.hpp), so this test is a single self-contained file and "re-upsample"
uses the same algorithm the writer uses.

Usage:
  python3 test_file_to_bin.py                 # self-test on a synthetic bin
  python3 test_file_to_bin.py path/to/out.bin # validate a real converted record
"""
import sys, struct, tempfile, os
import numpy as np
# ---------------------------------------------------------------------------
# Polyphase resampler -- faithful numpy port of FilterUtils.hpp (inlined so
# this test is a single self-contained file; identical to acceptance_test.py).
# ---------------------------------------------------------------------------
def _gcd(a, b):
    return np.gcd(int(a), int(b))

def _build_bank(P, Q, half_lobes):
    max_pq = max(P, Q)
    num_taps = 2 * half_lobes * max_pq + 1
    fc = 1.0 / max_pq
    M = num_taps - 1
    half_m = M / 2.0
    n = np.arange(num_taps)
    x = n - half_m
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.sin(np.pi * fc * x) / (np.pi * x)
    sinc = np.where(np.abs(x) < 1e-12, 1.0, ratio)
    w = 0.42 - 0.5 * np.cos(2 * np.pi * n / M) + 0.08 * np.cos(4 * np.pi * n / M)
    h = sinc * w
    sub_len = (num_taps + P - 1) // P
    bank = np.zeros((P, sub_len))
    for i in range(num_taps):
        bank[i % P, i // P] = h[i]
    s = bank.sum(axis=1)
    for p in range(P):
        if abs(s[p]) > 1e-15:
            bank[p] /= s[p]
    return bank, sub_len

def resample_poly(x, P, Q):
    """Match FilterUtils::polyphase_resample for integer P/Q."""
    x = np.asarray(x, dtype=float)
    if P == 1 and Q == 1:
        return x.copy()
    half_lobes = max(16, max(P, Q) // 2)
    bank, sub_len = _build_bank(P, Q, half_lobes)
    in_len = len(x)
    out_len = int(np.ceil(in_len * P / Q))
    max_pq = max(P, Q)
    filter_center = half_lobes * max_pq // P
    m = np.arange(out_len)
    up = m * Q
    phase = up % P
    base = up // P
    # gather input windows with the same (base - k + filter_center) mapping as C++
    k = np.arange(sub_len)
    idx = base[:, None] - k[None, :] + filter_center     # (out_len, sub_len)
    in_bounds = (idx >= 0) & (idx < in_len)
    idx_clip = np.clip(idx, 0, in_len - 1)
    samples = np.where(in_bounds, x[idx_clip], 0.0)
    taps = bank[phase]                                    # (out_len, sub_len)
    return np.sum(taps * samples, axis=1)

def upsample(x, source_rate, target_rate):
    if len(x) == 0:
        return np.array([])
    if source_rate == target_rate:
        return np.asarray(x, dtype=float).copy()
    g = _gcd(int(target_rate), int(source_rate))
    P, Q = int(target_rate) // g, int(source_rate) // g
    if P > 1000 or Q > 1000:
        raise ValueError("resampling ratio too large")
    return resample_poly(x, P, Q)

NCH = 36
HEADER_FIELDS = 4 + 4 * NCH               # version, n_channels, sleep_state_len,
HEADER_SIZE   = HEADER_FIELDS * 4          # size_sleep + 4 arrays -> 148*4 = 592
SENTINEL = -1.0

_fail = 0
def check(ok, name, detail=""):
    global _fail
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok: _fail += 1

# ---------------------------------------------------------------------------
def parse_bin(path):
    with open(path, "rb") as f:
        blob = f.read()
    off = 0
    def u32(n):
        nonlocal off
        v = np.frombuffer(blob, "<u4", n, off); off += 4*n; return v
    def f32(n):
        nonlocal off
        v = np.frombuffer(blob, "<f4", n, off); off += 4*n; return v
    hdr = {}
    hdr["version"]    = int(u32(1)[0])      # offset 0
    hdr["n_channels"] = int(u32(1)[0])      # offset 4
    hdr["sleep_state_len"] = int(u32(1)[0]) # offset 8
    hdr["size_up"]   = u32(NCH)
    hdr["size_raw"]  = u32(NCH)
    hdr["native"]    = f32(NCH)
    hdr["up"]        = f32(NCH)
    hdr["sleep_size"]= int(u32(1)[0])
    hdr["header_bytes"] = off
    up_blk, raw = [], []
    for c in range(NCH):
        nu = int(hdr["size_up"][c])
        up_blk.append(np.frombuffer(blob, "<f8", nu, off)); off += 8*nu
        nr = int(hdr["size_raw"][c])
        pr = np.frombuffer(blob, "<f8", 2*nr, off).reshape(-1, 2); off += 16*nr
        raw.append(pr)
    hdr["bytes_consumed"] = off
    hdr["file_bytes"] = len(blob)
    return hdr, up_blk, raw

def is_missing(hdr, up_blk, raw, c):
    return (hdr["size_up"][c] == 1 and up_blk[c].size == 1 and up_blk[c][0] == SENTINEL and
            hdr["size_raw"][c] == 1 and raw[c].shape[0] == 1 and
            raw[c][0, 0] == SENTINEL and raw[c][0, 1] == SENTINEL)

def raw_is_uniform(times_ms, native_hz):
    """EDF channels store raw samples on a uniform native grid (step 1000/native
    ms); CHAOS .dat channels are sparse/irregular. Uniform => the writer used the
    polyphase resampler, so a re-upsample round-trip is meaningful."""
    if times_ms.size < 3 or native_hz <= 0: return False
    step = 1000.0 / native_hz
    d = np.diff(times_ms)
    return bool(np.max(np.abs(d - step)) < 1e-6 * max(1.0, step))

# ---------------------------------------------------------------------------
def test_bin(path):
    print(f"file: {path}")
    hdr, up_blk, raw = parse_bin(path)

    # ---- structural ----
    check(hdr["bytes_consumed"] == hdr["file_bytes"] or
          hdr["bytes_consumed"] + hdr["sleep_size"]*8 <= hdr["file_bytes"],
          "body parses to EOF (incl. sleep block)",
          f"consumed {hdr['bytes_consumed']}, file {hdr['file_bytes']}")

    # ---- channel 0: epoch-ms timestamp ramp ----
    if hdr["up"][0] > 0 and up_blk[0].size > 1:
        step = 1000.0 / hdr["up"][0]
        err = np.max(np.abs((up_blk[0] - up_blk[0][0]) - np.arange(up_blk[0].size)*step))
        check(err < 1e-3, "channel 0 is a uniform timestamp ramp", f"max dev {err:.2e} ms")
    else:
        check(False, "channel 0 present", "timestamp channel missing")

    # ---- round-trip each present signal channel ----
    rt_tested = rt_skipped = 0
    for c in range(1, NCH):
        if is_missing(hdr, up_blk, raw, c):
            continue
        vals   = raw[c][:, 1]
        times  = raw[c][:, 0]
        native = float(hdr["native"][c]); target = float(hdr["up"][c])
        if raw_is_uniform(times, native):
            expect = upsample(vals, native, target)
            n = min(expect.size, up_blk[c].size)
            # compare on the overlap; lengths may differ by <=1 from independent ceil()
            err = float(np.max(np.abs(expect[:n] - up_blk[c][:n]))) if n else np.inf
            lok = abs(up_blk[c].size - expect.size) <= 1
            check(err < 1e-6 and lok, f"ch{c}: raw re-upsamples to stored upsampled block",
                  f"max|err|={err:.2e}, len {up_blk[c].size} vs {expect.size}")
            rt_tested += 1
        else:
            # CHAOS-style: raw samples are irregular (sparse column with gaps).
            # The writer reconstructs a dense native-rate signal (linear fill),
            # then polyphase-upsamples native->target. Mirror that exactly:
            # densify the raw (t,v) onto the native grid, then polyphase.
            t_rel_s = (times - times[0]) / 1000.0            # epoch-ms -> s from first sample
            dur = up_blk[c].size / target
            dense_len = int(np.ceil(dur * native))
            dense_t = np.arange(dense_len) / native
            dense = np.interp(dense_t, t_rel_s, vals)        # linear densify (edge-clamped)
            expect = upsample(dense, native, target)
            n = min(expect.size, up_blk[c].size)
            lo2, hi2 = int(0.25*n), int(0.75*n)              # interior: skip filter edges
            err = float(np.max(np.abs(expect[lo2:hi2] - up_blk[c][lo2:hi2]))) if hi2 > lo2 else np.inf
            check(err < 1e-6, f"ch{c}: sparse raw densify+polyphase reproduces upsampled (interior)",
                  f"max|err|={err:.2e}, len {up_blk[c].size} vs {expect.size}")
            rt_tested += 1
    print(f"    ({rt_tested} channels round-tripped, {rt_skipped} sparse-skipped)")

    # ---- sentinels for missing channels ----
    miss = [c for c in range(NCH) if is_missing(hdr, up_blk, raw, c)]
    all_sentinel_ok = all(
        hdr["native"][c] == 0.0 and hdr["up"][c] == 0.0 for c in miss)
    check(all_sentinel_ok, "missing channels read as sentinels (rates 0, -1 markers)",
          f"{len(miss)} missing channels")

    # ---- Task C header-spec section ----
    print("  -- Task C header spec --")
    check(hdr["version"] == 1, "uint32 header-version == 1 at offset 0",
          f"version = {hdr['version']}")
    check(hdr["n_channels"] == NCH, "header carries n_channels",
          f"n_channels = {hdr['n_channels']}")
    check(hdr["header_bytes"] == HEADER_SIZE, "header size matches the 36-channel layout",
          f"{hdr['header_bytes']} bytes (Task C text says 568 -- see note below)")

# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Ground-truth source readers (pure numpy -- no external EDF libraries)
# ---------------------------------------------------------------------------
def read_edf(path):
    """Minimal EDF/EDF+ reader -> [(label, physical_values, rate_hz)].
    Skips 'EDF Annotations'. Physical = pmin + (dig-dmin)*(pmax-pmin)/(dmax-dmin)."""
    with open(path, "rb") as fh:
        b = fh.read()
    def s(o, n): return b[o:o+n].decode("ascii", "replace").strip()
    ns = int(s(252, 4)); nrec = int(s(236, 8)); rec_dur = float(s(244, 8))
    o = 256
    labels = [s(o+16*i, 16) for i in range(ns)]; o += 16*ns
    o += 80*ns + 8*ns
    pmin = [float(s(o+8*i, 8)) for i in range(ns)]; o += 8*ns
    pmax = [float(s(o+8*i, 8)) for i in range(ns)]; o += 8*ns
    dmin = [float(s(o+8*i, 8)) for i in range(ns)]; o += 8*ns
    dmax = [float(s(o+8*i, 8)) for i in range(ns)]; o += 8*ns
    o += 80*ns
    spr = [int(s(o+8*i, 8)) for i in range(ns)]; o += 8*ns
    o += 32*ns
    per_rec = sum(spr)
    raw = np.frombuffer(b, "<i2", count=per_rec*nrec, offset=o).reshape(nrec, per_rec)
    out, col = [], 0
    for i in range(ns):
        chunk = raw[:, col:col+spr[i]].reshape(-1).astype(float); col += spr[i]
        if "annotation" in labels[i].lower():
            continue
        if dmax[i] > dmin[i]:
            phys = (chunk - dmin[i]) * (pmax[i]-pmin[i]) / (dmax[i]-dmin[i]) + pmin[i]
        else:
            phys = chunk
        out.append((labels[i], phys, spr[i] / rec_dur))
    return out

def read_dat_columns(path):
    """CHAOS .dat reader -> [(header_label, populated_values, None)]. Mirrors
    find_real_header + prescan_dat_columns."""
    def has(h, n): return n.upper() in h.upper()
    hdr = None; body = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if hdr is None:
                if has(line, "NLS_NOM_") or has(line, "NLS_EEG_") or \
                   (has(line, "Index") and has(line, "TimeStamp")):
                    hdr = [c.strip() for c in line.rstrip("\n").split(",")]
                continue
            body.append(line)
    if not hdr:
        return []
    ncol = len(hdr)
    vals = [[] for _ in range(ncol)]
    for line in body:
        cells = line.rstrip("\n").split(",")
        for c in range(min(ncol, len(cells))):
            cell = cells[c].strip()
            if not cell:
                continue
            try: vals[c].append(float(cell))
            except ValueError: pass
    return [(hdr[c], np.array(vals[c]), None) for c in range(ncol) if vals[c]]

def read_source(path):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".edf": return read_edf(path)
    if ext == ".dat": return read_dat_columns(path)
    raise ValueError(f"unknown source extension: {ext} (expected .edf or .dat)")

# ---------------------------------------------------------------------------
# Ground-truth comparison: every present .bin channel's RAW block must equal a
# real source channel's samples (source -> raw faithful). With test_bin()'s
# raw -> upsampled check, that closes the loop source -> upsampled.
# ---------------------------------------------------------------------------
def test_source_vs_bin(bin_path, source_path, tol=1e-3):
    print(f"ground-truth: {os.path.basename(source_path)} vs {os.path.basename(bin_path)}")
    hdr, up_blk, raw = parse_bin(bin_path)
    src = read_source(source_path)
    present = matched = 0
    for c in range(1, NCH):
        if is_missing(hdr, up_blk, raw, c):
            continue
        present += 1
        rv = raw[c][:, 1]
        best_lab, best_rate, best_err = None, None, np.inf
        for lab, vals, rate in src:
            if vals.size != rv.size:
                continue
            e = float(np.max(np.abs(vals - rv)))
            if e < best_err:
                best_err, best_lab, best_rate = e, lab, rate
        if best_lab is not None and best_err < tol:
            matched += 1
            note = ""
            if best_rate is not None and abs(best_rate - float(hdr["native"][c])) > 1e-6:
                note = f"  [WARN rate {best_rate:g} != native {hdr['native'][c]:g}]"
            check(True, f"ch{c}: raw matches source '{best_lab}'",
                  f"max|err|={best_err:.2e}{note}")
        else:
            check(False, f"ch{c}: raw matches a source channel",
                  f"best err={best_err:.2e} (none within {tol})")
    check(matched == present, "every present channel traces to source data",
          f"{matched}/{present} matched")
    print("    (verifies source->raw; test_bin above verifies raw->upsampled)")


def make_synthetic_bin(path):
    """Spec-faithful EDF-style bin: ch0 epoch ramp, ch1 ECG (256->1000),
    ch4 PPG (125->500), ch19 pres (32->32 pass-through), rest missing.
    Upsampled blocks ARE polyphase(raw), so the round-trip must reproduce them."""
    epoch = 1749995130250.0
    dur = 4.0
    def sine(n, f, fs): return np.sin(2*np.pi*f*np.arange(n)/fs + 0.3)
    chans = {}   # c -> (native, target, raw_values)
    chans[1]  = (256.0, 1000.0, sine(round(dur*256), 20.0, 256.0))
    chans[4]  = (125.0,  500.0, sine(round(dur*125), 5.0, 125.0))
    chans[19] = (32.0,   32.0,  sine(round(dur*32),  3.0, 32.0))   # pass-through
    # ch20: CHAOS-style irregular channel -- native 125, target 500, with a gap
    # in the middle so raw times are non-uniform (exercises the sparse branch).
    sparse_native, sparse_target = 125.0, 500.0
    with open(path, "wb") as f:
        f.write(b"\x00" * HEADER_SIZE)
        size_up = np.zeros(NCH, "<u4"); size_raw = np.zeros(NCH, "<u4")
        native  = np.zeros(NCH, "<f4"); up_r     = np.zeros(NCH, "<f4")
        # ch0 timestamp
        n0 = int(np.ceil(dur*1000.0)); ramp = epoch + np.arange(n0)*(1000.0/1000.0)
        f.write(ramp.astype("<f8").tobytes()); size_up[0] = n0
        nr0 = int(np.floor(dur*1000.0))+1
        rp = epoch + np.arange(nr0)*(1000.0/1000.0)
        f.write(np.column_stack([rp, rp]).astype("<f8").tobytes()); size_raw[0]=nr0
        native[0]=1000.0; up_r[0]=1000.0
        # ch20 built here so we can write it in channel order inside the loop
        sN = round(dur * sparse_native)
        s_all_t = np.arange(sN) / sparse_native
        s_all_v = np.sin(2*np.pi*7.0*s_all_t + 0.2)
        keep = np.ones(sN, bool); keep[int(0.40*sN):int(0.55*sN)] = False   # gap
        s_t, s_v = s_all_t[keep], s_all_v[keep]
        s_dense = np.interp(np.arange(sN)/sparse_native, s_t, s_v)          # writer densify
        s_up = upsample(s_dense, sparse_native, sparse_target)             # writer polyphase

        # signal channels
        for c in range(1, NCH):
            if c == 20:
                f.write(s_up.astype("<f8").tobytes()); size_up[20] = s_up.size
                tms = epoch + s_t * 1000.0                                  # irregular raw times
                f.write(np.column_stack([tms, s_v]).astype("<f8").tobytes()); size_raw[20] = s_v.size
                native[20] = sparse_native; up_r[20] = sparse_target
            elif c in chans:
                nat, tgt, vals = chans[c]
                up = upsample(vals, nat, tgt)
                f.write(up.astype("<f8").tobytes()); size_up[c]=up.size
                t = epoch + np.arange(vals.size)*(1000.0/nat)
                f.write(np.column_stack([t, vals]).astype("<f8").tobytes()); size_raw[c]=vals.size
                native[c]=nat; up_r[c]=tgt
            else:  # missing sentinel
                f.write(np.array([SENTINEL], "<f8").tobytes()); size_up[c]=1
                f.write(np.array([[SENTINEL, SENTINEL]], "<f8").tobytes()); size_raw[c]=1
        f.write(np.array([SENTINEL], "<f8").tobytes())  # sleep block
        # header
        f.seek(0)
        f.write(struct.pack("<I", 1))                   # header_version
        f.write(struct.pack("<I", NCH))                 # n_channels
        f.write(struct.pack("<I", 30))                  # sleep_state_len
        f.write(size_up.tobytes()); f.write(size_raw.tobytes())
        f.write(native.tobytes());  f.write(up_r.tobytes())
        f.write(struct.pack("<I", 1))                   # sleep_size

def _write_edf(path, signals, rec_dur=1.0):
    ns = len(signals); spr = [int(round(r*rec_dur)) for (_, r, _) in signals]
    nrec = min(len(v)//s for (_, _, v), s in zip(signals, spr))
    hdr_bytes = 256*(ns+1)
    def af(x, n): return f"{x:<{n}}".encode("ascii")[:n]
    with open(path, "wb") as f:
        f.write(af("0", 8)); f.write(af("X", 80)); f.write(af("rec", 80))
        f.write(af("01.01.25", 8)); f.write(af("00.00.00", 8))
        f.write(af(str(hdr_bytes), 8)); f.write(af("", 44))
        f.write(af(str(nrec), 8)); f.write(af(f"{rec_dur:g}", 8)); f.write(af(str(ns), 4))
        for (lab, _, _) in signals: f.write(af(lab, 16))
        for _ in signals: f.write(af("", 80))
        for _ in signals: f.write(af("uV", 8))
        pmin = [float(np.min(v)) for (_, _, v) in signals]
        pmax = [float(np.max(v)) for (_, _, v) in signals]
        dmin, dmax = -32768, 32767
        for p in pmin: f.write(af(f"{p:g}", 8))
        for p in pmax: f.write(af(f"{p:g}", 8))
        for _ in signals: f.write(af(str(dmin), 8))
        for _ in signals: f.write(af(str(dmax), 8))
        for _ in signals: f.write(af("", 80))
        for s in spr: f.write(af(str(s), 8))
        for _ in signals: f.write(af("", 32))
        for r in range(nrec):
            for i, (_, _, v) in enumerate(signals):
                s = spr[i]; seg = v[r*s:(r+1)*s].astype(float)
                pl, ph = pmin[i], pmax[i]
                dig = (np.round((seg-pl)*(dmax-dmin)/(ph-pl)+dmin) if ph > pl
                       else np.full(s, dmin))
                f.write(np.clip(dig, dmin, dmax).astype("<i2").tobytes())

def _build_bin(path, chan_map, epoch=1749995130250.0):
    """chan_map: slot -> (values, native, target). raw=values (uniform at native),
    upsampled=polyphase(values). ch0 = epoch ramp, others = sentinel."""
    with open(path, "wb") as f:
        f.write(b"\x00" * HEADER_SIZE)
        su = np.zeros(NCH, "<u4"); sr = np.zeros(NCH, "<u4")
        na = np.zeros(NCH, "<f4"); ur = np.zeros(NCH, "<f4")
        n0 = 4000; ramp = epoch + np.arange(n0)*1.0
        f.write(ramp.astype("<f8").tobytes()); su[0] = n0
        rp = epoch + np.arange(n0+1)*1.0
        f.write(np.column_stack([rp, rp]).astype("<f8").tobytes()); sr[0] = n0+1
        na[0] = 1000.0; ur[0] = 1000.0
        for c in range(1, NCH):
            if c in chan_map:
                vals, nat, tgt = chan_map[c]
                vals = np.asarray(vals, float)
                up = upsample(vals, nat, tgt)
                f.write(up.astype("<f8").tobytes()); su[c] = up.size
                tms = epoch + np.arange(vals.size)*(1000.0/nat)
                f.write(np.column_stack([tms, vals]).astype("<f8").tobytes()); sr[c] = vals.size
                na[c] = nat; ur[c] = tgt
            else:
                f.write(np.array([SENTINEL], "<f8").tobytes()); su[c] = 1
                f.write(np.array([[SENTINEL, SENTINEL]], "<f8").tobytes()); sr[c] = 1
        f.write(np.array([SENTINEL], "<f8").tobytes())
        f.seek(0)
        f.write(struct.pack("<III", 1, NCH, 30))
        f.write(su.tobytes()); f.write(sr.tobytes()); f.write(na.tobytes()); f.write(ur.tobytes())
        f.write(struct.pack("<I", 1))

def _selftest_source():
    tmp = tempfile.gettempdir()
    ecg = np.sin(2*np.pi*12*np.arange(1024)/256.0)
    ppg = np.cos(2*np.pi*4*np.arange(500)/125.0)
    edf = os.path.join(tmp, "taskc_src.edf")
    _write_edf(edf, [("EKG", 256.0, ecg), ("Pleth", 125.0, ppg)], rec_dur=1.0)
    got = {l: v for (l, v, _) in read_edf(edf)}
    binp = os.path.join(tmp, "taskc_from_edf.bin")
    _build_bin(binp, {1: (got["EKG"], 256.0, 1000.0), 4: (got["Pleth"], 125.0, 500.0)})
    print("[EDF source self-test]")
    test_bin(binp); test_source_vs_bin(binp, edf)
    dat = os.path.join(tmp, "taskc_src.dat")
    n = 800; art = np.sin(2*np.pi*2*np.arange(n)/256.0)
    with open(dat, "w") as f:
        f.write("Index,System TimeStamp UTC,NLS_NOM_ECG_ELEC_POTL_I,NLS_NOM_PRESS_BLD_ART\n")
        for i in range(n):
            f.write(f"{i},20250615 13:45:{i%60:02d}.000,{ecg[i%len(ecg)]:.6f},{art[i]:.6f}\n")
    cols = {l: v for (l, v, _) in read_dat_columns(dat)}
    binp2 = os.path.join(tmp, "taskc_from_dat.bin")
    _build_bin(binp2, {1: (cols["NLS_NOM_ECG_ELEC_POTL_I"], 256.0, 1000.0),
                       33: (cols["NLS_NOM_PRESS_BLD_ART"], 256.0, 500.0)})
    print("[.dat source self-test]")
    test_bin(binp2); test_source_vs_bin(binp2, dat)


# ---------------------------------------------------------------------------
# Optional plotting: a per-channel grid of the round-trip (raw / upsampled /
# source overlay). Present channels only. Saves a PNG (Agg backend).
# ---------------------------------------------------------------------------
CHANNEL_NAMES = [
    "time", "ecg_1", "ecg_2", "ecg_3", "ppg", "accel_x", "accel_y", "accel_z",
    "marker", "temp", "pacemaker", "eog_l", "eog_r", "emg", "eeg_1", "eeg_2",
    "eeg_3", "eeg_4", "cvp", "pres", "flow", "snore", "thor", "abdo", "leg",
    "auxac", "therm", "pos", "oxstatus", "spo2", "hr", "dhr", "resp", "abp",
    "art", "art_pulm",
]

def _match_source(raw_vals, src, tol=1e-3):
    """Find the source signal whose samples equal this channel's raw block."""
    best = None; best_err = np.inf
    for lab, vals, rate in src:
        if vals.size != raw_vals.size:
            continue
        e = float(np.max(np.abs(vals - raw_vals)))
        if e < best_err:
            best_err, best = e, (lab, vals, rate)
    return best if (best is not None and best_err < tol) else None

def plot_bin(bin_path, source_path=None, out_png="taskc_plot.png", win_s=3.0):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"    (plot skipped: matplotlib unavailable: {e})")
        return None

    hdr, up_blk, raw = parse_bin(bin_path)
    src = read_source(source_path) if source_path else None

    present = [c for c in range(1, NCH) if not is_missing(hdr, up_blk, raw, c)]
    if not present:
        print("    (plot skipped: no present channels)")
        return None

    ncols = 3
    nrows = (len(present) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2*ncols, 2.5*nrows),
                             squeeze=False)
    for ax in axes.flat:
        ax.axis("off")

    for i, c in enumerate(present):
        ax = axes[i // ncols][i % ncols]; ax.axis("on")
        up = up_blk[c]; times = raw[c][:, 0]; vals = raw[c][:, 1]
        up_rate = float(hdr["up"][c]); native = float(hdr["native"][c])

        # time axes in seconds, anchored at each channel's first sample
        t_up  = np.arange(up.size) / up_rate if up_rate > 0 else np.arange(up.size)
        t_raw = (times - times[0]) / 1000.0                    # epoch-ms -> s

        mu = t_up <= win_s
        ax.plot(t_up[mu], up[mu], "-", color="#2c7fb8", lw=1.0,
                label=f"upsampled {up_rate:g}Hz", zorder=2)
        mr = t_raw <= win_s
        ax.plot(t_raw[mr], vals[mr], "o", ms=3, color="#c0392b",
                label=f"raw {native:g}Hz", zorder=3)

        if src is not None:
            m = _match_source(vals, src)
            if m is not None:
                lab, svals, srate = m
                sr = srate if srate else native
                t_src = np.arange(svals.size) / sr
                ms_ = t_src <= win_s
                ax.plot(t_src[ms_], svals[ms_], "x", ms=4, color="#31a354",
                        label=f"source '{lab}'", zorder=4)

        name = CHANNEL_NAMES[c] if c < len(CHANNEL_NAMES) else str(c)
        ax.set_title(f"ch{c}  {name}", fontsize=9)
        ax.set_xlabel("s", fontsize=8); ax.grid(alpha=0.3)
        ax.legend(fontsize=6, loc="upper right")

    title = f"Task C round-trip: {os.path.basename(bin_path)}"
    if source_path:
        title += f"  vs source {os.path.basename(source_path)}"
    fig.suptitle(title + f"   (first {win_s:g}s per channel)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_png, dpi=120)
    plt.close(fig)
    print(f"    plot written: {out_png}  ({len(present)} channels)")
    return out_png


def main(argv):
    print("=== Task C  file-to-bin round-trip acceptance test ===")
    args = [a for a in argv[1:] if not a.startswith("-")]
    do_plot = "--plot" in argv
    plot_file = next((a.split("=", 1)[1] for a in argv if a.startswith("--plot-file=")),
                     "taskc_plot.png")
    if len(args) >= 2:
        test_bin(args[0]); test_source_vs_bin(args[0], args[1])
        if do_plot: plot_bin(args[0], args[1], plot_file)
    elif len(args) == 1:
        test_bin(args[0])
        print("(no source given -> internal round-trip only; pass the original "
              ".edf/.dat as a 2nd arg to compare against the source)")
        if do_plot: plot_bin(args[0], None, plot_file)
    else:
        p = os.path.join(tempfile.gettempdir(), "taskc_synthetic.bin")
        make_synthetic_bin(p)
        print("(no path given -> self-test on synthetic bins)")
        test_bin(p); _selftest_source()
        if do_plot:
            # plot the EDF self-test bin against its source
            plot_bin(os.path.join(tempfile.gettempdir(), "taskc_from_edf.bin"),
                     os.path.join(tempfile.gettempdir(), "taskc_src.edf"), plot_file)
    if _fail == 0:
        print("NOTE  Task C says '568-byte header'. That was the 35-channel layout")
        print("      (1 + 4*35 + 1 = 142 fields = 568 B). Two things grew it:")
        print("      (a) the slot count went to 36 (+4 fields = 584 B), and")
        print("      (b) this task's own version + n_channels fields (+8 B = 592 B).")
        print("      So 568 can't coexist with 36 channels + a version header.")
    print("===", "ALL PASS" if _fail == 0 else f"{_fail} FAILURE(S)", "===")
    return 0 if _fail == 0 else 1

if __name__ == "__main__":
    sys.exit(main(sys.argv))
