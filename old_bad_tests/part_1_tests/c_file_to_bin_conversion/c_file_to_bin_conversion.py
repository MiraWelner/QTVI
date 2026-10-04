#!/usr/bin/env python3
"""
round-trip check for the uniform 36-channel .bin against its original source (.dat or .edf)

Run it with no arguments and answer the prompts

Checks, per channel:
  * header parses (version=1, n_channels=36, 592-byte header),
  * MISSING channels read back as sentinels (size_up=1/-1, size_raw=1/(-1,-1), rates 0),
  * each PRESENT channel's stored RAW samples match a source column's values
    (source -> bin round-trip, within tolerance),
and optionally writes a per-channel plot: bin upsampled line + bin raw dots +
source overlay.
"""

import os
import struct
import sys

import numpy as np

NCH, HDR = 36, 592
CH_NAMES = ["timestamp","ecg_1","ecg_2","ecg_3","ppg","accel_x","accel_y","accel_z",
    "marker","temp","pacemaker","eog_l","eog_r","emg","eeg1","eeg2","eeg3","eeg4",
    "cvp","pres","flow","snore","thor","abdo","leg","auxac","therm","pos",
    "oxstatus","spo2","hr","dhr","resp","abp","art","art_pulm"]
SENT = -1.0
PLOT_NAME = "c_file_to_bin_conversion.png"   # written next to wherever you run this
_fail = 0
def check(ok, msg, detail=""):
    global _fail
    print(f"  [{'PASS' if ok else 'FAIL'}] {msg}" + (f"  --  {detail}" if detail else ""))
    if not ok: _fail += 1

# ---------------------------------------------------------------- prompts
def ask_path(prompt, exts):
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

# ---------------------------------------------------------------- bin reader
def parse_bin(path):
    with open(path, "rb") as f:
        head = f.read(HDR)
    ver, nch, sleep_len = struct.unpack_from("<III", head, 0)
    size_up  = np.frombuffer(head, "<u4", NCH, 12)
    size_raw = np.frombuffer(head, "<u4", NCH, 12 + 4*NCH)
    native   = np.frombuffer(head, "<f4", NCH, 12 + 8*NCH)
    up_rate  = np.frombuffer(head, "<f4", NCH, 12 + 12*NCH)
    size_sleep = struct.unpack_from("<I", head, 588)[0]
    # byte offset of each channel's data block
    offs, off = [], HDR
    for c in range(NCH):
        offs.append(off)
        off += int(size_up[c])*8 + int(size_raw[c])*16
    body_end = off + int(size_sleep)*8
    hdr = dict(ver=ver, nch=nch, sleep_len=sleep_len, size_up=size_up, size_raw=size_raw,
               native=native, up_rate=up_rate, size_sleep=size_sleep, offs=offs,
               body_end=body_end, file_size=os.path.getsize(path))
    return hdr

def read_channel(path, hdr, c, up_limit=None, raw_limit=None):
    nu = int(hdr["size_up"][c]); nr = int(hdr["size_raw"][c])
    if up_limit is not None: nu = min(nu, up_limit)
    if raw_limit is not None: nr = min(nr, raw_limit)
    with open(path, "rb") as f:
        f.seek(hdr["offs"][c])
        up = np.frombuffer(f.read(nu*8), "<f8")
        f.seek(hdr["offs"][c] + int(hdr["size_up"][c])*8)
        raw = np.frombuffer(f.read(nr*16), "<f8").reshape(-1, 2) if nr else np.empty((0,2))
    return up, raw

def is_missing(hdr, c):
    return int(hdr["size_up"][c]) == 1 and int(hdr["size_raw"][c]) == 1 and \
           float(hdr["native"][c]) == 0.0 and float(hdr["up_rate"][c]) == 0.0

# ---------------------------------------------------------------- source readers
def read_dat_columns(path):
    def has(h, n): return n.upper() in h.upper()
    hdr = None; body = []
    with open(path, "r", errors="replace") as f:
        for line in f:
            if hdr is None:
                if has(line,"NLS_NOM_") or has(line,"NLS_EEG_") or (has(line,"Index") and has(line,"TimeStamp")):
                    hdr = [c.strip() for c in line.rstrip("\n").split(",")]
                continue
            body.append(line)
    if not hdr: return []
    ncol = len(hdr); vals = [[] for _ in range(ncol)]
    for line in body:
        cells = line.rstrip("\n").split(",")
        for c in range(min(ncol, len(cells))):
            t = cells[c].strip()
            if not t: continue
            try: vals[c].append(float(t))
            except ValueError: pass
    return [(hdr[c], np.array(vals[c])) for c in range(ncol) if vals[c]]

def read_edf(path):
    with open(path, "rb") as fh: b = fh.read()
    def s(o,n): return b[o:o+n].decode("ascii","replace").strip()
    ns=int(s(252,4)); nrec=int(s(236,8)); dur=float(s(244,8)); o=256
    labels=[s(o+16*i,16) for i in range(ns)]; o+=16*ns; o+=80*ns+8*ns
    pmin=[float(s(o+8*i,8)) for i in range(ns)];o+=8*ns
    pmax=[float(s(o+8*i,8)) for i in range(ns)];o+=8*ns
    dmin=[float(s(o+8*i,8)) for i in range(ns)];o+=8*ns
    dmax=[float(s(o+8*i,8)) for i in range(ns)];o+=8*ns; o+=80*ns
    spr=[int(s(o+8*i,8)) for i in range(ns)];o+=8*ns; o+=32*ns
    per=sum(spr); raw=np.frombuffer(b,"<i2",count=per*nrec,offset=o).reshape(nrec,per)
    out=[]; col=0
    for i in range(ns):
        ch=raw[:,col:col+spr[i]].reshape(-1).astype(float); col+=spr[i]
        if "annotation" in labels[i].lower(): continue
        ph=(ch-dmin[i])*(pmax[i]-pmin[i])/(dmax[i]-dmin[i])+pmin[i] if dmax[i]>dmin[i] else ch
        out.append((labels[i], ph))
    return out

def read_source(path):
    ext=os.path.splitext(path)[1].lower()
    if ext==".dat": return read_dat_columns(path)
    if ext==".edf": return read_edf(path)
    raise ValueError(f"unknown source extension {ext}")

# ---------------------------------------------------------------- match & verify
def match_source(raw_vals, src, tol=1e-3):
    """Find source column whose populated values equal this channel's raw values."""
    best=None; besterr=np.inf
    for lab, vals in src:
        if vals.size != raw_vals.size: continue
        e=float(np.max(np.abs(vals-raw_vals))) if raw_vals.size else np.inf
        if e<besterr: besterr=e; best=(lab, vals, e)
    return best, besterr

def main():
    print("=== verify_bin: uniform 36-channel .bin round-trip check ===")
    binp = ask_path("  .bin file to verify         : ", exts={".bin"})
    srcp = ask_path("  original source (.dat/.edf) : ", exts={".dat", ".edf"})
    print()

    print(f"=== verify {os.path.basename(binp)} ===")
    hdr=parse_bin(binp)
    check(hdr["ver"]==1, "header version == 1", f"got {hdr['ver']}")
    check(hdr["nch"]==NCH, "n_channels == 36", f"got {hdr['nch']}")
    check(hdr["body_end"]==hdr["file_size"], "byte layout matches file size",
          f"walk end {hdr['body_end']} vs file {hdr['file_size']}")

    src = read_source(srcp)
    check(bool(src), "source yielded readable columns",
          f"{len(src)} column(s) from {os.path.basename(srcp)}")
    present=[c for c in range(1,NCH) if not is_missing(hdr,c)]
    missing=[c for c in range(1,NCH) if is_missing(hdr,c)]

    # sentinels
    sent_ok=True
    for c in missing:
        up,raw=read_channel(binp,hdr,c)
        if not (up.size==1 and up[0]==SENT and raw.shape[0]==1 and raw[0,0]==SENT and raw[0,1]==SENT):
            sent_ok=False
    check(sent_ok, "missing channels read as sentinels", f"{len(missing)} missing")

    # round-trip each present channel's raw block against the source
    if src:
        matched=0
        for c in present:
            _,raw=read_channel(binp,hdr,c)
            rv=raw[:,1]
            best,err=match_source(rv, src)
            if best and err<1e-3:
                matched+=1
                check(True, f"ch{c} {CH_NAMES[c]}: raw matches source '{best[0]}'", f"max|err|={err:.2e}")
            else:
                check(False, f"ch{c} {CH_NAMES[c]}: raw matches a source column",
                      f"best err={err:.2e}")
        check(matched==len(present), "every present channel traces to source data",
              f"{matched}/{len(present)}")

    make_plot(binp,hdr,present,src,PLOT_NAME)
    print("===", "ALL PASS" if _fail==0 else f"{_fail} FAILURE(S)", "===")
    return 0 if _fail==0 else 1

def make_plot(binp,hdr,present,src,out):
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception as e:
        print(f"    (plot skipped: {e})"); return
    win_s=6.0
    ncols=2; nrows=(len(present)+ncols-1)//ncols
    fig,axes=plt.subplots(nrows,ncols,figsize=(6.5*ncols,2.2*nrows),squeeze=False)
    for ax in axes.flat: ax.axis("off")
    for i,c in enumerate(present):
        ax=axes[i//ncols][i%ncols]; ax.axis("on")
        up_r=float(hdr["up_rate"][c]); nat=float(hdr["native"][c])
        nup=int(min(hdr["size_up"][c], int(win_s*up_r)+1))
        nraw=int(min(hdr["size_raw"][c], int(win_s*nat)+2))
        up,raw=read_channel(binp,hdr,c,up_limit=nup,raw_limit=nraw)
        t_up=np.arange(up.size)/up_r if up_r>0 else np.arange(up.size)
        ax.plot(t_up, up, "-", color="#2c7fb8", lw=0.8, label=f"bin upsampled {up_r:g}Hz")
        if raw.shape[0]:
            t0=raw[0,0]; t_raw=(raw[:,0]-t0)/1000.0 if raw[0,0]>1e9 else raw[:,0]
            m=t_raw<=win_s
            ax.plot(t_raw[m], raw[m,1], "o", ms=3, color="#c0392b", label=f"bin raw {nat:g}Hz")
            if src is not None:
                best,err=match_source(raw[:,1], src)
                if best and err<1e-3:
                    sv=best[1][:int(win_s*nat)+2]
                    ax.plot(np.arange(sv.size)/nat, sv, "x", ms=4, color="#31a354",
                            label=f"source '{best[0]}'")
        ax.set_title(f"ch{c} {CH_NAMES[c]}", fontsize=9); ax.set_xlabel("s",fontsize=8)
        ax.grid(alpha=.3); ax.legend(fontsize=6, loc="upper right")
    fig.suptitle(f"round-trip: {os.path.basename(binp)}  (first {win_s:g}s)", fontsize=11)
    fig.tight_layout(rect=[0,0,1,0.98]); fig.savefig(out, dpi=110); plt.close(fig)
    print(f"    plot written: {out}  ({len(present)} channels)")

if __name__=="__main__":
    sys.exit(main())
