"""
read_morphology_bin.py -- reader for the Section 4.6 morphology binaries
written by morphology_csv.hpp: <stem>_beats.bin and <stem>_templates.bin.

WHY THESE FILES EXIST IN BINARY ONLY. The beats file holds one waveform per
beat, so its size is (total beats) x (axis width). On one real record that was
33,384 beats x 3,160 samples = 105.5 million values. As CSV that is a ~1 GB text
file written one `ostream <<` at a time, which is why the writer looked like a
hang. Here the same values are raw little-endian doubles: ~840 MB, one pass, and
memory-mappable. The CSV keeps only the descriptor rows, which are what a person
reads.

Converting back to CSV, which the pipeline no longer writes:

    python read_morphology_bin.py --csv SUBJ_beats.bin

produces SUBJ_beats.csv byte-for-byte in the layout writeBeats() used to emit:
descriptor rows first, then one row per sample with one cell per beat. That file
is enormous by construction -- the reason it left the pipeline -- so generate it
only when something actually needs to read cells rather than slice arrays.

    from read_morphology_bin import read_beats, read_templates

    b = read_beats("SUBJ_beats.bin")
    ch1 = b["CH1"]
    ch1.frame            # pandas DataFrame, one row per beat, descriptors only
    ch1.waveforms        # (n_beats, width) float64, NaN where absent
    ch1.waveform(17)     # one beat's samples

    t = read_templates("SUBJ_templates.bin")
    t["CH1"].frame[["bin", "letter", "n_members", "beat_share"]]

pandas and numpy are optional. Without them you still get plain lists of dicts
and array.array rows -- the file format does not need either.

ON PADDING. The C++ side writes BeatRecord and TemplateRecord as raw structs, so
their compiler padding is part of the format. The struct format strings below
encode it explicitly ('2x' and '4x'), and morphology_csv.hpp has static_asserts
pinning sizeof() to 20 and 32 so a compiler that packs differently fails the
build rather than emitting a file this script would silently misparse.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, field
from typing import Any, BinaryIO

try:
    import numpy as _np
except ImportError:
    _np = None

try:
    import pandas as _pd
except ImportError:
    _pd = None


BEATS_MAGIC = b"DEXBEAT1"
TEMPLATES_MAGIC = b"DEXTMPL1"
SUPPORTED_VERSION = 1

# uint32 bin | 6x uint8 | 2 pad | int32 template_id | int32 r_col  == 20 bytes
_BEAT_FMT = "<I6B2xii"
_BEAT_SIZE = struct.calcsize(_BEAT_FMT)
assert _BEAT_SIZE == 20, _BEAT_SIZE

# uint32 bin | 7x uint8 | 1 pad | int32 x2 | uint32 n_members | double
# No padding before beat_share: n_members ends at offset 24, which is already
# 8-aligned, so the double follows immediately. Verified against sizeof()==32.
_TMPL_FMT = "<I7BxiiId"
_TMPL_SIZE = struct.calcsize(_TMPL_FMT)
assert _TMPL_SIZE == 32, _TMPL_SIZE


# --- code -> word, mirroring the enums the writer used ----------------------
# Kept as plain dicts rather than IntEnums so an unknown code round-trips as
# itself instead of raising: a file from a newer build should still be readable
# for every field this script does understand.
CATEGORY = {1: "pqrst", 2: "ectopic", 3: "noise"}
TUKEY_BEAT = {
    0: "not_eligible", 1: "kept", 2: "rej_rr_length",
    3: "rej_amplitude", 4: "rej_r_location", 5: "rej_wave_score",
}
TUKEY_AGG = {0: "not_eligible", 1: "kept", 2: "partial", 3: "removed"}
PREMATURE_BEAT = {0: "no", 1: "premature", 2: "vote"}
PREMATURE_AGG = {0: "no", 1: "mixed", 2: "vote", 3: "premature"}
CONFIRMED = {0: "presumed", 1: "confirmed", 2: "na"}
# annotation_types.hpp codes. 0 means unlabeled, which is NOT "unknown class":
# an unlabeled template is displayed as PQRST until an operator marks it.
LABEL = {
    0: "unlabeled", 1: "1) R Peak Noise", 2: "2) Minor Noise",
    3: "3) Blank.+Thresh.", 4: "4) PVC", 5: "5) PAC", 6: "6) Cond. Delay",
    7: "7) AF", 8: "8) SVT", 9: "9) VT", 10: "Benign Arr.",
    11: "Sig. Arr.", 12: "Other", 13: "Invert/Noninvert",
}

TEMPLATE_ID_UNSCORABLE = -2
TEMPLATE_ID_UNASSIGNED = -1


@dataclass
class Block:
    """One channel's records plus, for beats files, the sample matrix."""

    channel: str
    width: int
    records: list[dict[str, Any]] = field(default_factory=list)
    _samples: Any = None          # numpy (n, width) or list of array.array

    def __len__(self) -> int:
        return len(self.records)

    @property
    def waveforms(self):
        """(n_records, width) float64 with NaN where a sample is absent.

        NaN is meaningful, not missing data: a beat dropped by the baseline
        filter has descriptors -- including why it was dropped -- but no
        captured waveform, so its whole row is NaN. Compare against
        `template_id == -1` (unassigned) and `tukey` to tell the cases apart.
        """
        return self._samples

    def waveform(self, i: int):
        if self._samples is None:
            return None
        return self._samples[i]

    @property
    def frame(self):
        """pandas DataFrame of the descriptors, or a list of dicts without pandas."""
        if _pd is None:
            return self.records
        return _pd.DataFrame(self.records)


def _read_exact(f: BinaryIO, n: int) -> bytes:
    b = f.read(n)
    if len(b) != n:
        raise EOFError(f"expected {n} bytes, got {len(b)}")
    return b


def _read_header(f: BinaryIO, expect_magic: bytes) -> int:
    magic = _read_exact(f, 8)
    if magic != expect_magic:
        raise ValueError(
            f"not a morphology binary: magic {magic!r}, expected {expect_magic!r}"
        )
    (version, n_blocks) = struct.unpack("<II", _read_exact(f, 8))
    if version != SUPPORTED_VERSION:
        # A hard stop rather than a best-effort read. Guessing a record stride
        # is how a reader lands on the wrong offsets and produces plausible
        # numbers -- the failure this format's magic exists to prevent.
        raise ValueError(
            f"file version {version}; this reader handles {SUPPORTED_VERSION} only"
        )
    return n_blocks


def _read_block_header(f: BinaryIO) -> tuple[str, int, int]:
    (name_len,) = struct.unpack("<I", _read_exact(f, 4))
    channel = _read_exact(f, name_len).decode("ascii")
    (width, n_cols) = struct.unpack("<IQ", _read_exact(f, 12))
    return channel, width, n_cols


def _read_samples(f: BinaryIO, n_cols: int, width: int):
    """The bulk of the file: n_cols rows of `width` little-endian doubles."""
    if width == 0 or n_cols == 0:
        return None
    n = n_cols * width
    raw = _read_exact(f, n * 8)
    if _np is not None:
        return _np.frombuffer(raw, dtype="<f8", count=n).reshape(n_cols, width)
    import array
    a = array.array("d")
    a.frombytes(raw)
    return [a[i * width:(i + 1) * width] for i in range(n_cols)]


def read_beats(path: str) -> dict[str, Block]:
    """Read <stem>_beats.bin. Returns {channel: Block}, one record per beat."""
    out: dict[str, Block] = {}
    with open(path, "rb") as f:
        n_blocks = _read_header(f, BEATS_MAGIC)
        for _ in range(n_blocks):
            channel, width, n_cols = _read_block_header(f)
            blk = Block(channel=channel, width=width)

            # Records and samples are INTERLEAVED -- record, then that beat's
            # samples, then the next record. Read in that order; a reader that
            # slurps all records first lands mid-waveform.
            rows: list[dict[str, Any]] = []
            sample_rows = bytearray()
            for _ in range(n_cols):
                (bin_idx, cat, prem, tuk, conf, label, letter,
                 tid, r_col) = struct.unpack(_BEAT_FMT, _read_exact(f, _BEAT_SIZE))
                rows.append({
                    "bin": bin_idx,
                    "category": CATEGORY.get(cat, cat),
                    "premature": PREMATURE_BEAT.get(prem, prem),
                    "tukey": TUKEY_BEAT.get(tuk, tuk),
                    "confirmed": CONFIRMED.get(conf, conf),
                    "label_code": label,
                    "label": LABEL.get(label, f"code_{label}"),
                    "letter": chr(ord("A") + letter) if letter < 26 else letter,
                    "template_id": tid,
                    "r_col": r_col,
                })
                if width:
                    sample_rows += _read_exact(f, width * 8)

            blk.records = rows
            if width and sample_rows:
                if _np is not None:
                    blk._samples = _np.frombuffer(
                        bytes(sample_rows), dtype="<f8"
                    ).reshape(n_cols, width)
                else:
                    import array
                    a = array.array("d")
                    a.frombytes(bytes(sample_rows))
                    blk._samples = [a[i * width:(i + 1) * width]
                                    for i in range(n_cols)]
            out[channel] = blk
    return out


def read_templates(path: str) -> dict[str, Block]:
    """Read <stem>_templates.bin. One record per template column.

    Note which templates are ABSENT here by design: morphology_csv omits any
    template whose members were all Tukey-removed or all premature. Those beats
    are still in the beats file with their descriptors -- the templates file is
    the set of morphologies, not the set of everything the bank produced.
    """
    out: dict[str, Block] = {}
    with open(path, "rb") as f:
        n_blocks = _read_header(f, TEMPLATES_MAGIC)
        for _ in range(n_blocks):
            channel, width, n_cols = _read_block_header(f)
            blk = Block(channel=channel, width=width)

            rows: list[dict[str, Any]] = []
            sample_rows = bytearray()
            for _ in range(n_cols):
                (bin_idx, cat, prem, tuk, conf, label, letter, landmark,
                 tid, r_col, n_members, share) = struct.unpack(
                    _TMPL_FMT, _read_exact(f, _TMPL_SIZE))
                rows.append({
                    "bin": bin_idx,
                    "category": CATEGORY.get(cat, cat),
                    "premature": PREMATURE_AGG.get(prem, prem),
                    "tukey": TUKEY_AGG.get(tuk, tuk),
                    "confirmed": CONFIRMED.get(conf, conf),
                    "label_code": label,
                    "label": LABEL.get(label, f"code_{label}"),
                    "letter": chr(ord("A") + letter) if letter < 26 else letter,
                    "landmark_marked": bool(landmark),
                    "template_id": tid,
                    "r_col": r_col,
                    "n_members": n_members,
                    "beat_share": share,
                })
                if width:
                    sample_rows += _read_exact(f, width * 8)

            blk.records = rows
            if width and sample_rows:
                if _np is not None:
                    blk._samples = _np.frombuffer(
                        bytes(sample_rows), dtype="<f8"
                    ).reshape(n_cols, width)
                else:
                    import array
                    a = array.array("d")
                    a.frombytes(bytes(sample_rows))
                    blk._samples = [a[i * width:(i + 1) * width]
                                    for i in range(n_cols)]
            out[channel] = blk
    return out


def _fmt(v: float) -> str:
    """Match writeBeats()'s cell formatting: default ostream precision, and an
    EMPTY cell for NaN rather than the string "nan".

    That distinction mattered in the original and still does: a blank cell
    round-trips through every CSV reader as missing, which is what an absent
    sample is, while "nan" parses as a float in some readers and a string in
    others. Sample rows are dense enough that one surprise string in a numeric
    column poisons a whole column's dtype.
    """
    if v != v:            # NaN
        return ""
    return repr(float(v)) if v != int(v) else str(int(v))


def to_csv(path: str, out_path: str | None = None) -> str:
    """Convert a *_beats.bin (or *_templates.bin) to the CSV the pipeline used
    to write. Returns the path written.

    Streams row by row rather than building the table in memory: the sample rows
    are (beats x width) cells, and materialising 105 million strings to join
    them is how a converter runs out of memory on the file it exists to handle.
    """
    is_templates = path.endswith("_templates.bin")
    blocks = read_templates(path) if is_templates else read_beats(path)
    if out_path is None:
        out_path = path[:-len(".bin")] + ".csv"

    # Row order is writeBeats()/writeTemplates()' order exactly, so anything that
    # parsed the old files keeps working.
    beat_rows = [
        ("category", "category"), ("bin", "bin"), ("template", "_template"),
        ("premature", "premature"), ("tukey", "tukey"),
        ("confirmed", "confirmed"),
    ]
    tmpl_rows = beat_rows + [
        ("marking", "_marking"), ("n_members", "n_members"),
        ("beat_share", "beat_share"),
    ]
    rows = tmpl_rows if is_templates else beat_rows

    with open(out_path, "w", newline="") as f:
        for name, blk in blocks.items():
            if not blk.records:
                continue
            recs = blk.records
            # Derived columns the C++ writer composed rather than stored.
            for r in recs:
                tid = r["template_id"]
                if tid == TEMPLATE_ID_UNSCORABLE:
                    r["_template"] = "unscorable"
                elif tid == TEMPLATE_ID_UNASSIGNED:
                    r["_template"] = "unassigned"
                else:
                    cls = "PQRST" if r["label_code"] == 0 else r["label"]
                    r["_template"] = f"{cls}_{r['letter']}"
                if is_templates:
                    r["_marking"] = ("landmark" if r["landmark_marked"]
                                     else "class_only")

            f.write("channel," + ",".join([name] * len(recs)) + "\n")
            f.write("r_col," + ",".join(str(r["r_col"]) for r in recs) + "\n")
            for label, key in rows:
                # beat_share is the one non-integer descriptor, and the C++
                # writer produced it with std::to_string(double) -- six decimals,
                # always. Matched exactly so the converted file is
                # byte-identical to what writeTemplates() emits.
                if key == "beat_share":
                    cells = [f"{r[key]:.6f}" for r in recs]
                else:
                    cells = [str(r[key]) for r in recs]
                f.write(label + "," + ",".join(cells) + "\n")

            w = blk.waveforms
            if w is None:
                continue
            for s_idx in range(blk.width):
                cells = [_fmt(w[i][s_idx]) for i in range(len(recs))]
                f.write(str(s_idx) + "," + ",".join(cells) + "\n")
    return out_path


def _summarize(path: str) -> None:
    is_templates = path.endswith("_templates.bin")
    blocks = read_templates(path) if is_templates else read_beats(path)
    kind = "templates" if is_templates else "beats"
    print(f"{path}  ({kind})")
    for name, blk in blocks.items():
        print(f"  {name}: {len(blk)} records, width {blk.width}")
        if not blk.records:
            continue
        bins = {r["bin"] for r in blk.records}
        cats: dict[Any, int] = {}
        for r in blk.records:
            cats[r["category"]] = cats.get(r["category"], 0) + 1
        print(f"    bins {min(bins)}..{max(bins)}  categories {cats}")
        if is_templates:
            marked = sum(1 for r in blk.records if r["landmark_marked"])
            print(f"    landmark-marked {marked}/{len(blk)}")
        if blk.waveforms is not None and _np is not None:
            finite = int(_np.isfinite(blk.waveforms).any(axis=1).sum())
            print(f"    rows with any finite sample: {finite}/{len(blk)}")


if __name__ == "__main__":
    import sys
    args = [a for a in sys.argv[1:] if a != "--csv"]
    want_csv = "--csv" in sys.argv[1:]
    if not args:
        print(__doc__)
        print("usage: python read_morphology_bin.py [--csv] <file>_beats.bin "
              "| <file>_templates.bin ...")
        raise SystemExit(2)
    for p in args:
        if want_csv:
            out = to_csv(p)
            import os
            print(f"{p} -> {out}  ({os.path.getsize(out) / 1e6:.1f} MB)")
        else:
            _summarize(p)
