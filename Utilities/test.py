#!/usr/bin/env python3
"""
edf_to_bin.py

Converts an EDF file to the .bin / .inf format expected by the
De-Mazumder / Prucka-style analysis tools.

Output layout
─────────────
  .bin  – 64-bit doubles written in Fortran (column-major) order for a
          [numChannels × numSamples] matrix.  That is: all channels for
          sample 0, then all channels for sample 1, … (matches MATLAB's
          fwrite of a [C×N] matrix).
  .inf  – Plain-text header that mirrors the format produced by
          edf_to_bin3.m (Sampling Rate on line 6, Number of Channels on
          line 7, Number of Samples on line 8, then lead labels).

Features
────────
  • Resamples every channel to a target rate (default 1000 Hz).
  • Slices long recordings into 24-hour chunks.
  • Accepts an optional list of target signal names (default: ECG_1,
    ECG_2, ECG_3).  If none are found it falls back to every signal
    whose label contains "ECG" or "EKG", and finally to all signals.

Requirements
────────────
  pip install pyedflib numpy scipy
"""

import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pyedflib

# ──────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────

INPUT_DIR = r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\0_original_files\mesa_files"
OUTPUT_DIR = r"D:\USERS\MiraWelner\QTVI\QTVI-data-files\1_bin_mat_files\mesa_bin_deep"


# ──────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────


def pick_channels(signal_labels: list[str], targets: list[str] | None) -> list[int]:
    """Return indices of the channels to extract.

    Priority:
      1. Exact matches against *targets* (case-insensitive).
      2. Any label containing 'ECG' or 'EKG'.
      3. All channels.
    """
    labels_upper = [l.strip().upper() for l in signal_labels]

    if targets:
        targets_upper = [t.strip().upper() for t in targets]
        indices = [i for i, l in enumerate(labels_upper) if l in targets_upper]
        if indices:
            return indices

    # Fallback: anything ECG / EKG-like
    indices = [i for i, l in enumerate(labels_upper) if "ECG" in l or "EKG" in l]
    if indices:
        return indices

    # Last resort: everything
    return list(range(len(signal_labels)))


def resample_channel(data: np.ndarray, fs_orig: float, fs_new: float) -> np.ndarray:
    """Resample a 1-D signal from fs_orig to fs_new."""
    if fs_orig == fs_new:
        return data
    new_len = int(np.ceil(len(data) * fs_new / fs_orig))
    old_t = np.arange(len(data)) / fs_orig
    new_t = np.arange(new_len) / fs_new
    return np.interp(new_t, old_t, data)


# ──────────────────────────────────────────────
# Core conversion
# ──────────────────────────────────────────────


def edf_to_bin(
    edf_path: str,
    output_dir: str,
    target_fs: int = 1000,
    chunk_hours: int = 24,
    target_signals: list[str] | None = None,
    lead_labels: list[str] | None = None,
) -> None:
    """Read *edf_path* and write .bin / .inf chunk(s) into *output_dir*."""

    edf_path = Path(edf_path)
    if not edf_path.is_file():
        sys.exit(f"Error: EDF file not found: {edf_path}")

    reader = pyedflib.EdfReader(str(edf_path))
    try:
        n_signals = reader.signals_in_file
        signal_labels = [reader.getLabel(i) for i in range(n_signals)]
        channel_indices = pick_channels(signal_labels, target_signals)

        if not channel_indices:
            sys.exit("Error: no matching channels found in the EDF file.")

        n_ch = len(channel_indices)
        sample_rates = [reader.getSampleFrequency(i) for i in channel_indices]
        fs_orig = sample_rates[0]

        # Read & (optionally) resample each selected channel
        print(f"Reading {n_ch} channel(s) from {edf_path.name} …")
        channels: list[np.ndarray] = []
        for idx in channel_indices:
            sig = reader.readSignal(idx)
            fs_ch = reader.getSampleFrequency(idx)
            channels.append(resample_channel(sig, fs_ch, target_fs))

        # Trim to the shortest channel length (safety measure)
        min_len = min(len(c) for c in channels)
        data = np.vstack([c[:min_len] for c in channels])  # [C x N]

        # Downstream tools expect exactly 3 channels (Lead I, II, III).
        # If fewer are available, duplicate the first channel to fill.
        if data.shape[0] < 3:
            reps = int(np.ceil(3 / data.shape[0]))
            data = np.vstack([data] * reps)[:3, :]

        n_ch = data.shape[0]

        # Attempt to get the recording start date
        try:
            start_dt = reader.getStartdatetime()
            date_str = start_dt.strftime("%m/%d/%Y")
        except Exception:
            date_str = "01/01/2024"

        # Patient ID from filename
        patient_id = edf_path.stem
    finally:
        reader.close()

    # ── Chunk into 24-hour segments ──────────────────────
    total_samples = data.shape[1]
    chunk_samples = chunk_hours * 3600 * target_fs
    n_chunks = int(np.ceil(total_samples / chunk_samples))

    out_dir = Path(output_dir) / patient_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # Default lead labels — always 3 leads for downstream tools
    if lead_labels is None:
        lead_labels = ["Lead I", "Lead II", "Lead III"]

    for i in range(n_chunks):
        start = i * chunk_samples
        end = min((i + 1) * chunk_samples, total_samples)
        chunk = data[:, start:end]  # [C x chunk_len]

        base_name = f"{patient_id}_part{i + 1}_ECG_fs{target_fs}"

        # ── Write .bin ────────────────────────────────────
        # MATLAB fwrite of [C×N] writes column-major:
        #   ch0[0], ch1[0], ch2[0], ch0[1], ch1[1], ch2[1], …
        # Transpose to [N×C] then write in C order to match.
        bin_path = out_dir / f"{base_name}.bin"
        with open(bin_path, "wb") as f:
            f.write(chunk.T.astype(np.float64).tobytes(order="C"))

        # ── Write .inf ────────────────────────────────────
        inf_path = out_dir / f"{base_name}.inf"
        n_samples = chunk.shape[1]
        with open(inf_path, "w", newline="\r\n") as f:
            f.write(f"Patient = {patient_id}\n")
            f.write(f"Description = Resampled ECG Data\n")
            f.write(f"Export Date = \n")
            f.write(f"Number of Channel = {n_ch}\n")
            f.write(f"Points for Each Channel = {n_samples}\n")
            f.write(f"Data Sampling Rate = {target_fs} points/second\n")
            f.write(f"Start Time = {date_str} 12:00:00 AM\n")
            f.write(f"Stop Time = \n")
            f.write(f"Units: \n")
            f.write(f"Channel Number  Channel Label\n")
            for idx, lbl in enumerate(lead_labels[:n_ch], start=1):
                f.write(f"{idx}              {lbl}\n")

        print(f"  -> Saved: {base_name}  ({n_samples:,} samples)")

    print(f"\nDone — {n_chunks} chunk(s) written to {out_dir}")


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────


def main():
    input_dir = Path(INPUT_DIR)
    output_dir = Path(OUTPUT_DIR)

    edf_files = []
    for folder in sorted(input_dir.iterdir()):
        if not folder.is_dir() or folder.name.startswith("."):
            continue
        for f in folder.iterdir():
            if f.is_file() and f.suffix.lower() == ".edf":
                edf_files.append(f)
    edf_files.sort()
    if not edf_files:
        sys.exit(f"No EDF files found in subfolders of {input_dir}")

    print(f"Found {len(edf_files)} EDF file(s) in {input_dir}\n")

    for edf_path in edf_files:
        print(f"{'═' * 60}")
        print(f"Processing: {edf_path.name}")
        print(f"{'═' * 60}")
        try:
            edf_to_bin(
                edf_path=str(edf_path),
                output_dir=str(output_dir),
                target_fs=1000,
                chunk_hours=24,
                target_signals=["ECG_1", "ECG_2", "ECG_3"],
                lead_labels=None,
            )
        except Exception as e:
            print(f"  !! Error processing {edf_path.name}: {e}\n")

    print(f"\n{'═' * 60}")
    print("All files processed.")


if __name__ == "__main__":
    main()
