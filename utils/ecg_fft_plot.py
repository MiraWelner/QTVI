#!/usr/bin/env python3
"""ECG raw and ECG upsampled amplitude spectra: 0-100 Hz, plus 42-46 Hz and 55-65 Hz zooms.
    python ecg_fft_plot.py 7018053_20110922_01h15m55s.csv
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def spectrum(v, fs, pad=8):
    x = v - np.mean(v)
    w = np.hanning(x.size)
    nfft = 1 << int(np.ceil(np.log2(x.size * pad)))
    X = np.fft.rfft(x * w, nfft)
    return np.fft.rfftfreq(nfft, 1.0 / fs), 2.0 * np.abs(X) / w.sum()


path = Path(sys.argv[1])
d = pd.read_csv(path)
fig, ax = plt.subplots(2, 3, figsize=(19, 7))
for i, (sig, label) in enumerate([("ECG1_raw", "ECG raw"), ("ECG1_upsampled", "ECG upsampled")]):
    t = d[f"{sig}_t_ms"].dropna().to_numpy(float)
    v = d[f"{sig}_mv"].dropna().to_numpy(float)
    fs = 1000.0 / np.median(np.diff(t))
    f, a = spectrum(v, fs)
    for j, (lo, hi) in enumerate([(0, 100), (42, 46), (55, 65)]):
        m = (f >= lo) & (f <= hi)
        ax[i, j].semilogy(f[m], a[m] + 1e-12, lw=0.8 if j == 0 else 1.0)
        ax[i, j].set_xlim(lo, hi)
        ax[i, j].set_title(f"{label} ({fs:.0f} Hz)  -  {lo}-{hi} Hz", fontsize=10)
        ax[i, j].set_xlabel("frequency (Hz)", fontsize=9)
        ax[i, j].set_ylabel("amplitude (mV)", fontsize=9)
        ax[i, j].grid(alpha=0.3, which="both")
        if lo <= 60 <= hi:
            ax[i, j].axvline(60, color="crimson", ls="--", lw=0.8)
        if lo <= 44.52 <= hi:
            ax[i, j].axvline(44.52, color="darkorange", ls="--", lw=0.8)
        ax[i, j].tick_params(labelsize=8)
fig.suptitle(f"{path.name}: ECG amplitude spectra (Hann window); 60 Hz in red, 44.52 Hz in orange", fontsize=11)
fig.tight_layout()
out = path.with_name(path.stem + "_ecg_fft.png")
fig.savefig(out, dpi=110)
print(f"plot: {out}")
