import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv(
    "Z1011347_20210625T184500_20210626T024500.dat",
    usecols=["Index", "NLS_NOM_ECG_ELEC_POTL_I", "NLS_NOM_PULS_OXIM_PLETH"],
    header=2,
    nrows=4000,
)

ppg_timestamp_idx = pd.to_numeric(df.iloc[:, 0], errors="coerce")
ppg = pd.to_numeric(df.iloc[:, 2], errors="coerce")

valid = ppg_timestamp_idx.notna() & ppg.notna()
x_uniform = np.linspace(0, ppg_timestamp_idx[valid].max(), valid.sum())

num_empty = ppg.isna().sum()
num_nonempty = ppg.notna().sum()

print(f"Empty: {num_empty}")
print(f"Non-empty: {num_nonempty}")

fig, axs = plt.subplots(3, 1, figsize=(18, 6), sharex=True)

axs[0].scatter(ppg_timestamp_idx[valid], ppg[valid], s=1)
axs[0].set_ylabel("Timestamp: accurate")

axs[1].scatter(x_uniform, ppg[valid], s=1)
axs[1].set_ylabel("Timestamp: shifted")


axs[2].scatter(range(2000), df.iloc[:, 1], s=1)
axs[2].set_ylabel("ECG CHAN 1")

plt.tight_layout()
plt.savefig("plot.png", dpi=300, bbox_inches="tight")
