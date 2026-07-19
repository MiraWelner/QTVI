import numpy as np, pandas as pd, glob
import plot_templates as A

f = sorted(glob.glob("templates_healthy/*_template.csv"))[0]
df = pd.read_csv(f)

# RAW: first bin, straight from CSV, no processing
changed = df["bin_num"].ne(df["bin_num"].shift())
g0 = [x for _, x in df.groupby(changed.cumsum(), sort=False)][0].reset_index(drop=True)
yraw = pd.to_numeric(g0["ch1_Normalized_r"], errors="coerce").to_numpy()
def qrs_width(y):
    b = np.nanmedian(y); d = np.abs(y - b); r = np.nanargmax(d)
    half = d[r] / 2
    return int(np.sum(d[max(0,r-150):r+150] > half))
print("RAW bin0 QRS width (samples):", qrs_width(yraw), " R val:", round(np.nanmax(np.abs(yraw-np.nanmedian(yraw))),1))

# AFTER my pipeline
beats = A.collect_beats(df, "ch1", True)
print("beats kept:", len(beats))
print("first kept beat QRS width:", qrs_width(beats[0]["y"]))
mat, anch = A.align_r(beats)
print("aligned matrix beat0 QRS width:", qrs_width(mat[0][~np.isnan(mat[0])]))
