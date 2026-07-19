# QTVI script Pseudocode

Mira Welner

---

## 0. Config and Dataset collection

1. A file called `config.csv` in the same directory as the `.exe` file is loaded. It has 3 rows corresponding to the 3 datasets that this script is capable of loading: MESA, Bittium, and CHAOS. The columns are:

| Column                 | Type   | Meaning                                                                                                                                                                                                              |
| ---------------------- | ------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `data_type`            | string | Dataset key; row is used when this equals the selected dataset. Compared upper-cased.                                                                                                                                |
| `main_file_extention`  | string | File extension of main loaded file - .dat for CHAOS, .edf for others                                                                                                                                                 |
| `sleep_file_extention` | string | File extension of sleep state files - .XML for MESA                                                                                                                                                                  |
| `sleepstate_length`    | double | Length of time that a sleepstate takes up - 1/HZ                                                                                                                                                                     |
| `blanking_period`      | double | Time in milliseconds after the previous R peak in which a new R peak can NOT occur                                                                                                                                   |
| `threshold`            | double | The 'high' is the median R peak in the previous 10 second section, the 'low' is the median of the current 10s section. The threshold of height at which a peak is a possible R peak is $threshold \times (high-low)$ |
| `bin_size_minutes`     | double | Bin length for annealing and templating in minutes                                                                                                                                                                   |
| `original_file_path`   | string | Input location; if blank, prompted.                                                                                                                                                                                  |
| `output_folder`        | string | Parent output dir; subfolders derived from it. If blank, prompted.                                                                                                                                                   |

### Per-channel sample rates

Every channel contributes two columns: `<prefix>_raw_rate` (native Hz) and `<prefix>_upsampled_rate` (target Hz). A `0`/blank rate means the channel is absent. The channels and their column prefixes are:

| Column prefix     | Channel                        |
| ----------------- | ------------------------------ |
| `ecg`             | All 3 ECG leads (or 1 in MESA) |
| `ppg`             | Photoplethysmogram             |
| `cvp`             | Central venous pressure        |
| `pres`            |                                |
| `abp`             | Ambulatory Blood Pressure      |
| `art`             | Arterial Pressure              |
| `art_pulm`        | Arterial Pulmonary Pressure    |
| `accel`           | Accelerometer                  |
| `temp`            | Temperature                    |
| `marker`          |                                |
| `resp`            | Respiration                    |
| `pacemaker_event` | Pacemaker events               |
| `eeg`             | Electroencephalogram           |
| `eogl`            | Left electrooculogram          |
| `eogr`            | Right electrooculogram         |
| `emg`             | electromyography               |
| `flow`            |                                |
| `snore`           |                                |
| `thor`            |                                |
| `abdo`            |                                |
| `leg`             |                                |
| `auxac`           |                                |
| `therm`           |                                |
| `pos`             |                                |
| `oxstatus`        |                                |
| `spo2`            |                                |
| `hr`              | Heart rate                     |
| `dhr`             |                                |

If I do not put an explanation in the channel meaning it is because I'm unsure of the biological meaning

---

## 1. Noise markings

The first step for the user is to mark the noise and other phenomena on the noise marking GUI. The GUI is formatted as follows depending on dataset:

### GUI design

```
MESA
+------------+--------------------------------------------+
| controls   | [ ECG ampogram ---------------------------]|
|            | [ PPG ampogram -----] [  -----------------]|
|            | [ Sleep states ---------------------------]|
|            +--------------------------------------------+
|            |  ECG1(EKG)  <- markable                    |
|            |  PPG(Pleth) <- markable                    |
+------------+--------------------------------------------+

BITTIUM
+------------+--------------------------------------------+
| controls   | [ ECG ampogram ---------------------------]|
|            | [ Temperature ------] [  -----------------]|
|            | [ Marker ---] [ Resp -----] [ Pacer ------]|
|            +--------------------------------------------+
|            |  ECG1       <- markable                    |
|            |  ECG2       <- markable                    |
|            |  ECG3       <- markable                    |
|            |  Accel      <- markable                    |
+------------+--------------------------------------------+

CHAOS
+------------+--------------------------------------------+
| controls   | [ ECG ampogram ---------------------------]|
|            | [ PPG ampogram -----] [  -----------------]|
|            +--------------------------------------------+
|            |  ECG1       <- markable                    |
|            |  ECG2       <- markable                    |
|            |  ECG3       <- markable                    |
|            |  PPG        <- markable                    |
|            |  ABP        <- markable                    |
|            |  ART        <- markable                    |
|            |  ART_PULM   <- markable                    |
+------------+--------------------------------------------+
```

### Potential Markings

The markings that the user can make are, ordered by the number that represents them in the log and output files, are:

1. R Peak Noise — noise so significant that R peaks are not detected in the marked range. Nothing is labled as R peak noise in the output log because the log is beat by beat and R peak noise has no beats.

2. Minor noise — Indicitive of noise, but does not alter R peak detection in any way

3. Blank.+Thresh. — Labeling a region with this marker will create a popup box that allows you to edit the blanking and thresholds for that particular reason.

4. Invert/Noninvert — The signal in this region is analyzed as inverted signal, unless the checkbox for inverted signal is checked, in which case this region is uninverted.

5. PVC

6. PAC

7. Cond. Delay

8. AF

9. SVT

10. VT

11. Benign Arr.

12. Sig. Arr.

13. Other

### Output bin file

This will be found in `ouput_folder/noise_marking_output/[ID]_noise_markings.bin`

```
uint64  count                      # number of segments
repeat count times:
    double  start_sample
    double  end_sample
    double  start_sec              # = start_sample / sampleRate
    double  end_sec                # = end_sample   / sampleRate
    double  channel_code           # PPG=1, ECG1=2, ECG2=3, ECG3=4, ABP=5, ACCEL=6, ART=7, ART_PULM=8
    double  annotation_code        # The above mentioned noise markings
```

### Output CSV.

This will be found in `ouput_folder/noise_marking_output/[ID]_noise_markings.bin`

```
header:  start_sample,end_sample,start_sec,end_sec,label,marking_type
row:     one row per user created annotation
```

---

## 2. Annealing

1. Signal is split into segments of size `bin_size_minutes`

2. Noise is removed only if every markable channel that is present is marked with `1` ie. R peak noise in that particular location. An exception to this is acceleration data which doesn't effect what is removed.

3. The resulting fragments will be either `bin_size_minutes` at max, or whatever remains after the noise was removed form the fragment of `bin_size_minutes`

4. The clean fragments which are greater or equal to half the `bin_size_minutes` are kept, while ones less that this are merged into adjacent epochs. 

5. Adjacent epoch is determined by which half of the bin remains. If the midpoint of the bin is closer to the temporal location of the previous bin, it is merged with the previous bin. Else, it is merged with the subsequent bin.

6. If there is a too-small bin which should be merged with the previous bin on the start of the signal, or one which should be merged with the subsequent bin on the end of the signal, it is discarded.

The bin looks like:

```
# ---- header ----
uint64   nSegments                      # = n
if nSegments > 0:                        # NOTE: these 3 are written ONLY when n>0
    double   ppgSampleRate               # upsampled PPG rate
    double   ecgSampleRate               # upsampled ECG rate
    double   scoringEpochSec             # sleep-stage epoch length (s)
uint32   nChannels                       # = 36
float32  nativeRates[36]                 # per-channel native/raw rate

# ---- per segment  (nSegments of these) ----
    # index pairs (1-based sample ranges), each: uint64 count, then count x (uint64 first, uint64 second)
    uint64 cnt;  {uint64 first, uint64 second} x cnt      # ppg_bin_indexs
    uint64 cnt;  {uint64 first, uint64 second} x cnt      # ecg_bin_indexs

    # primary/secondary signal slices, each: uint64 sz, then sz x double
    uint64 sz;  double ppg[sz]
    uint64 sz;  double ecg1[sz]
    uint64 sz;  double ecg2[sz]
    uint64 sz;  double ecg3[sz]
    uint64 sz;  double sleep_stages[sz]

    # then all 36 pass-through channels, each as { upsampled slice, raw slice }:
    repeat ch = 0..35:
        uint64 sz;      double  upsampled[sz]             # proportional-index slice
                                                          #   missing channel -> sz=1, value -1.0 (sentinel)
        uint64 nPairs;  double  raw[2*nPairs]             # (t, v) interleaved, t = abs seconds
                                                          #   missing channel -> nPairs=1, pair (-1.0, -1.0)
```

This will be found in `ouput_folder/\annealed_output/[ID]_annealed.bin`

Notes

- PPG being first is an artifact of the code being based off Daniel's code which cared a lot about PPG
* Index pairs are **1-based, inclusive** `[first, second]` sample ranges on the
  channel's own grid; a segment may hold several pairs (redistributed fragments).
- The upsampled slices are stored as series, split by raw temporal location not taking into account the upsampling rate.
* The raw signal is stored as  `(t,v)` pairs whose `t` falls in the segment's ECG time window, timestamps in absolute seconds-from-recording-start. This is an artifact of when I was trying to mark timeskips, but I removed that functionality when I learned that CHAOS .dat files sometimes have issues where the signal is 2xsampling rate for a bit and absent for a bit, so it may be unnecessary now.
- Missing channels use sentinels (upsampled `[-1.0]`; raw `[(-1.0,-1.0)]`) so the
  schema stays uniform across all 36 slots.

---

## 3. R-peaks Detection

R-peak detection is conducted via weighted multi-algorithm consensus.

1. The signal is linearly detrended via subtracting the mean of each segment from every sample in the segment

2. The six algorithms used are:
   
   1. 3 passes of a Hamilton/Tompkins-style detector at sensitivity levels 0.1,0.2 and 0.4
      
      1. Bandpass filter (order-2 zero-phase `filtfilt`, ~0.05–100 Hz)
      
      2. Differentiate
         3- Each pass then thresholds that integrated signal at a **fraction of `max_h`** — **0.2**, **0.1**, and **0.4** — and localizes the accepted peaks back on the band-passed trace.
         4- The fraction is the sensitivity knob: **0.1** = most sensitive (more detections, more false positives), **0.4** = strictest (fewer, higher-confidence), **0.2** = middle.
      
      3. Square
      
      4. Moving sum integration (window ~7 samples, scaled by sample rate
      
      5. Median filter
      
      6. Remove filter/integration delay
      
      7. Estimate the peak level of integrated signal, 
      
      8. Each pass then thresholds that integrated signal at a fraction of the estimated peak level:  **0.2**, **0.1**, and **0.4** — and localizes the accepted peaks back on the band-passed trace.
   
   2. One pass of a Pan-Tompkins detector
      
      1. Bandpass 5–15 Hz
      
      2. derivative
      
      3. square → **moving-window integration of 150 ms**.
      
      4. Dual adaptive thresholds on the integrated and filtered signals (running signal/noise estimates, adaptation factor 0.125), with a **200 ms** minimum peak distance (refractory).
      
      5. T-wave discrimination**: when an RR interval is shorter than **360 ms**, it compares the candidate's slope (over a 75 ms window) to the previous QRS; if the slope ratio is below **0.5**, it's rejected as a T-wave rather than a beat.
      
      6. **Search-back**: if too long passes with no detection, it re-scans with a lowered threshold (200 ms margin; fast/slow recovery factors 0.25/0.75) to recover missed beats.
   
   3. One pass of a LMS/moving-integration detector
      
      1. Mean-subtract Butterworth bandpass (`filtfilt`) 
      
      2. Square the differenced signal → moving-window integration of ~175 ms.
      
      3. Adaptive threshold as a fraction of the running MWI maximum: fraction 0.25 with a forgetting/decay factor of 0.80.
      
      4- 250 ms refractory between accepted peaks.
   
   4. One pass of a squared-signal threshold detector.
      
      1. Square the raw ECG (no bandpass/integration).
      
      2. Static threshold** = mean + **2·standard deviation** of the squared signal.
      
      3. Runs `findpeaks` with a **minimum peak spacing ≈ half the median inter-beat interval**, then keeps only peaks above the threshold.
         5- The pan tompkins, LMS, and squared-signal detector are run through a peak-snapping algorithm
      
      4. Take the consecutive gaps between that detector's detections, `d[i] = idx[i] − idx[i−1]`, and compute their medians
      
      5. Create a window set to 1/6th of the median between-detection gaps
      
      3- Snap that detection to the location of that maximum, i.e. move `idx[i]` to the local peak apex; the original index is kept if nothing in the window is larger.
      
      4- Return the refined index list.
      
      5. In each of these windows — centered on a detection and spanning `[idx[i] − half_win, idx[i] + half_win]`, clamped to the signal bounds — find the sample with the maximum ECG amplitude.

3. After this, each sample that is selected by one of the algorithms as a 'peak' is weighted via the following:
   
   1. Hamilton/Tompkins algorithm with 0.2 threshold: **0.75**
   
   2. Hamilton/Tompkins algorithm with 0.1 threshold: **0.25**
   
   3. Hamilton/Tompkins algorithm with 0.4 threshold: **0.25**
   
   4. Pan-Tompkin: **0.25**
   
   5. LMS: **1.5**
   
   6. Squared-signal: **0.75**

4. A candidate is accepted as a true R-peak only if its summed weight reaches a consensus threshold of **2.4**

5. The same process above is run for the square of the signal and the absolute value of the signal

### Outputs

The binary file can be found at: `output_path\r_peak_finding_output\4011465_20110822_peak_locations_all_beats.bin`

```
*Header**
- `uint64 numBins`

**Per bin** (in this exact order):

1. **R-peak index arrays — 9** (3 ECG channels × 3 preprocessing methods).
   Each: `uint64 count` + `count × uint64` indices. Order:
   - `ch1.raw`, `ch1.squared`, `ch1.absval`
   - `ch2.raw`, `ch2.squared`, `ch2.absval`
   - `ch3.raw`, `ch3.squared`, `ch3.absval`

2. **PPG event index arrays — 2** (same `count + uint64[]` layout):
   - `ppgMaxAmps`, then `ppgMinAmps`

3. **Preprocessed signal traces — 6** (each `uint64 count` + `count × double`):
   - `ch1.squared`, `ch1.absval`, `ch2.squared`, `ch2.absval`, `ch3.squared`, `ch3.absval`
   - (kept on disk so template rebuilds don't recompute them)

4. **Noise flags — `uint8 flags[9]`** (9 raw bytes), one per channel×method,
   same order as the index arrays (ch1 raw/squared/absval, ch2 …, ch3 …);
   each is the `*_noisy` boolean.

5. **ECG–PPG pairs** — `uint64 numPairs`, then `numPairs × (int64 ppg, int64 ecg)` interleaved.
```

Notes:

All indices are 1-based on disk (writer adds 1 to every R-peak index,

In the pairs, an unpaired side (NaN / negative sentinel) is written as **`-1`** (not `+1`-shifted).

The CSV output can be found at: `csv_for_analysis/<stem>_peak_locations_all_beats.csv`

| Column             | Meaning                                             |
| ------------------ | --------------------------------------------------- |
| `beat`             | Row/beat index                                      |
| `<chan>_x`         | Time of detected peak, per channel                  |
| `<chan>_y`         | Amplitude of detected peak, per channel             |
| `blanking_<chan>`  | Blanking period active at detection                 |
| `threshold_<chan>` | Detection threshold active at detection             |
| `marked_<chan>`    | Annotation type code on that beat (0 = none)        |
| `post_<chan>`      | Post-arrhythmia tag (e.g. post-PVC), else 0         |
| `inverted_<chan>`  | 1 if channel polarity inverted at that beat, else 0 |

---

## 4. Templating

### ECG templating algorithm (R-aligned and Q-aligned)

- Pass 0 = R-aligned (`g_q_align = false`)
- Pass 1 = Q-aligned (`g_q_align = true`), built on top of the R-aligned result

**Shared steps (both passes)**

- Slice per beat: for each consecutive R-peak pair, cut `[R_i - 0.3*RR, R_i + 1.5*RR]` (RR = interval to the next R peak); beat length varies with RR
- Tukey outlier rejection (`[Q1 - k*IQR, Q3 + k*IQR]`, k = 1.5): drop beats whose length is an outlier, or whose R-column amplitude is an outlier
- Horizontal R-alignment: find the longest surviving RR, size the shared axis from it, right-shift every beat (NaN-padded) so its own R lands on a common anchor column
- Reference beat = the first beat whose length equals the median RR length
- Vertical DC alignment: measure each beat's PR-segment baseline (small window just before R), shift each beat vertically to match the reference beat's baseline
- Template = column-wise, NaN-skipping median across all aligned beats (NaN cells past a beat's real extent don't participate); per-sample std computed the same way (ddof = 1), only for the "raw" method

**Q-aligned pass -- extra step**

- Starting from the already R-aligned beats, find each beat's Q point: locate the steepest rising slope on the R-upstroke, walk left until slope drops below 20% of that max -> Q
  - Fallback: 10% of median RR before the upstroke onset, if no Q found
- Take the largest Q column across beats as the common target column
- Shift each beat right by `target - its_own_Q_column` (NaN-padded) so every beat's Q lands on the same column
- Template = column-wise median of these Q-shifted beats

### PPG / arterial channels (unaffected by the R/Q-align flag)

- Slicing driven by the same ECG R-peak pairs (`ch1.raw`), converted to each channel's own sample rate, with a fixed real-time pad (0.25-0.3s) so the first R always lands at a fixed column
- Within that window: find systolic peak (max) and foot (min before peak), then locate the 50%-upslope (half-height) crossing between them -- this is the horizontal alignment fiducial (well-localized on the steep upstroke, unlike the peak or foot)
- Vertical alignment matches each beat's foot-baseline to the reference beat's
- Template = column-wise, NaN-skipping median, same as ECG

## 5. Template Marking

The templates are loaded and the autodetection works as follows:

### Autodetection algorithms (FeatureMarks)

**R peak (fixed, auto-only)**

- `r_peak()`: finds argmax|v - baseline| over the first 3/4 of the beat window; baseline is the mean of that same window. Returns the index and whether the deflection is positive or negative, so downstream detectors can flip the trace to a consistent "upright" orientation
- `detect_r_peak()` wraps this for seeding

**ECG Human Placed Markers can be dragged, but they are automatically seeded to land somewhere in the plot**

- **Q begin**: finds the lowest point (trough) before R, then walks left from that trough until the first derivative goes non-negative (the trough's upslope onset)
- **P peak**: max value in the region before Q begin
- **S end**: finds the S trough after R, takes a baseline from a short window ~50-100 samples past R, then walks right from S until the signal recovers to 90% of the S-to-baseline depth
- **T begin**: first local max right after S end (restricted to the first 2/3 of the beat, picking the tallest candidate so an early ST bump isn't mistaken for the T wave), then walks left to its foot
- **T end**: same T-wave peak search, then walks right from the peak to its foot (must stay right of S end)

**ECG X markers**

- **P wave** = max within ±0.05s of the P marker
- **T end** = passthrough of the user's T-end marker
- **S end** = recovery knee within ±0.05s of the S marker
- **T peak** = max between T begin/T end
- **Q onset** = curvature knee of a cubic fit within ±0.05s of the Q marker
- **R wave** = argmax|v-baseline| over [Q begin, S end]

**PPG movable markers** (all operate on the pulse waveform)

- **Onset (foot)**: minimum before the systolic peak
- **Peak**: global max of the trace
- **End**: minimum after the peak
- **Dicrotic notch**: first interior local minimum between peak and end (falls back to 1/3 of the way from peak to end if none found)
- **T80**: point between foot and peak closest to 80% of the amplitude rise

Two bin files are output:

This can be found in `outputs_folder/template_outputs/[ID]_templates.bin`

| Field                                               | Meaning                                                                                                        |
| --------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `n_bins` (file header)                              | Total number of bins in the file                                                                               |
| `ch{1-3}_raw`                                       | ECG template (waveform + std + alignment point + avg_r_expand) built from the raw/filtered signal, per channel |
| `ch{1-3}_squared`                                   | Same, built from the squared signal (no std stored)                                                            |
| `ch{1-3}_absval`                                    | Same, built from the absolute-value signal (no std stored)                                                     |
| `ch{1-3}_unfiltered`                                | Same, built from the original unfiltered signal (no std stored)                                                |
| `ppgTemplate`, `ppgTemplate_std`                    | PPG template waveform and its per-sample std                                                                   |
| `abpTemplate`, `abpTemplate_std`                    | ABP background template and per-sample std (empty if channel absent)                                           |
| `artTemplate`, `artTemplate_std`                    | ART background template and per-sample std (empty if channel absent)                                           |
| `artPulmTemplate`, `artPulmTemplate_std`            | ART_PULM background template and per-sample std (empty if channel absent)                                      |
| `ch1_n_beats_raw, ch2_n_beats_raw, ch3_n_beats_raw` | Number of beat slices that survived drop rules and fed each channel's raw-method median                        |
| `ppg_n_beats`                                       | Number of beat slices that fed the PPG median                                                                  |
| `bad_segment`                                       | 1 if the whole bin is flagged unusable, else 0                                                                 |

This can be found in `outputs_folder/qtvi_marker_path/[ID]_template_mark.bin``

It starts out as `outputs_folder/qtvi_marker_path/[ID]_template_mark_r_align.bin`` but when the q align overwrites it it just becomes template_mark.

| Field                                                                                | Meaning                                                      |
| ------------------------------------------------------------------------------------ | ------------------------------------------------------------ |
| `index`                                                                              | Bin number                                                   |
| `bad_r_ch1/2/3`                                                                      | Whether R-detection is flagged bad on that ECG channel       |
| `ppg_issue`                                                                          | 0 = ok, 1 = bad, 2 = no PPG                                  |
| `p_peak_ch1-3, q_begin_ch1-3, r_peak_ch1-3, s_end_ch1-3, t_begin_ch1-3, t_end_ch1-3` | Sample index of each ECG marker, per channel (-1 = unmarked) |
| `ppg_onset, ppg_p50, ppg_peak, ppg_dicrotic, ppg_peak2, ppg_t80, ppg_end`            | Sample index of each PPG marker                              |
| `abp_issue`                                                                          | 0 = ok, 1 = bad, 2 = absent                                  |
| `abp_onset, abp_peak, abp_dicrotic, abp_peak2, abp_end`                              | Sample index of each ABP marker                              |
| `art_issue`                                                                          | 0 = ok, 1 = bad, 2 = absent                                  |
| `art_onset, art_peak, art_dicrotic, art_peak2, art_end`                              | Sample index of each ART marker                              |
| `art_pulm_issue`                                                                     | 0 = ok, 1 = bad, 2 = absent                                  |
| `art_pulm_onset, art_pulm_peak, art_pulm_dicrotic, art_pulm_peak2, art_pulm_end`     | Sample index of each ART_PULM marker                         |

And the CSV outputs are:

For every `<chan>` in:

- `ch1` — ECG channel 1

- `ch2` — ECG channel 2

- `ch3` — ECG channel 3

- `art_pulm` — Pulmonary artery

- `abp` — Arterial blood pressure

- `art` — Arterial line

- `ppg` — PPG

- 

`output_folder/csv_for_analysis/template.csv`

| Column pattern                                                                          | Meaning                                                                                     |
| --------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| `file_id`                                                                               | ID of subject                                                                               |
| `bin_num`                                                                               |                                                                                             |
| `x_ms`                                                                                  | Time in miliseconds                                                                         |
| `<chan>_raw_mv_{r/q}`                                                                   | Raw amplitude at this sample, per channel (ch1, ch2, ch3, ppg, abp, art, art_pulm)          |
| `<chan>_Normalized_{r/q}`                                                               | Normalized amplitude at this sample, per channel                                            |
| `<chan>_raw_std_{r/q}`                                                                  | Per-sample std of raw amplitude, per channel                                                |
| `<chan>_normalized_std_{r/q}`                                                           | Per-sample std of normalized amplitude, per channel                                         |
| `<pt>_ch{1-3}_location_{autodetect/user}_{r/q}`                                         | 1 if this row's sample index is that marker's position, else blank (per ECG marker/channel) |
| `<marker>_location_{autodetect/user}_{r/q}` (ppg/abp/art/art_pulm)                      | Same one-hot marker flag for pulse-channel markers                                          |
| `p_wave_ch{1-3}, q_onset_ch{1-3}, r_wave_ch{1-3}, t_peak_ch{1-3}` (auto only, `_{r/q}`) | One-hot flags for the autodetected computed ECG glyphs                                      |
| `ppg foot, p1` (auto only, `_{r/q}`)                                                    | One-hot flags for the autodetected computed PPG glyphs                                      |

# 
