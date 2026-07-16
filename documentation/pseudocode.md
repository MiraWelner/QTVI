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

2. There 3 passes of a Hamilton/Tompkins-style detector at sensitivity levels 0.1,0.2 and 0.4

3. One pass of a Pan-Tompkins detector

4. three passes of a Hamilton/Tompkins-style detector (`rpeakdetect`) at progressively different sensitivity thresholds (0.2, 0.1, 0.4 of a reference amplitude),

5. a Pan–Tompkins detector,

6. an LMS/moving-integration detector, and

7. a squared-signal threshold detector.

Each detector contributes candidate peak locations with a fixed reliability weight (0.75, 0.25, 0.25, 1.25, 1.5, 0.75 respectively), so more trusted algorithms count more heavily. Detections from the last three detectors are first refined by snapping each to the local ECG maximum within a window of ±(median RR)/6 samples, correcting small localization differences between algorithms. All candidate peaks are then pooled; detections from different algorithms falling within 8 samples of one another are treated as the same beat, coalesced to the sharper (higher-amplitude) sample, and their weights summed. A candidate is accepted as a true R-peak only if its summed weight reaches a consensus threshold of **2.4** — i.e. it must be supported by several independently-weighted detectors, not just one. This suppresses the idiosyncratic false positives of any single method while retaining beats that most methods agree on.

The individual detectors follow the standard QRS-enhancement pipeline — band-pass filtering to isolate the QRS frequency band, differentiation, squaring/rectification, and moving-window integration — differing in band edges, integration length, and thresholding. Representative parameters: Pan–Tompkins uses a 5–15 Hz band, 150 ms integration, adaptive dual signal/noise thresholds, a 200 ms refractory period, T-wave discrimination by slope for RR < 360 ms, and missed-beat search-back; `rpeakdetect` band-passes 0.05–100 Hz, integrates over ~7·fs/256 samples, and thresholds relative to the maximum integrated amplitude over the central 50% of the segment. In all cases the final peak position is refined on a cleaner version of the signal (the band-passed or raw trace) rather than on the squared/integrated envelope, so reported R-peak times coincide with the true QRS apex.

Detection runs per channel and per method; the pipeline additionally computes squared and absolute-value variants of the signal and detects on those, and pairs each ECG R-peak with the corresponding PPG pulse for cross-modal analysis. Detected peaks, their amplitudes (measured against a local PR-segment baseline), and the resulting RR intervals are exported per beat.

---

Two honest caveats for a paper:

- **"Consensus threshold 2.4" is a tuned constant**, not a probabilistic quantity — with the given weights it roughly means "agreement from the equivalent of the strong detectors, or several weaker ones." If you report it, state the weights alongside so the 2.4 is interpretable, and ideally note how it was tuned/validated.
- The **T-wave discrimination, refractory, and search-back parameters** (360/200 ms, RR bounds) are Pan–Tompkins-internal defaults; if your ECG sampling rate or population differs materially from the assumptions, those are the parameters a reviewer will ask whether you re-validated.

---

## 6. Templates — `<stem>_templates.bin` (+ `.partial.bin`) & snips `.csv`

Written by `template_io::write_template_binfile` / `write_snips_csv`. The **`.partial.bin`** is the provisional fast build the viewer opens; the canonical `.bin` is (re)written by `finalizeViewerJob`, then the partial is deleted.

### 6a. `.bin`

```
uint64  n_bins
repeat per bin:
    # 12 ECG method blocks, fixed order:
    #   ch1_raw, ch1_squared, ch1_absval, ch1_unfiltered,
    #   ch2_raw, ch2_squared, ch2_absval, ch2_unfiltered,
    #   ch3_raw, ch3_squared, ch3_absval, ch3_unfiltered
    for each block:
        uint64 sz;  double ecgTemplate[sz]
        uint64 sz;  double ecgTemplate_std[sz]      # sz=0 for non-raw methods
        double alignment_point
        double avg_r_expand

    # PPG template + std:
    uint64 sz;  double ppgTemplate[sz]
    uint64 sz;  double ppgTemplate_std[sz]

    # arterial background templates (foot-anchored, no std), sz=0 if absent:
    uint64 sz;  double abpTemplate[sz]
    uint64 sz;  double artTemplate[sz]
    uint64 sz;  double artPulmTemplate[sz]

    # per-channel surviving beat counts:
    uint64  ch1_n_beats_raw
    uint64  ch2_n_beats_raw
    uint64  ch3_n_beats_raw
    uint64  ppg_n_beats

    uint8   bad_segment
```

Only the ECG **raw** method writes a populated `std` vector; squared/absval/
unfiltered write `sz=0`. (Note: the header comment lists arterial std vectors;
the authoritative field order is the writer above — treat arterial std as
present-only-when-nonempty and verify against `write_template_binfile` if you
extend it.)

### 6b. snips `.csv` (`<stem>_template_snips.csv`, in `csv_for_analysis/`)

One **column per retained beat** ("snip"); rows are sample indices.

```
header:  {CHAN}_{n}, ...          # CHAN order: CH1,CH2,CH3,PPG,ABP,ART,ART_PULM; n = 0-based snip index
row s:   value at sample s for each snip column (blank if s >= that snip's length or NaN)
```

Bad segments (`bad_segment[b]`) are skipped. Column count = total non-bad snips
across all bins/channels; row count = longest snip.

---

## 7. Beat log — `<stem>_log.csv` (`beat_log.hpp`)

One row per beat index; 8 channels per group
(ecg1, ecg2, ecg3, ppg, abp, art, art_pulm, accel).

```
header:
  beat,
  ecg1_x,ecg1_y, ecg2_x,ecg2_y, ecg3_x,ecg3_y, ppg_x,ppg_y,
  abp_x,abp_y, art_x,art_y, art_pulm_x,art_pulm_y, accel_x,accel_y,
  blanking_ecg1,threshold_ecg1, ... blanking_accel,threshold_accel,
  marked_ecg1..marked_accel,        # annotation code covering this beat (per channel)
  post_ecg1..post_accel,            # code of the PREVIOUS eligible arrhythmia (post-beat)
  inverted_ecg1..inverted_accel     # inversion flag per channel
row: beat index, then the above fields (x = sample, y = amplitude)
```

`marked_*` / `post_*` use the annotation `code` field (§8). Fully-empty rows
(no beat on any channel) are skipped.

---

## 8. Annotation codes (`annotation_types.hpp`, `noise_types[13]`)

`codeFor(label)` returns the `code` field below; used by the noise `.bin`/`.csv` (§3) and the beat log's `marked_*`/`post_*` columns (§7).

```
label               code  suppressesDetection  postEligible
1) R Peak Noise      1    true
2) Minor Noise       2                          (kept in threshold stats)
3) Blank.+Thresh.    3    (param edit: drag sets blanking/threshold)
4) PVC               4*   -                     true
5) PAC               5*   -                     true
6) Cond. Delay       4                          -
7) AF                5                           true
8) SVT               6                           true
9) VT                7                           true
Benign Arr.          10
Sig. Arr.            11
Other                12
Invert/Noninvert     13   (invert edit: drag flips inversion)
```

`*` PVC/PAC set to 4/5 per the latest edit — note this collides with Cond.Delay(4)
and AF(5); renumber those if the export needs all codes distinct.

---

## 9. Fast/slow split & the two template passes (`post_process.hpp`)

```
prepareViewerJob(cfg, binFs):                 # FAST — everything to open the viewer
    reset alignment::g_q_align = false        # each subject starts R-aligned only
    anneal if stale -> <stem>annealed.bin
    if wave+templates fresh: return job(canonical, needsFinalize=false)
    peakResults = raw R-peaks (create_ecg_ppg_pairs_raw) OR load from disk
    fast        = buildTemplatesAndBeatsFast(peakResults, rates)   # raw/unfilt/PPG only
    write_template_binfile(<stem>_templates.partial.bin, fast.tmpl)
    return job(provisional, needsFinalize=true)

finalizeViewerJob(job):                       # SLOW — background thread, file writes only
    if needSqabsDetection: write_output_binfile + write_output_csvfile   # §5
    write_template_binfile(<stem>_templates.bin, job.tmpl)              # §6a canonical
    remove(<stem>_templates.partial.bin)
    write_snips_csv(<stem>_template_snips.csv, job.beats)              # §6b
    # NOTE: squared/absval detection + slow merge are skipped (abs/square ignored).

regenerateWithQAlign(job):                    # on the viewer's first "Finish and Next"
    alignment::g_q_align = true
    fast = buildTemplatesAndBeatsFast(job.peakResults, rates)   # now R- then Q-aligned
    write_template_binfile(<stem>_templates.qalign.partial.bin, fast.tmpl)
    job.viewerTemplatePath = that path         # viewer reopens on it
```

Template alignment itself (R-align, optional Q-align, column-wise **median** across beats) lives in the template-generation headers
(`build_templates` → `CreateEcgTemplates` → `alignment`), which are separate from
this backend snapshot.

---

# 10. Mathematical Algorithms

Notation: `x[n]` input sample; `fs` sample rate (Hz); all indices 0-based unless a `.bin` field is noted 1-based.

## 10.1 Filtering (`bandpass.hpp`)

**Butterworth biquad design (bilinear transform).** For cutoff `fc`, pre-warp `wc = tan(π·fc/fs)`. Pole angles `θk = π(2k+order+1)/(2·order)`. Each conjugate
pole pair becomes one second-order section; with `A=1, B=−2cos(θ)·wc, C=wc²` and `a0=A+B+C`:

```
lowpass  section:  b = [C, 2C, C]/a0
highpass section:  b = [A, −2A, A]/a0
both:              a = [1, (−2A+2C)/a0, (A−B+C)/a0]
```

Odd order adds one first-order section.

**Biquad application** — Direct Form II Transposed:

```
y[n]  = b0·x[n] + z1
z1    = b1·x[n] − a1·y[n] + z2
z2    = b2·x[n] − a2·y[n]
```

**`filtfilt` (zero-phase).** Reflect-pad both ends by `pad = 3·(#sections)` using `2·x[0] − x[pad−i]` (and mirror at the tail), filter forward through the SOS
cascade, reverse, filter again, reverse, then strip the padding. Net phase = 0,
magnitude squared.

**`bandpass_filtfilt(order, lo, hi, fs, x)`** = `filtfilt(LP_hi, filtfilt(HP_lo, x))`.

## 10.2 Derivatives & smoothing

**`diff2` (smoothed derivative).** Weighted blend of four first differences:

```
d[i] = (2·(x[i+1]−x[i]) + 2·(x[i]−x[i−1]) + (x[i−1]−x[i−2]) + (x[i+2]−x[i+1])) / 6
```

Savitzky–Golay-flavored; apply `nd` times for higher orders.

**`nanfastsmooth` (NaN-aware moving average).** Running sum `s[k]` and running
valid-count `np[k]` over a width-`w` window (even windows half-weight the two
boundary samples to stay centered); output `s[k]/np[k]`, set to NaN where `np[k] < w·(1−tol)`. `type∈{1,2,3}` applies 1/2/3 passes (rectangular →
triangular → ~Gaussian).

## 10.3 Outlier rejection

**`stdoutlier`.** On the first difference `d = diff(x)`, moving mean `μ[i] = movmean(d, W)` and global `σ = std(d)`. Flag `d[i]` if it leaves `[μ[i] − k·σ, μ[i] + k·σ]` (direction = lower/upper/both); both endpoints of a
flagged difference are marked.

**Tukey (template beat rejection, in the alignment layer).** Keep values inside `[Q1 − 1.5·IQR, Q3 + 1.5·IQR]`; applied to beat *length* and to R *amplitude* (and to PPG fiducial position at 3.0·IQR).

## 10.4 Generic peak finder (`PeakFinder.hpp` / `findpeaks`)

1. Collect local maxima (plateau → center index).
2. Sort tallest-first (ties → earlier index) — matches MATLAB `findpeaks`.
3. Greedy min-distance elimination: accept a candidate only if no already-accepted
   peak lies at distance `< minPeakDistance` (O(N log N) via a sorted set of
   accepted positions).
4. Return survivors in temporal order.

## 10.5 Single-detector R-peak paths

**`rpeakdetect` (Hamilton-style).** Prep (threshold-independent):
detrend → `bandpass_filtfilt(2, 0.05, 100 Hz)` → `diff` → square →
moving-sum integrate (kernel length `≈7·fs/256`) → `medfilt1(·,10)` → remove
filter delay `⌈win/2⌉`. Reference `max_h` = max of the integrated signal over the
middle 50% of the record. Apply (per threshold `t`): mark regions where
integrated `> t·max_h`; within each region the R-peak is the max of the *band-passed* signal. `JoinedRR` runs this at `t = 0.2, 0.1, 0.4`.

**`pan_tompkin` (Pan–Tompkins).** Bandpass 5–15 Hz (dedicated coeffs at 256 Hz,
else `butter(3,·)`), normalize by max|·|; derivative; square; moving-window
integrate over `0.150 s`. Dual **adaptive thresholds** updated as `THR ← α·peak + (1−α)·THR` with signal/noise running estimates
(`ADAPT_FAST=0.125`, `SEARCHBACK` reweights 0.25/0.75); refractory `0.200 s`; **T-wave discrimination** at RR `< 0.360 s` by comparing max slope in a `0.075 s` window (ratio `0.5`); missed-beat **search-back** when a gap exceeds `1.66·RR_average`, with RR bounds `[0.92, 1.16]·RR_average`.

**`ecgLms`.** Baseline removal (subtract mean); moving-window integration width `175 ms`; threshold fraction `0.25` of an adaptively-forgotten peak
(`ff = 0.80`); refractory `250 ms`. (Called from the ensemble with a scalar `filtfilt(5,12,·)` = constant gain `(5/12)²`.)

**`RRsimpleSquared`.** Square the signal; threshold `mean + 2·std`; `findpeaks` with `minDist`; keep peaks above threshold.

## 10.6 Ensemble R-peak detection (`JoinedRR`) — the production detector

```
detrend(ecg)
run 6 detectors, weights w = [0.75, 0.25, 0.25, 1.25, 1.5, 0.75]:
   0..2  rpeakdetect @ thresh 0.2 / 0.1 / 0.4   (shared prep)
   3     pan_tompkin
   4     ecgLms
   5     RRsimpleSquared(minDist = median_of_detector_median_RR / 2)
refine detectors 3..5: snap each detection to the local max within
                       ± median(RR)/6 samples  (RPeakfromRWave)
pool all detections as (position, weight)
merge neighbours within diff_range = 8 samples -> keep the taller sample,
      reassign the loser's weight to the winner
sum weights per surviving position
ACCEPT positions whose summed weight >= 2.4      # i.e. agreement across detectors
```

So a peak is kept when enough independently-weighted detectors agree on it.

## 10.7 PPG fiducials

**Foot / onset (`find_foot_pulseox`, intersecting-tangent).** Shift the pulse so `d = row − max(row)`; the foot is located by the max-derivative "rotation"
(intersecting-tangent) method — the point maximizing the perpendicular distance
from the chord after rotating the pulse onto its rising tangent. Robust to
baseline offset (only relative shape matters).

**Segmentation (`SegmentPPG`).** Smooth (`nanfastsmooth`); build an
above/below-baseline mask; `RunLength` gives contiguous runs; each `0`-run
contributes its **min** (valley), each `1`-run its **max** (peak). Valleys flagged
as outliers are recomputed as the min between the two adjacent non-outlier peaks.

## 10.8 Template construction (`build_templates` / `CreateEcgTemplates` + alignment layer)

1. **Slice** each beat as `[R_i − 0.3·RR_i, R_i + 1.5·RR_i]` (fractions `percent_interval_preceeding/following_rpeak`).
2. **Reject** with Tukey (§10.3) on length, then on R amplitude.
3. **Representative length** = **median** RR (middle of the sorted lengths — a
   length some beat actually has), used to pick the DC-alignment reference beat
   and window sizes.
4. **Pass 1 — R-align:** place every beat on a shared axis with R at the common
   column `R_anchor` (NaN-pad the rest).
5. **Pass 3 — PR-baseline DC align:** subtract each beat's PR-segment mean so
   baselines match the reference beat.
6. **Optional Pass 4 — Q-align** (2nd template pass): per beat find the R-upstroke
   onset (steepest rising slope walked back to 20% of its max), then Q = first
   local minimum within ~20 ms before that onset (fallback = 10% of median RR
   before the onset); shift each beat so all Qs land on a common column.
7. **Average = column-wise NaN-skipping MEDIAN** across aligned beats: `T[c] = median{ beat[c] : not NaN }`.
8. **Per-sample std** (optional, gray band): sample std with `ddof = 1`, `σ[c] = sqrt( Σ(beat[c]−mean)² / (n−1) )`, NaN-skipped.

## 10.9 Annealing (`anneal_handler`)

Split the recording into fixed-length bins; **excise** noisy regions (an ECG span
is excluded only when *every* ECG channel the user marked agrees it is noisy —
un-marked channels don't constrain; PPG is handled independently); then **redistribute** leftover clean fragments into neighbouring bins where they fit,
so each emitted bin is ~one bin-length of clean signal. All 40 pass-through
channels are preserved into the annealed `.bin` (§4).

## 10.10 Constants quick-reference

| quantity                        | value                          | where              |
| ------------------------------- | ------------------------------ | ------------------ |
| beat window before/after R      | 0.3 / 1.5 × RR                 | alignment          |
| ensemble detector weights       | 0.75,0.25,0.25,1.25,1.5,0.75   | JoinedRR           |
| ensemble accept threshold       | ≥ 2.4                          | JoinedRR           |
| ensemble merge distance         | 8 samples                      | JoinedRR           |
| rpeakdetect thresholds          | 0.2 / 0.1 / 0.4 × max_h        | JoinedRR           |
| Pan–Tompkins band               | 5–15 Hz                        | pan_tompkin        |
| Pan–Tompkins integration        | 150 ms                         | pan_tompkin        |
| Pan–Tompkins refractory         | 200 ms                         | pan_tompkin        |
| T-wave refractory / slope win   | 360 ms / 75 ms                 | pan_tompkin        |
| ecgLms integration / refractory | 175 ms / 250 ms                | ecgLms             |
| RRsimpleSquared threshold       | mean + 2·std                   | RRsimpleSquared    |
| Tukey fence                     | 1.5·IQR (3.0 for PPG position) | alignment          |
| template average                | column-wise median             | CreateEcgTemplates |

> The alignment layer (`alignment.hpp`, R/Q-align + Tukey) and the template-marking
> GUI (viewer, glyph auto-detectors) are **not in this project snapshot**; §10.8's
> alignment steps are described from the template-generation contract. The
> column-wise median average itself lives in `CreateEcgTemplates.hpp` (shown in §10.8).
