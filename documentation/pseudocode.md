# Cardiac Signal Electronic Record Processing Script Pseudocode

Mira Welner

---

## 0. Standalone File to .bin script

This code is a standalone script that is used to convert a dataset into a uniform .bin file.

The three datasets from which it can read are:

## <span style="color: orange;">MESA</span>

| Source | File Type | Channel      | Original Sampling Frequency | Data Precision | File Length         | Upsampled to |
| ------ | --------- | ------------ | --------------------------- | -------------- | ------------------- | ------------ |
| MESA   | .edf      | EKG          | 256 Hz                      | 16-bit int     | 8–12 hours (approx) | 1000 Hz      |
|        |           | EOG_L        | 256 Hz                      | 16-bit int     |                     | 1000 Hz      |
|        |           | EOG_R        | 256 Hz                      | 16-bit int     |                     | 1000 Hz      |
|        |           | EMG          | 265 Hz                      | 16-bit int     |                     | 1000 Hz      |
|        |           | EEG1         | 256 Hz                      | 16-bit int     |                     | 1000 Hz      |
|        |           | EEG2         | 256 Hz                      | 16-bit int     |                     | 1000 Hz      |
|        |           | EEG3         | 256 Hz                      | 16-bit int     |                     | 1000 Hz      |
|        |           | Pres         | 32 Hz                       | 16-bit int     |                     | 32 Hz        |
|        |           | Flow         | 32 Hz                       | 16-bit int     |                     | 32 Hz        |
|        |           | Snore        | 32 Hz                       | 16-bit int     |                     | 32 Hz        |
|        |           | Thor         | 32 Hz                       | 16-bit int     |                     | 32 Hz        |
|        |           | Abdo         | 32 Hz                       | 16-bit int     |                     | 32 Hz        |
|        |           | Leg          | 32 Hz                       | 16-bit int     |                     | 32 Hz        |
|        |           | Aux_Ac       | 32 Hz                       | 16-bit int     |                     | 32 Hz        |
|        |           | Therm        | 32 Hz                       | 16-bit int     |                     | 32 Hz        |
|        |           | Pos          | 32 Hz                       | 16-bit int     |                     | 32 Hz        |
|        |           | EKG_Off      | 1 Hz                        | 16-bit int     |                     |              |
|        |           | EOG-L_Off    | 1 Hz                        | 16-bit int     |                     |              |
|        |           | EOG-R_Off    | 1 Hz                        | 16-bit int     |                     |              |
|        |           | EMG_Off      | 1 Hz                        | 16-bit int     |                     |              |
|        |           | EEG1_Off     | 1 Hz                        | 16-bit int     |                     |              |
|        |           | EEG2_Off     | 1 Hz                        | 16-bit int     |                     |              |
|        |           | EEG3_Off     | 1 Hz                        | 16-bit int     |                     |              |
|        |           | Pleth        | 256 Hz                      | 16-bit int     |                     | 1000 Hz      |
|        |           | OxStatus     | 1 Hz                        | 16-bit int     |                     | 1 Hz         |
|        |           | SpO2         | 1 Hz                        | 16-bit int     |                     | 1 Hz         |
|        |           | HR           | 1 Hz                        | 16-bit int     |                     | 1 Hz         |
|        |           | DHR          | 256 Hz                      | 16-bit int     |                     | 1000 Hz      |
|        | .xml      | Sleep Stages | 1/30 Hz                     | ASCII (XML)    |                     | 1/30Hz       |

Notes:

- The _off channels are not included in the final .bin file because the GUI will automatically not display absent data

- I am unclear why there are 3 EEG channels but that definitely seems to be the case looking at the MESA data

- MESA data is stored in a folder, so it is automatically read from the folder

## <span style="color: blue;">Bittium</span>

| Source  | File Type | Channel         | Sampling Frequency | Data Precision | File Length       | Upsampled to |
| ------- | --------- | --------------- | ------------------ | -------------- | ----------------- | ------------ |
| Bittium | .edf      | ECG 1           | 500 Hz             | 16-bit int     | 5–7 days (approx) | 1000 Hz      |
|         |           | ECG 2           | 500 Hz             | 16-bit int     |                   |              |
|         |           | ECG 3           | 500 Hz             | 16-bit int     |                   |              |
|         |           | Accelerometer X | 25 Hz              | 16-bit int     |                   |              |
|         |           | Accelerometer Y | 25 Hz              | 16-bit int     |                   |              |
|         |           | Accelerometer Z | 25 Hz              | 16-bit int     |                   |              |
|         |           | Marker          | 1 Hz               | 16-bit int     |                   |              |
|         |           | DEV_Temperature | 1 Hz               | 16-bit int     |                   |              |
|         |           | Pacemaker Event | 8 Hz               | 16-bit int     |                   |              |

## <span style="color: green;">CHAOS</span>

| Source | File Type | Channel Name                | Channel Type | Sampling Frequency | Data Precision           | File Length | Upsampled to |
| ------ | --------- | --------------------------- | ------------ | ------------------ | ------------------------ | ----------- | ------------ |
| Chaos  | .dat      | Timestamp                   | Timestamp    | 500 Hz             |                          | 15 minutes  |              |
|        |           | NLS_NOM_ECG_ELEC_POTL_I     | ECG I        | 500 Hz             | ASCII loaded into double |             | 1000 Hz      |
|        |           | NLS_NOM_ECG_ELEC_POTL_II    | ECG II       | 500 Hz             | ASCII loaded into double |             | 1000 Hz      |
|        |           | NLS_NOM_ECG_ELEC_POTL_III   | ECG III      | 500 Hz             | ASCII loaded into double |             | 1000 Hz      |
|        |           | NLS_NOM_PULS_OXIM_PLETH     | PPG          | 125 Hz             | ASCII loaded into double |             | 1000 Hz      |
|        |           | NLS_EEG_NAMES_EEG_CHAN1_LBL | EEG I        |                    |                          |             |              |
|        |           | NLS_EEG_NAMES_EEG_CHAN2_LBL | EEG II       |                    |                          |             |              |
|        |           | NLS_EEG_NAMES_EEG_CHAN3_LBL | EEG III      |                    |                          |             |              |
|        |           | NLS_EEG_NAMES_EEG_CHAN4_LBL | EEG IIII     |                    |                          |             |              |
|        |           | NLS_NOM_PRESS_BLD_VEN_CENT  | CVP          | 125 Hz             | ASCII loaded into double |             | 1000 Hz      |
|        |           | NLS_NOM_PRESS_BLD_ART_ABP   | aBP          | 125 Hz             | ASCII loaded into double |             | 1000 Hz      |
|        |           | NLS_NOM_RESP                | Resp         | 62.5 Hz            | ASCII loaded into double |             | 500 Hz       |
|        |           | NLS_NOM_PRESS_BLD_ART       | Art          | 125 Hz             | ASCII loaded into double |             | 1000 Hz      |
|        |           | NLS_NOM_PRESS_BLD_ART_PULM  | Art Pulm     | 125 Hz             | ASCII loaded into double |             | 1000 Hz      |
|        |           | CSUM                        | Checksum     | 500Hz              |                          |             |              |

Notes:

- While there exist 4 EEG channels, in none of the CHAOS files could I find any actual examples of files that had EEG data

- Not every CHAOS file has all these channels, some are missing a few

- The timestamp used to be used to detect drops, but some CHAOS files have an error where the Hz is 2x for n seconds and then dropped for n seconds, so the timestamp is unreliable.

## Upsampling Algorithm

There are two upsampling algorithms used. One is for the MESA and Bittium files, which are large and thus need a faster algorithm:

- **Prototype low-pass filter (windowed sinc).**
  - Cutoff: `fc = 1 / max(P, Q)` (normalized frequency, i.e. cutoff at the more restrictive of the input Nyquist or the output Nyquist).
  - Length: `2 · halfLobes · max(P, Q) + 1` taps. Filter grows with the resample ratio so longer filters are used when the ratio is aggressive.
  - Impulse response: `h[n] = sinc(fc · (n − M/2)) · w[n]`, where `sinc` is the ideal low-pass impulse response and `w[n]` is a **Blackman window** (`0.42 − 0.5·cos + 0.08·cos`) that suppresses ringing.
  - Symmetric around the center, so it's a linear-phase FIR — no frequency-dependent time delay in the passband.
- **Polyphase decomposition.**
  - The single prototype filter of length `numTaps` is split into `P` sub-filters ("polyphase banks") by indexing: sub-filter `p` contains taps `h[p], h[p+P], h[p+2P], ...`.
  - Each sub-filter has length `subLen = ⌈numTaps / P⌉`.
  - Each sub-filter is then **DC-normalized** — its taps are divided by their sum so each sub-filter has gain 1 at DC. Preserves signal amplitude at low frequencies.
- **Output generation.** For each output sample at time-index `m`:
  - Compute the corresponding real-valued input position `m · Q / P`, split into an integer input base `baseInput` and a sub-filter selector `p = (m · Q) mod P`
  - Dot-product sub-filter `p` against a `subLen`-sample window of the input around `baseInput`. That single dot product is the output sample.
  - Because polyphase reindexes the math this way, the filter's zeros in the zero-stuffed intermediate signal are never actually multiplied or stored — each output sample costs `subLen` multiply-adds instead of `numTaps` multiply-adds against a mostly-zero buffer.
- **Group-delay compensation.** A `filterCenter` offset (`halfLobes · max(P, Q) / P`) shifts the output-to-input alignment so that a peak in the input lands at the correct time in the output. Prevents a systematic time-shift between channels resampled at different `P/Q` ratios.
- **Three-region processing.** The output run is split by where the filter window fits relative to the input's edges:
  - **Left boundary:** filter partially overhangs the start of the input. Missing samples treated as zero (zero-padding).
  - **Middle (main):** filter fits entirely inside the input. Fast path, no bounds checks per sample. Parallelized across CPU cores.
  - **Right boundary:** filter partially overhangs the end of the input. Same zero-padding as the left boundary.
- **Per-channel independence.** Every channel calls this function with its own `(P, Q)` — no shared clock is assumed between channels. This is what lets the pipeline handle datasets like MESA where signals in the same recording have wildly different native rates (256 Hz, 32 Hz, 1 Hz) and each ends up on its own upsampled grid.

CHAOS files are shorter so they have a much simpler algorithm

- **Row-grid setup.** The `.dat` has a single high-rate row clock (`row_rate`) shared across all columns. Each channel's real samples are extracted along with the **row index where each real sample sits**, producing a pair of arrays: `rawValues[k]` and `rawRowIdx[k]`. Time of the k-th real sample is `rawRowIdx[k] / row_rate`.
- **Double-rate-region correction (pre-resample).** Some stretches of the `.dat` are recorded at twice the nominal rate as a packing artifact, with the extra rows placed adjacent rather than at their true times. `fix_double_rate_regions` walks the row indices and rewrites them so:
  - Doubled-up rows get spread back onto the correct time grid.
  - Genuine gaps (missing/dropped rows) are preserved.
  - The upsampled and raw blocks that follow share one consistent clock.
- **Resampling itself: linear interpolation between real samples.** For each output sample at target-rate index `m` in `resample_from_sparse`:
  - Compute output time `t = m / targetRate`.
  - Find the pair of consecutive real samples that bracket `t`: `rawValues[k]` at time `t0 = rawRowIdx[k] / row_rate`, and `rawValues[k+1]` at time `t1 = rawRowIdx[k+1] / row_rate`.
  - Compute fractional position `f = (t − t0) / (t1 − t0)`.
  - Output value = `rawValues[k] · (1 − f) + rawValues[k+1] · f`.
  - A running index `k` advances monotonically through the input as `m` advances, so the whole resample is one linear pass, not a per-output search.
- **Edge handling.**
  - Before the first real sample: output holds the first real sample's value.
  - After the last real sample: output holds the last real sample's value.
  - Between real samples: linear interpolation (as above).
  - So gaps get filled with a straight line between the nearest observed values, not with zeros — appropriate for the physiological signals this format carries (a heart-rate value doesn't drop to zero just because a row is blank).
- **Raw block.** In parallel with the upsampled block, each real sample is written to the `.bin`'s raw section as a `(time, value)` pair, using the same row-index-derived time (`rawRowIdx[k] / row_rate`, post-correction). The `.dat`'s own wall-clock timestamps are unreliable and are **not** used for either block — row index is the single source of truth after `fix_double_rate_regions`.
- **Timestamp channel is still synthetic.** Channel 0 in the output `.bin` is `k / rate`, same as every other dataset — the CHAOS-specific `parse_dat_timestamps` is only used internally for detecting doubled regions, never written to the output.

## Output .bin file

The structure of a bin file which will be output by the standalone code to be read by the GUI is:

Each channels holds two arrays: 

**Raw Data** 

- Stored as `(x, y)` pairs. Each element is two float64s.
- `x` = time in seconds from the start of the recording.
- `y` = the sample value in physical units (mV for ECG, etc.).
- Length is `sizes_raw[i]` pairs.

**Upsampled array**

- Filtered and resampled to a uniform target rate (`up_rates[i]` Hz).
- Stored as a plain sequence of float64 values — no `x`, because timing is implicit: sample `k` is at time `k / up_rates[i]` seconds.
- Length is `sizes_up[i]` samples

**Header (bytes 0..583)**

| Byte range | Type         | Field              | Contents                                                         |
| ---------- | ------------ | ------------------ | ---------------------------------------------------------------- |
| 0..3       | uint32       | `sleep_state_len`  | Sleep-stage epoch length in seconds (30 for MESA, 0 for others)  |
| 4..147     | 36 × uint32  | `sizes_up[0..35]`  | Element count of each channel's upsampled array                  |
| 148..291   | 36 × uint32  | `sizes_raw[0..35]` | Pair count of each channel's raw array                           |
| 292..435   | 36 × float32 | `native_rates`     | Each channel's source rate in Hz. 0.0 = channel absent           |
| 436..579   | 36 × float32 | `up_rates`         | Each channel's upsampled target rate in Hz. 0.0 = channel absent |
| 580..583   | uint32       | `size_sleep`       | Sleep-stage epoch count appended at end of file (0 if not MESA)  |

**Body**

| Idx | Channel         | MESA up rate | Bittium up rate | CHAOS up rate |
| --- | --------------- | ------------ | --------------- | ------------- |
| 0   | Timestamp       | 1000 Hz      | 1000 Hz         | 1000 Hz       |
| 1   | ECG 1           | 1000 Hz      | 1000 Hz         | 1000 Hz       |
| 2   | ECG 2           |              | 1000 Hz         | 1000 Hz       |
| 3   | ECG 3           |              | 1000 Hz         | 1000 Hz       |
| 4   | PPG             | 1000 Hz      | 1000 Hz         | 1000 Hz       |
| 5   | Accelerometer X |              | 200 Hz          |               |
| 6   | Accelerometer Y |              | 200 Hz          |               |
| 7   | Accelerometer Z |              | 200 Hz          |               |
| 8   | Marker          |              | 1 Hz            |               |
| 9   | Temperature     |              | 1 Hz            |               |
| 10  | Pacemaker Event |              | 8 Hz            |               |
| 11  | EOG-L           | 1000 Hz      |                 |               |
| 12  | EOG-R           | 1000 Hz      |                 |               |
| 13  | EMG             | 1000 Hz      |                 |               |
| 14  | EEG 1           | 1000 Hz      |                 |               |
| 15  | EEG 2           | 1000 Hz      |                 |               |
| 16  | EEG 3           | 1000 Hz      |                 |               |
| 17  | EEG 4           |              |                 |               |
| 18  | CVP             |              |                 | 1000 Hz       |
| 19  | Pres            | 32 Hz        |                 |               |
| 20  | Flow            | 32 Hz        |                 |               |
| 21  | Snore           | 32 Hz        |                 |               |
| 22  | Thor            | 32 Hz        |                 |               |
| 23  | Abdo            | 32 Hz        |                 |               |
| 24  | Leg             | 32 Hz        |                 |               |
| 25  | Aux_AC          | 32 Hz        |                 |               |
| 26  | Therm           | 32 Hz        |                 |               |
| 27  | Pos             | 32 Hz        |                 |               |
| 28  | OxStatus        | 1 Hz         |                 |               |
| 29  | SpO2            | 1 Hz         |                 |               |
| 30  | HR              | 1 Hz         |                 |               |
| 31  | DHR             | 1000 Hz      |                 |               |
| 32  | Resp            |              |                 | 500 Hz        |
| 33  | aBP             |              |                 | 1000 Hz       |
| 34  | ART             |              |                 | 1000 Hz       |
| 35  | ART_PULM        |              |                 | 1000 Hz       |

**<span style="color: red;">TODO: split .bin files such that they are max 8 hours and longer bins are separated into multiple bins of max 8 hours</span>**

---

## 1. Config and Dataset collection

1. A file called `config.csv` in the same directory as the `.exe` file is loaded. It has 3 rows corresponding to MESA, Bittium, and CHAOS. The columns are:

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

**<span style="color: red;">TODO: Add to the config.csv an 'hours loaded at time' section - right now default 8 hours of the .bin are loaded into the gui at one time but if somebody has a small computer this can be hard so it should be configurable.</span>**

---

## 1. Noise markings

The first step for the user is to mark the noise and other phenomena on the noise marking GUI. The GUI is formatted as follows depending on dataset:

### GUI design

**MESA**

```
|-------------------------------------------------------|
| controls | [ ECG ampogram ---------------------------]|
|          | [ PPG ampogram ---------------------------]|
|          | [ Sleep states ---------------------------]|
|-------------------------------------------------------|
|           ECG1(EKG) - markable                        |
|           PPG(Pleth) - markable                       |
|-------------------------------------------------------|
```

**Bittium**

```
|---------------------------------------------------------|
| controls | [ ECG ampogram -----------------------------]|
|          | [ Temperature ------------------------------]|
|          | [ Marker----] [ Resp-------] [ Pacer -------]|
|---------------------------------------------------------|
|                 ECG1 -  markable                        |
|                 ECG2 -  markable                        |
|                 ECG3 -  markable                        |
|                 Accel - markable                        |
|---------------------------------------------------------|
```

**CHAOS**

```
|---------------------------------------------------------|
| controls | [ ECG ampogram -----------------------------]|
|          | [ PPG ampogram -----------------------------]|
|---------------------------------------------------------|
|                      ECG1 -- markable                   |
|                      ECG2 -- markable                   |
|                      ECG3 -- markable                   |
|                      PPG -- markable                    |
|                      ABP -- markable                    |
|                      ART -- markable                    |
|                   ART_PULM -- markable                  |
|---------------------------------------------------------|
```

### Potential Markings

The markings that the user can make are:

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

In the output bin file and the log file (which will be discussed later) each marking is represented by the number it corresponds to above. So for example, minor noise is marked with a 2.

### Peak Finding Algorithm

There is a faster and less accurate peak finding algorithm which is used for the GUI to mark the peaks - it is not what is eventually used to make the templates etc but it can be used to inform the user what should be marked as noise, pvc, etc, and the noise is excluded from the templates. The algorithms are:

### ECG R-peak detection

- **Scan for candidates.** For each sample index `i` in the visible range (visible range determined by user - you can select 1s, 3s, 10s, 30s, 1m, 2m 5m):
  - Let `t` be the time of the candidate peak, `s` be the sign (+1 upright -1 inverted) and `y` be the height of the candidate peak.
  - Only a candidate if the absolute value is greater than the prior sample and subsequent sample (local max), and it is not in a `no_peaks_in_region` region (such a region exists when marked as R peak noise by user)
  - Fetch reference statistics for this beat (see next section). Skip if not OK.
  - Remove candidate peaks below threshold
    - Threshold = `t(median of r peaks in 10 second region prior - median of signal in 10 second region prior)` where `t` is selected by user in config file.
    - User can change threshold region via markings.
  - Passing samples enter a preliminary list.
- **Reference statistics for a beat at time `t`**
  - The reference always comes from samples strictly earlier than `t` (causal windowing), so the detector can run on sliding chunks without needing the whole recording.
  - **Choose reference window.**
    - If `t` is inside a `withinSpan [aS, aE]`: use `[max(aS, t - W), t]` where `W = previous_seconds_to_train_on = 10 s`. No cache; the reference is bespoke to this annotated span. This is the "local reference" mode used for user-annotated regions — the detector adapts inside the annotation instead of reaching back to prior clean data.
    - Otherwise (frame mode): compute frame `m = floor(t / W)`, set `k = m - 1`, and walk backward past any frame that overlaps `do_not_learn_from_region`. Floor at 0. Reference window is `[k*W, (k+1)*W]`. Cache key `(k, sign)` — `RefStats` computed once per `(frame_index, sign)` is reused for every candidate in that frame.
  - Cache hit in frame mode returns immediately.
  - **Amplitude range:** collect sign-corrected amplitudes in the window, dropping samples inside `do_not_learn_from_region` (unless we're in a `withinSpan`, where `do_not_learn_from_region` is ignored inside). If exclusion left fewer than 2 samples, retry without exclusions. Return not-OK if `vMax_raw <= vMin_raw`.
  - **Reference peak collection:** scan the window for the same three-sample local-max pattern, gated by the local threshold formula `vMin_raw + threshold(t) * (vMax_raw - vMin_raw)`.
  - **Aggregate:**
    - `vMin` = median of sign-corrected amplitudes in the window (baseline).
    - `gateTop` = median amplitude of reference peaks, or `vMax_raw` if empty.
    - `meanRR` = mean interval between reference peaks; in frame mode, gaps inside `do_not_learn_from_region` are skipped. (Carried on `RefStats` for the cache key's frame; not consumed by the blanking pass below.)
  - Cache in frame mode and return.
- **Blanking pass.** Walk the preliminary list in time order. For each beat, compute `blank = blanking(t) / 1000`. If `blank > 0` and `t - lastT < blank`, drop the beat; otherwise keep it and set `lastT = t`.
- **Emit** only kept beats with `detStart <= t <= detEnd`.
- **How the two exclusion lists are used.**
  - `do_not_learn_from_region` protects the baseline — its spans are ignored when computing `vMin`, `gateTop`, and `meanRR` for the reference-window statistics.
  - `no_peaks_in_region` suppresses emissions — candidates whose time falls inside these spans are dropped from the scan.
  - The two are populated independently, so a region can be in one, both, or neither.
- **Time-varying parameters.** `threshold(t)`, `blanking(t)`, and `sgn(t)` are all functions of time, backed by `ParamIndex` for O(log n) piecewise-constant lookup. User-painted per-region overrides apply cleanly at boundaries.

### PPG / arterial peak detection

- **Scan for candidates.** For each sample index `k` in the scan range:
  - Compute three consecutive squared first-differences:
    - `d²[k-1] = (value[k] - value[k-1])²`
    - `d²[k] = (value[k+1] - value[k])²`
    - `d²[k+1] = (value[k+2] - value[k+1])²`
  - **Upstroke local-max test:** require `d²[k] >= d²[k-1]` and `d²[k] > d²[k+1]`. Skip if it fails.
  - **Apex walk:** from `tUp = time[k]`, march forward through subsequent samples for up to 0.2 s, tracking a running max `(apexIdx, apexVal)`. Stop as soon as the value drops. The reported peak time is `time[apexIdx]`, not `tUp` — the systolic peak arrives slightly after the fastest upstroke.
  - Skip if `time[apexIdx]` is inside `no_peaks_in_region`.
  - Fetch reference statistics for this beat. Skip if not OK.
  - **Upstroke-strength gate:** require `d²[k] >= upstrokeGate`.
  - **Amplitude gate:** require `apexVal >= vMin + threshold(t) * (gateTop - vMin)`.
  - Passing beats enter the preliminary list.
- **Reference statistics for a beat at time `t`**
  - The reference always comes from samples strictly earlier than `t` (causal windowing), so the detector can run on sliding chunks without needing the whole recording.
  - **Choose reference window.**
    - If `t` is inside a `withinSpan [aS, aE]`: use `[max(aS, t - W), t]`. No cache; local reference for user-annotated regions.
    - Otherwise (frame mode): frame `m = floor(t / W)`, walk `k` backward from `m - 1` past frames overlapping `do_not_learn_from_region`, floor at 0. Reference window `[k*W, (k+1)*W]`. Cache key = frame index — `RefStats` computed once per frame is reused for every candidate in that frame.
  - **Amplitude range:** as in ECG, but require at least 4 samples in the window; same exclusion-and-fallback logic. Return not-OK if `vMax_raw <= vMin_raw`.
  - **Squared-derivative distribution:** collect `d²` values across the window (respecting `do_not_learn_from_region` outside `withinSpans`, with fallback if empty). Sort and read the 90th-percentile value. Return not-OK if empty or if that value is `<= 0`.
    - `upstrokeGate = 0.6 * d²_at_90th_percentile`.
  - **Reference peak collection:** run `collectSystolic` over the reference window (same upstroke-max + apex-walk procedure, gated by the local threshold formula and the freshly-computed `upstrokeGate`). Then `cleanReferencePeaks` drops any that fall inside `do_not_learn_from_region`.
  - **Aggregate:**
    - `vMin` = median amplitude in the window.
    - `gateTop` = median amplitude of reference peaks, or `vMax_raw` if empty.
    - `meanRR` = mean interval between reference peaks. (Carried on `RefStatsD`; not consumed by the blanking pass below.)
  - Cache and return `{vMin, gateTop, meanRR, upstrokeGate}`.
- **Blanking pass.** Walk the preliminary list in time order. For each beat, compute `blank = blanking(t) / 1000`. If `blank > 0` and `t - lastT < blank`, drop the beat; otherwise keep it and set `lastT = t`.
- **Emit** only kept beats with `detStart <= t <= detEnd`.
- **How the two exclusion lists are used.**
  - `do_not_learn_from_region` protects the baseline — its spans are ignored when computing `vMin`, `gateTop`, `meanRR`, and the `d²` distribution that produces `upstrokeGate`.
  - `no_peaks_in_region` suppresses emissions — candidates whose time falls inside these spans are dropped from the scan.
  - The two are populated independently, so a region can be in one, both, or neither.
- **Time-varying parameters.** `threshold(t)` and `blanking(t)` are functions of time, backed by `ParamIndex` for O(log n) piecewise-constant lookup. User-painted per-region overrides apply cleanly at boundaries.

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

### Output Logfile

At the beginning of the program after choosing which dataset to analyze, the user is prompted for their initials. The log can be found in `log_[initials]/[ID]_log.csv`.

One purpose the logfile serves is to ensure that if a log has already been made of a file, when the GUI comes upon that file, it skips it over and prints that it skipped it, because it has already been analyzed by that particular person, who reported with their initials that they were the one who logged it.

- <span style="color: red;">**TODO:  right now the log file prevents a load even if the user only opened the file for a few seconds and did nothing. There should be a way to mark that the log is in progress and reload what has been logged without skipping over the file unless it is fully logged.**</span>**

Unlike the other binfiles and csvs, the log is only updated when a beat is viewed by the user. So if the user uses the ampogram to skip over many beats, they will simply not be logged. It works like this:

* The log file starts its life full of zeroes, to decrease the amount of memory used when writing to the log.

* The viewer automatically marks a part of the signal as an R peak or arterial peak

* `handle_data_plot()` is called, which logs all visible peaks in the log file.

* Visible peaks get added to a <map> object, keyed by global x location so that the subsequent times a peak is seen, it overwrites the previous wrote 

* ACCEL has no peaks so it is logged when ECG1 has a peak
  
  * **<span style="color: red;">TODO: acceleration data should be split into discrete bins and the log file should record when a state switches. However it doesn't easily split itself to discrete bins right now so further analysis must be done to figure out what the discrete states should be.</span>**

* Every 30 seconds, a `QTimer` calls `flushPending()`, which merges pending into the committed table, and for each channel with a non-empty pending map:
- For each channel with a non-empty pending map:
  
  - `ensureCapacity(merged.size())` grows `all_beats_in_log` in 5700-row chunks if the merged set is bigger than the current table.
  - **Every row's slot for that channel is zeroed out first** (`x=0, y=0, markType=0, postType=0, inverted=0`), then the merged map — which iterates in ascending time order — is written back in row order `0, 1, 2, …`.
  - The pending map for that channel is cleared.
  - Immediately after `flushPending()`, the timer callback also calls `writeCsv(...)`, rewriting the full CSV.

- It scans `all_beats_in_log`, skips any row where every channel's `x == 0.0` (fully empty), and writes the rest — so the file's row count tracks real beat count, not the vector's allocated capacity.

- Edits invalidate stale log entries via `removeInRange`.**
  
  - Whenever a marking is added/erased, or a threshold/blanking/invert override is applied, the code does a few things:
    - Erases any pending entries in `[t0, t1]` from the pending map.
    - Zeroes `x` (and the other fields) on any **committed** row whose beat time falls in `[t0, t1]` — this is the "empty slot" sentinel, so `writeCsv` will skip that row if nothing else is beating there.
    - No re-detection happens inside `removeInRange` itself — it only deletes. The next `handle_data_plot()` (which the caller always triggers right after) re-runs detection over the now-changed region and re-logs whatever peaks actually exist there now.

### Log file column reference

| Column                  | Type          | Meaning                                      | Source                                                                         | Notes                                                                                                                                         |
| ----------------------- | ------------- | -------------------------------------------- | ------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------- |
| `beat`                  | int           | Row index in the table                       | Row position in `all_beats_in_log`                                             | Not a stable ID — a beat's row can shift between flushes as the channel's merged map is rewritten in time order. Not aligned across channels. |
| `ecg1_x` … `art_pulm_x` | double (6 dp) | Global time of the beat (seconds)            | Detected peak's `x` + chunk offset                                             | `0.0` = empty slot (no beat in this row for this channel)                                                                                     |
| `ecg1_y` … `art_pulm_y` | double (6 dp) | Amplitude of the beat                        | Detected peak's `y`                                                            | Raw detector output, not gain-scaled                                                                                                          |
| `accel_x`               | double (6 dp) | Time of the corresponding ECG1 beat          | Anchored to ECG1's beat time, not a real ACCEL detection                       | ACCEL has no detector of its own                                                                                                              |
| `accel_y`               | double (6 dp) | Accel **X-axis** reading at that time        | `m_accelX[idx]` sampled at the ECG1 beat index                                 | Despite the name, this is not the accel Y-axis; Y/Z axes aren't logged at all                                                                 |
| `blanking_<chan>`       | double (6 dp) | Blanking value (ms) in effect for that beat  | Config default, or active per-region override                                  | One column per channel, incl. `blanking_accel`                                                                                                |
| `threshold_<chan>`      | double (6 dp) | Threshold value in effect for that beat      | Config default, or active per-region override                                  | One column per channel, incl. `threshold_accel`                                                                                               |
| `marked_<chan>`         | int           | Annotation-type code covering the beat       | `annotation_types` code, via `MarkSpanIndex::codeAt`                           | `0` = no annotation at this beat                                                                                                              |
| `marked_accel`          | int           | Always `0`                                   | Hardcoded placeholder                                                          | ACCEL isn't its own markable channel yet                                                                                                      |
| `post_<chan>`           | int           | Post-arrhythmia tag                          | Nonzero if this beat is the first one after an eligible AF/SVT/VT/PVC/PAC span | `0` = not a post-beat                                                                                                                         |
| `post_accel`            | int           | Always `0`                                   | Hardcoded placeholder                                                          | Concept doesn't apply to ACCEL                                                                                                                |
| `inverted_<chan>`       | int (0/1)     | Whether the channel is inverted at this beat | Lead-reversal checkbox XOR any active invert-override span                     | `1` = inverted, `0` = not                                                                                                                     |
| `inverted_accel`        | int           | Always `0`                                   | Hardcoded placeholder                                                          | ACCEL has no invert concept                                                                                                                   |

**Channel order used in every repeated-column group:** `ecg1, ecg2, ecg3, ppg, abp, art, art_pulm, accel`

**Column count:** 57 total = `beat` (1) + 8 channels × (`x, y, blanking, threshold, marked, post, inverted` = 7 fields each)

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
- <span style="color: red;">:**TODO For tracings that are too noisy that cannot be avoided, use 20% of the onset value instead of the minimum foot baseline to stay above the noise level.**</span>

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

### Normalization:

The ECG Normalization algorithm is as follows:

1. Calculate the total QRS vector magnitude for each beat to control for cardiac axis rotation:
    RS_peak(t) = abs(R_peak(t)) + abs(S_peak(t))
2. Find the global reference by taking the median across all bins for each individual:
    Global_Ref_person = median(RS_peak(t))
3. Normalize any amplitude feature like P, R, or T waves using the equation:
    Feature_peak_norm_abs = Feature_peak(t) / Global_Ref_person

The PPG normalization algorithm is as follows

1. First, calculate the local PI for each beat:
    PI(t) = ((systolic_peak(t) - diastolic_trough(t)) / abs(diastolic_trough(t))) * 100.
2. Find the global reference by taking the median PI across the entire recording for that individual:  Global_Ref_person = median(PI(t)).
3. Normalize your amplitude feature by converting it to its local baseline ratio first and then dividing by the global reference:  Feature_peak_norm_abs = Feature_Local_Ratio(t) / Global_Ref_person.



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
- **P50**: half way between onset and peak
- **Peak**: global max of the trace
- **End**: minimum after the peak
- **Dicrotic notch**: first interior local minimum between peak and end (falls back to 1/3 of the way from peak to end if none found)
- **Second Peak**: After the dicrotic notch there is sometimes a second peak before the foot/trough
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
| `file_id                                                                                | ID of subject                                                                               |
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

## CSV output — `<subject>_template.csv`

Per-sample trace CSV, written by `writeAlignedTemplateCsv`. One row per sample, all bins stacked. R-aligned and Q-aligned passes are merged into this single file with `_r`/`_q` suffixed columns (shared key columns appear once).

| Column pattern                                                                          | Meaning                                                                                     |
| --------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| `file_id, bin_num, x_ms`                                                                | Subject ID, bin number, sample time (shared, not suffixed)                                  |
| `<chan>_raw_mv_{r/q}`                                                                   | Raw amplitude at this sample, per channel (ch1, ch2, ch3, ppg, abp, art, art_pulm)          |
| `<chan>_Normalized_{r/q}`                                                               | Normalized amplitude at this sample, per channel                                            |
| `<chan>_raw_std_{r/q}`                                                                  | Per-sample std of raw amplitude, per channel                                                |
| `<chan>_normalized_std_{r/q}`                                                           | Per-sample std of normalized amplitude, per channel                                         |
| `<pt>_ch{1-3}_location_{autodetect/user}_{r/q}`                                         | 1 if this row's sample index is that marker's position, else blank (per ECG marker/channel) |
| `<marker>_location_{autodetect/user}_{r/q}` (ppg/abp/art/art_pulm)                      | Same one-hot marker flag for pulse-channel markers                                          |
| `p_wave_ch{1-3}, q_onset_ch{1-3}, r_wave_ch{1-3}, t_peak_ch{1-3}` (auto only, `_{r/q}`) | One-hot flags for the autodetected computed ECG glyphs                                      |
| `ppg foot, p1` (auto only, `_{r/q}`)                                                    | One-hot flags for the autodetected computed PPG glyphs                                      |
