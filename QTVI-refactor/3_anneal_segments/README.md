## Step 1: Initialization and batch setup

1. First, the script calls `readProps.m` to parse a settings file. It identifies where the raw data is and where the results should go

2. File discovery: it performs a recursive serach across the folder its in for files ending in `*annealedSegments.mat` 

3. Iterates through every file found, for every file found it triggers main processing engine `FindWaveBounds.m`

### Step 2: Data Preperation

1. Loading the Segment: the script loads the `.mat` file, which contains ecg and ppg. Right now it does NOT contain sleep states because they do not work.

2. It cleans up the signal, removing global means and handling `NaN` to ensure the algorithm doesn't crash.



## Step 3: The ECG "Ensemble" Detection

This takes place in `JoinedRR.m`

1. Using `paarfor` it runs 5 different detection attempts simultaniously
   
   1. The `pan_tompkin` algorithm is run in `pan_tompkin.m`
   
   2. The `rpeakdetect.m` algorithm run 3 times, once with standard threshold, once with low threshold (0.1) and once with high threshold (0.4)
   
   3. The least mean squares appreach implemented in `ecgLms.m`

2. Peak fusion: This collects all detected time stamps and uses weighted voting to determine where the peaks are. 
   
   1. If it is found by `ecgLms` ig gets 1.5 points
   
   2. If found by `pan_tompkin` it gets 1.25 points
   
   3. If found by otheres it gets 0.25 to 0.75

3. Only peaks that accumulate a total score $\ge 2.4$   

4. Then in `RPeakfromRWave.m` the ensemble gives a general area of a beat. This function looks at the raw signal within a tiny window to find the absolute maxium point - the true R-Peak

## Phase 4: PPG Pulse Segmentation

This takes place in `SegmentPPG.m` and analyzes the PPG simultaniously to the ECG

1. The PPG is smoothed using `nanfastsmooth` with a 0.25 window to remove dicrotic notches and high-frequency noise

2. Baseline tracking: A 1 second moving average `movmean` is calculated to find the center of the wave

3. Valley detection
   
   1. The code looks for sections where the signal stays below the moving average
   
   2. It uses `RunLength.m` to identify these low zones
   
   3. THe lowest point in each zone is marked as a PPG Valley which is the start of a pulse

4. `stdoutlier.m` is used to calculate time between valleys - if a valley appears to soon or late then it is flagged and removed or corrected



## Phase 5: Signal Synchronization

This takes place in `pairRtoPPGBeat.m` Now the code has 2 lists: a list of ECG Heartbeats and a list of PPG Pulses. It must link them

1. Time Windowing: For every ECG R-peak detected in Phase 3, the code looks forward in time.

2. Matching: It searches for a PPG Valley that occurs shortly after R-peak (since the blok takes time to travel from the heart to the sensor)

3. Conflict Resolution:
   
   1. If it finds two PPG valleys for one R-peak, it chooses the one with the smallest timing error.
   
   2. If it finds an R-peak with no pulse, it marks it as a dropped beat
   
   3. If it finds a ppg pulse with no r-peak ti labels it as an unpaired pulse



## Phase 6: Result Assembly and Storage

1. Data Packaging: the code creates a large structure containing
   
   1. The exact sample indices of every R-peak
   
   2. The exact sample indices of every PPG Valley and Peak
   
   3. The "Pairs" matrix showing which ECG beat belongs to which PPG pulse

2. The final results are written as a .mat file in the output directory

3. Memory cleared, next file loaded
