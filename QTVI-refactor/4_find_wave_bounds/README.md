## Find Wave Bounds

This script takes in a bin file with the following format:

Row 1: ECG Signal

Row 2: PPG Signal

Row 3: ECG sample rate 

Row 4: PPG sample rate

Row 5: ECG Bin indexes

Row 6: PPG Bin indexes

### Phase 1:

1) Configuration loading - the original MATLAB read a config.txt, the new C++ will read a config.csv file which identifies the input directory containing the files to be loaded. Originally mat, now bin

2) Recursively seraches the input directory for all files matching the pattern *annealedSegments.mat

3) For each file found, it checks if a corresponding wave_data.mat already exists in the output folder, if so, skip, unless force_run is set

### Phase 2:

1. Loads the cell array from the .mat file

2. Iterates through every segment (cell) in the array. EAch segment typically represents a few minutes of ECG/PPG data

### Phrase 3:

1. Looks for existing `r_peaks` in the input data. If missing, triggers detection pipeline
   
   1. Run 5 different R-R peak detection algorithms in parallel
      
      1. Algorithms 1-3 are modified Pan Tompkins based algorithm with 3 different sensitivity thresholds (standard, 0.1, and 0.4)
      
      2. Algorithm 4 is a full implmentation of the Pan Tompkins algorithm
      
      3. Algorithm 5 uses least median squares adaptive filtering to find peaks
   
   2. Each algorithm is assigned a weight form 0.25 to 1.5, detected peaks from all algorithms are combined
   
   3. Peaks within 2 samples of eachother are merged to the location of the highest amplitude
   
   4. A peak is only 'accepted' if the sum of weights from the algorithm that it found are $\ge 2.4$ 
   
   5. It then validates it by checking if the number of detected R peaks are physically plausible compared to the number of PPG pulses found. If it fails, the segment is flagged as noise
   
   ## Phase 4: PPG pulse segmentation
   
   1. Smooths the PPG signal using a window of 25% of the sampling rate to remove high-frequency noise
   
   2. Baseline clipping: calculates a 1-second moving mean. It identifies 'valleys' by looking at signals ares that fall below the moving mean
   
   3. Finds the local minima and maxima, calculates the amplitude and period for all of them
   
   4. Identifies beats with durations and amplitudes that deviate more than 2.5 std from the mean
   
   5. For gaps created by outliers, research the raw signal between the surrounding valid peaks to find the most likely missing pulse valley

2. 
3. 
   
   1. 


