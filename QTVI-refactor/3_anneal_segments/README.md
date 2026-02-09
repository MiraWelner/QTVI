## Step 1: Config and Environment Setup

First, parameters are defined

1.  `targetLength` = the desired length of each segment - this is currently hardcoded to 15 minutes

2. `min_size_minutes` = the absolute minimum length a segment can be to be considered valid

3. `expansion_seconds` = how much buffer to add around a noisy point to ensure the whole artifact is removed



### Step 2: Data Loading and Alignment

The code loads the raw data, which can include:

1) ECG

2) PPG

3) Sleep state

These are not super likely to all have the same sampling rate, thus the one with the lowest sampling rate is upsampled.



### Step 3: Noise Detection

First, it reads the start/stop times of known noise segments that were mark in the previous state

If no manual file exists, it runs a custom detection algorithm

1. First, it calculates the difference between the maxium and minimum values in a sliding window

2. It uses the Generalized Extreme Studentized Deviate to find spikes or flat-line dropouts ihn the PPG signal

3. It uses a 1D version of the DBSCAN clustering algorityhm. If a single sample is noisy, it marks a 15 seconds before and after as "noise" to ensure the segment is  fully excluded (in the original file this is set in the config.txt file)

4. 
