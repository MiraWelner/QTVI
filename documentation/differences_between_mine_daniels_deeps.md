My script was based off Daniel's QTVI code, because although it is less complete than Deep's code, I had an easier time reading it. However, after replicating Daniel's code, I changed it to be based off Deep's code. The core is still Daniel's code. This script describes the algorithmic differences between his code and mine.

This does not describe cosmetic differences in the GUI, just algorithmic differences.

## Filtering

* Neither Daniel's nor my code uses any filtering prior to the noise marking stage. The reason for lack of apparent baseline drift is likely due to the fact that the screen chooses min and max based on the min and max of the currently displayed window, which is small.

* I am unclear why Daniel's code filters noise better than mine.

## Templating:

- In Daniel's code, the templates are negative. You cannot see this because he doesn't display the y axis in his marking GUI. This is because of a bug (I think) where he shifts them, he accidentally shifts by the average peak, so the peaks are vertically aligned at -average peak. I removed this bug so they are aligned at the correct positive location (if average peak is positive - if for some reason it is negative then it'll be aligned at a negative y location)

- For throwing away data that does not align with the 'normal' morphology (in this case Daniel and I both used) Daniel used a Gaussian approach (mean ± N·σ) while I use Tukey ([Q1 − 1.5·IQR, Q3 + 1.5·IQR]). We both remove length, amplitude, and peak position.

- For a reason I don't understand, Daniel arbitrarily threw away all but the first 10% of the bin for ECG, but not PPG? I don't do this.

- 
