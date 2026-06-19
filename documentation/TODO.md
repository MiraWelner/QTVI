## Deep Feedback

### Now

1. Check if Mohamed's computer needs video driver update, also check display

### Later

1. Make the bin sizes in the QTVI markings variable based on how much the signal changes

2. Make a pacemaker detector by verifying if the STD = 0

3. This also helps you measure the noise of the system if the STD is slightly greater than 0.

4. Use an FFT to detect where there is noise, and make a dynamic filter, this prevents unnecessary messing up of the T wave.

## Team Feedback

1. Blanking thresholding not always working, sometimes very low peaks like t wave peaks selected  - i think range is huge or infinite due to specific bug?

2. Fix scaling is not intuitively designed

3. Since the max of threshold is determined by median of previous segment, changing thresholding can actually screw things up in later areas
