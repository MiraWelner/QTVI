## Deep Feedback

1. Make the bin sizes in the QTVI markings variable based on how much the signal changes

2. Make a pacemaker detector by verifying if the STD = 0

3. This also helps you measure the noise of the system if the STD is slightly greater than 0.

4. Use an FFT to detect where there is noise, and make a dynamic filter, this prevents unnecessary messing up of the T wave.

5. Make the file to bin read 

## Team Feedback

1. Blanking thresholding not always working, sometimes very low peaks like t wave peaks selected  - i think range is huge or infinite due to specific bug?

2. The logs ensure that a file that is previously worked on will not be automatically reloaded. However it doesn't allow you to come back to a file that you had previously been working on and were partway through.

3. Fix scaling is not intuitively designed

4. Since the max of threshold is determined by median of previous segment, changing thresholding can actually screw things up in later areas

5. Add ctrl-z to 'undo'

6. Weird scaling of the window on DK's computerf
