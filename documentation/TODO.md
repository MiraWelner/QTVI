1. Make a pacemaker detector by verifying if the STD = 0

2. This also helps you measure the noise of the system if the STD is slightly greater than 0.

3. Use an FFT to detect where there is noise, and make a dynamic filter, this prevents unnecessary messing up of the T wave.
4. Make the bin sizes in the QTVI markings variable based on how much the signal changes
