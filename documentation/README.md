# File to Uniform Bin

The executable file can be found in  `executables/exe_files/convert_uniform_bin.exe`

To run it, download all of the contents of `exe_files` and run the executable from within the folder, it needs the `.dll` files and such to work.

Double click on exe file, select whether you want to run MESA, BITTIUM, or CHAOS. BITTIUM files are not uploaded to the sample folders, but the CHAOS files can be found in `executables/chaos_raw` and the MESA in `executables/mesa_raw`

Select the raw files as input and wherever you want to put the bin files as the output folder

If you wish to change the rate being upsampled, it is in config.csv

# The Noise Marking GUI

The executable file can be found in  `executables/exe_files/qtvi_template_marking.exe`

To run it, download all of the contents of `exe_files` and run the executable from within the folder, it needs the `.dll` files and such to work.

Double click on exe file, select whether you want to run MESA, BITTIUM, or CHAOS. For input, navigate to where you put the output from the uniform bin converter.

To switch between what you are marking, you can use the 1-9 hotkeys

To mark, drag across the marking row, or click the start/stop buttons if you are marking a large region

The peak identifier is not the full algorithm, so if it does badly, that doesn't mean the actual peak finder will find garbage

The 'grid' button will create a grid that has thick lines every 0.2 seconds, and thin lines every 0.04 seconds. It scales based on the amount of time displayed in the window.

If you want to change the threshold for how high a R peak must be compared to the span between the median R peak and the minimum value the previous 10 seconds, click 'threshold' and a box will pop up. You put in the threshold value, between 0 and 1, and the highlighted region will be set to that value.

If you want to change the blanking period, click the 'blanking' button and do the same.

When you are done noise and feature marking, click 'save' and the processing will being automatically. Some of the processing will need to be completed, some will be done in parallel when you are doing the next step. It will still take about a minute for the processing to be done.

Then, an new screen will pop up where you will mark the templates. There will be nine markers which you will use to label the templates.

If the 9 markers are too crowded, you can un-check the box for the ppg or ecg markers

To mark a bad R, rightclick once, for a bad PPG, rightclick twice

The markers are:

1. P: ECG P wave peak

2. Q: Right before ECG Q wave

3. Tb: Before the T wave

4. Te: After the T wave

5. On: PPG Onset

6. 50: The 50% up the PPG peak

7. Pk: Peak of PPG

8. Dc: Dicrotic notch

9. En: End

When you are done, you can save the file and the next file will pop up for annotation.
