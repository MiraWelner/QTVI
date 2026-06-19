# File to Uniform Bin

The executable file can be found in  `executables/exe_files/convert_uniform_bin.exe`

To run it, download all of the contents of `exe_files` and run the executable from within the folder, it needs the `.dll` files and such to work.

Double click on exe file, select whether you want to run MESA, BITTIUM, or CHAOS. BITTIUM files are not uploaded to the sample folders, but the CHAOS files can be found in `executables/chaos_raw` and the MESA in `executables/mesa_raw`

Select the raw files as input and wherever you want to put the bin files as the output folder

If you wish to change the rate being upsampled, it is in config.csv

# The Noise Marking GUI

The executable file can be found in  the zipped folder `exe_files`. To run it, download`exe_files.zip` and unzip it. The executable is called `noise_marking_gui.exe` and you can run it by double clicking. Run the executable from within the folder, it needs the `.dll` files and such to work.

Once you run the executable, select whether you want to run MESA, BITTIUM, or CHAOS. 

If you leave the columns `original_file_path` and `output_folder` empty in the `config.csv` file, then you will be prompted to select the input and output folders. If you like, you can also put the absolute paths to the input and output folders in the appropriate cells.

To switch between the type of annoation you are marking, you can use the 1-9 hotkeys

To mark, drag and release. If you want to select a long segment that requires scrolling, you can instead use the mark button and place a start and stop line a the start and stop of the marking. To remove a marking, right click on it. 

The peak identifier is not the full algorithm, so if it does badly, that doesn't mean the actual peak finder will find garbage. However if you want to improve it for the log, which will help improve the algorithm, you can highlight a blanking and threshold region by clicking the blanking and threshold button. Then drag over the region and a popup will appear where you can change the blanking and thresold button between 0 and 1.

Blanking is the distance after the R peak in which another R peak can NOT be found. Threshold is the height required for an R peak. Removing noise helps the algorithm because it doesn't use regions marked as noise to determine blanking and threshold.

The 'grid' button will create a grid that has thick lines every 0.2 seconds, and thin lines every 0.04 seconds. It scales based on the amount of time displayed in the window.

When you are done noise and feature marking, click 'save' and the processing will being automatically. Some of the processing will need to be completed, some will be done in parallel when you are doing the next step. It will still take about a minute for the processing to be done.

Then, a new screen will pop up where you will mark the templates. There will be nine markers which you will use to label the templates.

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
