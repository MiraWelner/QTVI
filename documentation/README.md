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

If you are Mohamed, you can click the 'grid' button to get the grid :)

# The Template Marking GUI

The executable file can be found in  `executables/exe_files/qtvi_template_marking.exe`

To run it, download all of the contents of `exe_files` and run the executable from 

within the folder, it needs the `.dll` files and such to work.

Double click on exe file, select whether you want to run MESA, BITTIUM, or CHAOS. For input, navigate to where you put the output from the noise marking.

If the 9 markers are too crowded, you can un-check the box for the ppg or ecg markers

To mark a bad R, rightclick once, for a bad PPG, rightclick twice

The markers are:

1) P: ECG P wave peak

2) Q: Right before ECG Q wave

3) Tb: Before the T wave

4) Te: After the T wave

5) On: PPG Onset

6) 50: The 50% up the PPG peak

7) Pk: Peak of PPG

8) Dc: Dicrotic notch

9) En: End

    
