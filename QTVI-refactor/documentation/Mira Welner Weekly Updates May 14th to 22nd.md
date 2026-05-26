* Edits to the noise marking GUI
  
  * Adding peak marking and heartrate
    
    * The heartrate is now calculated by a minimum of 10s even if the window currently displayed is shorter than 10s
  
  * Add hotkeys (1-9) to the noise marking
  
  * Fixed a bug in the scaling that caused it not to reset when checkbox unchecked
  
  * Made the annealing, peak finding, and templating run in the background to a file after it has been marked

* Edited the Templating GUI
  
  * Added in a shift to the PPG so it is accurately displayed as to the right of the ECG 
    
    * This doesn't work if the PPG signal is shifted to the right due to recording error(?)
  
  * Added right and left to the end of the displayed PPG signal so you can better mark the onset and end
  
  * Added 9 markers that the user can place around the PPG/ECG

* CHAOS data combination
  
  * I made an error where I did not combine CHAOS datsets that had a gap - this was fixed and now even if there is a gap between two CHAOS files of the same Patient ID and Room number, it is combined up until 8 hours
  
  * Got 248 chaos files form the server and ran a script to convert them to bin files
