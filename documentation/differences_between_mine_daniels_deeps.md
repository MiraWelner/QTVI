## Template marking:

- In Daniel's code, the templates are negative. You cannot see this because he doesn't display the y axis in his marking GUI. This is because of a bug (I think) where he shifts them, he accidentally shifts by the average peak, so the peaks are vertically alligned at -average peak. I removed this bug so they are aligned at the correct positive location (if average peak is positive - if for some reason it is negative then it'll be aligned at a negative y location)
