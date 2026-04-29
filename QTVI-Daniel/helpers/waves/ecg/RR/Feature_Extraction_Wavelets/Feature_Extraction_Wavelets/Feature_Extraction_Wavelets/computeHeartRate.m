%Copyright 2017 The MathWorks, Inc.
function hr = computeHeartRate(tm)
beat =  diff(tm);
hr = mean(60./beat);
