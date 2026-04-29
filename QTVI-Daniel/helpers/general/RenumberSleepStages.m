function [sleepStages] = RenumberSleepStages(sleepStages)
    % according to documention below:
    % 'Sleep stages: hypnogram is broken into periods of 30 seconds
    % (epoch). 0 is wake stage, 1-4: sleep stage 1-4 and 5 is REM sleep'
    % https://github.com/nsrr/edf-editor-translator/wiki/Compumedics-Annotation-Format
    % This is not always the case as there are other numbers then 0-5.
    % Therefore this sets them all to 6 it also re orders the 0-5 numbers
    % and inverts them to make more visual sense between transitioning
    % sleep cycles
    
    sleepStages(sleepStages < 0 | sleepStages > 5) = 6; % unknown values
    sleepStages(sleepStages == 5) = inf; % temp move to re orginize
    sleepStages(sleepStages > 0 & sleepStages < 5) = sleepStages(sleepStages > 0 & sleepStages < 5) + 1; % move them up one
    sleepStages(sleepStages == inf) = 1; % temp move to re orginize
    sleepStages = sleepStages * -1;
end