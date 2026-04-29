function [outputwaveform, trend] = HautmanDetrend(waveform, Fs)
    %   takes long sgolay fit, then runs a basic smoothing function on the
    %   output. uses that output as a "trend" and subtracts that trend from
    %   the original waveform.
    trend = sgolayfilt(waveform, 3, (4 * Fs) + 1);
    trend = smoothdata(trend, 'movmean', Fs);
    outputwaveform = waveform - trend + mean(waveform);
end
