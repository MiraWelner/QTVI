function [regularityM] = regularity(data, sampleRate)
    % author: Michael Hautman

    winSize = sampleRate * 4;
    ThresSize = sampleRate * 40;

    %% detrend data
    [data, trend] = HautmanDetrend(data, sampleRate);

    %% filtering outliers to the variance
    threshold = sqrt(movvar(data, ThresSize));
    rollingSTD = sqrt(movvar(data, winSize));
    rollingSTD2 = sqrt(movvar(rollingSTD, winSize));

    regularityM = rollingSTD2 ./ threshold;

end
