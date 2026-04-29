function [waveformDDPlus, waveformDD, waveformD] = doubleDerive(waveform)
    waveformD = diff(waveform);
    waveformD = [waveformD' NaN];

    waveformDD = diff(waveformD);
    waveformDD = [waveformDD NaN];
    waveformDD = sgolayfilt(waveformDD, 3, 11); %small amount of smoothing done

    %Perform the switch for positive and negative first-derivative values
    waveformDDPlus = waveformDD .* (waveformD > 0 & waveformDD > 0);
    waveformDDPlus = waveformDDPlus.^2;
end
