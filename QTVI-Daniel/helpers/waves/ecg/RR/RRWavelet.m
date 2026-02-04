function [ecgQRSBandSqlocs,ecgQRSBandSqPks] = RRWavelet(ecg,MinPeakDistance)
    % from mathworks
    % Level 4 Details ==> 11 - 22 Hz 
    % Level 5 Details ==> 5 Hz- 11 Hz
    ecgQRSBand = extractQRSband(ecg,5,4,5);
    
    ecgQRSBandSq = ecgQRSBand.^2;
    [ecgQRSBandSqPks,ecgQRSBandSqlocs] = findpeaks(ecgQRSBandSq,'MinPeakHeight',mean(ecgQRSBandSq)+std(ecgQRSBandSq)*2,...
        'MinPeakDistance',MinPeakDistance); 
    
    ecgQRSBandSqlocs = RPeakfromRWave(ecg,ecgQRSBandSqlocs);
end

