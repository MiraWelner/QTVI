function [rridx,rramps] = RRsimpleSquared(ecg,minDist)
    ecgSigSq =ecg.^2;
    try
        [rramps,rridx] = findpeaks(ecgSigSq,'MinPeakHeight',mean(ecgSigSq)+std(ecgSigSq)*2,'MinPeakDistance',minDist);
    catch
        rridx = [];
    end
    rridx = rridx';
end