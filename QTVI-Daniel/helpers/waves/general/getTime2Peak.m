function [timeTP PPindex] = getTime2Peak(foot, peak, percent, bpwaveform, Fs)
    n = length(foot);
    timeTP = zeros(1, n);
    PPindex = zeros(1, n);
    time = 1 / Fs;

    for i = 1:n
        indexdis = foot(i):peak(i);
        minimum = min(bpwaveform(indexdis));
        maximum = max(bpwaveform(indexdis));
        [~, index] = min(abs(bpwaveform(indexdis) - (maximum - minimum) * percent - minimum));

        try
            PPindex(i) = index + foot(i) - 1;
            timeTP(i) = index * time;
        catch
            PPindex(i) = foot(i);
            timeTP(i) = 0;
        end

    end

end
