function [MPS MPSindex] = getMPS(forwardsDelim, backwardsDelim, bpwaveform, sf)
    % Most Positive Slope between given indices to search between;
    n = length(forwardsDelim);
    MPS = zeros(1, n);
    MPSindex = zeros(1, n);

    for i = 1:n
        indexdis = forwardsDelim(i):backwardsDelim(i);

        if (length(indexdis) > 3)
            [slp slpin] = RollingDeriv(bpwaveform(indexdis), 3);
            [MPS(i), tempi] = max(slp);
            MPSindex(i) = forwardsDelim(i) + slpin(tempi);
            MPS(i) = MPS(i) * sf;
        else
            slp = mean(diff(indexdis));
            slpin = backwardsDelim(i);
            [MPS(i), tempi] = max(slp);
            MPSindex(i) = forwardsDelim(i) + tempi;
            MPS(i) = MPS(i) * sf;
        end

    end

end
