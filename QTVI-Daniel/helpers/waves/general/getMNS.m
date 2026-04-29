function [MNS MNSindex] = getMNS(forwardsDelim, backwardsDelim, offset, bpwaveform, sf)

    % Most Negative Slope between given indices to search between;
    if (offset == 1)
        backwardsDelim = backwardsDelim(2:end);
        forwardsDelim = forwardsDelim(1:end - 1);
    end

    n = length(forwardsDelim);
    MNS = zeros(1, n);
    MNSindex = zeros(1, n);

    for i = 1:n

        try
            indexdis = forwardsDelim(i):backwardsDelim(i);
            [slp slpin] = RollingDeriv(bpwaveform(indexdis), 3);
            [MNS(i) tempi] = min(slp);
            MNSindex(i) = forwardsDelim(i) + tempi;
        catch

            if (isempty(indexdis))
                indexdis = forwardsDelim(i):forwardsDelim(i) + sf / 15;
            end

            slp = diff(bpwaveform(indexdis));
            MNS(i) = min(slp) * sf;
            MNSindex(i) = forwardsDelim(i);
        end

    end

end
