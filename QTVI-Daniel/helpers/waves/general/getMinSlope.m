function [MNS MNSindex] = getMinSLope(forwardsDelim, backwardsDelim, offset, bpwaveform, sf, pd)

    % Most Negative Slope between given indices to search between;
    if (offset == 1)
        backwardsDelim = backwardsDelim(2:end);

        if (~pd)
            forwardsDelim = forwardsDelim(1:end - 1);
        end

    end

    n = length(forwardsDelim);
    MNS = zeros(1, n);
    MNSindex = zeros(1, n);

    for i = 1:n

        if (~(i == n && pd))

            try
                indexdis = forwardsDelim(i):backwardsDelim(i);
                [slp slpin] = RollingDeriv(bpwaveform(indexdis), 3);
                [MNS(i) tempi] = min(slp);
                MNSindex(i) = forwardsDelim(i) + tempi;
            catch

                try

                    if (isempty(indexdis))
                        indexdis = forwardsDelim(i):forwardsDelim(i) + sf / 15;
                    end

                    slp = diff(bpwaveform(indexdis));
                    MNS(i) = min(slp) * sf;
                    MNSindex(i) = forwardsDelim(i);
                catch
                    slp = diff(bpwaveform(forwardsDelim(i):forwardsDelim(i) + 1));
                    MNS(i) = slp * sf;
                    MNSindex(i) = forwardsDelim(i);
                end

            end

        else

            try
                indexdis = forwardsDelim(i):min([length(bpwaveform) forwardsDelim(i) + 2 * sf]);
                [slp slpin] = RollingDeriv(bpwaveform(indexdis), 3);
                [MNS(i) tempi] = min(slp);
                MNSindex(i) = forwardsDelim(i) + tempi;
            catch

                if (forwardsDelim(i) ~= length(bpwaveform))
                    indexdis = forwardsDelim(i):forwardsDelim(i) + 1;
                    [MNS(i)] = diff(bpwaveform(indexdis));
                    MNSindex(i) = forwardsDelim(i) + 0;
                else
                    indexdis = forwardsDelim(i) - 1:forwardsDelim(i);
                    [MNS(i)] = diff(bpwaveform(indexdis));
                    MNSindex(i) = forwardsDelim(i) + 0;
                end

            end

        end

    end
