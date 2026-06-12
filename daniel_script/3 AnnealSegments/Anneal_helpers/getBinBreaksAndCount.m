function [bin_breaks, bin_count] = getBinBreaksAndCount(ppg, bin_size_samples, ppgSampleRate, min_bin_size_mins)
%% last bin will be targetLength + whatever is left if the remainder is <
    % min bin size. otherwise last bin is whatever length it equals.
    remainder = mod(length(ppg), bin_size_samples);
    remainder_mins = remainder / ppgSampleRate / 60; % in mins
    bin_count = floor(length(ppg) / bin_size_samples) + (remainder_mins > min_bin_size_mins);

    %% calculate inclusive index of bin endings. modify last index to be targetlength min+
    % whatever was left at end of data
    bin_breaks = ones(length(bin_count),2);
    prev = 1;
    for x = 1:bin_count
        next = prev+bin_size_samples;
        if next > length(ppg)
            bin_breaks(x,1) = prev;
            bin_breaks(x,2) = length(ppg);
        else
            bin_breaks(x,1) = prev;
            bin_breaks(x,2) = next;
        end
        prev = next+1;
    end

%     bin_breaks = (bin_size_samples + 1:bin_size_samples:length(ppg)); % add 1 for 1 indexing 
%     if length(bin_breaks) < bin_count
%        bin_breaks = [bin_breaks length(ppg)]; % add last bin
%     else
%         bin_breaks(end) = length(ppg); % make last bin include remainder
%     end
end

