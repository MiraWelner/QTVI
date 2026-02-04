function [markedContainers] = markForMovement(exclusions, bin_breaks, bin_count, sampleRate, bin_size_samples, min_bin_size_mins)
% markForMovement 
    %% for bins that have exclusions (need to be updated)
%% determine bins to update
    update_bins = unique(exclusions(:, 3));
    markedContainers = zeros(numel(update_bins) + 1, 4);
    count = 1;
 %% for bins that have exclusions (need to be updated)
    for i = 1:numel(update_bins)
        cur_bin = update_bins(i);
        bin_begin = bin_breaks(cur_bin) - bin_size_samples;
        bin_end = bin_breaks(cur_bin);
        bin_half = bin_end - (bin_size_samples / 2);
        exclusion_mask = exclusions(:, 3) == cur_bin;
        bin_exclusions = exclusions(exclusion_mask, 1:2);
        good = getExclusionIntervals(bin_begin, bin_end, bin_exclusions);
        good_length = diff(good')';
        good = good(good_length ~= 0, :);
        good_length = good_length(good_length ~= 0);        
        time_mask = ((good_length / sampleRate) / 60) < min_bin_size_mins;
        
        % determine the movement direction
        % movement_dir 1 = left 2 = right 0=none;
        % if max value is neg that means move left otherwize move right
        % smarkedContainerscial case for beginning and end bins begining always moves
        % right end always left        
        [m, movement_dir] = max((good - bin_half), [], 2);        
        if cur_bin == 1 % beginning bin
            movement_dir(m <= 0) = 2;
        elseif cur_bin == bin_count % end bin
            movement_dir(m >= 0) = 1;
        else
            movement_dir(m <= 0) = 1; % negative value
        end             
        movement_dir(~time_mask) = 0;
        move_flag = ones(length(movement_dir), 1);
        move_flag(~time_mask) = 0;
        
        potential_additions = [good movement_dir move_flag];
        
        % TEMPORARY? for the moment if the first or last bin has a good segment < min bin size
        % and lies on the left half or right half resmarkedContainersctively, ignore it.
        if cur_bin == 1 % beginning bin
            potential_additions = potential_additions(m>0,:);
        elseif cur_bin == bin_count % end bin
            potential_additions = potential_additions(m<0,:);
        end        
        newCount = size(potential_additions, 1);
        markedContainers(count: count + newCount - 1, :) = potential_additions;
        count = count + newCount;
    end
    markedContainers = markedContainers(1:count - 1, :);
end

