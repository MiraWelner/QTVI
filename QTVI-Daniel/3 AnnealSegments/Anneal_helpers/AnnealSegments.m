function [annealedSegments, final_bin_idx] = AnnealSegments(ppg, ppgSampleRate, ecg, ecgSampleRate, noiseSEG, scoring_epoch_size_sec, sleepStages, targetLength, rs, dbg_plot)

    %% changable
    min_exclusion_bin_size_seconds =  5;
    min_bin_size_mins = targetLength/2;

    %% determine ideal bin times
    %use 1:length instead of 0 for ecg to make calculating easier later
    ecg_time_seconds = (0:length(ecg)-1) / ecgSampleRate;
    ecg_time_seconds = ecg_time_seconds';
    po_time_seconds = (0:length(ppg) - 1) / ppgSampleRate;
    po_time_seconds = po_time_seconds';
    bin_size_samples = ppgSampleRate * 60 * targetLength;

    %% last bin will be targetLength + whatever is left if the remander is <
    % min bin size. otherwise last bin is whatever length it equals.
    remander = mod(length(ppg), bin_size_samples);
    remander_mins = remander / ppgSampleRate / 60; % in mins
    if remander_mins < min_bin_size_mins
        bin_count = floor(length(ppg) / bin_size_samples);
    else
        bin_count = ceil(length(ppg) / bin_size_samples);
    end

    %% calculate inclusive index of bin endings. modify last index to be targetlength min+
    % whatever was left at end of data
    bin_breaks = (bin_size_samples + 1:bin_size_samples:length(ppg)); % add 1 for 1 indexing
    
    if length(bin_breaks) < bin_count
       bin_breaks = [bin_breaks length(ppg)]; % add last bin
    else
        bin_breaks(end) = length(ppg); % make last bin include remander
    end

    bin_times_seconds = [0; po_time_seconds(bin_breaks)];
    bin_times_seconds = diff(bin_times_seconds);

    %% ignore exclusions < min_exclusion_bin_size_seconds
    exclusions_length_seconds = zeros(size(noiseSEG, 1), 1);
    greater_then_min = false(size(noiseSEG, 1), 1);

    parfor i = 1:size(noiseSEG, 1)
        exclusions_length_seconds(i) = noiseSEG(i, 2) - noiseSEG(i, 1);
        greater_then_min(i) = exclusions_length_seconds(i) >= min_exclusion_bin_size_seconds;
    end

    exclusions_seconds = noiseSEG(greater_then_min, :); % only take exclusions > min
    %exclusions_length_seconds = exclusions_length_seconds(gt_5_min);

    %% since exclusions can be marked at any time the marked time must be 
    % matched to the nearest index in the actual data
    flat_exclusions = reshape(exclusions_seconds', 1, []); % flatten for easier looping
    exclusions_indexs = zeros(numel(flat_exclusions), 1);

    parfor i = 1:numel(flat_exclusions)
        exclusions_indexs(i) = closest_idx(po_time_seconds, flat_exclusions(i));
    end

    %% determine which bin exclusion bound currently lies in
    exclusions_bin = RoundToClosestBin(bin_breaks, exclusions_indexs);

    % reshape to be 2xN instead of flat so they can be combined.
    exclusions_bin = reshape(exclusions_bin, 2, [])';
    exclusions_indexs = reshape(exclusions_indexs, 2, [])';
    
    %format: begin_idx, end_idx, begin_bin_idx, end_bin_idx
    exclusions = [exclusions_indexs, exclusions_bin];

    %% update exclusions which are overlapping bin breaks
    for i = 1:size(exclusions_indexs, 1)

        % occurs when an excluded zone overlaps bin breaks. When this happens,
        % split each portion into its respective bin
        if exclusions(i, 3) ~= exclusions(i, 4)
            temp_bin_end = exclusions(i, 2);
            temp_bin_num = exclusions(i, 4);

            % for bin # that are overlaped
            for bin = (exclusions(i, 3):1:exclusions(i, 4))

                if bin == exclusions(i, 3) % == first bin
                    exclusions(i, 2) = bin_breaks(bin);
                    exclusions(i, 4) = bin;
                elseif bin == temp_bin_num % == last bin
                    exclusions(end + 1, 1) = bin_breaks(bin - 1); % add one since bin breaks are inclusive
                    exclusions(end, 2) = temp_bin_end;
                    exclusions(end, 3) = bin;
                    exclusions(end, 4) = bin;
                else % middle bin
                    exclusions(end + 1, 1) = bin_breaks(bin - 1); % add one since bin breaks are inclusive
                    exclusions(end, 2) = bin_breaks(bin);
                    exclusions(end, 3) = bin;
                    exclusions(end, 4) = bin;
                end

            end

        end

    end

    % sort and drop last column as both bin indexs are the same now
    [~, order] = sort(exclusions(:, 1));
    exclusions = exclusions(order, 1:3);

    %% determine bins to update
    update_bins = unique(exclusions(:, 3));

    %% determine good bins beginings and ends
    good_bins = setdiff(1:bin_count, update_bins);
    % section_beg, section_end, move_direction(1 left, 2 right, 0 none),
    % move flag.
    good_sections = zeros(length(good_bins),4);

    for i = 1:numel(good_bins)
        cur_bin = good_bins(i);

        bin_begin = bin_breaks(cur_bin) - bin_size_samples;
        bin_end = bin_breaks(cur_bin);
        good_sections(i,1) = bin_begin;
        good_sections(i,2) = bin_end;
    end

    %% for bins that have exclusions (need to be updated)
    for i = 1:numel(update_bins)
        cur_bin = update_bins(i);

        bin_begin = bin_breaks(cur_bin) - bin_size_samples;
        bin_end = bin_breaks(cur_bin);
        bin_half = bin_end - (bin_size_samples / 2);

        exclusion_mask = exclusions(:, 3) == cur_bin;
        bin_exclusions = exclusions(exclusion_mask, 1:2);

        good = [bin_begin reshape(bin_exclusions', 1, []) bin_end];
        good = reshape(good, 2, [])';

        good_length = diff(good')';
        good = good(good_length ~= 0, :);
        good_length = good_length(good_length ~= 0);

        good_time_mins = (good_length / ppgSampleRate) / 60;
        time_mask = good_time_mins < min_bin_size_mins;
        
        % determine the movement direction
        % movement_dir 1 = left 2 = right 0=none;
        % if max value is neg that means move left otherwize move right
        % special case for beginning and end bins begining always moves
        % right end always left
        
        [m, movement_dir] = max((good - bin_half)');
        movement_dir = movement_dir';
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
        % and lies on the left half or right half respectively, ignore it.
        if cur_bin == 1 % beginning bin
            potential_additions = potential_additions(m>0,:);
        elseif cur_bin == bin_count % end bin
            potential_additions = potential_additions(m<0,:);
        end
        
        
        good_sections = [good_sections; potential_additions];
    end

    [~, order] = sort(good_sections(:, 1));
    good_sections = good_sections(order, :);

    %merge together index's which share a border and are both moving
    i = 1;

    while (i < size(good_sections, 1))

        if good_sections(i, 2) == good_sections(i + 1, 1) && good_sections(i, 4) ~= 0 && good_sections(i + 1, 4) ~= 0
            seg1_size_min = (good_sections(i, 2) - good_sections(i, 1)) / ppgSampleRate / 60;
            seg2_size_min = (good_sections(i + 1, 2) - good_sections(i + 1, 1)) / ppgSampleRate / 60;

            if seg1_size_min + seg2_size_min >= min_bin_size_mins
                good_sections(i, 3) = 0;
                good_sections(i, 4) = 0;
            else
                [~, idx] = max([seg1_size_min seg2_size_min]);
                idx = idx - 1;
                good_sections(i, 3) = good_sections(i + idx, 3);
                good_sections(i, 4) = 1;
            end

            good_sections(i, 2) = good_sections(i + 1, 2);

            good_sections(i + 1, :) = [];
        end

        i = i + 1;
    end

    q.po = [];
    q.ecg = [];
    q.sleepStage = [];
    q.scoring_epoch_size_sec = scoring_epoch_size_sec;
    final_bin_idx = {q};
    current_bin = 1;
    temp_bin = [];

    for i = 1:size(good_sections, 1)

        if good_sections(i, 4)
            %move left or right
            if good_sections(i, 3) == 1
                try
                    final_bin_idx{current_bin - 1}.po = [final_bin_idx{current_bin - 1}.po; good_sections(i, 1:2)];
                catch
                    final_bin_idx{current_bin}.po = [final_bin_idx{current_bin}.po; good_sections(i, 1:2)];
                end
            else
                temp_bin = [temp_bin; good_sections(i, 1:2)];
            end

        else
            final_bin_idx{current_bin}.po = [temp_bin; good_sections(i, 1:2)];
            current_bin = current_bin + 1;
            temp_bin = [];
        end

    end
    
    % remove duplicate indexs
    for i = 1:numel(final_bin_idx)
        if size(final_bin_idx{i}.po,1) > 1
            final_bin_idx{i}.po = MergeSegments(final_bin_idx{i}.po);
        end
    end
    
    for i = 1:numel(final_bin_idx)-1
        if final_bin_idx{i}.po(end,2) == final_bin_idx{i+1}.po(1,1)
            final_bin_idx{i}.po(end,2) = final_bin_idx{i}.po(end,2)-1;
        end
    end

    %match exclusion index to ecg
    for i = 1:numel(final_bin_idx)
        po_bin_indexs = final_bin_idx{i}.po;
        po_bin_times = (po_bin_indexs-1) / ppgSampleRate;
        flat = reshape(po_bin_times', 1, []);
        ecg_idx = zeros(numel(flat), 1);

        for q = 1:numel(flat)
            ecg_idx(q) = closest_idx(ecg_time_seconds, flat(q));
        end

        final_bin_idx{i}.ecg = reshape(ecg_idx, 2, [])';
    end
    
    sleep_stage_times = (1:length(sleepStages)) * scoring_epoch_size_sec;
    
%     for i = 1:numel(final_bin_idx)
%         time = (final_bin_idx{i}.po(2) - final_bin_idx{i}.po(1)) / POsampleRate;
%         
%         time
%     end
    
    %fill data from index's
    annealedSegments = cell(numel(final_bin_idx),1);
    for i = 1:numel(final_bin_idx)
        po_bin_indexs = final_bin_idx{i}.po;
        ecg_bin_indexs = final_bin_idx{i}.ecg;

        annealedSegments{i}.ppg_bin_indexs = final_bin_idx{i}.po;
        annealedSegments{i}.ecg_bin_indexs = final_bin_idx{i}.ecg;
        annealedSegments{i}.ppgSampleRate = ppgSampleRate;
        annealedSegments{i}.ecgSampleRate = ecgSampleRate;

        for w = 1:size(po_bin_indexs, 1)
            time = (po_bin_indexs(w, 1)-1:po_bin_indexs(w, 2)-1) / ppgSampleRate;
            if ~isempty(time)
                try
                    annealedSegments{i}.sleep_stages = [annealedSegments{i}.sleep_stages sleepStages(sleep_stage_times >= time(1) & sleep_stage_times <= time(end))'];
                catch
                    annealedSegments{i}.sleep_stages = [];
                    annealedSegments{i}.sleep_stages = [annealedSegments{i}.sleep_stages sleepStages(sleep_stage_times >= time(1) & sleep_stage_times <= time(end))'];
                end
            end
            
            if ~isempty(rs)
                mask = rs >= annealedSegments{i}.ecg_bin_indexs(w,1) & rs <= annealedSegments{i}.ecg_bin_indexs(w,2);
                rpeaks =  rs(mask);
                if ~isempty(rpeaks)
                    dif = rpeaks(1) - annealedSegments{i}.ecg_bin_indexs(w,1);
                    rpeaks = rpeaks - rpeaks(1) + 1 + dif;
                    try
                        annealedSegments{i}.r_peaks = [annealedSegments{i}.r_peaks rpeaks];
                    catch
                        annealedSegments{i}.r_peaks = [];
                        annealedSegments{i}.r_peaks = [annealedSegments{i}.r_peaks rpeaks];
                    end
                end
            end 
            
%             try
%                 annealedSegments{i}.ppg_bin_indexs_adjusted = [annealedSegments{i}.ppg_bin_indexs_adjusted; annealedSegments{i}.po_bin_indexs - annealedSegments{i}.ppg_bin_indexs_adjusted(end,2)];
%                 annealedSegments{i}.ecg_bin_indexs_adjusted = [annealedSegments{i}.ecg_bin_indexs_adjusted; annealedSegments{i}.ecg_bin_indexs - annealedSegments{i}.ecg_bin_indexs_adjusted(end,2)];
% 
%             catch
%                 annealedSegments{i}.ppg_bin_indexs_adjusted = po_bin_indexs(w,:);
%                 annealedSegments{i}.ecg_bin_indexs_adjusted = ecg_bin_indexs(w,:);
%             end
            
            try
                annealedSegments{i}.po = [annealedSegments{i}.po ppg(po_bin_indexs(w, 1):po_bin_indexs(w, 2))];
            catch
                annealedSegments{i}.po = [];
                annealedSegments{i}.po = [annealedSegments{i}.po ppg(po_bin_indexs(w, 1):po_bin_indexs(w, 2))];
            end

            try
                annealedSegments{i}.ecg = [annealedSegments{i}.ecg ecg(ecg_bin_indexs(w, 1):ecg_bin_indexs(w, 2))];
            catch
                annealedSegments{i}.ecg = [];
                annealedSegments{i}.ecg = [annealedSegments{i}.ecg ecg(ecg_bin_indexs(w, 1):ecg_bin_indexs(w, 2))];
            end

        end
    end

    if dbg_plot==1
        % for plotting
        fig = figure('visible','off','Name','Annealed Segments');
        height = 5;
        flat_exclusions = reshape(exclusions(:, 1:2)', 1, []);
        flat_exclusions(2, :) = 1;
        temp_breaks = bin_breaks;
        temp_breaks(2, :) = 0;

        all_breaks = [temp_breaks flat_exclusions];
        all_breaks = all_breaks';
        all_breaks(end + 1, 1) = 1;
        [~, order] = sort(all_breaks(:, 1));
        all_breaks = all_breaks(order, :);

        [~, ind] = unique(all_breaks(:, 1));
        duplicate_ind = setdiff(1:size(all_breaks, 1), ind);
        dup_mask = duplicate_ind(diff(duplicate_ind) == 1) - 1;

        all_breaks = all_breaks(ind, :);
        all_breaks(dup_mask, 2) = 1;
        total_size_mins = sum(bin_times_seconds) / 60;
        rectangle('Position', [0 0 total_size_mins height]);

        bin_breaks_mins = cumsum(bin_times_seconds) / 60;

        %% final bins
       for z = 1:size(final_bin_idx{i}, 1)

            for z = 1:size(final_bin_idx{i}, 1)
                bin_indexs = final_bin_idx{i};
                half = ((bin_indexs.po(z, 1) + bin_indexs.po(z, 2)) / 2) / ppgSampleRate / 60;
                text(half, height * 0.9, ['New Bin' newline num2str(i)], 'Color', 'g','HorizontalAlignment','Center');
            end

        end

        for i = 1:size(good_sections, 1)
            half = ((good_sections(i, 1) + good_sections(i, 2)) / 2) / ppgSampleRate / 60;
            % old bin
            text(half, height * 0.75, ['Old Bin' newline num2str(i)], 'Color', 'r','HorizontalAlignment','Center');
            % move direction
            if good_sections(i, 3) == 0
            	text(half, height * 0.8, 'Stay','HorizontalAlignment','Center');

            elseif good_sections(i, 3) == 1
            	text(half, height * 0.8, 'Left','HorizontalAlignment','Center');
            else
            	text(half, height * 0.8, 'Right','HorizontalAlignment','Center');
            end
        end

        for i = 1:bin_count
            text(bin_breaks_mins(i) - 17, height * 0.5, ['Ideal' newline num2str(i)]);
        end

        %
        %         for i = 1:size(exclusions_indexs, 1)
        %             b = exclusions_indexs(i, 1) / sampleRate / 60;
        %             e = exclusions_indexs(i, 2) / sampleRate / 60;
        %             patch([b e e b], [0 0 2.5 2.5], 'g', 'FaceAlpha', .3);
        %             text(b, 1.25, num2str(i));
        %         end

        for i = 1:size(exclusions, 1)
            b = exclusions(i, 1) / ppgSampleRate / 60;
            e = exclusions(i, 2) / ppgSampleRate / 60;
            patch([b e e b], [0 0 height height], 'red', 'FaceAlpha', .3);
        end

        dif = diff(all_breaks(:, 1));

        for i = 1:numel(dif)
            half = dif(i) / 2;
            x = all_breaks(i + 1, 1) - half - half / 2;
            x = x / ppgSampleRate / 60;
            if dif(i) / ppgSampleRate / 60 < 1
                continue
            end
            t = text(x, height, [num2str(dif(i) / ppgSampleRate / 60)  ' min']);
            set(t, 'Rotation', 45);
        end

        xlim([0, total_size_mins]);
        vline(bin_breaks_mins, 'k:');
        ylim([0, height + .75]);

        xticks(0:30:total_size_mins);
        xlabel('Time (Seconds)');
        title(['Results of annealing process' newline 'press space to continue...']);
        set(gca, 'YTick', []);     
        ShowDbgPlots({fig});
    end
end