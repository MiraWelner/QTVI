clear;
myNoiseDir = '/home/deeplab/Desktop/aha/exports/';
myEDFDir = '/hdd/data/mesa/EDF All';
exportDir = '/home/deeplab/Desktop/aha/';
myFiles = dir(fullfile(myNoiseDir,'*noiseSEG.mat')); 
sample_rate = 256;
targetLength = 30;



for k = 1:length(myFiles)
    baseFileName = myFiles(k).name;
    fullFileName = fullfile(myNoiseDir, baseFileName);
    [~, name, ~] = fileparts(fullFileName);
    name = name(1:end-9);
    fprintf(1, 'Reading %s\n', fullFileName);
    
    [hdr, data] = edfread(fullfile(myEDFDir,[name '_EDF/' name '.edf']), 'verbose',0,'targetSignals', [1, 23]); % 1 is ecg index, 23 pulse
    ecg = data(1, :);
    po = data(2, :);
    noiseSEG = load(fullFileName);
    noiseSEG = noiseSEG.noiseSEG;
    output = saveSegs(ecg,po,noiseSEG,sample_rate,targetLength,0);
    outfile = [fullfile(exportDir,name) '_output.mat'];
    fprintf(1, 'Saving %s\n', outfile);
    save(outfile,'output')
end

function [output] =  saveSegs(ecg,po,noiseSEG,sample_rate,targetLength,dbg)
    %% below copied from anneal segs
    min_exclusion_bin_size_seconds = 0;
    min_bin_size_mins = 1;

    ecg_time_seconds = (1:length(ecg)) / sample_rate;
    ecg_time_seconds = ecg_time_seconds';
    po_time_seconds = (0:length(po) - 1) / sample_rate;
    po_time_seconds = po_time_seconds';
    bin_size_samples = sample_rate * 60 * targetLength;

    %% last bin will be targetLength + whatever is left if the remander is <
    % min bin size. otherwise last bin is whatever length it equals.
    remander = mod(length(po), bin_size_samples);
    remander_mins = remander / 256 / 60; % in mins
    if remander_mins < min_bin_size_mins
        bin_count = floor(length(po) / bin_size_samples);
    else
        bin_count = ceil(length(po) / bin_size_samples);
    end

    %% calculate inclusive index of bin endings. modify last index to be 30min+
    % whatever was left at end of data
    bin_breaks = (bin_size_samples:bin_size_samples:length(po) - 1); % sub 1 for 1 indexing
    
    if length(bin_breaks) < bin_count
       bin_breaks = [bin_breaks length(po) - 1]; % add last bin
    else
        bin_breaks(end) = length(po) - 1; % make last bin include remander
    end
    bin_breaks = bin_breaks + 1; % add 1 to make 1 indexing work correct

    bin_times_seconds = [0; po_time_seconds(bin_breaks)];
    bin_times_seconds = diff(bin_times_seconds);

    %% ignore exclusions < min_exclusion_bin_size_seconds
    exclusions_length_seconds = zeros(size(noiseSEG, 1), 1);
    greater_then_min = false(size(noiseSEG, 1), 1);

    for i = 1:size(noiseSEG, 1)
        exclusions_length_seconds(i) = noiseSEG(i, 2) - noiseSEG(i, 1);
        greater_then_min(i) = exclusions_length_seconds(i) >= min_exclusion_bin_size_seconds;
    end

    exclusions_seconds = noiseSEG(greater_then_min, :); % only take exclusions > min
    %exclusions_length_seconds = exclusions_length_seconds(gt_5_min);

    %% since exclusions can be marked at any time the marked time must be 
    % matched to the nearest index in the actual data
    flat_exclusions = reshape(exclusions_seconds', 1, []); % flatten for easier looping
    exclusions_indexs = zeros(numel(flat_exclusions), 1);

    for i = 1:numel(flat_exclusions)
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

    
    %% above copied from anneal segs
    %% new save stuff
    exclusions_indexs = (exclusions_indexs / sample_rate) / 60;
    time_min = po_time_seconds / 60;
    time_min = time_min';
    
    if dbg == 1
        close all;
        plot(time_min,po);
        for i = 1:length(exclusions_indexs)
            x = [exclusions_indexs(i,1) exclusions_indexs(i,2) exclusions_indexs(i,2) exclusions_indexs(i,1)];
            y = [min(po) min(po) max(po) max(po)];
            p = patch(x,y,'r','LineStyle','none');
            alpha(p,0.15);
        end


        for i = targetLength:targetLength:time_min(end)
            vline(i);
        end
    end
    
    break_times = (targetLength:targetLength:time_min(end));
    begin_idx = 0;
    bin = 1;
    output = cell(length(targetLength:targetLength:time_min(end)),1);
    for i = targetLength:targetLength:time_min(end)
        end_idx = i*60*sample_rate;
        begin_idx = begin_idx *60*sample_rate;
        thirtyMinSection.ecg = ecg(begin_idx+1:end_idx+1);
        thirtyMinSection.po = po(begin_idx+1:end_idx+1);
        
        thirtyMinSection.excluded = exclusions(exclusions(:,3) == bin,1:2);
        thirtyMinSection.timeExcluded = sum(diff((thirtyMinSection.excluded /256/60)'));
        thirtyMinSection.timeGood = targetLength - thirtyMinSection.timeExcluded;
        thirtyMinSection.precentGood = 1 - thirtyMinSection.timeExcluded/targetLength;

        output{bin} = thirtyMinSection;
        begin_idx = i;
        bin = bin + 1;
    end
end