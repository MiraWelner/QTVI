clear;
close all;

props = readProps('config.txt');
features_path = props('windowedFeatures_input');
outputLoc = props('windowedFeatures_output');

window_size_sec = 300;


headers = [];

first = 0;
[analysisFiles] = windowFeaturesSetup(features_path);
time = 0;
parfor i = 1:size(analysisFiles,1)
    if isfile(fullfile(outputLoc, [analysisFiles{i,1} '_windowed_feature.mat']))
        continue
    end

    try
%     if strcmp(analysisFiles{i, 1}, '3013138_20111121')
        
    tStart = tic;

%     avg_time = time/i;
%     disp(['Avg Time (s): ' num2str(avg_time)]);

%     disp(['Est finish (min): ' num2str((avg_time*(size(analysisFiles,1)-i))/60)]);

    disp(join(['Saving analysis of ' analysisFiles{i, 1}]));
    beats_flattened = load(analysisFiles{i, 2});
    beats_flattened = beats_flattened.beats_flattened;
    beats_flattened = rmfield(beats_flattened,'ppg_wout_noise');
    beats_flattened = rmfield(beats_flattened,'ppg_flat_time_msec');
    
    time_sec = beats_flattened.corrected_time_sec/1000;
    d = diff(time_sec);

    for x = length(d):-1:1
       if sign(d(x)) >= 0
           break
       end
    end
    max_time = time_sec(x+1);
    window_count = ceil(max_time/window_size_sec);
    
    windowed_feature = [];
    headers = fieldnames(beats_flattened);
    len = length(headers);
    
    dropBlacklist = {'sleep_stages','idx_begin','idx_end','idx_foot','idx_pos_slope','idx_systolic','idx_neg_slope_b4','idx_neg_slope_after','idx_diastolic','idx_dnotch','proportional_pulse_amp','correct_idx_begin','edge_beat_mask'};
    % fill out window features
    intrestHeaders = {};
    for h = 1:len
        if ~any(strcmp(dropBlacklist,headers{h}))
            intrestHeaders{end+1} = headers{h};
            windowed_feature.([headers{h} '_n']) =nan(window_count,1);
            windowed_feature.([headers{h} '_nn']) =nan(window_count,1);   
            if strcmp('adjusted_sleep_state',headers{h})
                windowed_feature.([headers{h} '_sleep_order_compressed']) = cell(window_count,1);
                windowed_feature.([headers{h} '_sleep_order_compressed_n']) = cell(window_count,1);
                windowed_feature.([headers{h} '_0']) = nan(window_count,1);
                windowed_feature.([headers{h} '_1']) = nan(window_count,1);
                windowed_feature.([headers{h} '_2']) = nan(window_count,1);
                windowed_feature.([headers{h} '_3']) = nan(window_count,1);
                windowed_feature.([headers{h} '_4']) = nan(window_count,1);
                windowed_feature.([headers{h} '_5']) = nan(window_count,1);
                windowed_feature.([headers{h} '_other']) = nan(window_count,1);
            elseif any(strcmp({'corrected_time_sec','sec_to_first_onset_of_sleep','sec_from_last_onset_of_sleep'},headers{h}))
                windowed_feature.(headers{h}) = nan(window_count,1);
            elseif any(strcmp({'error_ppg_segmentation','review_bad_ppg_template','review_bad_r_template'},headers{h}))
                windowed_feature.([headers{h} '_n0']) = nan(window_count,1);
                windowed_feature.([headers{h} '_n1']) = nan(window_count,1);
            else
                windowed_feature.([headers{h} '_avg']) = nan(window_count,1);
                windowed_feature.([headers{h} '_med']) = nan(window_count,1);
                windowed_feature.([headers{h} '_max']) = nan(window_count,1);
                windowed_feature.([headers{h} '_min']) = nan(window_count,1);
                windowed_feature.([headers{h} '_q1']) = nan(window_count,1);
                windowed_feature.([headers{h} '_q3']) = nan(window_count,1);
                windowed_feature.([headers{h} '_slo']) = nan(window_count,1);
            end
        end
    end
    
    win_beg = 0;
    win_end = window_size_sec;
    len = length(intrestHeaders);
    for win = 1:window_count
%         disp(num2str(win));
        for h = 1:len
            mask = (time_sec >= win_beg) & (time_sec < win_end);
            tmp = beats_flattened.(intrestHeaders{h});
            slice = tmp(mask);
            
            if isempty(slice)
                continue
            end
            
            windowed_feature.([intrestHeaders{h} '_n'])(win) = numel(slice);
            windowed_feature.([intrestHeaders{h} '_nn'])(win) = sum(isnan(slice));
            
            if strcmp('corrected_time_sec',intrestHeaders{h})
                windowed_feature.(intrestHeaders{h})(win) = win_beg;
            elseif any(strcmp({'sec_to_first_onset_of_sleep','sec_from_last_onset_of_sleep'},intrestHeaders{h}))
                windowed_feature.(intrestHeaders{h})(win) = slice(1);
            elseif strcmp('adjusted_sleep_state',intrestHeaders{h})
                [B, N, ~] = RunLength(slice);
                windowed_feature.([intrestHeaders{h} '_sleep_order_compressed']){win} = [B];
                windowed_feature.([intrestHeaders{h} '_sleep_order_compressed_n']){win} = [N];
                windowed_feature.([intrestHeaders{h} '_0'])(win) = sum(slice==0);
                windowed_feature.([intrestHeaders{h} '_1'])(win) = sum(slice==1);
                windowed_feature.([intrestHeaders{h} '_2'])(win) = sum(slice==2);
                windowed_feature.([intrestHeaders{h} '_3'])(win) = sum(slice==3);
                windowed_feature.([intrestHeaders{h} '_4'])(win) = sum(slice==4);
                windowed_feature.([intrestHeaders{h} '_5'])(win) = sum(slice==5);
                windowed_feature.([intrestHeaders{h} '_other'])(win) = sum(~ismember(slice,[0,1,2,3,4,5]));
            elseif any(strcmp({'error_ppg_segmentation','review_bad_ppg_template','review_bad_r_template'},intrestHeaders{h}))

                windowed_feature.([intrestHeaders{h} '_n0'])(win) = sum(slice==0);
                windowed_feature.([intrestHeaders{h} '_n1'])(win) = sum(slice==1);

            else
                windowed_feature.([intrestHeaders{h} '_avg'])(win) = nanmean(slice);
                windowed_feature.([intrestHeaders{h} '_med'])(win) = nanmedian(slice);
                windowed_feature.([intrestHeaders{h} '_max'])(win) = nanmax(slice);
                windowed_feature.([intrestHeaders{h} '_min'])(win) = nanmin(slice);
                
                windowed_feature.([intrestHeaders{h} '_q1'])(win) = quantile(slice,0.25);
                windowed_feature.([intrestHeaders{h} '_q3'])(win) = quantile(slice,0.75);
                if length(slice) == 1
                    windowed_feature.([intrestHeaders{h} '_slo'])(win) = 0;
                else
                    p = polyfit(time_sec(mask), slice, 1);
                    windowed_feature.([intrestHeaders{h} '_slo'])(win) = p(1);
                end
            end
        end
        win_beg = win_end;
        win_end = win_end + window_size_sec;
    end
%     disp('Saving...');
    parsave(fullfile(outputLoc, [analysisFiles{i,1} '_windowed_feature']), windowed_feature);
% 
    catch
        disp([num2str(i) ' error']);
    end

end

function parsave(location, windowed_feature)
    save(location, 'windowed_feature');
end