function [flattened] = flatten_beat_idx(bin, processSegments)
    flattened = [];
    lens = cellfun(@(x) trylength(x),bin);
    beats = sum(lens);
    % make flattened array w/ proper length of beats
    idx1 = find((lens>0)==1,1);
    fields = fieldnames(bin{idx1});
    for y = 1:length(fields)
        if contains(fields{y},'sqi')
            if contains(fields{y},'labels')
                lbls = bin{idx1}.(fields{y});
                for z = 1:size(lbls,2)
                    flattened.(['sqi_' lbls{z}]) = nan(beats,1);
                end  
            end
            continue
        end
        flattened.(fields{y}) = nan(beats,1);
    end
    
    start = 0;
    prev = 1;
    for x = 1:length(bin)
        len = lens(x);
        if ~isempty(bin{x})
            if x > 1
                start = start + length(processSegments{x - 1}.po) - 1;
            end

            for y = 1:length(fields)
                if contains(fields{y},'idx')
                    flattened.(fields{y})(prev:prev+len-1) = bin{x}.(fields{y}) + start;
                elseif contains(fields{y},'labels')
                    continue
                elseif contains(fields{y},'sqi')
                    lbls = [];
                    for w = 1:length(fields)
                        if contains(fields{w},'labels')
                            lbls=bin{x}.(fields{w});
                            break
                        end
                    end
                    for z = 1:size(bin{x}.(fields{y}),2)
                        tmp = bin{x}.(fields{y});
                        flattened.(['sqi_' lbls{z}])(prev:prev+len-1) = tmp(:,z);
                    end
                else
                    flattened.(fields{y})(prev:prev+len-1) = bin{x}.(fields{y});
                end
            end
        end
        prev = prev + len;
    end

    begin_sec = flattened.idx_begin / processSegments{1}.ppgSampleRate;
    
    flattened.sec_valley_2_valley = [0; diff((flattened.idx_begin / processSegments{1}.ppgSampleRate))];
    flattened.sec_foot_2_foot = [0; diff((flattened.idx_foot / processSegments{1}.ppgSampleRate))];
    
    flattened.sec_tP_20_2_tP_20 = [0; diff(begin_sec + flattened.tP_20_x /processSegments{1}.ppgSampleRate)];
    flattened.sec_tP_50_2_tP_50 = [0; diff(begin_sec + flattened.tP_50_x / processSegments{1}.ppgSampleRate)];
    flattened.sec_tP_80_2_tP_80 = [0; diff(begin_sec + flattened.tP_80_x / processSegments{1}.ppgSampleRate)];
    flattened.sec_tP_20_inv_2_tP_20_inv = [0; diff(begin_sec + flattened.tP_20_x_inv / processSegments{1}.ppgSampleRate)];
    flattened.sec_tP_50_inv_2_tP_50_inv = [0; diff(begin_sec + flattened.tP_50_x_inv / processSegments{1}.ppgSampleRate)];
    flattened.sec_tP_80_inv_2_tP_80_inv = [0; diff(begin_sec + flattened.tP_80_x_inv / processSegments{1}.ppgSampleRate)];
    
    flattened.sec_pos_slope_2_pos_slope = [0; diff((flattened.idx_pos_slope / processSegments{1}.ppgSampleRate))];
    flattened.sec_systolic_2_systolic = [0; diff((flattened.idx_systolic / processSegments{1}.ppgSampleRate))];
    flattened.sec_neg_slope_b4_2_neg_slope_b4 = [0; (diff(flattened.idx_neg_slope_b4 / processSegments{1}.ppgSampleRate))];
    flattened.sec_neg_slope_after_2_neg_slope_after = [0; diff((flattened.idx_neg_slope_after / processSegments{1}.ppgSampleRate))];
    flattened.sec_diastolic_2_diastolic = [0; diff((flattened.idx_diastolic / processSegments{1}.ppgSampleRate))];
    flattened.sec_dnotch_2_dnotch = [0; diff((flattened.idx_dnotch / processSegments{1}.ppgSampleRate))];
    
    flattened.sec_tR_20_2_tR_20 = [0; diff(begin_sec + flattened.tR_20_x / processSegments{1}.ppgSampleRate)];
    flattened.sec_tR_50_2_tR_50 = [0; diff(begin_sec + flattened.tR_50_x / processSegments{1}.ppgSampleRate)];
    flattened.sec_tR_80_2_tR_80 = [0; diff(begin_sec + flattened.tR_80_x / processSegments{1}.ppgSampleRate)];
    flattened.sec_tR_20_inv_2_tR_20_inv = [0; diff(begin_sec + flattened.tR_20_x_inv / processSegments{1}.ppgSampleRate)];
    flattened.sec_tR_50_inv_2_tR_50_inv = [0; diff(begin_sec + flattened.tR_50_x_inv / processSegments{1}.ppgSampleRate)];
    flattened.sec_tR_80_inv_2_tR_80_inv = [0; diff(begin_sec + flattened.tR_80_x_inv / processSegments{1}.ppgSampleRate)];
    
    [B(:,1), B(:,2), ~] = RunLength(flattened.sleep_stages);
    sleepIDX = find(B(:,1) == 0,1);
    before_sleep_wake_count = B(sleepIDX,2);
    
    wakeIDX = find(B(:,1) ~= 0 & ~isnan(B(:,1)),1,'last') + 1;
    if wakeIDX > size(B,1)
        after_sleep_count = 0;
    else
        after_sleep_count = B(wakeIDX,2);
    end

    % sleep states go:
    % 0: Awake
    % -1: REM
    % -2 to -5:  NREM 1-5
    % nan: sleepstate unknown
    %
    % We want to group nrems as 1 and make 3 new classes with 5 total 'sleep states':
    % 0: awake before sleep
    % 1: NREM 1-4
    % 2: REM
    % 3: awake during sleep
    % 4: awake after sleep
    
    % 4=Awake after sleep without sleeping again. Set everything after last
    % wake to 4
    sleep_states = flattened.sleep_stages;
    if after_sleep_count > 0
        for i = length(sleep_states)-after_sleep_count+1:length(sleep_states)
            sleep_states(i) = 4;
        end
    end
    
    % Set all REM to 2
    sleep_states(sleep_states==-1) = 2;%2=REM
    
    % Set all NREM to 1
    sleep_states(sleep_states ~= 0 & sleep_states ~= 2 & sleep_states ~= 4) = 1;
    
    % Set all instances of awake to 3 so that awake during sleep is
    % labeled.
    sleep_states(sleep_states == 0) = 3;
    
    % Correct wake before sleep back to 0
    for i = 1:before_sleep_wake_count
        sleep_states(i) = 0;
    end
    
    flattened = correctTimes(flattened, processSegments, processSegments{1}.ppgSampleRate);
    flattened.ppg_flat_time_msec = (0:sum(cellfun(@(x) length(x.po), processSegments))) / processSegments{1}.ppgSampleRate;
    
    flattened.adjusted_sleep_state =  sleep_states;
    flattened.sec_to_first_onset_of_sleep = (flattened.corrected_time_sec - flattened.corrected_time_sec(before_sleep_wake_count)) * -1;
    flattened.sec_from_last_onset_of_sleep = (flattened.corrected_time_sec - flattened.corrected_time_sec(end - after_sleep_count));
    
end

function l = trylength(x)
    try
        l = length(x.area);
    catch
        l=0;
    end
end

function [flattened] = correctTimes(flattened,bins,samplingRate)
    begIdx = 1;
    flattened.correct_idx_begin = zeros(length(flattened.idx_begin),1);
    flattened.edge_beat_mask = zeros(length(flattened.idx_begin),1);
    for x = 1:length(bins)
        for y = 1:size(bins{x}.ppg_bin_indexs,1)

            %bins{x}.ppg_bin_indexs
            len = bins{x}.ppg_bin_indexs(y,2)-bins{x}.ppg_bin_indexs(y,1);
                    
            endIdx = begIdx + len;
            
            idx_mask = (flattened.idx_begin >= begIdx & flattened.idx_begin <= endIdx) ;
            flattened.correct_idx_begin(idx_mask) = flattened.idx_begin(idx_mask) - begIdx + bins{x}.ppg_bin_indexs(y,1);
            
            first = find(idx_mask==1,1);
            last = find(idx_mask==1,1,'last');

            flattened.edge_beat_mask(first) = 1;
            flattened.edge_beat_mask(last) = 1;

            begIdx = endIdx ;

        end
    end
    % remove first and last as they won't affect diff for sec_to_sec
    % measurements
    first = find(flattened.edge_beat_mask==1,1);
    last = find(flattened.edge_beat_mask==1,1,'last');
    flattened.edge_beat_mask(first) = 0;
    flattened.edge_beat_mask(last) = 0;
            
    flattened.corrected_time_sec = ((flattened.correct_idx_begin-1) / samplingRate);
end
