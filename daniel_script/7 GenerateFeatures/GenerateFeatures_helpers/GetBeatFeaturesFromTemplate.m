function [intraBeatMeasures, beats_total] = GetBeatFeaturesFromTemplate(sqi, threshold, ppg, sleepstates, pairs, ppgSamplingRate, ecgSamplingRate,dnotch_ratio_sp)

    if ~exist('dnotch_ratio_sp','var') == 1
        dnotch_ratio_sp = nan;
    end

    %% inner idx
    if size(pairs,1) >= 1
        idx_vallies = pairs(:,1);
    else
        idx_vallies = [];
    end
    if length(idx_vallies) > 1
        idx_begin = idx_vallies(1:end-1);
        idx_end = idx_vallies(2:end);
    %     idx_vallies = idx_vallies(sqi < threshold & sqi ~= 0);
    %     idx_begin = idx_begin(sqi < threshold & sqi ~= 0);
    %     idx_end = idx_end(sqi < threshold & sqi ~= 0);
    % 
    %     pairs = pairs(1:end-1,:);
    %     
    %     pairs = pairs(sqi < threshold & sqi ~= 0,:);
    %     len = length(sqi(sqi < threshold & sqi ~= 0));
%         idx_vallies = idx_vallies(sqi < threshold);
        idx_begin = idx_begin(sqi < threshold);
        idx_end = idx_end(sqi < threshold);

        pairs = pairs(1:end-1,:);

        pairs = pairs(sqi < threshold,:);
        len = length(sqi(sqi < threshold));

    else
        len=0;
        idx_begin = nan(len, 1);
        idx_end = nan(len, 1);
    end
    sleep_stages = nan(len, 1);
    idx_feets = nan(len, 1);
    
    tP_20_x = nan(len, 1);
    tP_50_x = nan(len, 1);
    tP_80_x = nan(len, 1);
    tP_20_x_inv = nan(len, 1);
    tP_50_x_inv = nan(len, 1);
    tP_80_x_inv = nan(len, 1);
    
    idx_pos_slopes = nan(len, 1);

    idx_systolic_peaks = nan(len, 1);
    
    idx_neg_slopes_before_dnotch = nan(len, 1);
    idx_dicrotic_notches = nan(len, 1);
    idx_diastolic_peaks = nan(len, 1);
    idx_neg_slopes_after_dnotch = nan(len, 1);
    
    tR_20_x = nan(len, 1);
    tR_50_x = nan(len, 1);
    tR_80_x = nan(len, 1);
    tR_20_x_inv = nan(len, 1);
    tR_50_x_inv = nan(len, 1);
    tR_80_x_inv = nan(len, 1);

    %% inner times
    msec_tP_50_2_first_valley = nan(len, 1);
    msec_tP_50_2_foot = nan(len, 1);
    
    msec_tP_50_2_tP_20 = nan(len, 1);
    msec_tP_50_2_tP_80 = nan(len, 1);
    msec_tP_50_2_tP_20_inv = nan(len, 1);
    msec_tP_50_2_tP_80_inv = nan(len, 1);

    msec_tP_50_2_pos_slope = nan(len, 1);
    msec_tP_50_2_systolic_peak = nan(len, 1);
    msec_tP_50_2_negslopes_pre_dnotch = nan(len, 1);
    msec_tP_50_2_dicrotic_notch = nan(len, 1);
    msec_tP_50_2_diastolic_peak = nan(len, 1);
    msec_tP_50_2_negslopes_post_dnotch = nan(len, 1);
    msec_tP_50_2_second_valley = nan(len, 1);

    msec_tP_50_2_tR_20 = nan(len, 1);
    msec_tP_50_2_tR_50 = nan(len, 1);
    msec_tP_50_2_tR_80 = nan(len, 1);
    msec_tP_50_2_tR_20_inv = nan(len, 1);
    msec_tP_50_2_tR_50_inv = nan(len, 1);
    msec_tP_50_2_tR_80_inv = nan(len, 1);

    
    %% inner ints
    msec_beat_length = nan(len, 1);
    msec_total_duration_tP_20 = nan(len, 1);
    msec_total_duration_tP_50 = nan(len, 1);
    msec_total_duration_tP_80 = nan(len, 1);
    msec_total_duration_tR_20 = nan(len, 1);
    msec_total_duration_tR_50 = nan(len, 1);
    msec_total_duration_tR_80 = nan(len, 1);
    
    %% amps raw
    proportional_pulse_amp = nan(len, 1);
    amp_raw_vallies = nan(len, 1);
    
    amp_raw_feets = nan(len, 1);
    amp_raw_tP_20 = nan(len, 1);
    amp_raw_tP_50 = nan(len, 1);
    amp_raw_tP_80 = nan(len, 1);
    amp_raw_tP_20_inv = nan(len, 1);
    amp_raw_tP_50_inv = nan(len, 1);
    amp_raw_tP_80_inv = nan(len, 1);
    
    amp_raw_pos_slopes = nan(len, 1);
    amp_raw_systolic_peaks = nan(len, 1);
    amp_raw_neg_slopes_pre_dnotch = nan(len, 1);
    amp_raw_dicrotic_notches = nan(len, 1);
    amp_raw_diastolic_peaks = nan(len, 1);
    amp_raw_neg_slopes_after_dnotch = nan(len, 1);
    
    amp_raw_tR_20 = nan(len, 1);
    amp_raw_tR_50 = nan(len, 1);
    amp_raw_tR_80 = nan(len, 1);
    amp_raw_tR_20_inv = nan(len, 1);
    amp_raw_tR_50_inv = nan(len, 1);
    amp_raw_tR_80_inv = nan(len, 1);
    
    %% amps baselined
    amp_baselined_feets = nan(len, 1);
    amp_baselined_tP_20 = nan(len, 1);
    amp_baselined_tP_50 = nan(len, 1);
    amp_baselined_tP_80 = nan(len, 1);
    amp_baselined_tP_20_inv = nan(len, 1);
    amp_baselined_tP_50_inv = nan(len, 1);
    amp_baselined_tP_80_inv = nan(len, 1);
    
    amp_baselined_pos_slopes = nan(len, 1);
    amp_baselined_systolic_peaks = nan(len, 1);
    amp_baselined_neg_slopes_pre_dnotch = nan(len, 1);
    amp_baselined_dicrotic_notches = nan(len, 1);
    amp_baselined_diastolic_peaks = nan(len, 1);
    amp_baselined_neg_slopes_after_dnotch = nan(len, 1);
    
    amp_baselined_tR_20 = nan(len, 1);
    amp_baselined_tR_50 = nan(len, 1);
    amp_baselined_tR_80 = nan(len, 1);
    amp_baselined_tR_20_inv = nan(len, 1);
    amp_baselined_tR_50_inv = nan(len, 1);
    amp_baselined_tR_80_inv = nan(len, 1);
    
    area_baselined = nan(len, 1);
    area = nan(len, 1);
    amp_delta_systolic = nan(len,1);
    abs_amp_foot = nan(len,1);
    abs_amp_peak = nan(len,1);


    
    %% R timings
    msec_R_2_first_valley = nan(len, 1);
    msec_R_2_foot = nan(len, 1);
    msec_R_2_tP_20 = nan(len, 1);
    msec_R_2_tP_50 = nan(len, 1);
    msec_R_2_tP_80 = nan(len, 1);
    msec_R_2_tP_20_inv = nan(len, 1);
    msec_R_2_tP_50_inv = nan(len, 1);
    msec_R_2_tP_80_inv = nan(len, 1);
    
    msec_R_2_pos_slope = nan(len, 1);
    msec_R_2_systolic_peak = nan(len, 1);
    msec_R_2_negslopes_pre_dnotch = nan(len, 1);
    msec_R_2_dicrotic_notch = nan(len, 1);
    
    msec_R_2_tR_20 = nan(len, 1);
    msec_R_2_tR_50 = nan(len, 1);
    msec_R_2_tR_80 = nan(len, 1);
    msec_R_2_tR_20_inv = nan(len, 1);
    msec_R_2_tR_50_inv = nan(len, 1);
    msec_R_2_tR_80_inv = nan(len, 1); 
    
    msec_R_2_diastolic_peak = nan(len, 1);
    msec_R_2_negslopes_post_dnotch = nan(len, 1);
    msec_R_2_second_valley = nan(len, 1);

    
    intraBeatMeasures = cell(len, 1);
    parfor i = 1:len
       
        beat_info = [];
        beat_info.id = i;

        beat_valley_one = idx_begin(i);
        beat_valley_two = idx_end(i);
        
        if beat_valley_two > length(ppg)
           beat_valley_two = length(ppg); 
        end
        

        %% find general beat info
        beat = ppg(beat_valley_one:beat_valley_two);
        if isempty(beat) || length(beat) < 4
           beat_info.sleep_stage = nan;        
           sleep_stages(i) = beat_info.sleep_stage;
           continue 
        end
        beat_info.sleep_stage = mode(sleepstates(beat_valley_one:beat_valley_two));
        sleep_stages(i) = beat_info.sleep_stage;

        beat_time_msec = (0:length(beat) - 1) / ppgSamplingRate * 1000;
        [~, beat_foot] = find_foot_pulseox(beat', 0);

        %% first half of wave up to peak

        % max amp
        [max_amp, max_peak] = max(beat(beat_foot:end));
        systolic_peak = beat_foot + max_peak - 1;
        
        % tP
        beat_zeroed = beat - beat(1);
        max_amp_zeroed = beat_zeroed(systolic_peak);
        dif = max_amp_zeroed - max_amp;
        max_tP20 = (max_amp_zeroed * .2) - dif;
        max_tP50 = (max_amp_zeroed * .5) - dif;
        max_tP80 = (max_amp_zeroed * .8) - dif;

        x = (1:length(beat(1:systolic_peak)))';
        y = beat(1:systolic_peak);
        
        
        if length(y) > 1
            l1 = [x y]';

            arr = [max_tP20, max_tP50, max_tP80];
            ttime_x = zeros(3, 1);
            ttime_y = zeros(3, 1);
            for q = 1:3
                l2 = [[1 length(beat(1:systolic_peak))]; [arr(q) arr(q)]];
                p = InterX(l1, l2);

                if isempty(p)
                    ttime_x(q) = nan;
                    ttime_y(q) = nan;
                    continue
                end
                %ttime_x(q) = round(p(1,end));
                ttime_y(q) = p(2,end);
                ttime_x(q) = p(1,end);
            end
            tP_20_X = ttime_x(1);
            tP_50_X = ttime_x(2);
            tP_80_X = ttime_x(3);
            
            tP_20_y = ttime_y(1);
            tP_50_y = ttime_y(2);
            tP_80_y = ttime_y(3);
        else
            tP_20_X = nan;
            tP_50_X = nan;
            tP_80_X = nan;
            
            tP_20_y = nan;
            tP_50_y = nan;
            tP_80_y = nan;
        end
        
        % tP-inv
        x = (1:length(beat(systolic_peak:end)))';
        y = beat(systolic_peak:end);
        
        if length(y) > 1
            l1 = [x y]';

            arr = [max_tP20, max_tP50, max_tP80];
            ttime_x = zeros(3, 1);
            ttime_y = zeros(3, 1);

            for q = 1:3
                l2 = [[1 length(y)]; [arr(q) arr(q)]];
                p = InterX(l1, l2);

                if isempty(p)
                    ttime_x(q) = nan;
                    ttime_y(q) = nan;
                    continue
                end
                %ttime_x(q) = round(p(1,end));
                ttime_y(q) = p(2,end);
                ttime_x(q) = p(1,end);
            end
            tP_20_inv_x = systolic_peak + ttime_x(1) - 1;
            tP_50_inv_x = systolic_peak + ttime_x(2) - 1;
            tP_80_inv_x = systolic_peak + ttime_x(3) - 1;
            
            tP_20_inv_y = ttime_y(1);
            tP_50_inv_y = ttime_y(2);
            tP_80_inv_y = ttime_y(3);
        else
            tP_20_inv_x  = nan;
            tP_50_inv_x  = nan;
            tP_80_inv_x  = nan;
            
            tP_20_inv_y = nan;
            tP_50_inv_y = nan;
            tP_80_inv_y = nan;
        end
        

        % pos slope
        % assumption: biggest slope between foot and peak
        [~, pos_slope_idx] = max(diff(beat(beat_foot:beat_foot + max_peak - 1)));
        if isempty(pos_slope_idx)
           pos_slope_idx = nan; 
        end
        max_positive_slope = beat_foot + pos_slope_idx - 1;

        % dnotch
        if isnan(dnotch_ratio_sp)
            dicrotic_notch = dumbDicrotic(beat);
        else
            dicrotic_notch = dumbDicrotic(beat,dnotch_ratio_sp);
        end

        if dicrotic_notch > length(beat) || dicrotic_notch <= systolic_peak || isnan(dicrotic_notch)
            dicrotic_notch = nan;
        end
        
        if isnan(dicrotic_notch)
            max_neg_slope_before_dicrotic_notch = nan;
            max_neg_slope_after_dicrotic_notch = nan;
            diastolic_peak = nan;
        else
            % b4 dnotch
            [~, b4dnotch] = max(-diff(beat(systolic_peak:dicrotic_notch)));
            if isempty(b4dnotch) || b4dnotch >= length(beat) || b4dnotch == 1
                b4dnotch = nan;
            end
            max_neg_slope_before_dicrotic_notch = systolic_peak + b4dnotch - 1;


            % after dnotch
            [~, afterNotch] = max(-diff(beat(dicrotic_notch:end)));
            if isempty(afterNotch) || afterNotch >= length(beat) || afterNotch == 1 
                afterNotch = nan;
            end
            max_neg_slope_after_dicrotic_notch = dicrotic_notch + afterNotch - 1;

            % diastolic_peak
            if isnan(max_neg_slope_after_dicrotic_notch)
                diastolic_peak = nan;
            else
                [~, dpeak] = max(beat(dicrotic_notch:max_neg_slope_after_dicrotic_notch));
                if dpeak == 1 || dpeak >= length(beat) || isempty(dpeak)
                    diastolic_peak = nan;
                    
                else
                    diastolic_peak = dicrotic_notch+dpeak-1;
                end
            end
        end
        
        % tR
        beat_zeroed = beat - beat(end);
        max_amp_zeroed = beat_zeroed(systolic_peak);
        dif = max_amp_zeroed - max_amp;
        max_tR20 = (max_amp_zeroed * .2) - dif;
        max_tR50 = (max_amp_zeroed * .5) - dif;
        max_tR80 = (max_amp_zeroed * .8) - dif;
        x = (1:length(beat(systolic_peak:end)))';
        y = beat(systolic_peak:end);
        
        if length(y) > 1
            l1 = [x y]';

            arr = [max_tR80, max_tR50, max_tR20];
            ttime_x = zeros(3, 1);
            ttime_y = zeros(3, 1);

            for q = 1:3
                l2 = [[1 length(y)]; [arr(q) arr(q)]];
                p = InterX(l1, l2);

                if isempty(p)
                    ttime_x(q) = nan;
                    ttime_y(q) = nan;
                    continue
                end
                %ttime_x(q) = round(p(1,end));
                ttime_y(q) = p(2,end);
                ttime_x(q) = p(1,end);
            end          
            tR_20_X = systolic_peak + ttime_x(1) - 1;
            tR_50_X = systolic_peak + ttime_x(2) - 1;
            tR_80_X = systolic_peak + ttime_x(3) - 1;
            
            tR_20_y = ttime_y(1);
            tR_50_y = ttime_y(2);
            tR_80_y = ttime_y(3);
        else
            tR_20_X  = nan;
            tR_50_X  = nan;
            tR_80_X  = nan;
            
            tR_20_y  = nan;
            tR_50_y  = nan;
            tR_80_y  = nan;
        end
        
        % tR-inv
        x = (1:length(beat(1:systolic_peak)))';
        y = beat(1:systolic_peak);
        
        if length(y) > 1
            l1 = [x y]';

            arr = [max_tR80, max_tR50, max_tR20];
            ttime_x = zeros(3, 1);
            ttime_y = zeros(3, 1);

            for q = 1:3
                l2 = [[1 length(y)]; [arr(q) arr(q)]];
                p = InterX(l1, l2);

                if isempty(p)
                    ttime_x(q) = nan;
                    ttime_y(q) = nan;
                    continue
                end
                %ttime_x(q) = round(p(1,end));
                ttime_y(q) = p(2,end);
                ttime_x(q) = p(1,end);
            end
            tR_20_inv_x = ttime_x(1);
            tR_50_inv_x = ttime_x(2);
            tR_80_inv_x = ttime_x(3);
            
            tR_20_inv_y = ttime_y(1);
            tR_50_inv_y = ttime_y(2);
            tR_80_inv_y = ttime_y(3);
        else
            tR_20_inv_x  = nan;
            tR_50_inv_x  = nan;
            tR_80_inv_x  = nan;
            
            tR_20_inv_y = nan;
            tR_50_inv_y = nan;
            tR_80_inv_y = nan;
        end

        
        beat_zeroed = beat - min(beat);
        area(i) = trapz(0:length(beat_zeroed)-1,beat_zeroed);
        try
            tmp = InterX([[1 length(beat)]; [beat(1) beat(end)]], [[systolic_peak systolic_peak]; [max(beat)+1 min(beat)-1]] );
            amp_delta_systolic(i) = abs(getVal(beat,systolic_peak) - tmp(1,end)) ;
            
            slope = (beat(1)-beat(end))/(0-length(beat)-1);
            b = beat(1)-(slope*0);
            ys = nan(length(beat),1);
            for iter = 1:length(beat)
                ys(iter) = (slope*iter)+b;
            end
            area_baselined(i) = trapz(1:length(beat),abs(beat-ys));
        catch
            amp_delta_systolic(i) = nan;
            area_baselined(i) = nan;
        end
        
        beat_norm = beat/min(beat);

        abs_amp_foot(i) = getVal(beat_norm,beat_foot);
        abs_amp_peak(i) = getVal(beat_norm,systolic_peak);

        beat_info.sampling_rate = ppgSamplingRate;
        beat_info.area_baselined = area_baselined(i);
        beat_info.area = area(i);
        beat_info.amp_delta_systolic = amp_delta_systolic(i);
        beat_info.abs_amp_foot = abs_amp_foot(i);
        beat_info.abs_amp_peak = abs_amp_peak(i);

        beat_info.min_amplitude_one = 1;
        beat_info.min_amplitude_two = length(beat);
        beat_info.beat_time_msec = beat_time_msec;
        beat_info.beat_foot = beat_foot;
        beat_info.systolic_peak = systolic_peak;
        beat_info.tP_20_x = tP_20_X;
        beat_info.tP_50_x = tP_50_X;
        beat_info.tP_80_x = tP_80_X;
        beat_info.tP_20_x_inv = tP_20_inv_x;
        beat_info.tP_50_x_inv = tP_50_inv_x;
        beat_info.tP_80_x_inv = tP_80_inv_x;
        beat_info.tP_20_y = tP_20_y;
        beat_info.tP_50_y = tP_50_y;
        beat_info.tP_80_y = tP_80_y;
        beat_info.tP_20_y_inv = tP_20_inv_y;
        beat_info.tP_50_y_inv = tP_50_inv_y;
        beat_info.tP_80_y_inv = tP_80_inv_y;
        beat_info.max_positive_slope = max_positive_slope;
        beat_info.dicrotic_notch = dicrotic_notch;
        beat_info.max_neg_slope_before_dicrotic_notch = max_neg_slope_before_dicrotic_notch;
        beat_info.max_neg_slope_after_dicrotic_notch = max_neg_slope_after_dicrotic_notch;
        beat_info.diastolic_peak = diastolic_peak;
        beat_info.tR_20_x = tR_20_X;
        beat_info.tR_50_x = tR_50_X;
        beat_info.tR_80_x = tR_80_X;
        beat_info.tR_20_x_inv = tR_20_inv_x;
        beat_info.tR_50_x_inv = tR_50_inv_x;
        beat_info.tR_80_x_inv = tR_80_inv_x;     
        beat_info.tR_20_y = tR_20_y;
        beat_info.tR_50_y = tR_50_y;
        beat_info.tR_80_y = tR_80_y;
        beat_info.tR_20_y_inv = tR_20_inv_y;
        beat_info.tR_50_y_inv = tR_50_inv_y;
        beat_info.tR_80_y_inv = tR_80_inv_y;  
        intraBeatMeasures{i} = beat_info;        

        %plt(beat,beat_info);
        
        %% index adjusted to beginning
        idx_feets(i) = beat_valley_one + beat_info.beat_foot - 1;
       
        idx_pos_slopes(i) = beat_valley_one + beat_info.max_positive_slope - 1;
        idx_systolic_peaks(i) =  beat_valley_one + beat_info.systolic_peak - 1;
        idx_neg_slopes_before_dnotch(i) =  beat_valley_one + beat_info.max_neg_slope_before_dicrotic_notch - 1;
        idx_dicrotic_notches(i) = beat_valley_one + beat_info.dicrotic_notch - 1;
        idx_diastolic_peaks(i) = beat_valley_one + beat_info.diastolic_peak - 1;
        idx_neg_slopes_after_dnotch(i) = beat_valley_one + beat_info.max_neg_slope_after_dicrotic_notch - 1;
        
        
                
        tP_20_x(i) = beat_info.tP_20_x - 1;
        tP_50_x(i) = beat_info.tP_50_x - 1;
        tP_80_x(i) = beat_info.tP_80_x - 1;
        tP_20_x_inv(i) = beat_info.tP_20_x_inv - 1;
        tP_50_x_inv(i) = beat_info.tP_50_x_inv - 1;
        tP_80_x_inv(i) = beat_info.tP_80_x_inv - 1;
        
        tR_20_x(i) = beat_info.tR_20_x - 1;
        tR_50_x(i) = beat_info.tR_50_x - 1;
        tR_80_x(i) = beat_info.tR_80_x - 1;
        tR_20_x_inv(i) = beat_info.tR_20_x_inv - 1;
        tR_50_x_inv(i) = beat_info.tR_50_x_inv - 1;
        tR_80_x_inv(i) = beat_info.tR_80_x_inv - 1;

        %% t-50 times
        t50 = (beat_info.tP_50_x/ppgSamplingRate)*1000;
        msec_tP_50_2_first_valley(i) =  beat_time_msec(1) - t50;
        msec_tP_50_2_foot(i) = (beat_info.beat_foot/ppgSamplingRate)*1000 - t50;
        
        msec_tP_50_2_tP_20(i) = (beat_info.tP_20_x/ppgSamplingRate)*1000 - t50;
        msec_tP_50_2_tP_80(i) = (beat_info.tP_80_x/ppgSamplingRate)*1000 - t50;
        msec_tP_50_2_tP_20_inv(i) = (beat_info.tP_20_x_inv/ppgSamplingRate)*1000 - t50;
        msec_tP_50_2_tP_80_inv(i) = (beat_info.tP_80_x_inv/ppgSamplingRate)*1000 - t50;
        
        msec_tP_50_2_pos_slope(i) = (beat_info.max_positive_slope/ppgSamplingRate)*1000 - t50;
        msec_tP_50_2_systolic_peak(i) = (beat_info.systolic_peak/ppgSamplingRate)*1000 - t50;
        msec_tP_50_2_negslopes_pre_dnotch(i) = (beat_info.max_neg_slope_before_dicrotic_notch/ppgSamplingRate)*1000 - t50;
        msec_tP_50_2_dicrotic_notch(i) = (beat_info.dicrotic_notch/ppgSamplingRate)*1000 - t50;
        msec_tP_50_2_diastolic_peak(i) = (beat_info.diastolic_peak/ppgSamplingRate)*1000 - t50;
        msec_tP_50_2_negslopes_post_dnotch(i) = (beat_info.max_neg_slope_after_dicrotic_notch/ppgSamplingRate)*1000 - t50;
%       msec_tP_50_2_second_valley(i) = (beat_info.min_amplitude_two/ppgSamplingRate)*1000 - t50;
%         if isnan(msec_tP_50_2_second_valley(i))
%            disp('why') 
%         end
                
        msec_tP_50_2_tR_20(i) = (beat_info.tR_20_x/ppgSamplingRate)*1000 - t50; 
        msec_tP_50_2_tR_50(i) = (beat_info.tR_50_x/ppgSamplingRate)*1000 - t50;
        msec_tP_50_2_tR_80(i) = (beat_info.tR_80_x/ppgSamplingRate)*1000 - t50;
        msec_tP_50_2_tR_20_inv(i) = (beat_info.tR_20_x_inv/ppgSamplingRate)*1000 - t50;
        msec_tP_50_2_tR_50_inv(i) = (beat_info.tR_50_x_inv/ppgSamplingRate)*1000 - t50;
        msec_tP_50_2_tR_80_inv(i) = (beat_info.tR_80_x_inv/ppgSamplingRate)*1000 - t50;

        %% inter beat durations
        msec_beat_length(i) = beat_time_msec(beat_info.min_amplitude_two);
        msec_total_duration_tP_20(i) = abs((beat_info.tP_20_x/ppgSamplingRate) - (beat_info.tP_20_x_inv/ppgSamplingRate))*1000;
        msec_total_duration_tP_50(i) = abs((beat_info.tP_50_x/ppgSamplingRate) - (beat_info.tP_50_x_inv/ppgSamplingRate))*1000;
        msec_total_duration_tP_80(i) = abs((beat_info.tP_80_x/ppgSamplingRate) - (beat_info.tP_80_x_inv/ppgSamplingRate))*1000;
        msec_total_duration_tR_20(i) = abs((beat_info.tR_20_x/ppgSamplingRate) - (beat_info.tR_20_x_inv/ppgSamplingRate))*1000;
        msec_total_duration_tR_50(i) = abs((beat_info.tR_50_x/ppgSamplingRate) - (beat_info.tR_50_x_inv/ppgSamplingRate))*1000;
        msec_total_duration_tR_80(i) = abs((beat_info.tR_80_x/ppgSamplingRate) - (beat_info.tR_80_x_inv/ppgSamplingRate))*1000;

        %% amplitudes
        %raw
        proportional_pulse_amp(i) = (beat(beat_info.systolic_peak) - getVal(beat,beat_info.diastolic_peak)) / beat(beat_info.systolic_peak);
        amp_raw_vallies(i) = getVal(beat,beat_info.min_amplitude_one);
        amp_raw_feets(i) = getVal(beat,beat_info.beat_foot);
        
        amp_raw_tP_20(i) = beat_info.tP_20_y;
        amp_raw_tP_50(i) = beat_info.tP_50_y;
        amp_raw_tP_80(i) = beat_info.tP_80_y;
        amp_raw_tP_20_inv(i) = beat_info.tP_20_y_inv;
        amp_raw_tP_50_inv(i) = beat_info.tP_50_y_inv;
        amp_raw_tP_80_inv(i) = beat_info.tP_80_y_inv;
        
        amp_raw_pos_slopes(i) = getVal(beat,beat_info.max_positive_slope);
        amp_raw_systolic_peaks(i) = getVal(beat,beat_info.systolic_peak);
        
        amp_raw_tR_20(i) = beat_info.tR_20_y;
        amp_raw_tR_50(i) = beat_info.tR_50_y;
        amp_raw_tR_80(i) = beat_info.tR_80_y;
        amp_raw_tR_20_inv(i) = beat_info.tR_20_y_inv;
        amp_raw_tR_50_inv(i) = beat_info.tR_50_y_inv;
        amp_raw_tR_80_inv(i) = beat_info.tR_80_y_inv;
        
        amp_raw_neg_slopes_pre_dnotch(i) = getVal(beat,beat_info.max_neg_slope_before_dicrotic_notch);
        amp_raw_dicrotic_notches(i) = getVal(beat,beat_info.dicrotic_notch);
        amp_raw_diastolic_peaks(i) = getVal(beat,beat_info.diastolic_peak);
        amp_raw_neg_slopes_after_dnotch(i) = getVal(beat,beat_info.max_neg_slope_after_dicrotic_notch);

        %baselined
        amp_baselined_feets(i) = getVal(beat,beat_info.beat_foot) - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_tP_20(i) = beat_info.tP_20_y - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_tP_50(i) = beat_info.tP_50_y - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_tP_80(i) = beat_info.tP_80_y -getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_tP_20_inv(i) = beat_info.tP_20_y_inv - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_tP_50_inv(i) = beat_info.tP_50_y_inv - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_tP_80_inv(i) = beat_info.tP_80_y_inv - getVal(beat,beat_info.min_amplitude_one);
        
        amp_baselined_pos_slopes(i) = getVal(beat,beat_info.max_positive_slope) - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_systolic_peaks(i) = getVal(beat,beat_info.systolic_peak) - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_neg_slopes_pre_dnotch(i) = getVal(beat,beat_info.max_neg_slope_before_dicrotic_notch) - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_dicrotic_notches(i) = getVal(beat,beat_info.dicrotic_notch) - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_diastolic_peaks(i) = getVal(beat,beat_info.diastolic_peak) - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_neg_slopes_after_dnotch(i) = getVal(beat,beat_info.max_neg_slope_after_dicrotic_notch) - getVal(beat,beat_info.min_amplitude_one);
        
        amp_baselined_tR_20(i) = beat_info.tR_20_y - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_tR_50(i) = beat_info.tR_50_y - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_tR_80(i) = beat_info.tR_80_y - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_tR_20_inv(i) = beat_info.tR_20_y_inv - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_tR_50_inv(i) = beat_info.tR_50_y_inv - getVal(beat,beat_info.min_amplitude_one);
        amp_baselined_tR_80_inv(i) = beat_info.tR_80_y_inv - getVal(beat,beat_info.min_amplitude_one);

        % R ints
        if pairs(i,2) ~= -1
            r = ((pairs(i,1)/ppgSamplingRate) * 1000) - ((pairs(i,2)/ecgSamplingRate) * 1000);
            msec_R_2_first_valley(i) = beat_time_msec(beat_info.min_amplitude_one)*1000 + r;
            msec_R_2_foot(i) = (beat_info.beat_foot/ppgSamplingRate)*1000 + r;
            msec_R_2_pos_slope(i) = (beat_info.max_positive_slope/ppgSamplingRate)*1000 + r;
            
            msec_R_2_tP_20(i) = (beat_info.tP_20_x/ppgSamplingRate)*1000 + r;
            msec_R_2_tP_50(i) = (beat_info.tP_50_x/ppgSamplingRate)*1000 + r;
            msec_R_2_tP_80(i) = (beat_info.tP_80_x/ppgSamplingRate)*1000 + r;
            msec_R_2_tP_20_inv(i) = (beat_info.tP_20_x_inv/ppgSamplingRate)*1000 + r;
            msec_R_2_tP_50_inv(i) = (beat_info.tP_50_x_inv/ppgSamplingRate)*1000 + r;
            msec_R_2_tP_80_inv(i) = (beat_info.tP_80_x_inv/ppgSamplingRate)*1000 + r;

            msec_R_2_systolic_peak(i) = (beat_info.systolic_peak/ppgSamplingRate)*1000 + r;
            msec_R_2_negslopes_pre_dnotch(i) = (beat_info.max_neg_slope_before_dicrotic_notch/ppgSamplingRate)*1000 + r;
            msec_R_2_dicrotic_notch(i) = (beat_info.dicrotic_notch/ppgSamplingRate)*1000 + r;
            msec_R_2_diastolic_peak(i) = (beat_info.diastolic_peak/ppgSamplingRate)*1000 + r;
            msec_R_2_negslopes_post_dnotch(i) = (beat_info.max_neg_slope_after_dicrotic_notch/ppgSamplingRate)*1000 + r;
            msec_R_2_second_valley(i) = (beat_info.min_amplitude_two/ppgSamplingRate)*1000 + r;
            
            msec_R_2_tR_20(i) = (beat_info.tR_20_x/ppgSamplingRate)*1000 + r;
            msec_R_2_tR_50(i) = (beat_info.tR_50_x/ppgSamplingRate)*1000 + r;
            msec_R_2_tR_80(i) = (beat_info.tR_80_x/ppgSamplingRate)*1000 + r;
            msec_R_2_tR_20_inv(i) = (beat_info.tR_20_x_inv/ppgSamplingRate)*1000 + r;
            msec_R_2_tR_50_inv(i) = (beat_info.tR_50_x_inv/ppgSamplingRate)*1000 + r;
            msec_R_2_tR_80_inv(i) = (beat_info.tR_80_x_inv/ppgSamplingRate)*1000 + r;
        end
    end

    beats_total.sleep_stages = sleep_stages;
    beats_total.area_baselined = area_baselined;
    beats_total.area = area;
    beats_total.amp_delta_systolic = amp_delta_systolic;
    beats_total.abs_amp_foot = abs_amp_foot;
    beats_total.abs_amp_peak = abs_amp_peak;
    
    beats_total.idx_begin = idx_begin;
    beats_total.idx_end = idx_end;

    beats_total.idx_foot = idx_feets;
    
    beats_total.tP_20_x = tP_20_x;
    beats_total.tP_50_x = tP_50_x;
    beats_total.tP_80_x = tP_80_x;
    beats_total.tP_20_x_inv = tP_20_x_inv;
    beats_total.tP_50_x_inv = tP_50_x_inv;
    beats_total.tP_80_x_inv = tP_80_x_inv;
    
    beats_total.idx_pos_slope = idx_pos_slopes;
    beats_total.idx_systolic = idx_systolic_peaks;
    beats_total.idx_neg_slope_b4 = idx_neg_slopes_before_dnotch;
    beats_total.idx_neg_slope_after = idx_neg_slopes_after_dnotch;
    beats_total.idx_diastolic = idx_diastolic_peaks;
    beats_total.idx_dnotch = idx_dicrotic_notches;
    
    beats_total.tR_20_x = tR_20_x;
    beats_total.tR_50_x = tR_50_x;
    beats_total.tR_80_x = tR_80_x;
    beats_total.tR_20_x_inv  = tR_20_x_inv ;
    beats_total.tR_50_x_inv  = tR_50_x_inv ;
    beats_total.tR_80_x_inv  = tR_80_x_inv ;
    
    beats_total.msec_beat_length = msec_beat_length;
    
    beats_total.msec_tP_50_2_first_valley = msec_tP_50_2_first_valley;
    beats_total.msec_tP_50_2_foot = msec_tP_50_2_foot;

    beats_total.msec_tP_50_2_tP_20 = msec_tP_50_2_tP_20;
    beats_total.msec_tP_50_2_tP_80 = msec_tP_50_2_tP_80;
    beats_total.msec_tP_50_2_tP_20_inv = msec_tP_50_2_tP_20_inv;
    beats_total.msec_tP_50_2_tP_80_inv = msec_tP_50_2_tP_80_inv;

    beats_total.msec_tP_50_2_pos_slope = msec_tP_50_2_pos_slope;
    beats_total.msec_tP_50_2_systolic_peak = msec_tP_50_2_systolic_peak;
    beats_total.msec_tP_50_2_negslopes_pre_dnotch = msec_tP_50_2_negslopes_pre_dnotch;
    beats_total.msec_tP_50_2_dicrotic_notch = msec_tP_50_2_dicrotic_notch;
    beats_total.msec_tP_50_2_diastolic_peak = msec_tP_50_2_diastolic_peak;
    beats_total.msec_tP_50_2_negslopes_post_dnotch = msec_tP_50_2_negslopes_post_dnotch;
    beats_total.msec_tP_50_2_second_valley = msec_tP_50_2_second_valley;

    beats_total.msec_tP_50_2_tR_20 = msec_tP_50_2_tR_20;
    beats_total.msec_tP_50_2_tR_50 = msec_tP_50_2_tR_50;
    beats_total.msec_tP_50_2_tR_80 = msec_tP_50_2_tR_80;
    beats_total.msec_tP_50_2_tR_20_inv = msec_tP_50_2_tR_20_inv;
    beats_total.msec_tP_50_2_tR_50_inv = msec_tP_50_2_tR_50_inv;
    beats_total.msec_tP_50_2_tR_80_inv = msec_tP_50_2_tR_80_inv;


    beats_total.msec_total_duration_20 = msec_total_duration_tP_20;
    beats_total.msec_total_duration_50 = msec_total_duration_tP_50;
    beats_total.msec_total_duration_80 = msec_total_duration_tP_80;
    beats_total.msec_total_duration_tR_20 = msec_total_duration_tR_20;
    beats_total.msec_total_duration_tR_50 = msec_total_duration_tR_50;
    beats_total.msec_total_duration_tR_80 = msec_total_duration_tR_80;

    beats_total.msec_R_2_first_valley = msec_R_2_first_valley;
    beats_total.msec_R_2_foot = msec_R_2_foot;
    beats_total.msec_R_2_tP_20 = msec_R_2_tP_20;
    beats_total.msec_R_2_tP_50 = msec_R_2_tP_50;
    beats_total.msec_R_2_tP_80 = msec_R_2_tP_80;
    beats_total.msec_R_2_tP_20_inv = msec_R_2_tP_20_inv;
    beats_total.msec_R_2_tP_50_inv = msec_R_2_tP_50_inv;
    beats_total.msec_R_2_tP_80_inv = msec_R_2_tP_80_inv;

    beats_total.msec_R_2_pos_slope = msec_R_2_pos_slope;
    beats_total.msec_R_2_systolic_peak = msec_R_2_systolic_peak;
    beats_total.msec_R_2_negslopes_pre_dnotch = msec_R_2_negslopes_pre_dnotch;
    beats_total.msec_R_2_dicrotic_notch = msec_R_2_dicrotic_notch;

    beats_total.msec_R_2_tR_20 = msec_R_2_tR_20;
    beats_total.msec_R_2_tR_50 = msec_R_2_tR_50;
    beats_total.msec_R_2_tR_80 = msec_R_2_tR_80;
    beats_total.msec_R_2_tR_20_inv = msec_R_2_tR_20_inv;
    beats_total.msec_R_2_tR_50_inv = msec_R_2_tR_50_inv;
    beats_total.msec_R_2_tR_80_inv = msec_R_2_tR_80_inv;

    beats_total.msec_R_2_diastolic_peak = msec_R_2_diastolic_peak;
    beats_total.msec_R_2_negslopes_post_dnotch = msec_R_2_negslopes_post_dnotch;
    beats_total.msec_R_2_second_valley = msec_R_2_second_valley;
    
    
    beats_total.amp_raw_vallies = amp_raw_vallies;
    beats_total.amp_raw_feets = amp_raw_feets;
    
    beats_total.amp_raw_tP_20 = amp_raw_tP_20;
    beats_total.amp_raw_tP_50 = amp_raw_tP_50;
    beats_total.amp_raw_tP_80 = amp_raw_tP_80;
    beats_total.amp_raw_tP_20_inv = amp_raw_tP_20_inv;
    beats_total.amp_raw_tP_50_inv = amp_raw_tP_50_inv;
    beats_total.amp_raw_tP_80_inv = amp_raw_tP_80_inv;
    
    beats_total.amp_raw_tR_20 = amp_raw_tR_20;
    beats_total.amp_raw_tR_50 = amp_raw_tR_50;
    beats_total.amp_raw_tR_80 = amp_raw_tR_80;
    beats_total.amp_raw_tR_20_inv = amp_raw_tR_20_inv;
    beats_total.amp_raw_tR_50_inv = amp_raw_tR_50_inv;
    beats_total.amp_raw_tR_80_inv = amp_raw_tR_80_inv;
    
    beats_total.amp_raw_pos_slopes = amp_raw_pos_slopes;
    beats_total.amp_raw_systolic_peaks = amp_raw_systolic_peaks;
    beats_total.amp_raw_neg_slopes_pre_dnotch = amp_raw_neg_slopes_pre_dnotch;
    beats_total.amp_raw_dicrotic_notches = amp_raw_dicrotic_notches;
    beats_total.amp_raw_diastolic_peaks = amp_raw_diastolic_peaks;
    beats_total.amp_raw_neg_slopes_after_dnotch = amp_raw_neg_slopes_after_dnotch;
    
    beats_total.amp_baselined_feets = amp_baselined_feets;
    
    beats_total.amp_baselined_tP_20 = amp_baselined_tP_20;
    beats_total.amp_baselined_tP_50 = amp_baselined_tP_50;
    beats_total.amp_baselined_tP_80 = amp_baselined_tP_80;
    beats_total.amp_baselined_tP_20_inv = amp_baselined_tP_20_inv;
    beats_total.amp_baselined_tP_50_inv = amp_baselined_tP_50_inv;
    beats_total.amp_baselined_tP_80_inv = amp_baselined_tP_80_inv;
    
    beats_total.amp_baselined_tR_20 = amp_baselined_tR_20;
    beats_total.amp_baselined_tR_50 = amp_baselined_tR_50;
    beats_total.amp_baselined_tR_80 = amp_baselined_tR_80;
    beats_total.amp_baselined_tR_20_inv = amp_baselined_tR_20_inv;
    beats_total.amp_baselined_tR_50_inv = amp_baselined_tR_50_inv;
    beats_total.amp_baselined_tR_80_inv = amp_baselined_tR_80_inv;
    
    beats_total.amp_baselined_pos_slopes = amp_baselined_pos_slopes;
    beats_total.amp_baselined_systolic_peaks = amp_baselined_systolic_peaks;
    beats_total.amp_baselined_neg_slopes_pre_dnotch = amp_baselined_neg_slopes_pre_dnotch;
    beats_total.amp_baselined_dicrotic_notches = amp_baselined_dicrotic_notches;
    beats_total.amp_baselined_diastolic_peaks = amp_baselined_diastolic_peaks;
    beats_total.amp_baselined_neg_slopes_after_dnotch = amp_baselined_neg_slopes_after_dnotch;
    
    beats_total.proportional_pulse_amp = proportional_pulse_amp;

end

function plt(beat,beat_info)
    close all;
    try
    figure('units','normalized','outerposition',[0 0 1 1]);

%     subplot(2, 1, 1);
%     plot(beat,'.-');
%     hold on;
% 
%     plot(beat_info.beat_foot, beat(beat_info.beat_foot), '.r')
%     text(beat_info.beat_foot, beat(beat_info.beat_foot), 'Foot')
%     plot(beat_info.min_amplitude_one, beat(beat_info.min_amplitude_one), '.r')
%     text(beat_info.min_amplitude_one, beat(beat_info.min_amplitude_one), 'Min Amp 1')
%     plot(beat_info.max_positive_slope, beat(beat_info.max_positive_slope), '.r')
%     text(beat_info.max_positive_slope, beat(beat_info.max_positive_slope), 'Pos Slope')
%     plot(beat_info.systolic_peak, beat(beat_info.systolic_peak), '.r')
%     text(beat_info.systolic_peak, beat(beat_info.systolic_peak), 'Systolic Peak')
%     plot(beat_info.dicrotic_notch, beat(beat_info.dicrotic_notch), '.r')
%     text(beat_info.dicrotic_notch, beat(beat_info.dicrotic_notch), 'Dicrotic Notch')
%     plot(beat_info.diastolic_peak, beat(beat_info.diastolic_peak), '.r')
%     text(beat_info.diastolic_peak, beat(beat_info.diastolic_peak), 'Diastolic Peak')
%     plot(beat_info.max_neg_slope_before_dicrotic_notch, beat(beat_info.max_neg_slope_before_dicrotic_notch), '.r')
%     text(beat_info.max_neg_slope_before_dicrotic_notch, beat(beat_info.max_neg_slope_before_dicrotic_notch), 'Slope before Dnotch')
%     plot(beat_info.min_amplitude_two, beat(beat_info.min_amplitude_two), '.r')
%     text(beat_info.min_amplitude_two, beat(beat_info.min_amplitude_two), 'Min Amp 2')
%     plot(beat_info.max_neg_slope_after_dicrotic_notch, beat(beat_info.max_neg_slope_after_dicrotic_notch), '.r')
%     text(beat_info.max_neg_slope_after_dicrotic_notch, beat(beat_info.max_neg_slope_after_dicrotic_notch), 'Slope after Dnotch')
% 
%     hline(beat(beat_info.tP_20_x),{'Color','r'},'tP 20',.5, {}, gca, true,[beat_info.tP_20_x beat_info.tP_20_x_inv]);
%     hline(beat(beat_info.tP_50_x),{'Color','r'},'tP 50',.5, {}, gca, true,[beat_info.tP_50_x beat_info.tP_50_x_inv]);
%     hline(beat(beat_info.tP_80_x),{'Color','r'},'tP 80',.5, {}, gca, true,[beat_info.tP_80_x beat_info.tP_80_x_inv]);
% 
%     plot(beat_info.tP_20_x, beat(beat_info.tP_20_x), '.r')
%     text(beat_info.tP_20_x, beat(beat_info.tP_20_x), 'tP 20')
%     plot(beat_info.tP_50_x, beat(beat_info.tP_50_x), '.r')
%     text(beat_info.tP_50_x, beat(beat_info.tP_50_x), 'tP 50')
%     plot(beat_info.tP_80_x, beat(beat_info.tP_80_x), '.r')
%     text(beat_info.tP_80_x, beat(beat_info.tP_80_x), 'tP 80')
% 
%     plot(beat_info.tP_20_x_inv, beat(beat_info.tP_20_x_inv), '.r')
%     text(beat_info.tP_20_x_inv, beat(beat_info.tP_20_x_inv), 'tP^-^1 20')
%     plot(beat_info.tP_50_x_inv, beat(beat_info.tP_50_x_inv), '.r')
%     text(beat_info.tP_50_x_inv, beat(beat_info.tP_50_x_inv), 'tP^-^1 50')
%     plot(beat_info.tP_80_x_inv, beat(beat_info.tP_80_x_inv), '.r')
%     text(beat_info.tP_80_x_inv, beat(beat_info.tP_80_x_inv), 'tP^-^1 80')
%     
%     hline(beat(beat_info.tR_20_x),{'Color','m'},'tR 20',.5, {}, gca, true,[beat_info.tR_20_x_inv beat_info.tR_20_x]);
%     hline(beat(beat_info.tR_50_x),{'Color','m'},'tR 50',.5, {}, gca, true,[beat_info.tR_50_x_inv beat_info.tR_50_x]);
%     hline(beat(beat_info.tR_80_x),{'Color','m'},'tR 80',.5, {}, gca, true,[beat_info.tR_80_x_inv beat_info.tR_80_x]);
% 
%     plot(beat_info.tR_20_x, beat(beat_info.tR_20_x), '.r')
%     text(beat_info.tR_20_x, beat(beat_info.tR_20_x), 'tR 20')
%     plot(beat_info.tR_50_x, beat(beat_info.tR_50_x), '.r')
%     text(beat_info.tR_50_x, beat(beat_info.tR_50_x), 'tR 50')
%     plot(beat_info.tR_80_x, beat(beat_info.tR_80_x), '.r')
%     text(beat_info.tR_80_x, beat(beat_info.tR_80_x), 'tR 80')
%     
%     plot(beat_info.tR_20_x_inv, beat(beat_info.tR_20_x_inv), '.r')
%     text(beat_info.tR_20_x_inv, beat(beat_info.tR_20_x_inv), 'tR^-^1 20')
%     plot(beat_info.tR_50_x_inv, beat(beat_info.tR_50_x_inv), '.r')
%     text(beat_info.tR_50_x_inv, beat(beat_info.tR_50_x_inv), 'tR^-^1 50')
%     plot(beat_info.tR_80_x_inv, beat(beat_info.tR_80_x_inv), '.r')
%     text(beat_info.tR_80_x_inv, beat(beat_info.tR_80_x_inv), 'tR^-^1 80')
%     
% 
%     vline(beat_info.systolic_peak,'g--');   
% 
%     subplot(2, 1, 2);
    plot(beat_info.beat_time_msec, beat,'.-');
    hold on;

    plot(beat_info.beat_time_msec(beat_info.beat_foot), beat(beat_info.beat_foot ), '.r')
    text(beat_info.beat_time_msec(beat_info.beat_foot), beat(beat_info.beat_foot ), 'Foot')
    plot(beat_info.beat_time_msec(beat_info.min_amplitude_one), beat(beat_info.min_amplitude_one), '.r')
    text(beat_info.beat_time_msec(beat_info.min_amplitude_one), beat(beat_info.min_amplitude_one), 'Min Amp 1')
    plot(beat_info.beat_time_msec(beat_info.max_positive_slope), beat(beat_info.max_positive_slope), '.r')
    text(beat_info.beat_time_msec(beat_info.max_positive_slope), beat(beat_info.max_positive_slope), 'Pos Slope')
    plot(beat_info.beat_time_msec(beat_info.systolic_peak), beat(beat_info.systolic_peak), '.r')
    text(beat_info.beat_time_msec(beat_info.systolic_peak), beat(beat_info.systolic_peak), 'Systolic Peak')
    plot(beat_info.beat_time_msec(beat_info.dicrotic_notch), beat(beat_info.dicrotic_notch), '.r')
    text(beat_info.beat_time_msec(beat_info.dicrotic_notch), beat(beat_info.dicrotic_notch), 'Dicrotic Notch')
    plot(beat_info.beat_time_msec(beat_info.diastolic_peak), beat(beat_info.diastolic_peak), '.r')
    text(beat_info.beat_time_msec(beat_info.diastolic_peak), beat(beat_info.diastolic_peak), 'Diastolic Peak')
    plot(beat_info.beat_time_msec(beat_info.max_neg_slope_before_dicrotic_notch), beat(beat_info.max_neg_slope_before_dicrotic_notch), '.r')
    text(beat_info.beat_time_msec(beat_info.max_neg_slope_before_dicrotic_notch), beat(beat_info.max_neg_slope_before_dicrotic_notch), 'Slope before Dnotch')
    plot(beat_info.beat_time_msec(beat_info.min_amplitude_two), beat(beat_info.min_amplitude_two), '.r')
    text(beat_info.beat_time_msec(beat_info.min_amplitude_two), beat(beat_info.min_amplitude_two), 'Min Amp 2')
    plot(beat_info.beat_time_msec(beat_info.max_neg_slope_after_dicrotic_notch), beat(beat_info.max_neg_slope_after_dicrotic_notch), '.r')
    text(beat_info.beat_time_msec(beat_info.max_neg_slope_after_dicrotic_notch), beat(beat_info.max_neg_slope_after_dicrotic_notch), 'Slope after Dnotch')

    hline(beat_info.tP_20_y,{'Color','r'},'tP 20',.5, {}, gca, true,[beat_info.tP_20_x/beat_info.sampling_rate*1000 beat_info.tP_20_x_inv/beat_info.sampling_rate*1000]);
    hline(beat_info.tP_50_y,{'Color','r'},'tP 50',.5, {}, gca, true,[beat_info.tP_50_x/beat_info.sampling_rate*1000 beat_info.tP_50_x_inv/beat_info.sampling_rate*1000]);
    hline(beat_info.tP_80_y,{'Color','r'},'tP 80',.5, {}, gca, true,[beat_info.tP_80_x/beat_info.sampling_rate*1000 beat_info.tP_80_x_inv/beat_info.sampling_rate*1000]);

    plot(beat_info.tP_20_x/beat_info.sampling_rate*1000, beat_info.tP_20_y, '.r')
    text(beat_info.tP_20_x/beat_info.sampling_rate*1000, beat_info.tP_20_y, 'tP 20')
    plot(beat_info.tP_50_x/beat_info.sampling_rate*1000, beat_info.tP_50_y, '.r')
    text(beat_info.tP_50_x/beat_info.sampling_rate*1000, beat_info.tP_50_y, 'tP 50')
    plot(beat_info.tP_80_x/beat_info.sampling_rate*1000, beat_info.tP_80_y, '.r')
    text(beat_info.tP_80_x/beat_info.sampling_rate*1000, beat_info.tP_80_y, 'tP 80')

    plot(beat_info.tP_20_x_inv/beat_info.sampling_rate*1000, beat_info.tP_20_y_inv, '.r')
    text(beat_info.tP_20_x_inv/beat_info.sampling_rate*1000, beat_info.tP_20_y_inv, 'tP^-^1 20')
    plot(beat_info.tP_50_x_inv/beat_info.sampling_rate*1000, beat_info.tP_50_y_inv, '.r')
    text(beat_info.tP_50_x_inv/beat_info.sampling_rate*1000, beat_info.tP_50_y_inv, 'tP^-^1 50')
    plot(beat_info.tP_80_x_inv/beat_info.sampling_rate*1000, beat_info.tP_80_y_inv, '.r')
    text(beat_info.tP_80_x_inv/beat_info.sampling_rate*1000, beat_info.tP_80_y_inv, 'tP^-^1 80')

    hline(beat_info.tR_20_y,{'Color','m'},'tR 20',.5, {}, gca, true,[beat_info.tR_20_x_inv/beat_info.sampling_rate*1000 beat_info.tR_20_x/beat_info.sampling_rate*1000]);
    hline(beat_info.tR_50_y,{'Color','m'},'tR 50',.5, {}, gca, true,[beat_info.tR_50_x_inv/beat_info.sampling_rate*1000 beat_info.tR_50_x/beat_info.sampling_rate*1000]);
    hline(beat_info.tR_80_y,{'Color','m'},'tR 80',.5, {}, gca, true,[beat_info.tR_80_x_inv/beat_info.sampling_rate*1000 beat_info.tR_80_x/beat_info.sampling_rate*1000]);

    plot(beat_info.tR_20_x/beat_info.sampling_rate*1000, beat_info.tR_20_y, '.r')
    text(beat_info.tR_20_x/beat_info.sampling_rate*1000, beat_info.tR_20_y, 'tR 20')
    plot(beat_info.tR_50_x/beat_info.sampling_rate*1000, beat_info.tR_50_y, '.r')
    text(beat_info.tR_50_x/beat_info.sampling_rate*1000, beat_info.tR_50_y, 'tR 50')
    plot(beat_info.tR_80_x/beat_info.sampling_rate*1000, beat_info.tR_80_y, '.r')
    text(beat_info.tR_80_x/beat_info.sampling_rate*1000, beat_info.tR_80_y, 'tR 80')

    plot(beat_info.tR_20_x_inv/beat_info.sampling_rate*1000, beat_info.tR_20_y_inv, '.r')
    text(beat_info.tR_20_x_inv/beat_info.sampling_rate*1000, beat_info.tR_20_y_inv, 'tR^-^1 20')
    plot(beat_info.tR_50_x_inv/beat_info.sampling_rate*1000, beat_info.tR_50_y_inv, '.r')
    text(beat_info.tR_50_x_inv/beat_info.sampling_rate*1000, beat_info.tR_50_y_inv, 'tR^-^1 50')
    plot(beat_info.tR_80_x_inv/beat_info.sampling_rate*1000, beat_info.tR_80_y_inv, '.r')
    text(beat_info.tR_80_x_inv/beat_info.sampling_rate*1000, beat_info.tR_80_y_inv, 'tR^-^1 80')

    vline(beat_info.beat_time_msec(beat_info.systolic_peak),'g--');  
    title('PPG beat');
    xlabel('Time (msec)');
    ylabel('mV');
    catch
    end
end

function c = setVal(arr,a,b)
    try
        c = arr(a) - b;
    catch 
        c = nan;
    end
end

function c = getVal(arr,x)
    try
        c = arr(x);
    catch 
        c = nan;
    end
end
