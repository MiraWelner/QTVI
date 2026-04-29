function [values] = PPGVallies(ppg, ppgSamplingRate,precentages)
    ppg_len = length(ppg);
    
    [start_val,start_idx] = min(ppg(1:ppgSamplingRate*2)); % min value in first 2 secs (arbitrary chosen)
    
%     close all;
%     plot(ppg); hold on;
%     set(gcf, 'Position', get(0, 'Screensize'));
    
    values = [];
    while start_idx < length(ppg)
        begin_index = start_idx;
        
        peak_value = -inf;
        peak_index = begin_index;
        end_index = begin_index;
        end_value = inf;
        
        values = [values start_idx];
        
%         v = vline(start_idx,{'--m','LineWidth',2});
%         p = [];
        
        %% find peak
        % search forward till max or end
        overshoot = 0;
        for i = begin_index + 1:ppg_len
            if peak_value < ppg(i)
                peak_index = i;
                peak_value = ppg(i);
%                 p(end+1) = vline(peak_index,'g');
            else
                overshoot = overshoot +1;
%                 p(end+1) = vline(i,'r');
            end
            if overshoot > round(ppgSamplingRate*.3)
%                 p(end+1) = vline(i,'r');
                break;
            end
        end

%        xlim([begin_index-300,begin_index+300])
        
        %% find valley
        overshoot = 0;
%         delete(p);
%         delete(v);
%         p = [];

%         v = vline(peak_index,{'--m','LineWidth',2});

        for i =  peak_index + 1:ppg_len
            if end_value >= ppg(i)
                end_index = i;
                end_value = ppg(i);
%                 p(end+1) = vline(end_index,'g');
            else
                overshoot = overshoot +1;
%                 p(end+1) = vline(i,'r');
            end
            if overshoot > round(ppgSamplingRate*.3)
%                 p(end+1) = vline(i,'r');
                break;
            end
        end
%         xlim([peak_index-300,peak_index+300])
        
        if peak_index == ppg_len
            start_idx = ppg_len;
        else
            start_idx = end_index;
        end
%         delete(p);
%         delete(v);
        
    end


%     if dbg_plot == 1
%         fig = figure('Name', 'PPG Valleys','visible','on');
%         plot(time,ppg);
%         hold on; 
%         plot(time(values),ppg(values),'o');
%         hold off;
%         ShowDbgPlots({fig});
%     end

end
