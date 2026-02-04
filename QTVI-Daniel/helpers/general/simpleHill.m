
function [peak_index] = simpleHill(beat)
    len = length(beat);
    
    [start_val,start_idx] = min(beat(1:round(len*.2)));
    peak_value = start_val;
    peak_index = start_idx;


    %% find peak
    % search forward till max or end
    overshoot = 0;
    for i = start_idx:len
        if peak_value < beat(i)
            peak_index = i;
            peak_value = beat(i);
        else
            overshoot = overshoot +1;
        end
        if overshoot > 3
            break;
        end
    end
        
end