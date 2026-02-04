function [value_array] = RoundToClosestBin(bins,value_array)
    processed = false(numel(value_array),1);
    for i= 1:numel(bins)
        value_array(value_array<=bins(i) & ~processed) = i;
        
        processed = value_array <= i; 
    end