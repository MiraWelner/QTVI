function [data] = splitOverLappingBins(indices, bins, breaks)
%splitOverLappingBins: Given the indices, pick up those that have start and
%end different. Use the breaks to split the indices so that the resulting
%indices all lie in the same bin.
%format: begin_idx, end_idx, begin_bin_idx, end_bin_idx
data = [indices, bins];
diff_idx = find(data(:, 3) ~= data(:, 4));
td = data(diff_idx, :);
for row = 1:size(td, 1)
    data_idx = diff_idx(row);
    bin_end = td(row, 2);
    last_bin_num = td(row, 4);
    for j = (td(row, 3): 1: td(row, 4))
        if j == td(row, 3)
            data(data_idx, :) = [td(row, 1) breaks(j) j j];                  
        elseif j == last_bin_num
            data(end + 1, :) = [breaks(j - 1) bin_end j j];
        else
            data(end + 1, :) = [breaks(j - 1) breaks(j) j j];
        end
    end
end
[~, order] = sort(data(:, 1));
data = data(order, 1 : 3);
