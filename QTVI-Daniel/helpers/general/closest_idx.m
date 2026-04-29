function [idx] = closest_idx(testArr, val)
    tmp = abs(testArr - val);
    [~, idx] = min(tmp);
    %val = testArr(idx);
