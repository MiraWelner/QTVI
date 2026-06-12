function [idxs] = ExpandIndexs(rr_indexes, subtractionleft,additionright,max_size)
    idxs = [rr_indexes(1:end - 1)' rr_indexes(2:end)'];

    idxs(:,1) = round(idxs(:,1) - subtractionleft);
    idxs(:,2) = round(idxs(:,2) + additionright);
    
    idxs(idxs < 1) = 1;
    idxs(idxs > max_size) = max_size;
end

