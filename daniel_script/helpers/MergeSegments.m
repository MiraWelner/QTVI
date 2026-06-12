function final = MergeSegments(segments)        
    final = sortrows(segments,1);
    index=1;
    while 1
        if index > size(final,1)
           break 
        end

        begidx = final(index,1);
        endidx = final(index,2);
        mask = (begidx <= final) & (final <= endidx);

        idxs = find(sum(mask,2)==2); % sections btw intervals
        idxs = idxs(idxs>index);
        for x = 1:length(idxs)                  
           final(idxs(x),:) = []; % remove section
           mask(idxs(x),:) = []; % remove section
           idxs = idxs-1;

        end

        flag = 0;
        idxs = find(mask(:,1)==1); % sections overlapping above intervals
        idxs = idxs(idxs>index);
        for x = 1:length(idxs) 
            if final(idxs(x),2) > endidx
                final(index,2) = final(idxs(x),2); % replace current end w/ new end
            end
            final(idxs(x),:) = [];
            idxs = idxs-1;

            flag = 1;
        end 

        if flag == 0
            index = index+1;                
        end

    end