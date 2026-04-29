function clusters = OneDDbscan(x,eps)
    if eps <1
        eps =1;
    end
    prev=0;
    s = sum(x);
    if s == 0
        clusters = [];
        return
    elseif s == 1
        z = find(x==1);
        clusters = [z z];
        return
    end
    
    clusters = [];
    l = length(x);
    for i = 1:l
        if x(i) == 1
            if i == l
                if prev ~= 0
                    clusters(end,2) = i;
                else
                    clusters(end+1,1) = i;
                    clusters(end,2) = i;

                end
            else
                if prev==0
                    clusters(end+1,1) = i;
                    prev = i;
                else
                    if i-prev <= eps
                        prev = i;
                    end
                end
            end
            
        elseif i == l
            if prev ~= 0
                clusters(end,2) = prev;
            end
        else
            if prev ~= 0
                if i-prev >= eps
                    clusters(end,2) = prev;
                    prev = 0;
                end
            end
        end
    end



end