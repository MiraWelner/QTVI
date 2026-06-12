function [pairs] = pairRtoPPGBeat(ecg, ppg, ecgSamplingRate, ppgSamplingRate, ecgRIndex, ppgMinAmps)
    ecgtime = (0:length(ecg) - 1) / ecgSamplingRate / 60;
    ppgtime = (0:length(ppg) - 1) / ppgSamplingRate / 60;

    % ppg vally idx, ecgR idx, outliers
    pairs = nan(length(ecgRIndex), 2);
    pairs(:, 2) = ecgRIndex;
    % pair every r to a ppg idx or multiple ppg idxs
    for i = 1:length(ecgRIndex)

        if i == 1
            beginidx = 1;
            endidx = ecgRIndex(i + 1);

        elseif i == length(ecgRIndex)
            beginidx = ecgRIndex(i - 1);
            endidx = length(ecgtime);
        else
            beginidx = ecgRIndex(i - 1);
            endidx = ecgRIndex(i + 1);
        end

        begtime = ecgtime(beginidx);
        enbtime = ecgtime(endidx);

        possible_ppg_parings = find(begtime <= ppgtime(ppgMinAmps) & enbtime >= ppgtime(ppgMinAmps));

        errors = zeros(length(possible_ppg_parings), 1);

        for x = 1:length(possible_ppg_parings)
            errors(x) = abs(ppgtime(ppgMinAmps(possible_ppg_parings(x))) - ecgtime(ecgRIndex(i)));
        end

        [~, minIdx] = min(errors);

        if isempty(minIdx)
            pairs(i, 1) = -1;
        else
            pairs(i, 1) = ppgMinAmps(possible_ppg_parings(minIdx));
        end

        %         v= vline(ecgRIndex(i));
        %         s = scatter(ppgMinAmps(possible_ppg_parings),ppg(ppgMinAmps(possible_ppg_parings)),'r');
        %
        %         if pairs(i,1) ~= -1
        %             plot([pairs(i,1),pairs(i,2)],[ppg(pairs(i,1)),ecg(pairs(i,2))],'c')
        %         end
        %         xlim([ecgRIndex(i)-1000 ecgRIndex(i)+1000]);
        %         delete(s);
        %         delete(v);

    end

    % add the ppg vallies that were not paired to rs
    unusedset = setdiff(ppgMinAmps, pairs(:, 1));

    for i = 1:length(unusedset)
        pairs(end + 1, 1) = unusedset(i);
        pairs(end, 2) = -1;
    end

    % after sorting if there is a -1 in pairs column 1 a r existed that couldn't pair
    % with a ppg,  if theres a -1 in column 2 a ppg existed that was not
    % paired with an r
    pairs = sortrows(pairs, 1);
    [B, N, BI] = RunLength(pairs(:, 1));

    % correct r's assigned to more then 1 ppg
    i = 1;
    while i <= length(N)

        if N(i) == 2 && B(i) ~= -1

            idxs = BI(i):BI(i) + N(i) - 1;
            p = pairs(idxs, :);
            
            if idxs(end) == length(pairs) % if last?
                pairs(idxs(end), 1) = -1;
                i = i + 1;
                continue
            elseif pairs(idxs(end) + 1, 2) == -1 % if next r is unpaired
                pairs(idxs(end) + 1, 2) = pairs(idxs(end), 2);
                pairs(idxs(end), :) = [];
                [B, N, BI] = RunLength(pairs(:, 1));
            else

                errors = zeros(length(p), 1);

                for x = 1:size(p,1)
                    errors(x) = abs(ppgtime(p(x, 1)) - ecgtime(p(x, 2)));
                end

                [~, minIdx] = min(errors);

                if minIdx == 1
                    begidx = p(end, 1);
                    if length(pairs) >= idxs(end) + 1; endidx = pairs(idxs(end) + 1, 1); else; endidx = length(ppg); end

                else
                    if idxs(1) - 1 > 1 && pairs(idxs(1) - 1, 1) > 1
                        begidx = pairs(idxs(1) - 1, 1); 
                    else
                        begidx = 1; 
                    end
                    endidx = p(end, 1);
                end

                ppgSeg = ppg(begidx:endidx);
                smoothed = nanfastsmooth(ppgSeg, length(ppgSeg) / 10, 3);
                [a, b] = findpeaks(smoothed);

                max_2 = maxk(a, 2);

                if length(max_2) > 1
                    [~, newp] = min(ppgSeg(b(1):b(2)));
                    newp = newp + b(1) - 1;
                else % nothing found, ecg had noise in it

                    if minIdx == 1
                        pairs(idxs(2), :) = pairs(idxs(2), 2);
                        [B, N, BI] = RunLength(pairs(:, 1));
                        i = i+1;
                        continue;
                    else
                        pairs(idxs(1), 1) = pairs(idxs(1), 2);
                        [B, N, BI] = RunLength(pairs(:, 1));
                        i = i+1;
                        continue;
                    end

                end

                if minIdx == 1
                    pairs(idxs(2), 1) = pairs(idxs(1), 1) + newp;
                    [B, N, BI] = RunLength(pairs(:, 1));
                else
                    pairs(idxs(1), 1) = begidx + newp - 1;
                    [B, N, BI] = RunLength(pairs(:, 1));
                end
            end

        elseif N(i) > 2 && B(i) ~= -1 % multiple rs for one ppg probably bad ppg.  update based on r
            idxs = BI(i):BI(i) + N(i) - 1;
            p = pairs(idxs, :);
            update = 0;
            for x = 1:length(p)
                if p(1,1) == p(x, 2) % if exact match between ppg and ecg
                    update = 1;
                    
                    for z = 1:length(p)
                        if z ~= x
                            pairs(idxs(z), 1) = pairs(idxs(z), 2);               
                        end
                    end
                    [B, N, BI] = RunLength(pairs(:, 1));
                    break;
                end
            end
            if update == 0
                error('implment me');
            end
        end

        i = i + 1;

        if i > length(N)
            break;
        end

    end

    
    %% final clean
    % not 100% sure what to do when r is not paired atm.  for now just
    % propagate the time they exist to the ppg i guess...
    pairs(pairs(:,1)==-1,1) = pairs(pairs(:,1)==-1,2);
    pairs = sortrows(pairs, 1);

    %plt(pairs, ecgRIndex,ppg, ppgMinAmps, ecg)
    %i = 1;
%     while true
%         
%         
%         low_outliers = stdoutlier(pairs(:,1), 3.5, 50, 'lower',0);
%         high_outliers = stdoutlier(pairs(:,1), 4, 50, 'upper',0);
% 
%         
%         
%         if i > length(pairs)-1
%             break;
%         end 
%     end
    

end

function plt(pairs, ecgRIndex,ppg, ppgMinAmps, ecg)
    close all;

    figure('units', 'normalized', 'outerposition', [0 0 1 1])

    plot(ecg, ':', 'Color', rgb('dimgrey'));
    hold on;
    plot(ppg, 'Color', rgb('white'));

    scatter(ecgRIndex, ecg(ecgRIndex), '*w');
    scatter(ppgMinAmps, ppg(ppgMinAmps), 'w');

    set(gca, 'Color', 'k')
    grid;

    for x = 1:size(pairs, 1)

        if pairs(x, 1) ~= -1 && pairs(x, 2) ~= -1
            plot([pairs(x, 1), pairs(x, 2)], [ppg(pairs(x, 1)), ecg(pairs(x, 2))], 'c')
        end

    end

end
