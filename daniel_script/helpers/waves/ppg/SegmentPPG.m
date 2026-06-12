function [ppgMinAmps, maxAmps] = SegmentPPG(ppg, sampleRate)

    %agressive smooth and tight moving mean
    ppg_smooth = nanfastsmooth(ppg, sampleRate * .25, 3);
    M = movmean(ppg_smooth, sampleRate);

    plow = ppg;
    plow(ppg_smooth > M) = nan;
    p_mask(~isnan(plow)) = 0;
    p_mask(isnan(plow)) = 1;
    [peakidx, vallyidx] = segBeats(ppg, p_mask);

    if abs(length(peakidx) - length(vallyidx)) > 1
        error('not expecting this')
    end

    ppg_outliers_peaks = stdoutlier(peakidx', 2.5, 100, 'both', 0);
    vallyoutlier_time = stdoutlier(vallyidx', 2.5, 100, 'both', 0);
    %vallyoutlier_amp = stdoutlier(ppg(vallyidx), 2.5, 100,'both', 0);

    if vallyidx(1) < peakidx(1)%vally peak vally peak ....
            pidx = 1;
    else % peak valley peak valley ....
            pidx = 2;
    end

    %
    %     close all;
    %     figure('units','normalized','outerposition',[0 0 1 1])
    %     plot(ppg,'Color',rgb('grey'));
    %     hold on;
    %     plot(ppg_smooth,':','Color',rgb('white'))
    %     plot(M,'Color',rgb('yellow'));
    %     plot(plow,'k--');
    %
    %     plot(peakidx, ppg(peakidx), '*w')
    %     plot(vallyidx, ppg(vallyidx), 'ow')
    %
    %     set(gca,'Color','k')

    for vidx = 1:length(vallyidx) - 1

        if vidx == 1 || vidx == length(vallyidx) - 1
            pidx = pidx + 1;
            continue
        end

        if vallyoutlier_time(vidx)
            %             v1 = vline(peakidx(pidx),'r');
            %             v2 = vline(vallyidx(vidx),'g');
            %             xlim([vallyidx(vidx)-1000,vallyidx(vidx)+1000]);
            vallyidx(vidx) = newVallyFromPeaks(ppg, vidx, pidx, peakidx, vallyidx, ppg_outliers_peaks, vallyoutlier_time);
            %             delete(v1);
            %             delete(v2);
        end

        pidx = pidx + 1;

    end

    %

    %   debug plot

    % %         vline(vallyidx(ppg_outliers_peaks),'--g');
    % %         vline(vallyidx(ppg_outliers_vallies),'--c');
    %
    %     h = [];
    %     h(end+1) = plot(nan,':','Color',rgb('white'));
    %     h(end+1) = plot(nan,'Color',rgb('yellow'));
    %
    %     legend(h,{'\color{white} Smooth PPG','\color{white} Moving Mean'});
    %     minOutliers = stdoutlier(vallyidx', 2.5, 10, 0);
    %     plot(vallyidx(minOutliers), ppg(vallyidx(minOutliers)), 'or')
    %
    %
    %
    %
    %     [B, N, BI] = RunLength(minOutliers);
    %     if ~isempty(N(N>3 & B == 1))
    %         idxs = find((N>3 & B == 1) == 1);
    %         for q = 1:length(idxs)
    %             idx1 = BI(idxs(q));
    %             idx2 = BI(idxs(q))+N(idxs(q))-1;
    %             section = vallyidx(idx1:idx2);
    %             out = PpgBestFitOfBeatTrain(vallyidx(1:idx1-1), vallyidx(idx2+1:end), section, 15);
    %             vallyidx = [vallyidx(1:idx1-1) out vallyidx(idx2+1:end)];
    %         end
    %     end

    ppgMinAmps = vallyidx;
    maxAmps = peakidx;
end

function [newidx] = newVallyFromPeaks(ppg, currVal, currPeak, peaks_idxs, vally_idxs, peakoutliers, vallyoutlier)
    beat = ppg(peaks_idxs(currPeak - 1):peaks_idxs(currPeak));
    [~, idx] = min(beat);

    if ~peakoutliers(currPeak - 1) &&~peakoutliers(currPeak)% no weird peaks

        if peaks_idxs(currPeak - 1) + idx - 1 == vally_idxs(currVal)
            newidx = vally_idxs(currVal);
        else
            newidx = peaks_idxs(currPeak - 1) + idx - 1;
        end

    else
        newidx = vally_idxs(currVal);
    end

end

function [dif] = oscillatingDiff(a, b)
    % diff skipping 1 each time.
    % assumtion is in time the two vectors would go [a,b,a,b,a,b,a,b...]
    % function gives back each ab pair diff, aka the diff on each odd index
    all = [a; b];
    all = sort(all);

    d = diff(all);
    dif = d(1:2:end);
end

function [peaks, vallies] = segBeats(ppg, mask)
    % find where 0 vs 1 and length of peroid
    [B(:, 1), B(:, 2), B(:, 3)] = RunLength(mask);

    % elements above and below line
    high = zeros(sum(B(:, 1)), 2);
    low = zeros(size(B, 1) - sum(B(:, 1)), 2);
    highidx = 1;
    lowidx = 1;

    for i = 1:size(B, 1)

        if B(i, 1) == 0
            low(lowidx, 1) = B(i, 3);
            low(lowidx, 2) = B(i, 3) + B(i, 2) - 1;
            lowidx = lowidx + 1;
        else
            high(highidx, 1) = B(i, 3);
            high(highidx, 2) = B(i, 3) + B(i, 2) - 1;
            highidx = highidx + 1;
        end

    end

    % find max
    peaks = zeros(size(high, 1), 1);

    for i = 1:size(high, 1)
        % max idx is indexed to begining of section (high(i,1)) so add that
        % to temp to move to proper position.
        [~, tmp] = max(ppg(high(i, 1):high(i, 2)));
        peaks(i) = tmp + high(i, 1);
    end

    peaks = peaks - 1;

    vallies = zeros(size(low, 1), 1);

    for i = 1:size(low, 1)
        % max idx is indexed to begining of section (high(i,1)) so add that
        % to temp to move to proper position.
        [~, tmp] = min(ppg(low(i, 1):low(i, 2)));
        vallies(i) = tmp + low(i, 1);
    end

    vallies = vallies - 1;

end

function [out] = PpgBestFitOfBeatTrain(sectionBefore, sectionAfter, dataToMatch, beats)

    try
        sectionBefore = sectionBefore(end - beats + 1:end);
    catch
    end

    try
        sectionAfter = sectionAfter(1:beats);
    catch
    end

    out = [];

    for i = 1:length(dataToMatch)
        medBefore = round(mean(diff(sectionBefore)));
        medAfter = round(mean(diff(sectionAfter)));

        endBefore = sectionBefore(end);
        beginAfter = sectionAfter(1);

        beforeRate = endBefore:medBefore:beginAfter;
        afterRate = endBefore:medAfter:beginAfter;

        if length(beforeRate) > length(afterRate)
            midLengths = round((beforeRate(1:numel(afterRate)) + afterRate) / 2);
        else
            midLengths = round((afterRate(1:numel(beforeRate)) + beforeRate) / 2);
        end

        v1 = vline(midLengths, 'r');

        if length(midLengths) < 2
            mid = endBefore + medBefore;
        else
            mid = midLengths(2);
        end

        errors = zeros(length(dataToMatch), 1);

        for z = 1:length(dataToMatch)
            errors(z) = abs(mid - dataToMatch(z));
        end

        [val, minIdx] = min(errors);
        out(end + 1) = dataToMatch(minIdx);

        sectionBefore(end + 1) = dataToMatch(minIdx);
        dataToMatch(dataToMatch <= minIdx) = [];

        delete(v1);
    end

end
