function [Regions] = findGoodData(verbose, Limit, patchSignal, filterLevel, review, varargin)
    % author: Michael Hautman
    % Identify largest contiguous area of usable data across multiple sets for useful
    % comparison
    % Verbose : boolean, displays all sets and where the Identified intervals
    % reside on those sets
    % Limit: integer: target length of intervals identified, in minutes
    % if set to 0, returns all good intervals regardless of length
    % patchSignal: boolean, if enabled, ignores small patches of noise
    % filterLevel : check IDgd for explanation of the filterLevel argument.
    % options 'strict' 'normal' 'loose', or the corresponding 2, 1, 0
    %
    % varargin inputs:
    % provide first the dataset and then it's sampling rate.
    % example:
    % findGoodData(1, 5, dataset1, samprate1, dataset2, samprate2, ... datasetN, samprateN)

    %see IDgd implementation for explanation of filter level

    numInputs = length(varargin);

    if (numInputs < 4 && mod(numInputs, 2) ~= 0)
        error("not enough inputs")
    end

    accepted = 0;

    while (~accepted)
        numSigs = numInputs / 2;

        for I = 1:numSigs
            [goodIntervals{I}] = IDgd(varargin{I * 2 - 1}, varargin{I * 2}, 0, Limit, filterLevel, patchSignal);
            goodTimes{I} = goodIntervals{I} / varargin{I * 2};
        end

        if (~review)
            accepted = 1;
        end

        %% Identifying regions of overlap between ranges identified as having good data
        overlapRegions = [];
        totalRegions = 1;

        template = goodTimes{1};

        for I = 1:numSigs - 1

            if (~isempty(overlapRegions))
                template = overlapRegions;
            end

            for d = 1:size(template, 1)
                tempStart = template(d, 1);
                tempEnd = template(d, 2);
                compare = goodTimes{I + 1};

                for c = 1:size(compare, 1)
                    compareStart = compare(c, 1);
                    compareEnd = compare(c, 2);

                    if (tempStart <= compareStart && ...
                            compareStart <= tempEnd)
                        overlapRegions(totalRegions, 1) = compareStart;

                        if (tempEnd <= compareEnd)
                            overlapRegions(totalRegions, 2) = tempEnd;
                        else
                            overlapRegions(totalRegions, 2) = compareEnd;
                        end

                        totalRegions = totalRegions + 1;
                    elseif (tempStart >= compareStart && ...
                            tempStart <= compareEnd)
                        overlapRegions(totalRegions, 1) = tempStart;

                        if (tempEnd <= compareEnd)
                            overlapRegions(totalRegions, 2) = tempEnd;
                        else
                            overlapRegions(totalRegions, 2) = compareEnd;
                        end

                        totalRegions = totalRegions + 1;
                    end

                end

            end

        end

        if (numSigs < 2)
            overlapRegions = goodTimes{1};
        end

        %% eliminating regions shorter than given limit
        badRegions = [];
        i = 1;
        CLimit = 2;

        if (Limit > 0)
            CLimit = Limit;
        end

        for I = 1:size(overlapRegions, 1)

            if (overlapRegions(I, 2) - overlapRegions(I, 1) < CLimit * 60)
                badRegions(i) = I;
                i = i + 1;
            end

        end

        overlapRegions(badRegions, :) = [];

        %% dividing regions into segments of the length given by `Limit` in minutes
        I = 1;
        totalRegions = size(overlapRegions, 1);

        while (I <= totalRegions && Limit > 0)

            try
                numSegs = floor((overlapRegions(I, 2) - overlapRegions(I, 1)) / (Limit * 60));
            catch
            end

            if (numSegs < 2)
                overlapRegions(I, 2) = overlapRegions(I, 1) + Limit * 60;
            else
                tempA = overlapRegions(1:I - 1, :);
                tempB = overlapRegions(I + 1:end, :);
                range = overlapRegions(I, :);

                for B = 0:numSegs - 1
                    tempA(end + 1, 1) = range(1) + B * Limit * 60;
                    tempA(end, 2) = range(1) + (B + 1) * Limit * 60 - 1;
                end

                totalRegions = totalRegions + numSegs - 1;
                overlapRegions = [tempA; tempB];
            end

            I = I + 1;
        end

        %% optional: displaying output
        if (verbose)

            if (exist('h'))
                close(h);
            end

            h = figure;
            hold on;

            for I = 1:numSigs
                timespace = linspace(0, length(varargin{I * 2 - 1}) / varargin{I * 2}, length(varargin{I * 2 - 1}));
                plot(timespace, varargin{I * 2 - 1});
            end

            for I = 1:size(overlapRegions, 1)
                line([overlapRegions(I, 1) overlapRegions(I, 2)], [0 0], 'LineWidth', 5, 'Color', 'g');
            end

            disp("Clean data regions are highlighted in green");
            disp('press enter to continue...');
            pause;
        end

        if (review)
            disp('Output acceptable? passing data regions are highlighted in green');
            disp('Enter y for yes, the program will continue ');
            disp('else, to re-attempt enter `s` for strict filtering, `n` for normal filtering, and `l` for loose filtering.');
            disp('Enter e for to exit');
            str = input('>> ', 's');

            if (str == 'y')
                accepted = 1;

                if (exist('h'))
                    close(h);
                end

            else

                if (str == 's')
                    filterLevel = 2;
                end

                if (str == 'n')
                    filterLevel = 1;
                end

                if (str == 'l')
                    filterLevel = 0;
                end

            end

            if (str == 'e')

                if (exist('h'))
                    close(h);
                end

                return;
            end

        else

            if (exist('h'))
                close(h);
            end

        end

        %% converting time based markers back into indices for each set
        for I = 1:numSigs
            Regions{I}.intervals = round(overlapRegions * varargin{I * 2});
        end

    end

end

function [goodSpans] = IDgd(data, sampleRate, verbose, Limit, filterLevel, patchSignal)
    % author: Michael Hautman
    % assumes good data will be defined as both niether extremely regular or
    % irregular, finds intervals of data greater than 'Limit' minutes long that meet
    % this requirement

    %algorithm, takes a rolling standard deviation, and then takes a rolling
    %standard deviation of that output. then it takes the ratio of that ouput
    %to a large rolling deviation of the set. Data is filtered based on this ratio.
    %Filter Settings:
    % "strict" : ratio less than .2 and greater than .002
    % "normal" : ratio less than .3 and greater than .003
    % "loose"  : ratio less than .5 and greater than .003

    if (isstring(filterLevel))

        switch (filterLevel)
            case "loose"
                high = .5;
                low = .003;
            case "normal"
                high = .3;
                low = .003;
            case "strict"
                high = .5;
                low = .003;
            otherwise
                high = .3;
                low = .003;
        end

    else

        switch (filterLevel)
            case 0
                high = .5;
                low = .003;
            case 1
                high = .3;
                low = .003;
            case 2
                high = .2;
                low = .003;
            otherwise
                high = .3;
                low = .003;
        end

    end

    winSize = sampleRate * 4;
    ThresSize = sampleRate * 40;
    targetLength = Limit * 60 * sampleRate;

    %% detrend data
    [data, trend] = HautmanDetrend(data, sampleRate);

    %% filtering outliers to the variance
    threshold = sqrt(movvar(data, ThresSize));
    rollingSTD = sqrt(movvar(data, winSize));
    rollingSTD2 = sqrt(movvar(rollingSTD, winSize));

    regularityM = rollingSTD2 ./ threshold;

    goodIndices = regularityM > low & regularityM < high;

    if (sum(goodIndices) < length(data))
        trans = diff(goodIndices);
        starting = find(trans == 1);
        ending = find(trans == -1);

        if (starting(1) - ending(1) > 0)
            starting = [1 starting'];

            if (length(starting) > length(ending))
                ending = [ending' length(data)];
                ending = ending';
            end

            starting = starting';
        end

        goodSpans = [starting ending];

        %% Segments broken by small interferences are merged
        if (Limit > 0)
            I = 1;

            while (I < size(goodSpans, 1) && patchSignal)

                if (goodSpans(I + 1, 1) - goodSpans(I, 2) < winSize / 2)
                    tempstart = goodSpans(I, 1);
                    goodSpans(I, :) = [];
                    goodSpans(I, 1) = tempstart;
                end

                I = I + 1;
            end

            badSpans = find((goodSpans(:, 2) - goodSpans(:, 1)) < targetLength);
            goodSpans(badSpans, :) = [];
        end

    else
        goodSpans = [1 length(data)];
    end

    if (verbose)
        %testSpans = rollingSTD(goodSpans(:,1):goodSpans(:,2));
        plot(data);
        hold on;
        plot(rollingSTD);

        ave = mean(data);

        for I = 1:size(goodSpans, 1)
            line([goodSpans(I, 1) goodSpans(I, 2)], [0 0], 'LineWidth', 5, 'Color', 'r');
        end

    end

end

function [forward, backward] = fixRelate(forward, backward)
    % Fixes relationships between two sets that are supposed to be in an
    % alternating pattern. No assurances that the filtered data is actually
    % correct

    %remove duplicates
    forward = unique(forward);
    backward = unique(backward);

    % fix relationships between sets that should have a distinct 'set1(N) set2(N) set1(N+1) set2(N+1)' pattern
    badb = [];
    badf = [];

    cntf = 0;
    cntb = 0;
    i = 1;

    while (length(backward) >= (i + cntb) && length(forward) >= (i + cntf))
        % for not less than back o
        while (length(backward) >= (i + cntb) && (~(forward(cntf + i) < backward(cntb + i))))
            badb(end + 1) = cntb + i;
            cntb = cntb + 1;
        end

        while (length(forward) >= (i + cntf + 1) && (~(forward(cntf + i + 1) > backward(cntb + i))))
            badf(end + 1) = cntf + i + 1;
            cntf = cntf + 1;
        end

        i = i + 1;
    end

    backward(badb) = [];
    forward(badf) = [];
    sizelimit = min([length(forward) length(backward)]);
    forward = forward(1:sizelimit);
    backward = backward(1:sizelimit);
end
