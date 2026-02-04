
function [out] = bestFitOfBeatTrain(sectionBefore, sectionAfter, dataToMatch, beats)
    % for each point in data to match check where the best matching r to
    % the 'beat train' is and output that
    try
        sectionBefore = sectionBefore(end-beats+1:end);
    catch
    end
    
    try
        sectionAfter = sectionAfter(1:beats);
    catch
    end
    
    
    
    
    out = [];
    
    if length(dataToMatch) == 2
        medBefore = round(mean(diff(sectionBefore)));
        medAfter = round(mean(diff(sectionAfter)));

        endBefore = sectionBefore(end);
        beginAfter = sectionAfter(1);

        beforeRate = endBefore:medBefore:beginAfter;
        afterRate = endBefore:medAfter:beginAfter;

        if length(beforeRate) > length(afterRate)
            mid = round((beforeRate(1:numel(afterRate))+afterRate)/2);
        else
            mid = round((afterRate(1:numel(beforeRate))+beforeRate)/2);
        end

        if length(mid) < 2
           mid = endBefore + medBefore;
        else
           mid = mid(2);
        end

        errors = zeros(length(dataToMatch),1);
        for z = 1:length(dataToMatch)
            errors(z) = abs(mid - dataToMatch(z));
        end
        [~,minIdx] = min(errors);
        out(end+1) =dataToMatch(minIdx);
        return 
    else
        for i = 1:length(dataToMatch)
            medBefore = round(mean(diff(sectionBefore)));
            medAfter = round(mean(diff(sectionAfter)));

            endBefore = sectionBefore(end);
            beginAfter = sectionAfter(1);

            beforeRate = endBefore:medBefore:beginAfter;
            afterRate = endBefore:medAfter:beginAfter;

            if length(beforeRate) > length(afterRate)
                mid = round((beforeRate(1:numel(afterRate))+afterRate)/2);
            else
                mid = round((afterRate(1:numel(beforeRate))+beforeRate)/2);
            end

            if length(mid) < 2
               mid = endBefore + medBefore;
            else
               mid = mid(2);
            end


            errors = zeros(length(dataToMatch),1);
            for z = 1:length(dataToMatch)
                errors(z) = abs(mid - dataToMatch(z));
            end

            [~,minIdx] = min(errors);
            out(end+1) = dataToMatch(minIdx);

            sectionBefore(end+1) = dataToMatch(minIdx);
            dataToMatch(dataToMatch<=minIdx) = [];

        end
    end
    out = unique(out);
end