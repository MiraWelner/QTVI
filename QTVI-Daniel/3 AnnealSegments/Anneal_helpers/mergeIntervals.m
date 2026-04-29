function [good_sections] = mergeIntervals(good_sections, sampleRate, minimumBinSize)
%mergeIntervals merge various good sections, based on sample rate, and the
%minimum valid intervals. 
%merge together index's which share a border and are both moving
i = 1;
    while (i < size(good_sections, 1))
        if good_sections(i, 2) == good_sections(i + 1, 1) && good_sections(i, 4) ~= 0 && good_sections(i + 1, 4) ~= 0
            seg1_size_min = (good_sections(i, 2) - good_sections(i, 1)) / sampleRate / 60;
            seg2_size_min = (good_sections(i + 1, 2) - good_sections(i + 1, 1)) / sampleRate / 60;

            if seg1_size_min + seg2_size_min >= minimumBinSize
                good_sections(i, 3) = 0;
                good_sections(i, 4) = 0;
            else
                [~, idx] = max([seg1_size_min seg2_size_min]);
                idx = idx - 1;
                good_sections(i, 3) = good_sections(i + idx, 3);
                good_sections(i, 4) = 1;
            end
            good_sections(i, 2) = good_sections(i + 1, 2);
            good_sections(i + 1, :) = [];
        end
        i = i + 1;
    end

end

