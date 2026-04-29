function [exclusions] = getExclusionIntervals(a, b, breaks)
% Given a collections of breaks and the overarching interval (a, b), return
% the intervals with in (a, b) that are included in breaks.

%since matlab reads arrays in column order, we need to convert the arrays
%to column order (take a transpose) and then reshape
exclusions = [a, reshape(breaks', 1, []), b]; % [a a1 a2 b1 b2 b]
exclusions = reshape(exclusions', 2, [])'; %[a a1; a2 b1; b2 b]
end

