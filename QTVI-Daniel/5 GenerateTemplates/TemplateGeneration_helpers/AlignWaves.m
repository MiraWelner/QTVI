function [alignedWaves, move_dist] = AlignWaves(waves, wave_alignment_points)
    % Author: Daniel Wendelken
    % AlignWaves - aligns matrix of waves to the specified point based on the
    % desired point in each wave.
    % Syntax:  [alignedWaves] = AlignWaves(waves, wave_alignment_points)
    %
    % Inputs:
    %    waves - NxN matrix where each row represents a different wave
    %    wave_alignment_points - 1xN matrix representing the point on each
    %    wave that should be aligned

    %
    % Outputs:
    %    alignedWaves - matrix with the waves aligned. matrix may be
    %    different dimensions in columns depending on if alignment requires
    %    more or less columns.

    wave_count = size(waves, 1);

    max_left = max(wave_alignment_points);

    lengths = zeros(wave_count, 1);

    for i = 1:wave_count
        lengths(i) = length(waves(i, ~isnan(waves(i, :))));
    end

    % determine greatest distance to right of peak
    max_right = max(lengths - wave_alignment_points);

    % determine distance and direction of movement (neg == left,pos==right)
    % for each wave
    move_dist = max_left - wave_alignment_points;

    alignedWaves = nan(wave_count, max_left + max_right);

    for I = 1:wave_count
        alignedWaves(I, 1 + move_dist(I):lengths(I) + move_dist(I)) = waves(I, ~isnan(waves(I, :)));
    end

    %% align peaks vert
    avg_peak = nanmean(alignedWaves(:, max_left));

    for I = 1:wave_count
        alignedWaves(I, :) = alignedWaves(I, :) - (avg_peak + alignedWaves(I, max_left));
    end

end
