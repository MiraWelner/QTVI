function plot_matlab_indices(mat_path, bins, output_dir, time_range)
% PLOT_RPEAKS  Plot ECG signal with R-peak markers and export R-peak indices to CSV.
%
%   plot_matlab_indices('path/to/wave_data.mat', [33 34], 'output_plots', [10 20])
%
%   mat_path   — path to the _wave_data.mat file
%   bins       — array of 0-based bin indices to plot
%   output_dir — folder where .png and .csv files will be saved (created if needed)
%   time_range — [start, end] in seconds (optional)

    if nargin < 3 || isempty(output_dir)
        output_dir = 'rpeaks_plots';
    end
    if nargin < 4
        time_range = [];
    end
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    data = load(mat_path);
    wave_data = data.wave_data;
    [~, file_id, ~] = fileparts(mat_path);
    file_id = strrep(file_id, '_wave_data', '');

    % CSV file for R-peak indices
    for k = 1:length(bins)

        bin_idx = bins(k);
        matlab_idx = bin_idx + 1;

        csv_path = fullfile(output_dir, sprintf('%s_bin%03d_rpeaks.csv', file_id, bin_idx));
        fid = fopen(csv_path, 'w');
        fprintf(fid, 'bin,r_peak_index\n');

        if matlab_idx > length(wave_data) || matlab_idx < 1
            fprintf('Bin %d out of range (max %d), skipping.\n', bin_idx, length(wave_data)-1);
            continue;
        end

        seg = wave_data{matlab_idx};
        ecg = seg.ecgSeg;
        fs  = seg.ecgSamplingRate;
        r   = seg.ecgRIndex;
        c_index_r = r-1;
        time = (0:length(ecg)-1) / fs;

        % Write R-peak indices to CSV
        for j = 1:length(c_index_r)
            fprintf(fid, '%d,%d\n', bin_idx, c_index_r(j));
        end

        % Filter by time range if specified
        if ~isempty(time_range)
            mask = (time >= time_range(1)) & (time <= time_range(2));
            if ~any(mask)
                fprintf('No data in range [%g, %g] for bin %d, skipping.\n', time_range(1), time_range(2), bin_idx);
                continue;
            end

            idx_start = find(mask, 1, 'first');
            idx_end = find(mask, 1, 'last');

            ecg = ecg(idx_start:idx_end);
            time = time(idx_start:idx_end);

            if ~isempty(r)
                r = r(r >= idx_start & r <= idx_end);
            end
        end

        fig = figure('Visible', 'off', 'Position', [70 70 1400 300]);

        plot(time, ecg, 'Color', [0.4 0.4 0.4], 'LineWidth', 0.5);
        hold on;

        if ~isempty(r)
            scatter((r-1)/fs, seg.ecgSeg(r), 50, 'r', 'v', 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end

        hold off;
        xlabel('Time (s)');
        ylabel('Amplitude');
        title(sprintf('%s — Bin %d | %d R-peaks | Fs = %g Hz', ...
            strrep(file_id,'_','\_'), bin_idx, length(r), fs));
        grid on;
        set(gca, 'GridAlpha', 0.15);

        if ~isempty(time)
            xlim([time(1) time(end)]);
        end

        legend({'ECG', 'R-peaks'}, 'Location', 'northeast');

        out_file = fullfile(output_dir, sprintf('%s_bin%03d_matlab.png', file_id, bin_idx));
        exportgraphics(fig, out_file, 'Resolution', 200);
        close(fig);
        fprintf('Saved: %s\n', out_file);
        fclose(fid);

    end

    fprintf('R-peak indices saved to: %s\n', csv_path);
end