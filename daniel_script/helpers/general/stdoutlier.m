function [outliers] = stdoutlier(data, multiplier, mean_window, direction, debug_plot)
    %% data == nx1 or 1xn matrix, Multiplier == std multipler (normally ~2.5 or 3)
    %% mean_window == length(data) * .02 normally good value
    x = movmean(diff(data), mean_window);
    s = std(diff(data));
    d = diff(data);
    upper_bound = x + s * multiplier;
    lower_bound = x - s * multiplier;
    if strcmp(direction,'lower')
        weridOnes = (d < lower_bound);
    elseif strcmp(direction,'upper')
        weridOnes = (d > upper_bound);
    else
        weridOnes = (d > upper_bound | d < lower_bound);
    end
    
    if debug_plot
        figure;
        plot(d, '.');
        hold on;
        plot(x);
        plot(upper_bound);
        plot(lower_bound);
        t = (1:length(d));
        plot(t(weridOnes == 1), d(weridOnes == 1), 'or');
        hold off;
    end

    outliers = zeros(length(data), 1);

    for i = 1:length(data) - 1

        if weridOnes(i) == 1
            outliers(i:i + 1) = 1;
        end

    end

    outliers = logical(outliers);
end