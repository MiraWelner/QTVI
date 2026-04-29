function [bool] = intervalIsOutlier(interval,location,data, multiplier, mean_window)


    x = movmean(diff(data), mean_window);
    d = diff(data);
    s = std(d);
    upper_bound = x + s * multiplier;
    lower_bound = x - s * multiplier;
    
%     figure;
%     plot(d, '.');
%     hold on;
%     plot(x);
%     plot(upper_bound);
%     plot(lower_bound);
%     t = (1:length(d));
%     plot(t(location), interval, 'or');
%     hold off;
    
    if interval >= lower_bound(location) && interval <= upper_bound(location)
        bool = 0;
    else
        bool = 1;
    end
end