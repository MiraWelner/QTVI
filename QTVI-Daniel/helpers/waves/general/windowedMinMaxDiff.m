function [amp] = windowedMinMaxDiff(data, sampling_rate, window_size_mins)
    % WindowedMinMaxDiff - take a windowed difference between the max and
    % min value of a windowed period defined in mins. Called an 'Ampogram'
    %
    % Syntax:  [amp] = WindowedMinMaxDiff(data, sampling_rate, window_size_mins)
    %
    % Inputs:
    %    data - data in column format (1xN)
    %    sampling_rate - sampling rate of data in hz
    %    window_size_mins - window size in minutes 
    %
    % Outputs:
    %    amp - matrix with the difference of the min and max values 
    %    contained within each window.
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    % Author: Daniel Wendelken
    % email address: wendeldr@ucmail.uc.edu
    % Spring 2019; Last revision: 1-23-2019
    desired_window_size_in_samples = (sampling_rate * 60) * window_size_mins;
    data_len_mins = length(data)/(desired_window_size_in_samples);
    windows = ceil(data_len_mins);
    
    window_breaks = (1:windows) * desired_window_size_in_samples;
    window_breaks(end) = length(data);
    window_breaks = [1 window_breaks];
    
    amp = zeros(windows,1);
    for i = 2:length(window_breaks)
        tmp = data(window_breaks(i-1):window_breaks(i));        %for full min
        %tmp = data(window_breaks(i-1):window_breaks(i-1)+2* sampling_rate);        %for first 2 sec
        amp(i-1) = max(tmp)-min(tmp);
    end
end
