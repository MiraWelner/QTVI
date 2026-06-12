% Rolling window functions
function [rwin] = rollingWindow(vector, winsize)
    vector = vector(:);
    vecsize = length(vector);
    rwin = NaN(winsize, vecsize);
    %     rwin = mat2cell(rwin, size(rwin,1), ones(size(rwin,2), 1));
    for i = 1:winsize
        tmp = vector(1:end - i + 1);
        rwin(i, i:end) = tmp;
    end

end
