function [aprDir aprDirIndex] = RollingDeriv(sig, ave)
    %Calculate averaged derivative
    sigdiff = diff(sig);
    aprDir = zeros(1, length(sig) - ave);

    for i = 1:length(sig) - ave
        aprDir(i) = sum(sigdiff(i:i + ave - 1)) / ave;
        aprDirIndex(i) = i + floor(ave / 2);
    end

end
