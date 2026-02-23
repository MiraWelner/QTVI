function [diff_matrix] = WaveDiff(waves)
    diff_matrix = zeros(size(waves,1))*1000;
    for i = 1:size(waves,1)
        for j = i+1:size(waves,1)
            diff_matrix(i,j) = nansum(waves(i,:)-waves(j,:))^2;
            %diff_matrix(j,i) = diff_matrix(i,j);
        end
    end
end