function [return_bins] = CombineBins(bins, bin_data)
    [~, n] = RunLength(bins);
    return_bins = cell(numel(n), 1);

    bin_idx = 1;
    i_diff = 0;

    for i = 1:length(bins)
        tmp = [];
        template_num = bins(i + i_diff);
        number_in_bin = n(bin_idx);

        if number_in_bin == 1
            tmp = bin_data{bin_idx};
            tmp.template_num = template_num;
            return_bins{bin_idx} = tmp;
        else
            z = i + i_diff;

            for q = z:z + number_in_bin - 1

                if isempty(return_bins{bin_idx})
                    tmp = bin_data{q};
                    tmp.template_num = template_num;
                    return_bins{bin_idx} = tmp;
                else
                    tmp = return_bins{bin_idx};
                    tmp.ppgSeg = [tmp.ppgSeg; bin_data{q}.ppgSeg];
                    idx = bin_data{q}.ppgMinAmps + tmp.ppgMinAmps(end) + 1;
                    tmp.ppgMinAmps = [tmp.ppgMinAmps; idx];

                    tmp.ppgSeg = [tmp.ecgSeg; bin_data{q}.ecgSeg];
                    idx = bin_data{q}.ecgRIndex + tmp.ecgRIndex(end) + 1;
                    tmp.ecgRIndex = [tmp.ecgRIndex; idx];

                    return_bins{bin_idx} = tmp;
                end

            end

            i_diff = i_diff + number_in_bin - 1;
        end

        bin_idx = bin_idx + 1;

        if i + i_diff >= length(bins)
            break
        end

    end

end
