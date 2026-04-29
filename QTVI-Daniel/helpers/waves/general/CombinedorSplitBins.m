function [errors, templates] = CombinedorSplitBins(bins, bin_templates, foot_locations, data, ecgSamplingRate)
    errors = inf(length(bins), 1);

    for i = 1:length(bins)
        errors(i) = CompareTemplateToBin(bins{i}.ppgSeg, bin_templates(bins{i}.template_num, :), bins{i}.ppgMinAmps, ...
            foot_locations(bins{i}.template_num), 4, 5, 6, 7);
    end

    for i = 1:length(bins)
        widened_rr_index = ExpandIndexs(bins{i}.ecgRIndex, ecgSamplingRate / 4, ecgSamplingRate / 4, numel(bins{i}.ecgSeg) - 1);
        bins{i}.ecg_template = EnsembleTemplate(bins{i}.ecgSeg, widened_rr_index, 2.5, 'ecg', 0);
    end

    templates = cell(length(bins), 1);

    for i = 1:length(bins)
        tmp = bin_templates(i, :);
        templates{i, 1} = tmp(~isnan(tmp));
        tmp = bins{i}.ecg_template;
        templates{i, 2} = tmp(~isnan(tmp));
    end

end
