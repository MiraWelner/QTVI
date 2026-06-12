function [template_matrix] = CreatePPGTemplates(bins, std_multiplier, dbg_plot)
    % calculate templates
    templates = cell(length(bins), 1);

    if dbg_plot == 1
        figs = cell(length(bins), 1);
    end

    for i = 1:length(bins)
        %[alignSegs] = FindAlignmentSegment(bins{i}.ppgSeg, bins{i}.ppgMinAmps,bins{i}.ppgMaxAmps);
        if bins{i}.bad_segment == 1
            templates{i} = [];
            continue
        end
        try
            if dbg_plot == 1
                [templates{i}, figs{i}] = EnsembleTemplate(bins{i}.ppgSeg, bins{i}.ppgMinAmps, std_multiplier, 'ppg', dbg_plot, ['Template #' num2str(i) ' of ' num2str(length(bins))]);
            else
                templates{i} = EnsembleTemplate(bins{i}.ppgSeg, bins{i}.ppgMinAmps, std_multiplier, 'ppg', dbg_plot, ['Template #' num2str(i) ' of ' num2str(length(bins))]);
            end
        catch
            continue
        end

    end

    % make templates all same size
    m = max(cellfun(@length, templates));
    template_matrix = nan(length(bins), m);

    for i = 1:length(bins)
        template_matrix(i, 1:length(templates{i})) = templates{i};
    end

    if dbg_plot == 1
        ShowDbgPlots(figs);
    end

end
