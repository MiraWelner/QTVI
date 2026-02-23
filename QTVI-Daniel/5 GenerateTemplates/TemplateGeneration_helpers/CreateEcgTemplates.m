function [ecgTemplates, ppg_alignment_point, avg_r_expand] = CreateEcgTemplates(bins, std_multiplier, dbg_plot)
    ecgTemplates = cell(length(bins), 1);
    ppg_alignment_point = zeros(length(bins), 1);
    avg_r_expand = zeros(length(bins), 1);
    if dbg_plot == 1
        figs = cell(length(bins), 1);
    end

    for i = 1:length(bins)
        reduced_size = round(length(bins{i}.ecgSeg) / 10);

        ecg = bins{i}.ecgSeg(1:reduced_size);
        r = bins{i}.ecgRIndex(bins{i}.ecgRIndex < reduced_size - 1);
        
        lens = floor(diff(r)/5);
        avg_r_expand(i) = round(median(lens));
        try
            if dbg_plot == 1
                [ecgTemplates{i}, figs{i}] = EnsembleTemplate(ecg, r, std_multiplier, 'ecg', dbg_plot, ['Template #' num2str(i) ' of ' num2str(length(bins))], lens);
            else
                ecgTemplates{i} = EnsembleTemplate(ecg, r, std_multiplier, 'ecg', dbg_plot, ['Template #' num2str(i) ' of ' num2str(length(bins))], lens);
            end


            ppg_alignment_point(i) = round(median(bins{i}.pairs(:,1) - bins{i}.pairs(:,2)));
        catch
            continue
        end
    end

    if dbg_plot == 1
        ShowDbgPlots(figs);
    end

end
