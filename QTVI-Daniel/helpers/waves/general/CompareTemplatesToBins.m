function [errors] = CompareTemplatesToBins(bins, bin_template_numbers, templates, template_foot_locations, mean_thresh, median_thresh, std_thresh, standard_error_thresh, dbg_plot)
    errors = inf(length(bins),1);
    if dbg_plot == 1
        figs = cell(length(bins),1);
    end
    for i = 1:length(bins)
        template_num = bin_template_numbers(i);
        foot_location = template_foot_locations(template_num);
        template = templates{template_num};
        wave = bins{i}.ppgSeg;
        wave_indexes = bins{i}.widened_RRindex;
        
        template_length = length(template);
        if template_length == 0
            return
        end

        wave_indexes_length = length(wave_indexes);

        template_array = nan(length(wave),1);
        for q = 1:wave_indexes_length
            val = wave(wave_indexes(q));
            difference = foot_location - val;
            tmp = template - difference;

            if q == 1
               left = min(foot_location, wave_indexes(q)); % -1 because matlab starts at 1...
            else
               left = min(foot_location, wave_indexes(q)-wave_indexes(q-1));
            end

            if q == wave_indexes_length
                right = min(template_length-foot_location, length(wave));
            else
                right = min(template_length-foot_location, wave_indexes(q+1)-wave_indexes(q));
            end

            beginpos = wave_indexes(q) - left + 1;
            endpos = wave_indexes(q) + right;

            if endpos - beginpos +1 < template_length || beginpos >= numel(wave) || endpos >= numel(wave)
               continue
            end

            template_array(beginpos:endpos) = tmp';
        end

        errors(i) = nansum((template_array - wave).^2);
        fig = figure('Name', ['Bin #' num2str(i) ' vs. Template #' num2str(template_num)],'visible','off');
        plot(wave');
        hold on;
        plot(template_array');
        title(['Bin #' num2str(i) ' vs. Template #' num2str(template_num)]);
        figs{i} = fig;
    end

    if dbg_plot == 1
        ShowDbgPlots(figs);
    end
end