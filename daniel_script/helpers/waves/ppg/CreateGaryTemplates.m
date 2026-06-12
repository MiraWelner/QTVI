
function [template_matrix,type] = CreateGaryTemplates(bins, windowlen, Fs)
    % calculate templates
    templates = cell(length(bins), 1);
    type = zeros(length(bins), 1);
    parfor i = 1:length(bins)
        anntime = bins{i}.ppgMinAmps;
        wave = PPGmedianfilter(bins{i}.ppgSeg, Fs, Fs);
        [t t2 v] = template_pleth(wave(1:min(windowlen, length(wave))), anntime(anntime < min(windowlen, length(wave))), 0, Fs);
        if v < 1  % Current template invalid 
            template = [];
            type(i) = v;
        else
            % Using t2 if available
            if v > 1
                type(i) = v;
                template = t2;
            else
                type(i) = v;
                template = t;
            end
        end
        templates{i} = template;
    end
    
    % make templates all same size
    m = max(cellfun(@length, templates));
    template_matrix = nan(length(bins), m);

    for i = 1:length(bins)
        template_matrix(i, 1:length(templates{i})) = templates{i};
    end
end