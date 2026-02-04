function [bin_numbers, bin_templates, foot_locations] = CombineTemplatesGraph(aligned_templates, diff_matrix, threshold_percent, dbg_plot)
    % right now this does nothing. Just outputs a template for bin number I
    % didn't have the time to revivisit "combining" the templates

    [row,cols] = find((diff_matrix <= threshold_percent & diff_matrix ~= 0) == 1);
    
    weights = diff_matrix((diff_matrix <= threshold_percent & diff_matrix ~= 0) == 1);
    idx = [cols row weights]; % order col row because we orig took difference columnwize
    unconnected_nodes = setdiff(1:size(aligned_templates,1), unique([unique(cols); unique(row)]));
    idx = [idx; [unconnected_nodes; unconnected_nodes; Inf(1, length(unconnected_nodes))]']; % add nodes that don't share an edge
    %G = digraph(idx(:,1),idx(:,2),idx(:,3));
    G = graph(idx(:,1),idx(:,2));
    %bin_numbers = conncomp(G);
    bin_numbers = 1:size(aligned_templates,1);

    bin_templates = cell(numel(unique(bin_numbers)),1);
    foot_locations = ones(numel(unique(bin_numbers)),1);
    bin_template_nans = nan(numel(unique(bin_numbers)), size(aligned_templates,2));
    for i = 1:numel(unique(bin_numbers))
        template = nanmean(aligned_templates(bin_numbers==i,:), 1);
        bin_template_nans(i,1:length(template)) = template;
        template = template(~isnan(template));
        [~,foot_locations(i)] = find_foot_pulseox(template, 0);
        bin_templates{i} = template;
    end
    
    % debugging
    if dbg_plot == 1
        fig = figure('Name', 'Template Combinations','visible','off');
        matvisual(diff_matrix,subplot(3,2,1));
        title('Diff Matrix');
        
        d = (diff_matrix <= threshold_percent & diff_matrix > 0);
        matvisual(d, subplot(3,2,2));
        title(['Diff where value <= ' num2str(threshold_percent)]);
        
        subplot(3,2,3);
%       plot(G,'LineWidth', G.Edges.Weight/4, 'ArrowSize', 5);
        plot(G);
        title('Diff Graph');
        
        subplot(3,2,4);
        plot(aligned_templates');
        legend('NumColumns',2);
        %PlotWaveNum(subplot(3,2,4), aligned_templates);
        title(['Original ' num2str(size(aligned_templates,1)) ' Templates']);
       
        GroupPlot(subplot(3,2,[5 6]), aligned_templates, bin_numbers);
        hold on;
        p=plot(subplot(3,2,[5 6]), bin_template_nans','c', 'LineWidth', 1);
        %PlotWaveNum(subplot(3,2,[5 6]), aligned_templates);
        title(['Final ' num2str(numel(unique(bin_numbers)))  ' bin groupings']);
        legend(p,'Templates');
        ShowDbgPlots({fig});

    end
end

function PlotWaveNum(ax, aligned_templates)
    hold(ax, 'on');
    str = sprintf('%i ',1:size(aligned_templates,1));
    strs = strsplit(str,' ');
    %legend(strs{1:size(aligned_templates,1)});
    good = ~isnan(aligned_templates(:,:));
    last_idx = arrayfun(@(x) find(good(x, :), 1, 'last'), 1:size(aligned_templates, 1));
    for i = 1:length(last_idx)
        text(last_idx(i), aligned_templates(i,last_idx(i)),[' ' strs(i)]);
    end
    hold(ax, 'off');
end

function GroupPlot(ax, data, bins)
    group_num = numel(unique(bins));
    c = linspecer(group_num);    
    hold(ax, 'on');
    for i = 1:group_num
        plot(data(bins==i,:)','Color',c(i,:));
    end
    hold(ax, 'off');
end
