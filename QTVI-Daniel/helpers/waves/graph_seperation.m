function graph_seperation(data, segment_idxs, precent)
    segment_idxs = [segment_idxs(1:end - 1)' segment_idxs(2:end)'];
    segment_lengths = segment_idxs(:, 2) - segment_idxs(:, 1) + 1;
    
    %% bad lengths
    med_length = median(segment_lengths);
    length_min = round(med_length - 2.5 * std(segment_lengths));
    length_max = round(med_length + 2.5 * std(segment_lengths));
    badlengths = (segment_lengths < length_min) | (segment_lengths > length_max);
    
    %% seperate segements from wave ignoring bad lengths
    good_segment_count = numel(badlengths(badlengths == 0));
    good_seg_lengths = segment_lengths(badlengths == 0);
    good_seg_ids = segment_idxs(badlengths == 0, :);
    good_segments = NaN(good_segment_count, max(good_seg_lengths) + 1);
    good_segment_peaks = NaN(good_segment_count, 2);
    
    for I = 1:size(good_segments, 1)
        good_segments(I, 1:good_seg_lengths(I)) = data(good_seg_ids(I, 1):good_seg_ids(I, 2));
        [good_segment_peaks(I, 1), good_segment_peaks(I, 2)] = max(good_segments(I, :));
    end

    alignedWaves = AlignWaves(good_segments, good_segment_peaks(:,2));
    
    diff_matrix = WaveDiff(alignedWaves);
   
    K = 20;
    a = 1;
    
    cnt = 10;
    error = nan(cnt,1);
    for groupNumber = 1:cnt
        clusteredLabels = gdlCluster(diff_matrix, groupNumber, K, a, true);
        
        tmp = nan(groupNumber,1);
        for x = 1:groupNumber
            tmp(groupNumber,1) = nansum(triu(WaveDiff(alignedWaves(clusteredLabels==x,:))),'All')^2;
        end
        error(groupNumber) = nansum(tmp);
    end
    
    
    idx = find(diff(error)>0, 1 );
    if idx >1
        idx = idx - 1;
    end
    idx
    clusteredLabels = gdlCluster(diff_matrix, idx, K, a, true);
    group_plot(alignedWaves, clusteredLabels)

    
    
%     %matvisual(diff_matrix*-1);
%     %matvisual(diff_matrix < .85);
%     [row,cols] = find((diff_matrix < precent) == 1);
%     idx = [cols row]; % order col row because we orig took difference column wize
%     
%     unconnected_nodes = setdiff(1:size(alignedWaves,1), 1:size(alignedWaves,1));
%     idx = [idx; [unconnected_nodes; unconnected_nodes]']; % add nodes that don't share an edge
%     G = graph(idx(:,1),idx(:,2));
%     p = plot(G);
%     bins = conncomp(G);
     %group_plot(alignedWaves, clusteredLabels)
    

end



function group_plot(data, bins)
    group_num = numel(unique(bins));
    c = linspecer(group_num);
    
    figure();
    hold on;

    for i = 1:group_num
        plot(data(bins==i,:)','Color',c(i,:));
    end

end