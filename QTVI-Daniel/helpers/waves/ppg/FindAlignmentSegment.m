function [seg_movement,template] = FindAlignmentSegment(wave,minAmps,maxAmps)
    segment_idxs = [minAmps(1:end - 1) minAmps(2:end)];

    lengths = zeros(length(segment_idxs),1);
    for i = 1:length(segment_idxs)
        lengths(i) = length(wave(minAmps(i):maxAmps(i)));
    end
    
    segs = nan(length(segment_idxs),max(lengths));
    t20to80lengths = zeros(length(segment_idxs),1);
    for i = 1:length(segment_idxs)
        % time 20% and %80 min/max is alignment segment
        d = wave(maxAmps(i))-wave(minAmps(i));
        
        
        t20 = d * .2;
        t80 = d * .8;
        
        % closest to t20
        tmp = abs(wave(minAmps(i):maxAmps(i)) - (wave(minAmps(i))+t20));
        [~, beg] = min(tmp);
        beg = minAmps(i)+beg-1;
        % closest to t80
        tmp = abs(wave(minAmps(i):maxAmps(i)) - (wave(minAmps(i))+t80));
        [~, endidx] = min(tmp);
        endidx = minAmps(i)+endidx-1;
        
        % segment t20 through t80
        t20to80lengths(i) = length(beg:endidx);
        segs(i,1:t20to80lengths(i)) = wave(beg:endidx);
    end
    
    
    outliers =  isoutlier(t20to80lengths);
    length_min = floor(mean(t20to80lengths(outliers == 0)));
    
    segs_zeroed = segs - segs(:,1);

%     plot(segs_zeroed(outliers==0,:)','b');
%     hold on;
    %plot(segs_zeroed(outliers==1,:)','r');

    %% by diff max
    dif = diff(good_segments');
    dif = dif';
    good_seg_diff_peak = NaN(good_segment_count, 2);
    for I = 1:good_segment_count
        [good_seg_diff_peak(I, 1), good_seg_diff_peak(I, 2)] = max(dif(I, :));
    end

    alignedSegments = AlignWaves(good_segments, good_seg_diff_peak(:,2));

    
    
    %% Create Template
    template = nanmedian(segs_zeroed(outliers==0,1:length_min));
%     plot(template,'c');


    seg_movement = zeros(length(segs_zeroed),2); 
    for i = 1:length(seg_movement)
        seg = segs(i,~isnan(segs(i,:)));
        [seg_movement(i,1),seg_movement(i,2)] = alignWithTemplate(template,seg);
    end

    %% beats
    
%     beats = nan(length(segment_idxs),max(lengths));
%     for i = 1:length(minAmps)-1
%         beats(i,1:lengths(i)) = wave(minAmps(i):maxAmps(i));
%     end
%         
%     max_left = abs(min(seg_movement(seg_movement(:,1) < 0,1)));
%     dataspace_max = max(lengths + abs(seg_movement(:,1)));
%     aligned_segments = nan(length(segment_idxs), dataspace_max);
%     for i = 1:size(seg_movement)
%         beat = beats(i,~isnan(beats(i,:)));
%         beg = 1 + max_left + seg_movement(i,1);
%             
%         aligned_segments(i,beg:beg+length(beat)-1) = beat + seg_movement(i,2);
%        
%     end
%     
    
    max_left = abs(min(seg_movement(seg_movement(:,1) < 0,1)));
    dataspace_max = max(lengths + abs(seg_movement(:,1)));
    aligned_segments = nan(length(segment_idxs), dataspace_max);
    for i = 1:size(seg_movement)
        seg = segs(i,~isnan(segs(i,:)));
        beg = 1 + max_left + seg_movement(i,1);
            
        aligned_segments(i,beg:beg+length(seg)-1) = seg + seg_movement(i,2);
       
    end
    
    
end


function [movex,movey] = alignWithTemplate(template,seg)
    tempMin = min(template);
    xtemplate = template -tempMin + 1;
    segMin = min(seg);
    xseg = seg -segMin + 1;

    [cor_seq,lags] = xcorr(xtemplate,xseg);
    [~,d]=max(cor_seq);

    movex = lags(d);  
    
    if movex < 0
        seg = seg(1+abs(movex):length(template)+abs(movex));
    else
        template = template(movex+1:movex+length(seg));
    end
    
    %% first pass
    % get estimated lowest error of segment by setting each point in
    % segment to height of each point in template and calculating error
    yerror_corse = zeros(length(template),length(seg));

    for t = 1:length(template)
       for x = 1:length(seg)
           difference = template(t) - seg(x);
           tempseg = seg + difference;
           yerror_corse(t,x) = sum((tempseg - template).^2);
       end
    end
    
    minMatrix = min(yerror_corse(:));
    [row,col] = find(yerror_corse==minMatrix);
    firstPass_diff = template(row) - seg(col);
    if length(firstPass_diff) > 1
       firstPass_diff = firstPass_diff(1); 
       row = row(1);
       col = col(1);
    end
    
    seg_firstPass = seg + firstPass_diff;

    
    if col == 1 
        upperrange = abs(seg_firstPass(1) - seg_firstPass(2));
        lowerange = 0;
    elseif col == length(seg_firstPass)
        upperrange = 0;
        lowerange = abs(seg_firstPass(end) - seg_firstPass(end-1));
    else
        upperrange = abs(seg_firstPass(col) - seg_firstPass(col+1));
        lowerange = abs(seg_firstPass(col) - seg_firstPass(col-1));
    end
    
    
    %% second pass
    % finer granularity
    steps_low = -linspace(0, lowerange, 5001);
    steps_high = linspace(0, upperrange, 5000);
    steps = [fliplr(steps_low(2:end)) steps_high];
    
    yerror_fine = zeros(length(steps),1);

    for x = 1:length(steps)
%         p = plot(seg_firstPass + steps(x),'.-m');
%         pause(0.00001);
%         delete(p)
        yerror_fine(x) = sum(((seg_firstPass + steps(x)) - template).^2);
    end
    
    [~,min2] = min(yerror_fine);
    movey = firstPass_diff+steps(min2);
end


