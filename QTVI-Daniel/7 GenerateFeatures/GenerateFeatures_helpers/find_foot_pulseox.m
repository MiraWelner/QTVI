function [val,idx] = find_foot_pulseox(data,dbg_plot)
    if isempty(data) || length(data) < 4
        idx = 1;
        val = 0;
        return
    end
    
    
    data = data - max(data,[],2);
    
    [~,m] = max(data,[],2);   
    a = cell(length(m),1);
    b = cell(length(m),1);
    for i = 1:length(m)
        if numel(data(~isnan(data(i,:)))) < 2
            a{i}=1;
            b{i}=1;
        else
        	[a{i},b{i}] = findpeaks(diff(data(i,:),1,2));
        end
    end
    
    points_of_intrest = cell(length(m),1);
    for i = 1:length(points_of_intrest)
        tmp = zeros(length(a{i}(b{i} <= m(i))),2);
        if isempty(tmp)
            tmp = zeros(1,2);
            tmp(:,1) = m(i);
            tmp(:,2) = 0;
            points_of_intrest{i} = tmp;
            continue
        end
        tmp(:,1) = b{i}(b{i} <= m(i));
        tmp(:,2) = a{i}(b{i} <= m(i));
        points_of_intrest{i} = tmp;
    end 
    
    sorted = cellfun(@(x) sortrows(x,'descend'),points_of_intrest,'UniformOutput', false);
    
    diffpeaks = zeros(length(sorted),2);
    overshoot = 0;
    for i = 1:length(sorted)
        wave = sorted{i};
        for x = 1:size(wave,1)
            curx = wave(x,1);
            cury = wave(x,2);
            if cury > diffpeaks(:,1)
                diffpeaks(:,1) = cury;
                diffpeaks(:,2) = curx;
                overshoot = 0;
            else
                overshoot = overshoot + 1;
            end
            if overshoot > 2
                break
            end
        end
    end
    
    for i = 1:size(diffpeaks,1)
        if diffpeaks(i,2) == 0
            [diffpeaks(i,1), diffpeaks(i,2)] = max(diff(data,1,2),[],2);
        end
    end

    beginPoints = [ones(length(data(:,1)),1) data(:,1)]; % start at 1 since matlab starts at 1...
    
    endPoints = zeros(size(diffpeaks,1),2);
    for i = 1:size(diffpeaks,1)
        if diffpeaks(i,2) < 1
            endPoints(i,1) = 1;
            endPoints(i,2) = data(i,1);
        else
            endPoints(i,1) = diffpeaks(i,2);
            endPoints(i,2) = data(i,diffpeaks(i,2));
        end
    end

    p1_prime = [1 0] - beginPoints; % start at 1 since matlab starts at 1...
    p2_prime = endPoints + p1_prime;

    moved = data + p1_prime(:,2);
    
    idx = zeros(size(moved,1),1);    
    val = zeros(size(moved,1),1);

    rotated = cell(size(moved,1),1);
%     theta = zeros(size(moved,1),1);
%     for i = 1:size(moved,1)
%         theta(i) = atand(p2_prime(i,1)/p2_prime(i,2));
%         R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
%         tmp =  R * [1:diffpeaks(i,2); moved(i,1:diffpeaks(i,2))];
%         rotated{i} = tmp(1,:);
%         [val(i),idx(i)] = max(rotated{i});
%     end
    
    theta = zeros(size(moved,1),1);
    for i = 1:size(moved,1)
        theta(i) = atand(p2_prime(i,1)/p2_prime(i,2));
        if theta(i) < 0
            theta(i) = theta(i) *-1;
        end
        R = [cosd(theta(i)) -sind(theta(i)); sind(theta(i)) cosd(theta(i))];
        tmp =  R * [1:diffpeaks(i,2); moved(i,1:diffpeaks(i,2))];
        rotated{i} = tmp(1,:);
        [val(i),idx(i)] = max(rotated{i});
    end

    if dbg_plot == 1
        figs = cell(size(diffpeaks,1),1);
        for i = 1:size(diffpeaks,1)
            fig = figure('Name', ['Find foot | #' num2str(i) ' of ' num2str(size(diffpeaks,1))],'NumberTitle', 'off', 'visible', 'off');
            plot(data(i,:),'b');
            hold on;
            plot(endPoints(i,1), endPoints(i,2),'.r');
            plot(data(i,1:endPoints(i,1))','Color',rgb('Cyan'));
            plot([1 endPoints(i,1)], [data(i,1) endPoints(i,2)], ':','Color', rgb('darkslategray'));
            plot(idx(i),data(i,idx(i)),'og');
            plot(rotated{i},'--','Color',rgb('dimgray'));
            vline(idx(i),'r',num2str(theta(i)));
            hline(0,'k');
            figs{i} = fig;
        end
        ShowDbgPlots(figs);
    end
    
end