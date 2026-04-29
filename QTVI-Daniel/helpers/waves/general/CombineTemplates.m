function [bin_template_number, templates] = CombineTemplates(aligned_templates)     
    % kinda a dumb way to do this but for now it works...
    
    
    bin_template_number = zeros(size(aligned_templates,1),1);
    for z = 1:size(aligned_templates,1)
        
        if bin_template_number(z) == 0
            bin_template_number(z) = z;
            t = aligned_templates(z,~isnan(aligned_templates(z,:)));
            arr = bin_template_number==0;
            for i = 1:length(arr)
                if arr(i) == 1
                    t2 = aligned_templates(i,~isnan(aligned_templates(i,:)));
                    val = GaryDiff(t,t2);
                    
                    if min(val) >= 90
                        bin_template_number(i) = z; % Excellent
                    else
                        if (median(val(1:3)) >= 80 && val(1) >= 50)|| (min(val) >= 70)
                            bin_template_number(i) = z; % Acceptable
                        end
                    end
                end
            end
            
        end
        
    end
    


    
    
    
%     
%     
%     
%     
%     if min(sqibuf) >= 90
%         typeMnemonic(j) = 'E'; % Excellent
%     else
%         if (median(sqibuf(1:3)) >= 80 && sqibuf(1) >= 50 && c4(j) < 0.3 )|| (min(sqibuf) >= 70 && c4(j) < 0.3)
%             typeMnemonic(j) = 'A'; % Acceptable
%         else
%             typeMnemonic(j) = 'Q'; % Unacceptable
%         end
%     end
% 
%     
%     
%     
%     
%     
%     
%     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%     diff_matrix = WaveDiff(aligned_templates);
%     
%     
%     maxSimilaritys = zeros(size(aligned_templates,1)-1,1);
%     y = zeros(size(aligned_templates,1)-1,1);
%     for x = 1:size(aligned_templates,1)-1
%        [maxSimilaritys(x),y(x)] = min(diff_matrix(x,x+1:end)); % get most similar template 
%     end
%     
%     
%     u = unique(y);
%     bin_template_number = inf(size(aligned_templates,1),1);
%     bin_template_number(end) = size(aligned_templates,1); % final template will always be self similar
% 
%     for i = 1:length(u)
%         set = maxSimilaritys(y==u(i));
%         med = median(set);
%         s = abs(std(set));
%         lower = med -(s*2);
%         upper = med +(s*2);
%         outliers = (set < lower) | (set > upper);
%         
%         z = 1;
%         arr = y==u(i);
%         for q = 1:length(bin_template_number)-1
%             if arr(q) == 1
%                if outliers(z) == 1
%                    bin_template_number(q) = q;
%                else
%                    bin_template_number(q) = u(i);
%                end
%                z = z+1;
%             end
%         end
%     end
%     
%     
%     
%     
%     
%     
% %     
% %     s = std(maxSimilaritys);
% %     x = [maxSimilaritys zeros(length(maxSimilaritys),1)];
% %     bin_template_number = dbscan(x, s/3, 1);
% %     
%     group_num = numel(unique(bin_template_number));
%     c = linspecer(group_num);    
%     hold on;
%     u = unique(bin_template_number);
%     for i = 1:group_num
%         plot(aligned_templates(bin_template_number==u(i),:)','Color',c(i,:));
%         hold on;
%     end
    
end













function PlotWaveNum(aligned_templates)
    %plot(aligned_templates');
    hold on;
    str = sprintf('%i ',1:size(aligned_templates,1));
    strs = strsplit(str,' ');
    %legend(strs{1:size(aligned_templates,1)});
    good = ~isnan(aligned_templates(:,:));
    last_idx = arrayfun(@(x) find(good(x, :), 1, 'last'), 1:size(aligned_templates, 1));
    for i = 1:length(last_idx)
        text(last_idx(i), aligned_templates(i,last_idx(i)),[' ' strs(i)]);
    end
end


