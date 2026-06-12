function [dicrotic_notch_idx] = dumbDicrotic(beat,sp_ratio)
    if isrow(beat)
        beat = beat';
    end
    try
        beat = nanfastsmooth(beat,15);
    catch
       dicrotic_notch_idx = nan; 
       return
    end
    
%     try
        [~,mx] = max(beat);
        x1 = mx;
%         [~,mn] = min(beat(mx:end));
%         x2 = x1+mn-1;
        x2 = length(beat);
        
        if x1 == x2 | x2-x1 < 2
            dicrotic_notch_idx = nan; 
            return
        end
        
        y1 = beat(x1);
        y2 = beat(x2);

        pmax = x1;
        pmin = x2;
        

        if exist('sp_ratio','var') == 1           
            endlen = length(beat)-pmax;
            tmp = round(endlen*sp_ratio);
            notch_estimate = pmax + tmp;
            p_min_dpdt_region = pmax + round((notch_estimate-pmax)*(2/3));
            if p_min_dpdt_region < pmax+2
                p_min_dpdt_region = pmax + round((pmin-pmax)/3);
            end
            init_EP = notch_estimate + round((length(beat)-notch_estimate)*(1/2));
            if init_EP <= p_min_dpdt_region
                init_EP = pmax + round(((pmin - pmax) * 3/4));
            end
        else
            p_min_dpdt_region = pmax + round((pmin-pmax)/3);
        	init_EP = pmax + round(((pmin - pmax) * 3/4));
        end 

        slice = beat(pmax:p_min_dpdt_region);
        [~,p_min_dpdt] = min(diff(slice));

        p_min_dpdt = p_min_dpdt + pmax;

        p_half = beat(pmax) - (beat(pmax) - beat(p_min_dpdt))/2;
        MINS = abs(slice-p_half);
        potential_sp = zeros(length(MINS),2);
        potential_sp(:,1) = (0:length(MINS)-1) + pmax;
        potential_sp(:,2) = MINS;
        [~,idx] = sort(potential_sp(:,2));
        potential_sp = potential_sp(idx,:);
    %     potential_sp(potential_sp == 0,:) = [];


        x = 1;
        found = 0;
        while x < length(potential_sp) && found == 0
            SP = potential_sp(x,1);

            transform = shear_transform(SP:init_EP,beat(SP:init_EP));
            if sum(transform < beat(SP:init_EP))/length(beat(SP:init_EP)) < .5
               found = 1 ;
            end
            x = x+1;
        end

        EP = init_EP;
        while EP > SP+1
            shearpressure = beat(SP:EP);
            time = SP:EP;

            c = polyfit([time(end),time(1)], [shearpressure(1), shearpressure(end)], 1);
            m = c(2)/(c(1));
            b = 0;
            shearline = zeros(length(shearpressure),1);
            for x=1:length(shearpressure)
               shearline(x)  = m * x + b;
            end

            norm_line = (shearline-min(shearline))./(max(shearline)-min(shearline)) ;
            norm_press = (shearpressure-min(shearpressure))./(max(shearpressure)-min(shearpressure)) ;
            norm_time = (time-min(time))./(max(time)-min(time));

            if orthoginal_dist_thresh(norm_time, norm_line, norm_press, 0.3)
                EP = EP-1;
            else
                break;
            end
        end

        transform = shear_transform(SP:EP,beat(SP:EP));
        [~, min_shear] = min(transform);
        min_shear = min_shear+SP-1;
        [~,start_diastolic_relax] = max(beat(min_shear:pmin));
        start_diastolic_relax = start_diastolic_relax + min_shear-1;
        [~,dicrotic_notch_idx] = min(beat(min_shear:start_diastolic_relax));
        dicrotic_notch_idx = dicrotic_notch_idx + min_shear-1;

%             figure();
%             plot(beat);
%             hold on;
%             hline(p_half);
%             plot(pmax,beat(pmax),'rd');
%             plot(pmin,beat(pmin),'cd');
%             plot(p_min_dpdt, beat(p_min_dpdt),'rd');
%             plot(SP, beat(SP),'yo');
%             plot(EP, beat(EP),'mo');
%             plot(SP:EP,beat(SP:EP),'r')
%             vline(min_shear)
%             vline(start_diastolic_relax)
%             
%             plot(dicrotic_notch_idx,beat(dicrotic_notch_idx),'g*')
%             
%             legend('beat','pmax','pmin','pmin_(dpdt)','SP');
    %     catch
    %         
    %         [~,mx] = max(beat);
    %         x1 = mx;
    %         [~,mn] = min(beat(mx:end));
    %         x2 = x1+mn-1;
    % 
    % 
    %         y1 = beat(x1);
    %         y2 = beat(x2);
    %         slope = (y2 - y1) / (x2 - x1);
    % 
    %         sub = zeros(1,length(beat(x1:x2)));
    % 
    %         for q = 1:length(beat(x1:x2))
    %             sub(q) = slope * q + y1;
    %         end
    % 
    %     %     figure;subplot(2,1,1);plot(beat);hold on;plot((0:length(sub)-1) + mx,sub);subplot(2,1,2);
    % 
    %         flat = beat(x1:x2) - sub;
    %         flat = -1*flat;
    % 
    %         try
    %             [B, N, ~] = RunLength(flat>=0);
    %             if B(1) == 0
    %                 xtmp = x1+N(1)-1;
    %                 [~,i] = min(diff(beat(xtmp:x2)));
    %                 xtmp = i+xtmp-1;
    %                 [~,tmp]=findpeaks(diff(beat(xtmp:x2)),'MinPeakDistance',25);
    %                 tmp(2) = tmp(1);
    %                 tmp(1) = 1;
    %                 x1=xtmp;
    %             else
    %                 [~,tmp]=findpeaks(-diff(beat(x1:x2)),'MinPeakDistance',25);
    % 
    %             end
    % 
    %             b = beat(x1+tmp(1)-1:x1+tmp(2)-1);
    % 
    %             arr = 1:length(b);
    %             slope = (b(end) - b(1)) / (arr(end) - arr(1));
    % 
    %             sub = zeros(1,length(b));
    % 
    %             for q = 1:length(b)
    %                 sub(q) = slope * q + y1;
    %             end
    % 
    %             f = b - sub;
    %             f = -1*f;
    % 
    % 
    %             [~,l] = max(f);
    %             dicrotic_notch_idx = x1+tmp(1)-1+l;
    % 
    %         catch
    %             loc = simpleHill(flat);
    %             dicrotic_notch_idx = x1 + loc - 1;
    %         end
    %         
    %     end
%     catch
%         dicrotic_notch_idx= nan; 
%     end

if isempty(dicrotic_notch_idx)
    dicrotic_notch_idx=nan;
end
%     plot(f);vline(tmp(1));vline(tmp(2));
%     subplot(2,1,1);vline(dicrotic_notch_idx);
end