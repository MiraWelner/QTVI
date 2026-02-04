function [dnotch] = template_plot(axis,ppg,ppgsamplingrate,ecg,ecgsamplingrate,alignment_point)
    lines = cell(4,1);
    hold on; 
    dnotch = nan;
    if ~isempty(ecg)
        yyaxis right;
        time = (0:length(ecg)-1);
        if abs(alignment_point)/ecgsamplingrate > .5
        else
            time = time-alignment_point;
        end
%         time = time/ecgsamplingrate;
        plot(axis,time,ecg,'Color',rgb('gray'),'Parent',axis,'HitTest','off');
    
        
    end
    
    if ~isempty(ppg)
        yyaxis left;
       % ppg_tmp = [nan(1, alignment_point) ppg];
%         time = (0:length(ppg)-1)/ppgsamplingrate;
        time = (0:length(ppg)-1);
        p = plot(axis,time,ppg,'Color',rgb('red'),'Parent',axis,'HitTest','off');
        
%         lines{1} = line([alignment_point alignment_point],[-10 10],'Parent',axis,'HitTest','off','Tag', 'Onset');
%         set(lines{1},'Color',[0 0 1],'LineWidth',1);
        
%         [~,max_idx] = max(ppg);
%         lines{2} = line([max_idx max_idx],[-10 10],'Parent',axis,'HitTest','off','Tag', 'Peak');
%         set(lines{2},'Color',[1 0 0],'LineWidth',1);
        
        dnotch = dumbDicrotic(ppg);
%         dnotch = dnotch/ppgsamplingrate;
       
        lines{3} = line([dnotch dnotch],[min(ppg) max(ppg)],'Parent',axis,'HitTest','off','Tag', 'Dicrotic');
        set(lines{3},'Color',rgb('blue'),'LineWidth',1);
        
%         idx = dumbEnd(ecg,dnotch, ppg, alignment_point );
%         lines{4} = line([idx idx],[-10 10],'Parent',axis,'HitTest','off','Tag', 'End');
%         set(lines{4},'Color',[128/255 0 128/255],'LineWidth',1);
        
        if min(ppg)~=max(ppg)
            ylim([min(ppg),max(ppg)]);
        end
    end
    hold off;
end

function [end_idx] = dumbEnd(ecg,dnotch,ppg,alignment_point)
    if isempty(ecg) || sum(isnan(ecg))
        end_idx = length(ppg)-alignment_point;
        return
    end

    try
        ecg = ecg(dnotch:end);
        [~,end_idx] = max(ecg);
        end_idx = dnotch + end_idx;
    catch
        end_idx = backup;
    end
    if end_idx > length(ppg)
        end_idx = length(ppg)-alignment_point;
    end
    if isempty(end_idx)
        end_idx = length(ppg)-alignment_point;
    end
        
end
