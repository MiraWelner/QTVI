function handles = updateCursor(handles)

    if (handles.setup)
        for i = 1:length(handles.cursors)
        	delete(handles.cursors{i});
        end
    else
        handles.setup = 1;
        handles.cursors = {};
    end

    view_end_in_hours = handles.endSeg / 3600;
    view_begin_in_hours = handles.startSeg / 3600;

    if handles.viewWidth <=30
       alpha_value = 1;
    elseif handles.viewWidth == 60
       alpha_value = .75;
    else
       alpha_value = .5;
    end
    
    handles.cursors = UpdateCursor(handles.ppg_amp_axes, handles.cursors, alpha_value, view_end_in_hours, view_begin_in_hours);
    handles.cursors = UpdateCursor(handles.ekg_amp_axes, handles.cursors, alpha_value, view_end_in_hours, view_begin_in_hours);
    handles.cursors = UpdateCursor(handles.sleep_state_axes, handles.cursors, alpha_value, view_end_in_hours, view_begin_in_hours);
    handles.cursors = UpdateCursor(handles.events_axes, handles.cursors, alpha_value, view_end_in_hours, view_begin_in_hours);   

end

function [cursors] = UpdateCursor(axes, cursors, alpha_value, view_end_in_hours, view_begin_in_hours)
    color = [132/255 142/255 255/255]; 
    hold(axes, 'on');
    y_min = axes.YLim(1);
    y_max = axes.YLim(2);
    cursors{end+1} = patch(axes, [view_begin_in_hours view_end_in_hours view_end_in_hours view_begin_in_hours], [y_min y_min y_max y_max], color, 'LineStyle', 'none', 'HitTest', 'off');
    alpha(cursors{end}, alpha_value);
    hold(axes, 'off');

end
