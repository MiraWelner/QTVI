function handles = updateWidth(handles)
    width = 0;

    switch handles.window_time_selection.Value
        case 1
            width = 10;
        case 2
            width = 30;
        case 3
            width = 60;
        case 4
            width = 300;
        case 5
            width = 600;
    end

    if (handles.viewWidth ~= width)
        handles.viewWidth = width;
        handles.CurrentPosition - width / 2;
        handles.startSeg = max([0, handles.CurrentPosition - width / 2]);
        handles.endSeg = min([handles.ppgTime_sec(end), handles.CurrentPosition + width / 2]);

        if (handles.startSeg == 0)
            handles.CurrentPosition = width / 2;
            handles.endSeg = width;

        elseif (handles.endSeg == handles.ppgTime_sec(end))
            handles.CurrentPosition = handles.ppgTime_sec(end) - width / 2;
            handles.startSeg = max([0, handles.CurrentPosition - width / 2]);

        end

        handles = updatePlots(handles);
        handles = updateScroll(handles);
    end
