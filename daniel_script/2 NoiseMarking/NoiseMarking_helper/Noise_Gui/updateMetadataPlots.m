function handles = updateMetadataPlots(handles, plot_type)
    if plot_type == 0
        InitPlots(handles);
    else
        if plot_type == 1
            InitPlots(handles);
        end
               
        UpdatePatches(handles.ekg_amp_axes, handles.GenExc.ind, plot_type, handles);
        UpdatePatches(handles.ppg_amp_axes, handles.GenExc.ind, plot_type, handles);
        UpdatePatches(handles.sleep_state_axes, handles.GenExc.ind, plot_type, handles);
        UpdatePatches(handles.events_axes, handles.GenExc.ind, plot_type, handles);
    end

end

function InitPlots(handles)
    PlotSleepStages(handles.sleep_state_axes, handles.time_sec, handles.sleepStages_timeIndex, handles.sleepStages);

    PlotScoredEvents(handles.events_axes, handles.time_sec, handles.scoredEvents);

    AmpogramPlot(handles.ekg_amp_axes, handles.ecgAmpogram, handles.time_sec, ...
        handles.ecgAmpogram_timeIndex, handles.ColorScheme(1, :), 'ECG Amp-O-Gram');

    AmpogramPlot(handles.ppg_amp_axes, handles.ppgAmpogram, handles.time_sec, ...
        handles.ppgAmpogram_timeIndex, handles.ColorScheme(2, :), 'PPG Amp-O-Gram');
end


function AmpogramPlot(axes, data, time_sec, time_idx, color, name)
    spacing_divisor = 6;

    set(axes, 'NextPlot', 'replacechildren');
    
    time_hrs = time_sec / 3600;

    min_ampogram = min(data);
    max_ampogram = max(data);
    min_ampogram_buffer = min_ampogram;
    max_ampogram_buffer = max_ampogram + max_ampogram / spacing_divisor;
    
    p = plot(axes, time_hrs(time_idx), data, 'Color', color, 'HitTest', 'off');
    set(axes, 'Xlim', [0, time_hrs(end)]);
    set(axes, 'YLim', [min_ampogram_buffer, max_ampogram_buffer]);
    text(axes, time_hrs(end) / 2, max_ampogram - max_ampogram / spacing_divisor, name, 'HorizontalAlignment', 'center', 'HitTest', 'off');
end

function UpdatePatches(axes, ind, plot_type, handles)
    if ~isempty(handles.GenExc.last)
        if (handles.GenExc.last(end) == 'l')
            color = handles.ColorScheme(3, :);
        elseif (handles.GenExc.last(end) == 'c')
            color = 'r';
        else
            color = handles.ColorScheme(4, :);
        end

        if plot_type == 1 % undo
            for i = 1:size(ind, 1)
                
                if (handles.GenExc.last(i) == 'l')
                    color = handles.ColorScheme(3, :);
                elseif (handles.GenExc.last(i) == 'c')
                    color = 'r';
                else
                    color = handles.ColorScheme(4, :);
                end

                beg = ind(i, 1) / 3600;
                last = ind(i, 2) / 3600;
                p = patch(axes, [beg last last beg], [axes.YLim(2) axes.YLim(2) axes.YLim(1) axes.YLim(1)], color, 'LineStyle', 'none', 'HitTest', 'off');
                alpha(p, .35);
            end
        else % just a new patch
            if (handles.reviewing)
               color =  [255/256 223/256 0/256];
            end
            beg = ind(end, 1) / 3600;
            last = ind(end, 2) / 3600;
            p = patch(axes, [beg last last beg], [axes.YLim(2) axes.YLim(2) axes.YLim(1) axes.YLim(1)], color, 'LineStyle', 'none', 'HitTest', 'off');
            alpha(p, .35);
        end
    end
end