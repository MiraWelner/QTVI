function handles = updateHighlightPlots(handles, nonnew)
    %% ppg amp plot
    if (~exist('nonnew'))
        nonnew = 0;
    end

    [handles.SetPlot(1), handles.patchesECG] = UpdateAmpogramPlot(handles.ekg_amp_axes, ...
        handles.ecg_ampogram, ...
        handles.ampogramtime, ...
        handles.time, ...
        handles.ColorScheme(1, :), ...
        handles.GenExc, ...
        nonnew, ...
        handles.patchesECG, ...
        handles.SetPlot(1));

    [handles.SetPlot(2), handles.patchesPPG] = UpdateAmpogramPlot(handles.ppg_amp_axes, ...
        handles.ppg_ampogram, ...
        handles.ampogramtime, ...
        handles.time, ...
        handles.ColorScheme(2, :), ...
        handles.GenExc, ...
        nonnew, ...
        handles.patchesPPG, ...
        handles.SetPlot(2));
    
    PlotSleepStages(handles.sleep_state_axes, handles.sleep_stage_time_hours, handles.sleepStages)
    
    for i = 1:length(handles.scoredEvents, 

end

function [setPlot, patches] = UpdateAmpogramPlot(axes, data, amp_time, time, color, genExc, nonnew, patches, setPlot)
    set(axes, 'NextPlot', 'replacechildren');

    % %% ecg ampogram
    min_ecg_ampogram = min(data);
    max_ecg_ampogram = max(data);
    min_ecg_ampogram_buffer = min_ecg_ampogram;
    max_ecg_ampogram_buffer = max_ecg_ampogram + max_ecg_ampogram / 5;

    if (setPlot == 0)
        p2 = plot(axes, (time(amp_time) / 60) / 60, data, 'Color', color);
        set(p2, 'HitTest', 'off');
        set(axes, 'xlim', [0, (time(end) / 60) / 60]);
        set(axes, 'YLim', [min_ecg_ampogram_buffer, max_ecg_ampogram_buffer]);
        text(axes, ((time(end) / 60) / 60) / 2, max_ecg_ampogram - max_ecg_ampogram / 5, 'ECG Amp-O-Gram', 'HorizontalAlignment', 'center');
        hold(axes, 'on');
        setPlot = 1;
        hold(axes, 'off');
        patches = [];
    else

        if (~isempty(genExc.ind) &&~nonnew)
            hold(axes, 'on');
            beg = genExc.ind(end, 1) / 60/60;
            last = genExc.ind(end, 2) / 60/60;
            upper = min_ecg_ampogram_buffer;
            lower = max_ecg_ampogram_buffer;

            if (genExc.last(end) == 'l')
                d2 = patch(axes, [beg last last beg], [upper upper lower lower], 'r', 'LineStyle', 'none');
            else
                d2 = patch(axes, [beg last last beg], [upper upper lower lower], 'g', 'LineStyle', 'none');
            end

            alpha(d2, .35);
            set(d2, 'HitTest', 'off');
            patches(end + 1) = d2;
            hold(axes, 'off');
        else
            p2 = plot(axes, (time(amp_time) / 60) / 60, data, 'Color', color);
            set(p2, 'HitTest', 'off');
            set(axes, 'xlim', [0, (time(end) / 60) / 60]);
            set(axes, 'YLim', [min_ecg_ampogram, max_ecg_ampogram]);
            hold(axes, 'on');
            text(axes, ((time(end) / 60) / 60) / 2, max_ecg_ampogram - max_ecg_ampogram / 5, 'ECG Amp-O-Gram', 'HorizontalAlignment', 'center');
            setPlot = 1;
            patches = [];

            for I = 1:size(genExc.ind, 1)
                beg = genExc.ind(I, 1) / 60/60;
                last = genExc.ind(I, 2) / 60/60;
                upper = min_ecg_ampogram_buffer;
                lower = max_ecg_ampogram_buffer;

                if (genExc.last(I) == 'l')
                    d2 = patch(axes, [beg last last beg], [upper upper lower lower], 'r', 'LineStyle', 'none');
                else
                    d2 = patch(axes, [beg last last beg], [upper upper lower lower], 'g', 'LineStyle', 'none');
                end

                alpha(d2, .35);
                set(d2, 'HitTest', 'off');
                patches(end + 1) = d2;
            end

            hold(axes, 'off');
        end

    end

end