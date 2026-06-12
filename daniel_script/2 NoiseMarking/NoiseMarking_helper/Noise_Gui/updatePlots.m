function handles = updatePlots(handles)
    %% ecg Plot
    beginInd = round(handles.startSeg * handles.ecgSampleingRate) + 1;
    endInd = round(handles.endSeg * handles.ecgSampleingRate) + 1;
    
    if endInd > length(handles.ecg)
        size = endInd-beginInd;
        endInd = length(handles.ecg);
        beginInd = endInd - size;     
        
    elseif beginInd < 1
        size = endInd-beginInd;
        beginInd = 1;
        endInd = beginInd + size;
        
    end
    set(handles.ekg_axes, 'NextPlot', 'replacechildren');
    p1 = plot(handles.ekg_axes, seconds(handles.ecgTime_sec(beginInd:endInd)), handles.ecg(beginInd:endInd), 'DurationTickFormat', 'hh:mm:ss', 'Color', handles.ColorScheme(4, :));
    set(p1, 'HitTest', 'off');
    ecgXMIN = seconds(handles.ppgTime_sec(beginInd));
    ecgXMAX = seconds(handles.ppgTime_sec(endInd));

    [ecgYMIN, ecgYMAX] = getMinMaxY('ecg', handles.data_y_scale, handles.ecg(beginInd:endInd), handles.ecg_mean, handles.ecg_std, handles.ecg_peak_median, handles.ecg_peak_std);
    text(handles.ekg_axes, (ecgXMIN + ecgXMAX) / 2, ecgYMAX - ecgYMAX / 5, 'ECG', 'HorizontalAlignment', 'center');
    set(handles.ekg_axes, 'xlim', [ecgXMIN ecgXMAX], 'ylim', [ecgYMIN ecgYMAX]);

    
    %% PO Plot
    set(handles.ppg_axes, 'NextPlot', 'replacechildren');
    beginInd = round(handles.startSeg *handles.ppgSamplingRate) + 1;
    endInd = round(handles.endSeg *handles.ppgSamplingRate) + 1;
    
    
    if endInd > length(handles.ecg)
        size = endInd-beginInd;
        endInd = length(handles.ecg);
        beginInd = endInd - size;        
    elseif beginInd < 1
        size = endInd-beginInd;
        beginInd = 1;
        endInd = beginInd + size;
        
    end
    
    p2 = plot(handles.ppg_axes, seconds(handles.ppgTime_sec(beginInd:endInd)), handles.ppg(beginInd:endInd), 'DurationTickFormat', 'hh:mm:ss', 'Color', handles.ColorScheme(3, :));
    set(p2, 'HitTest', 'off');
    ppgXMIN = seconds(handles.ppgTime_sec(beginInd));
    ppgXMAX = seconds(handles.ppgTime_sec(endInd));
   
    [ppgYMIN, ppgYMAX] = getMinMaxY('ppg',handles.data_y_scale, handles.ppg(beginInd:endInd), handles.ppg_mean, handles.ppg_std);
    text(handles.ppg_axes, (ppgXMIN + ppgXMAX) / 2, ppgYMAX - ppgYMAX / 5, 'PPG', 'HorizontalAlignment', 'center');
    set(handles.ppg_axes, 'xlim', [ppgXMIN ppgXMAX], 'ylim', [ppgYMIN ppgYMAX]);

    Highregions = checkExclusions(handles);

    if (~isempty(Highregions))

        for I = 1:length(Highregions)
            span = Highregions{I};
            hold(handles.ppg_axes, 'on');
            hold(handles.ekg_axes, 'on');
            beg = span(1);
            last = span(2);
            if (span(3) == 'l')
                d1 = patch(handles.ppg_axes, [beg last last beg], [ppgYMAX ppgYMAX ppgYMIN ppgYMIN], handles.ColorScheme(3, :), 'LineStyle', 'none');
            elseif (span(3) == 'r')
                d1 = patch(handles.ekg_axes, [beg last last beg], [ecgYMAX ecgYMAX ecgYMIN ecgYMIN], handles.ColorScheme(4, :), 'LineStyle', 'none');
            else
                d1 = patch(handles.ppg_axes, [beg last last beg], [ppgYMAX ppgYMAX ppgYMIN ppgYMIN], 'r', 'LineStyle', 'none');
                alpha(d1, .15);
                d1 = patch(handles.ekg_axes, [beg last last beg], [ecgYMAX ecgYMAX ecgYMIN ecgYMIN], 'r', 'LineStyle', 'none');
                alpha(d1, .15);
            end

            alpha(d1, .15);
            hold(handles.ppg_axes, 'off');
            hold(handles.ekg_axes, 'off');
        end

    end

end

function Highregions = checkExclusions(handles)
    spans = handles.GenExc.ind;
    beginInd = handles.startSeg;
    endInd = handles.endSeg;

    Highregions = {};
    numregions = 1;

    for B = 1:size(spans, 1)
        startSpan = spans(B, 1);
        endSpan = spans(B, 2);

        if (endSpan > beginInd && endSpan < endInd)

            if (startSpan > beginInd)
                Highregions{numregions} = [spans(B, :) handles.GenExc.last(B)];
                numregions = numregions + 1;
            else
                Highregions{numregions} = [beginInd spans(B, 2) handles.GenExc.last(B)];
                numregions = numregions + 1;
            end

        elseif (startSpan < endInd && startSpan > beginInd)
            Highregions{numregions} = [spans(B, 1) endInd handles.GenExc.last(B)];
        elseif (startSpan < beginInd && endSpan > endInd)
            Highregions{numregions} = [beginInd endInd handles.GenExc.last(B)];
        end

    end

end

function [YMIN,YMAX] = getMinMaxY(type, checkbox, data, zero, std, zero2, std2)
    switch checkbox
        case 0 % fixed
            if type == 'ecg'
                YMIN = zero - std * 2;
                YMAX = zero2 + std2;
            else
                YMIN = zero - std * 3;
                YMAX = zero + std * 3;
            end

        case 1 % Size to data
            YMIN = min(data);
            YMAX = max(data);
            dif = abs(YMAX - YMIN) / 4; 
            YMIN = YMIN - dif;
            YMAX = YMAX + dif;
    end
    if YMIN >= YMAX
        YMIN = -1;
        YMAX = 1;
    end
end
