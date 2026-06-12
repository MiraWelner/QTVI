function PlotScoredEvents(axes, time, lbled_events)
    %set(axes, 'NextPlot', 'replacechildren');

    cla(axes);
    for i = 1:length(lbled_events)
        e = lbled_events{i};
        line(axes,[e.Start_sec/3600 e.Start_sec/3600 + e.Duration_sec/3600], [0,0], 'LineWidth', 5, 'HitTest', 'off');
    end
    
    xlim(axes,[min(time)/3600 max(time)/3600])
end