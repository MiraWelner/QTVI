function PlotSleepStages(axis, time, sleepStages_timeIndex, sleepStages)
    %set(axes, 'NextPlot', 'replacechildren');

    time_hours = time/3600;
    time_hours = time_hours(sleepStages_timeIndex);
    marker_size = 4;
    uniq = unique(sleepStages);

    stair = stairs(axis, time_hours, sleepStages, 'Color', [230 230 230]/255);
    hold(axis,'on');
    wake = plot(axis, time_hours(sleepStages ==  0),sleepStages(sleepStages ==  0),'s','MarkerSize',marker_size,'MarkerEdgeColor','none','MarkerFaceColor', [0 0 0], 'HitTest', 'off'); % wake black
    rem = plot(axis, time_hours(sleepStages == -1),sleepStages(sleepStages == -1),'s','MarkerSize',marker_size,'MarkerEdgeColor','none','MarkerFaceColor', [0 1 0], 'HitTest', 'off'); % rem green
    s1 = plot(axis,time_hours(sleepStages == -2),sleepStages(sleepStages == -2),'s','MarkerSize',marker_size,'MarkerEdgeColor','none','MarkerFaceColor', [0 0 1], 'HitTest', 'off'); % sleep 1 blue
    s2 = plot(axis,time_hours(sleepStages == -3),sleepStages(sleepStages == -3),'s','MarkerSize',marker_size,'MarkerEdgeColor','none','MarkerFaceColor', [0 1 1], 'HitTest', 'off'); % sleep 2 cyan
    s3 = plot(axis,time_hours(sleepStages == -4),sleepStages(sleepStages == -4),'s','MarkerSize',marker_size,'MarkerEdgeColor','none','MarkerFaceColor', [1 0 0], 'HitTest', 'off'); % sleep 3 red
    s4 = plot(axis,time_hours(sleepStages == -5),sleepStages(sleepStages == -5),'s','MarkerSize',marker_size,'MarkerEdgeColor','none','MarkerFaceColor', [1 .5 .25], 'HitTest', 'off'); % sleep 4
    unknown = plot(axis,time_hours(sleepStages == -6),sleepStages(sleepStages == -6),'s','MarkerSize',marker_size,'MarkerEdgeColor','none','MarkerFaceColor', [1 0 1], 'HitTest', 'off'); % unknown magenta
%     tmp = time_hours >= 1.875 & time_hours <= 9.708333;
%     os = ones(sum(tmp),1);
%     unknown = plot(axis,time_hours(time_hours >= 1.875 & time_hours <= 9.708333),os*-3,'s','MarkerSize',marker_size,'MarkerEdgeColor','none','MarkerFaceColor', [1 0 1], 'HitTest', 'off'); % unknown magenta
    
    % suppress warning for extras
    warning('off','all');
    legend(axis,[wake,rem,s1,s2,s3,s4,unknown], {'Wake', 'REM', 'NREM 1', 'NREM 2' , 'NREM 3' , 'NREM 4', 'Unknown'},'Location','southwest','NumColumns',2,'AutoUpdate', 'off', 'HitTest', 'off');
    legend(axis,'boxoff');
    warning('on','all');

    ylim(axis,[min(uniq)-.5 , 1]);
    xlim(axis,[0,time(end)/3600]);
    yticks(axis,[]);
    set(axis, 'xticklabel', []);
    hold(axis,'off');

end