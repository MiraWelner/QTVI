% deep data:  MESA_ECG_results_20190720_stata13_v3
% clear;
close all;
analysisFiles = dir(fullfile("D:\deepLab\central data\MESA\MATS", '**/*.mat'));
% T = readtable('C:\Users\dan\OneDrive - University of Cincinnati\projects\000001-ecg_ppg_coupling\data\ecg\MESA_ECG_results_20190720_stata13_v3_times_only.csv');
figure('units','normalized','outerposition',[0 0 1 1]);
output = "C:\Users\dan\OneDrive - University of Cincinnati\projects\000001-ecg_ppg_coupling\data\processed_data\sleep_state_comparison_vs_qtvi";
for i = 1:length(analysisFiles)
    disp([num2str(i) ' of ' num2str(length(analysisFiles))]);
    load(fullfile(analysisFiles(i).folder, analysisFiles(i).name)); 
    disp(join(['Plotting Sleep for: ' analysisFiles(i).name]));
    id = str2num(analysisFiles(i).name(1:7));
    tmp = T(T.idno == id,:);
    tmp = sortrows(tmp,'ns');
    mintime = min(tmp.time);
    maxtime = max(tmp.time);
    

    
    axis = subplot(2,1,1);
    
    sleepStageTime_hrs = (0:length(sleepStages) - 1) / (1 / scoring_epoch_size_sec) / 3600;

    sleepStageSampleRate = 1/scoring_epoch_size_sec;
    
    max_length = numel(ecg);
    max_samplerate = ecgSamplingRate;
    time_sec = (0:max_length - 1) / max_samplerate;
    sleepStages_timeIndex = 1:(1 / sleepStageSampleRate * max_samplerate):numel(time_sec);
    sleepStages_timeIndex = sleepStages_timeIndex(1:length(sleepStages));
    clear B
    [B(:,1), B(:,2), ~] = RunLength(sleepStages);
    sleepIDX = find(B(:,1) == 0,1);
    before_sleep_wake_count = B(sleepIDX,2);
    
    wakeIDX = find(B(:,1) ~= 0 & ~isnan(B(:,1)),1,'last') + 1;
    if wakeIDX > size(B,1)
        after_sleep_count = 0;
    else
        after_sleep_count = B(wakeIDX,2);
    end

    % sleep states go:
    % 0: Awake
    % -1: REM
    % -2 to -5:  NREM 1-5
    % nan: sleepstate unknown
    %
    % We want to group nrems as 1 and make 3 new classes with 5 total 'sleep states':
    % 0: awake before sleep
    % 1: NREM 1-4
    % 2: REM
    % 3: awake during sleep
    % 4: awake after sleep
    
    % 4=Awake after sleep without sleeping again. Set everything after last
    % wake to 4
    collapsed_sleep_states = sleepStages;
    if after_sleep_count > 0
        for x = length(collapsed_sleep_states)-after_sleep_count+1:length(collapsed_sleep_states)
            collapsed_sleep_states(x) = 4;
        end
    end
    
    % Set all REM to 2
    collapsed_sleep_states(collapsed_sleep_states==-1) = 2;%2=REM
    
    % Set all NREM to 1
    collapsed_sleep_states(collapsed_sleep_states ~= 0 & collapsed_sleep_states ~= 2 & collapsed_sleep_states ~= 4) = 1;
    
    % Set all instances of awake to 3 so that awake during sleep is
    % labeled.
    collapsed_sleep_states(collapsed_sleep_states == 0) = 3;
    
    % Correct wake before sleep back to 0
    for x = 1:before_sleep_wake_count
        collapsed_sleep_states(x) = 0;
    end 
    collapsed_sleep_states = collapsed_sleep_states*-1;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% first plot
    
    time_hours = time_sec/3600;
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
    sgtitle(['ID: ', num2str(id)]);

    title('Original sleep states from edf');
    
    index_of_last  = find(collapsed_sleep_states == 0, 1, 'last');
    time_hours_sleep_adjusted = time_hours - time_hours(index_of_last);
    mask = time_hours_sleep_adjusted >= mintime & time_hours_sleep_adjusted <= maxtime;
    if sum(mask) == 0
        str = 'ID not in ECG dataset.';
        unknown = text(axis,0.45,0.9,str,'Units','Normalized' );
        warning('off','all');
        legend(axis,[wake,rem,s1,s2,s3,s4], {'Wake', 'REM', 'NREM 1', 'NREM 2' , 'NREM 3' , 'NREM 4'},'Location','southwest','NumColumns',2,'AutoUpdate', 'off', 'HitTest', 'off');
        legend(axis,'boxoff');
        warning('on','all');
    else
        os = ones(sum(mask),1);
        unknown = plot(axis,time_hours(mask),os,'s','MarkerSize',8,'MarkerEdgeColor','none','MarkerFaceColor', [1 0 0], 'HitTest', 'off'); % unknown magenta
        warning('off','all');
        legend(axis,[wake,rem,s1,s2,s3,s4,unknown], {'Wake', 'REM', 'NREM 1', 'NREM 2' , 'NREM 3' , 'NREM 4', 'Stat ECG Time'},'Location','southwest','NumColumns',2,'AutoUpdate', 'off', 'HitTest', 'off');
        legend(axis,'boxoff');
        warning('on','all');
    end
    % suppress warning for extras


    ylim(axis,[min(uniq)-1 , 2]);
    xlim(axis,[0,time_sec(end)/3600]);
    yticks(axis,[]);
%     set(axis, 'xticklabel', []);
    zeroidx = closest_idx(tmp.time,0);

    adj = tmp.time - min(tmp.time);
    t = time_hours(closest_idx(time_hours,adj(zeroidx)));
    idx = closest_idx(time_hours,t);
    
    vline(time_hours(index_of_last),{'color',rgb('gray'),'linestyle','--'});
    if ~isnan(adj)
        vline(time_hours(idx),{'color',rgb('red'),'linestyle','--'});
    end
    hold(axis,'off');
    a=1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%second plot   
    axis1 = subplot(2,1,2);
    adjusted_stair = collapsed_sleep_states;
    adjusted_stair(adjusted_stair==-3 | adjusted_stair==-4) = 0;  adjusted_stair(adjusted_stair==-2) = 1; adjusted_stair(adjusted_stair==-1) = -2;  adjusted_stair(adjusted_stair==1) = -1;
    stair = stairs(axis1, time_hours, adjusted_stair, 'Color', [230 230 230]/255);
    hold(axis1,'on');
    wakeb4sleep = plot(axis1, time_hours(collapsed_sleep_states ==  0),collapsed_sleep_states(collapsed_sleep_states ==  0),'s','MarkerSize',marker_size,'MarkerEdgeColor','none','MarkerFaceColor', [0 0 0], 'HitTest', 'off'); % wake black
    nrem = plot(axis1, time_hours(collapsed_sleep_states == -1),collapsed_sleep_states(collapsed_sleep_states == -1)-1,'s','MarkerSize',marker_size,'MarkerEdgeColor','none','MarkerFaceColor', [0 0 1], 'HitTest', 'off'); % rem green
    rem = plot(axis1,time_hours(collapsed_sleep_states == -2),collapsed_sleep_states(collapsed_sleep_states == -2)+1,'s','MarkerSize',marker_size,'MarkerEdgeColor','none','MarkerFaceColor', [0 1 0], 'HitTest', 'off'); % sleep 1 blue
    wakeduringsleep = plot(axis1,time_hours(collapsed_sleep_states == -3),collapsed_sleep_states(collapsed_sleep_states == -3)+3,'s','MarkerSize',marker_size,'MarkerEdgeColor','none','MarkerFaceColor', [128/255 128/255 128/255], 'HitTest', 'off'); % sleep 2 cyan
    wakeaftersleep = plot(axis1,time_hours(collapsed_sleep_states == -4),collapsed_sleep_states(collapsed_sleep_states == -4)+4,'s','MarkerSize',marker_size,'MarkerEdgeColor','none','MarkerFaceColor', [1 0 1], 'HitTest', 'off'); % sleep 3 red
    


    ylim(axis1,[min(uniq)-1 , 2]);
    xlim(axis1,[0,time_sec(end)/3600]);
    yticks(axis1,[]);
%     set(axis, 'xticklabel', []);

    linkaxes([axis,axis1],'x')
    title('Collapsed sleep states');
    if sum(mask) == 0
        str = 'ID not in ECG dataset.';
        unknown = text(axis1,0.45,0.9,str,'Units','Normalized' );
        warning('off','all');
        legend(axis1,[wakeb4sleep,nrem,rem,wakeduringsleep,wakeaftersleep], {'Wake Before Sleep (ss0)', 'NREM (ss1)', 'REM (ss2)', 'Wake During Sleep (ss3)' , 'Wake After Sleep (ss4)'},'Location','southwest','NumColumns',2,'AutoUpdate', 'off', 'HitTest', 'off');
        legend(axis1,'boxoff');
        warning('on','all');
    else
        os = ones(sum(mask),1);
        unknown = plot(axis1,time_hours(mask),os,'s','MarkerSize',8,'MarkerEdgeColor','none','MarkerFaceColor', [1 0 0], 'HitTest', 'off'); % unknown magenta
        warning('off','all');
        legend(axis1,[wakeb4sleep,nrem,rem,wakeduringsleep,wakeaftersleep,unknown],{'Wake Before Sleep (ss0)', 'NREM (ss1)', 'REM (ss2)', 'Wake During Sleep (ss3)' , 'Wake After Sleep (ss4)', 'Stat ECG Time'},'Location','southwest','NumColumns',2);
        legend(axis1,'boxoff');
        warning('on','all');
    end
    vline(time_hours(index_of_last),{'color',rgb('gray'),'linestyle','--'});
    if ~isnan(adj)
        vline(time_hours(idx),{'color',rgb('red'),'linestyle','--'});
    end
    hold(axis1,'off');
    
    x=1;

    saveas(gcf,fullfile(output,[analysisFiles(i).name(1:16) '.tif']));
    clf;
end