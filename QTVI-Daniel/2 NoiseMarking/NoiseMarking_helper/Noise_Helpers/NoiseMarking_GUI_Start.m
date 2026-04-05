% pulse annotate script
clear;
close all;


version = '0.1.0';
disp(['Version: ' version]);
% [analysisFiles, process_noise_manual, outputLoc] = PulseSetup('normal', '');
props = readProps('config.txt');
NoiseMarking_input = props('NoiseMarking_input'); 
NoiseMarking_output = props('NoiseMarking_output'); 
Skip_Existing = logical(str2num(props('Skip_Existing')));

analysisFiles = dir(fullfile(NoiseMarking_input, '**/*_ecg_ppg_sleep.mat'));

% [analysisFiles, process_noise_manual, outputLoc] = PulseSetup(edfFile, props);
disp('*********************************************************************');
completed = {};
for i = 1:size(analysisFiles, 1)

%     check_header = edfread(analysisFiles{i, 2}); % 1 is ecg index, 23 pulse
%     
%     if ~isequal(find(strcmp(check_header.label, 'Pleth')), 23) ||~isequal(find(strcmp(check_header.label, 'EKG')), 1)
%         cprintf('e rr', ['Error code #99: ', analysisFiles{i, 1}, ' marked as unusual (pleth/ekg not in proper place in .edf) continuing on to next file.' newline]);
%         disp(['____________________________________________________________________________________________________' newline]);
%         continue
%     end
% 
%     [edf_hdr, edf_data] = edfread(analysisFiles{i, 2}, 'verbose', 1, 'targetSignals', [1, 23]); % 1 is ecg index, 23 pulse
% 
%     ppg = edf_data(2, :);
%     ecg = edf_data(1, :);
%     ecgSamplingRate = edf_hdr.frequency(1); % 256 FOR MESA SET
%     ppgSamplingRate = edf_hdr.frequency(23); % 256 FOR MESA SET
% 
%     disp('Reading XML...');
%     try
%         [scoring_epoch_size_sec, sleepStages, scoredEvents] = ReadXML(analysisFiles{i, 3});
%     catch
%         cprintf('err', ['Error code #100: ', analysisFiles{i, 1}, '''s .xml could not be read. Marking error and continuing on to next file.']);
%         disp(['____________________________________________________________________________________________________' newline]);
%         continue
%     end
    name = analysisFiles(i).name;
    start_idx = regexp(name, '_ecg_ppg_sleep');
    start_idx = start_idx(1);
    name = name(1:start_idx-1);
    disp(join(['Beginning analysis of ' analysisFiles(i).name]));

    if isfile(fullfile(NoiseMarking_output, [name '_noise_markings.mat'])) && Skip_Existing
        disp([analysisFiles(i).name ' exists skipping because Skip_Existing = 1 in config.']);
        continue
    end
    
    load(fullfile(analysisFiles(i).folder, analysisFiles(i).name));
    if std(ppg) == 0
        disp([analysisFiles(i).name ' has only flat line ppg. Skipping']);
        continue
    end
    tStart = tic;
    disp('Generating prelim data...');
%     [sleepStages] = RenumberSleepStages(sleepStages);

    %prelim analysis
%     poSmooth = smoothdata(ppg, 'sgolay', round(ppgSamplingRate));
%     [syst, systLoc] = findpeaks(poSmooth, ppgSamplingRate);

    ppgTime_sec = (0:length(ppg)-1)/ppgSamplingRate;
    ecgTime_sec = (0:length(ecg)-1)/ecgSamplingRate;
%     ppgTime_sec = 0:1 / ppgSamplingRate:(length(poSmooth) / ppgSamplingRate - 1 / ppgSamplingRate);
%     ecgTime_sec = 0:1 / ecgSamplingRate:(length(ecg) / ecgSamplingRate - 1 / ecgSamplingRate);
    sleepStageTime_hrs = (0:length(sleepStages) - 1) / (1 / scoring_epoch_size_sec) / 3600;

    % ampograms
    window_size_mins = 1;

    ecgAmpogram = windowedMinMaxDiff(ecg, ecgSamplingRate, window_size_mins);
    ecgAmpogram_time_seconds = ((0:length(ecgAmpogram) - 1) * (60 * window_size_mins));
    [~, ~, ecgAmpogram_timeIndex] = intersect(ecgAmpogram_time_seconds, ecgTime_sec);

    ppgAmpogram = windowedMinMaxDiff(ppg, ppgSamplingRate, window_size_mins);
    ppgAmpogram_time_seconds = ((0:length(ppgAmpogram) - 1) * (60 * window_size_mins));
    [~, ~, ppgAmpogram_timeIndex] = intersect(ppgAmpogram_time_seconds, ppgTime_sec);
    
    
%     plot(ppgTime_sec,zeros(length(ppgTime_sec),1),'.'); hold on;
%     plot(ppgAmpogram_time_seconds,zeros(length(ppgAmpogram_time_seconds),1),'o')


    %% call gui for user to prune segs
    inputs.name = analysisFiles(i).name;

    inputs.ppg = ppg;
    inputs.ppgTime_sec = ppgTime_sec;
    inputs.ppgSamplingRate = ppgSamplingRate;
    inputs.ppgAmpogram = ppgAmpogram;
    inputs.ppgAmpogram_timeIndex = ppgAmpogram_timeIndex;

    inputs.ecg = ecg;
    inputs.ecgTime_sec = ecgTime_sec;
    inputs.ecgSamplingRate = ecgSamplingRate;
    inputs.ecgAmpogram = ecgAmpogram;
    inputs.ecgAmpogram_timeIndex = ecgAmpogram_timeIndex;

    inputs.sleepStageTime_hrs = sleepStageTime_hrs;
    inputs.sleepStageSampleRate = 1/scoring_epoch_size_sec;
    inputs.sleepStages = sleepStages;
    if exist('scoredEvents','var') == 0
        inputs.scoredEvents = {};
    else
        inputs.scoredEvents = scoredEvents;
    end
    inputs.rIndex = DeepR_and_PVC(ecg,ecgSamplingRate);
    inputs.reviewing = 0;
%     if process_noise_manual == -10
%         inputs.reviewing = 1;
%         inputs.noiseInfo = load(analysisFiles{i, 4});
%     end
    
    disp('Showing noise GUI...');

    prune_response = ReviewData(inputs);

    % if we closed and didn't click finalize then quit
    if prune_response.closeResponse == 0
        disp('Bye!');
        disp('Completed this session:');
        disp(completed);
        return;
    end
    noiseSEG = prune_response.GenExc.noiseExc;
    marktype = prune_response.GenExc.last;

    noise_markings = zeros(size(noiseSEG,1),5);
    for x = 1:size(noiseSEG,1)
        noise_markings(x,1) = noiseSEG(x,1);
        noise_markings(x,2) = noiseSEG(x,2);
        noise_markings(x,3) = noiseSEG(x,1)/ecgSamplingRate; % only works atm if sampling rates are same
        noise_markings(x,4) = noiseSEG(x,2)/ecgSamplingRate;
        noise_markings(x,5) = marktype(x);
    end

    data_description = 'First column and second column are times marked in seconds. Third and fourth are first and second column values divided by sampling rate, for closest actual sample round number to integer. The last column is the ascii letter designation of button clicked to mark. In this file l (108) = ppg noise | r(114)  = ecg noise | c(99) = both have noise';
    disp('Saving markings...');

    save(fullfile(NoiseMarking_output, [name '_noise_markings.mat']), 'noise_markings', 'data_description', 'version');
    completed{end+1} = analysisFiles(i).name;
    disp(['____________________________________________________________________________________________________' newline]);
    toc(tStart);

    disp('Completed this session:');
    disp(completed);
    
    
    
end

x=1;
