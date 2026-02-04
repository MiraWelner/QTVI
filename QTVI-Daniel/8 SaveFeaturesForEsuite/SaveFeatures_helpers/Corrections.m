% I believe I ran the feature generation with half finished 'Sleep state'
% code this corrects the files:
% *-adjusted_sleep_state.mat
% *-sec_to_first_onset_of_sleep.mat
% *-sec_from_last_onset_of_sleep.mat
clear;
clc;

patientFeatureDirectory = "/hdd/data/mesa/Manuscript 1/2020 November run/9 individual features";
mat_list = dirWithoutDots(fullfile(patientFeatureDirectory));

names = cellfun(@(x) x(1:16), {mat_list.name}, 'UniformOutput', false);

uniq_names = unique(names);

for x = 1:length(uniq_names)
    clearvars -except uniq_names x patientFeatureDirectory;

    disp(['Fixing sleep files for ' uniq_names{x}]);
    
    % load 4 files to be fixed. loaded into flattened struct so I can copy
    % function over to correct original code and not have to add it there
    flattened = [];
    sleep_stages = load(fullfile(patientFeatureDirectory, uniq_names{x},[uniq_names{x} '-sleep_stages.mat']));
    flattened.sleep_stages = sleep_stages.tmp;
%     adjusted_sleep_state = load(fullfile(patientFeatureDirectory, uniq_names{x},[uniq_names{x} '-adjusted_sleep_state.mat']));
%     flattened.adjusted_sleep_state = adjusted_sleep_state.tmp;
%     sec_to_first_onset_of_sleep = load(fullfile(patientFeatureDirectory, uniq_names{x},[uniq_names{x} '-sec_to_first_onset_of_sleep.mat']));
%     flattened.sec_to_first_onset_of_sleep = sec_to_first_onset_of_sleep.tmp;
%     sec_from_last_onset_of_sleep = load(fullfile(patientFeatureDirectory, uniq_names{x},[uniq_names{x} '-sec_from_last_onset_of_sleep.mat']));
%     flattened.sec_from_last_onset_of_sleep = sec_from_last_onset_of_sleep.tmp;
    
    

    corrected_time_sec = load(fullfile(patientFeatureDirectory, uniq_names{x},[uniq_names{x} '-corrected_time_msec.mat']));
    flattened.corrected_time_sec = corrected_time_sec.tmp;

    [B(:,1), B(:,2), ~] = RunLength(flattened.sleep_stages);
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
    sleep_states = flattened.sleep_stages;
    if after_sleep_count > 0
        for i = length(sleep_states)-after_sleep_count+1:length(sleep_states)
            sleep_states(i) = 4;
        end
    end
    
    % Set all REM to 2
    sleep_states(sleep_states==-1) = 2;%2=REM
    
    % Set all NREM to 1
    sleep_states(sleep_states ~= 0 & sleep_states ~= 2 & sleep_states ~= 4) = 1;
    
    % Set all instances of awake to 3 so that awake during sleep is
    % labeled.
    sleep_states(sleep_states == 0) = 3;
    
    % Correct wake before sleep back to 0
    for i = 1:before_sleep_wake_count
        sleep_states(i) = 0;
    end
    
    flattened.adjusted_sleep_state =  sleep_states;

    flattened.sec_to_first_onset_of_sleep = (flattened.corrected_time_sec - flattened.corrected_time_sec(before_sleep_wake_count)) * -1;
    flattened.sec_from_last_onset_of_sleep = (flattened.corrected_time_sec - flattened.corrected_time_sec(end - after_sleep_count));
    
    tmp = flattened.adjusted_sleep_state;
    save(fullfile(patientFeatureDirectory, uniq_names{x},[uniq_names{x} '-adjusted_sleep_state.mat']),'tmp');
    
    tmp = flattened.sec_to_first_onset_of_sleep;
    save(fullfile(patientFeatureDirectory, uniq_names{x},[uniq_names{x} '-sec_to_first_onset_of_sleep.mat']),'tmp');
    
    tmp = flattened.sec_from_last_onset_of_sleep;
    save(fullfile(patientFeatureDirectory, uniq_names{x},[uniq_names{x} '-sec_from_last_onset_of_sleep.mat']),'tmp');

    
    tmp = flattened.corrected_time_sec;
    save(fullfile(patientFeatureDirectory, uniq_names{x},[uniq_names{x} '-corrected_time_sec.mat']),'tmp');
end




