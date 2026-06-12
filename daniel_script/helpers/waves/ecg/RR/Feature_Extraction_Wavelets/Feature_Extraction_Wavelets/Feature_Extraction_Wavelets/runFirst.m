%Copyright 2017 The MathWorks, Inc.
%% Feature Extraction and Detection - Part 2
% Goal: : Isolate signal components by decomposing the signal into
% different frequency bands using Wavelet Transform


%% Steps
% Load and visualize the signal 
% Detect peaks in the signal using findpeaks function
% Decompose the signal into various subbands using Discrete Wavelet Transform
% Reconstruct the signal components with subbands that represent the QRS waveform


%% Load signals

load ecgSig; 

%% Visualize ECG signal and reference annotations
figure
hr = computeHeartRate(tm(expertAnns)); 
titleStr = sprintf('Heart Rate : %0.2f bpm (Reference)',hr);
helperTimeDomain(tm,ecg_Preprocessed,titleStr,[], 'b');
hold on
l1 = plot(tm(expertAnns),ecg_Preprocessed(expertAnns),'ro'); 
legend(l1,'R peaks -> Reference'); 
xlim([1,30]);


%% Detect R-peaks using a peak detection approach
ecgSigSq =ecg_Preprocessed.^2;
[ecgSqPks,ecgSqlocs] = findpeaks(ecgSigSq,'MinPeakHeight',0.5,...
     'MinPeakDistance',60); 
h1 = computeHeartRate(tm(ecgSqlocs)); 
titleStr = sprintf('Heart rate is :: %0.2f bpm, Reference is :%0.2f bpm',h1,hr);
 
figure; 
plot(tm,ecg_Preprocessed); 
hold on; 
h1 = plot(tm(expertAnns),ecg_Preprocessed(expertAnns),'ro'); 
hold on; 
h2 = plot(tm(ecgSqlocs),ecg_Preprocessed(ecgSqlocs),'k*'); 
legend([h1,h2],'R-wave Reference','R-wave from findpeaks'); axis tight; grid on;
xlim([0, 20]);
ylim([-3,3]);
title(titleStr);


%% 
% <Open qrsLims.bmp>
% Center of the QRS ranges from 5 Hz to 22 Hz
% We will extract this using the Discrete Wavelet Transform

%% %% Partition the signal into different bands

% Level 4 Details ==> 11 - 22 Hz 
% Level 5 Details ==> 5 Hz- 11 Hz
ecgQRSBand = extractQRSband(ecg_Preprocessed,5,4,5);
figure;
helperTimeDomain(tm,ecgQRSBand, 'Isolating QRS Wave (Level 4 and 5)',30,'b')

%% Visualize the results
ecgQRSBandSq = ecgQRSBand.^2;
[ecgQRSBandSqPks,ecgQRSBandSqlocs] = findpeaks(ecgQRSBandSq,'MinPeakHeight',0.35,...
    'MinPeakDistance',60); 
hw = computeHeartRate(tm(ecgQRSBandSqlocs)); 
fprintf('Avg. Heart rate obtained after Wavelet QRS band isolation:: %0.2f bpm, Reference value:%0.2f bpm \n',hw,hr);
figure; 
plot(tm,ecg_Preprocessed); 
hold on; 
h1 = plot(tm(expertAnns),ecg_Preprocessed(expertAnns),'ro'); 
hold on; 
h2 = plot(tm(ecgQRSBandSqlocs),ecg_Preprocessed(ecgQRSBandSqlocs),'k*'); 
title(sprintf('Reference heart rate:         %0.2f bpm \t\t\t\nWavelet based heart rate: %0.2f bpm',hr,hw));
legend([h1,h2],'R-wave Reference ', 'R-wave detected using wavelets'); axis tight; grid on;
xlim([0,100]);
ylim([-5,5]);
