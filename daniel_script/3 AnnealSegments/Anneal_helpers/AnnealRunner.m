clear;
close all;
props = readProps('config.txt');
AnnealSegments_AnnealType = str2num(props('AnnealSegments_AnnealType'));

AnnealSegments_use_handmarked_noise = str2num(props('AnnealSegments_use_handmarked_noise'));
AnnealSegments_noise_input = props('AnnealSegments_noise_input');
AnnealSegments_raw_mats = props('AnnealSegments_raw_mats');
AnnealSegments_output_path = props('AnnealSegments_output_path');
AnnealSegments_segment_size_mins = str2num(props('AnnealSegments_segment_size_mins'));
AnnealSegments_window_size_mins = str2num(props('AnnealSegments_ampo_window_size_mins'));
AnnealSegments_noise_expansion_seconds = str2num(props('AnnealSegments_noise_expansion_seconds'));
Skip_Existing = logical(str2num(props('Skip_Existing')));

%%
analysisFiles = SetupAnneal(AnnealSegments_raw_mats, AnnealSegments_noise_input, AnnealSegments_use_handmarked_noise);
time = 0;

for i = 1:size(analysisFiles, 1)
    name = analysisFiles{i,1};

    if Skip_Existing && isfile(fullfile(AnnealSegments_output_path, [name '_annealedSegments.mat']))
    	disp([name '_annealedSegments.mat exists skipping because Skip_Existing = 1 in config.']);
        continue
    end
%     if strcmp(name, '7017189_20120523') | strcmp(name, '7017170_20110630') 
%     else
%         continue
%     end
%     if isfile(fullfile(output_path, [analysisFiles{i, 1} '_annealedSegments.mat']))
%         update = 0;
%         load(fullfile(output_path, [analysisFiles{i, 1} '_annealedSegments.mat']))
%         for x = 1:length(annealedSegments)
%             if ~isfield(annealedSegments{x},'ppgSampleRate') | ~isfield(annealedSegments{x},'ppg_bin_indexs')
%                 update=1;
%             end
%         end
%         if ~update 
%             continue
%         end
%     end
%for i = randperm(length(1:size(analysisFiles, 1)))
%     if isfile(fullfile(output_path, [analysisFiles{i, 1} '_annealedSegments.mat']))
%         continue
%     end
   
%     tStart = tic;

%     disp(join(['Beginning analysis of ' analysisFiles{i, 1} ' | ' num2str(i) ' of ' num2str(size(analysisFiles,1))]));
%     avg_time = time/i;
%     disp(['Avg Time (s): ' num2str(avg_time)]);

%     disp(['Est finish (min): ' num2str((avg_time*(size(analysisFiles,1)-i))/60)]);
%     disp(join(['Output loc ' AnnealSegments_output_path]));


    data = load(analysisFiles{i, 3}, 'ecg', 'ecgSamplingRate', 'ppg', 'ppgSamplingRate', 'scoring_epoch_size_sec', 'sleepStages');
    ecg = data.ecg;
    ecgSamplingRate = data.ecgSamplingRate;
    ppg = data.ppg;
    ppgSamplingRate = data.ppgSamplingRate;
    scoring_epoch_size_sec = data.scoring_epoch_size_sec;
    sleepStages = data.sleepStages;
    
    rs = [];
    if ppgSamplingRate ~= 1
        % Noise Segments are in seconds. First element is beginning of noise
        % second is end of noise
        if AnnealSegments_use_handmarked_noise == 1
            % load from manual noise markings
            noiseSEG = load(analysisFiles{i, 2});
            if isfield(noiseSEG,'noiseSEG')
                noiseSEG = noiseSEG.noiseSEG;
            else
                tmp = noiseSEG.noise_markings(:,1:2);
                noiseSEG = tmp; % TODO: allow user to specify using only specific noise marking classes
            end
        else
            if size(analysisFiles,2) > 3
                if ~isempty(analysisFiles{i, 4})
                    rs = load(analysisFiles{i, 4});
                    try
                        rs = rs.rs;
                    catch
                        rs = rs.newRs;
                    end
                    noiseClass = load(analysisFiles{i, 2});
                    noiseClass = noiseClass.noiseClass;

                    badidxs = zeros(size(rs,1),1);
                    for j = 1:size(noiseClass,1)
                        mask = (rs(:,2) >= noiseClass{j,1}) & (rs(:,2) <= noiseClass{j,2});
                        badidxs(mask==1) = 1;
                    end
                    badidxs(rs(:,3)==2)=1; % 2==marked noise

                    rs = rs(~badidxs,1); % only get indexes where not noise
                else
                    rs = [];
                end
            else
                rs = [];
            end

            % simplistic noise markings taking out some of the spikes and areas
            % where wave is flat
            ppgTime_sec = 0:1 / ppgSamplingRate:(length(ppg) / ppgSamplingRate - 1 / ppgSamplingRate);

            % ampograms
            ppgAmpogram = windowedMinMaxDiff(ppg, ppgSamplingRate, AnnealSegments_window_size_mins);
            ppgAmpogram_time_seconds = ((0:length(ppgAmpogram) - 1) * (60 * AnnealSegments_window_size_mins));
            [~, ~, ppgAmpogram_timeIndex] = intersect(ppgAmpogram_time_seconds, ppgTime_sec);

            ppg_std = movstd(ppg,round(ppgSamplingRate/2));

            spikes_begin=find(isoutlier(ppgAmpogram,'gesd')==1);

            spikes = zeros(length(ppg),1);
            n = numel(ppg);
            for x = 1:length(spikes_begin)
                start = ppgAmpogram_timeIndex(spikes_begin(x));
                add = (ppgSamplingRate*60*AnnealSegments_window_size_mins);
                if start+add > n
                   add = n-start; 
                end
                spikes(start:start+add) = 1;
            end
%             yyaxis left
%             plot(ppg)
%             yyaxis right
%             plot(spikes);

            s = OneDDbscan(spikes==1, round(AnnealSegments_noise_expansion_seconds*ppgSamplingRate));

            SpikeNoiseSEGs = zeros(size(s,1),2);
            for x = 1:size(s,1)
                if ppgTime_sec(s(x,1))-AnnealSegments_noise_expansion_seconds > 0
                    SpikeNoiseSEGs(x,1) = ppgTime_sec(s(x,1))-AnnealSegments_noise_expansion_seconds; 
                else
                    SpikeNoiseSEGs(x,1) = 0; 
                end

               if ppgTime_sec(s(x,2))+AnnealSegments_noise_expansion_seconds < ppgTime_sec(end)
                    SpikeNoiseSEGs(x,2) = ppgTime_sec(s(x,2))+AnnealSegments_noise_expansion_seconds; 
                else
                    SpikeNoiseSEGs(x,2) = ppgTime_sec(end); 
                end

            end



            f = OneDDbscan(ppg_std==0, round(length(ppg_std)*0.01));
            FlatNoiseSEGs = zeros(size(f,1),2);
            for x = 1:size(f,1)
                if ppgTime_sec(f(x,1))-AnnealSegments_noise_expansion_seconds > 0
                    FlatNoiseSEGs(x,1) = ppgTime_sec(f(x,1))-AnnealSegments_noise_expansion_seconds; 
                else
                    FlatNoiseSEGs(x,1) = 0; 
                end

               if ppgTime_sec(f(x,2))+AnnealSegments_noise_expansion_seconds < ppgTime_sec(end)
                    FlatNoiseSEGs(x,2) = ppgTime_sec(f(x,2))+AnnealSegments_noise_expansion_seconds; 
                else
                    FlatNoiseSEGs(x,2) = ppgTime_sec(end); 
                end
            end
            tempSegs = [SpikeNoiseSEGs; FlatNoiseSEGs];

            tempSegs = sortrows(tempSegs,1);
            begidx = 1;
            index=1;
            while 1
                if index > size(tempSegs,1)
                   break 
                end

                begidx = tempSegs(index,1);
                endidx = tempSegs(index,2);
                mask = (begidx <= tempSegs) & (tempSegs <= endidx);

                idxs = find(sum(mask,2)==2); % sections btw intervals
                idxs = idxs(idxs>index);
                for x = 1:length(idxs)                  
                   tempSegs(idxs(x),:) = []; % remove section
                   mask(idxs(x),:) = []; % remove section
                   idxs = idxs-1;

                end

                flag = 0;
                idxs = find(mask(:,1)==1); % sections overlapping above intervals
                idxs = idxs(idxs>index);
                for x = 1:length(idxs) 
                    if tempSegs(idxs(x),2) > endidx
                        tempSegs(index,2) = tempSegs(idxs(x),2); % replace current end w/ new end
                    end
                    tempSegs(idxs(x),:) = [];
                    idxs = idxs-1;

                    flag = 1;
                end 

                if flag == 0
                    index = index+1;                
                end

            end


            noiseSEG = tempSegs;
            
%             Fs = 256;            % Sampling frequency                    
%             T = 1/Fs;             % Sampling period       
% 
%             w=256*1000;
%             
%                         L = w;             % Length of signal
%             t = (0:L-1)*T;        % Time vector
%             n=length(ppg);
%             nw=n-w+1;  
%             k=(1:w);
%             for x = 1:nw
%                   Y = fft(ppg(k));
%                   P2 = abs(Y/L);
%                     P1 = P2(1:L/2+1);
%                     P1(2:end-1) = 2*P1(2:end-1);
% 
%                     f = Fs*(0:(L/2))/L;
%                     subplot(2,1,1)
%                 plot(f,P1) 
%                 title('Single-Sided Amplitude Spectrum of X(t)')
%                 xlabel('f (Hz)')
%                 ylabel('|P1(f)|')
%                 xlim([0,15]);
%                 subplot(2,1,2)
%                 
%                 plot(ppg(k));
%                 
%                 k=k+w;
%                 clf;
%             end
            
            
%             close all;
%             plot(ppgTime_sec,ppg);
%             hold on;
%             stairs(ppgTime_sec(ppgAmpogram_timeIndex),ppgAmpogram,'r');
%             plot(ppgTime_sec(ppg_std==0),ppg(ppg_std==0),'oc');
%     
%     
%             for x = 1:size(SpikeNoiseSEGs,1)
%                 d2 = patch(gca, [SpikeNoiseSEGs(x,1) SpikeNoiseSEGs(x,2) SpikeNoiseSEGs(x,2) SpikeNoiseSEGs(x,1)], [5 5 -5 -5], 'r', 'LineStyle', 'none');
%                 alpha(d2, .35);
%             end
%             vline(ppgTime_sec(s(:,1)),{'Color','g'});
%             vline(ppgTime_sec(s(:,2)),{'Color','r'});
%             for x = 1:size(FlatNoiseSEGs,1)
%                 d2 = patch(gca, [FlatNoiseSEGs(x,1) FlatNoiseSEGs(x,2) FlatNoiseSEGs(x,2) FlatNoiseSEGs(x,1)], [5 5 -5 -5], 'g', 'LineStyle', 'none');
%                 alpha(d2, .35);
%             end
%             
%             
%             for x = 1:size(noiseSEG,1)
%                 d2 = patch(gca, [noiseSEG(x,1) noiseSEG(x,2) noiseSEG(x,2) noiseSEG(x,1)], [5 5 -5 -5], 'r', 'LineStyle', 'none');
%                 alpha(d2, .35);
%             end
            
            
        end

        
%         interactiveMatrixProfileVer3(ppg(1:256*5*60),256)
%         load('omnitemplate.mat')

%         bin_size_samples = ppgSamplingRate * 60 * 5;
%         [bin_breaks, bin_count] = getBinBreaksAndCount(ppg, bin_size_samples, ppgSamplingRate, 0);
%         profiles = cell(bin_count,1);
%         for z = 1:size(bin_breaks,1)
%             tic
%             temp = ppg(bin_breaks(z,1):bin_breaks(z,2));
%             profiles{z} = MatrixProfile(temp, ppgSamplingRate*2, 50);
%             toc
%             
%             
%         end
%         
%         prev = 0;
%         figure;
%         
%         for x = 1:length(profiles)
%             time = prev:prev+length(profiles{x})-1;
%             prev = prev+length(profiles{x})+1;
%             
%             plot(time,profiles{x});
%             hold on;
%         end
%         
        
        
        if AnnealSegments_AnnealType == 0
            annealedSegments = SimpleSegments(ppg, ppgSamplingRate, ecg, ecgSamplingRate, noiseSEG, scoring_epoch_size_sec, sleepStages, AnnealSegments_segment_size_mins);
        else
%         disp('Anealing segments...');
        	[annealedSegments, final_bin_idx] = AnnealSegments(ppg, ppgSamplingRate, ecg, ecgSamplingRate, noiseSEG, scoring_epoch_size_sec, sleepStages, AnnealSegments_segment_size_mins, rs, 1);
        end
%         for x = 1:length(annealedSegments)
% %            if std(annealedSegments{x}.po) == 0 || numel(annealedSegments{x}.po) < 256
% %               annealedSegments(x) = [];
% %            end
%         end
        parsave(annealedSegments,fullfile(AnnealSegments_output_path, [analysisFiles{i, 1} '_annealedSegments']));
        %movefile(analysisFiles{i, 2},finished_path);
    end

%     disp(['____________________________________________________________________________________________________' newline]);
%     time = time + toc(tStart);
%     toc(tStart);

end

function parsave(annealedSegments,path)
    save(path,'annealedSegments');
end
