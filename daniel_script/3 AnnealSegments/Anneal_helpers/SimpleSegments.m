function [annealedSegments] = SimpleSegments(ppg, ppgSampleRate, ecg, ecgSampleRate, noiseSEG, scoring_epoch_size_sec, sleepStages, targetLength)
    %targetLength given in minutes
    
    %% determine ideal bin times
    %use 1:length instead of 0 for ecg to make calculating easier later.
    ecg_time_seconds = ((0:length(ecg)-1) / ecgSampleRate)';
    po_time_seconds = ((0:length(ppg) - 1) / ppgSampleRate)';
    bin_size_samples = ppgSampleRate * 60 * targetLength;
    [bin_breaks, bin_count] = getBinBreaksAndCount(ppg, bin_size_samples, ppgSampleRate, 0);   
    
    % sleep data labels do not always match the length of the data. when
    % this is the case fill the extra time with unknown sleepstate.
    if scoring_epoch_size_sec * ppgSampleRate * length(sleepStages) < length(ppg)
        d = length(ppg) - scoring_epoch_size_sec * ppgSampleRate * length(sleepStages);
        addNum = ceil(d/(ppgSampleRate*scoring_epoch_size_sec));
        sleepStages = [sleepStages; nan(addNum,1)];
    end
    sleep_stages_propigated = RunLength(sleepStages, ones(length(sleepStages),1) * scoring_epoch_size_sec * ppgSampleRate);
    
    %fill data from index's
    annealedSegments = cell(numel(bin_count),1);
    for i = 1:bin_count
        beg_idx = bin_breaks(i,1);
        end_idx = bin_breaks(i,2);
        time = beg_idx:end_idx;
%         plot(time,ppg(beg_idx:end_idx))
%         annealedSegments{i}.index = i;

        annealedSegments{i}.ppg_bin_indexs = bin_breaks(i,:);
        annealedSegments{i}.ecg_bin_indexs = bin_breaks(i,:);
        annealedSegments{i}.ppgSampleRate = ppgSampleRate;
        annealedSegments{i}.ecgSampleRate = ecgSampleRate;          

        annealedSegments{i}.sleep_stages = sleep_stages_propigated(beg_idx:end_idx);
        annealedSegments{i}.po = ppg(beg_idx:end_idx);
        annealedSegments{i}.ecg = ecg(beg_idx:end_idx);
        
        mask = noiseSEG >= (beg_idx-1)/ppgSampleRate & noiseSEG <= (end_idx-1)/ppgSampleRate;
        annealedSegments{i}.ppg_noise_segs = [];
        if sum(mask(:,1))>=1 || sum(mask(:,2))>=1
           for x = 1:size(mask,1)
               if mask(x,1) == 1 && mask(x,2)==0
                   annealedSegments{i}.ppg_noise_segs = [annealedSegments{i}.ppg_noise_segs; [floor(noiseSEG(x,1)*ppgSampleRate) end_idx]];
               elseif mask(x,1) == 0 && mask(x,2)==1
                   annealedSegments{i}.ppg_noise_segs = [annealedSegments{i}.ppg_noise_segs; [beg_idx floor(noiseSEG(x,2)*ppgSampleRate)]];
               elseif mask(x,1) == 1 && mask(x,2)==1
                   annealedSegments{i}.ppg_noise_segs = [annealedSegments{i}.ppg_noise_segs; [floor(noiseSEG(x,1)*ppgSampleRate) floor(noiseSEG(x,2)*ppgSampleRate)]];
               end
           end
        end
    end
end