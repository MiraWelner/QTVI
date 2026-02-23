function [data] = FindWaveBounds_EKGandPPG(annealedSegments, dbg_plot, use_R_algorithms)

    %     if dbg_plot == 1
    %         figures = cell(length(annealedSegments),1);
    %     end

    % find rr index and associated poidx.
    data = cell(length(annealedSegments), 1);

    for i = 1:length(annealedSegments)
        %disp(['Section ' num2str(i) ' of '  num2str(length(annealedSegments))]);
        ppgSamplingRate = annealedSegments{i}.ppgSampleRate;
        ecgSamplingRate = annealedSegments{i}.ecgSampleRate;
        ecgSeg = annealedSegments{i}.ecg';
        rIsNoise = false;

        try
            ecgRIndex = annealedSegments{i}.r_peaks;
        catch
            ecgRIndex = [];
        end
        ppgSeg = annealedSegments{i}.po';
        try
        	[ppgMinAmps, ppgMaxAmps] = SegmentPPG(ppgSeg, ppgSamplingRate);
            data{i}.bad_segment = 0;
        catch
            ppgMinAmps = [];
            ppgMaxAmps = [];
            data{i}.bad_segment = 1;
        end
        %p = pspectrum(ecgSeg,ecgSamplingRate,'spectrogram','FrequencyResolution',.1,'Reassign',true);
        
%         if median(sum(p)) < 0.05
%             rIsNoise = true;
%         end
        if isempty(ecgRIndex) && use_R_algorithms
            if std(ecgSeg) == 0
                rIsNoise = true;
            end
        end
        
        
        if ~rIsNoise && isempty(ecgRIndex) && use_R_algorithms
            try
                ecgRIndex = JoinedRR(ecgSeg, ecgSamplingRate, 2);
                if length(ecgRIndex) < length(ppgMinAmps)/2 || length(ppgMinAmps) * 1.5 < length(ecgRIndex)
                    rIsNoise = true;
                end
            catch
               rIsNoise = true; 
            end
        end
%         if rIsNoise
%             disp(['Noisy ECG bin: ' num2str(i)]);
%         end

        if ~rIsNoise && ~isempty(ecgRIndex)
            try
                pairs = pairRtoPPGBeat(ecgSeg, ppgSeg, ecgSamplingRate, ppgSamplingRate, ecgRIndex, ppgMinAmps);
            catch
                if ~data{i}.bad_segment
                    ecgRIndex = [];
                    pairs = nan(length(ppgMinAmps),4);
                    pairs(:,1) = ppgMinAmps;
                    pairs(:,2) = -1;
                else
                    pairs = [];
                end
            end
        else
            if ~data{i}.bad_segment
                ecgRIndex = [];
                pairs = nan(length(ppgMinAmps),4);
                pairs(:,1) = ppgMinAmps;
                pairs(:,2) = -1;
            else
                pairs = [];
            end
        end
        
        segs.ecgSeg = ecgSeg;
        segs.ppgSeg = ppgSeg;
        segs.ecgRIndex = ecgRIndex;
        segs.ppgMinAmps = ppgMinAmps;
        segs.ppgMaxAmps = ppgMaxAmps;
        segs.pairs = pairs;

        data{i} = segs;
        data{i}.index = i;
        data{i}.ecgSamplingRate = ecgSamplingRate;
        data{i}.ppgSamplingRate = ppgSamplingRate;
        data{i}.ppg_bin_indexs = annealedSegments{i}.ppg_bin_indexs;
        data{i}.ecg_bin_indexs = annealedSegments{i}.ecg_bin_indexs;
        %
        %         if dbg_plot == 1
        %             subplot = @(m,n,p) subtightplot (m, n, p, [0 0]);
        %             %fig = figure('visible','off','Name',['RR index and PPG index | #' num2str(i) ' of ' num2str(length(annealedSegments))]);
        %             ax(1) = subplot(2,1,1);
        %             time = (0:length(segs.ecgSeg)-1)/ecgSamplingRate/60; plot(time,segs.ecgSeg);hold on; scatter(time(segs.ecgRIndex),segs.ecgSeg(segs.ecgRIndex));hold off;% plot(time(r),segs.ecgSeg(r),'o');
        %             %xlabel('Time (minutes)');
        %             ylabel('mV');
        %             title(['ECG vs ECG RR index from rrextract | #' num2str(i) ' of ' num2str(length(annealedSegments))]);
        %             grid;
        %             ax(2) = subplot(2,1,2);
        %             time = (0:length(segs.ppgSeg)-1)/ppgSamplingRate/60; plot(time,segs.ppgSeg);
        %             hold on;
        %             for z = 1:size(segs.ppgMinAmps,1)
        %                 line([time(segs.ppgMinAmps(z,1)) time(segs.ppgMinAmps(z,2))], [segs.ppgSeg(segs.ppgMinAmps(z,1)) segs.ppgSeg(segs.ppgMinAmps(z,2))], 'Color','c');
        %             end
        %             %scatter(time(segs.ppgMinAmps(:,1)),segs.ppgSeg(segs.ppgMinAmps(:,1)),'og');
        %             %scatter(time(segs.ppgMinAmps(:,2)),segs.ppgSeg(segs.ppgMinAmps(:,2)),'or');
        %             hold off;
        %             xlabel('Time (minutes)');
        %             ylabel('mV');
        %             %title(['PPG vs PPG index from inverted peak based on heartrate | #' num2str(i) ' of ' num2str(length(annealedSegments))]);
        %             linkaxes(ax,'x');
        %
        %             %figures{i} = fig;
        %             grid;
        %         end
    end

    %     if dbg_plot == 1
    %         ShowDbgPlots(figures);
    %     end

end

