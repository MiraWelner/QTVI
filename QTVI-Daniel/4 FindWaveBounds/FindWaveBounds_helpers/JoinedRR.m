function [rr] = JoinedRR(ecgSeg,ecgSamplingRate,diff_range)
    if std(ecgSeg) == 0
        return
    end
    
    algos = {{@rpeakdetect, {ecgSeg,ecgSamplingRate}}, ...
            {@rpeakdetect, {ecgSeg,ecgSamplingRate,0.1}}, ...
            {@rpeakdetect, {ecgSeg,ecgSamplingRate,0.4}},...
            {@pan_tompkin, {ecgSeg,ecgSamplingRate}}, ...
            {@ecgLms, {ecgSeg - mean(ecgSeg),ecgSamplingRate, 5,12,0}}};
            %{@SimpleRST,{ecgSeg,ecgSamplingRate}},...
            %{@getECGFeatures,{ecgSeg,ecgSamplingRate}},...
    %weights = [0.75,0.25,0.25,1.25,1.5,0.75,1]; % with wavelet

    weights = [0.75,0.25,0.25,1.25,1.5,0.75];
    weight_max = sum(weights);
        
    output = cell(length(algos),1);
    parfor i=1:length(algos)
        curr_alg = algos{i};
        try
            if length(curr_alg{2}) == 2
                output{i} = curr_alg{1}(curr_alg{2}{1}, curr_alg{2}{2});
            elseif length(curr_alg{2}) == 3
                output{i} = curr_alg{1}(curr_alg{2}{1},curr_alg{2}{2},curr_alg{2}{3});
            else
                output{i} = curr_alg{1}(curr_alg{2}{1},curr_alg{2}{2},curr_alg{2}{3},curr_alg{2}{4},curr_alg{2}{5});
            end
        catch
            output{i} = [];
%             disp(':(');
        end
    end
    
    potentialPeaks = sortedList(output,weights);
    median_amp = median(ecgSeg(potentialPeaks(:,1)));
    median_dists = zeros(length(output),1);
    for i = 1:length(output)
        median_dists(i) = nanmedian(diff(output{i}));
    end
    median_dist = nanmedian(median_dists);
    output{end+1} = RRsimpleSquared(ecgSeg,median_dist/2);

    
    
    parfor r = 4:6
        output{r} = RPeakfromRWave(ecgSeg,output{r});
    end
%     output{end+1} = RRWavelet(ecgSeg,median_dist/2);

         
   
    potentialPeaks = sortedList(output,weights);
    
    % shift peaks which are within range of specificed diff to largest value in
    % that range
    uniq = unique(potentialPeaks(:,1));
    mask = diff(uniq)<=diff_range;
    for i = 1:length(mask)
        if mask(i) == 1
            if ecgSeg(uniq(i)) > ecgSeg(uniq(i+1))
               potentialPeaks(potentialPeaks(:,1)==uniq(i+1),1) = uniq(i);
            else
               potentialPeaks(potentialPeaks(:,1)==uniq(i),1) = uniq(i+1);
            end
        end
    end
        
    uniq = unique(potentialPeaks(:,1));
    weighted_peaks = zeros(length(uniq),2);
    for i = 1:length(uniq)
        weighted_peaks(i,1) = uniq(i);
        weighted_peaks(i,2) = sum(potentialPeaks(potentialPeaks(:,1)==uniq(i),2));
    end
    
    rr = weighted_peaks(weighted_peaks(:,2) >= 2.4,1);
    
    
%     
%     ecgR_outliers = stdoutlier(rr', 2.5, round(length(rr)*.35), 'lower',0);
%     
%     if sum(ecgR_outliers) >1
%             close all;
% 
%             plot(ecgSeg);   hold on;
%             plot(rr(ecgR_outliers),ecgSeg(rr(ecgR_outliers)),'*r')
%             plot(rr(~ecgR_outliers),ecgSeg(rr(~ecgR_outliers)),'*g')
% 
%     end
 

            
%             
%             
%             outs = {};
%             [B, N, BI] = RunLength(ecgR_outliers);
%             
%             if ~isempty(N(N>=2 & B == 1))
%                 idxs = find((N>=2 & B == 1) == 1);
%                 for q = 1:length(idxs)
%                     idx1 = BI(idxs(q));
%                     idx2 = BI(idxs(q))+N(idxs(q))-1;
%                     section = ecgRIndex(idx1:idx2);
%                     out = bestFitOfBeatTrain(ecgRIndex(1:idx1-1), ecgRIndex(idx2+1:end), section, 15);
%                     outs{end+1} = [ecgRIndex(1:idx1-1); out';ecgRIndex(idx2+1:end)];
%                 end
%             end
%             
% 
% 
%             for z = 1:length(outs)
%                 tmp = ismember(ecgRIndex,outs{z});
%                 ecgRIndex = ecgRIndex(tmp);
%             end
    

    
%     close all;
%     time = 1:length(ecgSeg);
%     ax(1) = subplot(2,1,1);
%     axis tight;
%     plot(time,ecgSeg);
%     hold on;
%     c = distinguishable_colors(length(output));
%     for i = 1:length(output)
%         scatter(time(output{i}),ecgSeg(output{i}),'o','MarkerEdgeColor','none','MarkerFaceColor', c(i,:),'MarkerFaceAlpha',.4);
%     end
%     
%     for r = 1:length(rr)
%         text(weighted_peaks(r,1)+10, ecgSeg(weighted_peaks(r,1)),num2str(weighted_peaks(r,2)), 'Clipping', 'on','FontSize',5);
%     end
%     
%     legend('wave','rpeakdetect-norm','rpeakdetect-low','rpeakdetect-high','tompkin','deeps','simpleRSquared','wavelet');
%     
%     ax(2) = subplot(2,1,2);
% 
%     axis tight;
%     plot(time,ecgSeg);
%     hold on;
%     plot(time(rr),ecgSeg(rr),'o');
%     linkaxes(ax,'x');
    

end

function [lst] = sortedList(output,weight)
    lens = cellfun(@length,output);
    lst = zeros(sum(lens),2);
    prev = 1;
    for i = 1:length(output)
        lst(prev:prev+lens(i)-1,1) = output{i};
        lst(prev:prev+lens(i)-1,2) = weight(i);
        prev = prev+lens(i);
    end
    lst = sortrows(lst,1);
end