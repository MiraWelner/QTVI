    
function [ret] = GenerateFeatures(anneal_path, wave_path, template_path, output_path)
        [~, name, ~] = fileparts(anneal_path);
        name = regexp(name, '(\d+_\d+)', 'match');
        name = name{1};
% try
        window_length = 30;
        % value between 1 and 100 lower the value the more "different" a template vs beat can be. Set to inf to process all. 
        sqi_threshold_to_process = inf;
        
%         [~, name, ~] = fileparts(anneal_path);
%         name = regexp(name, '(\d+_\d+)', 'match');
%         name = name{1};
        
        annealedSegments = load(anneal_path);
        annealedSegments = annealedSegments.annealedSegments;
        
%         updated = 0;
%         for x = 1:length(annealedSegments)
%             if ~isfield(annealedSegments,'ppgSampleRate')
%                 annealedSegments{x}.ppgSampleRate = 256;
%                 updated=1;
%             end
%             
%             if ~isfield(annealedSegments,'ecgSampleRate')
%                 updated=1;
%                 annealedSegments{x}.ecgSampleRate = 256;
%             end
%             sleep_stages = annealedSegments{x}.sleep_stages;
%             if length(sleep_stages) < length(annealedSegments{x}.po)
%                 
%                 updated=1;
%             end
%         end
%         if updated
%             save(anneal_path,'annealedSegments');
%         end

        wave_data = load(wave_path);
        wave_data = wave_data.wave_data;

        if strcmpi(template_path,'')
            disp('Finding beat features...');
            bin_marks_mask = countBins(wave_data);           
        else
            template_info = load(template_path);
            template_info = template_info.template_info;

            disp('Finding beat features from templates...');
            
            if length(annealedSegments) ~= length(template_info)
                updated = 0;
                if length(annealedSegments) - length(template_info) == 1
                    if length(annealedSegments{1}.po) <= 258
                        tmp = [];
                        tmp.ecgTemplate = [];
                        tmp.ppgTemplate = [];
                        tmp.alignment_point = nan;
                        tmp.TemplateBad = 1;
                        tmp.bad_r_templates = 1;
                        tmp.bad_ppg_templates = 1;
                        tmp.Dicrotic = nan;
                        tmp.Onset = nan;
                        tmp.Peak = nan;
                        tmp.End = nan;
                        hold = template_info;
                        template_info = cell(length(hold)+1,1);
                        template_info{1} = tmp;
                        for x = 1:length(hold)
                            template_info{x+1} = hold{x};
                        end
                        save(template_path,'template_info'); 
                        updated = 1;
                    elseif length(annealedSegments{end}.po) <= 258
                        tmp = [];
                        tmp.ecgTemplate = [];
                        tmp.ppgTemplate = [];
                        tmp.alignment_point = nan;
                        tmp.TemplateBad = 1;
                        tmp.bad_r_templates = 1;
                        tmp.bad_ppg_templates = 1;
                        tmp.Dicrotic = nan;
                        tmp.Onset = nan;
                        tmp.Peak = nan;
                        tmp.End = nan;
                        hold = template_info;
                        template_info = cell(length(hold)+1,1);
                        template_info{end} = tmp;
                        for x = 1:length(hold)
                            template_info{x} = hold{x};
                        end
                        save(template_path,'template_info'); 
                        updated = 1;                   
                    end
                end
                if updated == 0
                    try
                        errors = load('errors.mat');
                    catch
                        errors.names = {};
                    end
                    errors.names{end+1} = name;
                    names = errors.names;
                    save('errors.mat','names')
                    ret= 0;
                    return
                end
            end
            
            bin_marks_mask = countBins(wave_data, template_info);
        end

        
        beats_per_bin_individual = cell(length(wave_data), 1);
        beats_in_bin = cell(length(wave_data), 1);

        for t = 1:length(wave_data)
            disp(['Section ' num2str(t) ' of '  num2str(length(wave_data))]);
            if(exist('template_info', 'var') )
                % don't change to double && needs to be single
%                 if bin_marks_mask.could_not_identify_PPG(t) | bin_marks_mask.ppg_template_manually_excluded(t)
%                     beats_in_bin{t}.sqi = -1;
%                     beats_in_bin{t}.error_ppg_segmentation = bin_marks_mask.could_not_identify_PPG(t);
%                     beats_in_bin{t}.review_bad_ppg_template = bin_marks_mask.ppg_template_manually_excluded(t);
%                     beats_in_bin{t}.review_bad_r_template = bin_marks_mask.ecg_template_manually_excluded(t);
%                     continue
%                 end
%                     error('not imp');
%                 else
                    try
                        ppgidxs = wave_data{t}.pairs(:,1)';
                    catch
                        ppgidxs = [];
                    end
                    [sqi, sqilabels] = PPG_SQI(annealedSegments{t}.po', ppgidxs, template_info{t}.ppgTemplate, window_length*annealedSegments{t}.ppgSampleRate, annealedSegments{t}.ppgSampleRate);
                    size_peak_2_end = template_info{t}.End - template_info{t}.Peak;
                    if isnan(size_peak_2_end) | size_peak_2_end < 1 % KEEP AS SINGLE DOUBLE FAILS
                        [beats_per_bin_individual{t}, beats_in_bin{t}] = GetBeatFeaturesFromTemplate(sqi(:,1),sqi_threshold_to_process,annealedSegments{t}.po',annealedSegments{t}.sleep_stages,wave_data{t}.pairs, annealedSegments{t}.ppgSampleRate,annealedSegments{t}.ecgSampleRate);
                    else
                        ratio_sp = abs(template_info{t}.Dicrotic-template_info{t}.Peak)/size_peak_2_end;
                        [beats_per_bin_individual{t}, beats_in_bin{t}] = GetBeatFeaturesFromTemplate(sqi(:,1),sqi_threshold_to_process,annealedSegments{t}.po',annealedSegments{t}.sleep_stages,wave_data{t}.pairs, annealedSegments{t}.ppgSampleRate,annealedSegments{t}.ecgSampleRate,ratio_sp);

                    end
                    
                    beats_in_bin{t}.sqi = sqi(sqi(:,1) < sqi_threshold_to_process,:);
                    beats_in_bin{t}.sqilabels = sqilabels;
                    beats_in_bin{t}.error_ppg_segmentation = bin_marks_mask.could_not_identify_PPG(t);
                    beats_in_bin{t}.review_bad_ppg_template = bin_marks_mask.ppg_template_manually_excluded(t);
                    beats_in_bin{t}.review_bad_r_template = bin_marks_mask.ecg_template_manually_excluded(t);
%                 end
            else
                beats_in_bin{t}.sqi = ones(length(wave_data{t}.ppgMinAmps)-1)*-1; % set everything to unknown acceptance
                beats_in_bin{t}.sqilabels = {'not_run'};
                [beats_per_bin_individual{t}, beats_in_bin{t}] = GetBeatFeaturesFromTemplate(zeros(size(wave_data{t}.pairs,1)-1,1),sqi_threshold_to_process,annealedSegments{t}.po',annealedSegments{t}.sleep_stages,wave_data{t}.pairs, annealedSegments{t}.ppgSampleRate,annealedSegments{t}.ecgSampleRate);
                beats_in_bin{t}.error_ppg_segmentation = bin_marks_mask.could_not_identify_PPG(t);
                beats_in_bin{t}.review_bad_ppg_template = bin_marks_mask.ppg_template_manually_excluded(t);
                beats_in_bin{t}.review_bad_r_template = bin_marks_mask.ecg_template_manually_excluded(t);
            end
        end
        disp('Flattening...');
        [beats_flattened] = flatten_beat_idx(beats_in_bin, annealedSegments);
                
        lens = cellfun(@(x) length(x.po), annealedSegments);
        ppg_wout_noise = zeros(sum(lens), 1);
        prev = 1;

        for t = 1:length(annealedSegments)
            len = length(annealedSegments{t}.po) - 1;
            ppg_wout_noise(prev:prev + len) = annealedSegments{t}.po;
            prev = prev + len;
        end

        beats_flattened.ppg_wout_noise = ppg_wout_noise;
%         plot(ppg_wout_noise);hold on;plot(beats_flattened.idx_systolic(~isnan(beats_flattened.idx_systolic)),ppg_wout_noise(beats_flattened.idx_systolic(~isnan(beats_flattened.idx_systolic))),'oc');
%         vline(length(annealedSegments{1}.po));
        disp('Saving...');
        save(fullfile(output_path, [name '_feature_output']), 'beats_flattened');

        ret = 1;
% catch e
%     try
%         errors = load('fatal.mat');
%     catch
%         errors.names = {};
%     end
%     errors.names{end+1} = name;
%     names = errors.names;
%     save('fatal.mat','names')
%     ret = 0;
%     return
% end

    
end

function goodBins = countBins(waveData,template)
    if exist('template','var')
        could_not_identify_PPG = logical(cellfun(@(x) int8(x.bad_segment),waveData));
        poor_ppg_or_ecg_template_manually_excluded = logical(cellfun(@(x) int8(x.TemplateBad),template));
        ecg_template_manually_excluded = logical(cellfun(@(x) int8(x.bad_r_templates),template));
        ppg_template_manually_excluded = logical(cellfun(@(x) int8(x.bad_ppg_templates),template));

%         c = logical(cellfun(@(x) length(x.ppgTemplate(~isnan(x.ppgTemplate))) > 1,template));
        goodBins.could_not_identify_PPG = could_not_identify_PPG;
        goodBins.poor_ppg_or_ecg_template_manually_excluded = poor_ppg_or_ecg_template_manually_excluded;
        goodBins.ecg_template_manually_excluded = ecg_template_manually_excluded;
        goodBins.ppg_template_manually_excluded = ppg_template_manually_excluded;

    else
       could_not_identify_PPG = logical(cellfun(@(x) int8(x.bad_segment),waveData));
       goodBins.could_not_identify_PPG = could_not_identify_PPG;
       goodBins.poor_ppg_or_ecg_template_manually_excluded = zeros(size(could_not_identify_PPG,1),1);
       goodBins.ecg_template_manually_excluded = zeros(size(could_not_identify_PPG,1),1);
       goodBins.ppg_template_manually_excluded = zeros(size(could_not_identify_PPG,1),1);
    end
end