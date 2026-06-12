function [ret] = GenerateTemplates(wave_data_path, output_path)
%     try
        std_multiplier = 2.5;
        threshold_percent = 15;

        [~, name, ~] = fileparts(wave_data_path);
        name = regexp(name, '(\d+_\d+)', 'match');
        name = name{1};

        wave_data = load(wave_data_path);
        wave_data = wave_data.wave_data;
        
        updated = 0;
        for i = 1:length(wave_data)
            %[alignSegs] = FindAlignmentSegment(bins{i}.ppgSeg, bins{i}.ppgMinAmps,bins{i}.ppgMaxAmps);
            if isfield(wave_data,'bad_segment') == 0
                wave_data{i}.bad_segment = 0;
                updated = 1;
            end
        end
        if updated == 1
            save(wave_data_path, 'wave_data');
        end

        [template_matrix] = CreatePPGTemplates(wave_data, std_multiplier, 0);
        
        template_good_mask = ones(size(template_matrix,1),1);
        
        for l = 1:size(template_matrix,1)
           if isempty(template_matrix(l,~isnan(template_matrix(l,:))))
              template_good_mask(l) = 0; 
           end
        end

%         template_good_mask = logical(template_good_mask);
        %[garys_template_matrix,type] = CreateGaryTemplates(wave_data, 60*256, 256);

        disp('Finding template feet...');
%         [~, idx] = find_foot_pulseox(template_matrix(template_good_mask,:), 0);
        [~, idx] = find_foot_pulseox(template_matrix, 0);

        %[~, gidx] = find_foot_pulseox(garys_template_matrix, 0);

        disp('Aligning templates...');
%         alignedTemplates = AlignWaves(template_matrix(template_good_mask,:), idx);
        alignedTemplates = AlignWaves(template_matrix, idx);

        %galignedTemplates = AlignWaves(garys_template_matrix, gidx);

%         disp('Combining similar templates...');
        template_Diffmatrix = WaveDiff(alignedTemplates);
    %     CombineTemplates(alignedTemplates);
    %     

        [bin_template_numbers, ppgTemplates, ~] = CombineTemplatesGraph(alignedTemplates, template_Diffmatrix, threshold_percent, 0);
        
%         for x =1:length(ppgTemplates)
%             hold on;
%             plot(ppgTemplates{x})
%         end
        
        % bins are all seperate atm this works but is commented out as each
        % has its own template
        %disp('Combining bins...');
        %[bins] = CombineBins(bin_template_numbers, wave_data);
        
        % not used atm wanted to use it to furhter rull
        [ecgTemplates, alignment_points, avg_r_expand] = CreateEcgTemplates(wave_data, std_multiplier, 0);

        template_info = cell(length(ecgTemplates), 1);

        for t = 1:length(ecgTemplates)
            if template_good_mask(t)
                info.index = t;
                info.ppg_bin_indexs = wave_data{t}.ppg_bin_indexs;
                info.ecg_bin_indexs = wave_data{t}.ecg_bin_indexs;
                info.ecgSamplingRate = wave_data{t}.ecgSamplingRate;
                info.ppgSamplingRate = wave_data{t}.ppgSamplingRate;
                info.bad_segment = wave_data{t}.bad_segment;

                
                info.ecgTemplate = ecgTemplates{t};
%                 info.ecgTemplate = [];
                info.ppgTemplate = ppgTemplates{bin_template_numbers(t)};
                info.alignment_point = alignment_points(t);
                info.avg_r_expand = avg_r_expand(t);
                template_info{t} = info;
            else
                info.index = t;
                info.ppg_bin_indexs = wave_data{t}.ppg_bin_indexs;
                info.ecg_bin_indexs = wave_data{t}.ecg_bin_indexs;
                info.ecgSamplingRate = wave_data{t}.ecgSamplingRate;
                info.ppgSamplingRate = wave_data{t}.ppgSamplingRate;
                info.bad_segment = wave_data{t}.bad_segment;
                
                info.ecgTemplate = [];
                info.ppgTemplate = [];
                info.alignment_point = nan;
                info.avg_r_expand = avg_r_expand(t);
                template_info{t} = info;
            end

        end

        disp('Saving...');
        save(fullfile(output_path, [name '_template_info']), 'template_info');
        ret = 1;
%     catch e
%         input.anneal_path = wave_data_path;
%         st = dbstack;
%         namestr = st.name;
%         LogError(namestr, output_path, input, e);
%         ret = 0;
%     end
end
