
function [analysisFiles] = SetupAnneal(data_path, noise_path, noise_flag)
    if noise_flag == 1
        data_list = dir(fullfile(data_path, '**/*_ecg_ppg_sleep.mat'));
        noise_list = dir(fullfile(noise_path, '**/*_noise*.mat'));

        tmp = cell(length(data_list), 2);
        patient_ids = cell(length(data_list), 1);
        for i = 1:length(data_list)
            tmp{i, 2} = fullfile(data_list(i).folder, data_list(i).name);
            [~, name, ~] = fileparts(tmp{i, 2});
            start_idx = regexp(name, '_ecg_ppg_sleep');
            start_idx = start_idx(1);
            tmp{i, 1} = name(1:start_idx-1);
            patient_ids{i, 1} = name(1:start_idx-1);
        end

        analysisFiles = cell(length(noise_list), 3);
        for i = 1:length(noise_list)
                noise = fullfile(noise_list(i).folder, noise_list(i).name);
                [~, name, ~] = fileparts(noise);

                start_idx = regexp(name, '_noise');
                start_idx = start_idx(1);
                
                % set name
                analysisFiles{i, 1} = name(1:start_idx-1);

                Index = count(patient_ids,  name(1:start_idx-1)); %analysisFiles{i, 1});
                Index = find(Index == 1);
                if numel(Index) ~= 0
                    if length(Index) > 1
                        Index = Index(1);
                    end
                    analysisFiles{i, 2} = noise;
                    analysisFiles{i, 3} = tmp{Index,2};
                end
        end
    elseif noise_flag == 0
        noise_list = {};
        data_list = dir(fullfile(data_path, '**/*_ecg_ppg_sleep.mat'));
        if noise_flag~=0
            rs_list = regexpdir(noise_path, '\d_golden_rs');
            n = regexpdir(noise_path, '\dnoiseclass_golden_rs');

            for x = 1:length(n)
                [~, name, ~] = fileparts(n{x});
                name = regexp(name, '(\d+_\d+)', 'match');
                name = name{1};
                noise_list{end+1} = name;
            end
        end
        tmp = cell(length(data_list), 2);
        patient_ids = cell(length(data_list), 1);
        for i = 1:length(data_list)
            tmp{i, 2} = fullfile(data_list(i).folder, data_list(i).name);
            [~, name, ~] = fileparts(tmp{i, 2});
            name = regexp(name, '(\d+_\d+)', 'match');
            tmp{i, 1} = name;
            patient_ids{i, 1} = name{1};
        end
        
        analysisFiles = cell(length(data_list), 4);
        for i = 1:length(data_list)
            data = data_list(i);
            name = data.name;
            name = regexp(name, '(\d+_\d+)', 'match');
            name = name{1};
            
            Index = count(noise_list,  name);
            Index = find(Index == 1);
            
            if ~isempty(Index)
                analysisFiles{i, 1} = name;

                analysisFiles{i, 2} = n{Index};
                analysisFiles{i, 3} = fullfile(data.folder,data.name);

                rs = rs_list(Index);
                rs = rs{1};
                analysisFiles{i, 4} = rs;
            else
                analysisFiles{i, 1} = name;

                analysisFiles{i, 3} = fullfile(data.folder,data.name);

            end
        end
%         
%         for i = 1:length(noise_list)
%             noise = noise_list(i);
%             noise = noise{1};
%             [~, name, ~] = fileparts(noise);
%             name = regexp(name, '(\d+_\d+)', 'match');
%             name = name{1};
%             
%             Index = count(patient_ids,  name);
%             Index = find(Index == 1);
%             
%             if ~isempty(Index)
%                 analysisFiles{i, 1} = name;
% 
%                 analysisFiles{i, 2} = noise;
%                 analysisFiles{i, 3} = tmp{Index,2};
% 
%                 rs = rs_list(i);
%                 rs = rs{1};
%                 analysisFiles{i, 4} = rs;
%             else
%                 disp(['No matchs for: ' name]);
%             end
%         end
%         
%         i=1;
%         while true
%             if i > size(analysisFiles,1)
%                 break
%             end
%             if isempty(analysisFiles{i,1})
%                 analysisFiles(i,:)= [];
%                 continue
%             end
%             i = i +1;
%         end
    else
        data_list = dir(fullfile(data_path, '**/*_ecg_ppg_sleep.mat'));
        
        tmp = cell(length(data_list), 2);
        patient_ids = cell(length(data_list), 1);
        for i = 1:length(data_list)
            tmp{i, 2} = fullfile(data_list(i).folder, data_list(i).name);
            [~, name, ~] = fileparts(tmp{i, 2});
            name = regexp(name, '(\d+_\d+)', 'match');
            tmp{i, 1} = name;
            patient_ids{i, 1} = name{1};
        end
        
        analysisFiles = cell(length(data_list), 3);
        for i = 1:length(data_list)
                noise = fullfile(data_list(i).folder, data_list(i).name);
                [~, name, ~] = fileparts(noise);
                name = regexp(name, '(\d+_\d+)', 'match');

                analysisFiles{i, 1} = name{1};
                
                Index = count(patient_ids,  name); %analysisFiles{i, 1});
                Index = find(Index == 1);
                if numel(Index) ~= 0
                    analysisFiles{i, 3} = tmp{Index,2};
                end
        end
        
        
    end

end