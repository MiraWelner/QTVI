function [analysisFiles, order] = SaveSetup(metadata_path, final_output_path, outputLoc)
    matched = import_matched_ids(metadata_path);
    file_list = dir(fullfile(final_output_path, '**/*feature_output.mat'));
    order = matched.Properties.VariableNames;
    
    analysisFiles = cell(length(file_list), size(matched,2));   
    for i = 1:length(file_list)
        analysisFiles{i, 2} = fullfile(file_list(i).folder, file_list(i).name);

        [path, name, ~] = fileparts(analysisFiles{i, 2});
        analysisFiles{i, 1} = name(1:end-17);

        s = split(analysisFiles{i, 1},'_');
        index = matched.id == str2num(s{1});
        index = find(index == 1);
        
        for x = 1:size(matched,2)
            analysisFiles{i, 2+x} = matched.(order{x})(index);
        end
    end

    
end