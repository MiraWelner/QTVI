function [analysisFiles] = windowFeaturesSetup(inputpath)

    file_list = dir(fullfile(inputpath, '**/*feature_output.mat'));

    analysisFiles = cell(length(file_list), 5);   
    for i = 1:length(file_list)
        analysisFiles{i, 2} = fullfile(file_list(i).folder, file_list(i).name);

        [~, name, ~] = fileparts(analysisFiles{i, 2});
        analysisFiles{i, 1} = name(1:end-15);
    end

    
end