function [analysisFiles] = TemplateSetup(guiless_path, outputLoc)
    template_info_list = dir(fullfile(guiless_path, '**/*_template_info.mat'));

    analysisFiles = cell(length(template_info_list), 4);   
    
    for i = 1:length(template_info_list)
        [~, name, ~] = fileparts(template_info_list(i).name);
        start_idx = regexp(name, '_template_info');
        start_idx = start_idx(1);
        name = name(1:start_idx-1);
        analysisFiles{i, 2} = fullfile(template_info_list(i).folder, template_info_list(i).name);
        analysisFiles{i, 1} = name;
        analysisFiles{i, 4} = outputLoc;
    end

    
end