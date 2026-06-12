function [analysisFiles] = SaveSetup( final_output_path, outputLoc)
%     matched = import_matched_ids(metadata_path);
    file_list = dir(fullfile(final_output_path, '**/*feature_output.mat'));

    analysisFiles = cell(length(file_list), 5);   
    for i = 1:length(file_list)
        analysisFiles{i, 2} = fullfile(file_list(i).folder, file_list(i).name);

        [path, name, ~] = fileparts(analysisFiles{i, 2});
        analysisFiles{i, 1} = name(1:end-15);

        analysisFiles{i, 5} = outputLoc;
%         s = split(analysisFiles{i, 1},'_');
%         index = matched.idno == str2num(s{1});
%         index = find(index == 1);
%         
%         analysisFiles{i, 3} = matched.idno(index);
%         analysisFiles{i, 4} = matched.mesaid(index);
        
%         analysisFiles{i, 5} = matched.Healthy(index);
%         analysisFiles{i, 6} = matched.MIDM(index);
%         analysisFiles{i, 7} = matched.MIDM1(index);
%         analysisFiles{i, 8} = matched.DMRX(index);
%         analysisFiles{i, 9} = matched.DMRX1(index);

    end

    
end