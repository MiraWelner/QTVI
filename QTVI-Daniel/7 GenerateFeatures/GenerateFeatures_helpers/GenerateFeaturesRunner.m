% pulse annotate script
clear;
% close all;

try
    parpool;
catch
end

props = readProps('config.txt');
GenerateFeatures_template_type = str2num(props('Gen_Features_use_manually_reviewed_templates'));
GenerateFeatures_annealed_path = props('GenerateFeatures_annealed_path');
GenerateFeatures_wavedata_path = props('GenerateFeatures_wavedata_path');
GenerateFeatures_template_path = props('GenerateFeatures_template_path');
GenerateFeatures_output = props('GenerateFeatures_output');
Skip_Existing = logical(str2num(props('Skip_Existing')));

[analysisFiles] = GenerateFeaturesSetup(GenerateFeatures_annealed_path, GenerateFeatures_wavedata_path, GenerateFeatures_template_path, GenerateFeatures_template_type);

disp('*********************************************************************');
time = 0;
success = 0;
fail = 0;

for i = 1:size(analysisFiles, 1)
    if isfile(fullfile(GenerateFeatures_output, [analysisFiles{i,1} '_feature_output.mat'])) && Skip_Existing
        disp([analysisFiles{i,1} '_feature_output.mat exists skipping because Skip_Existing = 1 in config.']);
       continue; 
    end

    tStart = tic;
    disp(join(['Beginning analysis of ' analysisFiles{i,1} ' | ' num2str(i) ' of ' num2str(size(analysisFiles,1))]));
    avg_time = time/i;
    disp(['Avg Time (s): ' num2str(avg_time)]);

    disp(['Est finish (min): ' num2str((avg_time*(size(analysisFiles,1)-i))/60)]);
    disp(join(['Output loc ' GenerateFeatures_output]));

    disp('Loading Data...')
    r = GenerateFeatures(analysisFiles{i, 2}, analysisFiles{i, 3}, analysisFiles{i, 4}, GenerateFeatures_output);
    if r
        success = success +1;
    else
        fail = fail +1;
    end
    
    disp(['____________________________________________________________________________________________________' newline]);
    time = time + toc(tStart);
    toc(tStart);
end
% %     disp(join(['Success: ' num2str(success/size(analysisFiles, 1)*100) ]));
% %         disp(join(['Fail: ' num2str(fail/size(analysisFiles, 1)*100) ]));