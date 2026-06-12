clear
props = readProps('config.txt');
GenerateTemplates_input = props('GenerateTemplates_input');
GenerateTemplates_output = props('GenerateTemplates_output');
Skip_Existing = logical(str2num(props('Skip_Existing')));

list = dir(fullfile(GenerateTemplates_input, '**/*wave_data.mat'));
time = 0;
parfor i = 1:size(list, 1)
    name = list(i).name;
    start_idx = regexp(name, '_wave_data');
    start_idx = start_idx(1);
    name = name(1:start_idx-1);
    if isfile(fullfile(GenerateTemplates_output, [name '_template_info.mat']))
        disp([name '_template_info.mat exists skipping because Skip_Existing = 1 in config.']);
        continue
    end
    file = fullfile(list(i).folder, list(i).name);

    tStart = tic;

    disp(join(['Beginning analysis of ' file ' | ' num2str(i) ' of ' num2str(length(list))]));
%     avg_time = time/i;
%     disp(['Avg Time (s): ' num2str(avg_time)]);

%     disp(['Est finish (min): ' num2str((avg_time*(length(list)-i))/60)]);
    disp(join(['Output loc ' GenerateTemplates_output]));
        
    r = GenerateTemplates(file, GenerateTemplates_output);
%     if r == 1 
%         movefile(file, '/hdd/data/mesa/Manuscript 1/2020 November run/3 wave data');
%     end
    
%     disp(['____________________________________________________________________________________________________' newline]);
%     time = time + toc(tStart);
end
% disp(time)