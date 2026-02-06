 clear
props = readProps('config.txt');
Annealed_segments = props('FWB_input_path');
output_path = props('FWB_output_path');
Skip_Existing = logical(str2num(props('Skip_Existing')));

list = dir(fullfile(Annealed_segments, '**/*annealedSegments.mat'));
list = flip(list);
time = 0;
for i = 1:size(list, 1)
    % set name   
    name = list(i).name;
    start_idx = regexp(name, '_annealedSegments');
    start_idx = start_idx(1);
    name = name(1:start_idx-1);
    if Skip_Existing && isfile(fullfile(output_path, [name '_wave_data.mat']))
        disp([name '_wave_data.mat exists skipping because Skip_Existing = 1 in config.']);
        continue
    end
    
    file = fullfile(list(i).folder, list(i).name);

    tStart = tic;

    disp(join(['Beginning analysis of ' file ' | ' num2str(i) ' of ' num2str(length(list))]));
    avg_time = time/i;
    disp(['Avg Time (s): ' num2str(avg_time)]);

    disp(['Est finish (min): ' num2str((avg_time*(length(list)-i))/60)]);
    disp(join(['Output loc ' output_path]));
    r = FindWaveBounds(file, output_path);
%     if r == 1 
%         movefile(file, '/hdd/data/mesa/Manuscript 1/2020 November run/2 Anneal Segs/');
%     end
    disp(['____________________________________________________________________________________________________' newline]);
    time = time + toc(tStart);

end
