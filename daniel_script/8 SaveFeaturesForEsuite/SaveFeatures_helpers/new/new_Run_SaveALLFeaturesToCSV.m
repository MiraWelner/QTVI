clear;
close all;

props = readProps('config.txt');
metadata_path = props('combinedFeatures_metadata');
features_path = props('combinedFeatures_input');
outputLoc = props('combinedFeatures_output');
output_file = fullfile(outputLoc, ['ppg_data_' datestr(now,'mm-dd-yyyy HH-MM') '.csv']);

FID = -1;
headers = [];

first = 0;
[analysisFiles,columnNames] = SaveSetup(metadata_path,features_path, outputLoc);
time = 0;
for i = 1:size(analysisFiles,1)
    tStart = tic;

    avg_time = time/i;
    disp(['Avg Time (s): ' num2str(avg_time)]);

    disp(['Est finish (min): ' num2str((avg_time*(size(analysisFiles,1)-i))/60)]);
    disp(join(['Output loc ' output_file]));
    disp(join(['Saving analysis of ' analysisFiles{i, 1}]));
    disp(['Remaining: ' num2str(size(analysisFiles,1)-i)]);
    
    windowed_feature = load(analysisFiles{i, 2});
    windowed_feature = windowed_feature.windowed_feature;
        
    if first == 0
        first = 1;
        FID = fopen(output_file,'w');
        headers = fieldnames(windowed_feature);
        m = length(headers);
        sz = zeros(m,2);
        
        l = '';
        for ii = 1:m
            sz(ii,:) = size(windowed_feature.(headers{ii}));   
            if ischar(windowed_feature.(headers{ii}))
                sz(ii,2) = 1;
            end
            l = [l,'',headers{ii},',',repmat(',',1,sz(ii,2)-1)];
        end
        l = [l(1:end-1),'\n'];

        fprintf(FID,l);
        output = toOutput(windowed_feature);
    else
        output = toOutput(windowed_feature);
    end
    

    dlmwrite(output_file,output,'-append','precision',15);
    disp(['____________________________________________________________________________________________________' newline]);
    time = time+toc(tStart);
    toc(tStart);
end

function output = toOutput(struct)
    headers = fieldnames(struct);
    output = nan(length(struct.area_n), length(headers));
    l = length(headers);
    for h = 1:l
        output(:,h) = struct.(headers{h});
    end
end
