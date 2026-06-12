clear;
input_dir = '/home/deeplab/Desktop/aha/30_min_bins_with_noise_data';
sample_rate = 256;

myFiles = dir(fullfile(input_dir,'*output.mat')); 
total = 0;
for k = 1:length(myFiles)
    disp('-----------------------------');
    baseFileName = myFiles(k).name;
    fullFileName = fullfile(input_dir, baseFileName);
    [~, name, ~] = fileparts(fullFileName);
    name = name(1:end-9);
    fprintf(1, 'Reading %s\n', fullFileName);
    out = load(fullFileName);
    out = out.output;
    
    % for 30 min seg in ppg
    for i = 1:length(out)
        total = total+1;
        %plot(out{i,1}.po)
        %plot(out{i,1}.ecg)
        disp([name ' 30 min bin # ' num2str(i) ' statistics:']);
        disp(['Time Excluded (mins) | ' num2str(out{i,1}.timeExcluded)]);
        disp(['Time Kept (mins) | ' num2str(out{i,1}.timeGood)]);
        disp(['% Good | ' num2str(out{i,1}.precentGood*100)]);
        disp('________________________________________________');
        
        
        %% Call what you want here %%
        
        
        
    end    
        
    
end

disp(['Total bins: ' num2str(total)])