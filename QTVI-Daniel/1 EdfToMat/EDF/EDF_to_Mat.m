props = readProps('../config.txt');
EdfToMat_input_path = props('EdfToMat_input_path'); %'/hdd/data/mesa/EDF All';
EdfToMat_output_path = props('EdfToMat_output_path'); %'/hdd/data/mesa/MATS';
Skip_Existing = logical(str2num(props('Skip_Existing')));

edf_list = dir(fullfile(EdfToMat_input_path, '**/*.edf'));
xml_list = dir(fullfile(EdfToMat_input_path, '**/*edf.XML'));

if isempty(xml_list)
    sleep_data_exists = 0;
else
    sleep_data_exists = 1;
end

for i = 1:length(edf_list)
    tic;
    disp([num2str(i) ' of ' num2str(length(edf_list))]);
    edf = fullfile(edf_list(i).folder, edf_list(i).name);
    if Skip_Existing && isfile(fullfile(EdfToMat_output_path, [edf_list(i).name(1:end-4) '_ecg_ppg_sleep.mat']))
        disp([edf_list(i).name(1:end-4) '_ecg_ppg_sleep.mat exists skipping because Skip_Existing = 1 in config.']);
        continue
    end
    
    if sleep_data_exists
        xml = fullfile(xml_list(i).folder, xml_list(i).name);
    else
        xml = '';
    end
    convert(edf, xml, EdfToMat_output_path);
    toc
end



function convert(edf_path,xml_path, output_path)
    [~,name,~] = fileparts(edf_path);
    

    % don't remake if it exists
    if isfile(fullfile(output_path, [name '_ecg_ppg_sleep.mat']))
        return
    end

    if strcmp(xml_path,'')
       [edf_hdr, edf_data] = edfread(edf_path, 'verbose', 0, 'targetSignals', {'EKG','PPG'});
    else
       [edf_hdr, edf_data] = edfread(edf_path, 'verbose', 0, 'targetSignals', {'EKG','Pleth'}); % 1 is ecg index, 23 pulse
    end
    ppg = edf_data(2, :);
    ecg = edf_data(1, :);

    ecgSamplingRate = edf_hdr.frequency(1); % 256 FOR MESA SET
    ppgSamplingRate = edf_hdr.frequency(2); % 256 FOR MESA SET
    try
        if strcmp(xml_path,'')
            scoring_epoch_size_sec = 30.0;
            sleepStages = length(edf_data) / scoring_epoch_size_sec ;
            sleepStages = zeros(fix(sleepStages / ecgSamplingRate), 1);
        else
            [scoring_epoch_size_sec, sleepStages, ~] = ReadXML(xml_path);
        end
        
    catch
        return
    end
    [sleepStages] = RenumberSleepStages(sleepStages);
    save(fullfile(output_path, [name '_ecg_ppg_sleep']), 'ppg','ecg','ecgSamplingRate','ppgSamplingRate','sleepStages','scoring_epoch_size_sec');
%     catch e
%         input.edf_path = edf_path;
%         input.output_path = output_path;
%         st = dbstack;
%         namestr = st.name;
%         LogError(namestr,output_path, input, e);
% 
%     end
end