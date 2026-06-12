props = readProps('path.txt');
alledfs_path = props('EDF_to_Mat_input_path'); %'/hdd/data/mesa/EDF All';
output_path = props('EDF_to_Mat_output_path'); %'/hdd/data/mesa/MATS';

edf_list = dir(fullfile(alledfs_path, '**/*.edf'));
xml_list = dir(fullfile(alledfs_path, '**/*edf.XML'));

for i = 1:length(edf_list)
    tic;
    disp([num2str(i) ' of ' num2str(length(edf_list))]);
    edf = fullfile(edf_list(i).folder, edf_list(i).name);
    xml = fullfile(xml_list(i).folder, xml_list(i).name);
    convert(edf, xml, output_path);
    toc
end



function convert(edf_path,xml_path, output_path)
    try 
        [~,name,~] = fileparts(edf_path);
        [edf_hdr, edf_data] = edfread(edf_path, 'verbose', 0, 'targetSignals', [1, 23]); % 1 is ecg index, 23 pulse
        ppg = edf_data(2, :);
        ecg = edf_data(1, :);

        ecgSamplingRate = edf_hdr.frequency(1); % 256 FOR MESA SET
        ppgSamplingRate = edf_hdr.frequency(23); % 256 FOR MESA SET
        [scoring_epoch_size_sec, sleepStages, ~] = ReadXML(xml_path);
        [sleepStages] = RenumberSleepStages(sleepStages);
        save(fullfile([output_path '/' name '_ecg_ppg_sleep']), 'ppg','ecg','ecgSamplingRate','ppgSamplingRate','sleepStages','scoring_epoch_size_sec');
    catch e
        input.edf_path = edf_path;
        input.output_path = output_path;
        st = dbstack;
        namestr = st.name;
        LogError(namestr,output_path, input, e);

    end
end