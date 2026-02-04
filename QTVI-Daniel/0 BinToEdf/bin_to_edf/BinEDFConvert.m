% Script by: Noah Campagna
% Date: 11/27/2020
% Bin/EDF Converter


close all;
clear;


disp('Bin/EDF Converter');
disp('___________________________________________________________________________________________');


% Choose to conversion direction

convertType = 1; % fixed for this version of code

props = readProps('config.txt');

% 
% flag = true;
% while(flag)
%     disp("Type 1 to merge a .bin and .inf file and convert it to a .edf file format.")
%     in = input("Type 2 to convert a .edf file into a .bin and a .inf file format.  \n: ");
%     clc
% 
%     if(in == 1)
%         convertType = 1;
%         flag = false;
%     elseif(in == 2)
%         convertType = 2;
%         flag = false;
%     else
%         fprintf("\n\n\n\n\n\n\n\nInvalid response. Please enter 1 or 2.");
%         fprintf('');
%     end
% end


%%
% Convert .bin and .inf to .edf
try
   
if(convertType == 1)
    
    % Choose files
%     input("Press enter to choose the data directory (.bin and .inf files)");
%     [baseName, folder] = uigetfile('..\*.bin');
%     binFolder = uigetdir();
    binFolder = props('BinToEdf_input_path');
    Skip_Existing = logical(str2num(props('Skip_Existing')));
%     clc
    
%     binFiles = dir(binFolder);
    binFiles = dir(fullfile(binFolder, '*.bin'));
    infFiles = dir(fullfile(binFolder, '*.inf'));

    binFileLoc = string([]);
    for x = 1:length(binFiles)
        binFileLoc(end+1) = fullfile(binFolder,binFiles(x).name);
    end
    
    fileOutputName = string([]);
    for x = 1:length(binFiles)
        fileOutputName(end+1) = strcat(binFiles(x).name(1:length(binFiles(x).name)-3), 'edf');
    end
    
%     fileOutputName = strcat(baseName(1:length(baseName)-3), 'edf');
%     binFileLoc = fullfile(folder, baseName);

%     input("Press enter to choose the .inf files directory.");
% %     [baseName, folder] = uigetfile('../*.inf');
%     infFolder = uigetdir();
%     clc
%     
    headFileLoc = string([]);
    for x = 1:length(infFiles)
        headFileLoc(end+1) = fullfile(binFolder, infFiles(x).name);
    end
    
%     input("Press enter to choose the ouput location.");
%     outputLoc = uigetdir;
%     clc
    outputLoc = props('BinToEdf_output_path'); 

    for b = 1:length(binFileLoc)
        
        infFileLoc = strrep(binFileLoc(b), '.bin', '.inf');
        % Open header, pick output location
        
        try
            fileID = fopen(infFileLoc);
        catch
            disp('.inf file not found.');
        end
            
           
        
        edfFnOut = char(fullfile(outputLoc, fileOutputName(b)));
        if(isfile(edfFnOut) && Skip_Existing)
            fprintf('%s already exists in this directory. Skipping because Skip_Existing = 1 in config.\n', fileOutputName);
            continue
        end
        fprintf('Converting file...');

        % Read main header file info
        head = textscan(fileID, '%s', 'delimiter', '\n');
        head = head{1};
        keys = {'Patient', 'Number of Channel', 'Points for Each Channel', 'Data Sampling Rate', 'Start Time', 'Channel Number'};
        for x = 1:length(keys)
            for y = 1:length(head)
                if(contains(head{y}, keys{x}))

                    switch x
                        case 1
                            line = strrep(head{y}, 'Patient = ', ''); %patient ID
                            patient_id = strrep(line, 'Patient = ', '');
                            break;

                        case 2
                            num_signals = strrep(head{y}, 'Number of Channel = ', ''); %num of signals
                            num_signals = str2double(num_signals);
                            break;

                        case 3
    %                         num_data_records = strrep(head{y}, 'Points for Each Channel = ', ''); %num data records
    %                         num_data_records = round(str2double(num_data_records) / 977);
                            break;

                        case 4
                            line = strrep(head{y}, 'Data Sampling Rate = ', '');
                            sampling_rate = '';
                            for z = 1:length(line)
                                if(line(z) ~= ' ')
                                    sampling_rate(end+1) = line(z);
                                else
                                    break
                                end
                            end
                            sampling_rate = str2double(sampling_rate);
                            break;

                        case 5
                            line = strrep(head{y}, 'Start Time = ', ''); 
                            try
                                recording_startdate = datetime(line, 'InputFormat', 'MM/dd/yyyy h:mm:ss a');
                            catch
                                line = strrep(line, 'PM', ''); 
                                line = strrep(line, 'AM', '');
                                recording_startdate = datetime(line);
                            end

                            recording_startdate.Format = 'dd-MM-yy'; %start date
                            recording_starttime = recording_startdate;
                            recording_startdate = char(recording_startdate);

                            recording_starttime.Format = 'HH:mm:ss'; %start time
                            recording_starttime = char(recording_starttime);
                        otherwise
                            channel = y + 1;  
                    end
                end
            end
        end

        % Read channel info from header file
        chan_names = {};
        chan_nums = {};
        for x = channel:length(head)
                name = [];
                number = head{x}(1);
                numFlag = true;
                nameFlag = true;
                space_loc = 0;
                charCount = 0;
            for y = 2:length(head{x})

                if((head{x}(y) ~= ' ') && (numFlag))
                    number = strcat(number, head{x}(y));
                    charCount = charCount + 1;
                elseif(head{x}(y) == ' ' && (nameFlag))
                    numFlag = false;
                    charCount = charCount + 1;
                else
                    name = strcat(name, head{x}(y));
                    if(head{x}(y) == ' ')
                        space_loc = y;
                    end
                    nameFlag = false;
                end
            end

            if(space_loc > 0)
                name1 = name(1:(space_loc - charCount-2));
                name2 = name((space_loc-charCount-1):end);
                space = ' ';
                name = [name1, space, name2];

                chan_names{end+1} = name;
            else

                if(strcmp(name, 'I'))
                    name = 'ECG I';
                elseif(strcmp(name, 'II') || strcmp(name, 'ECG-II') || strcmp(name, 'ECG II'))
                    name = 'EKG';
                elseif(strcmp(name, 'III'))
                    name = 'ECG III';
                end

                chan_names{end+1} = name;
            end
            chan_nums{end+1} = number;

        end

        % Hard-coded extra header info
        reserve_1 = ''; %reserve 1
        num_header_bytes = 256+256*num_signals;
        data_record_duration = 1; %data record duration, temporary
        header.edf_ver = 1;
        fclose(fileID);

        % Read binary file
        fileID = fopen(binFileLoc(b));
        binData = fread(fileID, [num_signals inf], 'double');
        fclose(fileID);

        % Cut off extra data for conversion
        % Makes total number of records per channel divisble by sampling rate
        remain = mod(length(binData),sampling_rate);
        binData = binData(1:num_signals, 1:(length(binData) - remain));
        num_data_records = length(binData)/sampling_rate;

        % Populate Header info
        header.edf_ver = '0';
        header.patient_id = patient_id;
        header.local_rec_id = '';
        header.recording_startdate = recording_startdate;
        header.recording_starttime = char(recording_starttime);
        header.num_header_bytes = num_header_bytes;
        header.reserve_1 = reserve_1;
        header.num_data_records = num_data_records;
        header.data_record_duration = data_record_duration;
        header.num_signals = num_signals;

        % Populate Signal Header info
        for s = 1 : num_signals
            signalHeader(s).signal_labels = chan_names{s};
            signalHeader(s).tranducer_type = '';
            if(strcmp(chan_names{s}, 'PPG'))
                signalHeader(s).physical_dimension = '';
            elseif(strcmp(chan_names{s}, 'SpO2'))
                signalHeader(s).physical_dimension = '%';
            elseif(strcmp(chan_names{s}, 'ABP') || strcmp(chan_names{s}, 'CVP') || strcmp(chan_names{s}, 'ART'))
                signalHeader(s).physical_dimension = 'mmHg';
            else
                signalHeader(s).physical_dimension = 'mV';
            end
            signalHeader(s).physical_min = max(max(binData));
            signalHeader(s).physical_max = min(min(-5));
            signalHeader(s).digital_min = -32768;
            signalHeader(s).digital_max = 32767;
            signalHeader(s).prefiltering = '';
            signalHeader(s).samples_in_record = sampling_rate;
            signalHeader(s).reserve_2 = ' ';
            signalCell{s} = binData(s,:);
        end

        % Write edf file
        status = blockEdfWrite(edfFnOut, header, signalHeader, signalCell);
        [header signalHeader signalCell] = blockEdfLoad(edfFnOut);
        fprintf('%s has been converted to .edf format.\n', binFiles(b).name);
    end
end
    
%%
% Covert .edf to .bin and .inf

if(convertType == 2)
    
    % Choose files
    input("Press enter to choose the .edf file location.");
    [baseName, folder] = uigetfile('../*.edf');
    clc
    edfFileLoc = fullfile(folder, baseName);
    
    [header signalHeader signalCell] = blockEdfLoad(edfFileLoc);
    
    input("Press enter to choose the ouput location.");
    outputLoc = uigetdir;
%     outputLoc = strcat(outputLoc, '\');
    infPath = fullfile(outputLoc,[header.patient_id '.inf']);
    binPath = fullfile(outputLoc,[header.patient_id '.bin']);
    clc
    
    % Create header file
    fileID = fopen(infPath,'wt');
    fprintf(fileID, 'Patient = ');
    fprintf(fileID, header.patient_id);
    fprintf(fileID, '\n');
    
    fprintf(fileID, 'Description = \n');
    
    fprintf(fileID, 'Export Date = \n');
    
    fprintf(fileID, 'Number of Channel = ');
    fprintf(fileID, num2str(header.num_signals));
    fprintf(fileID, '\n');
    
    fprintf(fileID, 'Points for Each Channel = ');
    fprintf(fileID, num2str(length(signalCell{1})));
    fprintf(fileID, '\n');
    
    fprintf(fileID, 'Data Sampling Rate = ');
    fprintf(fileID, num2str(length(signalCell{1})/header.num_data_records));
    fprintf(fileID, ' points/second\n');
    
    fprintf(fileID, 'Start Time = ');
    fprintf(fileID, header.recording_startdate);
    fprintf(fileID, ' ');
    fprintf(fileID, header.recording_starttime);
    fprintf(fileID, '\n');
    
    fprintf(fileID, 'Stop Time = \n');
    
    fprintf(fileID, 'Units: \n');
    
    fprintf(fileID, 'Channel Number  Channel Label\n');
    

    for x = 1 : length(signalHeader)
        fprintf(fileID, num2str(x));
        fprintf(fileID, '              ');
        fprintf(fileID, signalHeader(x).signal_labels);
        fprintf(fileID, '\n');
    end
    
    fclose(fileID);
    
    % Convert data to binary
    fileID = fopen(binPath,'w');
    data = cell2mat(signalCell);
    reformData = zeros(1, header.num_signals * length(data)); 
    count = 1;
    for y = 1 : length(data)
        for x = 1 : header.num_signals
            reformData(count) = data(y, x);
            count = count + 1;
        end
    end
    
    fwrite(fileID, reformData, 'double');
    fclose(fileID);
    fprintf('Success! Your file has been converted to a .bin and .inf file.');
end
catch
end