% mat_to_edf_converter.m
clear; clc;

matFolder = uigetdir(pwd, 'Select folder containing SHHS .mat files');
if matFolder == 0, error('No folder selected.'); end
matFiles = dir(fullfile(matFolder, '*.mat'));

for i = 1:length(matFiles)
    matFilename = fullfile(matFolder, matFiles(i).name);
    [~, baseName, ~] = fileparts(matFiles(i).name);
    edfFilename = fullfile(matFolder, [baseName '.edf']);
    
    try
        data = load(matFilename);
        write_edf_robust(edfFilename, data);
        fprintf('Converted: %s\n', [baseName '.edf']);
    catch ME
        fprintf('Error in %s: %s\n', matFiles(i).name, ME.message);
    end
end

function write_edf_robust(filename, data)
    % 1. Extract Signals
    if isfield(data, 'signal'), signals = data.signal;
    elseif isfield(data, 'val'), signals = data.val;
    else, fields = fieldnames(data); signals = data.(fields{1}); end
    
    if ~iscell(signals)
        if isstruct(signals)
            fn = fieldnames(signals);
            for j=1:length(fn), s{j} = signals.(fn{j}); end
            signals = s;
        else, signals = {signals}; end
    end
    
    ns = length(signals);
    fs = 256 * ones(1, ns); % Default sampling rate
    recordDuration = 1;
    numDataRecords = ceil(length(signals{1}) / fs(1));
    
    % 2. Open file in BINARY mode (CRITICAL)
    fid = fopen(filename, 'wb', 'ieee-le');
    if fid < 0, error('File open failed'); end
    
    % Helper to write fixed-length strings
    write_fix = @(str, len) fwrite(fid, sprintf(['%-', num2str(len), 's'], str(1:min(end,len))), 'char');

    % --- Main Header (256 bytes) ---
    write_fix('0', 8);           % Version
    write_fix('Unknown', 80);    % Patient ID
    write_fix('SHHS Data', 80);  % Recording ID
    write_fix(datestr(now, 'dd.mm.yy'), 8);
    write_fix(datestr(now, 'HH.MM.SS'), 8);
    write_fix(num2str(256 + ns * 256), 8); % Header bytes
    write_fix('', 44);           % Reserved
    write_fix(num2str(numDataRecords), 8);
    write_fix(num2str(recordDuration), 8);
    write_fix(num2str(ns), 4);
    
    % --- Signal Headers (ns * 256 bytes) ---
    % Label: Translated "Signal1" logic to "ECG"
    for j = 1:ns
        write_fix(sprintf('ECG%d', j), 16); 
    end
    for j = 1:ns, write_fix('', 80); end % Transducer
    for j = 1:ns, write_fix('uV', 8); end % Unit
    
    % Phys Min/Max and Dig Min/Max
    pMin = zeros(1, ns); pMax = zeros(1, ns);
    for j = 1:ns
        pMin(j) = min(signals{j}); pMax(j) = max(signals{j});
        if pMin(j) == pMax(j), pMax(j) = pMin(j) + 1; end
    end
    for j = 1:ns, write_fix(sprintf('%.4f', pMin(j)), 8); end
    for j = 1:ns, write_fix(sprintf('%.4f', pMax(j)), 8); end
    for j = 1:ns, write_fix('-32768', 8); end
    for j = 1:ns, write_fix('32767', 8); end
    for j = 1:ns, write_fix('', 80); end % Prefilter
    for j = 1:ns, write_fix(num2str(fs(j) * recordDuration), 8); end
    for j = 1:ns, write_fix('', 32); end % Reserved
    
    % --- Data Records ---
    for rec = 1:numDataRecords
        for j = 1:ns
            startIdx = (rec - 1) * fs(j) + 1;
            endIdx = rec * fs(j);
            seg = zeros(fs(j), 1);
            if startIdx <= length(signals{j})
                actual = signals{j}(startIdx:min(endIdx, length(signals{j})));
                seg(1:length(actual)) = actual;
            end
            % Scale to int16
            scale = 65535 / (pMax(j) - pMin(j));
            digitalValues = round((seg - pMin(j)) * scale - 32768);
            digitalValues = max(-32768, min(32767, digitalValues));
            fwrite(fid, digitalValues, 'int16');
        end
    end
    fclose(fid);
end
