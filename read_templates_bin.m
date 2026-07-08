function [bins, saecg] = read_templates_bin(path)
%READ_TEMPLATES_BIN  Load a QTVi *_templates.bin file.
%
%   [bins, saecg] = READ_TEMPLATES_BIN(PATH) returns:
%
%     bins   struct array, one element per bin, with fields:
%       bin_index                            (0-based index)
%       bad_segment                          (logical)
%       ch1_raw, ch1_squared, ch1_absval, ch1_unfiltered      struct
%       ch2_raw, ch2_squared, ch2_absval, ch2_unfiltered      struct
%       ch3_raw, ch3_squared, ch3_absval, ch3_unfiltered      struct
%       ppg               (Mx1 vector; empty if no PPG)
%       ppg_std           (Mx1 vector; empty if no PPG or not computed)
%
%     Each method struct (e.g. bins(k).ch1_raw) has:
%       ecgTemplate       (Nx1 vector)   averaged heartbeat waveform
%       ecgTemplate_std   (Nx1 vector)   per-sample std (empty for
%                                        squared/absval/unfiltered methods
%                                        -- only raw computes std)
%       alignment_point   (scalar)
%       avg_r_expand      (scalar)
%
%     saecg  struct with recording-wide SAECG averages. Fields:
%       ch1_raw, ch1_squared, ch1_absval, ch1_unfiltered      struct
%       ch2_raw, ch2_squared, ch2_absval, ch2_unfiltered      struct
%       ch3_raw, ch3_squared, ch3_absval, ch3_unfiltered      struct
%       ppg                                                   struct
%
%     Each SAECG struct has:
%       waveform          (Nx1 vector)
%       n_contributing    (scalar, count of bins contributing to the avg)
%
%   Layout (from template_io.hpp):
%     [uint64 n_bins]
%     For each bin:
%       12 method blocks (ch1_raw, ch1_squared, ch1_absval, ch1_unfiltered,
%                         ch2_raw, ..., ch3_unfiltered)
%         Each block:
%           [uint64 sz][sz doubles ecgTemplate]
%           [uint64 sz][sz doubles ecgTemplate_std]     (sz=0 for non-raw)
%           [double alignment_point]
%           [double avg_r_expand]
%       [uint64 sz][sz doubles ppgTemplate]
%       [uint64 sz][sz doubles ppgTemplate_std]
%       [uint8 bad_segment]
%     SAECG tail (13 averaged waveforms in fixed order):
%       [uint64 sz][sz doubles waveform]
%       [uint64 n_contributing]
%   Repeat SAECG block 13 times.
%
%   Example:
%     [bins, saecg] = read_templates_bin('Z1011347_..._15_templates.bin');
%     plot(bins(2).ch1_raw.ecgTemplate);           % bin 2's ch1 raw template
%     plot(saecg.ch1_raw.waveform);                % recording-wide avg for ch1
%     fprintf('SAECG ch1_raw n_contributing = %d\n', saecg.ch1_raw.n_contributing);

    f = fopen(path, 'rb');
    if f < 0
        error('read_templates_bin:cannotOpen', 'Cannot open file: %s', path);
    end
    cleanupObj = onCleanup(@() fclose(f));

    % Read n_bins.
    nBinsRaw = fread(f, 1, '*uint64');
    if isempty(nBinsRaw)
        error('read_templates_bin:truncated', 'File too short: %s', path);
    end
    nBins = double(nBinsRaw);

    METHOD_NAMES = {'ch1_raw', 'ch1_squared', 'ch1_absval', 'ch1_unfiltered', ...
                    'ch2_raw', 'ch2_squared', 'ch2_absval', 'ch2_unfiltered', ...
                    'ch3_raw', 'ch3_squared', 'ch3_absval', 'ch3_unfiltered'};

    % Preallocate bins struct array.
    bins = repmat(makeEmptyBin(METHOD_NAMES), nBins, 1);

    for k = 1:nBins
        bins(k).bin_index = k - 1;   % 0-based to match the C++ side

        % 12 ECG method blocks.
        for m = 1:numel(METHOD_NAMES)
            bins(k).(METHOD_NAMES{m}) = readMethodBlock(f);
        end

        % PPG template + std.
        bins(k).ppg     = readVecD(f);
        bins(k).ppg_std = readVecD(f);

        % bad_segment flag (uint8).
        bad = fread(f, 1, '*uint8');
        bins(k).bad_segment = ~isempty(bad) && bad ~= 0;
    end

    % SAECG tail: 13 averaged waveforms in the same order (12 ECG + PPG).
    saecg = struct();
    for m = 1:numel(METHOD_NAMES)
        saecg.(METHOD_NAMES{m}) = readAveraged(f);
    end
    saecg.ppg = readAveraged(f);
end

function s = makeEmptyBin(methodNames)
    s = struct('bin_index', 0, 'bad_segment', false, 'ppg', [], 'ppg_std', []);
    for m = 1:numel(methodNames)
        s.(methodNames{m}) = struct( ...
            'ecgTemplate', [], 'ecgTemplate_std', [], ...
            'alignment_point', 0, 'avg_r_expand', 0);
    end
end

function m = readMethodBlock(f)
    %READMETHODBLOCK  Read one ChannelMethodTemplate block.
    m.ecgTemplate     = readVecD(f);
    m.ecgTemplate_std = readVecD(f);
    ap = fread(f, 1, '*double');
    ar = fread(f, 1, '*double');
    if isempty(ap) || isempty(ar)
        error('read_templates_bin:truncated', 'File truncated in method block');
    end
    m.alignment_point = double(ap);
    m.avg_r_expand    = double(ar);
end

function v = readVecD(f)
    %READVECD  Read [uint64 sz][sz doubles].
    sz = fread(f, 1, '*uint64');
    if isempty(sz)
        error('read_templates_bin:truncated', 'File truncated at vector size');
    end
    n = double(sz);
    if n == 0
        v = [];
    else
        v = fread(f, n, '*double');
        if length(v) < n
            error('read_templates_bin:truncated', ...
                  'Vector truncated: expected %d doubles, got %d', n, length(v));
        end
    end
end

function a = readAveraged(f)
    %READAVERAGED  Read one AveragedTemplate: [uint64 sz][sz doubles][uint64 n_contributing].
    a.waveform = readVecD(f);
    nc = fread(f, 1, '*uint64');
    if isempty(nc)
        error('read_templates_bin:truncated', 'File truncated at n_contributing');
    end
    a.n_contributing = double(nc);
end