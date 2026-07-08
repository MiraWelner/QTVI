function T = read_template_markings(path)
%READ_TEMPLATE_MARKINGS  Load a QTVi *_template_markings.bin file.
%
%   T = READ_TEMPLATE_MARKINGS(PATH) returns a table with one row per bin.
%
%   Columns:
%     bin_index                            uint64
%     bad_r_ch1, bad_r_ch2, bad_r_ch3      logical
%     ppg_issue                            uint8   (0 = ok, 1 = bad, 2 = no ppg)
%     p_begin_ch1..3, q_begin_ch1..3       int32 (sample index; -1 = unmarked)
%     t_begin_ch1..3, t_end_ch1..3         int32 (sample index; -1 = unmarked)
%     ppg_onset, ppg_peak, ppg_dicrotic,
%     ppg_50, ppg_end                      int32 (sample index; -1 = unmarked)
%
%   Layout (must match writeTemplateMarkingsBin in TemplateBinIO.hpp):
%     uint64   numBins
%     per bin (80 bytes):
%       uint64  index                     (offset  0, size 8)
%       uint8   bad_r_ch1..3              (offset  8, size 3)
%       uint8   ppg_issue                 (offset 11, size 1)
%       int32   p_begin_ch1..3            (offset 12, size 12)
%       int32   q_begin_ch1..3            (offset 24, size 12)
%       int32   t_begin_ch1..3            (offset 36, size 12)
%       int32   t_end_ch1..3              (offset 48, size 12)
%       int32   ppg_onset, ppg_peak,
%               ppg_dicrotic, ppg_50,
%               ppg_end                   (offset 60, size 20)
%
%   This reader loads the entire file into a byte buffer and typecasts
%   fields at fixed offsets, rather than relying on fread positioning
%   through mixed integer types (which can misbehave in MATLAB on Windows).

    f = fopen(path, 'rb');
    if f < 0
        error('read_template_markings:cannotOpen', 'Cannot open file: %s', path);
    end
    cleanupObj = onCleanup(@() fclose(f));
    bytes = fread(f, inf, '*uint8');

    if length(bytes) < 8
        error('read_template_markings:truncated', ...
              'File too short (%d bytes) to contain nBins header', length(bytes));
    end

    numBins = typecast(bytes(1:8), 'uint64');
    n = double(numBins);
    if n == 0
        T = table();
        return
    end

    BIN_SIZE = 80;
    expectedTotal = 8 + n * BIN_SIZE;
    if length(bytes) < expectedTotal
        error('read_template_markings:truncated', ...
              'File has %d bytes; expected %d for %d bins', ...
              length(bytes), expectedTotal, n);
    end
    if length(bytes) > expectedTotal
        warning('read_template_markings:trailingBytes', ...
                '%d unexpected trailing bytes in %s', ...
                length(bytes) - expectedTotal, path);
    end

    bin_index    = zeros(n, 1, 'uint64');
    bad_r_ch1    = false(n, 1);
    bad_r_ch2    = false(n, 1);
    bad_r_ch3    = false(n, 1);
    ppg_issue    = zeros(n, 1, 'uint8');

    p_begin      = zeros(n, 3, 'int32');
    q_begin      = zeros(n, 3, 'int32');
    t_begin      = zeros(n, 3, 'int32');
    t_end_       = zeros(n, 3, 'int32');

    ppg_onset    = zeros(n, 1, 'int32');
    ppg_peak     = zeros(n, 1, 'int32');
    ppg_dicrotic = zeros(n, 1, 'int32');
    ppg_50       = zeros(n, 1, 'int32');
    ppg_end      = zeros(n, 1, 'int32');

    for i = 1:n
        base = 8 + (i - 1) * BIN_SIZE + 1;

        bin_index(i)    = typecast(bytes(base    : base +  7), 'uint64');
        bad_r_ch1(i)    = bytes(base +  8) ~= 0;
        bad_r_ch2(i)    = bytes(base +  9) ~= 0;
        bad_r_ch3(i)    = bytes(base + 10) ~= 0;
        ppg_issue(i)    = bytes(base + 11);

        p_begin(i, :)   = typecast(bytes(base + 12 : base + 23), 'int32');
        q_begin(i, :)   = typecast(bytes(base + 24 : base + 35), 'int32');
        t_begin(i, :)   = typecast(bytes(base + 36 : base + 47), 'int32');
        t_end_(i, :)    = typecast(bytes(base + 48 : base + 59), 'int32');

        ppg_onset(i)    = typecast(bytes(base + 60 : base + 63), 'int32');
        ppg_peak(i)     = typecast(bytes(base + 64 : base + 67), 'int32');
        ppg_dicrotic(i) = typecast(bytes(base + 68 : base + 71), 'int32');
        ppg_50(i)       = typecast(bytes(base + 72 : base + 75), 'int32');
        ppg_end(i)      = typecast(bytes(base + 76 : base + 79), 'int32');
    end

    T = table(bin_index, bad_r_ch1, bad_r_ch2, bad_r_ch3, ppg_issue, ...
              p_begin(:,1), p_begin(:,2), p_begin(:,3), ...
              q_begin(:,1), q_begin(:,2), q_begin(:,3), ...
              t_begin(:,1), t_begin(:,2), t_begin(:,3), ...
              t_end_(:,1),  t_end_(:,2),  t_end_(:,3),  ...
              ppg_onset, ppg_peak, ppg_dicrotic, ppg_50, ppg_end, ...
        'VariableNames', { ...
            'bin_index', 'bad_r_ch1', 'bad_r_ch2', 'bad_r_ch3', 'ppg_issue', ...
            'p_begin_ch1', 'p_begin_ch2', 'p_begin_ch3', ...
            'q_begin_ch1', 'q_begin_ch2', 'q_begin_ch3', ...
            't_begin_ch1', 't_begin_ch2', 't_begin_ch3', ...
            't_end_ch1',   't_end_ch2',   't_end_ch3',   ...
            'ppg_onset', 'ppg_peak', 'ppg_dicrotic', 'ppg_50', 'ppg_end'});
end