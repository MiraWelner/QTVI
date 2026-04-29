function [ret] = FindWaveBounds(anneal_path, output_path)
%     try
        [~, name, ~] = fileparts(anneal_path);
        start_idx = regexp(name, '_annealedSegments');
        start_idx = start_idx(1);
        name = name(1:start_idx-1);
        annealedSegments = load(anneal_path);
        annealedSegments = annealedSegments.annealedSegments;

        % ====================================================================
        % Realign each bin's data with the C++ pipeline's bin geometry,
        % then upsample 256 -> 1000 Hz.
        %
        % The MATLAB annealer wrote bins of 15361 samples placed back-to-
        % back at stride 15361. So each bin starts (n-1) samples too late
        % relative to C++'s clean stride of 15360 (= 60 s exactly), and
        % consecutive bins drift forward by ~3.9 ms per bin in absolute
        % time.
        %
        % Strategy (mirrors C++):
        %   1. Concatenate every bin back into the full recording at
        %      256 Hz. MATLAB's bins tile the source with no gap or
        %      overlap, so this exactly reconstructs the continuous
        %      signal.
        %   2. Resample the FULL signal once (256 -> 1000 Hz). This is
        %      what C++ does at file ingest. Doing it once means every
        %      bin sees real signal context on both sides instead of
        %      zero-padding -- which is what caused trailing-edge
        %      filter transients and spurious R-peaks at the end of
        %      each bin in earlier per-bin-resample versions of this
        %      function.
        %   3. Re-slice the upsampled signal at stride 60000 samples
        %      (= 60 s at 1000 Hz). Bin n now covers upsampled-sample
        %      indices [(n-1)*60000 + 1, n*60000] -- exactly C++'s
        %      geometry.
        %
        % Number of bins is preserved (or shrinks by at most 1 if the
        % final bin is partial after the realignment). Bin-index arrays
        % and sample rates are scaled to the new 1000-Hz coordinates.
        % ====================================================================
        orig_fs   = 256;
        target_fs = 1000;
        [p, q]    = rat(target_fs / orig_fs);   % 125 / 32
        scale     = target_fs / orig_fs;

        bin_size_target = target_fs * 60;        % 60000 samples per 60-s bin

        nSeg = length(annealedSegments);

        % --- Concatenate all bins (256 Hz) ---
        ecg_lens = zeros(nSeg, 1);
        ppg_lens = zeros(nSeg, 1);
        for i = 1:nSeg
            seg = annealedSegments{i};
            ecg_lens(i) = numel(seg.ecg);
            ppg_lens(i) = numel(seg.po);
        end
        ecg_concat = zeros(sum(ecg_lens), 1);
        ppg_concat = zeros(sum(ppg_lens), 1);
        ecg_off = 0;
        ppg_off = 0;
        for i = 1:nSeg
            seg = annealedSegments{i};
            n_e = ecg_lens(i);
            n_p = ppg_lens(i);
            ecg_concat(ecg_off + (1:n_e)) = seg.ecg(:);
            ppg_concat(ppg_off + (1:n_p)) = seg.po(:);
            ecg_off = ecg_off + n_e;
            ppg_off = ppg_off + n_p;
        end

        % --- Resample the FULL recording once (no per-bin transients) ---
        ecg_full_up = resample(ecg_concat, p, q);
        ppg_full_up = resample(ppg_concat, p, q);

        % --- Determine how many full 60-s bins fit at 1000 Hz ---
        n_full_ecg = floor(numel(ecg_full_up) / bin_size_target);
        n_full_ppg = floor(numel(ppg_full_up) / bin_size_target);
        n_new = min([nSeg, n_full_ecg, n_full_ppg]);

        % --- Re-slice the upsampled signal at the C++ stride ---
        new_segments = cell(1, n_new);
        for i = 1:n_new
            seg = annealedSegments{i};

            a = (i - 1) * bin_size_target + 1;
            b = i * bin_size_target;

            seg.ecg = ecg_full_up(a:b);
            seg.po  = ppg_full_up(a:b);

            % Scale sample indices to the new 1000-Hz rate.
            seg.ecg_bin_indexs = round(seg.ecg_bin_indexs * scale);
            seg.ppg_bin_indexs = round(seg.ppg_bin_indexs * scale);

            ecg_len = numel(seg.ecg);
            ppg_len = numel(seg.po);
            seg.ecg_bin_indexs(seg.ecg_bin_indexs > ecg_len) = ecg_len;
            seg.ppg_bin_indexs(seg.ppg_bin_indexs > ppg_len) = ppg_len;

            seg.ecgSampleRate = target_fs;
            seg.ppgSampleRate = target_fs;

            new_segments{i} = seg;
        end

        annealedSegments = new_segments;

        disp('Finding individual beats...');
        [wave_data] = FindWaveBounds_EKGandPPG(annealedSegments, 0, 1);
        disp('Saving...');
        save(fullfile(output_path, [name '_wave_data']), 'wave_data');
        ret = 1;
%     catch e
%         cprintf('err', ['Error in file: ', name, newline]);
%         input.anneal_path = anneal_path;
%         st = dbstack;
%         namestr = [st.name](http://st.name);
%         LogError(namestr, output_path, input, e);
%         ret = 0;
%     end
end