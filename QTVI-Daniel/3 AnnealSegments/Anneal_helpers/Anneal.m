function [ret] = Anneal(mat_path, noise_path, output_path)

    try
        [~, name, ~] = fileparts(mat_path);
        name = regexp(name, '(\d+_\d+)', 'match');
        load(mat_path, 'ecg', 'ecgSamplingRate', 'ppg', 'ppgSamplingRate', 'scoring_epoch_size_sec', 'sleepStages');

        % should make noise gui save as variables not struct.
        noiseSEG = load(noise_path);
        noiseSEG = noiseSEG.noiseSEG;

        [annealedSegments, ~] = AnnealSegments(ppg, ppgSamplingRate, ecg, ecgSamplingRate, noiseSEG, scoring_epoch_size_sec, sleepStages, 15, 0);
        save(fullfile(output_path, [name{1} '_annealedSegments']), 'annealedSegments');
        ret = 1;
    catch e
        input.mat_path = mat_path;
        input.noise_path = noise_path;
        st = dbstack;
        namestr = st.name;
        LogError(namestr, output_path, input, e);
        ret = 0;
    end

end
