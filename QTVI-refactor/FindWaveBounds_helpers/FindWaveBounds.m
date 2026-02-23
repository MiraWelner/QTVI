function [ret] = FindWaveBounds(anneal_path, output_path)
%     try
        [~, name, ~] = fileparts(anneal_path);
        
        start_idx = regexp(name, '_annealedSegments');
        start_idx = start_idx(1);
        name = name(1:start_idx-1);
        
        annealedSegments = load(anneal_path);
        annealedSegments = annealedSegments.annealedSegments;

        disp('Finding individual beats...');
        [wave_data] = FindWaveBounds_EKGandPPG(annealedSegments, 0, 1);
        disp('Saving...');
        save(fullfile(output_path, [name '_wave_data']), 'wave_data');

        ret = 1;
%     catch e
%         cprintf('err', ['Error in file: ', name, newline]);
%         input.anneal_path = anneal_path;
%         st = dbstack;
%         namestr = st.name;
%         LogError(namestr, output_path, input, e);
%         ret = 0;
%     end

end
