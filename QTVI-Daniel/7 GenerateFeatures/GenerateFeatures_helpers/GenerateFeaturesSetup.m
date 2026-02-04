function [analysisFiles] = GenerateFeaturesSetup(annealedseg_path, wave_data_path, template_path, type)

    if ~type
        annealedseg_list = dir(fullfile(annealedseg_path, '**/*_annealedSegments.mat'));
        anneal = parse_files_for_ids(annealedseg_list,'_annealedSegments');

        wavedata_list = dir(fullfile(wave_data_path, '**/*_wave_data.mat'));     
        wave = parse_files_for_ids(wavedata_list,'_wave_data');
        
        i1 = intersect(wave,anneal,'stable');
        smaller = 0;
        num = 0;
        if length(i1) < length(anneal)
            disp(['Annealed Segments list is smaller then others, defaulting to use only matched ids from: ' annealedseg_path]);
            disp(['*_annealedSegments.mat list length: ' num2str(length(anneal))]);
            disp(['*_wave_data.mat list length: ' num2str(length(wave))]);
            disp(['*_template_markings.mat list length: ' num2str(length(temps))]);
            smaller = 1;
            num = length(anneal);
        end
        if length(i1) < length(wave) && ~smaller
            disp(['Wave data list is smaller then others, defaulting to use only matched ids from: ' wavedata_list]);
             disp(['*_annealedSegments.mat list length: ' num2str(length(anneal))]);
            disp(['*_wave_data.mat list length: ' num2str(length(wave))]);
            disp(['*_template_markings.mat list length: ' num2str(length(temps))]);
            num = length(wave);
        end
        
        analysisFiles = cell(length(i1), 4);
        for i = 1:length(i1)
            uuid = i1{i};
            analysisFiles{i, 1} = uuid;
            analysisFiles{i, 2} = fullfile(annealedseg_path, [uuid '_annealedSegments.mat']);
            analysisFiles{i, 3} = fullfile(wave_data_path, [uuid '_wave_data.mat']);
            analysisFiles{i, 4} = fullfile(template_path, [uuid '_template_markings.mat']);
        end
    else
        template_list = dir(fullfile(template_path, '**/*_template_markings.mat'));
        if isempty(template_list) 
           error('Gen_Features_use_manually_reviewed_templates = 1 in config, expected but couldn''t find any *_template_markings.mat files in path set in config.');
        end
        temps = parse_files_for_ids(template_list,'_template_markings');

        annealedseg_list = dir(fullfile(annealedseg_path, '**/*_annealedSegments.mat'));
        anneal = parse_files_for_ids(annealedseg_list,'_annealedSegments');

        wavedata_list = dir(fullfile(wave_data_path, '**/*_wave_data.mat'));     
        wave = parse_files_for_ids(wavedata_list,'_wave_data');
        
        
        i1 = intersect(wave,anneal,'stable');
        smaller = 0;
        num = 0;
        if length(i1) < length(anneal)
            disp(['Annealed Segments list is smaller then others, defaulting to use only matched ids from: ' annealedseg_path]);
            disp(['*_annealedSegments.mat list length: ' num2str(length(anneal))]);
            disp(['*_wave_data.mat list length: ' num2str(length(wave))]);
            disp(['*_template_markings.mat list length: ' num2str(length(temps))]);
            smaller = 1;
            num = length(anneal);
        end
        if length(i1) < length(wave) && ~smaller
            disp(['Wave data list is smaller then others, defaulting to use only matched ids from: ' wavedata_list]);
             disp(['*_annealedSegments.mat list length: ' num2str(length(anneal))]);
            disp(['*_wave_data.mat list length: ' num2str(length(wave))]);
            disp(['*_template_markings.mat list length: ' num2str(length(temps))]);
            smaller = 1;
            num = length(wave);
        end
        i2 = intersect(i1,temps,'stable');

        if length(i2) < length(temps) && ~smaller
            disp(['Template list is smaller then others, defaulting to use only matched ids from: ' template_list]);
             disp(['*_annealedSegments.mat list length: ' num2str(length(anneal))]);
            disp(['*_wave_data.mat list length: ' num2str(length(wave))]);
            disp(['*_template_markings.mat list length: ' num2str(length(temps))]);
            smaller = 1;
            num = length(temps);
        end
        
        analysisFiles = cell(length(i2), 4);
        for i = 1:length(i2)
            uuid = i2{i};
            analysisFiles{i, 1} = uuid;
            analysisFiles{i, 2} = fullfile(annealedseg_path, [uuid '_annealedSegments.mat']);
            analysisFiles{i, 3} = fullfile(wave_data_path, [uuid '_wave_data.mat']);
            analysisFiles{i, 4} = fullfile(template_path, [uuid '_template_markings.mat']);
        end
    end
end

function lst = parse_files_for_ids(file_list, regstr)
    lst = {};
    for i = 1:length(file_list)
        [~, name, ~] = fileparts(file_list(i).name);
        start_idx = regexp(name, regstr);
        start_idx = start_idx(1);
        name = name(1:start_idx-1);
        lst{end+1} = name;
    end
end