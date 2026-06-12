% pulse annotate script
clear;
close all;
props = readProps('config.txt');
TemplateMarking_input = props('TemplateMarking_input');
TemplateMarking_output = props('TemplateMarking_output');
TemplateMarking_TemplatesToShow = str2num(props('TemplateMarking_TemplatesToShow'));
Skip_Existing = logical(str2num(props('Skip_Existing')));

%%
[analysisFiles] = TemplateSetup(TemplateMarking_input, TemplateMarking_output);
disp('*********************************************************************');

time = 0;
completed = 1;
for i = 1:size(analysisFiles, 1)

    if isfile(fullfile(TemplateMarking_output, [analysisFiles{i,1} '_template_markings.mat'])) && Skip_Existing
        disp([analysisFiles{i,1} '_template_markings.mat exists skipping because Skip_Existing = 1 in config.']);
        continue
    end
    tStart = tic;
    disp(['Showing ', analysisFiles{i,1}, ' | ' num2str(i) ' of ' num2str(size(analysisFiles,1)) ' | Remaining: ' num2str(size(analysisFiles,1)-i)]);
    avg_time = time/completed;
    disp(['Avg Time (s): ' num2str(avg_time)]);

    disp(['Est finish (min): ' num2str((avg_time*(length(analysisFiles)-i))/3600)]);
    disp(join(['Output loc ' TemplateMarking_output]));

    template_info = load(analysisFiles{i, 2});
    template_info = template_info.template_info;
    template_indexs = cell(length(template_info),1);
    if TemplateMarking_TemplatesToShow == inf || TemplateMarking_TemplatesToShow < 1
        [template_indexs, closearg] = template_viewer(template_info, analysisFiles{i,1}, analysisFiles{i,4});
        if closearg == 0
           disp('Bye!');
           return 
        end
    else
        runs = ceil(length(template_info)/TemplateMarking_TemplatesToShow);
        beginidx = find(mod(1:length(template_info),TemplateMarking_TemplatesToShow)==0);
        prev = 0;
        for r = 1:runs
            begin = prev+1;
            endidx = begin + TemplateMarking_TemplatesToShow-1;

            if endidx > length(template_info)
                endidx = length(template_info);
            end
            if begin > endidx
               begin = endidx; 
            end
           
            [template_indexs(begin:endidx), closearg] = template_viewer(template_info(begin:endidx), analysisFiles{i,1}, analysisFiles{i,4});

            if closearg == 0
               disp('Bye!');
               return 
            end
            prev = endidx;
        end
    end
    
    for t = 1:size(template_info)
        template_info{t}.TemplateBad = template_indexs{t}.TemplateBad;
        template_info{t}.bad_r_templates = template_indexs{t}.bad_r_templates; 
        template_info{t}.bad_ppg_templates = template_indexs{t}.bad_ppg_templates;
        if ~isnan(template_indexs{t}.Dicrotic)
            if template_indexs{t}.Dicrotic < 0
                template_info{t}.Dicrotic = 1;
            elseif template_indexs{t}.Dicrotic > length(template_info{t}.ppgTemplate)
                template_info{t}.Dicrotic = length(template_info{t}.ppgTemplate);
            else
                template_info{t}.Dicrotic = template_indexs{t}.Dicrotic;
            end
        else
            template_info{t}.Dicrotic = template_indexs{t}.Dicrotic;
        end

        if isempty(template_info{t}.ecgTemplate)
            template_info{t}.bad_r_templates = 1;
        end
        
        if isempty(template_info{t}.ppgTemplate)
            template_info{t}.bad_ppg_templates = 1;
            template_info{t}.TemplateBad = 1;
        end

        if ~template_info{t}.bad_ppg_templates
            template_info{t}.Onset = 0;
            [~, template_info{t}.Peak] = max(template_info{t}.ppgTemplate);
            template_info{t}.Peak = template_info{t}.Peak -1;
            template_info{t}.End = length(template_info{t}.ppgTemplate) - 1;
        else
            template_info{t}.Onset = nan;
            template_info{t}.Peak = nan;
            template_info{t}.End = nan;
        end
    end
    
    disp('Saving template info...');
    save(fullfile(analysisFiles{i,4}, [analysisFiles{i, 1} '_template_markings']), 'template_info');
    %movefile(analysisFiles{i, 2}, '/hdd/data/mesa/Manuscript 1/Handmarked Rs/4 template Generation');
    disp(['This review took ' num2str(toc(tStart)) ' sec']);
    %% finished
    disp(['____________________________________________________________________________________________________' newline]);
    time = time + toc(tStart);
    completed = completed+1;
end