% pulse annotate script
clear;
close all;

props = readProps('config.txt');
Prescreen_input = props('Prescreen_input');
Prescreen_output = props('Prescreen_output');



[analysisFiles] = TemplateSetup(Prescreen_input, Prescreen_output);
disp('*********************************************************************');

time = 0;
for i = 1:size(analysisFiles, 1)
    
    if isfile(fullfile(analysisFiles{i,4}, [analysisFiles{i, 1} '_template_markings']))
        continue
    end
    
    close all;
    tStart = tic;
    disp(['Showing ', analysisFiles{i,1}, ' | ' num2str(i) ' of ' num2str(length(analysisFiles))]);
    avg_time = time/i;
    disp(['Avg Time (s): ' num2str(avg_time)]);

    disp(['Est finish (min): ' num2str((avg_time*(length(analysisFiles)-i))/60)]);
    disp(join(['Output loc ' Prescreen_output]));

    template_info = load(analysisFiles{i, 2});
    template_info = template_info.template_info;
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    minimum= inf;
    maximum = -inf;
    for x = 1:length(template_info)
       minimum= min([template_info{x}.ppgTemplate, minimum]);
       maximum= max([template_info{x}.ppgTemplate, maximum]);
    end
    for x = 1:length(template_info)
        hold on;
        ma = max(template_info{x}.ppgTemplate);
        mi = min(template_info{x}.ppgTemplate);
        plot((template_info{x}.ppgTemplate-mi)./(ma-mi));
    end
    
     while true
        disp('Press Space to confirm ppg is ok otherwise press anything else...');
        try
            w = waitforbuttonpress;
        catch
            return
        end
        switch w
            case 1 % (keyboard press)
            key = get(gcf, 'currentcharacter');
            switch key
                case 32 % space
                    for t = 1:size(template_info)
                        template_info{t}.TemplateBad = 0;
                        template_info{t}.Onset = 1;
                        [~,idx] = max(template_info{x}.ppgTemplate);
                        template_info{t}.Peak = idx;
                        template_info{t}.Dicrotic = dumbDicrotic(template_info{x}.ppgTemplate);
                        template_info{t}.End = length(template_info{x}.ppgTemplate);
                    end
                    disp('Saving template info...');
                    save(fullfile(analysisFiles{i,4}, [analysisFiles{i, 1} '_template_markings']), 'template_info');
                    %movefile(analysisFiles{i, 2}, '/hdd/data/mesa/Manuscript 1/Handmarked Rs/4 template Generation');

                    break
                otherwise
                    break% Wait for a different command.
            end
        end
    end
    
    disp(['____________________________________________________________________________________________________' newline]);
    time = time + toc(tStart);
end