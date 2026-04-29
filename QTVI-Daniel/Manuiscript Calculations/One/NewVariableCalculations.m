% calculate new D/B variables, See "Summary of what needs to be done in the
% paper" email on 7/15/2021

files = dirWithoutDots('D:\data\deepLab\central data\MESA\Manuscript 1\9 individualFeatures');
output = "D:\data\deepLab\central data\MESA\Manuscript 1\newfeatures";
for f = 1:length(files)
    
    name = files(f).name;
%     if ~strcmp('7014457_20111102',name)
%         continue
%     end
%     if isfolder(fullfile(output,name,[name '-D.mat']))
%         continue
%     end
    begin = load(fullfile(files(f).folder, name, [name '-indx_intrabegin.mat']));
    begin = fillmissing(fillmissing(begin.tmp,'movmedian',20),'nearest');
    tp20 = load(fullfile(files(f).folder, name, [name '-indx_intratp20.mat']));
    tp20 = fillmissing(fillmissing(tp20.tmp,'movmedian',20),'constant',nanmean(tp20.tmp));
    r2tp20 = load(fullfile(files(f).folder, name, [name '-msec_r2tp20.mat']));
    r2tp20 = fillmissing(fillmissing(r2tp20.tmp,'movmedian',20),'constant',nanmean(r2tp20.tmp));
    time =  load(fullfile(files(f).folder, name, [name '-seco_time2sleep.mat']));
    time = sort(abs(time.tmp-time.tmp(1)));
    
    D = begin + tp20;
    D = diff(D); D = [D;D(end)]; % to make same length rep end. 1 beat wont throw things off.
    D = D/256; % to sec
    D = D*1000;
    
    A = r2tp20;
    
    B = A+D;
    

    parsave(A,B,D,time,name,files(f).folder);
    a=1;
end

function parsave(A,B,D,time,name,folder)
    output = "D:\data\deepLab\central data\MESA\Manuscript 1\newfeatures";
    if ~isfolder(fullfile(output,name))
       mkdir(fullfile(output,name)); 
    end
    tmp = A;
    save(fullfile(output,name,[name '-A.mat']),"tmp");
    tmp = B;
    save(fullfile(output,name,[name '-B.mat']),"tmp");
    tmp = D;
    save(fullfile(output,name,[name '-D.mat']),"tmp");   
    tmp = time;
    save(fullfile(output,name,[name '-time.mat']),"tmp");
    
    copyfile(fullfile(folder,name,[name '-slee_sleepstate.mat']),fullfile(output,name,[name '-slee_sleepstate.mat']));
    copyfile(fullfile(folder,name,[name '-seco_time2sleep.mat']),fullfile(output,name,[name '-seco_time2sleep.mat']));
    copyfile(fullfile(folder,name,[name '-seco_time2wake.mat']),fullfile(output,name,[name '-seco_time2wake.mat']));
end