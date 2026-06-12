clear
inputloc = 'D:\data\deepLab\central data\MESA\Manuscript 1\9 individualFeatures nanless';
e_len = 300;

idno = [];
epoch_num = [];
t2sleep = [];
t2wake = [];
raw_time = [];
epochclosest2sleep = [];
epochclosest2wake = [];
timefromepochbegin_tosleeponset_ = [];
timefromepochend_tosleeponset_se = [];
epoch_len_sec = [];
ss0 = [];
ss1 = [];
ss2 = [];
ss3 = [];
ss4 = [];



file_list = dirWithoutDots(inputloc);
for ff = 1:length(file_list)
    if ff < 10
        continue
    end
disp([num2str(ff) ' of ' num2str(length(file_list))]);
file = fullfile( file_list(ff).folder, file_list(ff).name, [file_list(ff).name '-seco_time2sleep.mat']);
tmp = load(file);
time2sleep = tmp.tmp*-1;

file = fullfile( file_list(ff).folder, file_list(ff).name, [file_list(ff).name '-seco_time2wake.mat']);
tmp = load(file);
time2wake = tmp.tmp;  % don't do negative to match whay time2sleep is

file = fullfile( file_list(ff).folder, file_list(ff).name, [file_list(ff).name '-slee_sleepstate.mat']);
tmp = load(file);
sleepstates = tmp.tmp;

file = fullfile(file_list(ff).folder, file_list(ff).name, [file_list(ff).name '-seco_time2sleep.mat']);
tmp = load(file);
time = tmp.tmp*-1;
time = time/3600;

epochtimes = [0];
%% find epoch closest to sleep
sleepidx = closest_idx(time2sleep, 0);
tmp = (time2sleep) - time2sleep(1);
sleeptime = tmp(sleepidx);
max_time = max(tmp);

next = e_len;

while (max_time-e_len) > next
    epochtimes = [epochtimes; next];
    next = next + e_len;
end

if next < max_time
    epochtimes = [epochtimes; next];
    epochtimes = [epochtimes; max_time];
else
    epochtimes = [epochtimes; max_time];
end
rtime = epochtimes;
eidx = find(epochtimes>=sleeptime,1);
epochtimes = epochtimes-epochtimes(eidx);
epochtimes = epochtimes/3600;
sleepidx = closest_idx(epochtimes, 0);

if sleepidx >= length(epochtimes)
    sleepidx = sleepidx-1;
end

%% find epoch closest to wake

waketimes = [0];
wakeidx = closest_idx(time2wake, 0);
% plot(time2wake);hold on;hline(0);scatter(wakeidx,time2wake(wakeidx));

tmp = (time2wake) - time2wake(1);
waketime = tmp(wakeidx);
max_time = max(tmp);

next = e_len;

while (max_time-e_len) > next
    waketimes = [waketimes; next];
    next = next + e_len;
end

if next < max_time
    waketimes = [waketimes; next];
    waketimes = [waketimes; max_time];
else
    waketimes = [waketimes; max_time];
end

eidx = find(waketimes>=waketime,1);
waketimes = waketimes-waketimes(eidx);
waketimes = waketimes/3600;
wakeidx = closest_idx(waketimes, 0);
if wakeidx >= length(waketimes)
    wakeidx = wakeidx-1;
end



for x = 1:length(epochtimes)-1
    b = epochtimes(x);
    e = epochtimes(x+1);
    mask = (time>=b) & (time<e);
    cut = sleepstates(mask);
    ss0 = [ss0; sum(cut==0)/length(cut)];
    ss1 = [ss1; sum(cut==1)/length(cut)];
    ss2 = [ss2; sum(cut==2)/length(cut)];
    ss3 = [ss3; sum(cut==3)/length(cut)];
    ss4 = [ss4; sum(cut==4)/length(cut)];
end




parts = split(file_list(ff).name, "_");
id = str2num(parts{1});


idno = [idno; zeros(length(epochtimes)-1,1) + id];
epoch_num = [epoch_num; (1:length(epochtimes)-1)'];
epoch_len_sec = [epoch_len_sec; diff(epochtimes*3600)];
raw_time = [raw_time; rtime(1:length(epochtimes)-1)];
t2sleep = [t2sleep; epochtimes(1:end-1)];
t2wake = [t2wake; waketimes(1:end-1)];

tmp = zeros(length(epochtimes)-1,1);
tmp(sleepidx) = 1;
epochclosest2sleep = [epochclosest2sleep; tmp];
tmp = zeros(length(epochtimes)-1,1);
tmp(wakeidx) = 1;
epochclosest2wake = [epochclosest2wake; tmp];

timefromepochbegin_tosleeponset_ = [timefromepochbegin_tosleeponset_; epochtimes(1:end-1)];
timefromepochend_tosleeponset_se = [timefromepochend_tosleeponset_se; epochtimes(2:end)];
% t = array2table([idno, epoch_num, epoch_len_sec, raw_time, t2sleep, t2wake, epochclosest2sleep, epochclosest2wake, timefromepochbegin_tosleeponset_, timefromepochend_tosleeponset_se, ss0, ss1, ss2, ss3, ss4],'VariableNames',{'idno','epoch_num','epoch_len_sec','recordtimesec','timetosleephr','timetowakehr','epochclosest2sleep','epochclosest2wake','timefromepochbegin_tosleeponset_','timefromepochend_tosleeponset_se', 'ss0', 'ss1', 'ss2', 'ss3', 'ss4'});

end

t = array2table([idno, epoch_num, epoch_len_sec, raw_time, t2sleep, t2wake, epochclosest2sleep, epochclosest2wake, timefromepochbegin_tosleeponset_, timefromepochend_tosleeponset_se, ss0, ss1, ss2, ss3, ss4],'VariableNames',{'idno','epoch_num','epoch_len_sec','recordtimesec','timetosleephr','timetowakehr','epochclosest2sleep','epochclosest2wake','timefromepochbegin_tosleeponset_','timefromepochend_tosleeponset_se', 'ss0', 'ss1', 'ss2', 'ss3', 'ss4'});
writetable(t,'C:\Users\dan\OneDrive - University of Cincinnati\projects\000001-ecg_ppg_coupling\src\run\v7\data\entropy results\times.csv')