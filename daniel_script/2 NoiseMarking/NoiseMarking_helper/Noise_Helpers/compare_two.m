clear;

files = {"D:\git\PPG\src\matlab\noise marking standalone\OneDrive_2020-10-25\Grace 5 Practice Files - Copy\g6018904_20110621_noise.MAT",
    "D:\git\PPG\src\matlab\noise marking standalone\OneDrive_2020-10-25 (1)\Andrew Practice Files - Copy\a6018904_20110621_noise.MAT"};

edf = "C:\Users\dan\Desktop\OneDrive_1_10-25-2020\6018904_20110621\6018904_20110621.edf";


[edf_hdr, edf_data] = edfread(edf, 'verbose', 1, 'targetSignals', [1, 23]); % 1 is ecg index, 23 pulse

ppg = edf_data(2, :);
ecg = edf_data(1, :);
ecgSamplingRate = edf_hdr.frequency(1); % 256 FOR MESA SET
ppgSamplingRate = edf_hdr.frequency(23); % 256 FOR MESA SET


ppgTime_sec = 0:1 / ppgSamplingRate:(length(ppg) / ppgSamplingRate - 1 / ppgSamplingRate);
ecgTime_sec = 0:1 / ecgSamplingRate:(length(ecg) / ecgSamplingRate - 1 / ecgSamplingRate);

%     a1=subplot(2,2,1);
%     plot(ecgTime_sec,ecg);
%     a2=subplot(2,2,2);
%     plot(ecgTime_sec,ecg);
% 
%     a3=subplot(2,2,3);
%     plot(ppgTime_sec,ppg);
%     a4=subplot(2,2,4);
%     plot(ppgTime_sec,ppg);
% 
%     linkaxes([a1 a2],'y');
%     linkaxes([a3 a4],'y');
%     linkaxes([a1 a2 a3 a4],'x');

noiseA = load(files{1});
noiseB = load(files{2});

close all;
stepsize = 500;

beg = 0;
end_idx = stepsize;
x=1;
while end_idx < ecgTime_sec(end)
    cla;
    mask = ecgTime_sec >= beg & ecgTime_sec <= end_idx;
    [ha, pos] = tight_subplot(2,2,[.01 .03],[.1 .01],[.01 .01]);
    a1 = ha(1);
    a2 = ha(2);
    a3 = ha(3);
    a4 = ha(4);
    
    axes(a1);
    plot(ecgTime_sec(mask),ecg(mask));
    axes(a2);
    plot(ecgTime_sec(mask),ecg(mask));

    axes(a3);
    plot(ppgTime_sec(mask),ppg(mask));
    axes(a4);
    plot(ppgTime_sec(mask),ppg(mask));

    linkaxes([a1 a2],'y');
    linkaxes([a3 a4],'y');
    linkaxes([a1 a2 a3 a4],'x');
    
    for n = 1:size(noiseA.noise_markings,2)
        beg = noiseA.noise_markings(n,1);
        last = noiseA.noise_markings(n,2);
        if (noiseA.noise_markings(n,5) == 'l')
            d1 = patch(a3, [beg last last beg],[a3.YLim(2) a3.YLim(2) a3.YLim(1) a3.YLim(1)],  'g', 'LineStyle', 'none');
            alpha(d1, .15);
        elseif (noiseA.noise_markings(n,5) == 'r')
            d1 = patch(a1, [beg last last beg], [a1.YLim(2) a1.YLim(2) a1.YLim(1) a1.YLim(1)], 'c', 'LineStyle', 'none');
            alpha(d1, .15);
        else
            d1 = patch(a3, [beg last last beg], [a3.YLim(2) a3.YLim(2) a3.YLim(1) a3.YLim(1)], 'r', 'LineStyle', 'none');
            alpha(d1, .15);
            d1 = patch(a1, [beg last last beg], [a1.YLim(2) a1.YLim(2) a1.YLim(1) a1.YLim(1)], 'r', 'LineStyle', 'none');
            alpha(d1, .15);
        end
    end
    
    for n = 1:size(noiseB.noise_markings,2)
        beg = noiseB.noise_markings(n,1);
        last = noiseB.noise_markings(n,2);
        if (noiseB.noise_markings(n,5) == 'l')
            d1 = patch(a4, [beg last last beg],[a4.YLim(2) a4.YLim(2) a4.YLim(1) a4.YLim(1)],  'g', 'LineStyle', 'none');
            alpha(d1, .15);
        elseif (noiseB.noise_markings(n,5) == 'r')
            d1 = patch(a2, [beg last last beg], [a2.YLim(2) a2.YLim(2) a2.YLim(1) a2.YLim(1)], 'c', 'LineStyle', 'none');
            alpha(d1, .15);
        else
            d1 = patch(a4, [beg last last beg], [a4.YLim(2) a4.YLim(2) a4.YLim(1) a4.YLim(1)], 'r', 'LineStyle', 'none');
            alpha(d1, .15);
            d1 = patch(a2, [beg last last beg], [a2.YLim(2) a2.YLim(2) a2.YLim(1) a2.YLim(1)], 'r', 'LineStyle', 'none');
            alpha(d1, .15);
        end
    end
    
    
    beg = end_idx;
    end_idx = stepsize * x;
    
    key = input('press button');
    
    x = x+1
end