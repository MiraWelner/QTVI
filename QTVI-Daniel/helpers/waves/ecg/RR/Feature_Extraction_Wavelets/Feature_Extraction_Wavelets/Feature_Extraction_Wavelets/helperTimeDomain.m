%Copyright 2017 The MathWorks, Inc.
function helperTimeDomain(t, sig,titleStr,plotDur, color)
if isempty(plotDur)
    plotDurSamp = numel(sig);
else
    plotDurSamp= plotDur*1000;
end

plot(t(1:plotDurSamp), sig(1:plotDurSamp)); axis tight
a = findall(gca,'Type','line');
switch color
    case('r')
        cvec = [1,0,0];
    case('b')
        cvec = [0,0,1];
    case('gray')
        cvec = [0.32,0.55,0.55];
 
    otherwise
        cvec = [0.5,0.5,0.5];
end
set(a,'Color',cvec);


title(titleStr);
grid on; box on;
xlabel('time (sec)');
ylabel('Amplitude');

