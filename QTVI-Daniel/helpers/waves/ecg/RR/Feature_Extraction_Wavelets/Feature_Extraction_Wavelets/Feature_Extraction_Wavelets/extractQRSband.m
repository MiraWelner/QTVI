%Copyright 2017 The MathWorks, Inc.
function ecgQRSBand = extractQRSband(ecgsig,level,stLevl,endLevel) %#codegen
wt = modwt(ecgsig,level);
wtrec = zeros(size(wt));
wtrec(stLevl:endLevel,:) = wt(stLevl:endLevel,:);
ecgQRSBand = imodwt(wtrec,'sym4');