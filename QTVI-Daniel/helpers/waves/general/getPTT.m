function [PTT] = getPTT(PPGfeature, ECGfeature, PPGtime, ECGtime)

    try
        PPGevent = PPGtime(PPGfeature);
        ECGevent = ECGtime(ECGfeature);
        [ECGevent2, PPGevent2] = fixRelate(ECGevent, PPGevent);
        PTT = PPGevent2 - ECGevent2;
    catch
        disp('oh fugg :DDDD');
    end

end
