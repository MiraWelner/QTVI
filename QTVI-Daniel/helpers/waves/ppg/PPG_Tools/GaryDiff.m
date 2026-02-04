function sqibuf = garyDiff(t,wave)
    d1 = t;
    d1 = (d1 - min(d1)) / (max(d1) - min(d1)) .* 100;
    [y1 pla1] = PLA(d1, 1, 1);

    complength = min(length(t), length(wave));

    cc = corrcoef(t(1:complength), wave(1:complength));
    c1= cc(1, 2);
    if (c1 < 0)
        c1 = 0;
    end
    subtype = int8(c1 * 100);
    %             SQI2: Linear resampling

    %             Calculate correlation coefficients based on the linear
    %             resampling (interp1)

    y = interp1(1:length(wave), wave, 1:(length(wave)-1) / (length(t) - 1):length(wave), 'spline');
    y(isnan(y)) = 0;
    cc = corrcoef(t, y);
    c2 = cc(1, 2);
    if (c2 < 0)
        c2 = 0;
    end
    chan = int8(c2 * 100);

    %             SQI3: Dynamic Time Warping

    %             Calculate correlation coefficients based on the dynamic time
    %             warping

    d2 = wave;

    % if beat too long, set SQI=0;
    if (length(d2) > length(d1) * 10)
        c3 = 0;
    else
        d2 = (d2 - min(d2)) / (max(d2) - min(d2)) .* 100;
        [y2 pla2] = PLA(d2, 1, 1);

        [w ta tb] = simmx_dtw(y1, pla1, y2, pla2);
        [p, q, Dm] = dp_dtw2(w, ta, tb);
        [ym1 ym2 yout1] = draw_dtw(y1, pla1, p, y2, pla2, q);
        cc = corrcoef(y1, ym2);
        c3 = cc(1, 2);
        if (c3 < 0)
            c3 = 0;
        end
    end
    num = int8(c3 * 100);



    %             SQI: Combined

    sqibuf = [subtype chan num];
end