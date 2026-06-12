% [template t2 valid] = PPG_SQI(wave,anntime,annot,template_ahead,windowlen,Fs)
%
% PPG_SQI.m - PPG SQI based on beat template correlation.
% (as an advice, the algorithm get 30 beats at each call and run in loop)
% by Qiao Li 30 Mar 2011
%
% PPG sampling frequency is Fs
%
% input:
%     wave:       PPG data;
%     anntime:    PPG annotation time (samples), read from ple annot file,
%                 But the ann-time is the OFFSET based on wave(1)
%     annot:      Annotation of beats, read from from ple annot file
%                 directly
%     template:   Last PPG beat template
%     windowlen:  length of window to calculate template(default: 30s)
%     Fs       :  sampling frequency (defatult: 125 to work with pervious code)
% output:
%     annot:      ppg sqi annotation
%                     annot.typeMnemonic: E - excellent beat;
%                                         A - acceptable beat;
%                                         Q - unacceptable beat
%                     annot.subtype: SQI based on Direct compare
%                     annot.chan:    SQI based on Linear resampling
%                     annot.num:     SQI based on Dynamic time warping
%                     annot.auxInfo: SQI based on Clipping detection
%     template:   Current PPG beat template
%     valid:      1 or greater for valid template,
%                 0 for invalid template
%
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information
%
% 12-01-2107 Modified by Giulia Da Poian: replace fixed samplig frequency
% (125) with a variable Fs

function [typeMnemonic, labels, template, valid] = PPG_SQI(wave, anntime, template, windowlen, Fs)
    warning('off','all');

    typeMnemonic = cell(length(anntime) - 1,1);
    if nargin < 5
        Fs = 256;
    end

    if nargin < 4
        windowlen = 30 * Fs;
    end

    if nargin < 3
        template = [];
    end

    if nargin < 2
        sprintf('Error: must provide wave, anntime');
        template = [];
        valid = 0;
        return;
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 19/09/2011 ADD baseline wander filter
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % get PPG template
%     if length(template) > 1
        v = 1;
        t2 = template;
        t = template;
%     else
%         [t t2 v] = template_pleth(wave(1:min(windowlen, length(wave))), anntime(anntime < min(windowlen, length(wave))), 0, Fs);
%     end

    if length(template) < 1 % && v < 1% Current template invalid && no previous template available
        template = t;
        valid = 0;
        typeMnemonic = cell2mat(typeMnemonic);
        if isempty(typeMnemonic)
            typeMnemonic = [0,0,0,0,0];
        end
        labels = {'mean_corr_dtw', 'corrcoff_direct', 'corrcoff_interp', 'dtw', 'frechet'};
        warning('on','all');
        return
%        error('shouldn''t happen');
% 
%         t = [];
    else
        wave = PPGmedianfilter(wave, Fs, Fs);
        % Using previous template as the template
        if v < 1
            t = template;
        end

        % Using t2 if available
        if v > 1
            t = t2;
        end

        % Calculate the PLA (Piecewise linear approximation (PLA)) of template for dynamic time warping
        d1 = t;
        d1 = (d1 - min(d1)) / (max(d1) - min(d1)) .* 100;
        [tmp, pla1] = PLA(d1, 1, 1);
%         plot(pla1,y1(pla1));
        % Main Loop
        parfor j = 1:length(anntime) - 1
            %             SQI1: Direct compare

            %             Calculate correlation coefficients based on the template
            %             length
            try
                y1=tmp;
                beatbegin = anntime(j);
                beatend = anntime(j + 1);
                % 07/11/2011 ADD max beat length <= 3s detection
                if beatend - beatbegin > 3 * Fs
                    beatend = beatbegin + 3 * Fs;
                end

                templatelength = length(t);
                complength = min(templatelength, beatend - beatbegin - 1);
                % check within bounds of wave for smallest length (either
                % template or beat)
                if beatbegin + complength - 1 > length(wave) || beatend > length(wave) || beatbegin < 1
                    typeMnemonic{j} = [0,0,0,0,0];
                    continue;
                end

                cc = corrcoef(t(1:complength), wave(beatbegin:beatbegin + complength - 1));
                c1 = cc(1, 2);


                if (c1 < 0)
                    c1 = 0;
                end
                
                %annot(currentb).subtype = int8(c1 * 100);
%                 subtype = int8(c1 * 100);
                %             SQI2: Linear resampling

                %             Calculate correlation coefficients based on the linear
                %             resampling (interp1)

                y = interp1(1:beatend - beatbegin, wave(beatbegin:beatend - 1), 1:(beatend - beatbegin - 1) / (templatelength - 1):(beatend - beatbegin), 'spline');
                y(isnan(y)) = 0;
                cc = corrcoef(t, y);
                c2 = cc(1, 2);

                if (c2 < 0)
                    c2 = 0;
                end

%                 chan = int8(c2 * 100);

                %             SQI3: Dynamic Time Warping

                %             Calculate correlation coefficients based on the dynamic time
                %             warping

                d2 = wave(beatbegin:beatend - 1);

                % if beat too long, set SQI=0;
                if (length(d2) > length(d1) * 10)
                    c3 = 0;
                else
                    % normalize wave
                    d2 = (d2 - min(d2)) / (max(d2) - min(d2)) .* 100;
                    [y2, pla2] = PLA(d2, 1, 1);
%                     plot(1:length(y2),y2);
%                     hold on;
%                      plot(pla1,y1(pla1),'.-');
%                      plot(pla2,y2(pla2),'.-');
                    [w, ta, tb] = simmx_dtw(y1, pla1, y2, pla2); % diff of pla slopes
                    [p, q, ~] = dp_dtw2(w, ta, tb); % dynamic programming shortest DTW path
                    [~, ym2 , ~] = draw_dtw(y1, pla1, p, y2, pla2, q); % create line representing difference
                    cc = corrcoef(y1, ym2); % calculate error btw template and line 
                    c3 = cc(1, 2);

                    if (c3 < 0)
                        c3 = 0;
                    end

                end
%                 
%                 plot(t);
%                 hold on;
%                 plot(y);
%                 legend('Template','Beat compare')
%                 text(.15,.15,['Corr Coff direct to shortest: ' num2str(c1)],'Units','normalized')
%                 text(.1,.1,['Corr Coff interp: ' num2str(c2)],'Units','normalized')
%                 text(.1,.1,['Corr Coff interp: ' num2str(c2)],'Units','normalized')
%                 text(.05,.05,['Corr Coff interp: ' num2str(c3)],'Units','normalized')
%                 num = int8(c3 * 100);

                %             SQI4: Clipping detection

%                 d2 = wave(beatbegin:beatend - 1);
%                 y = diff(d2);
%                 clipthreshold = 0.5;
%                 c4 = length(find(abs(y) > clipthreshold)) / length(y);

                %             SQI: Combined
                
                templatelength = length(t);
                                
                if ~iscolumn(y1)
                    y1=y1';
                end 
                if ~iscolumn(y2)
                    y2=y2';
                end 

                beatlen = beatend - beatbegin;
                if templatelength < beatlen
                    x = linspace(0,1,templatelength);
                    x_finer = linspace(0,1,beatlen);
                    z = interp1(x, y1/100, x_finer, 'spline');        
                    if ~iscolumn(z)
                        z=z';
                    end 
                    [cm, cSq] =DiscreteFrechetDist(z,y2/100);
                elseif templatelength > beatlen
                    x = linspace(0,1,beatlen);
                    x_finer = linspace(0,1,templatelength);
                    z = interp1(x, y2/100, x_finer, 'spline');

                    if ~iscolumn(z)
                        z=z';
                    end 
                    [cm, cSq] =DiscreteFrechetDist(y1/100,z);
                else
                    [cm, cSq] = DiscreteFrechetDist(y1/100, y2/100);
                end
                
                

                dif = sqrt(2)-(cm/2);
                Frechet=((dif * 1) / sqrt(2));
                typeMnemonic{j} = [mean([c1 c2 c3]), c1, c2, c3, Frechet];
                
% 
%                 if min(sqibuf) >= 90
%                     typeMnemonic(j) = 'A'; % Excellent %65
%                 else
% 
%                     if (median(sqibuf(1:3)) >= 80 && sqibuf(1) >= 50 && c4 < 0.3) || (min(sqibuf) >= 70 && c4 < 0.3)
%                         typeMnemonic(j) = 'B'; % Acceptable %66
%                     elseif (median(sqibuf(1:3)) >= 70 && sqibuf(1) >= 50)
%                         typeMnemonic(j) = 'C'; % bad % 67
%                     elseif (median(sqibuf(1:3)) > 60)
%                         typeMnemonic(j) = 'D'; % poor %68
%                     else
%                         typeMnemonic(j) = 'Z'; % Unacceptable %90
%                     end
% 
%                 end

            catch
                typeMnemonic{j} = [0,0,0,0,0];
            end

        end

    end

    template = t;
    valid = v;
    typeMnemonic = cell2mat(typeMnemonic);
    labels = {'mean_corr_dtw', 'corrcoff_direct', 'corrcoff_interp', 'dtw', 'frechet'};
    warning('on','all');
end
