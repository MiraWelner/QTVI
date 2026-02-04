function heartCandles(time,high,low,q1,q3,med,slope)
    width_fill=.95; %width of body as % of space btw time, change to draw body thicker or thinner
    
    % comment
    hold on
    timediff = nanmedian(diff(time));
    centers = timediff/2;
    centers = time+centers;
    
    width = (width_fill*timediff)/2;
    
    for i=1:length(time)
        if ~isnan(low(i))
            line([centers(i) centers(i)],[low(i) high(i)], 'color', 'k');

            if (sign(slope(i)) == 0) || (sign(slope(i)) == 1)
                x=[centers(i)-width centers(i)+width centers(i)+width centers(i)-width];
                y=[q1(i) q1(i) q3(i) q3(i)];
                patch(x,y,'g')
            else
                x=[centers(i)-width centers(i)+width centers(i)+width centers(i)-width];
                y=[q1(i) q1(i) q3(i) q3(i)];
                patch(x,y,'r')
            end
        end
    end 
    hold off
    
    % comment
end