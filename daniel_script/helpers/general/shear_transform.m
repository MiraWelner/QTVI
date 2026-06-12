function [transform] = shear_transform(time,slice)
    m= (slice(1)-slice(end))/(time(end)-time(1));
    c = -m*time(1);
    transform = zeros(length(slice),1);
    
    for x = 1:length(slice)
       transform(x) = slice(x) + m*x + c;
    end
%     sub = [];
%     for q = 1:length(slice)
%         sub(end+1) = m * q + 0;
% %     end
%     figure()
%     plot(slice);
%     hold on;
%     plot(transform);
%     z=1;
end