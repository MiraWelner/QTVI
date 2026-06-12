% 
% f = fields(beats_flattened_new);
% for x = 1:length(f)
%     tmp = beats_flattened_new.(f{x});
%    r = range(tmp);
%    m = nanmedian(tmp);
%    disp([f{x} ',' num2str(min(tmp)) ',' num2str(max(tmp)) ',' num2str(r) ',' num2str(m)]);
%     
% end




z=beats_flattened.msec_R_2_second_valley - beats_flattened.msec_R_2_first_valley;
z=(z*1000);
z+beats_flattened.msec_R_2_first_valley*-1;
fixed = z+beats_flattened.msec_R_2_first_valley*-1;



test = fixed(~isnan(fixed));
test2 = beats_flattened_new.msec_R_2_second_valley(~isnan(beats_flattened_new.msec_R_2_second_valley));




