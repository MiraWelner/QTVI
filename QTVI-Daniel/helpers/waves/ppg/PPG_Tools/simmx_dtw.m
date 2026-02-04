% function [w ta tb] = simmx_dtw(y1,pla1,y2,pla2)
%
% calculate a sim matrix between y1 and y2.
%

function [w, ta, tb] = simmx_dtw(y1,pla1,y2,pla2)

slope1=zeros(1,length(pla1)); % slope of peicewize line 1
ta=zeros(1,length(pla1)); % time axis of peicewize line 1
ta(1)=1;
for i=2:length(pla1)
    slope1(i)=(y1(pla1(i))-y1(pla1(i-1)))/(pla1(i)-pla1(i-1));
    ta(i)=pla1(i)-pla1(i-1);
end

slope2=zeros(1,length(pla2)); % slope of peicewize line 2
tb=zeros(1,length(pla2)); % time axis of peicewize line 2
tb(1)=1;
for i=2:length(pla2)
    slope2(i)=(y2(pla2(i))-y2(pla2(i-1)))/(pla2(i)-pla2(i-1));
    tb(i)=pla2(i)-pla2(i-1);
end

A1 = zeros(length(slope1),length(slope2));
for i=1:length(slope2)
    A1(:,i)=slope1;
end
B1 = zeros(length(slope1),length(slope2));
for j=1:length(slope1)
    B1(j,:)=slope2';
end

w = abs(A1-B1);