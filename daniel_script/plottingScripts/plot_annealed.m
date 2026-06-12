close all;
cellidx = 1;
prev=0;
figure;
hold on;
for x =1:size(annealedSegments{1}.ppg_bin_indexs,1)
    ppgidxs = annealedSegments{1}.ppg_bin_indexs(x,1):annealedSegments{1}.ppg_bin_indexs(x,2);
%    ppgtime = ppgidxs/annealedSegments{1}.ppgSampleRate;
%    ppgidxs = (ppgidxs-annealedSegments{1}.ppg_bin_indexs(x,1))+1;
%    ecgidxs = annealedSegments{1}.ecg_bin_indexs(x,1):annealedSegments{1}.ecg_bin_indexs(x,2);
%    ecgtime = ecgidxs/annealedSegments{1}.ecgSampleRate;
%    ecgidxs = (ecgidxs-annealedSegments{1}.ecg_bin_indexs(x,1))+1;
   arr = (prev:prev+length(ppgidxs)-1)+1;

   plot(annealedSegments{1}.po(arr),'b');
   plot(annealedSegments{1}.ecg(arr),'r');
   prev = arr(end);

end
