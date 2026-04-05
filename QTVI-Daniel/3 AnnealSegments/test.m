annealedSegments = load(analysisFiles{1,2});
annealedSegments = annealedSegments.annealedSegments;

for x = 1:length(annealedSegments)
    fprintf('Segment %d: PPG=%d, sleep=%d\n', x, length(annealedSegments{x}.po), length(annealedSegments{x}.sleep_stages));
end