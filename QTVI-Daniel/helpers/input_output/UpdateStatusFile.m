function UpdateStatusFile(location,analysisFiles)
    fid = fopen( fullfile(location,'pulse_analysis_status.csv'), 'wt' );
    fprintf(fid, 'Record,NoiseStatus,AnalysisStatus\n');
    for i = 1:size(analysisFiles,1)
      fprintf( fid, '%s,%i,%i\n', analysisFiles{i,1},analysisFiles{i,4},analysisFiles{i,5});
    end
    fclose(fid);
end

