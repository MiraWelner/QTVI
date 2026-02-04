function LogError(script, location, function_inputs, error)
    fid = fopen(fullfile(location, [script '_errors.csv']), 'at');
    %fprintf(fid, 'Record,NoiseStatus,AnalysisStatus\n');
    tbl = table2cell(struct2table( error.stack));
    str = '';
    for i = 1:size(tbl,1)
        str = [str tbl{i,1} sprintf('\t') tbl{i,2} sprintf('\t') num2str(tbl{i,3}) newline];
    end
    str = ['"' str(1:end-1) '"'];
    
    
    input_string = [];
    fn = fieldnames(function_inputs);
    for i = 1:length(fn)
        input_string = [input_string function_inputs.(fn{i}) ','];
            
    end
    input_string = ['"' input_string '"'];
    
    fprintf( fid, '%s%s,%s,%s\n', error.identifier, ...
                                  error.message, ...
                                  str,...
                                  input_string);
    fclose(fid);
end