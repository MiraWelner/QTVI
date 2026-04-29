function matched = import_matched_ids(filename)
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["idno", "mesaid"];
opts.VariableTypes = ["double", "double"];
%opts = setvaropts(opts, [3, 4, 5, 6, 7], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
matched = readtable(filename, opts);