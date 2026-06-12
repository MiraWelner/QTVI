function matched = import_matched_ids(filename)
opts = delimitedTextImportOptions("NumVariables", 63);
dataLines = [2, Inf];
% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["uid", "mesaid", "id", "mi", "rca", "ang", "ptca", "cbg", "othrev", "chf", "pvd", "strk", "tia", "dth", "chdh", "chda", "cvdh", "cvda", "revc", "chkdds", "cancer", "diab", "demen", "arsinf", "ourtd", "pneu", "copd", "asthma", "chawob", "orsysd", "hpfrct", "dvtpem", "onocvd", "mitt", "rcatt", "angtt", "ptcatt", "cbgtt", "othrevtt", "chftt", "pvdtt", "strktt", "tiatt", "dthtt", "chdhtt", "chdatt", "cvdhtt", "cvdatt", "revctt", "chkddstt", "cancertt", "diabtt", "dementt", "arsinftt", "ourtdtt", "pneutt", "copdtt", "asthmatt", "chawobtt", "orsysdtt", "hpfrcttt", "dvtpemtt", "onocvdtt"];
opts.VariableTypes = ["double", "double", "double", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["mi", "rca", "ang", "ptca", "cbg", "othrev", "chf", "pvd", "strk", "tia", "dth", "chdh", "chda", "cvdh", "cvda", "revc", "chkdds", "cancer", "diab", "demen", "arsinf", "ourtd", "pneu", "copd", "asthma", "chawob", "orsysd", "hpfrct", "dvtpem", "onocvd"], "EmptyFieldRule", "auto");

% Import the data
matched = readtable(filename, opts);



