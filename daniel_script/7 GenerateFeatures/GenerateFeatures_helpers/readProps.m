function [props] = readProps(inputFile)
%readProps read properties
%   Reads the properties from the inputFile and returns a struct
fid = fopen(inputFile);
entries  = [];
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if numel(line) == 0 || startsWith(line, '#')
        continue
    end
    [k, ~] = strsplit(line, '=');
    k = strtrim(k);
    entries{end + 1} = k;
end
props = containers.Map('KeyType', 'char', 'ValueType', 'char');
for i = 1:numel(entries)
    value = entries{i}{2};
    if isempty(entries{i}{2})
        value = '';
    end
    props(entries{i}{1}) = value;
end
fclose(fid);

