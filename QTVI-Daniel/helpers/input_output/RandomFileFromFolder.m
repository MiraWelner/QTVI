
function [file_name,  file_folder, full_path] = RandomFileFromFolder(folder, extension)
    extension = strrep(extension,'.','');
    FileList = dir(fullfile(folder, ['**/*.' extension]));
    Index    = randi([1 numel(FileList)]);
    full_path = fullfile(FileList(Index).folder,FileList(Index).name);
    file_folder = FileList(Index).folder;
    file_name = FileList(Index).name;
end
