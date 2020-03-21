function renameAllDogs(Prefix)


 [~, ProcPath] = DetermineLocalFolders(Prefix);

% Get all files in the current folder
dogFolder = [ProcPath, filesep, Prefix, '_', filesep, 'dogs', filesep];
files = dir([dogFolder, filesep, '*.tif';]);
% Loop through each
for id = 1:length(files)
    % Get the file name (minus the extension)
%     [a, f] = fileparts(files(id).name);
    oldName = files(id).name;
%     newName = ['DOG_2017-09-14-P2P-MCP-NoNLS-mCherry-doubledosage2_rerun_001_z01_ch01'];
    newName = ['prob', Prefix, '_', oldName(1:4), oldName(end-7:end)];
  
    %this renames files despite the name "movefile"'
    movefile([dogFolder, filesep, oldName], [dogFolder, filesep, newName]);
end


end