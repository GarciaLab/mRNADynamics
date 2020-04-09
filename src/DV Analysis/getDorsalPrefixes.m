function [dataFolder, resultsFolder, prefixes] = getDorsalPrefixes(DataType)

[dataFolder, resultsFolder] = getDorsalFolders;

if exist([resultsFolder,filesep,DataType,filesep,'prefixes.mat'], 'file')
    load([resultsFolder,filesep,DataType,filesep,'prefixes.mat'], 'prefixes');
else
    [~, prefixes, ~] = LoadMS2Sets(DataType, 'justPrefixes', 'noCompiledNuclei');
end

end