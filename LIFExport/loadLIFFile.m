%Loads LIF file and metadata
function [LIFDir, LIFIndex, LIFImages, LIFMeta] = loadLIFFile(Folder)
  %Leica confocal
  LIFDir = dir([Folder,filesep,'*.lif']);

  %Load the file using BioFormats
  %Figure out which one is not the FF
  LIFIndex = find(cellfun(@isempty, strfind({LIFDir.name}, 'FF')));

  %Load the data, this might cause problems with really large sets
  LIFImages = bfopen([Folder, filesep, LIFDir(LIFIndex).name]);

  %Extract the metadata for each series
  LIFMeta = LIFImages{:, 4};
end
