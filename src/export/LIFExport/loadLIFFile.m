%Loads LIF file and metadata
function [LIFImages, LIFMeta] = loadLIFFile(Folder)
  %Leica confocal
  LIFDir = dir([Folder,filesep,'*.lif']);

  %Load the file using BioFormats
   %Figure out which one is not the FF
  LIFIndex = find(~contains({LIFDir.name}, 'FF'));
  LIFPath = [Folder, filesep, LIFDir(LIFIndex).name];

  if isempty(LIFIndex)
      disp('couldn''t find movie file. try selecting it?');
      [file, path] = uigetfile(Folder);
      LIFPath = [path, file];
%       error('Only flat field LIF file found. No dataset LIF file found.\n%s', 'Check that your dataset LIF is present in the folder. Check that your full dataset name doesn''t contain the phrase ''FF''.')
  end

    LIFImages = bfopen(LIFPath);
    

  %Extract the metadata for each series
  LIFMeta = LIFImages{:, 4};
  
end