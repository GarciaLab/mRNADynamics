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


  
  %Load the data, this might cause problems with really large sets
  % Construct an empty Bio-Formats reader
    r = bfGetReader();
    % Decorate the reader with the Memoizer wrapper
    r = loci.formats.Memoizer(r);
   % Initialize the reader with an input file
    % If the call is longer than a minimal time, the initialized reader will
    % be cached in a file under the same directory as the initial file
    % name .large_file.bfmemo
    r.setId(LIFPath);
    LIFImages = bfopen(LIFPath);
    

  %Extract the metadata for each series
  LIFMeta = LIFImages{:, 4};
  
  r.close();
end