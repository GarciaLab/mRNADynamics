% Extracts time information from xml files
function [XMLFolder, SeriesFiles, SeriesFiles3] = getSeriesFiles(Folder)
  XMLFolder = Folder;
  SeriesFiles = dir([XMLFolder, filesep, '*Series*Properties.xml']);
  if isempty(SeriesFiles)
      XMLFolder = [Folder, filesep, 'MetaData'];
      SeriesFiles = dir([XMLFolder, filesep, '*Series*Properties.xml']);
      if isempty(SeriesFiles)
          error('XML MetaFiles could not be found. Did they get exported using the LAS software?')
      end
  end
  
  %Get the XML files in the metadata folder that do not contain the word "Properties"
  SeriesFilesTemp = dir([XMLFolder, filesep, '*Series*.xml']);
  SeriesFiles3 = [];
  for j = 1:length(SeriesFilesTemp)
    if mod(j,2) ~= 0
      SeriesFiles3 = [SeriesFiles3,SeriesFilesTemp(j)];
    end
  end
end
