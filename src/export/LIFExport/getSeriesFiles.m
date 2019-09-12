% Extracts time information from xml files
function [XMLFolder, seriesPropertiesXML, seriesXML] = getSeriesFiles(Folder)
  XMLFolder = Folder;
  seriesPropertiesXML = dir([XMLFolder, filesep, '*Series*Properties.xml']);
  if isempty(seriesPropertiesXML)
      XMLFolder = [Folder, filesep, 'MetaData'];
      seriesPropertiesXML = dir([XMLFolder, filesep, '*Properties.xml']);
      if isempty(seriesPropertiesXML)
          error('XML MetaFiles could not be found. Did they get exported using the LAS software?')
      end
  end
  
  %Get the XML files in the metadata folder that do not contain the word "Properties"
  SeriesFilesTemp = dir([XMLFolder, filesep, '*Series*.xml']);
  seriesXML = [];
  for j = 1:length(SeriesFilesTemp)
    if mod(j,2) ~= 0
      seriesXML = [seriesXML,SeriesFilesTemp(j)];
    end
  end
end
