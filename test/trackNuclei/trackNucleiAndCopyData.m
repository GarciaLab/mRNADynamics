function trackNucleiAndCopyData(Prefix, testPath, dynamicsResultsPath, preprocessedDataPath, processedDataPath)
  % Creates root folder if it does not exist
  if 7 ~= exist(testPath, 'dir')
    mkdir(testPath);
  end

  dynamicsResultsExpectedPath = createOrCleanTrackNucleiExpectedDataSubFolder(testPath, 'DynamicsResults', Prefix);
  preProcessedDataExpectedPath = createOrCleanTrackNucleiExpectedDataSubFolder(testPath, 'PreProcessedData', Prefix);
  processedDataExpectedPath = createOrCleanTrackNucleiExpectedDataSubFolder(testPath, 'ProcessedData', Prefix);
  
  TrackNuclei(Prefix);

  copyDataToExpectedData(Prefix, dynamicsResultsPath, dynamicsResultsExpectedPath, preprocessedDataPath,...
    preProcessedDataExpectedPath, processedDataPath, processedDataExpectedPath);
end

function expectedDataSubFolder = createOrCleanTrackNucleiExpectedDataSubFolder(testPath, subfolder, Prefix)
    % Adds underscore to the end of PreprocessedData subfolder
    if strcmpi('ProcessedData', subfolder)
      Prefix = [Prefix, '_'];
    end
  
  expectedDataSubFolder = [testPath, filesep, 'TrackNuclei', filesep, subfolder, filesep, Prefix];
  deleteDirectory(expectedDataSubFolder, Prefix);
  mkdir(expectedDataSubFolder);
end