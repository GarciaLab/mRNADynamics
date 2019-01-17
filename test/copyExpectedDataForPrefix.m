function copyExpectedDataForPrefix(Prefix, step)
  %Figure out the initial folders.
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];

  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');

  dynamicsResultsPath = getConfigValue(configValues, 'DropboxFolder');
  preprocessedDataPath = getConfigValue(configValues, 'PreProcPath');
  processedDataPath = getConfigValue(configValues, 'FISHPath');
    
  testPath = getConfigValue(configValues, 'TestPath');

  if 7 ~= exist(testPath, 'dir')
    error('Expected data for does not exist. Please run generateExpectedDataForPrefix.')
  end
  
  preprocessedDataExperimentPath = [preprocessedDataPath, filesep, Prefix];
  dynamicsResultsExperimentPath = [dynamicsResultsPath, filesep, Prefix];
  processedDataExperimentPath = [processedDataPath, filesep, Prefix, '_'];

  % Clean up previous runs
  deleteDirectory(dynamicsResultsExperimentPath, Prefix);
  deleteDirectory(preprocessedDataExperimentPath, Prefix);
  deleteDirectory(processedDataExperimentPath, Prefix);

  if strcmpi(step, 'ExportDataForLivemRNA')
    copyExpectedExportData(Prefix, testPath, dynamicsResultsExperimentPath, preprocessedDataExperimentPath);
  elseif strcmpi(step, 'filterMovieTifs')
    copyFilterMovieTifsData(Prefix, testPath, dynamicsResultsExperimentPath, preprocessedDataExperimentPath);
  elseif strcmpi(step, 'filterMovieWeka')
    copyFilterMovieWekaData(Prefix, testPath, dynamicsResultsExperimentPath, preprocessedDataExperimentPath, processedDataExperimentPath);
  elseif strcmpi(step, 'SegmentSpots')
    copyExpectedSegmentSpotsData(Prefix, testPath, dynamicsResultsExperimentPath);
  elseif strcmpi(step, 'TrackNuclei')
    copyExpectedTrackNucleiData(Prefix, testPath, dynamicsResultsExperimentPath, preprocessedDataExperimentPath, [], []);
  elseif strcmpi(step, 'TrackmRNADynamics')
    copyExpectedTrackmRNADynamicsData(Prefix, testPath, dynamicsResultsExperimentPath, preprocessedDataExperimentPath, [], []);
  else
    error('Pipeline step not recognized, data cannot be copied.')
  end

end
