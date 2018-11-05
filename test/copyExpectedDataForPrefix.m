function copyExpectedDataForPrefix(Prefix, step)
  %Figure out the initial folders.
  CONFIG_CSV_PATH = ['ComputerFolders.csv'];

  configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');

  dynamicsResultsPath = getConfigValue(configValues, 'DropboxFolder');
  preprocessedDataPath = getConfigValue(configValues, 'PreProcPath');
  processedDataPath = getConfigValue(configValues, 'FISHPath');
    
  codePath = getConfigValue(configValues, 'MS2CodePath');
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

  if strcmpi(step, 'ExportDataForFISH')
    copyExpectedExportData(Prefix, testPath, dynamicsResultsExperimentPath, preprocessedDataExperimentPath);
  elseif strcmpi(step, 'SegmentSpotsML')
    copyExpectedSegmentSpotsMLData(Prefix, testPath, dynamicsResultsExperimentPath, preprocessedDataExperimentPath,...
      processedDataExperimentPath);
  elseif strcmpi(step, 'TrackNuclei')
    copyExpectedTrackNucleiData(Prefix, testPath, dynamicsResultsExperimentPath, preprocessedDataExperimentPath,...
      processedDataExperimentPath);
  elseif strcmpi(step, 'TrackmRNADynamics')
    copyExpectedTrackmRNADynamicsData(Prefix, testPath, dynamicsResultsExperimentPath, preprocessedDataExperimentPath,...
      processedDataExperimentPath);
  else
    error('Pipeline step not recognized, data cannot be copied.')
  end

end
