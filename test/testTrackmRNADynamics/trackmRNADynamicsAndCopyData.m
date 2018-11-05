function trackmRNADynamicsAndCopyData(testCase, testPath, dynamicsResultsPath, preprocessedDataPath, processedDataPath)
  % Creates root folder if it does not exist
  if 7 ~= exist(testPath, 'dir')
    mkdir(testPath);
  end

  dynamicsResultsExpectedPath = createOrCleanTrackmRNADynamicsExpectedDataSubFolder(testPath, 'DynamicsResults', testCase.Prefix);
  preProcessedDataExpectedPath = createOrCleanTrackmRNADynamicsExpectedDataSubFolder(testPath, 'PreProcessedData', testCase.Prefix);
  processedDataExpectedPath = createOrCleanTrackmRNADynamicsExpectedDataSubFolder(testPath, 'ProcessedData', testCase.Prefix);
  
  TrackmRNADynamics(testCase.Prefix, testCase.Threshold1, testCase.Threshold2);

  copyDataToExpectedData(testCase.Prefix, dynamicsResultsPath, dynamicsResultsExpectedPath, preprocessedDataPath,...
    preProcessedDataExpectedPath, processedDataPath, processedDataExpectedPath);

  % Runs the tracking a second time and stores the results of it, so we can test the case of Particles.mat already existing.
  TrackmRNADynamics(testCase.Prefix, testCase.Threshold1, testCase.Threshold2);

  secondPassSubFolder = [testPath, filesep, 'TrackmRNADynamics', filesep, 'secondPass', filesep, 'DynamicsResults', filesep, testCase.Prefix];
  deleteDirectory(secondPassSubFolder, testCase.Prefix);
  mkdir(secondPassSubFolder);

  copyfile([dynamicsResultsPath, filesep, '*'], secondPassSubFolder);
  disp(['Expected data copied to folder ', secondPassSubFolder]);
end

function expectedDataSubFolder = createOrCleanTrackmRNADynamicsExpectedDataSubFolder(testPath, subfolder, Prefix)
    % Adds underscore to the end of PreprocessedData subfolder
    if strcmpi('ProcessedData', subfolder)
      Prefix = [Prefix, '_'];
    end
  
  expectedDataSubFolder = [testPath, filesep, 'TrackmRNADynamics', filesep, subfolder, filesep, Prefix];
  deleteDirectory(expectedDataSubFolder, Prefix);
  mkdir(expectedDataSubFolder);
end