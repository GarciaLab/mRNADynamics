function exportAndCopyData(testCase, testPath, dynamicsResultsPath, preprocessedDataPath)
  experimentTestRootPath = strcat(testPath, filesep, 'ExportDataForLivemRNA', filesep);
  expectedPreProcFolder = strcat(experimentTestRootPath, 'PreProcessedData', filesep, testCase.Prefix);
  expectedDynamicsResultsFolder = strcat(experimentTestRootPath, 'DynamicsResults', filesep, testCase.Prefix);

  deleteDirectory(expectedPreProcFolder, testCase.Prefix);
  deleteDirectory(expectedDynamicsResultsFolder, testCase.Prefix);

  mkdir(expectedPreProcFolder);
  mkdir(expectedDynamicsResultsFolder);

  if (~ isprop(testCase, 'PreferredFileName'))
    ExportDataForLivemRNA(testCase.Prefix, 'keepTifs');
  else
    testCase.initializeTestCase;
    ExportDataForLivemRNA(testCase.Prefix, testCase.PreferredFileName, 'keepTifs');
  end

  disp(['Copying expected data for Prefix ', testCase.Prefix]);
  copyfile([dynamicsResultsPath, filesep, '*'], expectedDynamicsResultsFolder);
  disp(['Expected data copied to folder ', expectedDynamicsResultsFolder]);
  copyfile([preprocessedDataPath, filesep, '*'], expectedPreProcFolder);
  disp(['Expected data copied to folder ', expectedPreProcFolder]);
end
