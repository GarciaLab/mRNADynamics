function copyExpectedExportData(Prefix, testPath, dynamicResultsPath, preprocessedDataPath)
  disp(['Copying expected data for ExportDataForFISH test case with Prefix', Prefix]);

  experimentTestRootPath = strcat(testPath, filesep, 'ExportDataForFISH', filesep);

  expectedDynamicsResultsFolder = strcat(experimentTestRootPath, 'DynamicsResults', filesep, Prefix);
  validateExpectedDataFolderExists(expectedDynamicsResultsFolder, Prefix)

  expectedPreProcFolder = strcat(experimentTestRootPath, 'PreProcessedData', filesep, Prefix);
  validateExpectedDataFolderExists(expectedPreProcFolder, Prefix)

  deleteDirectory(dynamicResultsPath, Prefix);
  mkdir(dynamicResultsPath);

  deleteDirectory(preprocessedDataPath, Prefix);
  mkdir(preprocessedDataPath);

  disp(['Copying expected data for Prefix ', Prefix]);
  copyfile([expectedDynamicsResultsFolder, filesep, '*'], dynamicResultsPath);
  disp(['Expected data copied to folder ', dynamicResultsPath]);
  copyfile([expectedPreProcFolder, filesep, '*'], preprocessedDataPath);
  disp(['Expected data copied to folder ', preprocessedDataPath]);
end
