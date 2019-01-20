function copyDataToExpectedData(Prefix, dynamicsResultsPath, dynamicsResultsExpectedPath, preprocessedDataPath, preProcessedDataExpectedPath, processedDataPath, processedDataExpectedPath)
  disp(['Copying expected data for Prefix ', Prefix]);
  
  copyfile([dynamicsResultsPath, filesep, '*'], dynamicsResultsExpectedPath);
  disp(['Expected data copied to folder ', dynamicsResultsPath]);

  copyfile([preprocessedDataPath, filesep, '*'], preProcessedDataExpectedPath);
  disp(['Expected data copied to folder ', preProcessedDataExpectedPath]);
end
