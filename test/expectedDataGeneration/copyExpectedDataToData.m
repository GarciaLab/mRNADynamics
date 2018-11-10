function copyExpectedDataToData(Prefix, dynamicsResultsPath, dynamicsResultsExpectedPath, preprocessedDataPath, preProcessedDataExpectedPath, processedDataPath, processedDataExpectedPath)
  disp(['Copying expected data for Prefix ', Prefix]);
  
  if ~isempty(dynamicsResultsExpectedPath) 
    copyfile([dynamicsResultsExpectedPath, filesep, '*'], dynamicsResultsPath);
    disp(['Expected data copied to folder ', dynamicsResultsPath]);
  end

  if ~isempty(preProcessedDataExpectedPath) 
    copyfile([preProcessedDataExpectedPath, filesep, '*'], preprocessedDataPath);
    disp(['Expected data copied to folder ', preprocessedDataPath]);
  end

  if ~isempty(processedDataExpectedPath) 
    copyfile([processedDataExpectedPath, filesep, '*'], processedDataPath);
    disp(['Expected data copied to folder ', processedDataPath]);
  end

end
