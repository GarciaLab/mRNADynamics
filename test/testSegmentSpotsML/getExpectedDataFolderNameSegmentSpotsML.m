function expectedDataFolder = getExpectedDataFolderNameSegmentSpotsML(testPath, step, subfolder, Prefix)
  % Adds underscore to the end of PreprocessedData subfolder
  if strcmpi('ProcessedData', subfolder)
    Prefix = [Prefix, '_'];
  end

  expectedDataFolder = [testPath, filesep, 'SegmentSpotsML', filesep, step, filesep, subfolder, filesep, Prefix];
end
