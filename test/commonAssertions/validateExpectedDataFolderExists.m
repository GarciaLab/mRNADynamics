function validateExpectedDataFolderExists(folder, Prefix)
  if 7 ~= exist(folder, 'dir')
    error(['Expected data for Prefix ', Prefix, ' does not exist. Please run generateExpectedDataForPrefix.'])
  end
end
