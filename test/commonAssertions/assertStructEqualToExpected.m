% Given the path of the dynamic results folder and the expected data path,
% asserts that two structs contained in a .mat file is equal to expected.
function assertStructEqualToExpected(testCase, containingFolder, expectedDataFolder, fileName)
  filePath = [containingFolder, filesep, fileName];
  testCase.assertTrue(exist(filePath, 'file') == 2);
  fileContents = load(filePath);

  expectedFilePath = [expectedDataFolder, filesep, fileName];
  expectedFileContents = load(expectedFilePath);

  fields = fieldnames(fileContents);
  fieldsInExpected = fieldnames(expectedFileContents);

  testCase.assertEqual(numel(fields), numel(fieldsInExpected));

  for i = 1:numel(fields)
    testCase.assertEqual(fileContents.(fields{i}), expectedFileContents.(fields{i}));
  end

end
