% Object to specify a preferred file to be used for testing experiments that have more
% than one file in the folder. This avoids asking for user prompt.
classdef PreferredFileForTest

  properties
    fileName;
  end

  methods 
    function obj = PreferredFileForTest(fileNameArg)
      obj.fileName = fileNameArg;
    end
  end

end
