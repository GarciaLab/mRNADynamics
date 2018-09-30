% Object to specify a preferred classifier to be used for testing SegmentSpotsML without requiring user input.
classdef ClassifierForTest

  properties
    classifierFolder;
    classifierPathCh1;
  end

  methods 
    function obj = ClassifierForTest(classifierFolder, classifierPathCh1)
      obj.classifierFolder = classifierFolder;
      obj.classifierPathCh1 = classifierPathCh1;
    end
  end

end
