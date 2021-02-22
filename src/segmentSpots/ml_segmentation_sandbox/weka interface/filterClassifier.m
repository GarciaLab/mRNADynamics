function filterClassifier(trainingFile, classifierFile)

maxDepth = 0;
nTrees = 200;

% trainingFile = 'S:\Armando\Dropbox\DorsalSyntheticsDropbox\training_data_and_classifiers\mRNA\veryspecialsimonclassifier_almost_3_filteredRanker.arff';

arffLoader = javaObject('weka.core.converters.ArffLoader'); %this constructs an object of the arffloader class
arffLoader.setFile(javaObject('java.io.File',trainingFile)); %construct an arff file object
trainingData= arffLoader.getDataSet;
trainingData.setClassIndex(trainingData.numAttributes - 1);
  

resampleFilter = weka.filters.supervised.instance.Resample;
resampleFilter.setOptions({'-B', '1.0','-S', '1', '-Z','100.0'});
resampleFilter.setInputFormat(trainingData);
trainingData = weka.filters.Filter.useFilter(trainingData, resampleFilter);

classifier =hr.irb.fastRandomForest.FastRandomForest;
options = {'-I', num2str(nTrees), '-threads', 12, '-K', '2', '-S',...
    '-1650757608', '-depth', num2str(maxDepth)};


classifier.setOptions(options);
classifier.buildClassifier(trainingData);
weka.core.SerializationHelper.write(classifierFile, classifier);


% writePath = 'E:\Armando\LivemRNA\classifier.model'; 
% readPath = 'S:\Armando\Dropbox\DorsalSyntheticsDropbox\training_data_and_classifiers\mRNA\veryspecialsimonclassifier_almost_3.model';
% classifier = weka.core.SerializationHelper.read(readPath);


