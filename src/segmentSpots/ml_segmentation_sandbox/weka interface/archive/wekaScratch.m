%% Weka Example 
%% initialize weka

javaaddpath('C:\Program Files\Weka-3-8-4\weka.jar','-end')
arffloader = javaObject('weka.core.converters.ArffLoader') %this constructs an object of the arffloader class
%% load up data

arff = javaObject('java.io.File','C:\Program Files\Weka-3-8-4\data\iris.arff') %construct an arff file object
arffloader.setFile(arff)
data = arffloader.getDataSet()
data.toString
%% make a classifier

classifier = javaObject('weka.classifiers.trees.J48')
data.setClassIndex(data.numAttributes - 1)
classifier.buildClassifier(data)
classifier.toString
%% check out the classifier

evaluation = javaObject('weka.classifiers.Evaluation', data)
evaluation.crossValidateModel(classifier, data, 10,...
 javaObject('java.util.Random', 1),...
 javaArray('java.lang.Object',0))
evaluation.toSummaryString

%% save arff file into matlab readable matrix

saver = javaObject('weka.core.converters.MatlabSaver')
saver.setFile(javaObject('java.io.File', 'E:\Armando\LivemRNA\mRNADynamics\src\segmentSpots\ml_segmentation_sandbox\iris.data'))
saver.setInstances(data)
saver.writeBatch
m = load('E:\Armando\LivemRNA\mRNADynamics\src\segmentSpots\ml_segmentation_sandbox\iris.data')
scatter(m(:, 3), m(:, 4), 20, m(:, 5))
%% reduce data

f = javaObject('weka.filters.unsupervised.attribute.Remove')
f.setAttributeIndices('1-2')
f.setInputFormat(data)
rD = javaMethod('useFilter','weka.filters.Filter', data, f)
classifier.buildClassifier(rD)
classifier.toString
%% turn reduced data into a matrix 

rM = zeros(rD.numInstances, rD.numAttributes);
for i = 1:rD.numInstances
 for j = 1:rD.numAttributes
 rM(i,j) = rD.instance(i - 1).value(j - 1);
 end
end

%% Store predictions for reduced dataset in a matrix

p = zeros(rD.numInstances, rD.numClasses);
for i = 1:rD.numInstances
    dist = classifier.distributionForInstance(rD.instance(i - 1));
 for j = 1:rD.numClasses
    p(i,j) = dist(j);
 end
end

%% Plot data using colors based on predicted probabilities:

scatter(rM(:,1), rM(:,2), 20, p)
%% Generating predictions for a grid of points  

[x, y] = meshgrid(1:.1:7, 0:.1:2.5);
x = x(:);
y = y(:);
gM = [x y];
save('E:\Armando\LivemRNA\mRNADynamics\src\segmentSpots\ml_segmentation_sandbox\grid.data','gM','-ascii')
l = javaObject('weka.core.converters.MatlabLoader')
l.setFile(javaObject('java.io.File','E:\Armando\LivemRNA\mRNADynamics\src\segmentSpots\ml_segmentation_sandbox\grid.data'))
gD = l.getDataSet
gD.insertAttributeAt(rD.attribute(2), 2)
gD.setClassIndex(2)
p = zeros(gD.numInstances, gD.numClasses);
for i = 1:gD.numInstances
 dist = c.distributionForInstance(gD.instance(i - 1));
 for j = 1:gD.numClasses
 p(i, j) = dist(j);
 end
end
scatter(gM(:,1), gM(:,2), 20, p)
%% Clustering and visualizing data  

f = javaObject('weka.filters.unsupervised.! ! attribute.Remove')
f.setAttributeIndices('last')
f.setInputFormat(d)
rD = javaMethod('useFilter','weka.filters.Filter', d, f)
c = javaObject('weka.clusterers.SimpleKMeans')
c.setNumClusters(3)
c.buildClusterer(rD)
c.toString
a = zeros(rD.numInstances, 1);
for i = 1:rD.numInstances
 a(i) = c.clusterInstance(rD.instance(i - 1));
end
scatter(m(:,3), m(:,4), 20, a)

%% Finding the most important predictors  

as = javaObject('weka.attributeSelection.! ! AttributeSelection')
s = javaObject('weka.attributeSelection.GreedyStepwise')
s.setSearchBackwards(true)
as.setSearch(s)
e = javaObject('weka.attributeSelection.!!! WrapperSubsetEval')
e.setClassifier(javaObject('weka.!!! classifiers.trees.J48'))
as.setEvaluator(e)
as.SelectAttributes(d)
as.toResultsString
%% 
% 
%% Build a classifier with attribute selection  

c = javaObject('weka.classifiers.meta.!!!AttributeSelectedClassifier')
c.setEvaluator(e)
c.setSearch(s)
c.setClassifier(javaObject('weka.!!!!classifiers.trees.J48'))
c.buildClassifier(d)
c.toString

e = javaObject('weka.classifiers.Evaluation', d)
e.crossValidateModel(c, d, 10,...
 javaObject('java.util.Random', 1),...
 javaArray('java.lang.Object',0))
e.toSummaryString
%% 
%