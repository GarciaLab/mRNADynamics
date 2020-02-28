function pMap = classifyImage(im, training, varargin)
%classifyImage Converts image to probability map.
%   pmap = classifyImage(im, training, varargin) Converts image to probability map.
%
%
%       IM        m x n numeric matrix
%       TRAINING training is the path to a training data set-
%                       e.g. 'E:\classifiers backup 9-20-19\classifiers\nuclear\'70nmgood.arff'
%
%
%   Notes:
%
%
%   Examples:
%
%       % example histone image
%       load('E:\DorsalSynthetics\Data\PreProcessedData\2020-01-21-1Dg-8D_EfEfEf_9\2020-01-21-1Dg-8D_EfEfEf_9_hisMat.Mat')
%       im = double(squeeze(hisMat(1, :, :)));
%       arffFile =  'E:\classifiers backup 9-20-19\classifiers\nuclear\70nmgood.arff';
%       pMap = classifyImage(im, arfffile)
%
%   See also GENERATEDOGSWEKA, GENERATEDOGS

%% PARAMS, OPTIONS

displayFigures = false;
tempPath = fileparts(training);
classifierPath = '';
reSc = false;

for i = 1:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

xDim=size(im, 2); yDim = size(im, 1);
ni = xDim*yDim;

%% helper functions
paren = @(x, varargin) x(varargin{:}); %emulates the syntax array(i)(j)

%% initialize weka
javaaddpath('C:\Program Files\Weka-3-8-4\weka.jar','-end');

%% load up training data
arffloader = javaObject('weka.core.converters.ArffLoader'); %this constructs an object of the arffloader class
arffloader.setFile(javaObject('java.io.File',training)); %construct an arff file object
trainingData= arffloader.getDataSet;

%normalize data to the max of the training set for better classification
if reSc
    v = trainingData.attributeToDoubleArray(0);
    im = (range(v)/median(v))* im./max(max(im));      
end

%% load or build classifier from training data
if ~isempty(classifierPath)
    classifier = weka.core.SerializationHelper.read(classifierPath);
else
    classifier = weka.classifiers.trees.RandomForest;
    trainingData.setClassIndex(trainingData.numAttributes - 1);
    options = {'-I' '200' '-K' '2' '-S' '-1650757608'};
    classifier.setOptions(options);
    classifier.buildClassifier(trainingData);
    
    if displayFigures
        classifier.toString
    end
    
end

%% generate test data by filtering the image

waitbarFigure = waitbar(0, 'filtering image');

testMatrix = zeros(ni, trainingData.numAttributes-1);
attributes = {};
nAtt = trainingData.numAttributes()-2;
for i = 0:nAtt
    attributes{i+1} = char(trainingData.attribute(i).name());
    attribute = trainingData.attribute(i);
    %might error with attributes with 2 sigmas. need to fix -AR
    if i
        filterType = regexp(attributes{i+1}, '.*(?=_)', 'match');
        sigma =  flip(num2str(sscanf(flip(attributes{i+1}), '%f')));
        fim = filterImage(im, filterType{1}, {sigma});
    else
        fim = im;
    end
    
    testMatrix(:,i+1) = fim(:);
    
     waitbar(i/nAtt, waitbarFigure);
    
end

close(waitbarFigure);

%now turn data matrix into weka arff
tempFile = [tempPath, filesep, 'temp.data'];
save(tempFile,var2str(testMatrix),'-ascii');
loader = javaObject("weka.core.converters.MatlabLoader");
loader.setFile(javaObject('java.io.File',tempFile));
testData = loader.getDataSet;
testData.insertAttributeAt(trainingData.classAttribute,trainingData.classIndex)
testData.setClassIndex(testData.numAttributes-1)

f = javaObject('weka.filters.unsupervised.attribute.Remove');
f.setInputFormat(trainingData)
testData = javaMethod('useFilter','weka.filters.Filter', testData, f); %match test header to training header

nInstances = testData.numInstances;

compatible = testData.equalHeaders(trainingData); %check compatability between arff and data
if ~compatible
    error(char(testData.equalHeadersMsg(trainingData)));
end

%this step is unbearably slow. >1min per frame. needs a massive overhaul
pLin = zeros(testData.numInstances, 1);
waitbarFigure = waitbar(0, 'classifying image');
set(waitbarFigure, 'units', 'normalized', 'position', [0.4, .15, .25,.1]);
for i = 1:nInstances
    
    pLin(i) = paren(classifier.distributionForInstance(testData.instance(i - 1)), 1);
    %     pLin(i) = paren(c.classifyInstance(rD.instance(i - 1)), 1);
    waitbar(i/nInstances, waitbarFigure);
    
end

close(waitbarFigure);

pMap = reshape(pLin, [yDim xDim]);

if displayFigures
    figure(1); imshowpair(im, pMap, 'montage');
end

end