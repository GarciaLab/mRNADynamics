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
%       pMap = classifyImage(im, arffFile)
%
%   See also GENERATEDOGSWEKA, GENERATEDOGS

%% PARAMS, OPTIONS

displayFigures = false;

if ischar(training), tempPath = fileparts(training);
else, tempPath = ''; end

classifierPath = '';
classifierObj = [];
reSc = false;
arffLoader = [];

for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

xDim=size(im, 2); yDim = size(im, 1); 
dim = length(size(im));
if dim==3
    zDim = size(im, 3);
end
ni = numel(im);

warning('off', 'MATLAB:Java:DuplicateClass');
warning('off', 'MATLAB:javaclasspath:jarAlreadySpecified');

%% initialize weka
javaaddpath('C:\Program Files\Weka-3-8-4\weka.jar','-end');

%% load up training data
if ischar(training)
    arffLoader = javaObject('weka.core.converters.ArffLoader'); %this constructs an object of  the arffloader class
    arffLoader.setFile(javaObject('java.io.File',training)); %construct an arff file object
    trainingData= arffLoader.getDataSet;
    trainingData.setClassIndex(trainingData.numAttributes - 1);
else
    trainingData = training;
end


%normalize data to the max of the training set for better classification
if reSc
    v = trainingData.attributeToDoubleArray(0);
    im = 5*im./max(max(im));     
end

%% load or build classifier from training data
if ~isempty(classifierObj)
    classifier = classifierObj;
elseif ~isempty(classifierPath) 
    classifier = weka.core.SerializationHelper.read(classifierPath);
else
    classifier = weka.classifiers.trees.RandomForest;
    options = {'-I', '64', '-K', '2', '-S', '-1650757608', '-depth', 20};
    javaaddpath('C:\Users\Armando\Desktop\fast random forest\fastrandomforest-2019.12.3.jar ')
   
%     classifier = javaObject('hr.irb.fastRandomForest.FastRandomForest');
%     options = {'-I', '20', '-threads', '1', '-K', '2', '-S', '-1650757608'};
%     
    classifier.setOptions(options);
    classifier.buildClassifier(trainingData);
    if displayFigures
        classifier.toString
    end
    
end

%% generate test data by filtering the image

% waitbarFigure = waitbar(0, 'filtering image');
% set(waitbarFigure, 'units', 'normalized', 'position', [0.4, .15, .25,.1]);

%create a new data set with the header copied from the training data. 

% testMatrix = zeros(ni, trainingData.numAttributes-1);

attributes = {};
nAtt = trainingData.numAttributes()-2;
for i = 0:nAtt+1
    attributes{i+1} = char(trainingData.attribute(i).name());
end

testData = arffLoader.getStructure;
for i = 1:ni 
    inst = weka.core.DenseInstance(nAtt+1);
    testData.add(inst);   
end

name = 'segment'; %hardcoded. change later. 


for i = 0:nAtt-1 %don't do the class label
    attributes{i+1} = char(trainingData.attribute(i).name());
    att = attributes{i+1};
    %might error with some attribtes. it'd be nice to write a filter that
    %removes attributes that aren't supported in matlab so the models and
    %classifiers will still work 
    
    disp(['Generating feature ', num2str((i+1)), '/', num2str(nAtt+1) , ': ', att]);
    
    if i > 0
        filterType = regexp(att, '.*(?=_\d)', 'match');
        sigmas = regexp(att, '(\d[.]\d)|(\d\d[.]\d)', 'match');
        filteredIm = filterImage(im, filterType{1}, sigmas);
    else
        filteredIm = im;
   end
    
%     testDataTemp = mat2instances(attributes, name, filteredIm(:));
       filteredIm = filteredIm(:);
     

       for k = 1:numel(filteredIm)
           
            testData = setInstancesVal(testData, k-1, i, filteredIm(k));
            
       end
        
    %     testMatrix(:,i+1) = filteredIm(:);
    
%     testData.add(filteredIm(:));
    
%      waitbar(i/nAtt, waitbarFigure);
    
end

% close(waitbarFigure);

%now turn data matrix into weka arff
% tempFile = [tempPath, filesep, 'temp.data'];
% save(tempFile,var2str(testMatrix),'-ascii');
% loader = javaObject("weka.core.converters.MatlabLoader");
% loader.setFile(javaObject('java.io.File',tempFile));
% testData = loader.getDataSet;
% testData.insertAttributeAt(trainingData.classAttribute,trainingData.classIndex);
testData.setClassIndex(testData.numAttributes-1);

% f = javaObject('weka.filters.unsupervised.attribute.Remove');
% f.setInputFormat(trainingData);
% testData = javaMethod('useFilter','weka.filters.Filter', testData, f); %match test header to training header

nInstances = testData.numInstances;

compatible = testData.equalHeaders(trainingData); %check compatability between arff and data
if ~compatible
    error(char(testData.equalHeadersMsg(trainingData)));
end

%this step is unbearably slow. ~1min per frame. needs a massive overhaul
pLin = zeros(testData.numInstances, 1);

% testData_constant = parallel.pool.Constant(testData);
% classifier_constant = parallel.pool.Constant(classifier);

for i = 1:nInstances
    
    pLin(i) = classifyInstance(i, testData, classifier);
    
end

if dim == 2
    pMap = reshape(pLin, [yDim xDim]);
elseif dim == 3
     pMap = reshape(pLin, [yDim xDim zDim]);
end

if displayFigures & dim==2
    figure(1); imshowpair(im, pMap, 'montage');
end

end

 function p = classifyInstance(ind, testData, classifier)
         k = ind - 1;
        inst = testData.instance(k);
        temp = classifier.distributionForInstance(inst);
        p = temp(1);
 end