function pMap = classifyImageMatlab(im, training, varargin)
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

numInstances = numel(im);

warning('off', 'MATLAB:Java:DuplicateClass');
warning('off', 'MATLAB:javaclasspath:jarAlreadySpecified');

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
    im = 5*im./max(max(im));
end

%% load or build classifier from training data
if ~isempty(classifierObj)
    classifier = classifierObj;
elseif ~isempty(classifierPath)
    classifier = uiopen;
else   
    [data,attributes,classIndex] = weka2matlab(trainingData);
    numAttributes = classIndex - 1;
    trainingResponse = data(:, classIndex);
    data = data(:, 1:numAttributes);
    
    classifier = TreeBagger(64,data,trainingResponse,'OOBPredictorImportance','Off', 'Method','classification');
    
    clear data; clear trainingResponse; clear trainingData; clear arffLoader;
end

%% generate test data by filtering the image


attributes = getAttributes(trainingData);
numAttributes = numel(classifier.PredictorNames);
testMatrix = zeros(numInstances, numAttributes);


for i = 1:numAttributes
    
    att = attributes{i};
    disp(['Generating feature ', num2str((i)), '/', num2str(numAttributes) , ': ', att]);
    
    
    if i > 1
        filterType = regexp(att, '.*(?=_\d)', 'match');
        sigmas = regexp(att, '(\d[.]\d)|(\d\d[.]\d)', 'match');
        filteredIm = filterImage(im, filterType{1}, sigmas);
    else
        filteredIm = im;
    end
    
   
    testMatrix(:,i) = filteredIm(:);
    
end

[~, pLin] = predict(classifier,testMatrix);

%pLin = cellfun(@str2num, mask)

if dim == 2
    pMap = reshape(pLin(:, 2), [yDim xDim]);
elseif dim == 3
    pMap = reshape(pLin, [yDim xDim zDim]);
end

if displayFigures & dim==2
    figure(1); imshowpair(im, pMap, 'montage');
end




end