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
addJavaPathsForLivemRNA()


if ischar(training), tempPath = fileparts(training);
else, tempPath = ''; end

classifier = [];
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
    scale = 5; %ar hardcoded temporarily til i figure out a non-sketch way to do this. do i standardize or normalize?
    im = scale*im./max(max(im));
end

%% load or build classifier from training data
if ischar(classifier)
    load(classifier, 'classifier');
elseif isempty(classifier)
  
    [classifier, trainingData] = loadClassifier(trainingData);
    
end

clear arffLoader;

%% generate test data by filtering the image


attributes = getAttributes(trainingData);
numAttributes = numel(classifier.PredictorNames);
testMatrix = zeros(numInstances, numAttributes);

usedFeatures = {};
ignoredFeatures = {};

for i = 1:numAttributes
    
    att = attributes{i};
    %     disp(['Generating feature ', num2str((i)), '/', num2str(numAttributes) , ': ', att]);
    
    
    if i > 1
        [filteredIm, successFlag]  = filterAttribute(att, im);
    else
        filteredIm = im;
        successFlag = true;
    end
    
    if successFlag
        testMatrix(:,i) = filteredIm(:);
        usedFeatures = [usedFeatures, att];
    else
        testMatrix(:,i) = nan(numel(im),1);
        ignoredFeatures = [ignoredFeatures, att];
    end
    
end

% disp(['features used: ', usedFeatures])
% disp(['features ignored: ', ignoredFeatures])

[~, pLin] = predict(classifier,testMatrix);

%pLin = cellfun(@str2num, mask)

if dim == 2
    pMap = reshape(pLin(:, 1), [yDim xDim]);
elseif dim == 3
    pMap = reshape(pLin(:,1), [yDim xDim zDim]);
end

if displayFigures 
    if dim==2
        figure(1); imshowpair(im, pMap, 'montage');
    elseif dim==3
        figure(1); imshowpair(max(im, [], 3),max(pMap, [], 3), 'montage');
end




end