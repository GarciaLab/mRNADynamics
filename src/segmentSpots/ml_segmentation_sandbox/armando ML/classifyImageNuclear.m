function pMap = classifyImageNuclear(im, training, varargin)
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
addJavaPathsForLivemRNA()
%% PARAMS, OPTIONS
displayFigures = false;
classifier = [];
classifierPath = '';
arffLoader = [];
shouldRescaleTrainingData = false;
useMatlabLoader = true;
parallelizeInstances = false;
classifyMethod = 'matlab'; %also accepts 'weka'
if ischar(training), tempPath = fileparts(training);
else, tempPath = ''; end

if ~displayFigures
    cleanupObj = onCleanup(@myCleanupFun);
end

for k = 1:2:(numel(varargin)-1)
    if k ~= numel(varargin)
        eval([varargin{k} '=varargin{k+1};']);
    end
end
clear varargin;

xDim=size(im, 2); yDim = size(im, 1);
dim = length(size(im));
if dim==3
    zDim = size(im, 3);
end

numInstances = numel(im);

warning('off', 'MATLAB:Java:DuplicateClass');
warning('off', 'MATLAB:javaclasspath:jarAlreadySpecified');

addJavaPathsForLivemRNA()

%% load up training data
if ischar(training)
    [trainingData, arffLoader] = loadArff(training);
else
    trainingData = training;
end


%normalize data to the max of the training set for better classification
if shouldRescaleTrainingData
    scale = 5; %ar hardcoded temporarily til i figure out a non-sketch way to do this. do i standardize or normalize?
    im = scale*im./max(max(im));
end

%% load or build classifier from training data

shouldLoadClassifier = isempty(classifier) && ~isempty(classifierPath);
shouldCreateClassifier = isempty(classifier) && isempty(classifierPath);

if shouldLoadClassifier
    switch classifyMethod
        case 'matlab'
            load(classifier, 'classifier');
        case 'weka'
            classifier = weka.core.SerializationHelper.read(classifierPath);
    end
elseif shouldCreateClassifier
    switch classifyMethod
        case 'matlab'
            [classifier, trainingData] = loadClassifier(trainingData);
        case 'weka'
            classifier = javaObject('hr.irb.fastRandomForest.FastRandomForest');
            options = {'-I', '20', '-threads', '1', '-K', '2', '-S', '-1650757608'};
            classifier.setOptions(options);
            classifier.buildClassifier(trainingData);
            if displayFigures, classifier.toString, end %#ok<*UNRCH>
    end
end
clear arffLoader;

%% generate test data by filtering the image


attributes = getAttributes(trainingData);
if strcmpi(classifyMethod, 'matlab')
    lastInd = numel(classifier.PredictorNames);
    testMatrix = zeros(numInstances, lastInd);
elseif strcmpi(classifyMethod, 'weka')
    if useMatlabLoader
        testMatrix = zeros(numInstances, trainingData.numAttributes-1);
        lastInd = numel(attributes) - 2;
    else
        testData = initializeEmptyDataSet(arffLoader, numInstances);
        testData.setClassIndex(testData.numAttributes-1);
        lastInd = numel(attributes) - 1;
    end
end


usedFeatures = {};
ignoredFeatures = {};

oneIndexingFlag = double(strcmpi(classifyMethod, 'matlab'));

for k = oneIndexingFlag:lastInd
    
    att = attributes{k + ~oneIndexingFlag};
    
    if k > oneIndexingFlag
        [filteredIm, successFlag]  = filterAttribute(att, im);
    else
        filteredIm = im;
        successFlag = true;
    end
    
    if strcmpi(classifyMethod, 'matlab') || useMatlabLoader
        
        if successFlag
            testMatrix(:,k+~oneIndexingFlag) = filteredIm(:);
            usedFeatures = [usedFeatures, att];
        else
            testMatrix(:,k+~oneIndexingFlag) = nan(numel(im),1);
            ignoredFeatures = [ignoredFeatures, att];
        end
        
    elseif strcmpi(classifyMethod, 'weka') && ~useMatlabLoader
        
        for k = 1:numel(filteredIm)
            testData = setInstancesVal(testData, k-1, k, filteredIm(k));
        end
        
    end %classifymethod conditional
    
end %feature loop 


if strcmpi(classifyMethod, 'matlab')
    
    [~, pLinear] = predict(classifier,testMatrix);
    
    pLinear = pLinear(:, 1);
    
    
elseif strcmpi(classifyMethod, 'weka')
    
    %now turn data matrix into weka arff
    
    if useMatlabLoader
        testData = mat2ascii2dataSet(testMatrix, tempPath, trainingData);
    end
    clear testMatrix; clear filteredIm;
    
    %check compatability between arff and data
    compatible = testData.equalHeaders(trainingData); 
    if ~compatible
        error(char(testData.equalHeadersMsg(trainingData)));
    end
    
    %this step is unbearably slow. ~1min per frame. needs a massive overhaul
    nInstances = testData.numInstances;
    pLinear = zeros(nInstances, 1);
    
    if parallelizeInstances
        
        testDataConstant = parallel.pool.Constant(testData);
        classifierConstant = parallel.pool.Constant(classifier);
        parfor k = 1:nInstances
            pLinear(k) = classifyInstance(k, testDataConstant.Value,...
                classifierConstant.Value);
        end
        
    else
        
        for k = 1:nInstances
            pLinear(k) = classifyInstance(k, testData, classifier);
        end
        
    end %parallelize instances conditional
        
end %classifymethod conditional

if dim == 2
    pMap = reshape(pLinear, [yDim xDim]);
    if displayFigures, figure; imshowpair(im, pMap, 'montage'); end
elseif dim == 3
    pMap = reshape(pLinear, [yDim xDim zDim]);
    if displayFigures
        figure; 
        imshowpair(max(im, [], 3),max(pMap, [], 3), 'montage'); 
    end
end

end


function probability = classifyInstance(instanceIndex, testData, classifier)

k = instanceIndex - 1;
instance = testData.instance(k);
temp = classifier.distributionForInstance(instance);
probability = temp(1); 

end

function testData = initializeEmptyDataSet(arffLoader, numInstances)

testData = arffLoader.getStructure;
numAttributes = testData.numAttributes;
for i = 1:numInstances
    inst = weka.core.DenseInstance(numAttributes+1);
    testData.add(inst);
end

end

function testData = mat2ascii2dataSet(mat, tempPath, trainingData)

tempFile = [tempPath, 'temp.data'];
save(tempFile,var2str(mat),'-ascii', '-double');
loader = weka.core.converters.MatlabLoader;
loader.setFile( java.io.File(tempFile) );

testData = loader.getDataSet;
testData.insertAttributeAt(trainingData.classAttribute,trainingData.classIndex);
testData.setClassIndex(testData.numAttributes-1);


f = weka.filters.unsupervised.attribute.Remove;
f.setInputFormat(trainingData);
%match test header to training header
testData = javaMethod('useFilter','weka.filters.Filter', testData, f); 

end