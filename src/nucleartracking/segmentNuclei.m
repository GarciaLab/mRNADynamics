function pMap = segmentNuclei(Prefix, varargin)

cleanupObj = onCleanup(@myCleanupFun);
warning('off', 'MATLAB:Java:DuplicateClass');
warning('off', 'MATLAB:javaclasspath:jarAlreadySpecified');

%% Initialize variables
displayFigures = false;
keepPool = false;
nWorkers = 1;
shouldRescaleTrainingData = false;
probabilityThreshold = .5;
classificationAlgorithm = 'TreeBagger';
NumPredictorsToSample = 2;
maxDepth = 20;
nTrees = 64;
hisMat = [];
shouldBalanceClasses = false; %resample to balance classes
cleanAttributes = false;
frameRange = [];
makeEllipses=false;
classifier = [];
classifyMethod = 'matlab';
tempPath = 'S:\livemRNATempPath\';
if ~exist(tempPath, 'dir')
    mkdir(tempPath);
end
matlabLoader = true;
parFrame = false;
parInstances = false;


%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

if strcmpi(classifyMethod, 'weka')
    classificationAlgorithm = 'FastRandomForest';
end


tic

addJavaPathsForLivemRNA()

%%
disp(['Segmenting nuclei on ', Prefix, '...']);

thisExperiment = liveExperiment(Prefix);

mlFolder = thisExperiment.MLFolder;

[trainingNameExt, trainingFolder] = uigetfile([mlFolder, filesep, '*.arff*']);
trainingFile = [trainingFolder, filesep, trainingNameExt];
[~ ,trainingName] = fileparts(trainingNameExt);


if isempty(frameRange)
    frameRange = [1, thisExperiment.nFrames];
end
if isempty(hisMat)
    hisMat = getHisMat(thisExperiment);
end


ySize = size(hisMat, 1);
xSize = size(hisMat, 2);
nFrames = size(hisMat, 3);

pMap = zeros(size(hisMat, 1), size(hisMat, 2), size(hisMat, 3));

%%
%only need to make the classifier from training data
%once and not every frame
if strcmpi(classifyMethod, 'matlab')
    
    trainingData = loadArff(trainingFile, 'balance', shouldBalanceClasses);
    arffLoader = [];
    if ~exist(classifier, 'var') || isempty(classifier)
        
        [classifier, trainingData] = loadClassifier(trainingData, 'cleanAttributes', cleanAttributes,...
            'NumPredictorsToSample', NumPredictorsToSample, 'nTrees', nTrees, 'nWorkers', nWorkers);
        compactClassifier = compact(classifier);
        suffix = strrep(strrep(char(datetime(now,'ConvertFrom','datenum')), ' ', '_'), ':', '-');
        save([trainingFolder, filesep, trainingName, '_', suffix '_classifier.mat'], 'compactClassifier', '-v7')
        
    end
    
elseif strcmpi(classifyMethod, 'weka')
    
    [trainingData, arffLoader] = loadArff(trainingFile, 'balance', shouldBalanceClasses);
    
    %remove the features we can't (currently) generate in matlab
    dim = 2;
    [~,attributes,~] = weka2matlab(trainingData);
    [~, ~, keepIndices, ~] = validateAttributes(attributes, dim);
    trainingData = cleanArff(trainingData, keepIndices);
    
    classifier =hr.irb.fastRandomForest.FastRandomForest;
    options = {'-I', num2str(nTrees), '-threads', num2str(nWorkers), '-K', '2', '-S', '-1650757608', '-depth', num2str(maxDepth)};
    
    switch classificationAlgorithm
        case 'FastRandomForest'
            %default
        case 'RandomForest'
            classifier = weka.classifiers.trees.RandomForest;
            options = {'-I', num2str(nTrees), '-K', '2', '-S', '-1650757608','-depth', num2str(maxDepth)};
        otherwise
            warning('classification algorithm not recognized. defaulting to fast random forest');
    end
    
    classifier.setOptions(options);
    classifier.buildClassifier(trainingData);
    suffix = strrep(strrep(char(datetime(now,'ConvertFrom','datenum')), ' ', '_'), ':', '-');
%     save([trainingFolder, filesep, trainingName, '_', suffix '.model'], 'classifier', '-v7.3')
    
end


%% make probability maps for each frame

if parFrame
    %parallel version
    startParallelPool(nWorkers, displayFigures, keepPool);
    hisMat = parallel.pool.Constant(hisMat);
    trainingData = parallel.pool.Constant(trainingData);
    classifier = parallel.pool.Constant(classifier);
    parfor f = 1:nFrames
        
        hisFrame = hisMat.Value(:, :, f);
        pMap(:, :, f) = classifyImageNuclear(hisFrame, trainingData.Value,...
            'tempPath', tempPath,...
            'shouldRescaleTrainingData', shouldRescaleTrainingData, 'classifierObj',...
            classifier.Value, 'arffLoader', arffLoader, 'matlabLoader', matlabLoader,...
            'classifyMethod', classifyMethod);
        
    end
    
else
    %non-parallel version
    deltaT = [];
    for f = 1:nFrames
        
        tic
        mean_dT = movmean(deltaT, [3, 0]);
        if f~=1, mean_dT = mean_dT(end); end
        
        if f~=1, tic, disp(['Making probability map for frame: ', num2str(f),...
                '. Estimated ', num2str(mean_dT*(nFrames-f)), ' minutes remaining.']); end
        
        hisFrame = hisMat(:, :, f);
        
        pMap(:, :, f) = classifyImageNuclear(hisFrame, trainingData,'tempPath', tempPath,...
            'shouldRescaleTrainingData', shouldRescaleTrainingData, 'classifier', classifier,...
            'arffLoader', arffLoader, 'matlabLoader', matlabLoader,...
            'parallelizeInstances', parInstances, 'displayFigures', displayFigures,...
            'classifyMethod', classifyMethod);
        
        deltaT(f)=toc/60;
        
    end
    
end

[~, ProcPath] = DetermineLocalFolders(Prefix);
procFolder = [ProcPath, filesep, Prefix, '_'];
mkdir(procFolder);
probHisFile = [procFolder, filesep, Prefix, '_probHis.mat'];

livemRNAImageMatSaver([procFolder, filesep, Prefix, '_probHis.mat'],...
            pMap);
        
%% Make ellipses from generated probability maps
if makeEllipses
    Ellipses = makeEllipses(pMap, probabilityThreshold);
    %track nuclei will complain if there are frames with no ellipses,
    %so we'll fake it for now.
    fakeFrame = Ellipses(~cellfun(@isempty, Ellipses));
    fakeFrame =  fakeFrame{1};
    Ellipses(cellfun(@isempty, Ellipses)) = {fakeFrame};
    save([thisExperiment.resultsFolder, 'Ellipses.mat'], 'Ellipses', '-v6');
    TrackNuclei(Prefix, 'retrack');
end

end