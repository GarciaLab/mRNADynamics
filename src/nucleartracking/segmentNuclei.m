function pMap = segmentNuclei(Prefix, varargin)

cleanupObj = onCleanup(@myCleanupFun);
warning('off', 'MATLAB:Java:DuplicateClass');
warning('off', 'MATLAB:javaclasspath:jarAlreadySpecified');

%% Initialize variables

%note that currently the parallelized code is buggy, so it's recommended to
%avoid it. 
displayFigures = false;
keepPool = false;
nWorkers = 1;
parFrame = false;
parInstances = false;


shouldRescaleTrainingData = false;
probabilityThreshold = .5;
classificationAlgorithm = 'TreeBagger';
NumPredictorsToSample = 2; % default value works well
maxDepth = 20; %RF tree height. generally default is fine
nTrees = 64; % generally fine.

%resample to balance classes. 
%currently produces poor results. not recommended
shouldBalanceClasses = false; 

%not recommended. 
cleanAttributes = false;
frameRange = [];

%untested. not recommended.
makeEllipses=false;

classifier = [];
classifyMethod = 'weka'; %matlab is faster. weka is more accurate
tempPath = 'S:\livemRNATempPath\';
if ~exist(tempPath, 'dir')
    mkdir(tempPath);
end
matlabLoader = true;



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

liveExperiment = LiveExperiment(Prefix);

mlFolder = liveExperiment.MLFolder;

[trainingNameExt, trainingFolder] = uigetfile([mlFolder, filesep, '*.arff*']);
trainingFile = [trainingFolder, filesep, trainingNameExt];
[~ ,trainingName] = fileparts(trainingNameExt);


if isempty(frameRange)
    frameRange = [1,liveExperiment.nFrames]; %#ok<*NASGU>
end

hisMat = getHisMat(liveExperiment);

pMap = zeros(size(hisMat, 1), size(hisMat, 2), frameRange(2)-frameRange(1)+1);

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
    
    classifier = javaObject('hr.irb.fastRandomForest.FastRandomForest');
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
    startParallelPool(nWorkers, displayFigures, keepPool); %#ok<*UNRCH>
    hisMat = parallel.pool.Constant(hisMat);
    trainingData = parallel.pool.Constant(trainingData);
    classifier = parallel.pool.Constant(classifier);
    parfor f = frameRange
        
        hisFrame = hisMat.Value(:, :, f);
        pMap(:, :, f) = classifyImageNuclear(hisFrame, trainingData.Value,...
            'tempPath', tempPath,...
            'shouldRescaleTrainingData', shouldRescaleTrainingData, 'classifier',...
            classifier.Value, 'arffLoader', arffLoader, 'matlabLoader', matlabLoader,...
            'classifyMethod', classifyMethod);
        
    end
    
else
    %non-parallel version
    deltaT = [];
    n = 0;
    for f = frameRange(1):frameRange(2)
        
        n = n + 1;
        tic
        mean_dT = movmean(deltaT, [3, 0]);
        if n~=1, mean_dT = mean_dT(end); end
        
        if n~=1, tic, disp(['Making probability map for frame: ', num2str(f),...
                '. Estimated ', num2str(mean_dT*(frameRange(2)-frameRange(1)-n)), ' minutes remaining.']); end
        
        hisFrame = hisMat(:, :, f);
        
        pMap(:, :, n) = classifyImageNuclear(hisFrame, trainingData,'tempPath', tempPath,...
            'shouldRescaleTrainingData', shouldRescaleTrainingData, 'classifier', classifier,...
            'arffLoader', arffLoader, 'matlabLoader', matlabLoader,...
            'parallelizeInstances', parInstances, 'displayFigures', displayFigures,...
            'classifyMethod', classifyMethod);
        
        deltaT(n)=toc/60;
        
    end
    
end

[~, ProcPath] = DetermineLocalFolders(Prefix);
procFolder = [ProcPath, filesep, Prefix, '_'];
mkdir(procFolder);
probHisFile = [procFolder, filesep, 'probHis.tif'];

imwrite(pMap(:, :, 1), probHisFile);
            
for f = 2:size(pMap, 3)

    imwrite(pMap(:, :, f), probHisFile, 'WriteMode', 'append');

end
            
        
%% Make ellipses from generated probability maps
if makeEllipses
    Ellipses = makeEllipses(pMap, probabilityThreshold);
    %track nuclei will complain if there are frames with no ellipses,
    %so we'll fake it for now.
    fakeFrame = Ellipses(~cellfun(@isempty, Ellipses));
    fakeFrame =  fakeFrame{1};
    Ellipses(cellfun(@isempty, Ellipses)) = {fakeFrame};
    save([liveExperiment.resultsFolder, 'Ellipses.mat'], 'Ellipses', '-v6');
    TrackNuclei(Prefix, 'retrack');
end

end