function filterMovieWeka(Prefix, varargin)

displayFigures = false;
keepPool = false;
nWorkers = 1;
optionalResults = '';
reSc = false;
algo = 'FastRandomForest';
maxDepth = 20;
nTrees = 64;
balance = false; %resample to balance classes
cleanAttributes = false;
frameRange = [];
ramDrive = 'R:\';
classifyWithMatlab = false;
classifyWithWeka = true;
matlabLoader = true;
parFrame = false;
parInstances = true;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

startParallelPool(nWorkers, displayFigures, keepPool);

warning('off', 'MATLAB:Java:DuplicateClass');
warning('off', 'MATLAB:javaclasspath:jarAlreadySpecified');
heapSize = java.lang.Runtime.getRuntime.maxMemory;
if heapSize < 1E10 && ~ ignoreMemoryCheck
    error('Please increase your Java heap memory allocation to at least 10GB (Home -> Preferences -> General -> Java Heap Memory.');
end

%%
disp(['Filtering ', Prefix, '...']);

cleanupObj = onCleanup(@myCleanupFun);

thisExperiment = liveExperiment(Prefix);

[~,ProcPath,DropboxFolder,~, PreProcPath,~, Prefix, ~,~,~,~,~, ~, ~, movieDatabase]...
    = readMovieDatabase(Prefix, optionalResults);
% 
% [~, ~, ~, ~, ~, ~, Channel1, Channel2, ~, ~, ~, ~, ~, ...
%     ~, ~, ~, ~, ~, ~, ~, Channel3, ~, ~] =...
%     getExperimentDataFromMovieDatabase(Prefix, movieDatabase);

% coats = getCoatChannel(Channel1, Channel2, Channel3);
coats = thisExperiment.spotChannel;

nCh = length(coats); 

if nCh > 1
    error('2spot2color not supported yet. Talk to AR.')
end

ch = coats(1); %assumes the experiment is _not_ 2spot2color

mlFolder = thisExperiment.MLFolder;


[trainingNameExt, trainingFolder] = uigetfile([mlFolder, filesep, '*.arff']);
trainingFile = [trainingFolder, filesep, trainingNameExt];
[~ ,trainingName] = fileparts(trainingNameExt);

% load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo');
FrameInfo = getFrameInfo(thisExperiment);

movieMat = getMovieMat(thisExperiment);

movieMat = double(squeeze(movieMat(:, :, :, :, ch)));

%need to change this later. will be loaded from computerfolders

nFrames = size(movieMat, 4);
nSlices = size(movieMat, 3);
yDim = size(movieMat, 1);
xDim = size(movieMat, 2);

pMap = zeros(yDim, xDim, nSlices, nFrames); % y x z f

%%
%only need to make the classifier from training data
%once and not every frame

arffLoader = javaObject('weka.core.converters.ArffLoader'); %this constructs an object of the arffloader class
arffLoader.setFile(javaObject('java.io.File',trainingFile)); %construct an arff file object
trainingData= arffLoader.getDataSet;
trainingData.setClassIndex(trainingData.numAttributes - 1);

%remove the features matlab we can't (currently) generate in matlab
if cleanAttributes
    dim = 3;
    [~,attributes,~] = weka2matlab(trainingData);
    [~, ~, keepIndices, ~] = validateAttributes(attributes, dim);
    trainingData = cleanArff(trainingData, keepIndices);
end

classifier = javaObject('hr.irb.fastRandomForest.FastRandomForest');
options = {'-I', num2str(nTrees), '-threads', num2str(nWorkers), '-K', '2', '-S', '-1650757608', '-depth', num2str(maxDepth)};

switch algo
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
save([trainingFolder, filesep, trainingName, '_', suffix '.model'], 'classifier', '-v7.3')

%%
if parFrame
    
    %parallel version
    movieMat = parallel.pool.Constant(movieMat);
    trainingData = parallel.pool.Constant(trainingData);
    classifier = parallel.pool.Constant(classifier);
    parfor f = 1:nFrames
        im = squeeze(movieMat.Value(:, :, :,f));
        pMap(:, :, f) = classifyImageMatlab(im, trainingData.Value,...
            'reSc', reSc, 'classifier', classifier.Value);
    end
    
else
    
    dT = [];
    wb = waitbar(0, 'Classifying frames');

    for f = 1:nFrames

        tic
        mean_dT = movmean(dT, [3, 0]);
        if f~=1, mean_dT = mean_dT(end); end

        if f~=1, tic, disp(['Making probability map for frame: ', num2str(f),...
                '. Estimated ', num2str(mean_dT*(nFrames-f)), ' minutes remaining.'])
        end
        im = squeeze(movieMat(:, :, :,f));
        pMap(:, :, :, f) = classifyImageWeka(im, trainingData,'tempPath',...
            ramDrive, 'reSc', reSc, 'classifier', classifier, 'par', parInstances);
        try waitbar(f/nFrames, wb); end
        dT(f)=toc/60;

    end
    
    try close(wb); end

end

mkdir([ProcPath, filesep, Prefix, '_']);

probSpotFile = [ProcPath, filesep, Prefix, '_', filesep, Prefix, '_probSpot.mat'];
if whos(var2str(pMap)).bytes < 2E9
    save(probSpotFile, 'pMap', '-v6')
else
newmatic(probSpotFile,true,...
            newmatic_variable('pMap', 'double', [ySize, xSize, zSize, nFrames], [ySize, xSize, 1, 1]));
end

      
end