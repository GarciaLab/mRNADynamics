function filterMovieWeka(Prefix, varargin)

displayFigures = false;
keepPool = false;
nWorkers = 1;
reSc = false;
algo = 'FastRandomForest';
maxDepth = 20;
nTrees = 64;

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

[~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

dataRoot = fileparts(PreProcPath);
mlFolder = [dataRoot, filesep, 'training_data_and_classifiers', filesep];

[trainingNameExt, trainingFolder] = uigetfile([mlFolder, filesep, '*.*']);
trainingFile = [trainingFolder, filesep, trainingNameExt];
[~ ,trainingName] = fileparts(trainingNameExt);

load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo');

movieMat = makeMovieMats(Prefix, PreProcPath, nWorkers, FrameInfo); %ch z t y x

movieMat = double(movieMat);

%need to change this later. will be loaded from computerfolders
ramDrive = 'R:\';

nFrames = size(movieMat, 3);
nSlices = size(movieMat, 2);
yDim = size(movieMat, 4);
xDim = size(movieMat, 5);

pMap = zeros(yDim, xDim, nSlices, nFrames); % y x z f

%%
%only need to make the classifier from training data
%once and not every frame

arffLoader = javaObject('weka.core.converters.ArffLoader'); %this constructs an object of the arffloader class
arffLoader.setFile(javaObject('java.io.File',trainingFile)); %construct an arff file object
trainingData= arffLoader.getDataSet;
trainingData.setClassIndex(trainingData.numAttributes - 1);

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
save([trainingFolder, filesep, trainingName, '_', suffix '.model'], 'classifier')

%%
dT = [];
wb = waitbar(0, 'Classifying frames');

ch = 1;
for f = 1:nFrames
    
    tic
    meandT = movmean(dT, [3, 0]);
    meandT = meandT(end);
    
    if f~=1, tic, disp(['Making probability map for frame: ', num2str(f),...
            '. Estimated ', num2str(meandT*(nFrames-f)), ' minutes remaining.'])
    end
    im = squeeze(movieMat(ch,:,f, :, :));
    pMap(:, :, :, f) = classifyImage(im, trainingData,'tempPath', ramDrive, 'reSc', reSc, 'classifierObj', classifier);
    waitbar(f/nFrames, wb);
    dT(f)=toc/60;
    
end

close(wb);

save([ProcPath, filesep, Prefix, filesep, Prefix, '_probSpot.mat'], 'pMap', '-v7.3', '-nocompression');
      
end