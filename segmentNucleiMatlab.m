function pMap = segmentNucleiMatlab(Prefix, varargin)

displayFigures = false;
keepPool = false;
nWorkers = 1;
reSc = false;
thresh = .5;
fish = false;
algo = 'TreeBagger';
NumPredictorsToSample = 2;
maxDepth = 20;
nTrees = 64;
hisMat = [];
classifier = [];


%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

tic

warning('off', 'MATLAB:Java:DuplicateClass');
warning('off', 'MATLAB:javaclasspath:jarAlreadySpecified');
javaaddpath('C:\Program Files\Weka-3-8-4\weka.jar','-end');

%%
disp(['Segmenting nuclei on ', Prefix, '...']);

[~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

dataRoot = fileparts(PreProcPath);
mlFolder = [dataRoot, filesep, 'training_data_and_classifiers', filesep];

[trainingNameExt, trainingFolder] = uigetfile([mlFolder, filesep, '*.arff*']);
trainingFile = [trainingFolder, filesep, trainingNameExt];
[~ ,trainingName] = fileparts(trainingNameExt);

load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo');

if isempty(hisMat)
    %     [~,hisMat] = makeMovieMats(Prefix, PreProcPath, nWorkers, FrameInfo, 'loadMovie', false);
    load([PreProcPath, filesep, Prefix, filesep, Prefix, '_hisMat.mat'], 'hisMat');
    hisMat = double(hisMat);
end

nFrames = size(hisMat, 3);

pMap = zeros(size(hisMat, 1), size(hisMat, 2), size(hisMat, 3));

%%
%only need to make the classifier from training data
%once and not every frame

arffLoader = javaObject('weka.core.converters.ArffLoader'); %this constructs an object of the arffloader class
arffLoader.setFile(javaObject('java.io.File',trainingFile)); %construct an arff file object
trainingData= arffLoader.getDataSet;
trainingData.setClassIndex(trainingData.numAttributes - 1);

if isempty(classifier)
    
    [trainingMat,~,classIndex] = weka2matlab(trainingData);
    numAttributes = classIndex - 1;
    trainingResponse = trainingMat(:, classIndex);
    trainingMat = trainingMat(:, 1:numAttributes);
    
    rng(1650757608);
    paroptions = statset('UseParallel',nWorkers>1);
    classifier = TreeBagger(64,trainingMat,trainingResponse,...
        'OOBPredictorImportance','Off', 'Method','classification',...
        'NumPredictorsToSample', NumPredictorsToSample,...
        'Reproducible', true, 'MinLeafSize', 1, 'Surrogate','On', 'Options',paroptions);
    
    suffix = strrep(strrep(char(datetime(now,'ConvertFrom','datenum')), ' ', '_'), ':', '-');
    save([trainingFolder, filesep, trainingName, '_', suffix '_classifier.mat'], 'classifier')
    
end

%% make probability maps for each frame

if nWorkers > 1
    %parallel version
    startParallelPool(nWorkers, displayFigures, keepPool);
    hisMat = parallel.pool.Constant(hisMat);
    trainingData = parallel.pool.Constant(trainingData);
    classifier = parallel.pool.Constant(classifier);
    parfor f = 1:nFrames
        im = squeeze(hisMat.Value(:, :, f));
        pMap(:, :, f) = classifyImageMatlab(im, trainingData.Value,...
            'reSc', reSc, 'classifierObj', classifier.Value);
    end
else
    %non-parallel version
    dT = [];
    wb = waitbar(0, 'Classifying frames');
    profile off
    profile on
    for f = 1:nFrames
        
        tic
        mean_dT = movmean(dT, [3, 0]);
        if f~=1, mean_dT = mean_dT(end); end
        
        if f~=1, tic, disp(['Making probability map for frame: ', num2str(f),...
                '. Estimated ', num2str(mean_dT*(nFrames-f)), ' minutes remaining.'])
        end
        im = squeeze(hisMat(:, :, f));
        pMap(:, :, f) = classifyImageMatlab(im, trainingData, 'reSc', reSc, 'classifierObj', classifier);
        try waitbar(f/nFrames, wb); end
        dT(f)=toc/60;
    end
    
end

profile off; profile viewer;
try close(wb); end

mkdir([ProcPath, filesep, Prefix, filesep, Prefix]);
save([ProcPath, filesep, Prefix, '_', filesep, Prefix, '_probHis.mat'], 'pMap', '-v7.3', '-nocompression');


%% Get Ellipses

if exist([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'] ,'file')
    
    ellipsePrompt = ('Ellipses.mat already exists. Do you want to overwrite?');
    ellipseAnswer = inputdlg(ellipsePrompt);
    
    if contains(ellipseAnswer,'y')
        
        %do morphology analysis to reduce our probability maps to a list of
        %ellipses
        
        Ellipses = makeEllipses(pMap, thresh);
        
        %track nuclei will complain if there are frames with no ellipses,
        %so we'll fake it for now.
        fakeFrame = Ellipses(~cellfun(@isempty, Ellipses));
        fakeFrame =  fakeFrame{1};
        
        Ellipses(cellfun(@isempty, Ellipses)) = {fakeFrame};
        
        save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses', '-v7.3', '-nocompression');
        
    end
    
end

%% Tracking
%Decide whether we need to re-track
userPrompt = 'Do you want to track nuclei now?';

trackAnswer = inputdlg(userPrompt);
if contains(trackAnswer,'n')
    disp('Ellipses saved. Per user input, not tracking. Exiting.')
else
    opts = {};  if fish, opts = [opts, 'markandfind']; end
    disp('Ellipses saved. Running TrackNuclei.')
    TrackNuclei(Prefix,'NoBulkShift','ExpandedSpaceTolerance', 1.5, 'retrack', 'nWorkers', 1, opts{:});
end


end