function pMap = segmentNuclei(Prefix, varargin)

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
balance = false; %resample to balance classes
cleanAttributes = false;
frameRange = [];
makeEllipses=false;
doTracking = false;
classifyWithMatlab = true;
classifyWithWeka = false;
ramDrive = 'R:\';
matlabLoader = true;


%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

if classifyWithWeka
    algo = 'FastRandomForest';
    classifyWithMatlab = false;
end


tic

warning('off', 'MATLAB:Java:DuplicateClass');
warning('off', 'MATLAB:javaclasspath:jarAlreadySpecified');
javaaddpath('C:\Program Files\Weka-3-8-4\weka.jar','-end');

%%
disp(['Segmenting nuclei on ', Prefix, '...']);

thisExperiment = liveExperiment(Prefix);

[~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);
% ProcPath = thisExperiment.procFolder;
% PreProcPath = thisExperiment.preFolder;
% DropboxFolder = thisExperiment.resultsFolder;
mlFolder = thisExperiment.MLFolder;

% dataRoot = fileparts(PreProcPath);
% mlFolder = [dataRoot, filesep, 'training_data_and_classifiers', filesep];

[trainingNameExt, trainingFolder] = uigetfile([mlFolder, filesep, '*.arff*']);
trainingFile = [trainingFolder, filesep, trainingNameExt];
[~ ,trainingName] = fileparts(trainingNameExt);

load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo');

if isempty(hisMat)
    
    hisFile = [PreProcPath, filesep, Prefix, filesep, Prefix, '_hisMat.mat'];
    hisMat = double(loadHisMat(hisFile, 'frameRange', frameRange));
    
end

ySize = size(hisMat, 1);
xSize = size(hisMat, 2);
nFrames = size(hisMat, 3);

pMap = zeros(size(hisMat, 1), size(hisMat, 2), size(hisMat, 3));

%%
%only need to make the classifier from training data
%once and not every frame
if classifyWithMatlab
    
    trainingData = loadArff(trainingFile, 'balance', balance);
    if isempty(classifier)
        
        [classifier, trainingData] = loadClassifier(trainingData, 'cleanAttributes', cleanAttributes,...
            'NumPredictorsToSample', NumPredictorsToSample, 'nTrees', nTrees);
        
        suffix = strrep(strrep(char(datetime(now,'ConvertFrom','datenum')), ' ', '_'), ':', '-');
        save([trainingFolder, filesep, trainingName, '_', suffix '_classifier.mat'], 'classifier', '-v7.3')
    end
    
elseif classifyWithWeka
    
    arffLoader = javaObject('weka.core.converters.ArffLoader'); %this constructs an object of the arffloader class
    arffLoader.setFile(javaObject('java.io.File',trainingFile)); %construct an arff file object
    trainingData= arffLoader.getDataSet;
    trainingData.setClassIndex(trainingData.numAttributes - 1);
    
    %remove the features matlab we can't (currently) generate in matlab
    dim = 2;
    [~,attributes,~] = weka2matlab(trainingData);
    [~, ~, keepIndices, ~] = validateAttributes(attributes, dim);
    trainingData = cleanArff(trainingData, keepIndices);
    
    classifier =hr.irb.fastRandomForest.FastRandomForest;
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
        if classifyWithMatlab
            pMap(:, :, f) = classifyImageMatlab(im, trainingData.Value,...
                'reSc', reSc, 'classifierObj', classifier.Value);
        elseif classifyWithWeka
            pMap(:, :, f) = classifyImageWeka(im, trainingData.Value,'tempPath', ramDrive,...
                'reSc', reSc, 'classifierObj', classifier.Value, 'arffLoader', arffLoader, 'matlabLoader', matlabLoader);
        end
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
        
        if classifyWithMatlab
            pMap(:, :, f) = classifyImageMatlab(im, trainingData, 'reSc', reSc, 'classifier', classifier);
            
        elseif classifyWithWeka
            pMap(:, :, f) = classifyImageWeka(im, trainingData,'tempPath', ramDrive,...
                'reSc', reSc, 'classifier', classifier, 'arffLoader', arffLoader, 'matlabLoader', matlabLoader);
            
        end
        try waitbar(f/nFrames, wb); end
        dT(f)=toc/60;
    end
    
end

profile off; profile viewer;

try close(wb); end

mkdir([ProcPath, filesep, Prefix, '_']);
probHisFile = [ProcPath, filesep, Prefix, '_', filesep, Prefix, '_probHis.mat'];
if whos(var2str(pMap)).bytes < 2E9
    save(probHisFile, 'pMap', '-v6')
else
    newmatic(probHisFile,true,...
        newmatic_variable('pMap', 'double', [ySize, xSize, nFrames], [ySize, xSize, 1]));
end

%% Get Ellipses
if makeEllipses
    
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
            
            save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses', '-v6');
            
        end
        
    end
end

if doTracking
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

end