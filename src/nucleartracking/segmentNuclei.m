function pMap = segmentNuclei(Prefix, varargin)

cleanupObj = onCleanup(@myCleanupFun);


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
classifier = [];
balance = false; %resample to balance classes
cleanAttributes = false;
frameRange = [];
makeEllipses=false;
doTracking = false;
classifyWithMatlab = true;
classifyWithWeka = false;
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

if classifyWithWeka
    classificationAlgorithm = 'FastRandomForest';
    classifyWithMatlab = false;
end


tic

warning('off', 'MATLAB:Java:DuplicateClass');
warning('off', 'MATLAB:javaclasspath:jarAlreadySpecified');

addJavaPathsForLivemRNA()

%%
disp(['Segmenting nuclei on ', Prefix, '...']);

thisExperiment = liveExperiment(Prefix);

[~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);
% ProcPath = thisExperiment.procFolder;
% PreProcPath = thisExperiment.preFolder;
% DropboxFolder = thisExperiment.resultsFolder;
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
        hisFrame = hisMat.Value(:, :, f);
        if classifyWithMatlab
            pMap(:, :, f) = classifyImageMatlab(hisFrame, trainingData.Value,...
                'shouldRescaleTrainingData', shouldRescaleTrainingData, 'classifierObj', classifier.Value);
        elseif classifyWithWeka
            pMap(:, :, f) = classifyImageWeka(hisFrame, trainingData.Value,...
                'tempPath', tempPath,...
                'shouldRescaleTrainingData', shouldRescaleTrainingData, 'classifierObj',...
                classifier.Value, 'arffLoader', arffLoader, 'matlabLoader', matlabLoader);
        end
    end
else
    %non-parallel version
    deltaT = [];
    profile off
    profile on
    for f = 1:nFrames
        
        tic
        mean_dT = movmean(deltaT, [3, 0]);
        if f~=1, mean_dT = mean_dT(end); end
        
        if f~=1, tic, disp(['Making probability map for frame: ', num2str(f),...
                '. Estimated ', num2str(mean_dT*(nFrames-f)), ' minutes remaining.'])
        end
        hisFrame = hisMat(:, :, f);
        
        if classifyWithMatlab
            pMap(:, :, f) = classifyImageMatlab(hisFrame, trainingData, 'reSc',...,
                shouldRescaleTrainingData, 'classifier', classifier, 'displayFigures', displayFigures);
            
        elseif classifyWithWeka
            pMap(:, :, f) = classifyImageWeka(hisFrame, trainingData,'tempPath', tempPath,...
                'reSc', shouldRescaleTrainingData, 'classifier', classifier,...
                'arffLoader', arffLoader, 'matlabLoader', matlabLoader, 'displayFigures', displayFigures);
            
        end
        deltaT(f)=toc/60;
    end
    
end

profile off; profile viewer;

mkdir([ProcPath, filesep, Prefix, '_']);
probHisFile = [ProcPath, filesep, Prefix, '_', filesep, Prefix, '_probHis.mat'];
if whos(var2str(pMap)).bytes < 2E9
    save(probHisFile, 'pMap', '-v6')
else
    newmatic(probHisFile,true,...
        newmatic_variable('pMap', 'double', [ySize, xSize, nFrames], [ySize, xSize, 1]));
end




%% Make ellipses from generated probability maps 
if makeEllipses
    if exist([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'] ,'file')
        ellipsePrompt = ('Ellipses.mat already exists. Do you want to overwrite?');
        ellipseAnswer = inputdlg(ellipsePrompt);
        if contains(ellipseAnswer,'y')
            %do morphology analysis to reduce our probability maps to a list of
            %ellipses
            Ellipses = makeEllipses(pMap, probabilityThreshold);
            %track nuclei will complain if there are frames with no ellipses,
            %so we'll fake it for now.
            fakeFrame = Ellipses(~cellfun(@isempty, Ellipses));
            fakeFrame =  fakeFrame{1};     
            Ellipses(cellfun(@isempty, Ellipses)) = {fakeFrame};
            save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses', '-v6');  
        end      
    end
end

end