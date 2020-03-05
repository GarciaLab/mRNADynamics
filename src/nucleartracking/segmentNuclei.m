function segmentNuclei(Prefix, varargin)

displayFigures = false;
keepPool = false;
nWorkers = 1;
reSc = false;
thresh = .5;
fish = false;
algo = 'FastRandomForest';
maxDepth = 20;
nTrees = 64;
matlabLoader = true;


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

%%
disp(['Segmenting nuclei on ', Prefix, '...']);

[~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

dataRoot = fileparts(PreProcPath);
mlFolder = [dataRoot, filesep, 'training_data_and_classifiers', filesep];

[trainingNameExt, trainingFolder] = uigetfile([mlFolder, filesep, '*.*']);
trainingFile = [trainingFolder, filesep, trainingNameExt];
[~ ,trainingName] = fileparts(trainingNameExt);

load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo');

[~,hisMat] = makeMovieMats(Prefix, PreProcPath, nWorkers, FrameInfo);

% load([PreProcPath, filesep, Prefix, filesep, Prefix, '_hisMat.mat'], 'hisMat');
hisMat = double(hisMat);

%need to change this later. will be loaded from computerfolders
ramDrive = 'R:\';

nFrames = size(hisMat, 1);

pMap = zeros(size(hisMat, 1), size(hisMat, 2), size(hisMat, 3));

%%
%only need to make the classifier from training data
%once and not every frame

arffLoader = javaObject('weka.core.converters.ArffLoader'); %this constructs an object of the arffloader class
arffLoader.setFile(javaObject('java.io.File',trainingFile)); %construct an arff file object
trainingData= arffLoader.getDataSet;
trainingData.setClassIndex(trainingData.numAttributes - 1);

%remove the features matlab we can't (currently) generate in matlab
dim = 2;
[~,attributes,~] = weka2matlab(trainingData);
[~, ~, keepIndices, ~] = validateAttributes(attributes, dim);
trainingData = cleanArff(trainingData, keepIndices);

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

for f = 1:nFrames
    
    tic
    mean_dT = movmean(dT, [3, 0]);
    if f~=1, mean_dT = mean_dT(end); end
    
    if f~=1, tic, disp(['Making probability map for frame: ', num2str(f),...
            '. Estimated ', num2str(mean_dT*(nFrames-f)), ' minutes remaining.'])
    end
    im = squeeze(hisMat(f, :, :));
    pMap(f, :, :) = classifyImageWeka(im, trainingData,'tempPath', ramDrive,...
        'reSc', reSc, 'classifierObj', classifier, 'arffLoader', arffLoader, 'matlabLoader', matlabLoader);
    try waitbar(f/nFrames, wb); end
    dT(f)=toc/60;
    
end

try close(wb); end

mkdir([ProcPath, filesep, Prefix, filesep, Prefix]);
save([ProcPath, filesep, Prefix, filesep, Prefix, '_probHis.mat'], 'pMap', '-v7.3', '-nocompression');


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