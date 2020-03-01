function segmentNuclei(Prefix, varargin)

displayFigures = false;
keepPool = false;
nWorkers = 1;
reSc = false;
thresh = .5;
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
    meandT = movmean(dT, [3, 0]);
    meandT = meandT(end);
    
    if f~=1, tic, disp(['Making probability map for frame: ', num2str(f),...
            '. Estimated ', num2str(meandT*(nFrames-f)), ' minutes remaining.'])
    end
    im = squeeze(hisMat(f, :, :));
    pMap(f, :, :) = classifyImage(im, trainingData,'tempPath', ramDrive, 'reSc', reSc, 'classifierObj', classifier);
    waitbar(f/nFrames, wb);
    dT(f)=toc/60;
    
end

close(wb);

save([ProcPath, filesep, Prefix, filesep, Prefix, '_probHis.mat'], 'pMap', '-v7.3', '-nocompression');


%% Get Ellipses

if exist([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'] ,'file')
    
    ellipsePrompt = ('Ellipses.mat already exists. Do you want to overwrite?');
    ellipseAnswer = inputdlg(ellipsePrompt);
    
    if contains(ellipseAnswer,'y')
        
        %do morphology analysis to reduce our probability maps to a list of
        %ellipses
        
        nuclearMask = pMap > thresh;
        
        for f = 1:nFrames
            
            mask = squeeze(nuclearMask(f, :, :)) > thresh;
            pim = squeeze(pMap(f, :, :));
            
            ellipseStats = regionprops(mask, pim, {...
                'WeightedCentroid',...
                'MajorAxisLength',...
                'MinorAxisLength',...
                'Orientation'});
            
            mjx = [ellipseStats.MajorAxisLength];
            mnx = [ellipseStats.MinorAxisLength];
            ori =  [ellipseStats.Orientation];
            wcs= zeros(length(ellipseStats), 2);
            
            for r = 1:length(ellipseStats)
                wcs(r, :) = ellipseStats(r).WeightedCentroid;
            end
            
            %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
            Ellipses{f} = [wcs(r,2),wcs(r,1),mjx(r),mnx(r),ori(r),0,0,0];
            
        end
        
        save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses', '-v7.3', '-nocompression');
        
    end
    
end

%%
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