function D = generateSyntheticNuclearData(Prefix, varargin)

nWorkers = 1;
displayFigures = false;
keepPool = true;


k = 2;
p = .3;
s = .2;

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

javaaddpath('C:\Program Files\Weka-3-8-4\weka.jar','-end');
%%
disp(['Generating synthetic data for ', Prefix, '...']);

[~, ProcPath, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

dataRoot = fileparts(PreProcPath);
mlFolder = [dataRoot, filesep, 'training_data_and_classifiers', filesep];

[trainingNameExt, trainingFolder] = uigetfile([mlFolder, filesep, '*.*']);
trainingFile = [trainingFolder, filesep, trainingNameExt];
[~ ,trainingName] = fileparts(trainingNameExt);

arffLoader = javaObject('weka.core.converters.ArffLoader'); %this constructs an object of the arffloader class
arffLoader.setFile(javaObject('java.io.File',trainingFile)); %construct an arff file object
trainingData= arffLoader.getDataSet;
trainingData.setClassIndex(trainingData.numAttributes - 1);

trainingMat = weka2matlab(trainingData);

A = trainingMat(:, 1:end-1);

if min(A(:)) ~= 0 || max(A(:)) ~= 1
    
    warning('data not normalized to [0, 1]. normalizing now.');
    A = normalize(A, 'range'); %this normalizes each column separately. 
    
end

trainingMat(:, 1:end-1) = A; 


D = munge(trainingMat, k, p, s);

end