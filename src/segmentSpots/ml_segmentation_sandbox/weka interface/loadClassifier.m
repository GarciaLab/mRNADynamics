function [classifier, trainingData] = loadClassifier(trainingData, varargin)

cleanAttributes = false;
nWorkers = 1;
NumPredictorsToSample = 2;
nTrees = 64;

for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end


if cleanAttributes
    [~,attributes,~] = weka2matlab(trainingData);
    %remove the features we can't (currently) generate in matlab
    dim = 2;
    [~, ~, keepIndices, ~] = validateAttributes(attributes, dim);
    trainingData = cleanArff(trainingData, keepIndices);
end

[trainingMat,~,classIndex] = weka2matlab(trainingData);
numAttributes = classIndex - 1;
trainingResponse = trainingMat(:, classIndex);
trainingMat = trainingMat(:, 1:numAttributes);

rng(1650757608);
if nWorkers > 1
    startParallelPool(nWorkers, false, true);
end
paroptions = statset('UseParallel',nWorkers>1);
classifier = TreeBagger(nTrees,trainingMat,trainingResponse,...
    'OOBPredictorImportance','Off', 'Method','classification',...
    'NumPredictorsToSample', NumPredictorsToSample,...
    'Reproducible', true, 'MinLeafSize', 1, 'Surrogate','On', 'Options',paroptions);

end