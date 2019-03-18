% Script to test feasibility of porting ML segmentation over to matlab
clear
close all
% specify file paths and test project
Prefix = '2018-10-25-embryo2_nc14_Linear unmixing';
dataRoot = 'E:\Liya\LivemRNA\Data\';

tp = 51;


preProcFiles = dir([dataRoot 'PreProcessedData/' Prefix '/*' sprintf('%03d',tp) '_z*ch02.tif']);
probProcFiles = dir([dataRoot 'ProcessedData/' Prefix '_/dogs/*' sprintf('%03d',tp) '_z*ch02.tif']);


% load raw and processed image stacks
rawStack = [];
probStack = [];
for i = 2:numel(preProcFiles)-1
    rawStack(:,:,i-1) = imread([preProcFiles(i).folder '/' preProcFiles(i).name]);
    probStack(:,:,i-1) = imread([probProcFiles(i).folder '/' probProcFiles(i).name]);
end
%%
% apply image filters to raw stack and use labelled stack to generate
% training data for random forest
binDist = round(bwdist(probStack>5000));
binDistVec = binDist(:);
% closeIDs = find(binDist==1);
training_array = probStack(:)>5000;
training_headers = {'labels'};   
sigmaVec = [1,3,5,7,9];
% feature_cell = {'Hessian_smallest','Hessian_largest','Difference_of_Gaussian',...
%     'Gaussian_blur','Structure_smallest','Structure_largest','Laplacian'};
feature_cell = {'Difference_of_Gaussian','Gaussian_blur','Laplacian'};

for i = 1:numel(feature_cell)
    feature = feature_cell{i};
    tic
    for j = 1:numel(sigmaVec)
        if strcmpi(feature,'Difference_of_Gaussian')
            ft_im = filterImage(rawStack, 'Difference_of_Gaussian', {sigmaVec(j), sigmaVec(j)*4});
        else
            ft_im = filterImage(rawStack, feature, {sigmaVec(j)});
        end
        training_array = [training_array ft_im(:)];
        training_headers = [training_headers{:} {[feature '_s' num2str(sigmaVec(j))]}];
    end
    toc
    disp(feature)
end


trainingTable = array2table(training_array, 'VariableNames',training_headers);

% create re-sampled version of table with balanced class instances
w0 = numel(trainingTable{:,1}) / sum(trainingTable{:,1}==0);
w1 = numel(trainingTable{:,1}) / sum(trainingTable{:,1}==1);
wt_vec = NaN(size(trainingTable{:,1}));
wt_vec(trainingTable{:,1}==1) = w1;
wt_vec(trainingTable{:,1}==0) = w0;

index_vec = 1:numel(trainingTable{:,1});
resamp_ids1 = find(trainingTable{:,1}==1);
binDistVec(resamp_ids1) = Inf;
resamp_ids2 = randsample(index_vec,sum(trainingTable{:,1}==1),true,1./binDistVec);

trainingTableResamp = trainingTable([resamp_ids1' resamp_ids2],:);
%% train radom forest classifier
if isempty(gcp('nocreate'))
    parpool(8);
end
paroptions = statset('UseParallel',true);

tic
rng(1)
Mdl = TreeBagger(15,trainingTableResamp,'labels','InBagFraction',1,'SampleWithReplacement',...
    'on','Options',paroptions,'OOBPrediction','On','PredictorSelection','curvature','OOBPredictorImportance','on');
toc


% examine relative importance of predictors
imp = Mdl.OOBPermutedPredictorDeltaError;
figure;
bar(imp);
title('Curvature Test');
ylabel('Predictor importance estimates');
xlabel('Predictors');
h = gca;
h.XTick = 1:numel(Mdl.PredictorNames);
h.XTickLabel = Mdl.PredictorNames;
h.XTickLabelRotation = 90;
h.TickLabelInterpreter = 'none';


% lets try training a stripped-down model
predictorsMinimized = {'Difference_of_Gaussian_s5','Gaussian_blur_s1','Gaussian_blur_s3','Gaussian_blur_s9','Laplacian_s7'};
predictorIndices = find(ismember(Mdl.PredictorNames,predictorsMinimized));

tic
rng(1)
MdlSimp = TreeBagger(15,trainingTableResamp(:,[1 predictorIndices]),'labels','InBagFraction',1,'SampleWithReplacement',...
    'on','Options',paroptions,'OOBPrediction','On','OOBPredictorImportance','on');
toc

figure;
oobErrorBaggedEnsemble = oobError(Mdl);
oobErrorBaggedEnsembleSimp = oobError(MdlSimp);
plot(oobErrorBaggedEnsemble)
hold on
plot(oobErrorBaggedEnsembleSimp)
legend('full','minimized')
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';

%%
tic
[mdl_fits, scores] = predict(MdlSimp,trainingTable);
toc

%%
predStack = reshape(scores(:,2),size(rawStack));
