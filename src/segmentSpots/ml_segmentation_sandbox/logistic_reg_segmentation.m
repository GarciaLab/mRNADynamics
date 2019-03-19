% Script to test feasibility of porting ML segmentation over to matlab
clear
close all
% specify file paths and test project
Prefix = '2018-10-25-embryo2_nc14_Linear unmixing';
dataRoot = 'E:\Liya\LivemRNA\Data\';

% hardcode path to sample set for now...
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
    for j = 1:numel(sigmaVec)
        if strcmpi(feature,'Difference_of_Gaussian')
            ft_im = filterImage(rawStack, 'Difference_of_Gaussian', {sigmaVec(j), sigmaVec(j)*4});
        else
            ft_im = filterImage(rawStack, feature, {sigmaVec(j)});
        end
        training_array = [training_array ft_im(:)];
        training_headers = [training_headers{:} {[feature '_s' num2str(sigmaVec(j))]}];
    end
end

trainingTable = array2table(training_array, 'VariableNames',training_headers);

% create re-sampled version of table with balanced class instances
w0 = numel(trainingTable{:,1}) / sum(trainingTable{:,1}==0);
w1 = numel(trainingTable{:,1}) / sum(trainingTable{:,1}==1);
index_vec = 1:numel(trainingTable{:,1});
resamp_ids1 = find(trainingTable{:,1}==1);
binDistVec(resamp_ids1) = Inf;
resamp_ids2 = randsample(index_vec,sum(trainingTable{:,1}==1),true,1./binDistVec);
trainingTableResamp = trainingTable([resamp_ids1' resamp_ids2],:);
%%% train radom forest classifier
if isempty(gcp('nocreate'))
    parpool(8);
end
%%%
% perform feature selection to determine minimal model
% First fit full model using all features
X = trainingTableResamp{:,2:end};
Y = categorical(trainingTableResamp{:,1});
[betaFull,devFull,statsFull] = mnrfit(X,Y);
% Perform forward feature selection
% maxdev = chi2inv(.90,1);     
opt = statset('display','iter','TolFun',.01*devFull,'TolTypeFun','abs');
inmodel = sequentialfs(@critfun,X,Y,'cv','none','nullmodel',true,'options',opt,...
                       'direction','forward');
[betaMin,devMin,statsMin] = mnrfit(X(:,inmodel),Y);  

%%% Label full stack 
pihat = mnrval(betaMin,trainingTable{:,[false inmodel]});
pdStackSpot = reshape(pihat(:,2),size(rawStack));
pdStackNotSpot = reshape(pihat(:,1),size(rawStack));
maxRaw = max(rawStack,[],3) / max(max(max(rawStack,[],3)));
% pdStackNotSpot = reshape(pihat(:,1),size(rawStack));
lbRGB = cat(3,mat2gray(max(pdStackSpot,[],3)+maxRaw),maxRaw,maxRaw);
qc_fig = figure;
imshow(lbRGB);


%% convenience functions
function dev = critfun(X,Y)
    model = fitglm(X,Y,'Distribution','binomial');
    dev = model.Deviance;
end