% Script to test feasibility of porting ML segmentation over to matlab
clc; 
clear
close all
imtool close all;
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

%% Experiment with labeling data
fontSize = 16;
labelCell = cell(1,size(rawStack,3));
[yq, xq] = meshgrid(1:size(rawStack,1),1:size(rawStack,2));
xq = xq(:);
yq = yq(:);

% Move through Z-stack and label spot and background regions
zc = 1;
exit_flag = 0;
while ~exit_flag
    grayImage = mat2gray(rawStack(:,:,zc));
    % generate masks to show existing labels
    labelMat = labelCell{zc};
    notSpotMask = zeros(size(grayImage));
    spotMask = zeros(size(grayImage));
    if ~isempty(labelMat)        
        spotMask(labelMat(labelMat(:,1)==1,2)) = .5;        
        notSpotMask(labelMat(labelMat(:,1)==0,2)) = .5;
    end
    % generate RGB
    RGB = cat(3,mat2gray(grayImage+spotMask),mat2gray(grayImage+notSpotMask),grayImage);
    % generate figure
    label_fig = figure;
    axis on
    imshow(RGB);
    title(['Class labels (z slice: ' sprintf('%02d',zc) ')'])
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%     message = sprintf('Highlight "spot" and "notSpot" regions');
%     uiwait(msgbox(message));
    % User signs name here.
    ax = gca;
    hFH = drawfreehand(ax,'FaceAlpha',.2,'LineWidth',1);        
    
    % Get the xy coordinates of where they drew.
    xy = hFH.Position;
    in = inpolygon(xq,yq,xy(:,2),xy(:,1));
    delete(hFH);
    titleBarCaption = 'Enter Region Type';
    promptMessage = sprintf('');
    button = questdlg(promptMessage, titleBarCaption, 'Spot', 'NotSpot', 'Cancel');
    
    
    if strcmpi(button,'Spot')
        indices = find(in);
        lb_mat = [ones(size(indices)), indices];
        labelMat = vertcat(labelMat,lb_mat);
    elseif strcmpi(button,'NotSpot')
        indices = find(in);
        lb_mat = [zeros(size(indices)), indices];
        labelMat = vertcat(labelMat,lb_mat);
    end    
    labelCell{zc} = labelMat;
    
%     promptMessage = sprintf('Press:\n "Enter" to label another region\n "a/z" to move up or down in Z stack\n "m/n" to move to next/previous frame \n "x" to save and exit');
    promptMessage = sprintf('Specify next action');
    titleBarCaption = 'Continue?';    
    promptMessage = sprintf('');
    button = questdlg(promptMessage, titleBarCaption, 'Continue','Exit','Exit');
    if strcmpi(button,'Exit')
        close all
        exit_flag = 1;
        continue
    end
    
    titleBarCaption = 'Where to?';    
    promptMessage = sprintf('');
    button = questdlg(promptMessage, titleBarCaption, 'SameZ', 'Up', 'Down','Down');

    if strcmpi(button,'Up')
        zc = max([1,zc-1]);
    elseif strcmpi(button,'Down')
        zc = min([size(rawStack,3),zc+1]);       
    else
        continue
    end
end



%%

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