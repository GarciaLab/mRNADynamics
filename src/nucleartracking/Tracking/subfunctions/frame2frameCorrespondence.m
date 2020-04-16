function [ mapping, nextNucleiXY, score, inverse_mapping, varargout ] =...
    frame2frameCorrespondence(FrameInfo, hisMat, frameNumber1, frameNumber2, ...
    nucleiFrame1, nucleusDiameter, precision, nucleiFrame2, manualMapping, ...
    shifts, ExpandedSpaceTolerance, NoBulkShift )
%FRAME2FRAMECORRESPONDENCE Summary of this function goes here
%   Detailed explanation goes here

%% 0. Initialize Parameters
space_resolution = getDefaultParameters(FrameInfo,'space resolution');
time_resolution = getDefaultParameters(FrameInfo,'time resolution');
maxShiftCorrection = getDefaultParameters(FrameInfo,'max Shift Correction', 'trackToTheNextFrame')*nucleusDiameter/space_resolution;
maxNucleusStep = getDefaultParameters(FrameInfo,'max Interphase Displacement')...
    *time_resolution/60*nucleusDiameter/space_resolution * ExpandedSpaceTolerance;
mapping = zeros(size(nucleiFrame1,1),1);

frame1 = double(hisMat(:, :, frameNumber1));
frame2 = double(hisMat(:, :, frameNumber2));

xDim = size(frame1, 2);
yDim = size(frame1, 1);

%% 1. Find overall movements between the two frames
if NoBulkShift
    vx = zeros(size(frame1));
    vy = zeros(size(frame2));
    varargout{1}(:,:,1) = vx;
    varargout{1}(:,:,2) = vy;
else
    if ~exist('shifts','var') || isempty(shifts)
        [x,y,vx,vy] = interpolatedShift(FrameInfo, frame1, frame2, maxShiftCorrection, [], [], precision);
        varargout{1}(:,:,1) = vx;
        varargout{1}(:,:,2) = vy;
    else
        vx = shifts(:,:,1);
        vy = shifts(:,:,2);
        
        varargout{1} = shifts;
    end
end

%[x,y,currentNucleiShiftedXY] = interpolatedShift(frame1, frame2, maxShiftCorrection, [], [], precision,nucleiFrame1);


%% 2. If not provided, find nuclei in the second frame
if ~exist('nucleiFrame2','var') || isempty(nucleiFrame2)
    [nextNucleiXY, intensity] = findNuclei(FrameInfo,hisMat, frameNumber2, nucleusDiameter, true(size(frame1)), [],[1 1 1 1 1]);
else
    nextNucleiXY = nucleiFrame2;
    intensity = zeros(size(nextNucleiXY,1),1);
end
inverse_mapping = zeros(size(nextNucleiXY,1),1);

%% 3. Update the position to correct for the shift

currentNucleiShiftedXY = zeros(size(nucleiFrame1));
nNuclei1 = size(nucleiFrame1, 1);

for j = 1:nNuclei1
    
    nucleusPosition = nucleiFrame1(j, :);
    %round and bracket the positions so they can be...
    %used as array indices
    xSub = min(max(round(nucleusPosition(2)), 1), xDim);
    ySub = min(max(round(nucleusPosition(1)), 1), yDim);
    
    currentNucleiShiftedXY(j,:) = nucleusPosition +...
        ...
        [vx( ySub, xSub) , vy( ySub, xSub)];
    
end


%% 4. Find the distance between each nuclei pair and build correspondences.
% Compute the distances between all nuclei
remainingNuclei1 = nucleiFrame1;
remainingNuclei2 = nextNucleiXY;
shiftedXY1 = currentNucleiShiftedXY;

indices1 = 1:size(remainingNuclei1,1);
indices2 = 1:size(remainingNuclei2,1);

mappedNuc1 = [];
mappedNuc2 = [];


% manualMapping = [];

%% 5. If manual data is provided, set them first.
if exist('manualMapping','var') && ~isempty(manualMapping)
    ind0 = find(manualMapping < 0); % Values below 0 enforce nuclei not to be mapped to anything on the next frame.
    ind1 = find(manualMapping > 0); % Values above 0 enforce nuclei to be mapped to a certain nucleus on the next frame.
    ind2 = manualMapping(manualMapping > 0); % Values to which nuclei ind1 were mapped to.
    
    
    mapping(ind1) = ind2;
    inverse_mapping(ind2) = ind1;
    
    mappedNuc1 = [mappedNuc1; remainingNuclei1(ind0,:)];
    mappedNuc1 = [mappedNuc1; remainingNuclei1(ind1,:)];
    
    mappedNuc2 = [mappedNuc2; remainingNuclei2(ind2,:)];
    indices1(ind0) = [];
    indices1(ind1) = [];
    indices2(ind2) = [];
    
    remainingNuclei1(ind0,:) = [];
    remainingNuclei1(ind1,:) = [];
    remainingNuclei2(ind2,:) = [];
    shiftedXY1 = remainingNuclei1;
    
end

keep_looping = true;
while keep_looping
    [dist] = pdist2(shiftedXY1,remainingNuclei2,'euclidean');
    
    % Build correspondences
    
    % First attribute nuclei that are mutually the closest
    [mc,ic] = min(dist,[],1);
    [mr,ir] = min(dist,[],2);
    
    indc = sub2ind(size(dist),ic,1:numel(ic));
    indr = sub2ind(size(dist),1:numel(ir),ir');
    
    indc(mc>maxNucleusStep) = [];
    indr(mr>maxNucleusStep) = [];
    
    Mc = zeros(size(dist));
    Mc(indc) = 1;
    Mr = zeros(size(dist));
    Mr(indr) = 1;
    
    attributedIndices = false(length(remainingNuclei1),1);
    ind = find(Mc & Mr);
    [rowInd,colInd] = ind2sub(size(dist),ind);
    for j = 1:numel(rowInd)
        mapping(indices1(rowInd(j))) = indices2(colInd(j));
        inverse_mapping(indices2(colInd(j))) = indices1(rowInd(j));
        attributedIndices(rowInd(j)) = true;
    end
    mappedNuc1 = [mappedNuc1; remainingNuclei1(rowInd,:)];
    mappedNuc2 = [mappedNuc2; remainingNuclei2(colInd,:)];
    
    indices1(rowInd) = [];
    indices2(colInd) = [];
    
    remainingNuclei1(rowInd,:) = [];
    remainingNuclei2(colInd,:) = [];
    
    try % If the number of mapped nuclei is not sufficient tpaps will fail.
        warning('off')
        st = tpaps(mappedNuc1', mappedNuc2',0.85);
        shiftedXY1 = fnval(st,remainingNuclei1')';
        warning('on')
    catch
        shiftedXY1 = remainingNuclei1;
    end
    
    
    if isempty(remainingNuclei1) || isempty(remainingNuclei2) || ~any(attributedIndices)
        keep_looping = false;
    end
    
end % while loop

if nargout > 2
    
    score = nan(size(nucleiFrame1,1),1);
    %projectedXY = fnval(st,nucleiFrame1');
    %dist = sqrt(sum( (projectedXY(:,mapping~=0)'-nextNucleiXY(mapping(mapping~=0),:)),2));
    
    %score(mapping~=0) = mvnpdf(dist , 0, maxNucleusStep).*intensity(mapping(mapping~=0));
    
    if nargout > 4
        
        %         naive_score = nan(size(nucleiFrame1,1),1);
        %         projectedXY = fnval(st,nucleiFrame1');
        %         dist = sqrt(sum( (projectedXY(mapping~=0,:)-nextNucleiXY(mapping(mapping~=0),:)),2));
        %         naive_score(mapping~=0) = mvnpdf(dist , 0, maxNucleusStep).*intensity(mapping(mapping~=0));
        
    end
end

end % function

