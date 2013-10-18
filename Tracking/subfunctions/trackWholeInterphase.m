function [ nuclei, varargout ] = trackWholeInterphase( names, startingFrame, previousMitosisInd, nextMitosisInd, nucleusDiameter, embryoMask, xy, mapping, nuclei, shifts, varargin )
%TRACKWHOLEINTERPHASE Summary of this function goes here
%   Detailed explanation goes here

if numel(varargin) > 0
    h_waitbar_tracking = varargin{1};
end

startingFrame = previousMitosisInd;

totalNumberOfFrames = numel(names);
numberOfFrames = nextMitosisInd-previousMitosisInd+1;
time_resolution = getDefaultParameters('time resolution');
space_resolution = getDefaultParameters('space resolution');
edgeClearance = getDefaultParameters('edge clearance')*nucleusDiameter/space_resolution;
img = imread(names{startingFrame});
marginBeforeMitosis = ceil(getDefaultParameters('increased precision before mitosis')/time_resolution);
marginAfterMitosis = ceil(getDefaultParameters('increased precision after mitosis')/time_resolution);

if exist('xy','var') && ~isempty(xy)
    skip_segmentation = true;
    XY = xy(previousMitosisInd:nextMitosisInd);
else
    skip_segmentation = false;
    XY = cell(numberOfFrames,1);
    XY{startingFrame-previousMitosisInd+1} = findNuclei(names,startingFrame,nucleusDiameter, embryoMask, [],[1 1 1 1 1]);
end
if ~exist('shifts','var') || isempty(shifts)
    shifts = cell(totalNumberOfFrames-1,1);
end

numberOfNuclei = size(XY{startingFrame-previousMitosisInd+1},1);

% initialize array
% nuclei = struct('position',nan(totalNumberOfFrames,2),'indXY',mat2cell([zeros(startingFrame-1,numberOfNuclei); 1:numberOfNuclei; zeros(totalNumberOfFrames-startingFrame,numberOfNuclei)],totalNumberOfFrames,ones(numberOfNuclei,1)),'score',nan(totalNumberOfFrames,1),'P',[],'D',[],'E',[],'approved',0);
% 
% 
% for j = 1:numberOfNuclei
%     
%     nuclei(j).position(startingFrame,:) = XY{startingFrame-previousMitosisInd+1}(j,:);
%     
% end

ind = true(length(XY{startingFrame-previousMitosisInd+1}),1);
nucleiIndices = nan(size(XY{1},1),1);
for j = 1:size(XY{1},1);
    for jj = 1:numel(nuclei)
        if nuclei(jj).indXY(previousMitosisInd) == j
            nucleiIndices(j) = jj;
            break
        end
    end
end
if any(isnan(nucleiIndices))
    error(' NAN ')
end
    %nucleiIndices = 1:numel(nuclei); % Temporar vector that contains the mapping between 'XY' rows and 'nuclei' elements, i.e. XY{frame}(j,:) = nuclei(nucleiIndices(j)).position(frame,:);
targetNumber = numel(nuclei);
%% Track forwards
for j = 1:(nextMitosisInd-startingFrame)%[]%1:(nextMitosisInd-startingFrame)
    
    currentFrameNumber = startingFrame+j-1; % Frame that was just analyzed before.
    currentFrameInd = startingFrame-previousMitosisInd+j; % Index in the XY cell array.
    newFrameNumber = startingFrame+j; % Number of the new frame to analyze (corresponds to current frame number +1 because we're tracking forwards).
    newFrameInd = startingFrame-previousMitosisInd+j+1; % Index in the XY cell array.
    
    %Added by HG
    [currentFrameNumber,currentFrameInd,newFrameNumber,newFrameInd];
    
    if skip_segmentation
        [mapping{currentFrameInd},dummy1,dummy2,ind, shifts{currentFrameNumber}] = frame2frameCorrespondence(names,currentFrameNumber,newFrameNumber,XY{currentFrameInd},nucleusDiameter,1,XY{newFrameInd},mapping{currentFrameNumber},shifts{currentFrameNumber});%, embryoMask, targetNumber, [1 1 1]);
    else
        [mapping{currentFrameInd},XY{newFrameInd},dummy,ind, shifts{currentFrameNumber}] = frame2frameCorrespondence(names,currentFrameNumber,newFrameNumber,XY{currentFrameInd},nucleusDiameter,1,[],shifts{currentFrameNumber});%, embryoMask, [],[1 1 0]);
    end
    
    % Put the output in the nuclei structure.
    nucInd = zeros(numel(ind),1);
    for jj = 1:numel(mapping{currentFrameInd})
        if mapping{currentFrameInd}(jj) == 0
            continue;
        end
        
        nuclei(nucleiIndices(jj)).indXY(startingFrame+j) = mapping{currentFrameInd}(jj);
        nuclei(nucleiIndices(jj)).position(startingFrame+j,:) = XY{newFrameInd}(mapping{currentFrameInd}(jj),:);

        %nuclei(nucleiIndices(jj)).score(startingFrame+j) = score(jj);
        nucInd(mapping{currentFrameInd}(jj)) = nucleiIndices(jj);
    end
    
    % Add the nuclei that were found in the next frame but that weren't
    % mapped to any current nucleus.
    orphanNuclei = find(ind == 0);
    indToDelete = [];
    for jj = 1:numel(orphanNuclei)
        xpos = XY{newFrameInd}(orphanNuclei(jj),1);
        ypos = XY{newFrameInd}(orphanNuclei(jj),2);
        
        %if xpos >= 1+edgeClearance && xpos <= size(img,1)-edgeClearance && ypos >= 1+edgeClearance && ypos <= size(img,2)-edgeClearance
            
            IND = numel(nuclei)+1;
            nuclei(IND).position = nan(totalNumberOfFrames,2);
            nuclei(IND).indXY = zeros(totalNumberOfFrames,1);
            nuclei(IND).position(startingFrame+j,:) = XY{newFrameInd}(orphanNuclei(jj),:);
            nuclei(IND).indXY(startingFrame+j) = orphanNuclei(jj);
            %            nuclei(IND).score(startingFrame+j) = score(orphanNuclei(jj));
            nucInd(orphanNuclei(jj)) = IND;
%         else
%             indToDelete = [indToDelete jj];
%         end
    end
    XY{newFrameInd}(orphanNuclei(indToDelete),:) = [];
    nucInd(orphanNuclei(indToDelete)) = [];
    nucleiIndices = nucInd;
    
    
    waitbar((j+startingFrame-1)/(totalNumberOfFrames-1), h_waitbar_tracking, ['Tracking progress : processing frames ' num2str(j+startingFrame-1) ' and ' num2str(j+startingFrame) ' out of ' num2str(totalNumberOfFrames) '...']);
end


varargout{1} = XY;
varargout{2} = shifts;

end

