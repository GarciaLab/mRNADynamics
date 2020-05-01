function [ nuclei, varargout ] = trackWholeInterphase(FrameInfo, hisMat, startingFrame, ...
    previousMitosisInd, nextMitosisInd, nucleusDiameter, embryoMask, xy, ...
    mapping, nuclei, shifts, varargin )
%TRACKWHOLEINTERPHASE Summary of this function goes here
%   Detailed explanation goes here

if numel(varargin) > 0
    h_waitbar_tracking = varargin{1};
    if numel(varargin) > 1
        ExpandedSpaceTolerance = varargin{2};
        NoBulkShift = varargin{3};
    end
end

startingFrame = previousMitosisInd;

totalNumberOfFrames = size(hisMat, 3);
numberOfFrames = nextMitosisInd-previousMitosisInd+1;
time_resolution = getDefaultParameters(FrameInfo,'time resolution');
space_resolution = getDefaultParameters(FrameInfo,'space resolution');
edgeClearance = getDefaultParameters(FrameInfo,'edge clearance')*nucleusDiameter/space_resolution;
img = hisMat(:, :, startingFrame);
marginBeforeMitosis = ceil(getDefaultParameters(FrameInfo,'increased precision before mitosis')/time_resolution);
marginAfterMitosis = ceil(getDefaultParameters(FrameInfo,'increased precision after mitosis')/time_resolution);
%This looks like it decides to use existing data that was passed to it, and
%if none was provided makes it itself. 
if exist('xy','var') && ~isempty(xy)
    
    skip_segmentation = true;
    
    xyInterphase = xy(previousMitosisInd:nextMitosisInd);
    
else
    skip_segmentation = false;
    xyInterphase = cell(numberOfFrames,1);
    xyInterphase{startingFrame-previousMitosisInd+1} = findNuclei(FrameInfo,hisMat,startingFrame,...
        nucleusDiameter, embryoMask, [],[1 1 1 1 1]);
end

%the code down below will break if frames are empty. 
%let's fix that now. 
emptyFrames = find(cellfun(@isempty, xyInterphase));

% fh = @(x) fillmissing(x, 'previous');
for f = 1:length(emptyFrames)
    xyInterphase{emptyFrames(f)} = [0 0];
end


if ~exist('shifts','var') || isempty(shifts)
    shifts = cell(totalNumberOfFrames-1,1);
end
%this startingFrame variable appears to be a holdover from an older version of the
%code. This like invariably returns size(xyInterphase{1},1);
numberOfNuclei = size(xyInterphase{startingFrame-previousMitosisInd+1},1);


ind = true(length(xyInterphase{startingFrame-previousMitosisInd+1}),1);
nucleiIndices = nan(size(xyInterphase{1},1),1);
for j = 1:size(xyInterphase{1},1)
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


%% Track forwards
for j = 1:(nextMitosisInd-startingFrame)
    
    currentFrameNumber = startingFrame+j-1; % Frame that was just analyzed before.
    currentFrameInd = startingFrame-previousMitosisInd+j; % Index in the XY cell array.
    newFrameNumber = startingFrame+j; % Number of the new frame to analyze (corresponds to current frame number +1 because we're tracking forwards).
    newFrameInd = startingFrame-previousMitosisInd+j+1; % Index in the XY cell array.
        
    if skip_segmentation
        [mapping{currentFrameInd},~,~,ind, shifts{currentFrameNumber}] =...
            frame2frameCorrespondence(FrameInfo, hisMat, currentFrameNumber,...
            newFrameNumber,xyInterphase{currentFrameInd},nucleusDiameter,1,xyInterphase{newFrameInd},...
            mapping{currentFrameNumber},shifts{currentFrameNumber}, ...
            ExpandedSpaceTolerance, NoBulkShift);
    else
        [mapping{currentFrameInd},xyInterphase{newFrameInd},dummy,ind, shifts{currentFrameNumber}]...
            = frame2frameCorrespondence(FrameInfo,hisMat,currentFrameNumber,...
            newFrameNumber,xyInterphase{currentFrameInd},nucleusDiameter,1,[],...
            shifts{currentFrameNumber}, ExpandedSpaceTolerance, NoBulkShift);
    end
    
    % Put the output in the nuclei structure.
    nucInd = zeros(numel(ind),1);
    for jj = 1:numel(mapping{currentFrameInd})
        if mapping{currentFrameInd}(jj) == 0
            continue;
        end
        
        nuclei(nucleiIndices(jj) ).indXY(startingFrame+j) = mapping{currentFrameInd}(jj);
        nuclei(nucleiIndices(jj) ).position(startingFrame+j,:) = xyInterphase{newFrameInd}(mapping{currentFrameInd}(jj),:);

        %nuclei(nucleiIndices(jj)).score(startingFrame+j) = score(jj);
        nucInd(mapping{currentFrameInd}(jj)) = nucleiIndices(jj);
    end
    
    % Add the nuclei that were found in the next frame but that weren't
    % mapped to any current nucleus.
    orphanNuclei = find(ind == 0);
    indToDelete = [];
    for jj = 1:numel(orphanNuclei)
        xpos = xyInterphase{newFrameInd}(orphanNuclei(jj),1);
        ypos = xyInterphase{newFrameInd}(orphanNuclei(jj),2);
        
        %if xpos >= 1+edgeClearance && xpos <= size(img,1)-edgeClearance && ypos >= 1+edgeClearance && ypos <= size(img,2)-edgeClearance
            
            IND = numel(nuclei)+1;
            nuclei(IND).position = nan(totalNumberOfFrames,2);
            nuclei(IND).indXY = zeros(totalNumberOfFrames,1);
            nuclei(IND).position(startingFrame+j,:) = xyInterphase{newFrameInd}(orphanNuclei(jj),:);
            nuclei(IND).indXY(startingFrame+j) = orphanNuclei(jj);
            %            nuclei(IND).score(startingFrame+j) = score(orphanNuclei(jj));
            nucInd(orphanNuclei(jj)) = IND;
%         else
%             indToDelete = [indToDelete jj];
%         end
    end
    xyInterphase{newFrameInd}(orphanNuclei(indToDelete),:) = [];
    nucInd(orphanNuclei(indToDelete)) = [];
    nucleiIndices = nucInd;
    
    
    waitbar((j+startingFrame-1)/(totalNumberOfFrames-1), h_waitbar_tracking,...
        ['Tracking progress : processing frames ' num2str(j+startingFrame-1),...
        ' and ' num2str(j+startingFrame) ' out of ' num2str(totalNumberOfFrames) '...']);
    
end


varargout{1} = xyInterphase;
varargout{2} = shifts;

end

