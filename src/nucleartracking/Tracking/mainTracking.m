function [ nuclei, varargout ] = mainTracking(FrameInfo, hisMovie, varargin)
%MAINTRACKING Run the full segmentation and tracking of a movie.
% Optionally, call again with more arguments to enforce manual
% corrections.
%
%
%   The function can be used in two ways :
%
%       1. [nuclei, centers] = mainTracking(names, 'indMitosis', indices_mitosis, 'embryoMask', embryo_mask)
%
%       Inputs :
%       - names: is a cell array containing the name of every frame in
%       order.
%       - indices_mitosis : nx2 array with each row defining the beginning
%       and ending of every mitosis. More precisely, if indices_mitosis =
%       [5   9
%        18 23]
%       two mitosis exist and during the first of the two, the first nuclei
%       divide between frame 5 and 6 and the last ones divide between 8 and
%       9. If not provided, the program will try to detect them on its own.
%       - embryoMask : logical array of the same size as the images that
%       is 'true' for values inside the embryo and 'false' outside. When no
%       edge of the embryo is visible, the mask can be set with
%       'true(size(one_image_from_the_movie))'. Everything outside the
%       embryo is simply ignored. If not provided,
%       the program will try to detect it on its own.
%
%       Outputs :
%       -nuclei : a structure containing the results from the tracking. The
%       structure contains the following fields, where n is the number of
%       frames in the movie :
%           - position : an nx2 array defining the position of the nucleus
%           on every frame. Rows are filled with NaNs if the nucleus is not
%           present on a frame.
%           - indXY : an nx1 array defining the indices of the nucleus in
%           the 'centers' output (see output n°2).
%           - P : index of the parent of this nucleus in the nuclei
%           structure.
%           - D, E : indices of the daughters for this nucleus in the nuclei
%           structure.
%
%
%        2. [nuclei] = mainTracking(names, 'centers', centers, 'mapping', map,...)
%
%        This usage is the same as before, but instead of segmenting the
%        nuclei from the images, the centers array is used. When correcting
%        tracking, the 'mapping' argument defining the matches to enforce
%        is provided.
%
%       'centers' has to be an nx1 cell array with each cell being an mx2
%       array with the x-y coordinates of the m nuclei present on the
%       corresponding frame.
%       'mapping' has to be an (n-1)x1 cell array with each cell being either
%       empty, if nothing has to be changed, or an mx2 array if a
%       correction has to be done. E.g. if map{2}(5,:) = [6 9], this means
%       that on frame 2, the nucleus with position 'centers(5,:)' is
%       dividing and its daughters are the nuclei with positions
%       'centers(6,:)' and 'centers(9,:)' on frame 3. If
%       map{2}(5,:) = [6 -1], nucleus 5 on frame 2 is forced to be mapped
%       nucleus 6 on frame 3 only. The -1 impairs the program to assume a
%       division happened and map it to a second nucleus. On the other hand
%       map{2}(5,:) = [6 0] means that the nucleus number 6 on frame 3
%       corresponds to the nucleus number 5 on frame 2 and, because of the
%       0, if these frames are within a division, the program could decide
%       that the nucleus divided and assign a second nucleus to the
%       mapping.
%
%
%       All the parameters used by the different functions are defined in
%       the function 'getDefaultParameters'. Typically, the pixel size,
%       nuclei diameters at each nuclear cycle, etc. are defined there.
%
%



h_waitbar_initialization = waitbar(0,'Initializing data...');

%%  0. Parse inputs
global provided_time_resolution
global provided_space_resolution
global provided_LoGratio
provided_time_resolution = NaN;
provided_space_resolution = NaN;
provided_LoGratio = NaN;
segmentationOnly = false;

% Edited to include multithresh by GM on 1/9/19
useMultithresh=false;
% Edited to include multithresh by GM on 1/9/19

% default values
ExpandedSpaceTolerance = 1;
NoBulkShift = 0;
 
% Edited to include multithresh by GM on 1/9/19
for j=1:2:numel(varargin)
    switch lower(varargin{j})
        case {'expandedspacetolerance'}
            ExpandedSpaceTolerance = varargin{j + 1};
        case {'nobulkshift'}
            NoBulkShift = varargin{j + 1};
        case {'indmitosis'}
            indMitosis = varargin{j+1};
        case {'centers'}
            centers = varargin{j+1};
        case {'mapping'}
            mapping = varargin{j+1};
        case {'embryomask'}
            embryoMask = varargin{j+1};
        case {'shifts'}
            shifts = varargin{j+1};
        case {'approved'}
            approvedCenters = varargin{j+1};
        case {'interpolated shifts'}
            interpolatedShifts = varargin{j+1};
        case {'space resolution'}
            provided_space_resolution = varargin{j+1};
        case {'time resolution'}
            provided_time_resolution = varargin{j+1};
        case {'logratio'}
            provided_LoGratio = varargin{j+1};
        case {'usemultithresh'}
            useMultithresh = varargin{j+1}; 
        case {'segmentationonly', 'segmentation only'}
            if ~islogical(segmentationOnly)
                segmentationOnly = false;
                warning('Unable to interpret the value of the ''Segmentation Only'' parameter. Value should be a logical. The value given was ignored.')
            else
                segmentationOnly = varargin{j+1};
            end
        case {'data structure', 'datastructure', 'data'}
            data = varargin{j+1};
            try
                provided_time_resolution = data.time_resolution;
            catch
                error('Unable to load the variable ''time_resolution'' from the provided data structure.')
            end
            try
                provided_space_resolution = data.space_resolution;
            catch
                error('Unable to load the variable ''space_resolution'' from the provided data structure.')
            end
            try
                interpolatedShifts = data.interpolatedShifts;
            catch
                error('Unable to load the variable ''interpolatedShifts'' from the provided data structure.')
            end
            try
                indMitosis = data.indMitosis;
            catch
                error('Unable to load the variable ''indMitosis'' from the provided data structure.')
            end
            try
                embryoMask = data.embryoMask;
            catch
                error('Unable to load the variable ''embryoMask'' from the provided data structure.')
            end
            try
                shifts = data.shifts;
            catch
                error('Unable to load the variable ''shifts'' from the provided data structure.')
            end
    end
end
% Edited to include multithresh by GM on 1/9/19


if isnan(provided_space_resolution)
    clear -global provided_space_resolution
end
if isnan(provided_time_resolution)
    clear -global provided_time_resolution
end
if isnan(provided_LoGratio)
    clear -global provided_LoGratio
end

%%


nFrames = size(hisMovie, 3);
originalNFrames = nFrames;

%If we don't have nc14 we'll fool the code into thinking that the last
    %frame of the movie was nc14
if isnan(indMitosis(end,1))
    
    indMitosis(end,1)= nFrames-2;
    indMitosis(end,2)= nFrames-1;
    
%     if nFrames - indMitosis(end, end) < 4
%         for k = 1:5
%             hisMovie(:, :, end+1) = zeros(size(hisMovie, 1), size(hisMovie, 2));
%             if exist('centers', 'var')
%                 centers(end+1) = centers(nFrames);
%             end
%         end
%         
%         indMitosis(end,1) = indMitosis(end,1) + 5;
%         indMitosis(end,2) = indMitosis(end,2) + 5;

%         nFrames = size(hisMovie, 3);
%     end 
%     
%     if exist('centers', 'var') && size(centers, 1) > nFrames
%         while size(centers, 1) > nFrames
%            centers(end) = []; 
%         end
%     end
%        
end
%%

waitbar(0.2,h_waitbar_initialization);

waitbar(0.3,h_waitbar_initialization)

% Mapping can't be provided without centers, otherwise ambiguous solution
% may arise.
if (~exist('centers','var') || isempty(centers) ) && exist('mapping','var')
    error('Provide centers that correspond to the provided mapping.')
end

if ~exist('mapping','var') || isempty(mapping)
    mapping = cell(nFrames-1,1);
end

if ~exist('interpolatedShifts','var') || isempty(interpolatedShifts)
    interpolatedShifts = cell(size(mapping));
end

% Get default parameters
time_resolution = getDefaultParameters(FrameInfo,'time resolution');
space_resolution = getDefaultParameters(FrameInfo,'space resolution');



%% 1. Start the segmentation if necessary.

waitbar(0.35,h_waitbar_initialization);

% Measure the average shift of the nuclei from one frame to the next.
if ~exist('shifts','var') || isempty('shifts')
    shifts = measureAllShifts(hisMovie,h_waitbar_initialization);
end

waitbar(0.8,h_waitbar_initialization)

% If not provided, find the mitosis.
if ~exist('indMitosis','var') || isempty('indMitosis')
    % Indices of putative mitosis:
    %     [indMitosis, dummy] = findMitosis(names,shifts);
    % Reshape in a correctly sized structure:
    tmp = zeros(numel(indMitosis,2));
    for j = 1:numel(indMitosis)
        tmp(j,1) = max(1,indMitosis(j)-4);
        tmp(j,2) = min(nFrames,indMitosis(j)+4);
    end
    indMitosis = tmp;
    indMitosis = [0 1; indMitosis; nFrames 0];
else
    % If the movie doesn't start within a mitosis, the first row has to be
    % [0,1] for the program to know it is an interphase.
    if indMitosis(1,1) > 1
        indMitosis = [0 1; indMitosis];
    end
    % Similarly, if the movie doesn't end in a mitosis, the last row has to
    % be [numberOfFrames 0] for the program to know it is an interphase.
    if all(indMitosis(end,:) < nFrames)
        indMitosis = [indMitosis; nFrames 0];
    end
end

waitbar(0.9,h_waitbar_initialization);

% Define the starting points for the tracking (at the middle of the
% interphase).
trackingStartingPoints = choseFramesToStartTracking(FrameInfo,indMitosis,nFrames);

% Determine the nuclear cycle of the interphases:
nucCyc = zeros(numel(trackingStartingPoints),1);
for j = 1:numel(trackingStartingPoints)
    nucCyc(end+1-j) = max(14-j+1,10);
end
% Associate the diameter of the cells at those nuclear cycles.
for j = 1:numel(trackingStartingPoints)
    diameters(j) = getDefaultParameters(FrameInfo,['d' num2str(nucCyc(j))]);
end
if ~exist('diameters','var') || isempty(diameters)
    % diameters is empty if there is no interphase in the movie (so all the
    %frames are from a single interphase). Assume the nucleus diameter
    %is between nc13 and nc14.
    diameters = 0.5*(getDefaultParameters(FrameInfo,'d13')+...
        getDefaultParameters(FrameInfo,'d14'));
end

close(h_waitbar_initialization);

if ~exist('centers','var') || isempty(centers)
    h_waitbar_segmentation = waitbar(0,'Segmentation...');
    xy = cell(nFrames,1);
    mitosisStartingPoint = [];
    if useMultithresh
        embryoMask = getMultithreshEmbryoMask(FrameInfo, names, diameters);
    end
    for j = 1:(size(indMitosis,1)-1)
        
        % Segment interphases
        
        segment = true;
        if size(indMitosis,1) == 1 % there is only one phase (mitosis or interphase) in the movie
            if indMitosis(1,1) ~= 1 % the only phase is an interphase
                first = 1;
                last = nFrames;
            else
                segment = false;
            end
        else
            first = max(indMitosis(j,2),1);
            if ~isempty(xy{first})
                first = first+1;
            end
            last = min(indMitosis(j+1,1),nFrames);
            if ~isempty(xy{last})
                last = last-1;
            end
        end
        if segment
            if ~useMultithresh
                xy(first:last) = segmentFrames(FrameInfo, hisMovie,first,last,diameters(j),embryoMask,h_waitbar_segmentation);
            else
                xy(first:last) = segmentFrames(FrameInfo, hisMovie,first,last,diameters(j),embryoMask,h_waitbar_segmentation, 'useMultithresh');
            end
        end
        
        % Segment mitosis
        
        segment = true;
        if all(indMitosis(j,:) ~= 0)
            if size(indMitosis,1) == 1 % there is only one phase (mitosis or interphase) in the movie
                if indMitosis(1,1) == 1 % the only phase is a mitosis
                    first = indMitosis(1,1);
                    last = indMitosis(1,2);
                else % the only phase is an interphase and was already segmented
                    segment = false;
                end
            else
                first = max(indMitosis(j,1),1);
                if ~isempty(xy{first})
                    first = first+1;
                end
                last = min(indMitosis(j,2),nFrames);
                if ~isempty(xy{last})
                    last = last-1;
                end
            end
            mitosisStartingPoint(end+1) = first;
            if j == 1 || size(indMitosis,1)
                D = diameters(j);
            else
                D = 0.5*sum(diameters(j-1:j));
            end
            
            if ~useMultithresh
                xy(first:last) = segmentFrames(FrameInfo,hisMovie,first,last,D,embryoMask,h_waitbar_segmentation);
            else
                xy(first:last) = segmentFrames(FrameInfo,hisMovie,first,last,D,embryoMask,h_waitbar_segmentation, 'useMultithresh');
            end

        end
        
    end
    
    
    close(h_waitbar_segmentation)
    

    %If the xy contains only one or zero nuclei then there's probably something
    %wrong. In that case just copy the information from the previous good
    %frame.
    
    %Edit AR 7/28/2020- this is commented out because it wreaks havoc
    %downstream in the tracking code for unknown reasons.
%     xy = fillEmptyXYFrames(xy);
    
else
    xy = centers;
end








% Initialize output
numberOfNuclei = size(xy{1},1);

% initialize array
nuclei = struct('position',nan(nFrames,2),'indXY',...
    mat2cell([1:numberOfNuclei; zeros(nFrames-1,numberOfNuclei)],...
    nFrames,ones(numberOfNuclei,1)),'P',[],'D',[],'E',[],'approved',0);

for j = 1:numel(nuclei)
    nuclei(j).position(1,:) = xy{1}(j,:);
end

if segmentationOnly
    varargout{1} = xy;
    return;
end



%% 2. Start the tracking

h_waitbar_tracking = waitbar(0,'Tracking...');

firstFrameIsAnInterphase = indMitosis(1,1) == 0;
numbMitosis = sum(all(indMitosis ~= 0,2));
numbInterphase = numbMitosis;
if any(indMitosis(1,:)==0) && any(indMitosis(end,:)==0)
    numbInterphase = numbInterphase + 1;
else if ~any(indMitosis(1,:)==0) && ~any(indMitosis(end,:)==0)
        numbInterphase = numbInterphase - 1;
    end
end
numberOfPhases = numbMitosis + numbInterphase;
if firstFrameIsAnInterphase
    phaseIsAMitosis = rem(0:(numberOfPhases-1),2);
else
    phaseIsAMitosis = rem(1:numberOfPhases,2);
end
if firstFrameIsAnInterphase
    offset = 1;
else
    offset = 0;
end
%Define the frames in each phase.
for j = 1:numberOfPhases
    %Gives the index of the most recent mitosis
    index = round((j+offset)/2);
    if phaseIsAMitosis(j)
        %Frames are from the start to finish of mitosis
        breakUpsFrames(j,:) = indMitosis(index,:);
    else
        %Frames are from the end of previous mitosis to the start of the next
        %mitosis.
        breakUpsFrames(j,:) = [indMitosis(index,2) indMitosis(index+1,1)];
    end
end
indInterphases = find(~phaseIsAMitosis);
diameters = zeros(numberOfPhases,1);

if phaseIsAMitosis(end)
    nuclearCycle = 13;
    diameters(end) = 0.5*(getDefaultParameters(FrameInfo,'d13')+getDefaultParameters(FrameInfo,'d14'));
else
    nuclearCycle = 14;
end
for j = numel(indInterphases):-1:1
    
    diameters(indInterphases(j)) = getDefaultParameters(FrameInfo,['d' num2str(max(nuclearCycle,10))]);
    previousNuclearCycle = nuclearCycle;
    nucCyc(j) = nuclearCycle;
    nuclearCycleFloored = max(nuclearCycle-1,10);
    if indInterphases(j) > 1
        diameters(indInterphases(j)-1) = 0.5 *...
            (getDefaultParameters(FrameInfo,['d' num2str(nuclearCycleFloored)])...
            + getDefaultParameters(FrameInfo, ['d' num2str(previousNuclearCycle)])) ;
    end
    nuclearCycle = nuclearCycle-1;
    
end
for j = 1:numberOfPhases
    %Sets the range of this phase based on the break up frames from above
    %(again checking for index out of bounds errors in advance)
    
    first = max(breakUpsFrames(j,1),1);
    last = min(breakUpsFrames(j,2),nFrames);
    
    if phaseIsAMitosis(j)
        if firstFrameIsAnInterphase
            fprintf(['Processing mitosis between nuclear cycle ' num2str(nucCyc(0.5*j)) ' and ' num2str(nucCyc(0.5*j)+1) '... ']);
        else
            try
                fprintf(['Processing mitosis between nuclear cycle ' num2str(nucCyc(0.5*(j+1))-1) ' and ' num2str(nucCyc(0.5*(j+1))) '... ']);
            catch
                warning('weirdness with how we specify nc13-only datasets in MovieDatabase')
            end
        end
        
%         try
            [ xy(first:last), mapping(first:last-1), nuclei ] =...
                trackMitosis(FrameInfo, hisMovie, first, last, shifts,...
                diameters(j), embryoMask, xy(first:last),...
                mapping(first:last-1), nuclei, h_waitbar_tracking );
%         catch
%             warning('couldn''t track this mitosis')
%         end
        
        fprintf('Done!\n')
        
        
    else
        
        
        if firstFrameIsAnInterphase
            fprintf(['Processing nuclear cycle ' num2str(nucCyc(0.5*(j+1))) '... ']);
        else
            fprintf(['Processing nuclear cycle ' num2str(nucCyc(0.5*j)) '... ']);
        end
        
%         try
            [nuclei, ~, interpolatedShifts] = trackWholeInterphase(FrameInfo,hisMovie,...
                trackingStartingPoints(1),first,last,diameters(j), embryoMask, ...
                xy, mapping,nuclei, interpolatedShifts, h_waitbar_tracking, ...
                ExpandedSpaceTolerance, NoBulkShift);
%         catch
%             warning('skipping this interphase');
%         end
        
        fprintf('Done!\n')
        
        
    end
end

close(h_waitbar_tracking)

if nargout > 1
    varargout{1} = xy;
    
    if nargout > 2
        if ~exist('approvedNuclei','var') || isempty(approvedNuclei)
            approvedNuclei = false(numel(nuclei),1);
        else
            approvedNuclei = false(numel(nuclei),1);
            for j = 1:numel(nuclei)
                indXY = nuclei(j).indXY;
                flags = false(sum(indXY>0),1);
                frames = find(indXY>0);
                indXY = indXY(indXY>0);
                for jj = 1:numel(flags)
                    if approvedCenters{frames(jj)}(indXY(jj))
                        flags(jj) = true;
                    end
                end
                if all(flags)
                    approvedNuclei(j) = true;
                else if ~any(flags)
                        approvedNuclei(j) = false;
                    else
                        'error'
                    end
                end
            end
            
        end
        
        if nargout > 3
            
            data.names = '';
            data.time_resolution = getDefaultParameters(FrameInfo,'time resolution');
            data.space_resolution = getDefaultParameters(FrameInfo,'space resolution');
            data.interpolatedShifts = interpolatedShifts;
            data.shifts = shifts;
            data.indMitosis = indMitosis;
            data.embryoMask = embryoMask;
            varargout{3} = data;
            
        end
    end
    
end


end
