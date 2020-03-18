function TrackNuclei(Prefix,varargin)
% TrackNuclei(Prefix, [Options])
%
% DESCRIPTION
% Use exported and projected nuclear images to segment and track nuclei over time.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% [Options]: See below.
%
% OPTIONS
% 'ExpandedSpaceTolerance': A multiplier for how how far away two nuclei of
% adjacent frames can be in order for them to be the same nuclei. It's
% recommended to set this to 1.5 if you use NoBulkShift.
% 'bulkShift': Runs the nuclear tracking with accounting for the bulk
% shift between frames (greatly reduces runtime).
% 'retrack': retrack
% 'integrate': integrate nuclear fluorescence
% 'segmentBetter': segment the nuclei well.
%
% OUTPUT
% '*_lin.mat' : Nuclei with lineages
% 'Ellipses.mat' : Just nuclei
%
% Author (contact): Hernan Garcia (hggarcia@berkeley.edu)
% Created: 01/01/2013 ish.
% Last Updated: 1/23/2019
%
% Documented by: Armando Reimer (areimer@berkeley.edu)
%
%


cleanupObj = onCleanup(@myCleanupFun);


disp(['Tracking nuclei on ', Prefix, '...']);




[stitchSchnitz, ExpandedSpaceTolerance, NoBulkShift,...
    retrack, nWorkers, track, noBreak, noStitch,...
    markandfind, fish, intFlag, chooseHis, segmentBetter]...
    = DetermineTrackNucleiOptions(varargin{:});


if segmentBetter
    if ~retrack
        resegmentAllFrames(Prefix);
    end
    retrack = true;
end

thisExperiment = liveExperiment(Prefix);

FrameInfo = getFrameInfo(thisExperiment);

[~, ProcPath, DropboxFolder, ~, PreProcPath] =...
    DetermineLocalFolders(Prefix);

anaphaseFrames = thisExperiment.anaphaseFrames';
nc9 = anaphaseFrames(1);
nc10 = anaphaseFrames(2);
nc11 = anaphaseFrames(3);
nc12 = anaphaseFrames(4);
nc13 = anaphaseFrames(5);
nc14 = anaphaseFrames(6);

if ~exist('nc9','var')
    error('Cannot find nuclear cycle values. Were they defined in MovieDatabase?')
end

%Do we need to convert any NaN chars into doubles?
if strcmpi(nc14,'nan')
    nc14=nan;
end
if strcmpi(nc13,'nan')
    nc13=nan;
end
if strcmpi(nc12,'nan')
    nc12=nan;
end
if strcmpi(nc11,'nan')
    nc11=nan;
end
if strcmpi(nc10,'nan')
    nc10=nan;
end
if strcmpi(nc9,'nan')
    nc9=nan;
end

%This checks whether all ncs have been defined
if length(anaphaseFrames)~=6
    error('Check the nc frames in the MovieDatabase entry. Some might be missing')
end

if (length(find(isnan(anaphaseFrames)))==length(anaphaseFrames))||(length(anaphaseFrames)<6)
    error('Have the ncs been defined in MovieDatabase?')
end

%Now do the nuclear segmentation and lineage tracking. This should be put
%into an independent function.

if chooseHis
    
    [hisFile, hisPath] = uigetfile([ProcPath, filesep, Prefix,'_',filesep,'*.mat']);
    hisStruct = load([hisPath, hisFile]);
    hisField = fieldnames(hisStruct);
    hisMat = hisStruct.(hisField);
    
else
    
    hisMat =  loadHisMat([PreProcPath, filesep, Prefix, filesep, Prefix, '_hisMat.mat']);
    
end


%Pull the mitosis information from ncs.
anaphaseFrames=anaphaseFrames(anaphaseFrames~=0);

%Note that I'm adding a range of two frames frames before and after the
%determines mitosis
indMit=[anaphaseFrames'-2,anaphaseFrames'+2];

%Make sure no indices are negative. This could happen is the nuclear cycle
%started at frame 1, for example.
indMit(indMit<1)=1;

%Check whether nc14 occurred very close to the end of the movie. For those
%frames we'll move the boundary for tracking purposes
nFrames = length(FrameInfo);
indMit(indMit>=nFrames)=indMit(indMit>=nFrames)-2;

%If indMit is empty this means that we didn't see any mitosis. Deal with
%this by assigning the first frames to it
if isempty(indMit)
    indMit=[1,2];
    %If we don't have nc14 we'll fool the code into thinking that the last
    %frame of the movie was nc14
elseif isnan(indMit(end,1))
    indMit(end,1)=nFrames-6;
    indMit(end,2)=nFrames-5;
end


%Embryo mask
ImageTemp=squeeze(hisMat(:, :, 1));
embryo_mask=true(size(ImageTemp));
clear ImageTemp



%Get information about the spatial and temporal resolution
settingArguments{1}='time resolution';
settingArguments{2}=median(diff([FrameInfo.Time]));     %Median separation between frames (in seconds)
settingArguments{3}='space resolution';
settingArguments{4}=FrameInfo(1).PixelSize;

schnitzFileExists = exist([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'file');

if schnitzFileExists & ~retrack &track
    answer=input('Previous tracking detected. Proceed with retracking? (y/n):','s');
    if strcmpi(answer,'y')
        retrack = true;
    elseif strcmpi(answer, 'n')
        retrack = false;
    else
        %do nothing
    end
end



%Do the tracking for the first time
if ~retrack
    
    
    if track
        [nuclei, centers, ~, dataStructure] = ...
            mainTracking(FrameInfo, hisMat,'indMitosis',indMit,'embryoMask', embryo_mask,...
            settingArguments{:}, 'ExpandedSpaceTolerance', ExpandedSpaceTolerance, ...
            'NoBulkShift', NoBulkShift);
    else
        [nuclei, centers] = ...
            mainTracking(FrameInfo, hisMat,'indMitosis',indMit,'embryoMask', embryo_mask,...
            settingArguments{:}, 'ExpandedSpaceTolerance', ExpandedSpaceTolerance, ...
            'NoBulkShift', NoBulkShift, 'segmentationonly', true);
    end
    % names is a cell array containing the names of all frames in the movie in order.
    % indMitosis is an nx2 array containing the first and last frame of mitosis in every row.
    % embryoMask is the possible mask of the embryo. If no embryo edge is visible,
    % true(size(an_image_from_the_movie)) can be given as input.
    % Convert the results to compatible structures and save them
    %Put circles on the nuclei
    [Ellipses] = putCirclesOnNuclei(FrameInfo,centers,nFrames,indMit);
    %Convert nuclei structure into schnitzcell structure
    [schnitzcells] = convertNucleiToSchnitzcells(nuclei);
    
else
    %Do re-tracking
    disp 'Re-tracking...'
    warning('Re-tracking.') %This code needs to be able to handle approved schnitz still
    
    %Load the Ellipses and re-generate the centers
    load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses')
    centers = updateCentersFromEllipses(FrameInfo, Ellipses);
    
    %Load the dataStructure to speed up retracking if it exists
    if exist([ProcPath,filesep,Prefix,'_',filesep,'dataStructure.mat'], 'file')
        load([ProcPath,filesep,Prefix,'_',filesep,'dataStructure.mat'])
    elseif exist([ProcPath,filesep,Prefix,'_',filesep,'TrackingDataStructure.mat'], 'file')
        load([ProcPath,filesep,Prefix,'_',filesep,'TrackingDataStructure.mat'])
    end
    
    % look for a settings file in the Raw Data folder.
    if exist([PreProcPath,filesep,Prefix,filesep,Prefix,'-AcquisitionSettings.mat'],'file')
        load([PreProcPath,filesep,Prefix,filesep,Prefix,'-AcquisitionSettings.mat']);
        fields = fieldnames(AcquisitionSettings);
        settingArguments = cell(1,2*length(fields));
        for i=1:length(fields)
            settingArguments{2*i-1} = fields{i}; % field names must match the parameter names in mainTracking
            settingArguments{2*i} = getfield(AcquisitionSettings, fields{i});
            % alter argument names if necessary
            if strcmp(settingArguments{2*i-1}, 'time_resolution')
                settingArguments{2*i-1} = 'time resolution';
            elseif strcmp(settingArguments{2*i-1}, 'space_resolution')
                settingArguments{2*i-1} = 'space resolution';
            end
        end
    else
        settingArguments = {};
    end
    
    
    %Re-run the tracking
    if exist('dataStructure', 'var')
        %Edit the names in dataStructure to match the current folder setup
        dataStructure.names='';
        
        [nuclei, centers, ~, dataStructure] = mainTracking(...
            FrameInfo, hisMat,'indMitosis',indMit,'embryoMask', embryo_mask,...
            'centers',centers,'dataStructure',dataStructure, settingArguments{:}, ...
            'ExpandedSpaceTolerance', ExpandedSpaceTolerance, ...
            'NoBulkShift', NoBulkShift);
    else
        [nuclei, centers, ~, dataStructure] = mainTracking(...
            FrameInfo, hisMat,'indMitosis',indMit,'embryoMask', embryo_mask,...
            'centers',centers, settingArguments{:}, 'ExpandedSpaceTolerance', ...
            ExpandedSpaceTolerance, 'NoBulkShift', NoBulkShift);
    end
    
    %Put circles on the nuclei
    %     [Ellipses] = putCirclesOnNuclei(FrameInfo,centers,nFrames,indMit);
    %Convert nuclei structure into schnitzcell structure
    [schnitzcells] = convertNucleiToSchnitzcells(nuclei);
end

%Add the radius information to the schnitz
for i=1:length(schnitzcells)
    for j=1:length(schnitzcells(i).frames)
        schnitzcells(i).len(:)=...
            single(mean(Ellipses{schnitzcells(i).frames(j)}(schnitzcells(i).cellno(j),3:4)));
    end
end

expandedAnaphaseFrames = [zeros(1,8),thisExperiment.anaphaseFrames'];

for s = 1:length(schnitzcells)
    midFrame = ceil(length(schnitzcells(s).frames)/2);
    dif = double(schnitzcells(s).frames(midFrame)) - expandedAnaphaseFrames;
    cycle = find(dif>0, 1, 'last' );
    schnitzcells(s).cycle = uint8(cycle);
end

schnitzcells = addRelativeTimeToSchnitzcells(schnitzcells, FrameInfo, expandedAnaphaseFrames);

%perform some quality control
schnitzcells = filterSchnitz(schnitzcells, [thisExperiment.yDim, thisExperiment.xDim]);

%Save everything at this point. It will be overwritten later, but it's
%useful for debugging purposes if there's a bug in the code below.
if ~exist([DropboxFolder,filesep,Prefix], 'dir')
    mkdir([DropboxFolder,filesep,Prefix]);
end

if whos(var2str(Ellipses)).bytes < 2E9
    
    save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses', '-v6');
else
    save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses', '-v7.3');
    
end

if whos(var2str(schnitzcells)).bytes < 2E9
    save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells', '-v6');
else
    save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells', '-v7.3');
end

if ~exist([ProcPath,filesep,Prefix,'_'], 'dir')
    mkdir([ProcPath,filesep,Prefix,'_']);
end
if exist('dataStructure', 'var')
    save([ProcPath,filesep,Prefix,'_',filesep,'dataStructure.mat'],'dataStructure');
end

%Extract the nuclear fluorescence values if we're in the right experiment
%type
if intFlag
    schnitzcells = integrateSchnitzFluo(Prefix, schnitzcells, FrameInfo, PreProcPath);
end

if fish
    schnitzcells = rmfield(schnitzcells, {'P', 'E', 'D'});
end


if track & ~noBreak
    [schnitzcells, Ellipses] = breakUpSchnitzesAtMitoses(schnitzcells, Ellipses, expandedAnaphaseFrames, nFrames);
    save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses');
    if (whos(var2str(schnitzcells)).bytes < 2E9)
        save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells', '-v6');
    else
        save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells', '-v7.3', '-nocompression');
    end
end

% Stitch the schnitzcells using Simon's code
if ~noStitch
    disp('stitching schnitzes')
    %     try
    StitchSchnitz(Prefix, nWorkers);
    %     catch
    %         disp('failed to stitch schnitzes')
    %     end
end



end