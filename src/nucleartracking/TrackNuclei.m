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
% 'mixedPolaritySegmentation': different segmentation method that works
% better when there are nuclei of mixed polarity (some dark, some bright)
% 'adjustNuclearContours': fit ellipses tightly around nuclei (done by
% adjusting regular segmentation)
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

postTrackingSettings = struct; 

[ExpandedSpaceTolerance, NoBulkShift,...
    retrack, nWorkers, track, postTrackingSettings.noBreak,...
    postTrackingSettings.noStitch,...
    markandfind, postTrackingSettings.fish,...
    postTrackingSettings.intFlag, chooseHis,...
    mixedPolaritySegmentation, min_rad_um,...
             max_rad_um, sigmaK_um, mu, nIterSnakes,...
             postTrackingSettings.doAdjustNuclearContours, radiusScale]...
    = DetermineTrackNucleiOptions(varargin{:});


postTrackingSettings.track = track;
postTrackingSettings.shouldConvertToAP = true;


liveExperiment = LiveExperiment(Prefix);

FrameInfo = getFrameInfo(liveExperiment);

ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;
PreProcPath = liveExperiment.preFolder;


ellipsesFile = [DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'];
schnitzcellsFile = [DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat']; 

anaphaseFrames = getAnaphaseFrames(liveExperiment); 

if iscolumn(anaphaseFrames)
    anaphaseFrames = anaphaseFrames';
end

%This checks whether all ncs have been defined
if length(anaphaseFrames)~=6
    error('Check the nc frames in the MovieDatabase entry. Some might be missing')
end

if length( find(isnan(anaphaseFrames))) ==...
        length(anaphaseFrames) || length(anaphaseFrames) < 6
    error('Have the ncs been defined in MovieDatabase or anaphaseFrames.mat?')
end

%Now do the nuclear segmentation and lineage tracking. This should be put
%into an independent function.

if chooseHis
    
    [hisFile, hisPath] = uigetfile([ProcPath, filesep, Prefix,'_',filesep,'*.*']);
    hisStruct = load([hisPath, hisFile]);
    hisField = fieldnames(hisStruct);
    hisMat = hisStruct.(hisField{1});
    
else
    
    hisMat =  getHisMat(liveExperiment);
    
end

if mixedPolaritySegmentation
    if ~retrack
        resegmentAllFrames(Prefix, 'hisMat', hisMat,...
            'min_rad_um', min_rad_um,...
            'max_rad_um', max_rad_um,'sigmaK_um',sigmaK_um,'mu', mu,...
            'nInterSnakes',nIterSnakes);
    end
    retrack = true;
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

expandedAnaphaseFrames = [zeros(1,8),liveExperiment.anaphaseFrames'];

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

if schnitzFileExists && ~retrack && track
    answer=input('Previous tracking detected. Erase existing segmentation? (y/n):','s');
    if strcmpi(answer,'n')
        retrack = true;
    elseif strcmpi(answer, 'y')
        retrack = false;
    end
end



%Perform both segmentation and tracking with the "MainTracking" script
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
    [Ellipses] = putCirclesOnNuclei(FrameInfo,centers,nFrames,indMit, radiusScale);
    
else
    %Retrack: Use "MainTracking" for tracking but not segmentation. 
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
    
    clear dataStructure;
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
    
end

disp('Finished main tracking.'); 

%Convert nuclei structure into schnitzcell structure
[schnitzcells] = convertNucleiToSchnitzcells(nuclei);




%Save everything at this point. It will be overwritten later, but it's
%useful for debugging purposes if there's a bug in the code below.
if ~exist([DropboxFolder,filesep,Prefix], 'dir')
    mkdir([DropboxFolder,filesep,Prefix]);
end


save2(ellipsesFile, Ellipses); 
save2(schnitzcellsFile, schnitzcells); 

if ~exist([ProcPath,filesep,Prefix,'_'], 'dir')
    mkdir([ProcPath,filesep,Prefix,'_']);
end
if exist('dataStructure', 'var')
    save([ProcPath,filesep,Prefix,'_',filesep,'dataStructure.mat'],'dataStructure');
end

%perform very important stuff subsequent to tracking proper
performPostNuclearTracking(Prefix,...
    expandedAnaphaseFrames, nWorkers, schnitzcellsFile,...
    ellipsesFile, postTrackingSettings)

end