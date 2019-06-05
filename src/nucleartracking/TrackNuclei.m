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
% 'StitchSchnitz' : Run the schnitzcells fixing code by Simon
% 'ExpandedSpaceTolerance': A multiplier for how how far away two nuclei of
% adjacent frames can be in order for them to be the same nuclei. It's
% recommended to set this to 1.5 if you use NoBulkShift. 
% 'NoBulkShift': Runs the nuclear tracking without accounting for the bulk
% shift between frames (greatly speeds up runtime). It's recommended you
% set ExpandedSpaceTolerance to 1.5 if you use this. 
% 'retrack': retrack
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

disp(['Tracking nuclei on ', Prefix, '...']);

[SkipStitchSchnitz, ExpandedSpaceTolerance, NoBulkShift, retrack] = DetermineTrackNucleiOptions(varargin{:});

%added this bulkshift requirement because i'm not sure it'll work without
%it
if NoBulkShift
    workers = 8;
    startParallelPool(workers, 0,0);
end



[~, ProcPath, DropboxFolder, ~, PreProcPath, ~, ~] = DetermineLocalFolders(Prefix);

CONFIG_CSV_PATH = 'ComputerFolders.csv';
configValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');    
DefaultDropboxFolder = getConfigValue(configValues, 'DropboxFolder');

%Determine division times
%Load the information about the nc from moviedatabase file
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

%If Channel2 was left empty, it would contain a NaN, which will cause
%problems below. In that case, replace it by an empty string.
if isnan(Channel2{1})
    Channel2{1}='';
end

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
ncCheck=[nc9,nc10,nc11,nc12,nc13,nc14];
if length(ncCheck)~=6
    error('Check the nc frames in the MovieDatabase entry. Some might be missing')
end

ncs=[nc9,nc10,nc11,nc12,nc13,nc14];

if (length(find(isnan(ncs)))==length(ncs))||(length(ncs)<6)
    error('Have the ncs been defined in MovieDatabase?')
end

%Now do the nuclear segmentation and lineage tracking. This should be put
%into an independent function.

%Create the cell array with the names.
D=dir([PreProcPath,filesep,Prefix,filesep,Prefix,'-His_*.tif']);
names = cell(length(D), 1);
for i=1:length(D)
    names{i}=[PreProcPath,filesep,Prefix,filesep,D(i).name];
end


%Pull the mitosis information from ncs.
ncs=ncs(ncs~=0);

%Note that I'm adding a range of two frames frames before and after the
%determines mitosis
indMit=[ncs'-2,ncs'+2];

%Make sure no indices are negative. This could happen is the nuclear cycle
%started at frame 1, for example.
indMit(indMit<1)=1;

%Check whether nc14 occurred very close to the end of the movie. For those
%frames we'll move the boundary for tracking purposes
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo')
numFrames = length(FrameInfo);
indMit(indMit>=numFrames)=indMit(indMit>=numFrames)-2;

%If indMit is empty this means that we didn't see any mitosis. Deal with
%this by assigning the first frames to it
if isempty(indMit)
    indMit=[1,2];
%If we don't have nc14 we'll fool the code into thinking that the last
%frame of the movie was nc14
elseif isnan(indMit(end,1))
    indMit(end,1)=numFrames-6;
    indMit(end,2)=numFrames-5;
end    


%Embryo mask
ImageTemp=imread(names{1});
embryo_mask=true(size(ImageTemp));
clear ImageTemp


  
%Get information about the spatial and temporal resolution
settingArguments{1}='time resolution';
settingArguments{2}=median(diff([FrameInfo.Time]));     %Median separation between frames (in seconds)
settingArguments{3}='space resolution';
settingArguments{4}=FrameInfo(1).PixelSize;

schnitzFileExists = exist([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'file');

if schnitzFileExists & ~retrack
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
    
    
    [nuclei, centers, ~, dataStructure] = ...
        mainTracking(FrameInfo, names,'indMitosis',indMit,'embryoMask', embryo_mask,...
        settingArguments{:}, 'ExpandedSpaceTolerance', ExpandedSpaceTolerance, ...
        'NoBulkShift', NoBulkShift);
    % names is a cell array containing the names of all frames in the movie in order.
    % indMitosis is an nx2 array containing the first and last frame of mitosis in every row.
    % embryoMask is the possible mask of the embryo. If no embryo edge is visible,
    % true(size(an_image_from_the_movie)) can be given as input.

    % Convert the results to compatible structures and save them
    %Put circles on the nuclei
    [Ellipses] = putCirclesOnNuclei(FrameInfo,centers,names,indMit);
    %Convert nuclei structure into schnitzcell structure
    [schnitzcells] = convertNucleiToSchnitzcells(nuclei); 
else
    %Do re-tracking
    disp 'Re-tracking...'
    warning('Re-tracking.') %This code needs to be able to handle approved schnitz still
  
    %Load the Ellipses and re-generate the centers
    load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses')
    %centers = updateCentersFromEllipses(Ellipses, centers);
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
        dataStructure.names=names;
        
        [nuclei, centers, ~, dataStructure] = mainTracking(...
            FrameInfo, names,'indMitosis',indMit,'embryoMask', embryo_mask,...
            'centers',centers,'dataStructure',dataStructure, settingArguments{:}, ...
            'ExpandedSpaceTolerance', ExpandedSpaceTolerance, ...
        'NoBulkShift', NoBulkShift);
    else
        [nuclei, centers, ~, dataStructure] = mainTracking(...
            FrameInfo, names,'indMitosis',indMit,'embryoMask', embryo_mask,...
            'centers',centers, settingArguments{:}, 'ExpandedSpaceTolerance', ...
            ExpandedSpaceTolerance, 'NoBulkShift', NoBulkShift);
    end

    %Put circles on the nuclei
    [Ellipses] = putCirclesOnNuclei(FrameInfo,centers,names,indMit);
    %Convert nuclei structure into schnitzcell structure
    [schnitzcells] = convertNucleiToSchnitzcells(nuclei); 
end

%Add the radius information to the schnitz
for i=1:length(schnitzcells)
    for j=1:length(schnitzcells(i).frames)
        schnitzcells(i).len(:)=...
            mean(Ellipses{schnitzcells(i).frames(j)}(schnitzcells(i).cellno(j),3:4));
    end
end


%Save everything at this point. It will be overwritten later, but it's
%useful for debugging purposes if there's a bug in the code below.
mkdir([DropboxFolder,filesep,Prefix])
save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses')
save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells')
if ~exist([ProcPath,filesep,Prefix,'_'], 'dir')
    mkdir([ProcPath,filesep,Prefix,'_'])
end
save([ProcPath,filesep,Prefix,'_',filesep,'dataStructure.mat'],'dataStructure')



%Extract the nuclear fluorescence values if we're in the right experiment
%type
if strcmpi(ExperimentType,'inputoutput')||strcmpi(ExperimentType,'input')
    
    
    %Parse the channel information for the different experiment types
    if strcmpi(ExperimentType,'inputoutput')
        InputChannelTemp=strfind({lower(Channel1{1}),lower(Channel2{1})},'mcp');
        if isempty(InputChannelTemp)
            InputChannelTemp=strfind({lower(Channel1{1}),lower(Channel2{1})},'pp7');
                if isempty(InputChannelTemp)
                    InputChannelTemp=strfind({lower(Channel1{1}),lower(Channel2{1})},'lambdan');
                end    
        end
        InputChannelTemp=cellfun(@isempty,InputChannelTemp);
        InputChannel = find(InputChannelTemp);
        
    elseif strcmpi(ExperimentType,'input')
        %Parse the information from the different channels
        Channels={Channel1{1},Channel2{1}};

        %Histone channel.
        histoneChannel=find(~cellfun(@isempty,strfind(lower(Channels),'his')));
        if isempty(histoneChannel)
            histoneChannel=0;
        end

        %Input channels
        InputChannel=~cellfun(@isempty,Channels);
        if histoneChannel
            InputChannel(histoneChannel)=0;
        end
        InputChannel=find(InputChannel);
    end
    
    
    %Create the circle that we'll use as the mask
    IntegrationRadius=2;       %Radius of the integration region in um
    IntegrationRadius=floor(IntegrationRadius/FrameInfo(1).PixelSize); %Radius of the integration in pixels
    if ~mod(IntegrationRadius,2)
        IntegrationRadius=IntegrationRadius+1;
    end
    Circle=false(3*IntegrationRadius,3*IntegrationRadius);
    Circle=MidpointCircle(Circle,IntegrationRadius,1.5*IntegrationRadius+0.5,...
        1.5*IntegrationRadius+0.5,1);
    
    %Initialize fields
    schnitzcells(1).Fluo = [];
    schnitzcells(1).Mask = [];
    
    if sum(InputChannel)
                
        %Extract the fluroescence of each schnitz, for each channel,
        %at each time point
         
        %Get the image dimensions
        PixelsPerLine=FrameInfo(1).PixelsPerLine;
        LinesPerFrame=FrameInfo(1).LinesPerFrame;
        %Number of z-slices
        NumberSlices=FrameInfo(1).NumberSlices;
        NumberSlices2 = NumberSlices + 2;
              
        %Generate reference frame for edge detection
        refFrame = ones(LinesPerFrame,PixelsPerLine,NumberSlices2);
        convRef = convn(refFrame, Circle, 'same');
        edgeMask = convRef~=sum(Circle(:));
        
        for ChN=1:length(InputChannel)
        
            h=waitbar(0,['Extracting nuclear fluorescence for channel ',num2str(ChN)]);
            
                if strcmpi(ExperimentType,'inputoutput') || (strcmpi(ExperimentType,'input')&&(length(InputChannel))==1)
                    nameSuffix=['_ch',iIndex(InputChannel,2)];
                elseif strcmpi(ExperimentType,'input')&&(length(InputChannel))>1
                    nameSuffix=['_ch',iIndex(InputChannel(ChN),2)];
                else
                    error('Not sure what happened here. Talk to YJK or SA.')
                end
            % NL: should parallelize this
          
            tempSchnitz = schnitzcells;
            for CurrentFrame=1:numFrames

                waitbar(CurrentFrame/numFrames,h);
% 
%                 %Initialize the image
                Image=zeros(LinesPerFrame,PixelsPerLine,NumberSlices2);
%                 %AR 1/12/18 moved this outside the for loops
%                 Image=zeros(LinesPerFrame,PixelsPerLine,NumberSlices2);
                %Load the z-stack for this frame        
                for CurrentZ=1:NumberSlices2   %Note that I need to add the two extra slices manually
                    Image(:,:,CurrentZ)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),nameSuffix,'.tif']);
                end
                
                % NL: rewriting this extraction step as a convolution                            
%                 convImage = convn(Image, Circle, 'same');    
                convImage = imfilter(Image, double(Circle), 'same');
                convImage(edgeMask) = NaN;
%                 [yDim, xDim] = size(convImage);
%                 tempFluo = []; tempMask = [];
                for j=1:length(tempSchnitz)
                    CurrentIndex=find(tempSchnitz(j).frames==CurrentFrame);
                    cenx=min(max(1,round(tempSchnitz(j).cenx(CurrentIndex))),PixelsPerLine);
                    ceny=min(max(1,round(tempSchnitz(j).ceny(CurrentIndex))),LinesPerFrame);
%                     tempFluo(j) = convImage(ceny,cenx,:);
%                      tempFluo(j).Fluo(CurrentIndex,1:NumberSlices2,ChN) = convImage(ceny,cenx,:);
                    tempSchnitz(j).Fluo(CurrentIndex,1:NumberSlices2,ChN) = convImage(ceny,cenx,:);
%                      tempMask(j).Mask = Circle;
                    tempSchnitz(j).Mask=Circle;                    
                end
                
%                 schnitzcells.Mask = tempMask;
                
%                 NL: commented out old version
%                 for j=1:length(schnitzcells)
% 
%                     %HG: Note that I'm calling a function here so that I can
%                     %debug the parfor loop above. Ideally, I would have 
%                     %parfor loop over images, not schnitzes within an image.
%                     %However, I couldn't quite figure out how to do that.
%                     %schnitzcells(j)=ExtractNuclearFluorescence(schnitzcells(j),...
%                     %    CurrentFrame,...
%                     %    Image,LinesPerFrame,PixelsPerLine,NumberSlices2,Circle,IntegrationRadius,InputChannel(ChN));
%                     schnitzcells(j)= ExtractNuclearFluorescence(schnitzcells(j),...
%                         CurrentFrame,...
%                         Image,LinesPerFrame,PixelsPerLine,NumberSlices2,Circle,IntegrationRadius,ChN);
%                 end
            end
            schnitzcells = tempSchnitz;
%             schnitzcells.Fluo = tempFluo;
%             schnitzcells.Mask = tempMask;
            
        close(h)
        end
    else
        error('Input channel not recognized. Check correct definition in MovieDatabase')
    end
end
%Save the information
%Now save
mkdir([DropboxFolder,filesep,Prefix])
save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses')

if strcmpi(ExperimentType,'inputoutput')||strcmpi(ExperimentType,'input')
    %Change the name of the Circle variable to make it more understandable when
    %loaded independently
    IntegrationArea=Circle;
    save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells','IntegrationArea')
else
    save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells')
end

if ~exist([ProcPath,filesep,Prefix,'_'], 'dir')
    mkdir([ProcPath,filesep,Prefix,'_'])
end
save([ProcPath,filesep,Prefix,'_',filesep,'dataStructure.mat'],'dataStructure')


%Stitch the schnitzcells using Simon's code
if ~SkipStitchSchnitz
    disp 'Skipping StitchSchnitz'
    StitchSchnitz(Prefix)
end

