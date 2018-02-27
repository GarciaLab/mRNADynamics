function [Particles, Spots, SpotFilter, schnitzcells] = CheckParticleTracking(varargin)
% function [Particles, Spots, SpotFilter, schnitzcells] = CheckParticleTracking(varargin)
%
% DESCRIPTION
% The point of this function is to check the tracking of particles. The
% logic of this function should be similar to schnitzcells: We want to be
% able to correct both the segmentation and tracking.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
%
% OPTIONS
% NoSort : Flag to sort or not particles according to their starting frame
% ForCompileAll : Flag to just save the data. This is good for CompileAll
% speedmode : Flag to plot only ellipses for current particle & save time
% sistermode : Decide whether you want to do sister chromatid analysis
% nc, NC : Only look at particles that show up in nc13 
%    % Currently this only starts at nc13...not restrict you to nc13 Added by Emma
%    % Also, this option shows you the max projection. 
%
% CONTROLS
% Frame specific:
% . , Move a frame forward/backward
% > < Move five frames forward/backward
% ; ' Move to the next empty frame within a particle
% a z Move up/down in Z
% j Jump to a specified frame
% g b Increase/decrease histone channel contrast
% !   Change the contranst in gfp channel
% 
% 
% Particle specific:
% m Move to the next particle
% n Move to the previous particle
% k Jump to a specified particle by inputting particle number
% \ Jump to a specified particle by clicking
% c Connect two existing particle traces. This will join the current
%  particle's trace to the clicked particle's trace. 
% d Separate traces forward. A new particle is created at the current frame
%  and this particle is disconnected from the current nucleus. If this is
%  done on a particle with only one frame then
%  it disconnects it from its nucleus.
%
% Disconnect backwards??
%
% q Cycle between approved status: green - approved; yellow - approved but
%  with conditions (drift of nucleus, for example)
% w Disapprove a trace
% p Identify a particle. It will also tell you the particle associated with
%  the clicked nucleus.
% e Approve/Disapprove a frame within a trace
% u Move a particle detected with Threshold2 into the our structure. 
% i Move a particle detected with Threshold2 into the our structure and
%  connect it to the current particle. This is a combination of "u" and
%  "c". %AR 1/15/18: This is currently deprecated. Talk to HG if you need
%  this function.
% [ Add a spot that was not recognized originally by segmentSpots to the
% 	current particle. Note that the command forces ZoomMode. To toggle, use
%   'o' or '+' depending on whether you're adding to an existing trace or creating a new
%    trace, respectively.
% 
% 
% 
% Nuclear tracking specific:
% l Split a nucleus and select one or two daughter nuclei or stop the
%  lineage. Usage:
%       Click on one new nucleus + ENTER: Continue the schnitz with that nucleus.
%       Click on the current nucleus + ENTER: Split the schnitz. This time
%           point will be the first frame of the new schnitz.
%       Click on two nuclei: Split the current nucleus into two daughter
%       nuclei.
%       Click on the same nucleus twice: Split the current nucleus, but
%       with only one daughter nucleus.
% 2 set parent of current nucleus
% p Find the particle associated with the clicked nucleus. It will also tell
%  you the closest particle associated you clicked on.
% 
% 
% General:
% 8 Change channels
% t Show/hide particles from the second threshold
% s Save the current Particles structure
% x Save and exit
% h Show non-approved particles yellow or dissapproved particlesz
% y Input the frame/nc information again. This only works in the absence of
%  the histone channel
% r Reorder the particles according to initial frame
% f Redo tracking. It only gets done on the non-approved particles.
% o Zoom in/out around the particle's first frame.
% + Zoom anywhere button. Click with the mouse after hitting this. 
% -/= Change the zoom factor when in zoom mode.
% 0 Enter debug mode to fix things manually
% ~ Switch figure 1 from a single plane image to a z or time projection. 
%
% OUTPUT
% Particles: A modified Particles
% Spots: A modified Spots
% SpotFilter: A modified SpotFilter
% schnitzcells: A modified schnitzcells
%
% Author (contact): Hernan Garcia (hgarcia@berkeley.edu)
% Created: 
% Last Updated: 1/13/2018

close all


warning('off','MATLAB:nargchk:deprecated')


%% Information about about folders

%Get the folders
[~,~,DefaultDropboxFolder,~,PreProcPath]=...
    DetermineLocalFolders;

%Also get the computer name. We'll use this later

%Find out which computer this is. That will determine the folder structure.
[ret, name] = system('hostname');  
if ret ~= 0 
   if ispc  
      name = lower(getenv('COMPUTERNAME'));  
   else  
      name = lower(getenv('HOSTNAME'));  
   end  
end  

if isempty(varargin)
    DataFolder=uigetdir(DefaultDropboxFolder,'Select data set to analyze');
else
    [~,~,DropboxFolder,~,~]=...
        DetermineLocalFolders(varargin{1});
    DataFolder=[DropboxFolder,filesep,varargin{1}];
end
    
%Flag to sort or not particles according to their starting frame
NoSort=0;
%Flag to just save the data. This is good for CompileAll
ForCompileAll=0;
%Flag to plot only ellipses for current particle & save time
SpeedMode = 0;
%Decide whether you want to do sister chromatid analysis
SisterMode = 0;
%Decide whether you want to only see nc13
ncRange = 0;
% This is for the projection mode
projectionMode = 'None (Default)';

Prefix = varargin{1};
if length(varargin)>1
    for i=2:length(varargin)
        if strcmpi(varargin{i},'NoSort')
            NoSort=1;
        elseif strcmpi(varargin{i},'ForCompileAll')
            ForCompileAll=1;
        elseif strcmpi(varargin{i}, 'speedmode')
            SpeedMode = 1;
        elseif strcmpi(varargin{i}, 'sistermode')
            SisterMode = 1;
        elseif strcmpi(varargin{i},'nc') % checking for the desired nc range
            ncRange = 1;
            NC = varargin{i+1};
            % startNC and endNC will be the varibale names that have the start of the nc(s) of interest
            if length((varargin{i+1})) == 2
                startNC = ['nc' num2str(varargin{i+1}(1))];
                endNC = ['nc' num2str(varargin{i+1}(2) +1)];% Not including the next nc
            else
                startNC = ['nc' num2str(varargin{i+1})]; 
                endNC = ['nc' num2str(varargin{i+1} + 1)]; % Not including the next nc
            end
        end
    end
end

%%

FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];

%Now get the actual folders
[~,~,DropboxFolder,~,PreProcPath]=...
    DetermineLocalFolders(FilePrefix(1:end-1));

load([DataFolder,filesep,'Particles.mat'])
load([DataFolder,filesep,'Spots.mat'])

%Check that FrameInfo exists
if exist([DataFolder,filesep,'FrameInfo.mat'], 'file')
    load([DataFolder,filesep,'FrameInfo.mat'])
else
    warning('No FrameInfo.mat found. Trying to continue')
    %Adding frame information
    DHis=dir([PreProcPath,filesep,FilePrefix(1:end-1),filesep,'*His*.tif']);
    FrameInfo(length(DHis)).nc=[];
    %Adding information
   Dz=dir([PreProcPath,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'*001*.tif']);
    NumberSlices=length(Dz)-1;
    for i=1:numFrames
        FrameInfo(i).NumberSlices=NumberSlices;
    end
end
xSize = FrameInfo(1).PixelsPerLine;
ySize = FrameInfo(1).LinesPerFrame;
numFrames =length(FrameInfo);

%See how  many frames we have and adjust the index size of the files to
%load accordingly
if numFrames<1E3
    NDigits=3;
elseif numFrames<1E4
    NDigits=4;
else
    error('No more than 10,000 frames supported. Change this in the code')
end



%Some parameters:
SnippetEdge=13;     %Size of the snippets generated by Michael's code in pixels.
storedTimeProjection = []; % Don't need to wait for timeProjection to finish each time its called

%Create the particle array. This is done so that we can support multiple
%channels. Also figure out the number of channels
if iscell(Particles)
    NChannels=length(Particles);
else
    Particles={Particles};
    Spots={Spots};
    SpotFilter={SpotFilter};
    NChannels=1;
end


%Add FramesApproved where necessary
for NCh=1:NChannels
    if ~isfield(Particles{NCh},'FrameApproved')
        for i=1:length(Particles{NCh})
            Particles{NCh}(i).FrameApproved=true(size(Particles{NCh}(i).Frame));
        end
    else
        for i=1:length(Particles{NCh})
            if isempty(Particles{NCh}(i).FrameApproved)
                Particles{NCh}(i).FrameApproved=true(size(Particles{NCh}(i).Frame));
            end
        end
    end
end



%Check if we have the histone channel and we have done the nuclear
%segmentation.
if exist([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
        FilePrefix(1:end-1),'-His_',iIndex(1,NDigits),'.tif'], 'file')||...
        exist([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
        FilePrefix(1:end-1),'_His_',iIndex(1,NDigits),'.tif'], 'file')
    %(MT, 2018-02-11) Added support for lattice imaging with bad histone 
    %channel, maybe temporary - FIX LATER
    if exist([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,'Ellipses.mat'])
        load([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,'Ellipses.mat'])
        UseHistoneOverlay=1;
    else
        warning('Ellipses.mat does not exist. Proceeding as though there is no Histone channel. If you expect a Histone channel, there is something wrong.')
        UseHistoneOverlay=0;
    end
else
    UseHistoneOverlay=0;
end

%Check that we have the nuclear tracking done using schnitzcells
if exist([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'], 'file')
    UseSchnitz=1;
else
    UseSchnitz=0;
end

% we name the variable DataFolderColumnValue to avoid shadowing previously defined DataFolder var, which is actually a subfolder inside dropbox
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolderColumnValue, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

if exist([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'file')
    
    for i=1:numFrames
        if i<nc9
            FrameInfo(i).nc=8;
        elseif (i>=nc9)&&(i<nc10)
            FrameInfo(i).nc=9;
        elseif (i>=nc10)&&(i<nc11)
            FrameInfo(i).nc=10;
        elseif (i>=nc11)&&(i<=nc12)
            FrameInfo(i).nc=11;
        elseif (i>=nc12)&&(i<=nc13)
            FrameInfo(i).nc=12;
        elseif (i>=nc13)&&(i<=nc14)
            FrameInfo(i).nc=13;
        elseif i>=nc14
            FrameInfo(i).nc=14;
        end
    end
else
    warning('No nuclear marker channel may result in strange behavior.');
    
    for i=1:numFrames
        if i<nc9
            FrameInfo(i).nc=8;
        elseif (i>=nc9)&&(i<nc10)
            FrameInfo(i).nc=9;
        elseif (i>=nc10)&&(i<nc11)
            FrameInfo(i).nc=10;
        elseif (i>=nc11)&&(i<=nc12)
            FrameInfo(i).nc=11;
        elseif (i>=nc12)&&(i<=nc13)
            FrameInfo(i).nc=12;
        elseif (i>=nc13)&&(i<=nc14)
            FrameInfo(i).nc=13;
        elseif i>=nc14
            FrameInfo(i).nc=14;
        end
    end
end

save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo') %this is here so that a user will still get an updated
%frameinfo.mat even if they abort checkparticletracking without saving (to
%prevent issues with compileparticles)

%Check if we have already determined nc
if (~isfield(FrameInfo,'nc'))&&(~UseHistoneOverlay)
    %FrameInfo=DetermineNC(fad,Particles,FrameInfo);  AR 3/14/16: This
    %script seems to have disappeared. 

elseif UseSchnitz
    load([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'])
    
    %Remove the schnitz fields that can give us problems potentially if
    %present. I don't know how this came to be, but it's for fields that
    %are not all that relevant. The fields are: approved, ang
    if isfield(schnitzcells,'approved')
        schnitzcells=rmfield(schnitzcells,'approved');
    end
    if isfield(schnitzcells,'ang')
        schnitzcells=rmfield(schnitzcells,'ang');
    end
end

% %Load the DoG images. Necessary for particle addition/subtraction
% dog = [];
% num_frames = numFrames;
% zSize = FrameInfo(1).NumberSlices + 2;
% OutputFolder1=[FISHPath,filesep,Prefix,'_',filesep,'dogs'];
% for current_frame = 1:num_frames
%     for i = 1:zSize   
%         dog(:,:,i,current_frame) = double(imread([OutputFolder1, filesep,'DOG_',Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),'.tif']));
%     end
% end

%Order particles by the earliest frame they appear at. This makes the
%tracking a lot easier!
if ~NoSort
    for ChN=1:NChannels
        for i=1:length(Particles{ChN})
            FirstFrame(i)=Particles{ChN}(i).Frame(1);
        end
        [~,Permutations]=sort(FirstFrame);
        Particles{ChN}=Particles{ChN}(Permutations);
        clear FirstFrame
    end
end

%Some flags and initial parameters
ShowThreshold2=0;                    %Whether to show particles below the threshold
HideApprovedFlag=0;
ParticleToFollow=[];
ZSlices=FrameInfo(1).NumberSlices+2; %Note that the blank slices are included
CurrentZ=round(ZSlices/2);          
ManualZFlag=0;
CurrentParticle=1;
PreviousParticle=1;
CurrentFrameWithinParticle=1;
CurrentChannel=1;
PreviousChannel=CurrentChannel;
CurrentFrame=Particles{1}(1).Frame(1);
DisplayRange=[];
DisplayRangeSpot=[];
ZoomMode=0;
GlobalZoomMode=0;
ZoomRange=50;
nameSuffix='';

%Set up the default contrast settings for the MCP channel depending on the
%microscope that was used used
if strcmpi(FrameInfo(1).FileMode,'dspin')
    %For spinning disk, we set the contrast to the maximum and minimum
    minContrast=[];
    maxContrast=[];
else
    %For all other microscopes, we have a default. HG is not sure this will
    %actually work well beyond Leica SP8.
    minContrast = 0; % Default contrast settings for gfp channel
    maxContrast = 80;
end

% Changing the intial frames and particle if justNC13
if ncRange
    
    if strcmpi('nc15',endNC) 
        lastNCFrame = numFrames;
    else
        lastNCFrame = eval(endNC)-1; % This will not include the 1st frame of the next NC
    end
    firstNCFrame = eval(startNC);
    particlesInRange = particlesWithinFrames(Prefix,firstNCFrame,lastNCFrame);
    CurrentParticle = particlesInRange(1);
    CurrentFrame = Particles{1}(CurrentParticle).Frame(1);
%     ncRangeFigure = figure();
%     set(gcf,'units', 'normalized', 'position',[0.35, 0.55, .2, .33])
%     uicontrol('Parent', ncRangeFigure, 'Style', 'text','String','Implementing justNC13','Units','normalized', 'Position', [0.25 0.5 0.5 0.35])
%     uicontrol('Parent', ncRangeFigure, 'Style', 'text','String',['nc13 : ' num2str(firstNCFrame)],'Units','normalized', 'Position', [0.25 0.40 0.5 0.35])
%     uicontrol('Parent', ncRangeFigure, 'Style', 'text','String',['Number of Particles: ' num2str(length(particlesInRange))],'Units','normalized','Position', [0.25 0.30 0.5 0.35])
%     uicontrol('Parent', ncRangeFigure, 'Style', 'text','String',['Particles in range: ' num2str(particlesInRange)],'Units','normalized','Position', [0.25 0.20 0.5 0.35])
    disp(['nc range: ' num2str(NC)])
    disp(['start frame: ' num2str(firstNCFrame)])
    disp(['end frame: ' num2str(lastNCFrame)])
    disp(['Particles in range: ' num2str(particlesInRange)])
    disp(['Number of Particles: ' num2str(length(particlesInRange))])
end

%Define the windows
Overlay=figure;
if UseHistoneOverlay 
    HisOverlayFig=figure;   
end
TraceFig=figure;
SnippetFig=figure;
ZProfileFig=figure;
% SisterFig=figure;
% SisterFig2 = figure;
% SisterFig3 = figure;
Gaussian = figure;
RawData = figure;


cc=1;

% Create the approved field if it does not exist
for NCh=1:NChannels
    if ~isfield(Particles{NCh},'Approved')
        for i=1:length(Particles{NCh})
            Particles{NCh}(i).Approved=0;
        end
    end
end    
    
%See if we just want to save the data
if ForCompileAll
    cc='x';
end


%Figure out channel-specific information
if NChannels==1
    if contains(Channel1{1}, 'MCP') || contains(Channel1{1}, 'PCP')
        nameSuffix=['_ch',iIndex(1,2)];
        coatChannel = 1;
    elseif contains(Channel2{1}, 'MCP') || contains(Channel2{1}, 'PCP')
        nameSuffix=['_ch',iIndex(2,2)];
        coatChannel = 2;
    end
elseif strcmpi(ExperimentType,'2spot2color')
        %We are assuming that channels 1 and 2 are assigned to coat
        %proteins. We should do a better job with this.
        coatChannels = [1,2];
        coatChannel=coatChannels(1);
else
    error('Experiment type not recognized')
end



%This flag allows the code to directly pass a command without waiting for
%the user to press a key or click on the figure
SkipWaitForButtonPress=[];


while (cc~='x')
    
    %Update the name suffic
    if strcmpi(ExperimentType,'2spot2color')
        nameSuffix=['_ch',iIndex(coatChannel,2)];
    end
    
    EllipseHandle=[];
    EllipseHandleYellow=[];
    EllipseHandleBlue=[];
    EllipseHandleWhite=[];
    EllipseHandleGreen=[];

    numParticles = length(Particles{CurrentChannel});
    
    %Get the coordinates of all the spots in this frame
    [x,y,z]=SpotsXYZ(Spots{CurrentChannel}(CurrentFrame));
        
    %If the approved field does not exist create it
    if ~isfield(Particles{CurrentChannel},'Approved')
        for i=1:numParticles
            Particles{CurrentChannel}(i).Approved=0;
        end
    end
    ApprovedParticles=[Particles{CurrentChannel}.Approved];
    
    %Pull out the right particle if it exists in this frame
    CurrentParticleIndex=...
        Particles{CurrentChannel}(CurrentParticle).Index(Particles{CurrentChannel}(CurrentParticle).Frame==...
        CurrentFrame);
    %This is the position of the current particle
    xTrace=x(CurrentParticleIndex);
    yTrace=y(CurrentParticleIndex);
    
    %These are the positions of all the approved, disapproved and
    %unflagged particles
    
    %Approved particles
    IndexApprovedParticles=[];
    for i=1:numParticles
        if sum(Particles{CurrentChannel}(i).Frame==CurrentFrame)&&...
                sum(Particles{CurrentChannel}(i).Approved==1)
            IndexApprovedParticles=[IndexApprovedParticles,...
                Particles{CurrentChannel}(i).Index(Particles{CurrentChannel}(i).Frame==CurrentFrame)];
        end
    end
    xApproved=x(IndexApprovedParticles);
    yApproved=y(IndexApprovedParticles);

    %Disapproved particles
    IndexDisapprovedParticles=[];
    for i=1:numParticles
        if sum(Particles{CurrentChannel}(i).Frame==CurrentFrame)&&sum(Particles{CurrentChannel}(i).Approved==-1)
            IndexDisapprovedParticles=[IndexDisapprovedParticles,...
                Particles{CurrentChannel}(i).Index(Particles{CurrentChannel}(i).Frame==CurrentFrame)];
        end
    end
    xDisapproved=x(IndexDisapprovedParticles);
    yDisapproved=y(IndexDisapprovedParticles);
    
    %Non-flagged particles (these are particles that have not been
    %processed)
    IndexNonFlaggedParticles=[];
    for i=1:numParticles
        if sum(Particles{CurrentChannel}(i).Frame==CurrentFrame)&&...
                ~(sum(Particles{CurrentChannel}(i).Approved==-1)||sum(Particles{CurrentChannel}(i).Approved==1))
            IndexNonFlaggedParticles=[IndexNonFlaggedParticles,...
                Particles{CurrentChannel}(i).Index(Particles{CurrentChannel}(i).Frame==CurrentFrame)];
        end
    end
    xNonFlagged=x(IndexNonFlaggedParticles);
    yNonFlagged=y(IndexNonFlaggedParticles);
    
    if (~isempty(xTrace))&&(~ManualZFlag)
        CurrentZ=z(CurrentParticleIndex);
        CurrentZIndex=find(...
                Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).z==...
                CurrentZ);
        ManualZFlag=0;
    end
        
    if NChannels==1 % inputoutput mode can also be in this case, changed CurrentChannel to the coatChannel (YJK : 1/15/2018)
        if strcmpi(projectionMode,'None (Default)')
            Image=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                FilePrefix,iIndex(CurrentFrame,NDigits),'_z',iIndex(CurrentZ,2),nameSuffix,'.tif']);
        elseif strcmpi(projectionMode,'Max Z')
            [Image,~] = zProjections(Prefix, coatChannel, CurrentFrame, ZSlices, NDigits,DropboxFolder,PreProcPath);
        elseif strcmpi(projectionMode,'Median Z')
            [~,Image] = zProjections(Prefix, coatChannel, CurrentFrame, ZSlices, NDigits,DropboxFolder,PreProcPath);
        elseif strcmpi(projectionMode,'Max Z and Time')
            if isempty(storedTimeProjection)
                if ncRange
                    Image = timeProjection(Prefix, coatChannel,'nc',NC);
                    storedTimeProjection = Image;
                else
                    Image = timeProjection(Prefix, CurrentChannel);
                    storedTimeProjection = Image;
                end
            else
                Image = storedTimeProjection;
            end
        end

    elseif NChannels>1
        Image=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
            FilePrefix,iIndex(CurrentFrame,NDigits),'_z',iIndex(CurrentZ,2),...
            nameSuffix,'.tif']);
    else
        error('ExperimentType and/or channel not supported.')          
    end
    
    figure(Overlay)
%     imshow(Image,[minContrast maxContrast],'Border','Tight') %AR
%     1/13/2018 this contrast setting does not work for dim particles.
    imshow(Image,DisplayRangeSpot,'Border','Tight')
    hold on
    %Show all particles in regular mode
    if ~SpeedMode
        plot(xNonFlagged,yNonFlagged,'ow')
        plot(xApproved,yApproved,'ob')
        plot(xDisapproved,yDisapproved,'^r')
        %plot(x, y, 'sw')
    end
    %Always show current particle
    plot(xTrace,yTrace,'og')
    hold off

    if isfield(FrameInfo, 'nc')
    set(gcf,'Name',['Particle: ',num2str(CurrentParticle),'/',num2str(numParticles),...
        ', Frame: ',num2str(CurrentFrame),'/',num2str(numFrames),...
        ', Z: ',num2str(CurrentZ),'/',num2str(ZSlices),' nc: ', num2str(FrameInfo(CurrentFrame).nc),...
        ', Ch: ',num2str(CurrentChannel)])
    end
    if UseSchnitz
        %Show all the nuclei in regular mode
        if ~SpeedMode
            hold on
            EllipseHandle=notEllipse(Ellipses{CurrentFrame}(:,3),...
                Ellipses{CurrentFrame}(:,4),...
                Ellipses{CurrentFrame}(:,5),...
                Ellipses{CurrentFrame}(:,1)+1,...
                Ellipses{CurrentFrame}(:,2)+1,'r',50);
            hold off
            
            
            %Show the ones that have been approved
            
            hold on
            schnitzCellNo=[];
            for i=1:numParticles
                if Particles{CurrentChannel}(i).Approved==1
                    schnitzIndex=find((schnitzcells(Particles{CurrentChannel}(i).Nucleus).frames)==CurrentFrame);
                    schnitzCellNo=[schnitzCellNo,schnitzcells(Particles{CurrentChannel}(i).Nucleus).cellno(schnitzIndex)];
                end
            end
            
            EllipseHandleBlue=notEllipse(Ellipses{CurrentFrame}(schnitzCellNo,3),...
                Ellipses{CurrentFrame}(schnitzCellNo,4),...
                Ellipses{CurrentFrame}(schnitzCellNo,5),...
                Ellipses{CurrentFrame}(schnitzCellNo,1)+1,...
                Ellipses{CurrentFrame}(schnitzCellNo,2)+1,'b',50);                       
            hold off
        end
        
        %Show the corresponding nucleus
        if ~isempty(Particles{CurrentChannel}(CurrentParticle).Nucleus)
            SchnitzIndex=find(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames==CurrentFrame);
            NucleusIndex=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).cellno(SchnitzIndex);

            if ~isempty(NucleusIndex)
                hold on
                EllipseHandleGreen=ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                    Ellipses{CurrentFrame}(NucleusIndex,4),...
                    Ellipses{CurrentFrame}(NucleusIndex,5),...
                    Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                    Ellipses{CurrentFrame}(NucleusIndex,2)+1);
                set(EllipseHandleGreen,'Color','g')
                hold off
            else
                %('Error: Particle without an associated nucleus?')
            end
        
        
        
            %Show the daughter nuclei if applicable
            DaughterE=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).E;
            DaughterD=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).D;

            
            if DaughterE~=0
                SchnitzIndex=find(schnitzcells(DaughterE).frames==CurrentFrame);
                NucleusIndex=schnitzcells(DaughterE).cellno(SchnitzIndex);

                if ~isempty(NucleusIndex)
                    hold on
                    EllipseHandleWhite=[EllipseHandleWhite,ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                        Ellipses{CurrentFrame}(NucleusIndex,4),...
                        Ellipses{CurrentFrame}(NucleusIndex,5),...
                        Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                        Ellipses{CurrentFrame}(NucleusIndex,2)+1)];
                    
                    hold off
                else
                    %('Error: Particle without an associated nucleus?')
                end
            end

            if DaughterD~=0
                SchnitzIndex=find(schnitzcells(DaughterD).frames==CurrentFrame);
                NucleusIndex=schnitzcells(DaughterD).cellno(SchnitzIndex);

                if ~isempty(NucleusIndex)
                    hold on
                    EllipseHandleWhite=[EllipseHandleWhite,ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                        Ellipses{CurrentFrame}(NucleusIndex,4),...
                        Ellipses{CurrentFrame}(NucleusIndex,5),...
                        Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                        Ellipses{CurrentFrame}(NucleusIndex,2)+1)];
                    hold off
                else
                    %('Error: Particle without an associated nucleus?')
                end
            end
            
            if ~isempty(EllipseHandleWhite)
                set(EllipseHandleWhite,'Color','w')
            end
            
            %Show the mother nucleus if applicable
            Mother=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).P;

            if Mother~=0
                SchnitzIndex=find(schnitzcells(Mother).frames==CurrentFrame);
                NucleusIndex=schnitzcells(Mother).cellno(SchnitzIndex);

                if ~isempty(NucleusIndex)
                    hold on
                    EllipseHandleYellow=ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                        Ellipses{CurrentFrame}(NucleusIndex,4),...
                        Ellipses{CurrentFrame}(NucleusIndex,5),...
                        Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                        Ellipses{CurrentFrame}(NucleusIndex,2)+1);
                    set(EllipseHandleYellow,'Color','y')
                    hold off
                else
                    %('Error: Particle without an associated nucleus?')
                end
            end
            
        else
            if UseHistoneOverlay
                warning('This particle does not have an associated nucleus.');
            end
        end  
    end
    

    
    if ApprovedParticles(CurrentParticle)==1
        set(gcf,'Color','g')
    elseif ApprovedParticles(CurrentParticle)==-1
        set(gcf,'Color','r')
    elseif ApprovedParticles(CurrentParticle)==2
        set(gcf,'Color','y')
    else
        set(gcf,'Color','default')
    end
    
    %Show the particles that were under threshold 2.
    if ShowThreshold2
        %Get the positions of all the spots in this frame
        [x2,y2]=SpotsXYZ(Spots{CurrentChannel}(CurrentFrame));
        %Filter those that were under threshold 2.
        CurrentSpotFilter=...
            ~logical(SpotFilter{CurrentChannel}(CurrentFrame,~isnan(SpotFilter{CurrentChannel}(CurrentFrame,:))));
        x2=x2(CurrentSpotFilter);
        y2=y2(CurrentSpotFilter);
        
        hold on
        plot(x2,y2,'sr')
        hold off
    end
    
    if ZoomMode
        %Find the closest frame
        [~,MinIndex]=min((Particles{CurrentChannel}(CurrentParticle).Frame-CurrentFrame).^2);
        if length(MinIndex)>1
            MinIndex=MinIndex(1);
        end
        [xForZoom,yForZoom]=...
            SpotsXYZ(Spots{CurrentChannel}(Particles{CurrentChannel}(CurrentParticle).Frame(MinIndex)));
        
        xForZoom=xForZoom(Particles{CurrentChannel}(CurrentParticle).Index(MinIndex));
        yForZoom=yForZoom(Particles{CurrentChannel}(CurrentParticle).Index(MinIndex));
       
 
        xlim([xForZoom-ZoomRange,xForZoom+ZoomRange])
        ylim([yForZoom-ZoomRange/2,yForZoom+ZoomRange/2])
    end
    
    if GlobalZoomMode       
        xlim([xForZoom-ZoomRange,xForZoom+ZoomRange])
        ylim([yForZoom-ZoomRange/2,yForZoom+ZoomRange/2])
    end
        
   
    if UseHistoneOverlay
        figure(HisOverlayFig)
        try
            ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,NDigits),'.tif']);
        catch %Had to do this for KITP
            ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                FilePrefix(1:end-1),'_His_',iIndex(CurrentFrame,NDigits),'.tif']);
        end

        if isempty(DisplayRange)
            HisOverlayImage=cat(3,mat2gray(ImageHis),mat2gray(Image),zeros(size(Image)));
        else
            HisOverlayImage=cat(3,mat2gray(ImageHis,double(DisplayRange)),mat2gray(Image),zeros(size(Image)));
        end
        imshow(HisOverlayImage,[],'Border','Tight')

       
        hold on
        if ~SpeedMode
            plot(xNonFlagged,yNonFlagged,'ow')
            plot(xApproved,yApproved,'ob')
        end
        plot(xTrace,yTrace,'og')
        hold off
        
        if ShowThreshold2
            hold on
            plot(x2,y2,'sw')
            hold off
        end
  
        
        if UseSchnitz
            
            copyobj(EllipseHandle,gca)
            copyobj(EllipseHandleBlue,gca)
            copyobj(EllipseHandleGreen,gca)
            copyobj(EllipseHandleWhite,gca)
            copyobj(EllipseHandleYellow,gca)
            
        end
    
        set(gcf,'Name',['Particle: ',num2str(CurrentParticle),'/',num2str(numParticles),...
            ', Frame: ',num2str(CurrentFrame),'/',num2str(numFrames),...
            ', Z: ',num2str(CurrentZ),'/',num2str(ZSlices),' nc: ', num2str(FrameInfo(CurrentFrame).nc),...
            ' Ch: ',num2str(CurrentChannel)])
        
        if ZoomMode || GlobalZoomMode
            xlim([xForZoom-ZoomRange,xForZoom+ZoomRange])
            ylim([yForZoom-ZoomRange/2,yForZoom+ZoomRange/2])
        end
 
    end
    
%     %AR 7/14/16: Need to fill in the details.
%     figure(SisterFig)
%     %plot sister 1 versus time
%     hold on
%     %plot sister 2 versus time
%     figure(SisterFig2)
%     %plot distance versus intensity
%     figure(SisterFig3)
%     %plot distance versus time
    if ~isempty(xTrace)
        MaxZIndex=find(...
            Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).z==...
            Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).brightestZ);
        CurrentZIndex=find(...
            Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).z==...
            CurrentZ);
        if isempty(CurrentZIndex)
%             warning('This particle has a gap in its z-profile. This is
%             highly suspect.'); %this if statement should only happen
%             between two spots, not past the PSF boundaries
        end
    end
    
    %Check to see if spots structure contains multi-slice fields
    multi_slice_flag = isfield(Spots{CurrentChannel}(CurrentFrame).Fits...
        (CurrentParticleIndex),'IntegralZ');
    
    figure(SnippetFig)
    if  ~isempty(xTrace) && ~isempty(CurrentZIndex)  
        %Get the snippet and the mask, and overlay them  
        %(MT, 2018-02-12): lattice data could use this, changed CurrentChannel to coatChannel 
        FullSlice=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(CurrentFrame,3)...
            ,'_z' iIndex(CurrentZ,2) '_ch' iIndex(coatChannel,2) '.tif']);
        xSpot = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).xDoG(CurrentZIndex);
        ySpot = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).yDoG(CurrentZIndex);
        if isfield(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex), 'snippet_size')
            snippet_size = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).snippet_size;
        %(MT, 2018-02-12): Hacky fix to get this to run with lattice data -
        %FIX LATER
        elseif strcmpi(ExperimentType,'lattice')
            snippet_size = 13;
        else
            snippet_size = 7;
        end
        if isempty(snippet_size)
                snippet_size = 7;
        end
%         CurrentSnippet=mat2gray(...
%             Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).Snippet{CurrentZIndex});
        
%         snippet_size = 15;        
        CurrentSnippet = FullSlice(max(1,ySpot-snippet_size):min(ySize,ySpot+snippet_size),...
                                max(1,xSpot-snippet_size):min(xSize,xSpot+snippet_size));
        imSnippet = mat2gray(CurrentSnippet);
        SnippetEdge=size(CurrentSnippet,1);     
        IntegrationRadius = 6; % this appears to be hard-coded into IdentifySingleSpot
        [xGrid, yGrid] = meshgrid(1:SnippetEdge,1:SnippetEdge);
        rGrid = sqrt((xGrid-ceil(SnippetEdge/2)).^2 + (yGrid-ceil(SnippetEdge/2)).^2); 
        SnippetMask = rGrid <= IntegrationRadius;
        IntegrationArea=bwperim(SnippetMask);
        
        SnippetOverlay=cat(3,IntegrationArea/2 + ...
            +imSnippet,imSnippet,imSnippet);
    
        imshow(SnippetOverlay,...
            [],'Border','Tight','InitialMagnification',1000)

        hold on
        
        SnippetX=(SnippetEdge-1)/2+1-...
            (Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).xDoG(CurrentZIndex)-...
            Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).xFit(CurrentZIndex));
        SnippetY=(SnippetEdge-1)/2+1-...
            (Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).yDoG(CurrentZIndex)-...
            Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).yFit(CurrentZIndex));
        hold off
    else
        imshow(zeros(SnippetEdge))
    end
    DoubleSnippet = double(CurrentSnippet);
    [mesh_y,mesh_x] = meshgrid(1:size(DoubleSnippet,2), 1:size(DoubleSnippet,1));

    % Single gaussian function: In future this should be a standalone
    % function file to ensure consistency with function used for fitting
    singleGaussian = @(params) (params(1).*...
        exp(-(...
        (((cos(params(7)))^2 / (2*params(3)^2) ) + ((sin(params(7)))^2 / 2*params(5)^2))  .* (mesh_x-params(2)).^2 ...
        - 2*((-sin(2*params(7)) / (4*params(3)^2) ) + (sin(2*params(7)) / 4*params(5)^2)) .* (mesh_x-params(2)).*(mesh_y-params(4))...
        + (((sin(params(7)))^2 / (2*params(3)^2) ) + ((cos(params(7)))^2 / 2*params(5)^2)).* (mesh_y-params(4)).^2 ...
            )))...
        + params(6) - DoubleSnippet;

    if ~isempty(xTrace) && ~isempty(CurrentZIndex)
        figure(Gaussian)
            %Get the snippet and the mask, and overlay them
%             domain=Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).rawSpot{CurrentZIndex};
%             snip = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).Snippet{CurrentZIndex};
            if isfield(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex),'gaussParams')
                gaussParams = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).gaussParams;
                if ~isempty(gaussParams)
                    gaussParams= gaussParams{CurrentZIndex};
                    gauss = singleGaussian(gaussParams);
                else
                    gauss = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).gaussSpot{CurrentZIndex};
                end
            elseif isfield(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex), 'gaussSpot')
                gauss = Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).gaussSpot{CurrentZIndex};
            else
                error('No Gaussian Fit Params or Gauss Snippet Found. Try Re-running segmentSpots')
            end

            surf(gauss + double(CurrentSnippet));
            title('Gaussian fit')
            set(gcf,'units', 'normalized', 'position',[0.815, 0.15, .2/3*2, .33/3*2]);
            zlimit = max(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).CentralIntensity);
            zlim([0, zlimit]);        
        figure(RawData)
%             surf(domain{1}, domain{2}, CurrentSnippet)
            surf(CurrentSnippet)
            title('Raw data');
            set(gcf,'units', 'normalized', 'position',[.6, 0.15, .2/3*2, .33/3*2]);
            zlim([0, zlimit]);
    else
        figure(Gaussian)
        clf;
        figure(RawData)
        clf;
    end

        
    figure(ZProfileFig)
    if ~isempty(xTrace)

        %Get the z DoG profile
        % check to see if Spots contains flag indicating type of
        % integration used
        title_string = '(max intensity)';
        if isfield(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex),'IntegralZ')
            if Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).IntegralZ
                ZProfile=Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).FixedAreaIntensity;
                title_string = '(raw integral)';
                IntegralZ_flag = 1;
            else
                ZProfile=Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).CentralIntensity;
                IntegralZ_flag = 0;
            end
        else
            ZProfile=Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).CentralIntensity;
            IntegralZ_flag = 0;
        end
        MaxZ=Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).brightestZ;
        
        plot(Spots{CurrentChannel}(CurrentFrame).Fits(CurrentParticleIndex).z,...
            ZProfile,'.-k');
        hold on
        if ~isempty(CurrentZIndex)
            plot(CurrentZ,ZProfile(CurrentZIndex),'ob')
        else
            plot(CurrentZ,CurrentZ,'or')
        end
        hold off
        title(['Z profile ' title_string],'FontSize',6)
    end

    figure(TraceFig)
    if ~strcmpi(ExperimentType,'inputoutput')
        %Only update the trace information if we have switched particles
        if (CurrentParticle~=PreviousParticle)||~exist('AmpIntegral', 'var')||(CurrentChannel~=PreviousChannel) 
            PreviousParticle=CurrentParticle;
            [Frames,AmpIntegral,GaussIntegral,AmpIntegral3,AmpIntegral5]=PlotParticleTrace(CurrentParticle,Particles{CurrentChannel},Spots{CurrentChannel});
        end       
       
        p1 = plot(Frames(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
            AmpIntegral(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-k');           
        hold on
        p2 = plot(Frames(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
            AmpIntegral3(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-','Color','green');                   
        p3 = plot(Frames(Particles{CurrentChannel}(CurrentParticle).FrameApproved),...
            GaussIntegral(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-','Color','blue');                   
        plot(Frames(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),AmpIntegral(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r');
        plot(Frames(Frames==CurrentFrame),AmpIntegral(Frames==CurrentFrame),'ob');
        legend([p1,p2,p3],'1 Slice Raw','3 Slice Raw','1 Slice Gaussian Fit')
        hold off
        try
            xlim([min(Frames)-1,max(Frames)+1]);
        catch
%             error('Not sure what happened here. Problem with trace fig x lim. Talk to AR if you see this, please.');
        end
        xlabel('frame')        
        ylabel('integrated intensity (a.u.)')
    else
        %Only update the trace information if we have switched particles
        if (CurrentParticle~=PreviousParticle)||~exist('Amp', 'var')||(CurrentChannel~=PreviousChannel)
            PreviousParticle=CurrentParticle;
            [Frames,Amp]=PlotParticleTrace(CurrentParticle,Particles{CurrentChannel},Spots{CurrentChannel});
        end
        clf('reset');
        %we'll plot the spot intensity first on the left axis.
        yyaxis left
        hold on
        plot(Frames(Particles{CurrentChannel}(CurrentParticle).FrameApproved),Amp(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-b')
        xlabel('frame')
        ylabel('transcript intensity (a.u.)')     
        plot(Frames(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),Amp(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.k')
        plot(Frames(Frames==CurrentFrame),Amp(Frames==CurrentFrame),'ob')
        hold off

        
        %now we'll plot the input protein intensity on the right-hand axis.
        yyaxis right
        plot(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames,...
            max(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).Fluo,[],2),'r.-')      
         try
            xlim([min(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames),max(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames)])
        catch
%             error('Not sure what happened here. Problem with trace fig x lim. Talk to AR if you see this, please.');
         end
        ylabel('input protein intensity (a.u.)');
        hold on
        plot(Frames(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),Amp(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r')
        hold off
    end

    %(MT, 2018-02-11) Added support for lattice imaging, maybe 
    %temporary - FIX LATER
    if strcmpi(ExperimentType,'lattice')
        FigureTitle={['Particle: ',num2str(CurrentParticle),'/',num2str(numParticles)],...
            ['Frame: ',num2str(CurrentFrame),'/',num2str(numFrames),')'],...
            ['Z: ',num2str(CurrentZ),'/',num2str(ZSlices),', Ch: ',num2str(CurrentChannel)]};
    else
        FigureTitle={['Particle: ',num2str(CurrentParticle),'/',num2str(numParticles)],...
            ['Frame: ',num2str(CurrentFrame),'/',num2str(numFrames), ' (nc',num2str(FrameInfo(CurrentFrame).nc),')'],...
            ['Z: ',num2str(CurrentZ),'/',num2str(ZSlices),', Ch: ',num2str(CurrentChannel)]};
    end
    
    if HideApprovedFlag==1
        FigureTitle=[FigureTitle,', Showing non-flagged particles'];
    elseif HideApprovedFlag==2
        FigureTitle=[FigureTitle,', Showing disapproved particles'];
    end
    title(FigureTitle)
    
    %Define the windows
    figure(Overlay)
    set(gcf,'units', 'normalized', 'position',[0.01, .55, .33, .33]);
    if UseHistoneOverlay 
        figure(HisOverlayFig)
        set(gcf,'units', 'normalized', 'position',[0.01, 0.1, .33, .33]);
    end
    figure(TraceFig);
    set(gcf,'units', 'normalized', 'position',[0.35, 0.55, .5, .3]);
    figure(SnippetFig);
    set(gcf,'units', 'normalized', 'position',[0.355, 0.15, .2/2, .33/2]);
    figure(ZProfileFig);
    set(gcf,'units', 'normalized', 'position',[0.47, 0.15, .2/2, .33/2]);
    
    figure(Overlay)
    if isempty(SkipWaitForButtonPress)
        ct=waitforbuttonpress;
        cc=get(Overlay,'currentcharacter');
        cm=get(gca,'CurrentPoint');
    else
        cc=SkipWaitForButtonPress;
        SkipWaitForButtonPress=[];
    end
        
    
    if cc=='.' & (CurrentFrame < length({Spots{1}.Fits})) %Move forward one frame
        CurrentFrame=CurrentFrame+1;
        ManualZFlag=0;
        %DisplayRange=[];
    elseif (cc==',')&(CurrentFrame>1) %Move backward one frame
        CurrentFrame=CurrentFrame-1;
        ManualZFlag=0;
        %DisplayRange=[];
    elseif (cc=='>')&(CurrentFrame+5<length({Spots{1}.Fits})) %Move forward five frames
        CurrentFrame=CurrentFrame+5;
        ManualZFlag=0;
        %DisplayRange=[];
    elseif (cc=='<')&(CurrentFrame-5>1) %Move backward five frames
        CurrentFrame=CurrentFrame-5;
        ManualZFlag=0;
        %DisplayRange=[];
    elseif (cc=='''')&(CurrentFrame<length({Spots{1}.Fits})) %Move to the next skipped frame
                                                             %within the particle
        
        %This is the total frame range possible for this particle. Note
        %that we could still want to add spots at the beginning or end of
        %this range.
        FrameRange=Particles{CurrentChannel}(CurrentParticle).Frame(1):...
            Particles{CurrentChannel}(CurrentParticle).Frame(end);
        %Frames in FrameRange that were not in this particle.
        SkippedFrames=FrameRange(~ismember(FrameRange,...
            Particles{CurrentChannel}(CurrentParticle).Frame));
        %Find the next skipped frame and set it up
        CurrentFrame=min(SkippedFrames(SkippedFrames>CurrentFrame));
        
        %If there is no next empty frame, jump to the last frame of this
        %particle
        if isempty(CurrentFrame)
            CurrentFrame=Particles{CurrentChannel}(CurrentParticle).Frame(end);
        end
    elseif (cc==';')&(CurrentFrame>1) %Move to the previous skipped frame
                                      %within the particle
                                      
        %This is the total frame range possible for this particle. Note
        %that we could still want to add spots at the beginning or end of
        %this range.
        FrameRange=Particles{CurrentChannel}(CurrentParticle).Frame(1):...
            Particles{CurrentChannel}(CurrentParticle).Frame(end);
        %Frames in FrameRange that were not in this particle.
        SkippedFrames=FrameRange(~ismember(FrameRange,...
            Particles{CurrentChannel}(CurrentParticle).Frame));
        %Find the next skipped frame and set it up
        CurrentFrame=max(SkippedFrames(SkippedFrames<CurrentFrame));
        
        %If there is no next empty frame, jump to the last frame of this
        %particle
        if isempty(CurrentFrame)
            CurrentFrame=Particles{CurrentChannel}(CurrentParticle).Frame(1);
        end
    elseif (cc=='a')&(CurrentZ<ZSlices) %Move up in Z
        CurrentZ=CurrentZ+1;
        ManualZFlag=1;
    elseif (cc=='z')&(CurrentZ>1) %Move down in Z
        CurrentZ=CurrentZ-1;
        ManualZFlag=1;
    elseif cc=='j'
        try
            iJump= inputdlg('Frame to jump to:',...
                'Move to frame');
            iJump=str2double(iJump{1});           
        catch
            iJump=CurrentFrame;
        end
        if (floor(iJump)>0)&&(iJump<length({Spots{1}.Fits}))
            CurrentFrame=iJump;
            ManualZFlag=0;
        end
        DisplayRange=[];
    elseif cc=='k'
        try
            ParticleJump = inputdlg('Particle to jump to:',...
                'Move to particle');
            ParticleJump=str2double(ParticleJump{1});
        catch
            ParticleJump=CurrentParticle;
        end
        if floor(ParticleJump)>0 && ParticleJump<=numParticles
            CurrentParticle=ParticleJump;
            CurrentFrame=Particles{CurrentChannel}(CurrentParticle).Frame(1);
            ManualZFlag=0;
        end
        
        DisplayRange=[];
    elseif cc=='g'      %Increase histone channel contrast
        if isempty(DisplayRange)
            DisplayRange=[min(min(ImageHis)),max(max(ImageHis))/1.5];
        else
            DisplayRange=[DisplayRange(1),DisplayRange(2)/1.5];
        end
        
    elseif cc=='b'      %Decrease histone channel contrast
        DisplayRange=[min(min(ImageHis)),max(max(ImageHis))*1.5];
    
    elseif cc=='#' %remove a spot from Spots and erase its frame in Particles
        %Check that we're in zoom mode. If not, set it up.
        if ~(ZoomMode || GlobalZoomMode)
            disp('You need to be in Zoom Mode to do this. You can switch using ''o'' or ''+''. Run the ''#'' command again.')
        else           
            %delete from particles
            del = 0;
            choice = questdlg('Are you sure you want to delete this spot? This can''t be undone.', ...
                '', 'Delete spot','Cancel','Cancel');
            switch choice
                case 'Delete spot'
                    disp 'Deleting spot.'
                    del = 1;
                case 'Cancel'
                    disp 'Spot deletion cancelled.'
                    del = 0;
            end
            
            if del
                ind = Particles{CurrentChannel}(CurrentParticle).Index(CurrentFrameWithinParticle);
                onlyFrame = length(Particles{CurrentChannel}(CurrentParticle).Frame) == 1;
                if onlyFrame
                    Particles{CurrentChannel}(CurrentParticle) = [];
                    numParticles = numParticles - 1;
                else
                    particleFields = fieldnames(Particles{CurrentChannel});
                    for i = 1:numel(particleFields)
                        if ~strcmpi(particleFields{i},'Nucleus') && ~strcmpi(particleFields{i},'Approved')
                            Particles{CurrentChannel}(CurrentParticle).(particleFields{i})(CurrentFrameWithinParticle) = [];
                        end
                    end
                end
                    %and this part changes the the index of other particles
                    %in the frame. 
                    for i=1:length(Particles{CurrentChannel})
                        for j = 1:length(Particles{CurrentChannel}(i).Frame)
                            if Particles{CurrentChannel}(i).Frame(j) == CurrentFrame
                                if Particles{CurrentChannel}(i).Index(j) > ind
                                    Particles{CurrentChannel}(i).Index(j) = Particles{CurrentChannel}(i).Index(j) - 1;
                                end
                            end
                        end
                    end
                %and this part deletes from the spots structure.
                CurrentSpot = CurrentParticleIndex; %renaming this to make it clear what it actually is
                Spots{CurrentChannel}(CurrentFrame).Fits(CurrentSpot)= [];
                %now delete from spotfilter
               spotRow = SpotFilter{CurrentChannel}(CurrentFrame,:);
               spotRow(CurrentSpot) = [];
               spotRow(end+1) = NaN;
               try
                SpotFilter{CurrentChannel}(CurrentFrame,:) = spotRow;
               catch
                   error('There probably wasn''t a spot in the frame you were trying to delete.')
               end
                if onlyFrame
                    %switch to another particle just to avoid any potential weirdness with
                    %checkparticletracking refreshing. simpler version of the
                    %'m' button
                    NextParticle = CurrentParticle+1;
                    if NextParticle>numParticles
                        NextParticle=NextParticle-2; %go backwards one particle if the deleted particle was the last. 
                    end
                    CurrentParticle=NextParticle;
                    CurrentFrame=Particles{CurrentChannel}(CurrentParticle).Frame(1);
                    ParticleToFollow=[];
                    DisplayRange=[];
                    disp 'Spot deleted successfully. Trace figures will refresh after switching particles.' 
                end
            end
            ZoomMode=0;
            GlobalZoomMode=0;
        end
    
%     note to ar: a potentially simpler version of this button deletes
%     particle frame but not spot. implement that with a different button.
%           
    elseif cc=='[' %Add particle and all of its shadows to Spots.
        
        %Check that we're in zoom mode. If not, set it up.
        if ~(ZoomMode || GlobalZoomMode)
            disp('You need to be in Zoom Mode to do this. You can switch using ''o'' or ''+''. Run the ''['' command again.')
        else
            %Click on the region we're going to fit in order to find a new
            %spot and add that spot to the current particle
            
            %Check that this particle doesn't already have a spot assigned
            %in this frame
            if sum(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame) &&  ~GlobalZoomMode
                warning('There is a spot assigned to this particle in this frame already.')
            else
                
                [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'r', 'linewidth',1);
                ConnectPositionx = round(ConnectPositionx);
                ConnectPositiony = round(ConnectPositiony);
                
                pixelSize = FrameInfo(1).PixelSize*1000; %nm
                LinesPerFrame = FrameInfo(1).LinesPerFrame;
                PixelsPerLine = FrameInfo(1).PixelsPerLine;
                snippet_size = 2*(floor(1300/(2*pixelSize))) + 1; % nm. note that this is forced to be odd
                
                % check that the clicked particle isn't too close to the
                % edge of the frame
                if (ConnectPositionx > snippet_size/2) && (ConnectPositionx + snippet_size/2 < PixelsPerLine)...
                        && (ConnectPositiony > snippet_size/2) && (ConnectPositiony + snippet_size/2 < LinesPerFrame)
                    SpotsIndex = length(Spots{CurrentChannel}(CurrentFrame).Fits)+1;
                    breakflag = 0;
                    parfor i = 2:ZSlices-1
                        spotsIm=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                             FilePrefix,iIndex(CurrentFrame,NDigits),'_z',iIndex(i,2),nameSuffix,'.tif']);  
                        Threshold = min(min(spotsIm));
                        dog = spotsIm;
                        im_thresh = dog >= Threshold;
                        [im_label, ~] = bwlabel(im_thresh);
                        microscope = FrameInfo(1).FileMode;
                        show_status = 0;
                        fig = [];
                        k = 1; %This is supposed to be the index for the partiles in an image.
                               %However, this image only contains one particle
                        neighborhood = round(1300 / pixelSize); %nm
                        %Get the information about the spot on this z-slice
                        temp_particles{i} = identifySingleSpot(k, spotsIm, im_label, dog, neighborhood, snippet_size, ...
                            pixelSize, show_status, fig, microscope, [1, ConnectPositionx, ConnectPositiony], [], '' );
                    end

                    for i = 2:ZSlices-1
                        if ~isempty(temp_particles{i})
                            %Copy the information stored on temp_particles into the
                            %Spots structure                            
                            if ~isempty(temp_particles{i}{1})
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).FixedAreaIntensity(i)=...
                                temp_particles{i}{1}{1};
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).xFit(i)=...
                                    temp_particles{i}{1}{2};
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).yFit(i)=...
                                    temp_particles{i}{1}{3};
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).Offset(i)=...
                                    temp_particles{i}{1}{4};
%                                 Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).Snippet{i}=...
%                                     temp_particles{i}{1}{5};
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).Area{i}=...
                                    temp_particles{i}{1}{6};
%                                     Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).xFitWidth{i}=...
%                                         temp_particles{i}{1}{7};
%                                     Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).yFitWidth{i}=...
%                                         temp_particles{i}{1}{8};
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).yDoG(i)=...
                                    temp_particles{i}{1}{9};
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).xDoG(i)=...
                                    temp_particles{i}{1}{10};
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).GaussianIntensity(i)=...
                                    temp_particles{i}{1}{11};
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).CentralIntensity(i)=...
                                    temp_particles{i}{1}{12};
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).DOGIntensity(i)=...
                                    temp_particles{i}{1}{13};
%                                     Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).snippet_mask{i}=...
%                                         NaN;
%                                     Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).SisterDistance(i)=... 
%                                         temp_particles{i}{1}{17};
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).ConfidenceIntervals{i}=...
                                    temp_particles{i}{1}{19};
%                                 Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).gaussSpot{i}=...
%                                     temp_particles{i}{1}{20};
%                                 Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).rawSpot{i}=...
%                                     temp_particles{i}{1}{21};
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).z(i)=...
                                    i;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).gaussParams{i}=...
                                    temp_particles{i}{1}{22};
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).discardThis=...
                                    0;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).frame=...
                                    CurrentFrame;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).r=...
                                    0;
                            else
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).FixedAreaIntensity(i)=...
                                    nan;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).xFit(i)=...
                                    nan;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).yFit(i)=...
                                    nan;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).Offset(i)=...
                                    nan;
%                                     Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).Snippet{i}=...
%                                         nan;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).Area{i}=...
                                    nan;
%                                     Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).xFitWidth{i}=...
%                                         nan;
%                                     Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).yFitWidth{i}=...
%                                         nan;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).yDoG(i)=...
                                    nan;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).xDoG(i)=...
                                    nan;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).GaussianIntensity(i)=...
                                    nan;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).CentralIntensity(i)=...
                                    nan;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).DOGIntensity(i)=...
                                    nan;
%                                     Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).snippet_mask{i}=...
%                                         nan;
%                                     Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).SisterDistance(i)=... 
%                                         nan;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).ConfidenceIntervals{i}=...
                                    nan;
%                                 Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).gaussSpot{i}=...
%                                     nan;
%                                 Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).rawSpot{i}=...
%                                     nan;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).z(i)=...
                                    i;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).gaussParams{i}=...
                                    nan;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).discardThis=...
                                    0;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).frame=...
                                    CurrentFrame;
                                Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).r=...
                                    0;  
                            end
%                             catch
%                                 error('This spot might be too close to the image boundary');
%                             end

                        else
                            disp('No spot added. Did you click too close to the image boundary?')
                            breakflag = 1;
                            break
                        end                                  
                    end
                                                    
                    if ~breakflag
                        %pull intensity value from particle snippets
                        RawIntensityVec = [Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).FixedAreaIntensity];            
                        CentralIntensityVec = [Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).CentralIntensity]; 
                        %find slice with brightest pixel
                        [~, MaxIndexCentral] = max(CentralIntensityVec);            
                        z_vec = [Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).z]; %convenience vector
%                         %Determine which centering method was used for this
%                         if isfield(Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex),'IntegralZ')
%                             IntegralZ_flag = Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).IntegralZ;
%                         else % assume set was generated prior to inclusion of integral option
%                             IntegralZ_flag = 0;
%                         end                        
                        % calculate concenience vectors
                        z_grid = min(z_vec):max(z_vec);
                        z_raw_values = NaN(size(z_grid));            
                        z_raw_values(ismember(z_grid,z_vec)) = RawIntensityVec;
                        if ~IntegralZ_flag                
                            CentralZ = z_vec(MaxIndexCentral); 
                            ZStackIndex = MaxIndexCentral;
                        else
                            % Convolve with gaussian filter to find "best" center
                            g = [-1 0 1];
                            gaussFilter = exp(-g .^ 2 / (2 ));
                            RawRefVec = conv(gaussFilter,z_raw_values);
                            RawRefVec = RawRefVec(2:end-1);
                            RawRefVec(1) = NaN;
                            RawRefVec(end) = NaN;
                            RawRefVec = RawRefVec(ismember(z_grid,z_vec));
                            [~, MaxIndexIntegral] = max(RawRefVec);
                            CentralZ = z_vec(MaxIndexIntegral);               
                            ZStackIndex = MaxIndexIntegral;
                        end            
                        
                        % if there are insufficient slices, these metrics will register as NaNs
                        RawIntegral3 = mean(z_raw_values(ismember(z_grid,CentralZ-1:CentralZ+1)));
                        RawIntegral5 = mean(z_raw_values(ismember(z_grid,CentralZ-2:CentralZ+2)));
                        if multi_slice_flag
                            Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).FixedAreaIntensity3 = RawIntegral3;
                            Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).FixedAreaIntensity5 = RawIntegral5;
                            Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).IntegralZ = IntegralZ_flag; 
                            Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).snippet_size = snippet_size;
                        end
                        Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).brightestZ = CentralZ;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).r = NaN;
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).discardThis = 0;
%                                                                         
                        %Add this to SpotFilter, which tells the code that this spot is
                        %above the threshold. First, check whether the
                        %dimensions of SpotFilter need to be altered. If so, pad it with NaNs
                        if size(SpotFilter{CurrentChannel},2)>SpotsIndex
                            SpotFilter{CurrentChannel}(CurrentFrame,SpotsIndex)=1;
                        else
                            %Pad with NaNs
                            SpotFilter{CurrentChannel}(:,end:SpotsIndex)=NaN;
                            SpotFilter{CurrentChannel}(CurrentFrame,SpotsIndex)=1;                
                        end

%                     if ~breakflag
%                         %Find the maximum Z-plane
%                         [~, max_index] = max(Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).CentralIntensity);
%                         Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).brightestZ = Spots{CurrentChannel}(CurrentFrame).Fits(SpotsIndex).z(max_index);
% 
%                         %Add this to SpotFilter, which tells the code that this spot is
%                         %above the threshold. First, check whether the
%                         %dimensions of SpotFilter need to be altered. If so, pad it with NaNs
%                         if size(SpotFilter{CurrentChannel},2)>SpotsIndex
%                             SpotFilter{CurrentChannel}(CurrentFrame,SpotsIndex)=1;
%                         else
%                             %Pad with NaNs
%                             SpotFilter{CurrentChannel}(:,end:SpotsIndex)=NaN;
%                             SpotFilter{CurrentChannel}(CurrentFrame,SpotsIndex)=1;                
%                         end
% 
                        %Turn this spot into a new particle. This is the equivalent of
                        %the 'u' command.
                        [SpotFilter{CurrentChannel},Particles{CurrentChannel}]=...
                            TransferParticle(Spots{CurrentChannel},...
                            SpotFilter{CurrentChannel},Particles{CurrentChannel},...
                            CurrentFrame,SpotsIndex);
                        numParticles = numParticles + 1;

                        %Connect this particle to the CurrentParticle. This is
                        %the equivalent of running the 'c' command. 
                        if ~GlobalZoomMode
                            Particles{CurrentChannel}=...
                            JoinParticleTraces(CurrentParticle,...
                            numParticles,Particles{CurrentChannel});
                        else
                            disp('Re-run TrackmRNADynamics to associate this particle with a nucleus and trace.')
                        end

                        %Finally, force the code to recalculate the fluorescence trace
                        %for this particle
                        PreviousParticle=0;
                        disp('Spot addded to the current particle.')                    
                    else 
                        warning('You clicked too close to the edge. A spot can''t be added here.');
                    end
                end
            end
        end

    elseif cc=='r'
        %Order particles by the earliest frame they appear at. This makes the
        %tracking a lot easier!
        clear FirstFrame
        for i=1:numParticles
            FirstFrame(i)=Particles{CurrentChannel}(i).Frame(1);
        end
        [~,Permutations]=sort(FirstFrame);
        Particles{CurrentChannel}=Particles{CurrentChannel}(Permutations);
        
    elseif cc=='f'
        Answer=lower(input('Are you sure you want to redo the tracking?  (y/n) ','s'));
        if Answer=='y'
            warning('HG: Not clear that this feature will work with the multiple channels')
            
            %We need to save the data
            save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo')
            if UseHistoneOverlay
                save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2', '-v7.3')
                save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'],'schnitzcells', '-v7.3')
            else
                save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2', '-v7.3')            
            end
        disp('Particles saved.')
        if NChannels==1
            Particles=Particles{1};
        end
            
           [Particles,schnitzcells,fad,fad2]=TrackmRNADynamics(FilePrefix(1:end-1),...
               Threshold1,Threshold2); 
        if NChannels==1
            Particles={Particles};
        end
           %Check the FrameApproved field
            for i=1:numParticles
                if isempty(Particles{CurrentChannel}(i).FrameApproved)
                    Particles{CurrentChannel}(i).FrameApproved=true(size(Particles{CurrentChannel}(i).Frame));
                end
            end
        end
    elseif cc=='c'
        PreviousParticle=0;
        if ~sum(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)
            %ConnectPosition=ginput(1);
            
            [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'b', 'linewidth',1);
            ConnectPosition = [ConnectPositionx,ConnectPositiony];
            
            if ~isempty(ConnectPosition)
                [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,Spots{CurrentChannel},Particles{CurrentChannel});
               
                
                %Check that the clicked particle doesn't exist in a previous
                %frame, that there is no overlap of frames. If it does
                %exist in a previous frame we will have to disconnect it.
                if sum(Particles{CurrentChannel}(ParticleOutput).Frame<CurrentFrame)
                    %Disconnect the clicked particle
                    Particles{CurrentChannel}=SeparateParticleTraces(ParticleOutput,CurrentFrame,Particles{CurrentChannel});
                    ParticleOutput=ParticleOutput+1;
                    
                    %If the current particle has an index larger than that
                    %of the clicked particle (ParticleOutput) we also need to
                    %move the index of the current Particle by one.
                    if ParticleOutput<CurrentParticle
                    	CurrentParticle=CurrentParticle+1;
                    end
                end
                
                %Check that there is no overlap. If so, split current particle
                overlap=0;
                for i=1:length(Particles{CurrentChannel}(ParticleOutput).Frame)
                    for j=1:length(Particles{CurrentChannel}(CurrentParticle).Frame)
                        if Particles{CurrentChannel}(ParticleOutput).Frame(i)==Particles{CurrentChannel}(CurrentParticle).Frame(j)
                            overlap=1;
                        end
                    end
                end
               
                if overlap
                    %Disconnect the clicked particle
                    Particles{CurrentChannel}=SeparateParticleTraces(CurrentParticle,CurrentFrame,Particles{CurrentChannel});
                    
                    %If the clicked particle has an index larger than that
                    %of the current particle we also need to
                    %move the index of the clicked particle by one.
                    if ParticleOutput>CurrentParticle
                    	ParticleOutput=ParticleOutput+1;
                    end
                end
                
                
                    
                Particles{CurrentChannel}=JoinParticleTraces(CurrentParticle,ParticleOutput,Particles{CurrentChannel});
                %Deals with the indexing changing because of the removal of
                %the old particle.
                 if ParticleOutput<CurrentParticle
                     CurrentParticle=CurrentParticle-1;
                 end
                %Sort the frames within the particle. This is useful if we
                %connected to a particle that came before.
                [SortedFrame,Permutations]=sort(Particles{CurrentChannel}(CurrentParticle).Frame);
                Particles{CurrentChannel}(CurrentParticle).Frame=Particles{CurrentChannel}(CurrentParticle).Frame(Permutations);
                Particles{CurrentChannel}(CurrentParticle).Index=Particles{CurrentChannel}(CurrentParticle).Index(Permutations);
                Particles{CurrentChannel}(CurrentParticle).FrameApproved=Particles{CurrentChannel}(CurrentParticle).FrameApproved(Permutations);                

            end
            
        else
            ConnectPosition=ginput(1);
            [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,Spots{CurrentChannel},Particles{CurrentChannel});
            
            %If it's an independent particle swap it with the frame in the
            %current particle
            if (length(Particles{CurrentChannel}(ParticleOutput).Frame)==1)&&...
                    (sum(Particles{CurrentChannel}(ParticleOutput).Frame==CurrentFrame)==1)
                
                ParticleTemp=Particles{CurrentChannel}(ParticleOutput);
                
                %Copy the particle out
                Particles{CurrentChannel}(ParticleOutput).Index=...
                    Particles{CurrentChannel}(CurrentParticle).Index(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame);

                %Copy the new particle in
                Particles{CurrentChannel}(CurrentParticle).Index(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)=...
                    ParticleTemp.Index;
                Particles{CurrentChannel}(CurrentParticle).FrameApproved(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)=1;
            else
                disp('Cannnot connect to two particles!')
            end
            
            
        end
    elseif cc=='p' %Identify a particle. It will also tell you the particle associated with
                   %  the clicked nucleus.
        [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'b', 'linewidth',1);
        ConnectPosition = [ConnectPositionx,ConnectPositiony];
        if ~isempty(ConnectPosition)
            %Find the closest particle
            [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,Spots{CurrentChannel},Particles{CurrentChannel});
            disp(['Clicked particle: ',num2str(ParticleOutput)]);
            
            if UseHistoneOverlay
                %Find the closest nucleus
                NewNuclei=ConnectPosition;

                %Find which schnitz this corresponds to
                SchnitzSuspect=[];
                xPosSuspect=[];
                yPosSuspect=[];
                for j=1:length(schnitzcells)
                    if sum(schnitzcells(j).frames==CurrentFrame)
                        SchnitzSuspect=[SchnitzSuspect,j];
                        xPosSuspect=[xPosSuspect,...
                            schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))];
                        yPosSuspect=[yPosSuspect,...
                            schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))];
                    end
                end

                %Find the closest one to the point where we clicked
                Distance=sqrt((NewNuclei(1)-xPosSuspect).^2+(NewNuclei(2)-yPosSuspect).^2);
                [MinValue,ClosestNucleusIndex]=min(Distance);

                ClickedSchnitz=SchnitzSuspect(ClosestNucleusIndex);

                %Now, find its associated particle
                for i=1:numParticles
                    if ~isempty(Particles{CurrentChannel}(i).Nucleus)
                        AssignedNuclei(i)=Particles{CurrentChannel}(i).Nucleus;
                    else
                        AssignedNuclei(i)=nan;
                    end
                end
                AssociatedParticle=find(AssignedNuclei==ClickedSchnitz);

                if isempty(AssociatedParticle)
                    disp(['Nucleus ',num2str(ClickedSchnitz),' does not have an associated particle'])
                else
                    disp(['Particle ',num2str(AssociatedParticle),' is associate with nucleus ',...
                        num2str(ClickedSchnitz)])
                end
            end
        end
        
    elseif cc=='\' %Moves to clicked particle.
        [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'b', 'linewidth',1);
        ConnectPosition = [ConnectPositionx,ConnectPositiony];
        if ~isempty(ConnectPosition)
            %Find the closest particle
            [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,Spots{CurrentChannel},Particles{CurrentChannel});
            disp(['Clicked particle: ',num2str(ParticleOutput)]);
            try
                ParticleJump=ParticleOutput;
            end
            if (floor(ParticleJump)>0)&&(ParticleJump<=numParticles)
                CurrentParticle=ParticleJump;
                CurrentFrame=Particles{CurrentChannel}(CurrentParticle).Frame(1);
                ManualZFlag=0;
            end
            
            if UseHistoneOverlay
                %Find the closest nucleus
                NewNuclei=ConnectPosition;

                %Find which schnitz this corresponds to
                SchnitzSuspect=[];
                xPosSuspect=[];
                yPosSuspect=[];
                for j=1:length(schnitzcells)
                    if sum(schnitzcells(j).frames==CurrentFrame)
                        SchnitzSuspect=[SchnitzSuspect,j];
                        xPosSuspect=[xPosSuspect,...
                            schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))];
                        yPosSuspect=[yPosSuspect,...
                            schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))];
                    end
                end

                %Find the closest one to the point where we clicked
                Distance=sqrt((NewNuclei(1)-xPosSuspect).^2+(NewNuclei(2)-yPosSuspect).^2);
                [MinValue,ClosestNucleusIndex]=min(Distance);

                ClickedSchnitz=SchnitzSuspect(ClosestNucleusIndex);

                %Now, find its associated particle
                for i=1:numParticles
                    if ~isempty(Particles{CurrentChannel}(i).Nucleus)
                        AssignedNuclei(i)=Particles{CurrentChannel}(i).Nucleus;
                    else
                        AssignedNuclei(i)=nan;
                    end
                end
                AssociatedParticle=find(AssignedNuclei==ClickedSchnitz);

                if isempty(AssociatedParticle)
                    disp(['Nucleus ',num2str(ClickedSchnitz),' does not have an associated particle'])
                else
                    disp(['Particle ',num2str(AssociatedParticle),' is associated with nucleus ',...
                        num2str(ClickedSchnitz)])
                end
            end
            
        end    
    elseif cc=='u'
        [x2,y2,z2]=SpotsXYZ(Spots{CurrentChannel}(CurrentFrame));
        if ~isempty(x2)
            ClickedSpot=ginput(1);

            UnfilterSpot(Spots{CurrentChannel},SpotFilter{CurrentChannel},...
                ClickedSpot,Particles{CurrentChannel},CurrentFrame)
         end
    elseif cc=='i'
%         PreviousParticle=0;
%          [x2,y2]=fad2xyzFit(CurrentFrame,fad2(CurrentChannel), 'addMargin'); 
%          if ~isempty(x2)
%             fad2Position=ginput(1);
%             if (~isempty(fad2Position))
%                 [fad(CurrentChannel),fad2(CurrentChannel),Particles{CurrentChannel}]=...
%                     Integratefad2Particle(fad(CurrentChannel),fad2(CurrentChannel),fad2Position,Particles{CurrentChannel},CurrentFrame);
%             end
%          end
% 
%          
%          if (~sum(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame))&(~isempty(fad2Position))
%             ConnectPosition=fad2Position;
%             [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,fad(CurrentChannel),Particles{CurrentChannel});
%             
%             %Check that the clicked particle doesn't exist in a previous
%             %frame, that there is no overlap of frames.  Maybe I can have
%             %those in a different color.
% 
%            
%             if sum(Particles{CurrentChannel}(ParticleOutput).Frame<CurrentFrame)
%                 disp(['Target particle (',num2str(ParticleOutput),') is already in a previous frame!']);
%             else
%                 Particles=JoinParticleTraces(CurrentParticle,ParticleOutput,Particles{CurrentChannel});
%                 %Do this in case the clicked particle comes before the current
%                 %particle in the structure
%                 if ParticleOutput<CurrentParticle
%                     CurrentParticle=ParticleOutput;
%                 end
%                 %
%                 %Sort the frames within the particle. This is useful if we
%                 %connected to a particle that came before.
%                 
%                 %There is an error here, cell contents reference from a non-cell array object.
%                 [SortedFrame,Permutations]=sort(Particles{CurrentChannel}(CurrentParticle).Frame);
%                 Particles{CurrentChannel}(CurrentParticle).Frame=Particles{CurrentChannel}(CurrentParticle).Frame(Permutations);
%                 Particles{CurrentChannel}(CurrentParticle).Index=Particles{CurrentChannel}(CurrentParticle).Index(Permutations);
%                 Particles{CurrentChannel}(CurrentParticle).FrameApproved=Particles{CurrentChannel}(CurrentParticle).FrameApproved(Permutations);
%                 
%                 if UseHistoneOverlay
%                     %Check for consistency within schnitzcell 
%                     [Particles{CurrentChannel},schnitzcells]=CheckSchnitzLineage(Particles{CurrentChannel},CurrentParticle,schnitzcells,CurrentFrame,...
%                         Overlay);
%                 end  
%             end 
%         else
%             disp('Cannnot connect to two particles!')
%         end

        warning(' AR 1/15/18: This is currently deprecated. Talk to HG if you need this function.')

    
     elseif cc=='d'  %d Separate traces forward at the current frame.
         %The separated particle (the trace following current frame) won't have a nucleus assigned!
         PreviousParticle=0;
        %Check that the particle does actually exist in this frame
        if ~(Particles{CurrentChannel}(CurrentParticle).Frame(1)==CurrentFrame)
            if sum(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)
                Particles{CurrentChannel}=SeparateParticleTraces(CurrentParticle,CurrentFrame,Particles{CurrentChannel});
            end
        elseif length(Particles{CurrentChannel}(CurrentParticle).Frame)==1
            Particles{CurrentChannel}(CurrentParticle).Nucleus=[];
        else
            disp('Cannot divide a trace at the first time point')
        end
        
    elseif cc=='v'  %d Separate traces forward at the current frame.
         %The separated particle (the trace following current frame) won't have a nucleus assigned!
         PreviousParticle=0;
        %Check that the particle does actually exist in this frame
        if ~(Particles{CurrentChannel}(CurrentParticle).Frame(1)==CurrentFrame)
            if sum(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)
                Particles{CurrentChannel}=SeparateParticleTraces(CurrentParticle,CurrentFrame,Particles{CurrentChannel});
            end
        elseif length(Particles{CurrentChannel}(CurrentParticle).Frame)==1
            Particles{CurrentChannel}(CurrentParticle).Nucleus=[];
        else
            disp('Cannot divide a trace at the first time point')
        end        
        
        
        
        
    elseif cc=='q'      %Approve a trace
        if Particles{CurrentChannel}(CurrentParticle).Approved==1
            Particles{CurrentChannel}(CurrentParticle).Approved=2;
        elseif Particles{CurrentChannel}(CurrentParticle).Approved==0
            Particles{CurrentChannel}(CurrentParticle).Approved=1;
        elseif Particles{CurrentChannel}(CurrentParticle).Approved==2
            Particles{CurrentChannel}(CurrentParticle).Approved=0;
        end
    elseif cc=='w'      %Disapprove a trace
        if Particles{CurrentChannel}(CurrentParticle).Approved==-1
            Particles{CurrentChannel}(CurrentParticle).Approved=0;
        else
            Particles{CurrentChannel}(CurrentParticle).Approved=-1;
        end    
        
    elseif cc=='s'
        
        %If we only have one channel bring Particles back to the legacy
        %format without any cells
        if NChannels==1
            Particles=Particles{1};
            Spots=Spots{1};
            SpotFilter=SpotFilter{1};
        end

        
        
        save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo')
        if UseHistoneOverlay
            save([DataFolder,filesep,'Particles.mat'],'Particles','SpotFilter','Threshold1','Threshold2', '-v7.3')
            save([DataFolder,filesep,'Spots.mat'],'Spots', '-v7.3') %CS20170912 necessary for saving Spots.mat if >2GB
            save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'],'schnitzcells', '-v7.3')
        else
            save([DataFolder,filesep,'Particles.mat'],'Particles','SpotFilter','Threshold1','Threshold2', '-v7.3')            
            save([DataFolder,filesep,'Spots.mat'],'Spots','-v7.3') %CS20170912 necessary for saving Spots.mat if >2GB
        end
        disp('Particles saved.')
        if NChannels==1
            Particles={Particles};
            Spots = {Spots};
            SpotFilter = {SpotFilter};           
        end
        
    elseif cc=='t'
        ShowThreshold2=~ShowThreshold2;
    elseif (cc=='y')&(~UseHistoneOverlay)
            FrameInfo=DetermineNC(fad,Particles{CurrentChannel},FrameInfo);
    elseif cc=='h'
        if HideApprovedFlag==0
            HideApprovedFlag=1;         %Show only non-approved traces
        elseif HideApprovedFlag==1
            HideApprovedFlag=2;         %Show only yellow and red traces
        elseif HideApprovedFlag==2
            HideApprovedFlag=0;
        end

        %HideApprovedFlag=~HideApprovedFlag;
    elseif cc=='o'
        if ~GlobalZoomMode
            ZoomMode=~ZoomMode;
        else
            disp('Try again after exiting global zoom mode by hitting ''+''') 
        end
    
    elseif cc=='+'
        if ~ZoomMode
            if ~GlobalZoomMode
                [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'r', 'linewidth',1);
                xForZoom = round(ConnectPositionx);
                yForZoom = round(ConnectPositiony);
            else
            end
            GlobalZoomMode =~ GlobalZoomMode;
        else
            disp('Try again after exiting zoom mode by hitting ''o''')
        end
    
    elseif (cc=='m')&(CurrentParticle<numParticles)
        
        NextParticle=CurrentParticle+1;
        
        if NextParticle>numParticles
            NextParticle=numParticles;
        end
        
        
        %Mode 1 - skip approved or flagged traces
        while (HideApprovedFlag)==1&&(NextParticle<numParticles)&&...
                ((Particles{CurrentChannel}(NextParticle).Approved==1)||(Particles{CurrentChannel}(NextParticle).Approved==-1)||...
                (Particles{CurrentChannel}(NextParticle).Approved==2))
                NextParticle=NextParticle+1;
        end
        
        %Mode 2 - skip approved traces
        while ((HideApprovedFlag)==2)&&(NextParticle<numParticles)&&...
                ((Particles{CurrentChannel}(NextParticle).Approved==1)||(Particles{CurrentChannel}(NextParticle).Approved==2))
                NextParticle=NextParticle+1;
        end
        
        
        CurrentParticle=NextParticle;
%         if ~isempty(Particles(CurrentParticle).Nucleus)
%             CurrentFrame=schnitzcells(Particles(CurrentParticle).Nucleus).frames(1)-1;
%         else
%             CurrentFrame=Particles(CurrentParticle).Frame(1);
%         end
        CurrentFrame=Particles{CurrentChannel}(CurrentParticle).Frame(1);

        ParticleToFollow=[];
        DisplayRange=[];
        
        msg = Particles{CurrentChannel}(CurrentParticle).Frame(find(diff(Particles{CurrentChannel}(CurrentParticle).Frame)>1));
        
        if ~isempty(msg)
%             disp('Missing frames:') %AR 12/3/17- Not sure what this
%             message is trying to say, so I am silencing it for now. 
%             msg        
        else 
            %do nothing
        end
            
    elseif (cc=='n')&(CurrentParticle>1)
        Approved=(find([Particles{CurrentChannel}.Approved]));
        %NotApproved=(find(~[Particles.Approved]));
        
        NextParticle=CurrentParticle-1;
        
        
        
        %Mode 1 - show non-flagged traces
        while (HideApprovedFlag)==1 && (NextParticle>1) &&...
                ((Particles{CurrentChannel}(NextParticle).Approved==1) || (Particles{CurrentChannel}(NextParticle).Approved==-1) ||...
                (Particles{CurrentChannel}(NextParticle).Approved==2))
            NextParticle=NextParticle-1;
            if NextParticle<1
                NextParticle=1;
            end
        end
        
        
        %Mode 2 - show disapproved traces
        while ((HideApprovedFlag)==2)&&(NextParticle>1)&&...
                ((Particles{CurrentChannel}(NextParticle).Approved==1)||(Particles{CurrentChannel}(NextParticle).Approved==2))
            NextParticle=NextParticle-1;
            if NextParticle<1
                NextParticle=1;
            end
        end
        
        
        if NextParticle<1
            NextParticle=CurrentParticle;
        end
        
        CurrentParticle=NextParticle;
        
%         if ~isempty(Particles(CurrentParticle).Nucleus)
%             CurrentFrame=schnitzcells(Particles(CurrentParticle).Nucleus).frames(1)-1;
%         else
%             CurrentFrame=Particles(CurrentParticle).Frame(1);
%         end
        CurrentFrame=Particles{CurrentChannel}(CurrentParticle).Frame(1);

        
        ParticleToFollow=[];
        
        DisplayRange=[];
        
    elseif cc=='e'
        Particles{CurrentChannel}(CurrentParticle).FrameApproved(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)=...
            ~Particles{CurrentChannel}(CurrentParticle).FrameApproved(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame);
    
    
    
    %Schnitzcells specific
    
    elseif cc=='l' %Split a nucleus and select one or two daughter cells or stop the lineage
        PreviousParticle=0;
        
        disp('Select one/two daughter nuclei and press ENTER or just press ENTER to terminate lineage')
        NewNuclei=ginput(2);
        
        if isempty(NewNuclei)
            warning('Write the disconnect schnitz part')
        else
            [NClicks,~]=size(NewNuclei);
        
            ClickedSchnitz=[];
            %Find the nuclei/schnitz that we clicked on
            for i=1:NClicks
                %Find which schnitz this corresponds to
                SchnitzSuspect=[];
                xPosSuspect=[];
                yPosSuspect=[];
                for j=1:length(schnitzcells)
                    if sum(schnitzcells(j).frames==CurrentFrame)
                        SchnitzSuspect=[SchnitzSuspect,j];
                        if (~isempty(schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))))&...
                                (~isempty(schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))))
                            xPosSuspect=[xPosSuspect,...
                                schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))];
                            yPosSuspect=[yPosSuspect,...
                                schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))];
                        else
                            xPosSuspect=[xPosSuspect,inf];
                            yPosSuspect=[yPosSuspect,inf];
                        end
                    end
                end

                %Find the closest one to the point where we clicked
                Distance=sqrt((NewNuclei(i,1)-xPosSuspect).^2+(NewNuclei(i,2)-yPosSuspect).^2);
                [MinValue,ClosestNucleusIndex]=min(Distance);

                ClickedSchnitz(i)=SchnitzSuspect(ClosestNucleusIndex);
            end
            
            %Now look at the different cases
            
            
%       Click on one nucleus + ENTER: Continue the schnitz with that nucleus.
%       Click on two nuclei: Split the current nucleus into two daughter
%       nuclei.
%       Click on the same nucleus twice: Split the current nucleus, but
%       with only one daughter nucleus.
            
            
            if length(ClickedSchnitz)==1
                if Particles{CurrentChannel}(CurrentParticle).Nucleus==ClickedSchnitz %Split the lineage
                    [Particles{CurrentChannel},schnitzcells]=SplitSchnitz(Particles{CurrentChannel},schnitzcells,...
                            CurrentFrame,...
                            CurrentParticle);
                else
                    try
                        [Particles{CurrentChannel},schnitzcells]=...
                            JoinSchnitz(Particles{CurrentChannel},schnitzcells,Particles{CurrentChannel}(CurrentParticle).Nucleus,...  
                            ClickedSchnitz,CurrentFrame);
                    catch
                        disp('Error in JoinSchnitz')
                    end
                end
            elseif length(ClickedSchnitz)==2
                if ClickedSchnitz(1)~=ClickedSchnitz(2)
                    [Particles{CurrentChannel},schnitzcells]=SplitSchnitzDaughters(Particles{CurrentChannel},schnitzcells,...
                            CurrentFrame,...
                            Particles{CurrentChannel}(CurrentParticle).Nucleus,ClickedSchnitz(1),ClickedSchnitz(2));
                else
                    [Particles{CurrentChannel},schnitzcells]=SplitSchnitzDaughters(Particles{CurrentChannel},schnitzcells,...
                        CurrentFrame,...
                        Particles{CurrentChannel}(CurrentParticle).Nucleus,ClickedSchnitz(1),0);
                end
            else
                error('Too many cells selected')
            end
        end
        
        %Check if there are any issues with cellno getting lost in the
        %schnitz
        for SchnitzN=1:length(schnitzcells)
            if length(schnitzcells(SchnitzN).frames)~=length(schnitzcells(SchnitzN).cellno)
                warning(['Problem with schnitz ',num2str(SchnitzN)])
                keyboard
            end
        end
        
    elseif cc=='2' %2 set parent of current nucleus
        disp('Select the mother nucleus or press enter to delete mother information')
        NewNuclei=ginput(1);
        
        if isempty(NewNuclei)
            error('Write the disconnect schnitz part')
        else
            [NClicks,~]=size(NewNuclei);
        
            ClickedSchnitz=[];
            %Find the nuclei/schnitz that we clicked on
            for i=1:NClicks
                %Find which schnitz this corresponds to
                SchnitzSuspect=[];
                xPosSuspect=[];
                yPosSuspect=[];
                for j=1:length(schnitzcells)
                    if sum(schnitzcells(j).frames==CurrentFrame)
                        SchnitzSuspect=[SchnitzSuspect,j];
                        if (~isempty(schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))))&...
                                (~isempty(schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))))
                            xPosSuspect=[xPosSuspect,...
                                schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))];
                            yPosSuspect=[yPosSuspect,...
                                schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))];
                        else
                            xPosSuspect=[xPosSuspect,inf];
                            yPosSuspect=[yPosSuspect,inf];
                        end
                    end
                end

                %Find the closest one to the point where we clicked
                Distance=sqrt((NewNuclei(i,1)-xPosSuspect).^2+(NewNuclei(i,2)-yPosSuspect).^2);
                [MinValue,ClosestNucleusIndex]=min(Distance);

                ClickedSchnitz(i)=SchnitzSuspect(ClosestNucleusIndex);
            end
            
            %Now look at the different cases. Note that I don't have a good
            %way to fix the parent nucleus itself. This might be a bad idea
            %after all
            
       
            
            if length(ClickedSchnitz)==1
                schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).P=ClickedSchnitz;

            elseif isempty(ClickedSchnitz)
                schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).P=0;
            end                
        end    
        
    elseif (cc=='8')&(NChannels>1)      %Switch channels
        %Update the channel number
        PreviousChannel=CurrentChannel;
        CurrentChannel=CurrentChannel+1;
        if CurrentChannel>NChannels
            CurrentChannel=1;
        end
        
        %Update the coatChannel
        coatChannel=coatChannels(CurrentChannel);
        
        %Do we have a histone channel? If so, we can find the particle in
        %the next channel corresponding to this nucleus.
        if UseHistoneOverlay
            %If a particle is associated with this same nucleus in the new
            %channel then change to it
            AssignedNucleusPreviousChannel=Particles{PreviousChannel}(CurrentParticle).Nucleus;

            %Now, find its associated particle
            AssignedNucleusNewChannel=[];
            for i=1:numParticles
                if ~isempty(Particles{CurrentChannel}(i).Nucleus)
                    AssignedNucleusNewChannel(i)=Particles{CurrentChannel}(i).Nucleus;
                else
                    AssignedNucleusNewChannel(i)=nan;
                end
            end

            if ~isempty(find(AssignedNucleusNewChannel==AssignedNucleusPreviousChannel))
                CurrentParticle=find(AssignedNucleusNewChannel==AssignedNucleusPreviousChannel);
            end
        
        %If we don't have a histone channel, go for the same particle
        %number in the new channel or for the last particle
        elseif numParticles<CurrentParticle
            CurrentParticle=numParticles;
        end
    
    elseif cc=='~'      %Switch projection mode
        projectionMode = chooseProjection;
        disp(['projectionMode : ' projectionMode])
    
    elseif cc=='!' %Increase contrast in the Overlay figure 
       if isempty(DisplayRangeSpot)
            DisplayRangeSpot=[min(min(Image)),max(max(Image))/1.5];
        else
            DisplayRangeSpot=[DisplayRangeSpot(1),DisplayRangeSpot(2)/1.5];
       end
    elseif cc=='@'      %Decrease spot channel contrast
        DisplayRangeSpot=[min(min(Image)),max(max(Image))*1.5];
    elseif cc=='0'      %Debugging mode
        keyboard;
    end
        
end


save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo')

%If we only have one channel bring Particles back to the legacy
%format without any cells
if NChannels==1
    Particles=Particles{1};
    Spots=Spots{1};
    SpotFilter=SpotFilter{1};
end


if UseHistoneOverlay
    save([DataFolder,filesep,'Particles.mat'],'Particles','SpotFilter','Threshold1','Threshold2', '-v7.3')
    save([DataFolder,filesep,'Spots.mat'],'Spots', '-v7.3') 
    save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'],'schnitzcells')
else
    save([DataFolder,filesep,'Particles.mat'],'Particles','SpotFilter','Threshold1','Threshold2', '-v7.3')            
    save([DataFolder,filesep,'Spots.mat'],'Spots','-v7.3')
end
close all
disp('Particles saved.')
disp(['(Left off at Particle #', num2str(CurrentParticle), ')'])
%% Extra stuff that is useful in debug mode

%Reset approve status of all approved particles in a certain nc
% nc=13;
% 
% %Determine the start and end frame of the nc
% if nc==14
%     disp('Do this')
% else
%     eval(['ncStart=nc',num2str(nc),';']);
%     eval(['ncEnd=nc',num2str(nc+1),';']);
% end
% 
% for i=1:length(Particles)
%     if Particles(i).Approved==1
%     
%         if (min(Particles(i).Frame(Particles(i).FrameApproved))>=ncStart)&...
%                 (min(Particles(i).Frame(Particles(i).FrameApproved))<ncEnd)
%             Particles(i).Approved=0;
%         end
%     end
% end


%This is if the reference from schnitz to nucleus is screwed up. It
%manifests itself in that cenx and ceny of the schnitz will not coincide
%with the actual ellipse. In this case we delete the schnitz (emptying it)
%and the nuclear reference and use '1' to recreate the schnitz.

% ParticlesToDecouple=[314:316,318,320:325,328,329,331:333,335,337,339,341,342,344:346]
% 
% for ParticleToDecouple=ParticlesToDecouple
% 
%     SchnitzToDecouple=Particles(ParticleToDecouple).Nucleus;
% 
%     %Delete the schnitz by emptying it
%     schnitzcells(SchnitzToDecouple).P=[];   
%     schnitzcells(SchnitzToDecouple).E=[];
%     schnitzcells(SchnitzToDecouple).D=[];
%     schnitzcells(SchnitzToDecouple).frames=[];
%     schnitzcells(SchnitzToDecouple).cenx=[];
%     schnitzcells(SchnitzToDecouple).ceny=[];
%     schnitzcells(SchnitzToDecouple).len=[];
%     schnitzcells(SchnitzToDecouple).cellno=[];
% 
%     %Decouple particle and schnitz
%     Particles(ParticleToDecouple).Nucleus=[];
% end
