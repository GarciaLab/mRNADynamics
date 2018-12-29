function CompileParticles(varargin)
% CompileParticles(varargin)
%
% DESCRIPTION
% This function puts together all the information we have about particles.
%
% ARGUMENTS
% varargin: A cell in which the first element is the prefix string of the data set
%           to analyze. Subsequent elements can be the options below.
%
% OPTIONS
% 'ForceAP': Force AP detection even if it's there already.
%
% 'SkipTraces': Don't output the individual traces.
%
% 'SkipFluctuations': Don't generate the plots of the correlation of signal
%                   and offset.
%
% 'SkipFits': Don't do the fits
%
% 'SkipMovie': Don't do the movie
%
% 'SkipAll': Skip all that can be skipped
%
% 'ApproveAll': Approves all particles. This is useful if we want to do a
%             quick check of, for example, the AP profiles
%
% 'MinParticles', N: Set the threshold for the minimum number of particles (N) per
%               AP bin for compilation. Default is 4.
%
% 'MinTime', M: %Require particles to exist for time M or else discard.
%               Default is 1.
%
% 'ROI', ROI1, ROI2: For Region of Interest (ROI) data. Assume that the ROI is top half of the imaging window.
%           Note that the origin is the left top of the image.
%           ROI1 and ROI2 are the y-position of the ROI rectangle. ROI1 is the
%           lower boundary of ROI and ROI2 is the upper boundary of non-ROI
%           since there is almost always scattering in the middle
% 'intArea': Change the area (in pixels) of integration used in offset calculations
% 'noHist': Force the code to assume there's no nuclear channel.
% 'doSingleFits': Generate single trace fits. Added by EL&AR 
% Author (contact): Hernan Garcia (hggarcia@berkeley.edu)
% Created:
% Last Updated: 6/17/17 (AR)
%
% Documented by: Hernan Garcia (hggarcia@berkeley.edu)

close all;

%% INITIALIZE ALL SAVED VARIABLES
% Please initialize any new variables you have added and want to save!!!!
APFilter  =  {};
APbinArea = [];
APbinID = [];
AllTracesAP  =  {};
AllTracesVector = {};
CompiledParticles = {};
DVFilter = {};
DVbinArea = [];
DVbinID = [];
ElapsedTime = [];
EllipsePos = {};
EllipsesFilteredPos = [];
EllipsesOnAP = {};
EllipsesOnDV = {};
FilteredParticlesPos = [];
MaxAPIndex = [];
MaxCyto = [];
MaxDVIndex = [];
MaxFrame = {};
MeanCyto = [];
MeanOffsetVector = [];
MeanSlopeVectorAP = {};
MeanSlopeVectorDV = {};
MeanVectorAP = {};
MeanVectorAP_ROI = {};
MeanVectorAP_nonROI = {};
MeanVectorAll = {};
MeanVectorAllAP = {};
MeanVectorAllDV = {};
MeanVectorAnterior = {};
MeanVectorDV = {};
MeanVectorDV_ROI = {};
MeanVectorDV_nonROI = {};
MedianCyto = [];
MinAPIndex = [];
MinDVIndex = [];
NEllipsesAP = [];
NEllipsesDV = [];
NOffsetParticles = [];
NParticlesAP = {};
NParticlesAP_ROI = {};
NParticlesAP_nonROI = {};
NParticlesAll = {};
NParticlesDV = {};
NParticlesDV_ROI = {};
NParticlesDV_nonROI = {};
NSlopeAP = {};
NSlopeDV = {};
NewCyclePos = [];
OnRatioAP = {};
OnRatioDV = {};
ParticleCountAP = {};
ParticleCountDV = {};
ParticleCountProbAP = {};
ParticleCountProbDV = {};
SDCyto = [];
SDOffsetVector = [];
SDSlopeVectorAP = {};
SDSlopeVectorDV = {};
SDVectorAP = {};
SDVectorAP_ROI = {};
SDVectorAP_nonROI = {};
SDVectorAll = {};
SDVectorDV = {};
SDVectorDV_ROI = {};
SDVectorDV_nonROI = {};
SEVectorAllAP = {};
SEVectorAllDV = {};
StemLoopEnd = ''; 
TotalEllipsesAP = [];
TotalEllipsesDV = [];
fittedLineEquations = [];
rateOnAP = [];
timeOnOnAP = [];
rateOnAPCell = [];
timeOnOnAPCell = [];

nc9 = [];
nc10 = [];
nc11 = [];
nc12 = [];
nc13 = [];
nc14 = [];
ncFilter = [];
ncFilterID = [];
%%

%Information about about folders
[~,~,DefaultDropboxFolder,~,~]=...
    DetermineLocalFolders;

%Look at the input parameters and use defaults if missing
Prefix='';
ForceAP=0;      %Force AP detection even if it's already there
SkipTraces=0;   %Do not output the individual traces.
SkipFluctuations=0;  %Do not generate the plots of correlations of fluctuations and offset
SkipFits=0;         %Do not generate the fit output (but still does the fit)
SkipMovie=0;        %Do not generate the movie
SkipAll=0;          %Do not do other things 
ApproveAll=0;       %Only use manually approved particles
MinParticles=4;
minTime = 1;
ROI=0; % No ROI
intArea = 109; %default for 220nm x 220nm zoom.
noHist = 0; 
doSingleFits = 0;



% Checking Varargin 
if isempty(varargin)%looks for the folder to analyze
    FolderTemp=uigetdir(DefaultDropboxFolder,'Select folder with data to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
else
    Prefix=varargin{1};
    for i=2:length(varargin)
        if strcmpi(varargin{i},'ForceAP')
            ForceAP=1;
        elseif strcmpi(varargin{i},'SkipTraces')
            SkipTraces=1;
        elseif strcmpi(varargin{i},'SkipFluctuations')
            SkipFluctuations=1;
        elseif strcmpi(varargin{i},'SkipFits')
            SkipFits=1;
        elseif strcmpi(varargin{i},'SkipMovie')
            SkipMovie=1;
        elseif strcmpi(varargin{i},'SkipAll')
            SkipTraces=1;
            SkipFluctuations=1;
            SkipFits=1;
            SkipMovie=1;
            SkipAll=1;
        elseif strcmpi(varargin{i},'doSingleFits')
            doSingleFits=1;
        elseif strcmpi(varargin{i},'ApproveAll')
            ApproveAll=1;
            disp('Approved')
        elseif strcmpi(varargin{i},'noHist')
            noHist = 1;
        elseif strcmp(varargin{i},'MinParticles')
            if ~isnumeric(varargin{i+1})
                error('Wrong input parameters. After ''MinParticles'' you should input the desired minimum number of particles per approved AP bin')
            else
                MinParticles=varargin{i+1};
            end
        elseif strcmpi(varargin{i},'intArea')
            if ~isnumeric(varargin{i+1})
                error('Wrong input parameters. After ''intArea'' you should input the desired number of pixels for intensity integration')
            else
                intArea=varargin{i+1};
            end
        elseif strcmpi(varargin{i},'MinTime')
            if ~isnumeric(varargin{i+1})
                error('Wrong input parameters. After ''MinTime'' you should input the desired minimum number of frames per particle.')
            else
                minTime=varargin{i+1};
            end
        elseif strcmpi(varargin{i},'ROI')
            ROI = 1;
            if ~isnumeric(varargin{i+1})||~isnumeric(varargin{i+2})
                error('Wrong input parameters. After ''ROI'' you should input the y-threshold of ROI ')
            else
                ROI1=varargin{i+1};
                ROI2=varargin{i+2};
            end
        end
    end
    
end

FilePrefix=[Prefix,'_'];

%What type of experiment are we dealing with? Get this out of MovieDatabase
[~,~,DropboxFolder,~, PreProcPath,...
    ~, ~, ~, ~, ~,~] = readMovieDatabase(Prefix);


% refactor in progress, we should replace readMovieDatabase with getExperimentDataFromMovieDatabase
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoopEnd, APResolution,...
    Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)


%Load all the information
load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'Spots.mat'])

if isempty(Particles)
    SkipTraces=1;
    SkipFluctuations=1;
    SkipFits=1;
    SkipMovie=1;
end

%Check that FrameInfo exists
if exist([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'file')
    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
else
    warning('No FrameInfo.mat found. Trying to continue')
    %Adding frame information
    DHis=dir([PreProcPath,filesep,FilePrefix(1:end-1),filesep,'*His*.tif']);
    if ~isempty(DHis)
        FrameInfo(length(DHis)).nc=[];
    end
    
    %Adding information
    
    Dz=dir([PreProcPath,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'*001*.tif']);
    NumberSlices=length(Dz)-1;
    
    for i=1:numFrames
        FrameInfo(i).NumberSlices=NumberSlices;
    end
end

numFrames = length(FrameInfo);



%See how  many frames we have and adjust the index size of the files to
%load accordingly
if numFrames<1E3
    NDigits=3;
elseif numFrames<1E4
    NDigits=4;
else
    error('No more than 10,000 frames currently supported.')
end

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



%Delete the files in folder where we'll write again.
if ~SkipTraces
    delete([DropboxFolder,filesep,Prefix,filesep,'ParticleTraces',filesep,'*.*'])
end
if ~SkipFits
    delete([DropboxFolder,filesep,Prefix,filesep,'Fits',filesep,'*.*'])
end



%See if we had any lineage/nuclear information
if exist([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'file') && ~noHist
    load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
    HistoneChannel=1;
else
    disp('No lineage / nuclear information found. Proceeding without it.');
    HistoneChannel=0;
end


if sum(~cellfun(@isempty,strfind({lower(Channel1{1}),lower(Channel2{1})},'mcherry'))|...
        ~cellfun(@isempty,strfind({lower(Channel1{1}),lower(Channel2{1})},'his')))
    
    % nothing to do, we already have ncs from abobe. refactor this if.
else
    warning('Warning: lack of histone channel may result in strange behavior.');
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



% % Read in which end the stem loops are at, if this information is available
% % (ES 2014-03-20)
% StemLoopEndColumn = find(strcmp(XLSRaw(1, :), 'StemLoopEnd'));
% if ~isempty(StemLoopEndColumn)
%     StemLoopEnd = XLSRaw{XLSEntry, StemLoopEndColumn};
% else
%     StemLoopEnd = '';
% end
%

NewCyclePos=[nc9,nc10,nc11,nc12,nc13,nc14];
NewCyclePos=NewCyclePos(~(NewCyclePos==0));
NewCyclePos=NewCyclePos(~isnan(NewCyclePos));



%Add the APPosition to Particles if they don't exist yet. Do this only if
%we took AP data. Otherwise just add XY.

if strcmpi(ExperimentAxis,'AP')
    if (~isfield(Particles{1},'APpos')) || ForceAP
        if HistoneChannel
            AddParticlePosition(Prefix);
        else
            AddParticlePosition(Prefix,'SkipAlignment')
        end
    else
        disp('Using saved AP information')
    end
elseif strcmpi(ExperimentAxis,'dv') && exist([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'],'file')
    AddParticlePosition(Prefix);
elseif strcmpi(ExperimentAxis,'dv') || strcmpi(ExperimentAxis,'NoAP')
    AddParticlePosition(Prefix,'NoAP');
else
    error('Experiment axis not recognized in MovieDatabase')
end

%Reload Particles.mat
load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
%Create the particle array. This is done so that we can support multiple
%channels. Also figure out the number of channels
if iscell(Particles)
    NChannels=length(Particles);
else
    Particles={Particles};
    NChannels=1;
end


if HistoneChannel
    load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'])
end


%Folders for reports
warning('off', 'MATLAB:MKDIR:DirectoryExists');
if strcmpi(ExperimentAxis,'AP')
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APMovie'])
end
mkdir([DropboxFolder,filesep,Prefix,filesep,'ParticleTraces'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Offset'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Fits'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Probabilities'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Various']);


%% Put together CompiledParticles

CompiledParticles = cell(NChannels,1);
%Approve all particles if the mode has been selected
if ApproveAll
    for ChN=1:NChannels
        %Check that the approved field is present. If not include
        %it. This can occur if CheckParticleTracking is not run
        %first.
        if ~isfield(Particles{ChN},'Approved')
            for i=1:length(Particles{ChN})
                Particles{ChN}(i).Approved=0;
            end
        end
        
        for i=1:length(Particles{ChN})
            %Make sure the particle has an associated nucleus if we are in
            %HistoneChannel mode
            if HistoneChannel
                if ~isempty(Particles{ChN}(i).Nucleus)
                    %If a particle has been explicitly rejected then don't
                    %approve it!
                    if Particles{ChN}(i).Approved~=-1
                        Particles{ChN}(i).Approved=1;
                    end
                end
            else
                Particles{ChN}(i).Approved=1;
            end
        end
    end
end


if HistoneChannel
    for ChN=1:NChannels
        %Get information about which particle is associated to which nucleus in
        %schnitzcells
        for i=1:length(Particles{ChN})
            if ~isempty(Particles{ChN}(i).Nucleus)
                AssignedNuclei{ChN}(i)=Particles{ChN}(i).Nucleus;
            else
                AssignedNuclei{ChN}(i)=nan;
            end
        end
    end
end


%First, figure out the AP position of each of the nuclei.
if strcmpi(ExperimentAxis, 'AP') || strcmpi(ExperimentAxis, 'DV')
    %Load the AP detection information
    load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
    %Angle between the x-axis and the AP-axis
    if exist('coordPZoom', 'var')
        %APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
        %Changed for DV compatibility
        APAngle=atan2((coordPZoom(2)-coordAZoom(2)),(coordPZoom(1)-coordAZoom(1)));
    else
        error('coordPZoom not defined. Was AddParticlePosition.m run?')
    end
    APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
end

if HistoneChannel&&strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')
    %The information in Ellipses is
    %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
    for i=1:length(Ellipses)
        for j=1:size(Ellipses{i},1)
            
            %Angle between the x-axis and the particle using the A position as a
            %zero
            %Angles=atan((Ellipses{i}(j,2)-coordAZoom(2))./(Ellipses{i}(j,1)-coordAZoom(1)));
            %Changed for DV compatibility
            Angles=atan2((Ellipses{i}(j,2)-coordAZoom(2)),(Ellipses{i}(j,1)-coordAZoom(1)));
            
            %Distance between the points and the A point
            Distances=sqrt((coordAZoom(2)-Ellipses{i}(j,2)).^2+(coordAZoom(1)-Ellipses{i}(j,1)).^2);
            APPositions=Distances.*cos(Angles-APAngle);
            EllipsePos{i}(j)=APPositions/APLength;
            
            %Added DV compatibility
            DVPositions=Distances.*sin(Angles-APAngle);
            EllipsePos_DV{i}(j)=DVPositions;
        end
    end
end


%Get the actual time corresponding to each frame
if isfield(FrameInfo,'FileMode')
    if strcmp(FrameInfo(end).FileMode,'TIF')
        for j=1:numFrames
            ElapsedTime(j)=etime(datevec(FrameInfo(j).TimeString),datevec(FrameInfo(1).TimeString));
        end
    elseif strcmp(FrameInfo(end).FileMode,'LSM')||strcmp(FrameInfo(end).FileMode,'LSMExport')||...
            strcmp(FrameInfo(end).FileMode,'LIFExport')||strcmp(FrameInfo(end).FileMode,'LAT')
        for j=1:numFrames
            ElapsedTime(j)=FrameInfo(j).Time-FrameInfo(1).Time;
        end
    else
        error('File mode not supported. Cannot extract time information. Include format in ExportDataForFISH.m')
    end
else
    warning('No FileMode information found. Assuming that this is TIF from the 2-photon.')
    for j=1:numFrames
        ElapsedTime(j)=etime(datevec(FrameInfo(j).TimeString),datevec(FrameInfo(1).TimeString));
    end
end

ElapsedTime=ElapsedTime/60;     %Time is in minutes


%Some parameters
MinAPArea=12500;%700;    %Minimum area in pixels in order to consider an AP bin as valid. AR 3/15/16: This should be recalculated in microns


if strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')
    
    APbinID=0:APResolution:1;
    %Create an image for the different AP bins
    APPosImage=zeros(FrameInfo(1).LinesPerFrame,FrameInfo(1).PixelsPerLine);
    [Rows,Columns]=size(APPosImage);
    
    for i=1:Rows
        for j=1:Columns
            try
                Angle=atan((i-coordAZoom(2))./(j-coordAZoom(1)));
            catch
                Angle=atan2((i-coordAZoom(2)),(j-coordAZoom(1)));
            end
            Distance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
            APPosition=Distance.*cos(Angle-APAngle);
            APPosImage(i,j)=APPosition/APLength;
        end
    end
    
    
    APPosBinImage=zeros(size(APPosImage));
    for i=1:(length(APbinID)-1)
        FilteredMask=(APbinID(i)<=APPosImage)&(APbinID(i+1)>APPosImage);
        APPosBinImage=APPosBinImage+FilteredMask*i;
    end
    
    
    %Calculate the area in pixels corresponding to each AP bin. We will use
    %this to get rid of small AP bins in the image and also to calculate
    %probabilities of nuclei being active.
    APbinArea = zeros(length(APbinID),1);
    %Calculate ther areas of the AP bins
    for i=1:length(APbinID)
        APbinArea(i)=sum(sum(APPosBinImage==i));
    end
    %Get the median of the non-zero areas
    MedianArea=median(APbinArea(APbinArea>0));
    %Only keep the bins with an area of at least 70% of the median
    APbinArea(APbinArea<MedianArea*0.7)=nan;
    
end

%Generate DV bin
if strcmpi(ExperimentAxis,'DV')
       
    DVbinID=linspace(-800,0,51); %JAKE: Would change to DV resolution later
    %Create an image for the different DV bins
    DVPosImage=zeros(FrameInfo(1).LinesPerFrame,FrameInfo(1).PixelsPerLine);
    [Rows,Columns]=size(DVPosImage);

    for i=1:Rows
        for j=1:Columns
            Angle=atan2((i-coordAZoom(2)),(j-coordAZoom(1)));   
            
            Distance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
            DVPosition=Distance.*sin(Angle-APAngle);
            DVPosImage(i,j)=DVPosition;
        end
    end


    DVPosBinImage=zeros(size(DVPosImage));
    for i=1:(length(DVbinID)-1)
        FilteredMask=(DVbinID(i)<=DVPosImage)&(DVbinID(i+1)>DVPosImage);
        DVPosBinImage=DVPosBinImage+FilteredMask*i;
    end


    %Calculate the area in pixels corresponding to each AP bin. We will use
    %this to get rid of small AP bins in the image and also to calculate
    %probabilities of nuclei being active.
    DVbinArea = zeros(length(DVbinID),1);
    %Calculate ther areas of the AP bins
    for i=1:length(DVbinID)
        DVbinArea(i)=sum(sum(DVPosBinImage==i));
    end
    %Get the median of the non-zero areas
    MedianArea=median(DVbinArea(DVbinArea>0));
    %Only keep the bins with an area of at least 70% of the median
    DVbinArea(DVbinArea<MedianArea*0.7)=nan;
    
end


%Now get the particle information for those that were approved
h=waitbar(0,'Compiling traces');
for ChN=1:NChannels
    k=1;
    for i=1:length(Particles{ChN})
        waitbar(i/length(Particles{ChN}),h)
        if (Particles{ChN}(i).Approved==1)
            
            for NCh=1:NChannels
                if ~isfield(Particles{NCh},'FrameApproved')
                    for i=1:length(Particles{NCh})
                        Particles{NCh}(i).FrameApproved=true(size(Particles{NCh}(i).Frame));
                    end
                end
            end
            
            %Which frames were approved manually?
            FrameFilter=Particles{ChN}(i).FrameApproved;
            %What is the first frame that was found, regardless of the column
            %condition?
            FirstFrame=Particles{ChN}(i).Frame(min(find(Particles{ChN}(i).FrameApproved)));
            
            %Check that for the remaining frames we got a good z-profile
            %                 for j=1:length(Particles{ChN}(i).Frame)
            %                     ZProfile=fad(ChN).channels(Particles{ChN}(i).Frame(j)).fits.shadowsDog{Particles{ChN}(i).Index(j)};
            %                     [Dummy,ZMax]=max(ZProfile);
            %                     if (ZMax==1)|(ZMax==length(ZProfile))
            %                         FrameFilter(j)=0;
            %                     end
            %                 end
            
            %Should I only keep traces of a certain length? We also just keep
            %the ones that have a real schnitz associated with them
            AnalyzeThisParticle=1;      %Flag to see if this particle should be analyzed.
            
            if HistoneChannel
                if ~((sum(FrameFilter)>0)&...
                        (~isempty(schnitzcells(Particles{ChN}(i).Nucleus).frames)))
                    AnalyzeThisParticle=0;
                end
            elseif ~(sum(FrameFilter)>0)
                AnalyzeThisParticle=0;
            elseif length(Particles{ChN}(i)) <  minTime
                AnalyzeThisParticle=0;
            end
            
            
            
            %See if this particle is in one of the approved AP bins
            try
                if strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')
                    CurrentAPbin=max(find(APbinID<mean(Particles{ChN}(i).APpos(FrameFilter))));
                    if isnan(APbinArea(CurrentAPbin))
                        AnalyzeThisParticle=0;
                    end
                end
            catch
                error(['You probably need to re-run AddParticlePosition again. If that ',...
                    'doesn''t fix things, talk to HG.'])
            end
            
            
            
            if AnalyzeThisParticle
                
                %Reference to the original Particles index
                CompiledParticles{ChN}(k).OriginalParticle=i;
                
                %Copy the filtered information
                CompiledParticles{ChN}(k).Frame=Particles{ChN}(i).Frame(FrameFilter);
                CompiledParticles{ChN}(k).Index=Particles{ChN}(i).Index(FrameFilter);
                CompiledParticles{ChN}(k).xPos=Particles{ChN}(i).xPos(FrameFilter);
                CompiledParticles{ChN}(k).yPos=Particles{ChN}(i).yPos(FrameFilter);
                %(MT, 2018-02-11) Hacky fix to get lattice to run - FIX LATER
                %CompiledParticles{ChN}(k).DVpos=Particles{ChN}(i).DVpos(FrameFilter);
                CompiledParticles{ChN}(k).FrameApproved = Particles{ChN}(i).FrameApproved;
                
                if strcmpi(ExperimentAxis,'AP')
                    CompiledParticles{ChN}(k).APpos=Particles{ChN}(i).APpos(FrameFilter);
                    
                    %Determine the particles average and median AP position
                    CompiledParticles{ChN}(k).MeanAP=mean(Particles{ChN}(i).APpos(FrameFilter));
                    CompiledParticles{ChN}(k).MedianAP=median(Particles{ChN}(i).APpos(FrameFilter));
                elseif strcmpi(ExperimentAxis,'DV')%&isfield(Particles,'APpos')
                    %AP information:
                    CompiledParticles{ChN}(k).APpos=Particles{ChN}(i).APpos(FrameFilter);
                    CompiledParticles{ChN}(k).MeanAP=mean(Particles{ChN}(i).APpos(FrameFilter));
                    CompiledParticles{ChN}(k).MedianAP=median(Particles{ChN}(i).APpos(FrameFilter));
                    %DV information:
                    CompiledParticles{ChN}(k).DVpos=Particles{ChN}(i).DVpos(FrameFilter);
                    CompiledParticles{ChN}(k).MeanDV=mean(Particles{ChN}(i).DVpos(FrameFilter));
                    CompiledParticles{ChN}(k).MedianDV=median(Particles{ChN}(i).DVpos(FrameFilter));
                    
                end
                
                %If we have the histone channel we will actually replace the AP
                %position by the position of the nucleus where the particle was
                %found. If there is no nucleus (like when a particle survives
                %past the nuclear division) we will still use the actual particle
                %position.
                if HistoneChannel&strcmpi(ExperimentAxis,'AP')
                    %Save the original particle position
                    CompiledParticles{ChN}(k).APposParticle=CompiledParticles{ChN}(k).APpos;
                    
                    FramesToCheck=schnitzcells(Particles{ChN}(i).Nucleus).frames(...
                        ismember(schnitzcells(Particles{ChN}(i).Nucleus).frames,Particles{ChN}(i).Frame(FrameFilter)));
                    EllipsesToCheck=schnitzcells(Particles{ChN}(i).Nucleus).cellno(...
                        ismember(schnitzcells(Particles{ChN}(i).Nucleus).frames,Particles{ChN}(i).Frame(FrameFilter)));
                    
                    for j=1:length(FramesToCheck)
                        IndexToChange=find(CompiledParticles{ChN}(k).Frame==FramesToCheck(j));
                        CompiledParticles{ChN}(k).APPos(IndexToChange)=EllipsePos{FramesToCheck(j)}(EllipsesToCheck(j));
                    end
                end
                
                %First frame it was detected at
                CompiledParticles{ChN}(k).FirstFrame=FirstFrame;
                CompiledParticles{ChN}(k).Approved=Particles{ChN}(i).Approved;
                
                %Copy the fit results if they are there
                if isfield(Particles,'Fit')
                    CompiledParticles{ChN}(k).Fit=Particles{ChN}(i).Fit;
                end
                
                
                %Extract information from Spots about fluorescence and background
                [Frame,AmpIntegral, AmpIntegral3, AmpIntegral5, AmpGaussian,...
                    Off, ErrorIntegral,ErrorGauss,optFit1, FitType, ErrorIntegral3, ErrorIntegral5,backGround3]...
                    = GetParticleTrace(k,CompiledParticles{ChN},Spots{ChN});
                CompiledParticles{ChN}(k).Fluo= AmpIntegral;
                CompiledParticles{ChN}(k).Fluo3= AmpIntegral3;
                CompiledParticles{ChN}(k).Fluo5= AmpIntegral5;
                CompiledParticles{ChN}(k).FluoGauss= AmpGaussian;
                CompiledParticles{ChN}(k).Off=Off;
                CompiledParticles{ChN}(k).FluoError=ErrorIntegral(1); % SEANCHANGED
                CompiledParticles{ChN}(k).optFit1=optFit1;
                CompiledParticles{ChN}(k).FitType=FitType;
                
                
                %Determine the nc where this particle was born
                try
                    CompiledParticles{ChN}(k).nc=FrameInfo(CompiledParticles{ChN}(k).Frame(1)).nc;
                catch
                end
                
                if HistoneChannel
                    CompiledParticles{ChN}(k).Nucleus=Particles{ChN}(i).Nucleus;
                    %We have two fields with the same information:
                    %"Nucleus" and "schnitz". In future versions we'll get rid of
                    %"Nucleus"
                    CompiledParticles{ChN}(k).schnitz=Particles{ChN}(i).Nucleus;
                    
                    %Save lineage information in terms of particles
                    if ~isempty(schnitzcells(Particles{ChN}(i).Nucleus).P)
                        if isempty(find(AssignedNuclei{ChN}==schnitzcells(Particles{ChN}(i).Nucleus).P))
                            CompiledParticles{ChN}(k).PParticle=0;
                        else
                            CompiledParticles{ChN}(k).PParticle=find(AssignedNuclei{ChN}==schnitzcells(Particles{ChN}(i).Nucleus).P);
                        end
                    else
                        CompiledParticles{ChN}(k).PParticle=[];
                    end
                    
                    if ~isempty(schnitzcells(Particles{ChN}(i).Nucleus).D)
                        if isempty(find(AssignedNuclei{ChN}==schnitzcells(Particles{ChN}(i).Nucleus).D))
                            CompiledParticles{ChN}(k).DParticle=0;
                        else
                            CompiledParticles{ChN}(k).DParticle=find(AssignedNuclei{ChN}==schnitzcells(Particles{ChN}(i).Nucleus).D);
                        end
                    else
                        CompiledParticles{ChN}(k).DParticle=[];
                    end
                    
                    if ~isempty(schnitzcells(Particles{ChN}(i).Nucleus).E)
                        if isempty(find(AssignedNuclei{ChN}==schnitzcells(Particles{ChN}(i).Nucleus).E))
                            CompiledParticles{ChN}(k).EParticle=0;
                        else
                            CompiledParticles{ChN}(k).EParticle=find(AssignedNuclei{ChN}==schnitzcells(Particles{ChN}(i).Nucleus).E);
                        end
                    else
                        CompiledParticles{ChN}(k).EParticle=[];
                    end
                    
                    %Save information about the nucleus birth and death
                    CompiledParticles{ChN}(k).NucStart=schnitzcells(Particles{ChN}(i).Nucleus).frames(1);
                    CompiledParticles{ChN}(k).NucEnd=schnitzcells(Particles{ChN}(i).Nucleus).frames(end);
                    
                    
                end
                
                
                
                %Plot and save this trace together with its offset value
                
                if ~SkipTraces
                    if ~isnan(nc9)|~isnan(nc10)|~isnan(nc11)|~isnan(nc12)|~isnan(nc13)|~isnan(nc14)
                        %ncFilterID just tells you the identity of the different
                        %filters stored in the cell ncFilter
                        ncFilterID=[];
                        if nc9~=0
                            ncFilterID=9;
                        end
                        if nc10~=0
                            ncFilterID=[ncFilterID,10];
                        end
                        if nc11~=0
                            ncFilterID=[ncFilterID,11];
                        end
                        if nc12~=0
                            ncFilterID=[ncFilterID,12];
                        end
                        if nc13~=0
                            ncFilterID=[ncFilterID,13];
                        end
                        if nc14~=0
                            ncFilterID=[ncFilterID,14];
                        end
                        %Add the first nc
                        ncFilterID=[min(ncFilterID)-1,ncFilterID];
                        
                        
                        %Create the filter
                        % for ChN=1:NChannels
                        
                        if isempty(CompiledParticles)==1
                            error(['No compiled particles found in channel ',num2str(ChN),'. Did you mean to run the code with ApproveAll?'])
                        end
                        
                        ncFilter=false(length(CompiledParticles{ChN})...
                            ,length(ncFilterID)); %AR 6/16/17: I think multi-channel data might require this to be a cell? Something for the future.
                        
                        for i=1:length(CompiledParticles{ChN})
                            %Sometimes CompiledParticles{1}(i).nc is empty. This is because of some
                            %problem with FrameInfo! In that case we'll pull the information out of
                            %the XLS file.
                            if ~isfield(CompiledParticles{ChN}(i), 'nc')
                                CompiledParticles{ChN}(i).nc = [];
                            end
                            if ~isempty(CompiledParticles{ChN}(i).nc)
                                ncFilter(i,find(CompiledParticles{ChN}(i).nc==ncFilterID))=true;
                            else
                                ncsFound=find(CompiledParticles{ChN}(i).Frame(1)>=[nc9,nc10,nc11,nc12,nc13,nc14]);
                                if ncsFound(end)==1
                                    CompiledParticles{ChN}(i).nc=9;
                                    ncFilter(i,ncFilterID==9)=true;
                                elseif ncsFound(end)==2
                                    CompiledParticles{ChN}(i).nc=10;
                                    ncFilter(i,ncFilterID==10)=true;
                                elseif ncsFound(end)==3
                                    CompiledParticles{ChN}(i).nc=11;
                                    ncFilter(i,ncFilterID==11)=true;
                                elseif ncsFound(end)==4
                                    CompiledParticles{ChN}(i).nc=12;
                                    ncFilter(i,ncFilterID==12)=true;
                                elseif ncsFound(end)==5
                                    CompiledParticles{ChN}(i).nc=13;
                                    ncFilter(i,ncFilterID==13)=true;
                                elseif ncsFound(end)==6
                                    CompiledParticles{ChN}(i).nc=14;
                                    ncFilter(i,ncFilterID==14)=true;
                                end
                            end
                        end
                        %end
                    end
                    
                    figure(2)
                    left_color = [213,108,85]/255;
                    right_color = [0, 0, 0]/255;
                    set(gcf,'defaultAxesColorOrder',[left_color; right_color]);
                    %Size of the snippet for each frame
                    SnippetSize=31; %AR 3/15/16: Why is this 31?
                    %Width of the screen
                    ScreenWidth=get( 0, 'ScreenSize' );
                    ScreenWidth=ScreenWidth(3);
                    
                    %Figure out the arrangement of snippets
                    NFrames=length(CompiledParticles{ChN}(k).Frame);
                    
                    NRows=ceil(NFrames/ScreenWidth*SnippetSize)*2;
                    NCols=max([2,ceil(NFrames/NRows)]);
                    
                    %Actual total number of rows
                    TotalRows=14;
                    
                    subplot(TotalRows,NCols,[1:((TotalRows-NRows)*NCols)])
                    
                    
                    
                    %Top left plot
                    FilterMatrix=zeros((TotalRows-NRows),NCols);
                    FilterMatrix(:,1:ceil(NCols/2))=1;
                    subplot(TotalRows,NCols,find(FilterMatrix'))
                    
                    %                     yyaxis left
                    errorbar(ElapsedTime(CompiledParticles{ChN}(k).Frame),...
                        CompiledParticles{ChN}(k).Fluo,ones(size(CompiledParticles{ChN}(k).Fluo))*...
                        CompiledParticles{ChN}(k).FluoError,...
                        '.-r');
                    ylabel('fluorescence (au)')
                    hold on
                    
                    %                     yyaxis right
                    plot(ElapsedTime(CompiledParticles{ChN}(k).Frame),...
                        CompiledParticles{ChN}(k).Off*intArea,'.-g');
                    if ~isempty(CompiledParticles{ChN}(k).optFit1)
                        
                        if strcmp(CompiledParticles{ChN}(k).FitType,'spline')
                            SplineValues=ppval(CompiledParticles{ChN}(k).optFit1,double(CompiledParticles{ChN}(k).Frame));
                        elseif strcmp(CompiledParticles{ChN}(k).FitType,'mean')
                            SplineValues=ones(size(CompiledParticles{ChN}(k).Frame))*CompiledParticles{ChN}(k).optFit1;
                        elseif strcmp(CompiledParticles{ChN}(k).FitType,'line')
                            SplineValues=polyval(CompiledParticles{ChN}(k).optFit1,CompiledParticles{ChN}(k).Frame);
                        end
                        
                        %                         yyaxis right
                        plot(ElapsedTime(CompiledParticles{ChN}(k).Frame),SplineValues*intArea,'-b')
                        try
                            title(['particle ',num2str(k),'(',num2str(i),'), nc',num2str(CompiledParticles{ChN}(k).nc),', Ch: ',num2str(ChN)])
                        catch
                        end
                    else
                        title(['particle ',num2str(k),'(',num2str(i),'), nc',num2str(CompiledParticles{1}(k).nc),', Ch: ',num2str(ChN),...
                            ' - WARNING: No offset fit'])
                    end
                    hold off
                    legend({'Particle','Offset','Offset fit'},'Location','Best')
                    xlabel('time (min)')
                    axis square
                    set(gca, 'Position', get(gca, 'OuterPosition') - ...
                        get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
                    drawnow
                    
                    
                    
                    
                    %Top right plot
                    if HistoneChannel
                        subplot(TotalRows,NCols,find(~FilterMatrix'))
                        
                        if length(CompiledParticles{ChN}(k).Frame)>1
                            colormap(jet(128));
                            cmap=colormap;
                            
                            ColorTime=[];
                            for j=1:length(CompiledParticles{ChN}(k).Frame)
                                ColorTime(j,:)= cmap(round((j-1)/...
                                    (length(CompiledParticles{ChN}(k).Frame)-1)*127+1),:);
                            end
                        else
                            ColorTime=[];
                            ColorTime(1,:)=[1,0,0];
                        end
                        
                        
                        
                        hold on
                        for j=1:length(CompiledParticles{ChN}(k).Frame)
                            PosSchnitz=find((schnitzcells(CompiledParticles{ChN}(k).Nucleus).frames)==...
                                CompiledParticles{ChN}(k).Frame(j));
                            PosEllipse=schnitzcells(CompiledParticles{ChN}(k).Nucleus).cellno(PosSchnitz);
                            CurrEllipse=Ellipses{CompiledParticles{ChN}(k).Frame(j)}(PosEllipse,:);
                            
                            if ~isempty(CurrEllipse)
                                EllipseHandle=ellipse(CurrEllipse(3),...
                                    CurrEllipse(4),...
                                    CurrEllipse(5),...
                                    0,0,[],[],gca);
                                set(EllipseHandle,'color',ColorTime(j,:))
                                plot(CompiledParticles{ChN}(k).xPos(j)-CurrEllipse(1),...
                                    CompiledParticles{ChN}(k).yPos(j)-CurrEllipse(2),'o','color',...
                                    ColorTime(j,:))
                            else
                                PosSchnitz=length(schnitzcells(CompiledParticles{ChN}(k).Nucleus).frames);
                                PosEllipse=schnitzcells(CompiledParticles{ChN}(k).Nucleus).cellno(PosSchnitz);
                                CurrEllipse=...
                                    Ellipses{schnitzcells(CompiledParticles{ChN}(k).Nucleus).frames(PosSchnitz)}(PosEllipse,:);
                                
                                
                                
                                plot(CompiledParticles{ChN}(k).xPos(j)-CurrEllipse(1),...
                                    CompiledParticles{ChN}(k).yPos(j)-CurrEllipse(2),'o','color',...
                                    ColorTime(j,:))
                            end
                            
                            
                            
                        end
                        hold off
                        box on
                        xlabel('x position (pixels)')
                        ylabel('y position (pixels)')
                        axis square
                        set(gca, 'Position', get(gca, 'OuterPosition') - ...
                            get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
                        set(gca, 'YDir', 'reverse')
                        BarHandle = colorbar;
                        set(BarHandle,'YTick',[])
                        
                        %%%%%AR 6/16/17: This still doesn't work in 2017. We need to find a
                        %%%%%replacement.
                        %                         try
                        %                             if ~isempty(cbfreeze(BarHandle))
                        %                                 BarHandle=cbfreeze(BarHandle);
                        %                             else
                        %                                 warning('Issue with cbfreeze.m. Skipping it. The color bar will not reflect time appropriately.')
                        %                             end
                        %                             ylabel(BarHandle,'Time')
                        %                         catch
                        %                             warning('Issue with cbfreeze.m. Skipping it. The color bar will not reflect time appropriately. This is an issue of Matlab 2014b.')
                        %                         end
                        %%%%%%%%%%%%%%%%
                        if strcmpi(ExperimentAxis,'AP')
                            title(['Mean AP: ',num2str(CompiledParticles{ChN}(k).MeanAP)])
                        end
                        drawnow
                    end
                    
                    
                    %Snippets
                    for j=1:NFrames
                        subplot(TotalRows,NCols,(TotalRows-NRows)*NCols+j)
                        spotFrame = CompiledParticles{ChN}(k).Frame(j);
                        [x,y,z]=SpotsXYZ(Spots{ChN}(spotFrame));
                        if ~isempty(x)
                            xTrace=x(CompiledParticles{ChN}(k).Index(j));
                            yTrace=y(CompiledParticles{ChN}(k).Index(j));
                            zTrace=z(CompiledParticles{ChN}(k).Index(j));
                            Image=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                                FilePrefix,iIndex(CompiledParticles{ChN}(k).Frame(j),NDigits),'_z',iIndex(zTrace,2),...
                                '_ch',iIndex(ChN,2),'.tif']);
                            [ImRows,ImCols]=size(Image);
                            ImageSnippet=zeros(SnippetSize,SnippetSize);
                            yRange=round(yTrace)+[-(SnippetSize-1)/2:(SnippetSize-1)/2];
                            yFilter=(yRange>0)&(yRange<=ImRows);
                            xRange=round(xTrace)+[-(SnippetSize-1)/2:(SnippetSize-1)/2];
                            xFilter=(xRange>0)&(xRange<=ImCols);
                            ImageSnippet(yFilter,xFilter)=Image(yRange(yFilter),...
                                xRange(xFilter));
                            imshow(ImageSnippet,[],'Border','Tight','InitialMagnification',200)
                            set(gca, 'Position', get(gca, 'OuterPosition') - ...
                                get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
                            if HistoneChannel
                                %Plot the corresponding nucleus
                                CurrentSchnitz=schnitzcells(CompiledParticles{ChN}(k).Nucleus);
                                if sum((CurrentSchnitz.frames)==CompiledParticles{ChN}(k).Frame(j))==1
                                    hold on
                                    EllipseNumber=CurrentSchnitz.cellno(...
                                        (CurrentSchnitz.frames)==CompiledParticles{ChN}(k).Frame(j));
                                    CurrEllipse=Ellipses{CompiledParticles{ChN}(k).Frame(j)}(EllipseNumber,:);
                                    EllipseHandle=ellipse(CurrEllipse(3),...
                                        CurrEllipse(4),...
                                        CurrEllipse(5),...
                                        CurrEllipse(1)-xTrace+(SnippetSize-1)/2,...
                                        CurrEllipse(2)-yTrace+(SnippetSize-1)/2,[],[],gca);
                                    %set(EllipseHandle,'color',ColorTime(j,:))
                                    set(EllipseHandle,'color','g')
                                    hold off
                                end
                            end
                        end
                    end
                    set(gcf,'Position',[1,41,1280,684])
                    
                    drawnow
                    if isfield(CompiledParticles{ChN}(k), 'nc')
                        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'ParticleTraces',filesep,iIndex(k,NDigits),...
                            '(',num2str(i),')-nc',...
                            num2str(CompiledParticles{ChN}(k).nc),'_ch',iIndex(ChN,2),'.tif'])
                    end
                    close(2)
                end
                
                
                k=k+1;
                
            end
        end
    end
end
close(h)


%% ROI option
% This option is separating the CompiledParticles defined above into
% CompiledParticles_ROI and COmpiledParticles_nonROI
% written by YJK on 10/24/2017
for ChN=1:NChannels
    if ROI
        % separate the CompileParticles into CompiledParticles_ROI and
        % CompiledParticles_nonROI using the Threshold (y position)
        t=1;
        s=1;
        
        % Use the ROI1 and ROI2 to split the Particles
        for ParticleIndex=1:length(CompiledParticles{ChN})
            if nanmean(CompiledParticles{ChN}(ParticleIndex).yPos) < ROI1
                CompiledParticles_ROI{ChN}(t)=CompiledParticles{ChN}(ParticleIndex);
                t=t+1;
            elseif nanmean(CompiledParticles{ChN}(ParticleIndex).yPos) > ROI2
                CompiledParticles_nonROI{ChN}(s)=CompiledParticles{ChN}(ParticleIndex);
                s=s+1;
            end
        end
    end
end
%% Create filters

% APFilter needs to be changed for ROI option, so that we can have two
% filters for each CompiledParticles (ROI and nonROI) (YJK on 10/24/2017)

%nc filters:

[ncFilterID, ncFilter, APFilter, APFilter_ROI, APFilter_nonROI, ...
    DVFilter, DVFilter_ROI, DVFilter_nonROI] = createFilters(nc9, nc10, nc11, nc12, ...
    nc13, nc14, NChannels, CompiledParticles, ExperimentAxis, ROI, APFilter,...
    APFilter_ROI, APFilter_nonROI, APbinID, DVbinID, DVFilter, DVFilter_ROI, ...
    DVFilter_nonROI, CompiledParticles_ROI, CompiledParticles_nonROI);

%% Binning and averaging data

% Here I need to think about how to bin ROI and non-ROI particles.
% One way is having a function that splits the CompiledParticles into CompiledParticles_ROI
% and CompiledParticles_nonROI according to its y-position (Threshold). Then, pass them through
% the AverageTraces (YJK on 10/22/2017)

[AllTracesVector, AllTracesAP, AllTracesDV, MeanVectorAP_ROI, ...
    SDVectorAP_ROI, NParticlesAP_ROI, MeanVectorAP_nonROI, SDVectorAP_nonROI, ...
    NParticlesAP_nonROI, MeanVectorAP, SDVectorAP, NParticlesAP, MeanVectorDV_ROI, ...
    SDVectorDV_ROI, NParticlesDV_ROI, MeanVectorDV_nonROI, SDVectorDV_nonROI, ...
    NParticlesDV_nonROI, MeanVectorDV, SDVectorDV, NParticlesDV, ...
    MeanVectorAnterior, MeanVectorAll, SDVectorAll, NParticlesAll] =...
    binAndAverage(NChannels, CompiledParticles, FrameInfo, ExperimentAxis, ...
    APFilter, ROI, CompiledParticles_ROI, APFilter_ROI, APFilter_nonROI, NewCyclePos, MaxFrame, ...
    MeanVectorAP_ROI, SDVectorAP_ROI, NParticlesAP_ROI, MeanVectorAP_nonROI, ...
    SDVectorAP_nonROI, NParticlesAP_nonROI, MeanVectorAP, SDVectorAP, ...
    NParticlesAP, MeanVectorDV_ROI, SDVectorDV_ROI, NParticlesDV_ROI, ...
    MeanVectorDV_nonROI, SDVectorDV_nonROI, NParticlesDV_nonROI, ...
    MeanVectorDV, SDVectorDV, NParticlesDV, MeanVectorAnterior, DVFilter_ROI, ...
    DVFilter_nonROI);


%% Instantaneous rate of change

[CompiledParticles, MeanSlopeVectorAP, SDSlopeVectorAP, NSlopeAP]...
    = instantRateOfChange(CompiledParticles, ElapsedTime, ExperimentAxis, APFilter);

%% Integrating each particle

CompiledParticles = integrateParticles(NChannels, ElapsedTime, CompiledParticles);


%% Information about the cytoplasm
%If the nuclear masks are present then use them. Otherwise just calculate
%the median of the images as a function of time

%HG on 8/6/16: Why was this commented out? Did I do this?


% if NChannels==1
%
%     if HistoneChannel&strcmpi(ExperimentAxis,'AP')
%         [MeanCyto,SDCyto,MedianCyto,MaxCyto]=GetCytoMCP(Prefix);
%     else
%         MeanCyto=[];
%         SDCyto=[];
%         MaxCyto=[];
%
%         h=waitbar(0,'Calculating the median cyto intentisy');
%         for i=1:numFrames
%             waitbar(i/numFrames,h)
%             for j=1:FrameInfo(1).NumberSlices
%                 Image(:,:,j)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(i,3),'_z',iIndex(j,2),'.tif']);
%             end
%             ImageMax=max(Image,[],3);
%             MedianCyto(i)=median(double(ImageMax(:)));
%         end
%         close(h)
%     end
% else
MeanCyto=[];
SDCyto=[];
MaxCyto=[];
MedianCyto=[];
% end




%% Offset and fluctuations

[MeanOffsetVector, SDOffsetVector, NOffsetParticles] = offsetAndFlux(NChannels, ...
    SkipFluctuations, ncFilter, ElapsedTime, CompiledParticles, DropboxFolder, ...
    Prefix, ExperimentAxis, intArea, MeanVectorAll, SDVectorAll);

%% Rate of mRNA production

CompiledParticles = mRNAProdRate(NChannels, CompiledParticles, ...
    ncFilter, ElapsedTime, SkipFits, DropboxFolder, Prefix);


%% First frames

firstFramesThing(NChannels, HistoneChannel, SkipAll, nc14, ...
    CompiledParticles, DropboxFolder, Prefix, ElapsedTime);

%% AP position of particle vs nucleus

if HistoneChannel&&strcmpi(ExperimentAxis,'AP')
    CompiledParticles = APPosParticleVsNucleus(NChannels, ...
        CompiledParticles, schnitzcells, EllipsePos, DropboxFolder, Prefix);
end


%% Fitting shapes` to single traces (includes time on and initial rate of loading)
% This section of code is will fit line segments piece-wise to the
% single traces. fittedLineEquations correspond to the stored fitted lines
% of the particles, where the indexing is as follows:
% fittedLineEquations(particleNumber) with fields: Coefficients,
% ErrorEstimation, and frameIndex, which are described below. This currently
% does not support more than one channel. Please contact Emma to work on
% implementing it for two channels.
if doSingleFits
    [CompiledParticles, fittedLineEquations] = fitShapesToTraces(Prefix, ...
        Particles, schnitzcells, FrameInfo, ElapsedTime, CompiledParticles);
end

%% Probability of being on

%I'm going to measure the probability of a nucleus having detectable
%expressiona as a function of time and AP. In order to do this I'll use
%Particles that have both the Approved flag set to 1 and 2. However, I'll
%also check that the nuclei are not too close to the edges.

%NOTE: I need a way to go back and check the nuclei that weren't on. Maybe
%I should move this to Check particles


%Create an image that is partitioned according to the AP bins. We will use
%this to calculate the area per AP bin.

if HistoneChannel&&strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')
    [NEllipsesAP, MeanVectorAllAP, SEVectorAllAP, EllipsesFilteredPos, ...
        FilteredParticlesPos, OnRatioAP, ParticleCountAP, ParticleCountProbAP, ...
        EllipsesOnAP, rateOnAP, rateOnAPCell, timeOnOnAP, timeOnOnAPCell, TotalEllipsesAP]...
        = APProbOn(NChannels, Particles, schnitzcells, ...
        CompiledParticles, Ellipses, APbinID, FrameInfo, ElapsedTime, DropboxFolder, ...
        Prefix);
end

% DV version
if HistoneChannel&&strcmpi(ExperimentAxis,'DV') %JAKE: Need to change this later
    [ParticleNuclei, ParticleFrames, CompiledParticleNuclei, CompiledParticleFrames, ...
        NEllipsesDV, MeanVectorAllDV, SEVectorAllDV, OnRatioDV, ParticleCountDV, ...
        ParticleCountProbDV, TotalEllipsesDV, EllipsesOnDV, EllipsesFilteredPos, ...
        FilteredParticlesPos] = DVProbOn(NChannels, ...
        Particles, schnitzcells, CompiledParticles, Ellipses, FrameInfo, ...
        DropboxFolder, Prefix, ElapsedTime, DVbinID);
end



%% Calculation of particle speed
calcParticleSpeeds(NChannels, Particles, ...
    Spots, ElapsedTime, schnitzcells, Ellipses) % this doesn't appear to do anything...

%% Movie of AP profile

%I want to make a movie of the average fluorescence as a function of AP as
%a function of time. In order to make life easier I'll just export to a
%folder. I can then load everything in ImageJ.

if ~SkipMovie&&strcmpi(ExperimentAxis,'AP')
    APProfileMovie(MeanVectorAP, NParticlesAP, MinParticles, numFrames, ...
        APbinID, SDVectorAP, FrameInfo, ElapsedTime, DropboxFolder, Prefix)
end

%% Checking correlations with position and expression
%
% 1+1;
%
% i=100
%
% for j=1:length(CompiledParticles{1}(i).Frame)
%
%     CurrentFrame=CompiledParticles{1}(i).Frame(j);
%     EllipsePos=find((schnitzcells(CompiledParticles{1}(i).Nucleus).frames)==CurrentFrame);
%     CurrentEllipse=Ellipses{CompiledParticles{1}(i).Frame};
%     CurrentEllipse=CurrentEllipse(EllipsePos,:);
%
%     Position(j)=(CompiledParticles{1}(i).xPos(j)-CurrentEllipse(1))^2+...
%         (CompiledParticles{1}(i).yPos(j)-CurrentEllipse(2))^2;
% end
%
% figure(9)
% plot(  Position ,CompiledParticles{1}(i).Fluo,'.k')
%
%


%% Input-output

%Compile the nuclear fluorescence information if we have the appropriate
%experiment type and axis
%Simon: This script is not required for Dl-Venus experiments in DV so I
%added the ExperimentAxis condition.
% ROI option added (YJK on 10/27/2017) :
% When the data is acquired in ROI mode, I added the ROI option,
% two thresholds to be incorporated into CompileNuclearProtein
if strcmpi(ExperimentType,'inputoutput') && (strcmpi(ExperimentAxis,'AP')||strcmpi(ExperimentAxis,'DV'))
    if ~ROI
        CompileNuclearProtein(Prefix)
    else
        CompileNuclearProtein(Prefix,'ROI',ROI1,ROI2)
    end
end



%% Save everything

%Now save all the information
save([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],...
    'APFilter', 'APbinArea', 'APbinID', 'AllTracesAP',...
    'AllTracesVector', 'CompiledParticles', 'DVFilter', 'DVbinArea',...
    'DVbinID', 'ElapsedTime', 'EllipsePos', 'EllipsesFilteredPos',...
    'EllipsesOnAP', 'EllipsesOnDV', 'FilteredParticlesPos', 'MaxAPIndex',...
    'MaxCyto', 'MaxDVIndex', 'MaxFrame', 'MeanCyto',...
    'MeanOffsetVector', 'MeanSlopeVectorAP', 'MeanSlopeVectorDV', 'MeanVectorAP',...
    'MeanVectorAP_ROI', 'MeanVectorAP_nonROI', 'MeanVectorAll', 'MeanVectorAllAP',...
    'MeanVectorAllDV', 'MeanVectorAnterior', 'MeanVectorDV', 'MeanVectorDV_ROI',...
    'MeanVectorDV_nonROI', 'MedianCyto', 'MinAPIndex', 'MinDVIndex',...
    'NEllipsesAP', 'NEllipsesDV', 'NOffsetParticles', 'NParticlesAP',...
    'NParticlesAP_ROI', 'NParticlesAP_nonROI', 'NParticlesAll', 'NParticlesDV',...
    'NParticlesDV_ROI', 'NParticlesDV_nonROI', 'NSlopeAP', 'NSlopeDV',...
    'NewCyclePos', 'OnRatioAP', 'OnRatioDV', 'ParticleCountAP',...
    'ParticleCountDV', 'ParticleCountProbAP', 'ParticleCountProbDV', 'Prefix',...
    'SDCyto', 'SDOffsetVector', 'SDSlopeVectorAP', 'SDSlopeVectorDV',...
    'SDVectorAP', 'SDVectorAP_ROI', 'SDVectorAP_nonROI', 'SDVectorAll',...
    'SDVectorDV', 'SDVectorDV_ROI', 'SDVectorDV_nonROI', 'SEVectorAllAP',...
    'SEVectorAllDV', 'StemLoopEnd', 'TotalEllipsesAP', 'TotalEllipsesDV',...
    'fittedLineEquations', 'rateOnAP', 'timeOnOnAP','rateOnAPCell', 'timeOnOnAPCell','nc10', 'nc11', 'nc12',...
    'nc13', 'nc14', 'nc9', 'ncFilter',...
    'ncFilterID','-v7.3');
 
end