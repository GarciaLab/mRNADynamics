function CompiledNuclei = CompileNuclearProteinFlexFluo(Prefix, varargin)
% date created: 9/25/20
% author: Gabriella Martini
% date last modified: 10/7/20


%This code gives me the input
%varargin Variable length input argument list.
%allows any number of arguments to a function.  The variable
%varargin is a cell array containing the optional arguments to the
%function.  varargin must be declared as the last input argument
%and collects all the inputs from that point onwards.

% ROI (Region of Interest, illumination) option is added on 10/27/2017 by YJK
% ex ) CompileNuclearProtein(Prefix,'ROI','direction',[ROI1 ROI1])
% 1st varargin : 'ROI' or 'nonROI'. This is for file saving name.
% 2nd varargin : the direction, either 'x', or 'y'.
% 3rd varargin : the boundaries of ROI.
% i.e.) [x1_ROI x2_ROI] (or [y1_ROI y2_ROI] )

% Note. The main change is in the 'APbin filter' and 'binning and
% averaging' part of the code.



%initialize variables
savedVariables = {};
MeanCytoAPProfile = [];
SDCytoAPProfile = [];
SECytoAPProfile = [];
IntegrationArea = NaN;
MeanVectorAll_ROI = [];
SDVectorAll_ROI = [];
NParticlesAll_ROI = [];
MeanVectorAll_nonROI = [];
SDVectorAll_nonROI = [];
NParticlesAll_nonROI = [];
MaxFrame = [];
MeanVectorAP_ROI = [];
SDVectorAP_ROI = [];
NParticlesAP_ROI = [];
MeanVectorAP_nonROI = [];
SDVectorAP_nonROI = [];
NParticlesAP_nonROI = [];
CompiledNuclei_ROI = [];
CompiledNuclei_nonROI = [];
MeanVectorAP_ROI = [];
MinDVIndex = [];
MaxDVIndex = [];
SDVectorDV_ROI = [];
NParticlesDV_ROI = [];
MeanVectorDV_nonROI = [];
SDVectorDV_nonROI = [];
NParticlesDV_nonROI = [];
DVbinID = [];
DVFilter = [];
MeanVectorDV = [];
SDVectorDV = [];
NParticlesDV = [];
MeanVectorAP = [];
SDVectorAP = [];
NParticlesAP = [];
CompiledNuclei = [];


methods = {'max', 'gaussiansmooth'};
if ~isempty(varargin)
    x = 1;
    while x <= length(varargin)
        switch varargin{x}

            case{'radii'}
                radii = varargin{x+1};
                x = x+1;
            case{'methods'}
                methods = varargin{x+1};
                x = x+1;
            otherwise
                error('invalid arguments passed')
            end
            x = x+1;
    end
end


%This function will add fluorescence information to each schnitz.

close all



[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath,...
    configValues, movieDatabasePath, movieDatabaseFolder, movieDatabase]=...
    DetermineLocalFolders(Prefix);


liveExperiment = LiveExperiment(Prefix);

%Load all the information
load([DropboxFolder,filesep,Prefix,'\Ellipses.mat'], 'Ellipses')
load([DropboxFolder,filesep,Prefix,'\FrameInfo.mat'], 'FrameInfo')
DetermineNucleiEndFrames(Prefix);
schnitzcells = CalculateNuclearMovement(Prefix);

FrameInfo = getFrameInfo(liveExperiment);
Ellipses = getEllipses(liveExperiment);
nc_info = [liveExperiment.nc9, liveExperiment.nc10, liveExperiment.nc11,...
    liveExperiment.nc12, liveExperiment.nc13, liveExperiment.nc14, length(FrameInfo)];
PixelSize = liveExperiment.pixelSize_um;
nucleusDiameters = zeros(1, 6);
for nc=9:14
    nucleusDiameters(nc-8) = getDefaultParameters(FrameInfo,[['d', num2str(nc)]])/PixelSize; % in pixels
end
xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
zDim = liveExperiment.zDim;
ChN = liveExperiment.inputChannels(1);
if length(liveExperiment.inputChannels) > 1
    error('multiple protein challenges currently unsupported.')
end
% start - Added by GM 9/25/20
if exist([DropboxFolder,filesep,Prefix,'\GridDivision.mat'], 'file')
    load([DropboxFolder,filesep,Prefix,'\GridDivision.mat'], 'GridDivision')
end
% end - Added by GM 9/25/20

load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
    

numFrames = length(FrameInfo);

%See if FrameInfo has information about the number of input channels. This
%is not fully implemented yet. If no information is found, then assume we
%have only one input channel.
if isfield(FrameInfo,'NChInput')
    NChannels=FrameInfo(1).NChInput;
else
    NChannels=1;
end

% refactor in progress, we should replace readMovieDatabase with getExperimentDataFromMovieDatabase
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
    Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, movieDatabase);

%Pre-calculating ExperimentAxis boolean for faster use in later if statements
ExperimentAxisIsNoAP = strcmpi(ExperimentAxis, 'NoAP');


fullEmbryo = exist([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'], 'file');





%% 



%Add the APPosition to Particles if they don't exist yet
if (~isfield(schnitzcells,'APpos'))&&(strcmpi(ExperimentAxis,'AP')||strcmpi(ExperimentAxis,'DV')) && fullEmbryo
    %First, run this to get the alignment between the zoom-in and zoom-out
    %images:
    %     AddParticlePosition(Prefix)
    %Now, add the nuclear position
    AddNuclearPosition(Prefix);
    schnitzcells = getSchnitzcells(liveExperiment);
end

if ~isfield(schnitzcells, 'HistoneFluo')
    schnitzcells = integrateHistoneFluo(Prefix, schnitzcells, FrameInfo);
end

%% Figure out which DVpositions are represented in the image
%Angle between the x-axis and the AP-axis
APAngle=atan2((coordPZoom(2)-coordAZoom(2)),(coordPZoom(1)-coordAZoom(1)));
APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);

ZoomImageSize = [FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine];
Distances = [sqrt((coordAZoom(2)-1).^2+(coordAZoom(1)-1).^2),...
    sqrt((coordAZoom(2)-ZoomImageSize(1)).^2+(coordAZoom(1)-ZoomImageSize(2)).^2),...
    sqrt((coordAZoom(2)-1).^2+(coordAZoom(1)-ZoomImageSize(2)).^2),...
    sqrt((coordAZoom(2)-ZoomImageSize(1)).^2+(coordAZoom(1)-1).^2)];
Angles=[atan2((1-coordAZoom(2)),...
         (-coordAZoom(1))),...
         atan2((ZoomImageSize(1)-coordAZoom(2)),...
         (ZoomImageSize(2)-coordAZoom(1))),...
         atan2((1-coordAZoom(2)),...
         (ZoomImageSize(2)-coordAZoom(1))),...
         atan2((ZoomImageSize(1)-coordAZoom(2)),...
         (1-coordAZoom(1)))];
DVlimits = Distances.*sin(Angles-APAngle)/APLength;

DVBinMinIdx = round(min(min(DVlimits))/APResolution)-1;
DVBinMin = DVBinMinIdx*APResolution;
DVBinMaxIdx = round(max(max(DVlimits))/APResolution)+1;
DVBinMax = DVBinMaxIdx*APResolution;
DVbinID = DVBinMin:APResolution:DVBinMax;

% for i=1:length(schnitzcells)
%     %Angle between the x-axis and the particle using the A position as a
%     %zero
%     Angles=atan2((double(schnitzcells(i).ceny)-coordAZoom(2)),...
%         (double(schnitzcells(i).cenx)-coordAZoom(1)));
% 
%     %Distance between the points and the A point
%     Distances=sqrt((coordAZoom(2)-double(schnitzcells(i).ceny)).^2+(coordAZoom(1)-double(schnitzcells(i).cenx)).^2);
%     APPositions=Distances.*cos(Angles-APAngle);
%     schnitzcells(i).APpos=APPositions/APLength;
% 
%     %Determine the distance perpendicular to the AP axis. This is a
%     %proxy for a DV axis.
% 
%     schnitzcells(i).DVpos=Distances.*sin(Angles-APAngle);
% 
% end

%% 




% Put together CompiledNuclei

%Get the actual time corresponding to each frame

for j=1:numFrames
    ElapsedTime(j)=FrameInfo(j).Time-FrameInfo(1).Time;%Finds the elapsed time by subtracting each time point by the initial time point
end


ElapsedTime=ElapsedTime/60;     %Time is in minutes

%If there is no Approved field then create it
if ~isfield(schnitzcells,'Approved')
    for i=1:length(schnitzcells)
        schnitzcells(i).Approved=true;
    end
end

%If there is no FrameApproved field then create it
if ~isfield(schnitzcells,'FrameApproved')
    for i=1:length(schnitzcells)
        schnitzcells(i).FrameApproved=true(size(schnitzcells(i).frames));
    end
end

schnitzcells = filterSchnitz(schnitzcells, [ liveExperiment.yDim, liveExperiment.xDim]);
[schnitzcells, Ellipses] = correctSchnitzCellErrors(schnitzcells, Ellipses);




%% 

%Now get the nuclear information for those that were approved
NZSlices=size(schnitzcells(1).Fluo,2);
%CompiledNuclei(length(schnitzcells))=struct;
FluoLabels = {};
if ~exist('radii', 'var')
    radii = [];
    SCfieldnames = fieldnames(schnitzcells);
    for i=1:length(SCfieldnames)
        fn = SCfieldnames{i};
        if regexp(fn, '^Fluo[_0-9]*um')
            FluoLabels{end+1} = fn;
            r = strrep(fn, 'Fluo', '');
            r = strrep(r, 'um', '');
            r = strrep(r, '_', '.');
            radii(end + 1) = str2double(r);
        end
    end
else
    
    for i=1:length(radii)
        r = radii(i);
        FluoLabels{i} = strrep(['Fluo', num2str(r), 'um'], '.', '_');
    end
end

MinCycle = find(nc_info > 0, 1) + 8;

h=waitbar(0,'Compiling nuclear traces');
k=1;
for i=1:length(schnitzcells)
    
    
    waitbar(i/length(schnitzcells),h)
    
    if ((schnitzcells(i).Approved) & (schnitzcells(i).cycle >= MinCycle))
        %Which frames were approved manually?
        FrameFilter=schnitzcells(i).FrameApproved;
        
        %Check that for the remaining frames we got a good z-profile
        for j=1:length(schnitzcells(i).frames)
            [MaxValue,MaxPos]=max(schnitzcells(i).Fluo(j,:));
            if NZSlices<3
                if (MaxPos==2)||(MaxPos==NZSclices-1)
                    FrameFilter(j)=false;
                end
            end
        end
        
        if sum(FrameFilter)
           
            %Copy the filtered information
%             if isfield(schnitzcells, 'P')
%                 CompiledNuclei(k).P=schnitzcells(i).P;
%                 CompiledNuclei(k).E=schnitzcells(i).E;
%                 CompiledNuclei(k).D=schnitzcells(i).D;
%             end
            CompiledNuclei(k).Frames=uint16(schnitzcells(i).frames(FrameFilter));
            CompiledNuclei(k).xPos=single(schnitzcells(i).cenx(FrameFilter));
            CompiledNuclei(k).yPos=single(schnitzcells(i).ceny(FrameFilter));
            CompiledNuclei(k).Radius=single(schnitzcells(i).len(FrameFilter));
            CompiledNuclei(k).cellno=uint16(schnitzcells(i).cellno(FrameFilter));
            CompiledNuclei(k).nc=schnitzcells(i).cycle;
           
            %Save the information about the original schnitz
            CompiledNuclei(k).schnitz=uint16(i);
            
            if ~ExperimentAxisIsNoAP && fullEmbryo
                CompiledNuclei(k).MeanDV=single(mean(schnitzcells(i).DVpos(FrameFilter)));
                CompiledNuclei(k).MedianDV=single(median(schnitzcells(i).DVpos(FrameFilter)));
                CompiledNuclei(k).MeanAP=single(mean(schnitzcells(i).APpos(FrameFilter)));
                CompiledNuclei(k).MedianAP=single(median(schnitzcells(i).APpos(FrameFilter)));
            end
            
            CompiledNuclei(k).APpos = double(schnitzcells(i).APpos(FrameFilter));
            CompiledNuclei(k).DVpos = double(schnitzcells(i).DVpos(FrameFilter));
            if isfield(schnitzcells, 'anaphaseFrame')
                CompiledNuclei(k).timeSinceAnaphase = double(schnitzcells(i).timeSinceAnaphase(FrameFilter));
                CompiledNuclei(k).anaphaseFrame=schnitzcells(i).anaphaseFrame;
                if ~isempty(schnitzcells(i).inferredAnaphaseFrame) 
                    if ~isnan(schnitzcells(i).inferredAnaphaseFrame)
                        CompiledNuclei(k).inferredAnaphaseFrame=logical(schnitzcells(i).inferredAnaphaseFrame);
                    else
                        CompiledNuclei(k).inferredAnaphaseFrame = false;
                    end
                else
                    CompiledNuclei(k).inferredAnaphaseFrame = false;
                end
            else
                CompiledNuclei(k).timeSinceAnaphase = NaN;
                CompiledNuclei(k).anaphaseFrame=NaN;
                CompiledNuclei(k).inferredAnaphaseFrame=false;
            end
            if isfield(schnitzcells, 'containsFirstFrameOfCycle')
                CompiledNuclei(k).containsLastFrameOfCycle = schnitzcells(i).containsFirstFrameOfCycle;
            end
            if isfield(schnitzcells, 'containsLastFrameOfCycle')
                CompiledNuclei(k).containsLastFrameOfCycle = schnitzcells(i).containsLastFrameOfCycle;
            end
            
            CompiledNuclei(k).VelocityInfo = schnitzcells(i).VelocityInfo;
            for j=1:length(radii)
                FluoLabel = FluoLabels{j};
                for l=1:length(methods)
                    method = methods{l};
                    var1 = [FluoLabel, '_', method, '_TimeTrace'];
                    var2 = [FluoLabel, '_', method, '_FrameApproved'];
                    var3 = [FluoLabel, '_', method, '_Z'];
                    [CompiledNuclei(k).(var1), CompiledNuclei(k).(var2),...
                        CompiledNuclei(k).(var3)] = CalculateNuclearFluorescence(schnitzcells(i).(FluoLabel)(FrameFilter,:), method);
                    CompiledNuclei(k).(var1) = CompiledNuclei(k).(var1)/(pi*radii(j)^2);
                end
                
            end
            CompiledNuclei(k).HistoneFluo = schnitzcells(i).HistoneFluo;
            CompiledNuclei(k).Flag1 = 0;
            CompiledNuclei(k).Flag2 = 0;
            CompiledNuclei(k).Flag3 = 0;
            CompiledNuclei(k).Flag4 = 0;
            CompiledNuclei(k).Flag5 = 0;
            CompiledNuclei(k).Flag6 = 0;
            CompiledNuclei(k).Flag7 = 0;
            if isfield(schnitzcells, 'Flag')
                if schnitzcells(i).Flag > 0
                    FlagVar = ['Flag', num2str(schnitzcells(i).Flag)];
                    CompiledNuclei(k).(FlagVar) = 1;
                end
            end
            
            
           
          
            
            %If there was only one time point and multiple channels,
            %squeeze can lead to a weird shape of the matrix
            
            k=k+1;
        end
    end
end
close(h)

%% ROI option
% This option is sorting the CompiledNuclei that are within the ROI region
% that is defined by direction, [ROI1 ROI2]
% written by YJK on 6/6/2019



%% Create AP and nc filters
% 
% %nc filters:
% 
% %ncFilterID just tells you the identity of the different
% %filters stored in the cell ncFilter
% ncFilterID=[];
% if nc9~=0
%     ncFilterID=9;
% end
% if nc10~=0
%     ncFilterID=[ncFilterID,10];
% end
% if nc11~=0
%     ncFilterID=[ncFilterID,11];
% end
% if nc12~=0
%     ncFilterID=[ncFilterID,12];
% end
% if nc13~=0
%     ncFilterID=[ncFilterID,13];
% end
% if nc14~=0
%     ncFilterID=[ncFilterID,14];
% end
% %Add the first nc
% ncFilterID=[min(ncFilterID)-1,ncFilterID];
% 
% 
% %Create the filter
% ncFilter=false(length(CompiledNuclei),length(ncFilterID));
% for i=1:length(CompiledNuclei)
%     if ~isempty(CompiledNuclei(i).Frames)
%         
%         if ~isempty(CompiledNuclei(i).nc)
%             ncFilter(i,find(CompiledNuclei(i).nc==ncFilterID))=true;
%         else
%             ncsFound=find(CompiledNuclei(i).Frames(1)>=[nc9,nc10,nc11,nc12,nc13,nc14]);
%             if ncsFound(end)==1
%                 CompiledNuclei(i).nc=9;
%                 ncFilter(i,ncFilterID==9)=true;
%             elseif ncsFound(end)==2
%                 CompiledNuclei(i).nc=10;
%                 ncFilter(i,ncFilterID==10)=true;
%             elseif ncsFound(end)==3
%                 CompiledNuclei(i).nc=11;
%                 ncFilter(i,ncFilterID==11)=true;
%             elseif ncsFound(end)==4
%                 CompiledNuclei(i).nc=12;
%                 ncFilter(i,ncFilterID==12)=true;
%             elseif ncsFound(end)==5
%                 CompiledNuclei(i).nc=13;
%                 ncFilter(i,ncFilterID==13)=true;
%             elseif ncsFound(end)==6
%                 CompiledNuclei(i).nc=14;
%                 ncFilter(i,ncFilterID==14)=true;
%             end
%             
%         end
%     end
% end
% 
% 
% 
% if strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV') && fullEmbryo
%     %AP filters:
%     
%     %Divide the AP axis into boxes of a certain AP size. We'll see which
%     %particle falls where.
%     
%     APbinID=0:APResolution:1;
%     
%     APFilter=false(length(CompiledNuclei),length(APbinID));
%     
%     for i=1:length(CompiledNuclei)
%         APFilter(i,max(find(APbinID<=CompiledNuclei(i).MeanAP)))=true;
%     end
% else
%     APbinID = [];
%     APFilter = [];
% end
% 
% 
% 
% %% Information about the cytoplasmic fluroescence
% %If the nuclear masks are present then use them. Otherwise just calculate
% %the median of the images as a function of time
% MeanCyto=[];
% SDCyto=[];
% MaxCyto=[];
% MedianCyto = [];
% if strcmp(ExperimentAxis,'AP')|| strcmp(ExperimentAxis,'DV')
%     [MeanCyto,SDCyto,MedianCyto,MaxCyto,...
%         MeanCytoAPProfile,SDCytoAPProfile,SECytoAPProfile]=GetCytoMCP(Prefix);
% if false
% else
%     MeanCyto=[];
%     SDCyto=[];
%     MaxCyto=[];
%     if (~isempty(strfind(Channel1{1}, 'Bcd')))
%         nameSuffix=['_ch',iIndex(1,2)];
%     else
%         nameSuffix=['_ch',iIndex(2,2)]; %assuming input channel is channel two. obviously this is going to be wrong in general.
%     end
%
%     MedianCyto = zeros(numFrames);
%     h=waitbar(0,'Calculating the median cyto intentisy');
%     for frame=1:numFrames
%         waitbar(i/numFrames,h)
%         for z=1:FrameInfo(1).NumberSlices
%             Image(:,:,z)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(i,3),'_z',iIndex(j,2),nameSuffix,'.tif']);
%         end
%         ImageMax=max(Image,[],3);
%         MedianCyto(frame)=median(double(ImageMax(:)));
%     end
%     close(h)
% end


%% Binning and averaging data
% 
% if strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV') && fullEmbryo
%     
%     %Get the data for the individual particles in a matrix that has the frame
%     %number and the particle number as dimensions. Also, get a vector that
%     %reports the mean position.
%     [AllTracesVector,AllTracesAP, AllTracesDV]=AllTracesNuclei(FrameInfo,CompiledNuclei);
%     
%     
%     
%     %Mean plot for different AP positions
%     
%     %Figure out the AP range to use
%     MinAPIndex=1;%min(find(sum(APFilter)));
%     MaxAPIndex=size(APFilter,2);%max(find(sum(APFilter)));
%     
%     %Get the corresponding mean information
%     k=1;
%     for ap=MinAPIndex:MaxAPIndex
%         [MeanVectorAPTemp,SDVectorAPTemp,NParticlesAPTemp]=AverageTracesNuclei(FrameInfo,...
%             CompiledNuclei(APFilter(:,ap)),NChannels);
%         MeanVectorAPCell{k}=MeanVectorAPTemp';
%         SDVectorAPCell{k}=SDVectorAPTemp';
%         NParticlesAPCell{k}=NParticlesAPTemp';
%         k=k+1;
%     end
%     
%     %Turn the information into useful structures
%     if NChannels>1
%         for ch=1:NChannels
%             for ap=MinAPIndex:MaxAPIndex
%                 MeanVectorAPCell2{ch,ap}=MeanVectorAPCell{ap}{ch};
%                 SDVectorAPCell2{ch,ap}=SDVectorAPCell{ap}{ch};
%                 NParticlesAPCell2{ch,ap}=NParticlesAPCell{ap}{ch};
%             end
%         end
%         
%         for ch=1:NChannels
%             MeanVectorAP{ch}=cell2mat({MeanVectorAPCell2{ch,:}}')';
%             SDVectorAP{ch}=cell2mat({SDVectorAPCell2{ch,:}}')';
%             NParticlesAP{ch}=cell2mat({NParticlesAPCell2{ch,:}}')';
%         end
%     else
%         for ap=MinAPIndex:MaxAPIndex
%             MeanVectorAPCell2{ap}=double(MeanVectorAPCell{ap});
%             SDVectorAPCell2{ap}=double(SDVectorAPCell{ap});
%             NParticlesAPCell2{ap}=double(NParticlesAPCell{ap});
%         end
%         
%         MeanVectorAP=cell2mat(MeanVectorAPCell2);
%         SDVectorAP=cell2mat(SDVectorAPCell2);
%         NParticlesAP=cell2mat(NParticlesAPCell2);
%     end
%     
%     
% elseif strcmpi(ExperimentAxis,'NoAP')
%     %Get the data for the individual particles in a matrix that has the frame
%     %number and the particle number as dimensions. Also, get a vector that
%     %reports the mean position.
%     [AllTracesVector,AllTracesAP, AllTracesDV]=AllTracesNuclei(FrameInfo,CompiledNuclei,'NoAP');
% end
% 
% 
% 
% 
% %Calculate the mean for all of them
% [MeanVectorAll,SDVectorAll,NParticlesAll]=AverageTracesNuclei(FrameInfo,CompiledNuclei);
% 

% %Now find the different maxima in each nc
%
% MaxFrame=[];
% for i=1:length(NewCyclePos)
%     if i==1
%         [~,MaxIndex]=max(MeanVectorAll(1:NewCyclePos(1)));
%         MaxFrame=[MaxFrame,MaxIndex];
%     elseif i<=length(NewCyclePos)
%         [~,MaxIndex]=max(MeanVectorAll(NewCyclePos(i-1):NewCyclePos(i)));
%         MaxFrame=[MaxFrame,NewCyclePos(i-1)+MaxIndex-1];
%     end
% end
%
% [~,MaxIndex]=max(MeanVectorAll(NewCyclePos(i):end));
% MaxFrame=[MaxFrame,NewCyclePos(i)+MaxIndex-1];

%% add some things to schnitzcells for convenience
% 
% 
% expandedAnaphaseFrames = [zeros(1,8),liveExperiment.anaphaseFrames'];
% 
% for s = 1:length(schnitzcells)
%     midFrame = ceil(length(schnitzcells(s).frames)/2);
%     dif = double(schnitzcells(s).frames(midFrame)) - expandedAnaphaseFrames;
%     cycle = find(dif>0, 1, 'last' );
%     schnitzcells(s).cycle = uint8(cycle);
% end
% 
% schnitzcells = addRelativeTimeToSchnitzcells(schnitzcells, FrameInfo, expandedAnaphaseFrames);
% 


% 
% %% Save everything
% 
% if ~fullEmbryo
%     
%     savedVariables = [savedVariables,...
%         'CompiledNuclei','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
%         'nc12','nc13','nc14','ncFilterID','ncFilter',...
%         'MeanVectorAll','SDVectorAll','NParticlesAll',...
%         'MaxFrame',...
%         'MeanCyto','SDCyto','MedianCyto','MaxCyto',...
%         'IntegrationArea'];
% else
%     savedVariables = [savedVariables,...
%         'CompiledNuclei','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
%         'nc12','nc13','nc14','ncFilterID','ncFilter','APbinID','APFilter',...
%         'MeanVectorAP','SDVectorAP','NParticlesAP',...
%         'MeanVectorAll','SDVectorAll','NParticlesAll',...
%         'MaxFrame',...
%         'AllTracesVector','AllTracesAP',...
%         'MeanCyto','SDCyto','MedianCyto','MaxCyto',...
%         'MeanCytoAPProfile','SDCytoAPProfile','SECytoAPProfile',...
%         'IntegrationArea'...n
%         'DVbinID','DVFilter','MeanVectorDV','SDVectorDV','NParticlesDV',...
%         'MinDVIndex','MaxDVIndex', 'AllTracesDV'];
%     
%     
% end

try
    save([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'],...
        'CompiledNuclei','-v6');
catch
    save([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'],...
        'CompiledNuclei','-v7.3', '-nocompression');
end
save2([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], Ellipses); 
if whos(var2str(schnitzcells)).bytes < 2E9
    save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells', '-v6')
else
    save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells', '-v7.3', '-nocompression')
    
end
