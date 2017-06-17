function CompileParticles(varargin)
% CompileParticles(varargin)
%
% DESCRIPTION
% This function puts together all the information we have about particles.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
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
%               AP bin for compilation
%
% 'MinTime', M: %Require particles to exist for time M or else discard 
%
%
% Author (contact): Hernan Garcia (hggarcia@berkeley.edu)
% Created: 
% Last Updated: 6/14/17 (AR)
%
% Commented by: Hernan Garcia (hggarcia@berkeley.edu)

close all

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
ApproveAll=0;       %Only use manually approved particles
MinParticles=4;     
minTime = 1;        

if isempty(varargin)%looks for the folder to analyze
    FolderTemp=uigetdir(DefaultDropboxFolder,'Select folder with data to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
else
    Prefix=varargin{1};
    for i=2:length(varargin)
        if strcmp(varargin{i},'ForceAP')
            ForceAP=1;
        elseif strcmp(varargin{i},'SkipTraces')
            SkipTraces=1;
        elseif strcmp(varargin{i},'SkipFluctuations')
            SkipFluctuations=1;
        elseif strcmp(varargin{i},'SkipFits')    
            SkipFits=1;
        elseif strcmp(varargin{i},'SkipMovie')    
            SkipMovie=1;
        elseif strcmp(varargin{i},'SkipAll')        
            SkipTraces=1;
            SkipFluctuations=1;
            SkipFits=1;
            SkipMovie=1;
        elseif strcmp(varargin{i},'ApproveAll')    
            ApproveAll=1;
        elseif strcmp(varargin{i},'MinParticles')
            if ~isnumeric(varargin{i+1})
                error('Wrong input parameters. After ''MinParticles'' you should input the desired minimum number of particles per approved AP bin')
            else
                MinParticles=varargin{i+1};
            end
        elseif strcmp(varargin{i},'MinTime')
            if ~isnumeric(varargin{i+1})
                error('Wrong input parameters. After ''MinTime'' you should input the desired minimum number of frames per particle.')
            else
                minTime=varargin{i+1};
            end
        end
    end

end
      
FilePrefix=[Prefix,'_'];

%What type of experiment are we dealing with? Get this out of
%MovieDatabase.xlsx
[~,~,DropboxFolder,~, PreProcPath,...
    ~, ~, ~, ~, ~,~] = readMovieDatabase(Prefix);

%Note that some of this information is redundant given what we get out of
%readMovieDatabase above. We'll have to integrate this better.
[~,XLSTxt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
ExperimentTypeColumn=find(strcmp(XLSRaw(1,:),'ExperimentType'));
ExperimentAxisColumn=find(strcmp(XLSRaw(1,:),'ExperimentAxis'));
APResolutionColumn = find(strcmp(XLSRaw(1,:),'APResolution'));

DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
        if isempty(PrefixRow)
            error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
        end
    end
        
if isempty(PrefixRow)
    error('Entry not found in MovieDatabase.xlsx')
end

ExperimentType=XLSRaw{PrefixRow,ExperimentTypeColumn};
ExperimentAxis=XLSRaw{PrefixRow,ExperimentAxisColumn};
APResolution = XLSRaw{PrefixRow,APResolutionColumn};



%Load all the information
load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'Spots.mat'])
%Check that FrameInfo exists
if exist([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
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
    
    for i=1:length(FrameInfo)
        FrameInfo(i).NumberSlices=NumberSlices;
    end
    
end



%See how  many frames we have and adjust the index size of the files to
%load accordingly
if length(FrameInfo)<1E3
    NDigits=3;
elseif length(FrameInfo)<1E4
    NDigits=4;
else
    error('No more than 10,000 frames supported. Change this in the code')
end





%Create the particle array. This is done so that we can support multiple
%channels. Also figure out the number of channels
if iscell(Particles)
    NChannels=length(Particles);
else
    Particles={Particles};
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
if exist([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
    HistoneChannel=1;
else
    disp('No lineage / nuclear information found. Proceeding without it.');
    HistoneChannel=0;
end




%Load the information about the nc from the XLS file
[~,Txt]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
XLSHeaders=Txt(1,:);
Txt=Txt(2:end,:);

%Find the different columns.

%Convert the prefix into the string used in the XLS file
Dashes=findstr(FilePrefix(1:end-1),'-');



%Determine division times
%Load the information about the nc from the XLS file

[~,Txt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
XLSHeaders=Txt(1,:);
Txt=Txt(2:end,:);

ExperimentTypeColumn=find(strcmp(XLSRaw(1,:),'ExperimentType'));
ExperimentAxisColumn=find(strcmp(XLSRaw(1,:),'ExperimentAxis'));

DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
if isempty(PrefixRow)
    PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
    end
end

ExperimentType=XLSRaw{PrefixRow,ExperimentTypeColumn};
ExperimentAxis=XLSRaw{PrefixRow,ExperimentAxisColumn};

%Find the different columns.
DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
nc9Column=find(strcmp(XLSRaw(1,:),'nc9'));
nc10Column=find(strcmp(XLSRaw(1,:),'nc10'));
nc11Column=find(strcmp(XLSRaw(1,:),'nc11'));
nc12Column=find(strcmp(XLSRaw(1,:),'nc12'));
nc13Column=find(strcmp(XLSRaw(1,:),'nc13'));
nc14Column=find(strcmp(XLSRaw(1,:),'nc14'));
CFColumn=find(strcmp(XLSRaw(1,:),'CF'));
Channel1Column=find(strcmp(XLSRaw(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSRaw(1,:),'Channel2'));

%Get the information
Channel1=XLSTxt(PrefixRow,Channel1Column);
Channel2=XLSTxt(PrefixRow,Channel2Column);

%Convert the prefix into the string used in the XLS file
Dashes=findstr(Prefix,'-');

%Find the corresponding entry in the XLS file
if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))&(isempty(findstr(Prefix,'Bcd1x')))
    XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
        [Date,'\BcdGFP-HisRFP']));
else
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
    if isempty(XLSEntry)
        XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
            [Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
        if isempty(XLSEntry)
            error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
        end
    end
    
end


if sum(~cellfun(@isempty,strfind({lower(Channel1{1}),lower(Channel2{1})},'mcherry'))|...
    ~cellfun(@isempty,strfind({lower(Channel1{1}),lower(Channel2{1})},'his')))
    nc9=cell2mat(XLSRaw(XLSEntry,nc9Column));
    nc10=cell2mat(XLSRaw(XLSEntry,nc10Column));
    nc11=cell2mat(XLSRaw(XLSEntry,nc11Column));
    nc12=cell2mat(XLSRaw(XLSEntry,nc12Column));
    nc13=cell2mat(XLSRaw(XLSEntry,nc13Column));
    nc14=cell2mat(XLSRaw(XLSEntry,nc14Column));
    %This is in case the last column for CF is all nan and is not part of
    %the Num matrix
    if ~isempty(CFColumn)    
        CF=cell2mat(XLSRaw(XLSEntry,CFColumn));
    else
        CF=nan;
    end
else
    warning('Warning: lack of histone channel may result in strange behavior.');
    
    nc9=cell2mat(XLSRaw(XLSEntry,nc9Column));
    nc10=cell2mat(XLSRaw(XLSEntry,nc10Column));
    nc11=cell2mat(XLSRaw(XLSEntry,nc11Column));
    nc12=cell2mat(XLSRaw(XLSEntry,nc12Column));
    nc13=cell2mat(XLSRaw(XLSEntry,nc13Column));
    nc14=cell2mat(XLSRaw(XLSEntry,nc14Column));
    %This is in case the last column for CF is all nan and is not part of
    %the Num matrix
    if ~isempty(CFColumn)    
        CF=cell2mat(XLSRaw(XLSEntry,CFColumn));
    else
        CF=nan;
    end
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



% Read in which end the stem loops are at, if this information is available
% (ES 2014-03-20)
StemLoopEndColumn = find(strcmp(XLSRaw(1, :), 'StemLoopEnd'));
if ~isempty(StemLoopEndColumn)
    StemLoopEnd = XLSRaw{XLSEntry, StemLoopEndColumn};
else
    StemLoopEnd = '';
end


NewCyclePos=[nc9,nc10,nc11,nc12,nc13,nc14];
NewCyclePos=NewCyclePos(~(NewCyclePos==0));
NewCyclePos=NewCyclePos(~isnan(NewCyclePos));



%Add the APPosition to Particles if they don't exist yet. Do this only if
%we took AP data. Otherwise just add XY.

if strcmp(ExperimentAxis,'AP')
    if (~isfield(Particles{1},'APpos')) || ForceAP
        if HistoneChannel
            AddParticlePosition(Prefix);
        else
            AddParticlePosition(Prefix,'SkipAlignment')
        end

    else
        display('Using saved AP information')
    end
elseif strcmp(lower(ExperimentAxis),'dv')& exist([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
    AddParticlePosition(Prefix);
elseif strcmpi(ExperimentAxis,'dv')
    AddParticlePosition(Prefix,'NoAP');
elseif strcmp(ExperimentAxis,'NoAP')
    AddParticlePosition(Prefix,'NoAP');
else
    error('Experiment axis not recognized in MovieDatabase.XLSX')
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
if strcmp(ExperimentAxis,'AP')
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APMovie'])
end
mkdir([DropboxFolder,filesep,Prefix,filesep,'ParticleTraces'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Offset'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Fits'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Probabilities'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Various']);


%% Put together CompiledParticles

CompiledParticles = {};
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
if strcmp(ExperimentAxis, 'AP')
    %Load the AP detection information
    load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
    %Angle between the x-axis and the AP-axis
    if exist('coordPZoom', 'var')
        APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
    else
        error('coordPZoom not defined. Was AddParticlePosition.m run?')
    end
    APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
end

if HistoneChannel&&strcmp(ExperimentAxis,'AP')
    %The information in Ellipses is
    %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
    for i=1:length(Ellipses)
        for j=1:size(Ellipses{i},1)

            %Angle between the x-axis and the particle using the A position as a
            %zero
            Angles=atan((Ellipses{i}(j,2)-coordAZoom(2))./(Ellipses{i}(j,1)-coordAZoom(1)));

            %Distance between the points and the A point
            Distances=sqrt((coordAZoom(2)-Ellipses{i}(j,2)).^2+(coordAZoom(1)-Ellipses{i}(j,1)).^2);
            APPositions=Distances.*cos(Angles-APAngle);
            EllipsePos{i}(j)=APPositions/APLength;
        end
    end
end


%Get the actual time corresponding to each frame
if isfield(FrameInfo,'FileMode')
    if strcmp(FrameInfo(end).FileMode,'TIF')
        for j=1:length(FrameInfo)
            ElapsedTime(j)=etime(datevec(FrameInfo(j).TimeString),datevec(FrameInfo(1).TimeString));
        end
    elseif strcmp(FrameInfo(end).FileMode,'LSM')||strcmp(FrameInfo(end).FileMode,'LSMExport')||...
            strcmp(FrameInfo(end).FileMode,'LIFExport')||strcmp(FrameInfo(end).FileMode,'LAT')
        for j=1:length(FrameInfo)
            ElapsedTime(j)=FrameInfo(j).Time-FrameInfo(1).Time;
        end
    else
        error('File mode not supported. Cannot extract time information. Include format in ExportDataForFISH.m')
    end
else
    warning('No FileMode information found. Assuming that this is TIF from the 2-photon.')
    for j=1:length(FrameInfo)
        ElapsedTime(j)=etime(datevec(FrameInfo(j).TimeString),datevec(FrameInfo(1).TimeString));
    end
end
    
ElapsedTime=ElapsedTime/60;     %Time is in minutes
    

%Some parameters
IntArea=500;%190        %Area of integration. AR 3/15/16: This should be recalculated in microns
MinAPArea=12500;%700;    %Minimum area in pixels in order to consider an AP bin as valid. AR 3/15/16: This should be recalculated in microns


if strcmp(ExperimentAxis,'AP')
       
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
                        Particles{NCh}(i).FrameApproved=logical(ones(size(Particles{NCh}(i).Frame)));
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
            if strcmp(ExperimentAxis,'AP')
                CurrentAPbin=max(find(APbinID<mean(Particles{ChN}(i).APpos(FrameFilter))));
                if isnan(APbinArea(CurrentAPbin))
                    AnalyzeThisParticle=0;
                end
            end



            if AnalyzeThisParticle

                %Reference to the original Particles index
                CompiledParticles{ChN}(k).OriginalParticle=i;            

                %Copy the filtered information
                CompiledParticles{ChN}(k).Frame=Particles{ChN}(i).Frame(FrameFilter);
                CompiledParticles{ChN}(k).Index=Particles{ChN}(i).Index(FrameFilter);
                CompiledParticles{ChN}(k).xPos=Particles{ChN}(i).xPos(FrameFilter);
                CompiledParticles{ChN}(k).yPos=Particles{ChN}(i).yPos(FrameFilter);
                CompiledParticles{ChN}(k).FrameApproved = Particles{ChN}(i).FrameApproved;

                if strcmp(ExperimentAxis,'AP')
                    CompiledParticles{ChN}(k).APpos=Particles{ChN}(i).APpos(FrameFilter);

                    %Determine the particles average and median AP position
                    CompiledParticles{ChN}(k).MeanAP=mean(Particles{ChN}(i).APpos(FrameFilter));
                    CompiledParticles{ChN}(k).MedianAP=median(Particles{ChN}(i).APpos(FrameFilter));
                elseif strcmp(ExperimentAxis,'DV')&isfield(Particles,'APpos')
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
                if HistoneChannel&strcmp(ExperimentAxis,'AP')
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
                [Frame,AmpIntegral,AmpGaussian,Off,...
                 ErrorIntegral,ErrorGauss,optFit1, FitType, noIntensityFlag]...
                 = GetParticleTrace(k,CompiledParticles{ChN},Spots);
                CompiledParticles{ChN}(k).Fluo= AmpIntegral;
                CompiledParticles{ChN}(k).FluoIntegral = AmpIntegral;
                CompiledParticles{ChN}(k).Off=Off;
                CompiledParticles{ChN}(k).FluoError=ErrorIntegral;
                CompiledParticles{ChN}(k).optFit1=optFit1;
                CompiledParticles{ChN}(k).FitType=FitType;
                CompiledParticles{ChN}(k).noIntensityFlag = noIntensityFlag;
                
                
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
                    if ~isnan(nc9)||~isnan(nc10)||~isnan(nc11)||~isnan(nc12)||~isnan(nc13)||~isnan(nc14)
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
                        for ChN=1:NChannels
                            if isempty(CompiledParticles)==1
                                error(['No compiled particles found in channel ',num2str(ChN),'. Did you mean to run the code with ApproveAll?'])
                            end
                            ncFilter=logical(zeros(length(CompiledParticles{ChN})...
                                ,length(ncFilterID))); %AR 6/16/17: I think multi-channel data might require this to be a cell? Something for the future.
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
                        end
                    end
               
                    figure(2)
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

                    errorbar(ElapsedTime(CompiledParticles{ChN}(k).Frame),...
                        CompiledParticles{ChN}(k).Fluo,ones(size(CompiledParticles{ChN}(k).Fluo))*...
                        CompiledParticles{ChN}(k).FluoError,...
                        '.-r');
                    hold on


                    plot(ElapsedTime(CompiledParticles{ChN}(k).Frame),...
                        CompiledParticles{ChN}(k).Off*IntArea,'.-g');
                    if ~isempty(CompiledParticles{ChN}(k).optFit1)

                        if strcmp(CompiledParticles{ChN}(k).FitType,'spline')
                            SplineValues=ppval(CompiledParticles{ChN}(k).optFit1,double(CompiledParticles{ChN}(k).Frame));
                        elseif strcmp(CompiledParticles{ChN}(k).FitType,'mean')
                            SplineValues=ones(size(CompiledParticles{ChN}(k).Frame))*CompiledParticles{ChN}(k).optFit1;
                        elseif strcmp(CompiledParticles{ChN}(k).FitType,'line')
                            SplineValues=polyval(CompiledParticles{ChN}(k).optFit1,CompiledParticles{ChN}(k).Frame);     
                        end

                        plot([ElapsedTime(CompiledParticles{ChN}(k).Frame)],SplineValues*IntArea,'-k')
                        try
                            title(['Particle ',num2str(k),'(',num2str(i),'), nc',num2str(CompiledParticles{ChN}(k).nc),', Ch: ',num2str(ChN)])
                        catch
                        end
                    else
                        title(['Particle ',num2str(k),'(',num2str(i),'), nc',num2str(CompiledParticles{1}(k).nc),', Ch: ',num2str(ChN),...
                            ' - WARNING: No offset fit'])
                    end
                    hold off
                    %legend({'Particle','Offset','Offset fit'},'Location','Best')
                    xlabel('Time (min)')
                    ylabel('Fluorescence (au)')
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
                                    0,0);
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
                        if strcmp(ExperimentAxis,'AP')
                            title(['Mean AP: ',num2str(CompiledParticles{ChN}(k).MeanAP)])
                        end
                        drawnow
                    end


                    %Snippets
                    for j=1:NFrames
                        subplot(TotalRows,NCols,(TotalRows-NRows)*NCols+j)
                        spotFrame = CompiledParticles{ChN}(k).Frame(j);
                        [x,y,z]=SpotsXYZ(Spots(spotFrame)); 
                        if ~isempty(x)
                            xTrace=x(CompiledParticles{ChN}(k).Index(j));
                            yTrace=y(CompiledParticles{ChN}(k).Index(j));
                            zTrace=z(CompiledParticles{ChN}(k).Index(j));
                        

                        if NChannels==1
                            Image=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                                FilePrefix,iIndex(CompiledParticles{ChN}(k).Frame(j),NDigits),'_z',iIndex(zTrace,2),...
                                '.tif']);
                        else
                            Image=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                                FilePrefix,iIndex(CompiledParticles{ChN}(k).Frame(j),NDigits),'_z',iIndex(zTrace,2),...
                                '_ch',iIndex(ChN,2),'.tif']);
                        end
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
                                    CurrEllipse(2)-yTrace+(SnippetSize-1)/2);
                                %set(EllipseHandle,'color',ColorTime(j,:))
                                set(EllipseHandle,'color','g')
                                hold off
                            end
                        end

                        end
                    end
                    set(gcf,'Position',[1,41,1280,684])  

                    drawnow
                    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'ParticleTraces',filesep,iIndex(k,NDigits),...
                        '(',num2str(i),')-nc',...
                        num2str(CompiledParticles{ChN}(k).nc),'_ch',iIndex(ChN,2),'.tif'])
                    close(2)
                end


                k=k+1;

            end
        end
    end
end
close(h) 





%% Create filters

%nc filters:

if ~isnan(nc9)||~isnan(nc10)||~isnan(nc11)||~isnan(nc12)||~isnan(nc13)||~isnan(nc14)
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
    for ChN=1:NChannels
        if isempty(CompiledParticles)==1
            error(['No compiled particles found in channel ',num2str(ChN),'. Did you mean to run the code with ApproveAll?'])
        end


        ncFilter=logical(zeros(length(CompiledParticles{ChN})...
            ,length(ncFilterID))); %AR 6/16/17: I think multi-channel data might require this to be a cell? Something for the future.
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





    %AP filters:
    if strcmp(ExperimentAxis,'AP')
        %Divide the AP axis into boxes of a certain AP size. We'll see which
        %particle falls where.

        APFilter{ChN}=logical(zeros(length(CompiledParticles{ChN}),length(APbinID)));
        for i=1:length(CompiledParticles{ChN})
            APFilter{ChN}(i,max(find(APbinID<=CompiledParticles{ChN}(i).MeanAP)))=1;
        end
    end
    end
end



%% Binning and averaging data

for ChN=1:NChannels
    %Get the data for the individual particles in a matrix that has the frame
    %number and the particle number as dimensions. Also, get a vector that
    %reports the mean AP position.
    [AllTracesVector{ChN},AllTracesAP{ChN}]=AllTraces(FrameInfo,CompiledParticles{ChN},'NoAP');

    if strcmp(ExperimentAxis,'AP')
        %Mean plot for different AP positions

        %Figure out the AP range to use
        MinAPIndex=1;%min(find(sum(APFilter)));
        MaxAPIndex=size(APFilter{ChN},2);%max(find(sum(APFilter)));

        %Get the corresponding mean information
        k=1;
        for i=MinAPIndex:MaxAPIndex
            [MeanVectorAPTemp,SDVectorAPTemp,NParticlesAPTemp]=AverageTraces(FrameInfo,...
                CompiledParticles{ChN}(APFilter{ChN}(:,i)));
            MeanVectorAPCell{k}=MeanVectorAPTemp';
            SDVectorAPCell{k}=SDVectorAPTemp';
            NParticlesAPCell{k}=NParticlesAPTemp';
            k=k+1;
        end
        MeanVectorAP{ChN}=cell2mat(MeanVectorAPCell);
        SDVectorAP{ChN}=cell2mat(SDVectorAPCell);
        NParticlesAP{ChN}=cell2mat(NParticlesAPCell);
        
        %Calculate the mean for only anterior particles
        try
            MeanVectorAPAnterior{ChN} = MeanVectorAP{ChN}(:,5:15); %Only average particles within window of 10% to 35% w/ 2.5% AP resolution. P2P expression is relatively flat here.
            MeanVectorAnterior{ChN} = nanmean(MeanVectorAPAnterior{ChN},2);
        catch
            %That didn't work
        end
        
    end

    %Calculate the mean for all of them
    [MeanVectorAll{ChN},SDVectorAll{ChN},NParticlesAll{ChN}]=AverageTraces(FrameInfo,CompiledParticles{ChN});
    
    %Now find the different maxima in each nc

    MaxFrame{ChN}=[];
    for i=1:length(NewCyclePos)
        if i==1
            [~,MaxIndex]=max(MeanVectorAll{ChN}(1:NewCyclePos(1)));
            MaxFrame{ChN}=[MaxFrame{ChN},MaxIndex];
        elseif i<=length(NewCyclePos)
            [~,MaxIndex]=max(MeanVectorAll{ChN}(NewCyclePos(i-1):NewCyclePos(i)));
            MaxFrame{ChN}=[MaxFrame{ChN},NewCyclePos(i-1)+MaxIndex-1];
        end
    end
    [~,MaxIndex]=max(MeanVectorAll{ChN}(NewCyclePos(i):end));
    if ~isempty(NewCyclePos)        %Why is this empty sometimes?
                                    %I think this only occurs with suboptimal
                                    %data
        MaxFrame{ChN}=[MaxFrame{ChN},NewCyclePos(i)+MaxIndex-1];
    end
end




%% Instantaneous rate of change
FrameWindow=5;

for ChN=1:NChannels
    %Calculate the derivative as a function of time for each particle
    for j=1:length(CompiledParticles{ChN})
        [CompiledParticles{ChN}(j).SlopeTrace,...
            CompiledParticles{ChN}(j).SDSlopeTrace]=...
            TraceDerivative(CompiledParticles{ChN}(j),...
            ElapsedTime,FrameWindow);
    end


    if strcmp(ExperimentAxis,'AP')
        %Calculate the average slope over an AP window
        MeanSlopeVectorAP{ChN}=nan(size(MeanVectorAP{ChN}));
        SDSlopeVectorAP{ChN}=nan(size(MeanVectorAP{ChN}));
        NSlopeAP{ChN}=nan(size(MeanVectorAP{ChN}));


        APBins=find(sum(APFilter{ChN}));



        for j=1:length(APBins)
            TraceCell=cell(length(ElapsedTime),1);

            ParticlesToAverage=find(APFilter{ChN}(:,APBins(j)));

            for k=1:length(ParticlesToAverage)

                for m=1:length(CompiledParticles{ChN}(ParticlesToAverage(k)).SlopeTrace)

                    TraceCell{CompiledParticles{ChN}(ParticlesToAverage(k)).Frame(m)}=...
                        [TraceCell{CompiledParticles{ChN}(ParticlesToAverage(k)).Frame(m)},...
                        CompiledParticles{ChN}(ParticlesToAverage(k)).SlopeTrace(m)];
                end
            end

            %Get rid of the nan in certain time points
            TraceCell=cellfun(@(x) x(~isnan(x)),TraceCell,'UniformOutput',false);

            MeanTrace=cellfun(@mean,TraceCell,'UniformOutput',false);
            SDTrace=cellfun(@std,TraceCell,'UniformOutput',false);
            NParticlesTrace=cellfun(@length,TraceCell,'UniformOutput',false);

            MeanSlopeVectorAP{ChN}(:,APBins(j))=[MeanTrace{:}];
            SDSlopeVectorAP{ChN}(:,APBins(j))=[SDTrace{:}];
            NSlopeAP{ChN}(:,APBins(j))=[NParticlesTrace{:}];
        end
    end
end



%% Integrating each particle

for ChN=1:NChannels

    %In order to take this seriously I need to come up with a way to deal with
    %the particles that survive during mitosis.

    for i=1:length(CompiledParticles{ChN})
        if length(ElapsedTime(CompiledParticles{ChN}(i).Frame))>1
            CompiledParticles{ChN}(i).TotalmRNA=trapz(ElapsedTime(CompiledParticles{ChN}(i).Frame),CompiledParticles{ChN}(i).Fluo);

            %Estimate the error
            if length(CompiledParticles{ChN}(i).Frame)==2
                CompiledParticles{ChN}(i).TotalmRNAError=(ElapsedTime(CompiledParticles{ChN}(i).Frame(2))-...
                    ElapsedTime(CompiledParticles{ChN}(i).Frame(1)))/2*...
                    CompiledParticles{ChN}(i).FluoError*sqrt(2);
            else

                ErrorTemp=[];
                %Calculate the error of the inner points
                for j=2:(length(CompiledParticles{ChN}(i).Frame)-1)
                     ErrorTemp(j)=(ElapsedTime(CompiledParticles{ChN}(i).Frame(j+1))-...
                         ElapsedTime(CompiledParticles{ChN}(i).Frame(j-1)))*...
                         CompiledParticles{ChN}(i).FluoError;
                end

                %Calculate the error of the outer points
                ErrorTemp(1)=(ElapsedTime(CompiledParticles{ChN}(i).Frame(2))-...
                    ElapsedTime(CompiledParticles{ChN}(i).Frame(1)))/2*...
                    CompiledParticles{ChN}(i).FluoError;

                ErrorTemp(length(CompiledParticles{ChN}(i).Frame))=...
                    (ElapsedTime(CompiledParticles{ChN}(i).Frame(end))-...
                    ElapsedTime(CompiledParticles{ChN}(i).Frame(end-1)))/2*...
                    CompiledParticles{ChN}(i).FluoError;

                %Now, add it all up
                CompiledParticles{ChN}(i).TotalmRNAError=sqrt(sum(ErrorTemp.^2));


            end
        else
            CompiledParticles{ChN}(i).TotalmRNA=[];
        end
    end
end


%% Information about the cytoplasm
%If the nuclear masks are present then use them. Otherwise just calculate
%the median of the images as a function of time

%HG on 8/6/16: Why was this commented out? Did I do this?


% if NChannels==1
% 
%     if HistoneChannel&strcmp(ExperimentAxis,'AP')
%         [MeanCyto,SDCyto,MedianCyto,MaxCyto]=GetCytoMCP(Prefix);
%     else
%         MeanCyto=[];
%         SDCyto=[];
%         MaxCyto=[];
% 
%         h=waitbar(0,'Calculating the median cyto intentisy');
%         for i=1:length(FrameInfo)
%             waitbar(i/length(FrameInfo),h)
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

if NChannels==1

    %Is there a correlation between the fluctuations coming from the offset and
    %those observed in the traces? In order to figure this out I'll fit the
    %nc13 and nc14 intensity data with splines and compute the deviations with
    %respect to them. I'll look into different versions of the data such as
    %with and without the offset subtracted


    if ~SkipFluctuations

        IntArea=109;


        FilteredParticles=find(ncFilter(:,end)|ncFilter(:,end-1));

        OffsetFluct=[];
        DataRawFluct=[];
        %DataOldFluct=[];
        DataSplineFluct=[];

        for j=1:length(FilteredParticles)

            try
                %Deviation from offset with respect to spline
                optFit = adaptiveSplineFit(double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]),...
                        double([CompiledParticles{1}(FilteredParticles(j)).Off*IntArea]),5);
                SplineValues=ppval(optFit,double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]));    



                %Deviation of the raw data, without background subtraction, with
                %respect to a spline.
                DataFitRaw = adaptiveSplineFit(double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]),...
                    double(CompiledParticles{1}(FilteredParticles(j)).FluoRaw),10);
                DataFitRawValues=ppval(DataFitRaw,double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]));


%                 %Deviation of the raw data minus the actual offset with respect to a
%                 %spline
%                 DataFitOld = adaptiveSplineFit(double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]),...
%                     double(CompiledParticles{1}(FilteredParticles(j)).FluoOld),10);
%                 DataSplineValuesOld=ppval(DataFitOld,double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]));



                %Deviation of the raw data minues the spline offset
                DataFit = adaptiveSplineFit(double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]),...
                    double(CompiledParticles{1}(FilteredParticles(j)).Fluo),10);
                DataSplineValues=ppval(DataFit,double([ElapsedTime(CompiledParticles{1}(FilteredParticles(j)).Frame)]));



                %Put all the data together for the plot
                OffsetFluct=[OffsetFluct,CompiledParticles{1}(FilteredParticles(j)).Off*IntArea-SplineValues];
                DataRawFluct=[DataRawFluct,double(CompiledParticles{1}(FilteredParticles(j)).FluoRaw)-DataFitRawValues];
                %DataOldFluct=[DataOldFluct,double(CompiledParticles{1}(FilteredParticles(j)).FluoOld)-DataSplineValuesOld];
                DataSplineFluct=[DataSplineFluct,double(CompiledParticles{1}(FilteredParticles(j)).Fluo)-DataSplineValues];
            end
        end


        xRange=linspace(-4500,4500);

        figure(4)
        plot(OffsetFluct,DataRawFluct,'.k')
        xlabel('Offset fluctuation')
        ylabel('Fluctuations in raw data')
        axis square
        xlim([-4500,4500])
        ylim([-4500,4500])
        R = corrcoef(OffsetFluct,DataRawFluct);
        title(['Correlation: ',num2str(R(2,1))])
        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations',filesep,'Fluct-OffsetVsRawData.tif'])

%         figure(5)
%         plot(OffsetFluct,DataOldFluct,'.k')
%         xlabel('Offset fluctuation')
%         ylabel('Fluctuations with instantaneous offset subtraction')
%         axis square
%         xlim([-4500,4500])
%         ylim([-4500,4500])
%         R = corrcoef(OffsetFluct,DataOldFluct)
%         title(['Correlation: ',num2str(R(2,1))])
%         saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations',filesep,'Fluct-OffsetVsInstData.tif'])



        figure(6)
        plot(OffsetFluct,DataSplineFluct,'.k')
        xlabel('Offset fluctuation')
        ylabel('Fluctuations with spline offset subtraction')
        axis square
        xlim([-4500,4500])
        ylim([-4500,4500])
        R = corrcoef(OffsetFluct,DataSplineFluct);
        title(['Correlation: ',num2str(R(2,1))])
        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations',filesep,'Fluct-OffsetVsSplineData.tif'])
    end



    %Look at the offset of each particle. Do they all look the same? This is
    %only for AP for now
    if strcmp(ExperimentAxis,'AP')
        figure(7)
        subplot(1,2,1)
        MaxAP=0;
        MinAP=inf;
        hold all
        for i=1:length(CompiledParticles{1})
            if sum(CompiledParticles{1}(i).Frame==MaxFrame{1}(end-1))
                MaxAP=max([CompiledParticles{1}(i).MeanAP,MaxAP]);
                MinAP=min([CompiledParticles{1}(i).MeanAP,MinAP]);
                FramePos=find(CompiledParticles{1}(i).Frame==MaxFrame{1}(end-1));
                plot(CompiledParticles{1}(i).MeanAP,CompiledParticles{1}(i).Off(FramePos),'.k')
            end
        end
        hold off
        if MinAP<MaxAP
            xlim([MinAP*0.8,MaxAP*1.2])
        end
        xlabel('AP position')
        ylabel('Offset fluorescence')
        title('Offset at maximum in nc13')
        axis square

        subplot(1,2,2)
        MaxAP=0;
        MinAP=inf;
        hold all
        for i=1:length(CompiledParticles)
            if sum(CompiledParticles{1}(i).Frame==MaxFrame{1}(end))
                MaxAP=max([CompiledParticles{1}(i).MeanAP,MaxAP]);
                MinAP=min([CompiledParticles{1}(i).MeanAP,MinAP]);
                FramePos=find(CompiledParticles{1}(i).Frame==MaxFrame{1}(end));
                plot(CompiledParticles{1}(i).MeanAP,CompiledParticles{1}(i).Off(FramePos),'.k')
            end
        end
        hold off
        if MinAP < Inf && MaxAP > 0
            % ES 2014-09-12: This change is for cases in which no spots were
            % detected during the final time point.
            xlim([MinAP*0.8,MaxAP*1.2]);
        end
        xlabel('AP position')
        ylabel('Offset fluorescence')
        title('Offset at maximum in nc14')
        axis square
        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Offset',filesep,'OffsetVsAP.tif'])
    end


    %Average over all time points

    %Average and SD over each time point. In order to do this we'll generate a
    %cell array with all the values for a given time point

    OffsetCell=cell(length(FrameInfo),1);


    for i=1:length(CompiledParticles{1})
        for j=1:length(CompiledParticles{1}(i).Frame)
            OffsetCell{CompiledParticles{1}(i).Frame(j)}=[OffsetCell{CompiledParticles{1}(i).Frame(j)},...
                CompiledParticles{1}(i).Off(j)];
        end
    end


    MeanOffsetTrace=cellfun(@nanmean,OffsetCell,'UniformOutput',false);
    SDOffsetTrace=cellfun(@nanstd,OffsetCell,'UniformOutput',false);
    NParticlesOffsetTrace=cellfun(@length,OffsetCell,'UniformOutput',false);


    MeanOffsetVector=[MeanOffsetTrace{:}];
    SDOffsetVector=[SDOffsetTrace{:}];
    NOffsetParticles=[NParticlesOffsetTrace{:}];


    if strcmp(ExperimentAxis,'AP')
        figure(8)
        IntArea=109;
        errorbar(1:length(MeanOffsetVector),MeanOffsetVector*IntArea,...
            SDOffsetVector*IntArea,'.-r')
        hold on
        errorbar(1:length(MeanVectorAll{1}),MeanVectorAll{1},...
            SDVectorAll{1},'.-k')
        hold off
        xlabel('Frame')
        ylabel('Fluorescence (au)')
        xlim([0,length(MeanOffsetVector)*1.1])
        legend('Offset','Spots','Location','Best')
        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Offset',filesep,'OffsetAndFluoTime.tif'])


        figure(9)
        errorbar(1:length(MeanOffsetVector),MeanOffsetVector*IntArea-min(MeanOffsetVector*IntArea)+...
            min(MeanVectorAll{1}),...
            SDOffsetVector*IntArea,'.-r')
        hold on
        errorbar(1:length(MeanVectorAll{1}),MeanVectorAll{1},...
            SDVectorAll{1},'.-k')
        hold off
        xlabel('Frame')
        ylabel('Fluorescence (au)')
        xlim([0,length(MeanOffsetVector)*1.1])
        legend('Offset (displaced)','Spots','Location','Best')
        axis square
        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Offset',filesep,'OffsetAndFluoTime-Displaced.tif'])
    end
else
    MeanOffsetVector=[];
    SDOffsetVector=[];
    NOffsetParticles=[];
end

%% Rate of mRNA production

for ChN=1:NChannels

    %Plot the results from fitting the individual traces. Notice that I'm
    %redoing the fits here using the range given by FrameRange. This is because
    %the fluorescence values were calculated slightly differently such that the
    %offset could vary.

    if isfield(CompiledParticles{ChN},'Fit')

        %First, find the maximum in intensity over traces for each cycle
        nc13Max=max([CompiledParticles{ChN}(ncFilter(:,end-1)).Fluo]);
        nc14Max=max([CompiledParticles{ChN}(ncFilter(:,end)).Fluo]);

        figure(10)
        clf

        for i=1:length(CompiledParticles{ChN})
            if ~isempty(CompiledParticles{ChN}(i).Fit)

                %Redo the fit and obtain the parameters

                FrameRange=CompiledParticles{ChN}(i).Fit.FrameRange;

                FramesRangeFilter=ismember(CompiledParticles{ChN}(i).Frame,[FrameRange(1):FrameRange(2)]);

                [a, b, sigma_a, sigma_b] = york_fit(ElapsedTime(CompiledParticles{ChN}(i).Frame(FramesRangeFilter)),...
                            CompiledParticles{ChN}(i).Fluo(FramesRangeFilter),...
                            ones(1,sum(FramesRangeFilter))*mean(diff(ElapsedTime))/2,...
                            ones(1,sum(FramesRangeFilter))*CompiledParticles{ChN}(i).FluoError);


                CompiledParticles{ChN}(i).Fit.Intercept=a;
                CompiledParticles{ChN}(i).Fit.SDIntercept=sigma_a;

                CompiledParticles{ChN}(i).Fit.Slope=b;
                CompiledParticles{ChN}(i).Fit.SDSlope=sigma_b;



                if ncFilter(i,end-1)
                    StartFrame=nc13;
                    EndFrame=nc14;
                    MaxFluo=nc13Max;
                elseif ncFilter(i,end)
                    StartFrame=nc14;
                    EndFrame=length(ElapsedTime);
                    MaxFluo=nc14Max;
                end




                if ~SkipFits



                    xRange=linspace(ElapsedTime(FrameRange(1)),...
                        ElapsedTime(FrameRange(end)));


                    plot(xRange,b*xRange+a,'-k','LineWidth',3)

                    hold on
                    errorbar(ElapsedTime(CompiledParticles{ChN}(i).Frame),CompiledParticles{ChN}(i).Fluo,...
                        ones(1,length(CompiledParticles{ChN}(i).Frame))*CompiledParticles{ChN}(i).FluoError,'.-r')
                    hold off
                    xlabel('Time (min)')
                    ylabel('Fluorescence (au)')
                    title(['Compiled particle ',num2str(i),', nc',num2str(CompiledParticles{ChN}(i).nc)])

                    xlim([ElapsedTime(StartFrame),ElapsedTime(EndFrame)])
                    ylim([0,MaxFluo])

                    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Fits',filesep,'Fit',iIndex(i,3),'-nc',...
                        num2str(CompiledParticles{1}(i).nc),'_ch',iIndex(ChN,2),'.tif'])
                end
            end
        end

        close(10)
    end
end
     

%% First frames

for ChN=1:NChannels

    %How does the first frame from schnitzcells compare to the general one set
    %by just looking at the movie? Do this only for nuclei where we have
    %identified the parent nucleus.

    if HistoneChannel

        figure(11)
        clf
        xRange=linspace(13.5,14.5);
        plot(xRange,ones(size(xRange))*nc14,'-k')
        hold on
        for i=1:length(CompiledParticles{ChN})
            if (CompiledParticles{ChN}(i).nc==14)&(CompiledParticles{ChN}(i).PParticle>0) %#ok<*AND2>
                plot(14*(1+(rand(1)-0.5)/100),CompiledParticles{ChN}(i).NucStart,'.k')
            end
        end

        xRange=linspace(12.5,13.5);
        plot(xRange,ones(size(xRange))*nc13,'-k')
        hold on
        for i=1:length(CompiledParticles{ChN})
            if (CompiledParticles{ChN}(i).nc==13)&(CompiledParticles{ChN}(i).PParticle>0)
                plot(13*(1+(rand(1)-0.5)/100),CompiledParticles{ChN}(i).NucStart,'.k')
            end
        end
        hold off
        try
            ylim([nc13-5,nc14+5])
        end
        set(gca,'XTick',[13,14])
        xlabel('nc')
        ylabel('Frame')
        title('Division time set by eye vs. actual division times of nuclei')
        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Various',filesep,'DivisionTimes_ch',...
            iIndex(ChN,2),'.tif'])



        %Is there any correlation between the first frame and the time of division?


        figure(12)
        subplot(1,2,1)
        hold on
        for i=1:length(CompiledParticles{ChN})
            if (CompiledParticles{ChN}(i).nc==13)&(CompiledParticles{ChN}(i).PParticle>0)
                plot(CompiledParticles{ChN}(i).NucStart*(1+(rand(1)-0.5)/100),CompiledParticles{ChN}(i).FirstFrame,'.k')
            end
        end
        hold off
        axis square
        xlabel('Nuclear birth (frame)')
        ylabel('First particle frame')
        title('nc13')


        subplot(1,2,2)
        hold on
        for i=1:length(CompiledParticles{ChN})
            if (CompiledParticles{ChN}(i).nc==14)&(CompiledParticles{ChN}(i).PParticle>0)
                plot(CompiledParticles{ChN}(i).NucStart*(1+(rand(1)-0.5)/100),CompiledParticles{ChN}(i).FirstFrame,'.k')
            end
        end
        hold off
        axis square
        title('nc14')
        xlabel('Nuclear birth (frame)')
        ylabel('First particle frame')
        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Various',filesep,'DivisionVsFirstFrame_ch',...
            iIndex(ChN,2),'.tif'])





        %First frame and AP position
        if strcmp(ExperimentAxis,'AP')
            figure(13)
            clf
            hold on
            for i=1:length(CompiledParticles{ChN})
                plot(CompiledParticles{ChN}(i).MeanAP,....
                    ElapsedTime(CompiledParticles{ChN}(i).FirstFrame)-...
                    ElapsedTime(nc14),'.k')
            end
            hold off
            box on
            xlabel('AP position (x/L)')
            ylabel('Particle first frame (min)')
    %         if length(ElapsedTime) > nc14+20
    %             ylim([0,ElapsedTime(nc14+20)-ElapsedTime(nc14)])
    %         elseif (ElapsedTime(end) - ElapsedTime(nc14))>0
    %             ylim([0, ElapsedTime(end) - ElapsedTime(nc14)])
    %         end
            % ES 2014-01-05 Testing early-nc14-only movies
            saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Various',filesep,'FirstFrameVsAP_ch',...
                iIndex(ChN,2),'.tif'])
        end
    end
end

%% AP position of particle vs nucleus

if HistoneChannel&&strcmp(ExperimentAxis,'AP')

    for ChN=1:NChannels
        %How different are the AP positions of the nuclei to the particles as a
        %function of time? Let's save the information about nuclear position in the
        %CompiledParticles structure.
        for i=1:length(CompiledParticles{ChN})
            for j=1:length(CompiledParticles{ChN}(i).Frame)
                CurrentNucleus=CompiledParticles{ChN}(i).Nucleus;
                CurrentFrame=CompiledParticles{ChN}(i).Frame(j);

                CurrentEllipse=schnitzcells(CurrentNucleus).cellno(...
                    find((schnitzcells(CurrentNucleus).frames)==CurrentFrame));

                if ~isempty(CurrentEllipse)
                    CompiledParticles{ChN}(i).NuclearAP(j)=EllipsePos{CurrentFrame}(CurrentEllipse);
                else
                    CompiledParticles{ChN}(i).NuclearAP(j)=nan;
                end
            end
        end

        figure(14)
        clf
        hold all
        for i=1:length(CompiledParticles{ChN})
            plot(CompiledParticles{ChN}(i).APpos-CompiledParticles{ChN}(i).NuclearAP)
        end
        hold off
        box on
        xlabel('Frame')
        ylabel('AP difference between spot and nucleus (x/L)')
        saveas(gca,[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'APNucVsParticle_ch',...
            iIndex(ChN,2),'.tif'])
    end
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

if HistoneChannel&&strcmp(ExperimentAxis,'AP')

    for ChN=1:NChannels
        %I'll use the Ellipses structure to count nuclei. This is because
        %schnitzcells sometimes misses things at the edges.


        %First, add the corresponding Ellipse number to each frame of the
        %particles. Also save the information in a cell array. This should make
        %searching easier.
        clear ParticleNuclei
        clear ParticleFrames
        for i=1:length(Particles{ChN})
            if ~isempty(Particles{ChN}(i).Nucleus)
                ParticleNuclei{i}=...
                    schnitzcells(Particles{ChN}(i).Nucleus).cellno(ismember(schnitzcells(Particles{ChN}(i).Nucleus).frames,...
                    Particles{ChN}(i).Frame));
                ParticleFrames{i}=...
                    schnitzcells(Particles{ChN}(i).Nucleus).frames(ismember(schnitzcells(Particles{ChN}(i).Nucleus).frames,...
                    Particles{ChN}(i).Frame));
            else
                ParticleNuclei{i}=[];
                ParticleFrames{i}=[];
            end
        end

        %Do the analogous for CompiledParticles. We'll use this to estimate
        %fluorescence per ALL nuclei
        for i=1:length(CompiledParticles{ChN})
            clear CompiledParticleNuclei
            clear CompiledParticleFrames
            if ~isempty(Particles{ChN}(i).Nucleus)
                CompiledParticleNuclei{i}=...
                    schnitzcells(CompiledParticles{ChN}(i).Nucleus).cellno(ismember(schnitzcells(CompiledParticles{ChN}(i).Nucleus).frames,...
                    CompiledParticles{ChN}(i).Frame));
                CompiledParticleFrames{i}=...
                    schnitzcells(CompiledParticles{ChN}(i).Nucleus).frames(ismember(schnitzcells(CompiledParticles{ChN}(i).Nucleus).frames,...
                    CompiledParticles{ChN}(i).Frame));
            else
                CompiledParticleNuclei{i}=[];
                CompiledParticleFrames{i}=[];
            end
        end


        EdgeWidth=10;
        %For each frame find the number of ellipses that are outside of an area
        %delimited from the edge of the image.
        %The information in Ellipses is
        %(x, y, a, b, theta, maxcontourvalue, time, particle_id)

        %Initialize matrices where we will store the information of number of
        %particles vs. AP vs. time
        NEllipsesAP=zeros(length(Ellipses),length(APbinID));
        NParticlesEllipsesAP{ChN}=zeros(length(Ellipses),length(APbinID));
        %Fluorescence per all of nuclei
        MeanVectorAllAP{ChN}=nan(length(Ellipses),length(APbinID));
        SEVectorAllAP{ChN}=nan(length(Ellipses),length(APbinID));
        



        for i=1:length(Ellipses)
            CurrentEllipses=Ellipses{i};

            Radius=max(CurrentEllipses(:,3:4)')';

            EllipseFilter=(CurrentEllipses(:,1)-Radius-EdgeWidth>0)&...
                (CurrentEllipses(:,2)-Radius-EdgeWidth>0)&...
                (CurrentEllipses(:,1)+Radius+EdgeWidth<FrameInfo(1).PixelsPerLine)&...
                (CurrentEllipses(:,2)+Radius+EdgeWidth<FrameInfo(1).LinesPerFrame);


            %Figure out which particles are in this frame and have an approved
            %flag of 1 or 2. Note that we haven't yet checked if their
            %corresponding Ellipses have been approved.
            CurrentParticlesIndex=cellfun(@(x) find(x==i),ParticleFrames,...
                'UniformOutput',false);
            CurrentParticlesFilter=~cellfun(@isempty,CurrentParticlesIndex);
            ParticlesToCheck=find(CurrentParticlesFilter);


            %Find which of the particles in this frame are related to filtered
            %ellipses and save their corresonding AP information.

            FilteredParticlesPos=[];
            for j=1:length(ParticlesToCheck)
                if EllipseFilter(ParticleNuclei{ParticlesToCheck(j)}(CurrentParticlesIndex{ParticlesToCheck(j)}))&...
                        (Particles{ChN}(ParticlesToCheck(j)).Approved==1 | Particles{ChN}(ParticlesToCheck(j)).Approved==2)
                    FilteredParticlesPos=[FilteredParticlesPos,...
                        EllipsePos{i}(ParticleNuclei{ParticlesToCheck(j)}(CurrentParticlesIndex{ParticlesToCheck(j)}))];
                end
            end

            %Count the number of filtered ellipses per AP bin
            EllipsesFilteredPos{i}=EllipsePos{i}(EllipseFilter);
            for j=1:length(EllipsesFilteredPos{i})
                NEllipsesAP(i,max(find(APbinID<=EllipsesFilteredPos{i}(j))))=...
                    NEllipsesAP(i,max(find(APbinID<=EllipsesFilteredPos{i}(j))))+1;
            end



            %Count the number of filtered particles per AP bin.
            for j=1:length(FilteredParticlesPos)
                NParticlesEllipsesAP{ChN}(i,max(find(APbinID<=FilteredParticlesPos(j))))=...
                    NParticlesEllipsesAP{ChN}(i,max(find(APbinID<=FilteredParticlesPos(j))))+1;
            end




            EllipsesFiltered{ChN}{i}=Ellipses{i}(EllipseFilter,:);

            NEllipsesFiltered{ChN}(i)=sum(EllipseFilter);


            %Calculate the fluorescence per ellipse. Here we'll draw the fluorescence from
            %CompiledParticles just to make sure that everything has been
            %quantified correctly.
            CurrentCompiledParticlesIndex=cellfun(@(x) find(x==i),CompiledParticleFrames,...
                'UniformOutput',false);
            CurrentCompiledParticlesFilter=~cellfun(@isempty,CurrentCompiledParticlesIndex);
            CompiledParticlesToCheck=find(CurrentCompiledParticlesFilter);

            %Find the fluorescence of each set of particles and the AP
            %positions of their Ellipses
            FluorescenceCompiledParticles=[];
            ErrorFluorescenceCompiledParticles=[];
            FilteredCompiledParticlesPos=[];
            for j=1:length(CompiledParticlesToCheck)
                if EllipseFilter(CompiledParticleNuclei{CompiledParticlesToCheck(j)}(CurrentCompiledParticlesIndex{CompiledParticlesToCheck(j)}))&...
                        (CompiledParticles{ChN}(CompiledParticlesToCheck(j)).Approved==1 | CompiledParticles{ChN}(CompiledParticlesToCheck(j)).Approved==2)
                    FilteredCompiledParticlesPos=[FilteredCompiledParticlesPos,...
                        EllipsePos{i}(CompiledParticleNuclei{CompiledParticlesToCheck(j)}(CurrentCompiledParticlesIndex{CompiledParticlesToCheck(j)}))];

                    FluorescenceCompiledParticles=[FluorescenceCompiledParticles,...
                        CompiledParticles{ChN}(CompiledParticlesToCheck(j)).Fluo(...
                            CurrentCompiledParticlesIndex{CompiledParticlesToCheck(j)})];

                    ErrorFluorescenceCompiledParticles=[ErrorFluorescenceCompiledParticles,...
                        CompiledParticles{ChN}(CompiledParticlesToCheck(j)).FluoError];
                end
            end

            %Sum the fluorescence values and divide by the number of ellipses
            for j=1:length(APbinID)
                if ~isnan(APbinArea(j))
                    APFilterTemp=(APbinID(j)<=FilteredCompiledParticlesPos)&...
                        (FilteredCompiledParticlesPos<APbinID(j+1));

                    MeanVectorAllAP{ChN}(i,j)=sum(FluorescenceCompiledParticles(APFilterTemp))/NEllipsesAP(i,j);
                    SEVectorAllAP{ChN}(i,j)=sqrt(sum((ErrorFluorescenceCompiledParticles(APFilterTemp)/NEllipsesAP(i,j)).^2));
                end
            end
        end


        %Minimum number of nuclei to actually do the calculation
        MinNuclei=3;
        MinAPIndexProb=min(find(sum(NEllipsesAP>=MinNuclei)));
        MaxAPIndexProb=max(find(sum(NEllipsesAP>=MinNuclei)));

        MinNucleiFilter=NEllipsesAP>=MinNuclei;

        %Calculate the ratio
        OnRatioAP{ChN}=NParticlesEllipsesAP{ChN}./NEllipsesAP;
        %Filter out the elements that correspond to a number of nuclei below our
        %limit of MinNuclei
        OnRatioAP{ChN}(~MinNucleiFilter)=nan;
        OnRatioAP{ChN}=reshape(OnRatioAP{ChN},size(MinNucleiFilter));


        if MaxAPIndexProb>MinAPIndexProb
            colormap(jet(128));
            cmap=colormap;

            Color=cmap(round((APbinID(MinAPIndexProb:MaxAPIndexProb)-...
                APbinID(MinAPIndexProb))/...
                (APbinID(MaxAPIndexProb)-APbinID(MinAPIndexProb))*127)+1,:);
            figure(15)
            clf
            PlotHandle=[];
            hold on
            for j=MinAPIndexProb:MaxAPIndexProb
                PlotHandle=[PlotHandle,...
                    plot(ElapsedTime,OnRatioAP{ChN}(:,j),'color',Color(j-MinAPIndexProb+1,:))];
            end
            hold off
            xlabel('Time (min)')
            ylabel('Fraction of on nuclei')
            h = colorbar;
            caxis([APbinID(MinAPIndexProb),APbinID(MaxAPIndexProb)])
            ylabel(h,'AP Position (x/L)')
            StandardFigure(PlotHandle,gca)
            xlim([0,ElapsedTime(end)])
            ylim([0,1.01])
            saveas(gca,[DropboxFolder,filesep,Prefix,filesep,'Probabilities',filesep,'ProbVsTimeVsAP.tif'])
        end




        %Now I want to compute the probability of nuclei being on in at least one
        %frame over the whole nc. This is a little bit tricky because I haven't
        %checked the tracking of the nuclei without particles. As a result, those
        %lineages are not complete and we could possibly overcount off nuclei if we
        %just looked at the schnitzcells we have right now. I'll fix this later,
        %but for now I'm going calculate the number of on nuclei per AP bin, where
        %I'll actually check how much area corresponds to that AP bin.





        %What's the probability of a nucleus being on for at least one frame in an nc?
        %This will be an array with rows for nc13 and nc14 and columns
        %corresponding to each AP bin.
        %This only works if we trust the tracking within one nc
        ParticleCountAP{ChN}=zeros(3,length(APbinID));

        try
            for i=1:length(Particles{ChN})
                %See if the particle has either the flag 1 or 2
                if (Particles{ChN}(i).Approved==1)|(Particles{ChN}(i).Approved==2)

                    %Determine the nc so we can add to the right position of ParticleCountAP
                    if (FrameInfo(min(Particles{ChN}(i).Frame(Particles{ChN}(i).FrameApproved))).nc)>=12
                        CurrentNC=FrameInfo(min(Particles{ChN}(i).Frame(Particles{ChN}(i).FrameApproved))).nc;

                        %Now determine the AP bin this particle is on. We'll use
                        %the nucleus positioning. Also, in order to count it we
                        %need to make sure that it fulfills the same criteria we
                        %used to count NEllipsesAP
                        EllipseAPPosTemp=[];
                        EllipsePosXTemp=[];
                        EllipsePosYTemp=[];
                        RadiusTemp=[];
                        xPosTemp=[];
                        yPosTemp=[];
                        for j=1:length(schnitzcells(Particles{ChN}(i).Nucleus).frames)
                            CurrentEllipse=Ellipses{schnitzcells(Particles{ChN}(i).Nucleus).frames(j)}(schnitzcells(Particles{ChN}(i).Nucleus).cellno(j),:);
                            RadiusTemp=[RadiusTemp,max(CurrentEllipse(3:4))];

                            %Determine the AP position
                            EllipseAPPosTemp=[EllipseAPPosTemp,...
                                EllipsePos{schnitzcells(Particles{ChN}(i).Nucleus).frames(j)}(schnitzcells(Particles{ChN}(i).Nucleus).cellno(j))];

                            %Determine the x and y positions
                            xPosTemp=[xPosTemp,CurrentEllipse(1)];
                            yPosTemp=[yPosTemp,CurrentEllipse(2)];
                        end

                        %Make sure that this nucleus was inside the limits at all
                        %points

                        if sum((xPosTemp-RadiusTemp-EdgeWidth>0)&...
                            (yPosTemp-RadiusTemp-EdgeWidth>0)&...
                            (xPosTemp+RadiusTemp+EdgeWidth<FrameInfo(1).PixelsPerLine)&...
                            (yPosTemp+RadiusTemp+EdgeWidth<FrameInfo(1).LinesPerFrame))==...
                            length(xPosTemp)

                            MeanAP=mean(EllipseAPPosTemp);

                            ParticleCountAP{ChN}(CurrentNC-11,max(find(APbinID<=MeanAP)))=...
                                ParticleCountAP{ChN}(CurrentNC-11,max(find(APbinID<=MeanAP)))+1;
                        end
                    end
                end
            end
        end

        %Calculate the probability using the mean number of nuclei per AP bin
        %in each nc. Note that we look at a reduced range within the nc to
        %reduce variability in counting at mitosis.
        ParticleCountProbAP{ChN}(:,1)=ParticleCountAP{ChN}(1,:)./mean(NEllipsesAP(nc12+5:nc13-5,:));
        ParticleCountProbAP{ChN}(:,2)=ParticleCountAP{ChN}(2,:)./mean(NEllipsesAP(nc13+5:nc14-5,:));
        ParticleCountProbAP{ChN}(:,3)=ParticleCountAP{ChN}(3,:)./...
            mean(NEllipsesAP(max(1,nc14-5):length(FrameInfo)-5,:));
        % ES 2014-01-08: accounting for movies started fewer than 5 frames before
        % mitosis 13

        figure(16)   
        plot(APbinID,ParticleCountProbAP{ChN}(:,1),'.-b')
        hold on
        plot(APbinID,ParticleCountProbAP{ChN}(:,2),'.-k')
        plot(APbinID,ParticleCountProbAP{ChN}(:,3),'.-r')
        hold off
        %ylim([0,max(ParticleCountAP(1,:))*2*1.1])
        xlabel('AP position (x/L)')
        ylabel('Active nuclei')
        title('Number active nuclei')
        saveas(gca,[DropboxFolder,filesep,Prefix,filesep,'Probabilities',filesep,'ProbVsAP_ch',...
            iIndex(ChN,2),'.tif'])


        %Use the alternative approach I used for the movies. We are going to
        %look at each nucleus towards the end of each nc and ask if they
        %correspond to an on or off particle in any frame previous to that one.


        %We'll go for 2.5 minutes before the next mitosis. I might relate this
        %later to the elongation time as a way to say that these particles
        %won't contribute to the total amount of mRNA produced anyway.
        FramesBack=ceil(2.5/mean(diff(ElapsedTime)));

        TotalEllipsesAP=zeros(length(APbinID),3);
        EllipsesOnAP{ChN}=zeros(length(APbinID),3);
        for nc=12:14

            %Figure out which frame we'll look at
            if nc==14
                FrameToUse=length(FrameInfo)-FramesBack;
            else
                FrameToUse=eval(['nc',num2str(nc+1)])-FramesBack;
            end

            if FrameToUse>0
                %Filter ellipses that are within the image
                CurrentEllipses=Ellipses{FrameToUse};

                Radius=max(CurrentEllipses(:,3:4)')';

                EllipseFilter=(CurrentEllipses(:,1)-Radius-EdgeWidth>0)&...
                    (CurrentEllipses(:,2)-Radius-EdgeWidth>0)&...
                    (CurrentEllipses(:,1)+Radius+EdgeWidth<FrameInfo(1).PixelsPerLine)&...
                    (CurrentEllipses(:,2)+Radius+EdgeWidth<FrameInfo(1).LinesPerFrame);

                %Check if the filtered ellipses had an associated particle
                EllipsesToCheck=find(EllipseFilter);

                for j=1:length(EllipsesToCheck)
                    %Find which AP bind we're in
                    CurrentAPbin=max(find(APbinID<EllipsePos{FrameToUse}(EllipsesToCheck(j))));
                    %Count the total amount of ellipses in the right AP bin
                    TotalEllipsesAP(CurrentAPbin,nc-11)=TotalEllipsesAP(CurrentAPbin,nc-11)+1;



                    %Find the schnitz this corresponds to
                    for k=1:length(schnitzcells)

                        IndexToUse=find((schnitzcells(k).frames)==FrameToUse);
                        if ~isempty(IndexToUse)

                            %Check this schnitz for consistency with cellno.
                            %Otherwise fix it. I obtained the fixing code from
                            %TrackmRNADynamicsV2.m
                            if length(schnitzcells(k).frames)~=length(schnitzcells(k).cellno)
                               %If there number of frames is different from the number of
                               %cellno then use the cenx and ceny to find the cellno in
                               %Ellipses an repopulate this schnitz
                               if (length(schnitzcells(k).frames)==length(schnitzcells(k).cenx))&...
                                       (length(schnitzcells(k).frames)==length(schnitzcells(k).ceny))
                                   for m=1:length(schnitzcells(k).frames)
                                        %The information in Ellipses is
                                        %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
                                        MaxDistance=2;  %Maximum pixel distance to identify an
                                                        %ellipse with a schnitz
                                        Distances=sqrt((Ellipses{schnitzcells(k).frames(m)}(:,1)-...
                                           schnitzcells(k).cenx(m)).^2+...
                                           (Ellipses{schnitzcells(k).frames(m)}(:,2)-...
                                           schnitzcells(k).ceny(m)).^2);
                                        [MinValue,MinIndex]=min(Distances);

                                       %Make sure no other schnitz is associated to this
                                       %ellipse
                                       EllipseFoundElsewhere=0;
                                       for n=[1:k-1,k+1:length(schnitzcells)]
                                           %Only consider it if the schnitzcell is also valid!
                                           if (length(schnitzcells(n).frames)==length(schnitzcells(n).cellno))
                                               if sum(schnitzcells(n).frames==schnitzcells(k).frames(m))
                                                  IndexToCheck=find(schnitzcells(n).frames==schnitzcells(k).frames(m));

                                                  %The schnitz I'm comparing to
                                                  %might also be screwed up.
                                                  %I'd have to compare its cenx
                                                  %and ceny to be sure
                                                  try
                                                      if schnitzcells(k).cellno(IndexToCheck)==MinIndex
                                                            error('duplicated schnitz?')
                                                            DistancesK=sqrt((Ellipses{schnitzcells(k).frames(m)}(:,1)-...
                                                               schnitzcells(n).cenx(IndexToCheck)).^2+...
                                                               (Ellipses{schnitzcells(k).frames(m)}(:,2)-...
                                                               schnitzcells(n).ceny(IndexToCheck)).^2);          

                                                            [MinValueK,MinIndexK]=min(DistancesK);
                                                            if MinValue<MinValueK
                                                                schnitzcells(n).cellno(IndexToCheck)=[];                                        
                                                            end
                                                        end
                                                  end
                                               end
                                           end
                                        end

                                        if ~EllipseFoundElsewhere
                                           schnitzcells(k).cellno(m)=MinIndex;
                                        else
                                           MinValue
                                           1+1; error('What to do here?')
                                        end
                                   end

                               else
                                   error('Cannnot rescue schnitz')
                               end
                            end

                        end

                        if schnitzcells(k).cellno(IndexToUse)==EllipsesToCheck(j)
                            %Now see if there is an associated particle with it
                            for m=1:length(CompiledParticles{ChN})
                                if CompiledParticles{ChN}(m).Nucleus==k
                                    [j,k,m];
                                    EllipsesOnAP{ChN}(CurrentAPbin,nc-11)=EllipsesOnAP{ChN}(CurrentAPbin,nc-11)+1;
                                end                        
                            end
                        end
                   end
                end
            end
        end

        figure(17)
        plot(APbinID,EllipsesOnAP{ChN}(:,1)./TotalEllipsesAP(:,1),'.-b')
        hold on
        plot(APbinID,EllipsesOnAP{ChN}(:,2)./TotalEllipsesAP(:,2),'.-k')
        plot(APbinID,EllipsesOnAP{ChN}(:,3)./TotalEllipsesAP(:,3),'.-r')
        hold off
        title('Fraction active nuclei')
        xlabel('AP (x/L)')
        ylabel('Fraction')
        legend('nc12', 'nc13', 'nc14')
    end
end




%% Movie of AP profile

%I want to make a movie of the average fluorescence as a function of AP as
%a function of time. In order to make life easier I'll just export to a
%folder. I can then load everything in ImageJ.

if ~SkipMovie&&strcmp(ExperimentAxis,'AP')
    
    for ChN=1:NChannels

        figure(17)

        MaxValue=max(max(MeanVectorAP{ChN}));
        NParticlesAPFilter=NParticlesAP{ChN}>=MinParticles;

        for i=1:length(FrameInfo)
            PlotHandle=errorbar(APbinID(NParticlesAPFilter(i,:)),...
                MeanVectorAP{ChN}(i,NParticlesAPFilter(i,:)),SDVectorAP{ChN}(i,NParticlesAPFilter(i,:)),'.-k');
            hold on
            PlotHandle=[PlotHandle,errorbar(APbinID(NParticlesAPFilter(i,:)),...
                MeanVectorAP{ChN}(i,NParticlesAPFilter(i,:)),...
                SDVectorAP{ChN}(i,NParticlesAPFilter(i,:))./sqrt(NParticlesAP{ChN}(i,NParticlesAPFilter(i,:))),'-k')];
            hold off
            xlim([0.1,0.8])
            ylim([0,MaxValue])
            xlabel('AP position (x/L)')
            ylabel('Mean fluorescence')

            if exist(['nc',num2str(FrameInfo(i).nc)])
                if eval(['nc',num2str(FrameInfo(i).nc)])>0
                    title(['nc',num2str(FrameInfo(i).nc),'. Time into nc: ',num2str(round((ElapsedTime(i)-...
                        ElapsedTime(eval(['nc',num2str(FrameInfo(i).nc)])))*10)/10),' min. Total time: ',...
                        num2str(round(ElapsedTime(i)*10)/10),' min (Frame ',num2str(i),').'])
                else
                    title(['nc',num2str(FrameInfo(i).nc),'. Total time: ',...
                        num2str(round(ElapsedTime(i)*10)/10),' min (Frame ',num2str(i),').'])
                end
            else
                title(['nc',num2str(FrameInfo(i).nc),'. Total time: ',...
                    num2str(round(ElapsedTime(i)*10)/10),' min (Frame ',num2str(i),').'])
            end

            StandardFigure(PlotHandle,gca)
            saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'APMovie',filesep,iIndex(i,3),'_ch',iIndex(ChN,2),'.tif']);   
        end
    end
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
if strcmp(lower(ExperimentType),'inputoutput') && strcmp(lower(ExperimentAxis),'ap')
    CompileNuclearProtein(Prefix)
end



%% Save everything
if HistoneChannel&strcmp(ExperimentAxis,'AP')

    %If we have only one channel get rid of all the cells
    if NChannels==1
        CompiledParticles=CompiledParticles{1};
        APFilter=APFilter{1};
        MeanVectorAP=MeanVectorAP{1};
        SDVectorAP=SDVectorAP{1};
        NParticlesAP=NParticlesAP{1};
        MeanVectorAll=MeanVectorAll{1};
        MeanVectorAnterior = MeanVectorAnterior{1};
        SDVectorAll=SDVectorAll{1};
        NParticlesAll=NParticlesAll{1};
        MaxFrame=MaxFrame{1};
        AllTracesVector=AllTracesVector{1};
        AllTracesAP=AllTracesAP{1};
        MeanSlopeVectorAP=MeanSlopeVectorAP{1};
        SDSlopeVectorAP=SDSlopeVectorAP{1};
        NSlopeAP=NSlopeAP{1};
        ParticleCountAP={1};
        OnRatioAP=OnRatioAP{1};
        ParticleCountProbAP=ParticleCountProbAP{1};
        EllipsesOnAP=EllipsesOnAP{1};
        MeanVectorAllAP=MeanVectorAllAP{1};
        SEVectorAllAP=SEVectorAllAP{1};    
        SNR = MeanVectorAll ./ (SDVectorAll*109); %AR 7/12/16: Rough estimate. Should do better.
    end
    
    save([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],...
        'CompiledParticles','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
        'nc12','nc13','nc14','ncFilterID','StemLoopEnd','ncFilter','APbinID','APFilter',...
        'MeanVectorAP','SDVectorAP','NParticlesAP','MeanVectorAll','MeanVectorAnterior',...
        'SDVectorAll','NParticlesAll','MaxFrame','MinAPIndex','MaxAPIndex',...
        'AllTracesVector','AllTracesAP','MeanCyto','SDCyto','MedianCyto','MaxCyto',...
        'MeanOffsetVector','SDOffsetVector','NOffsetParticles',...
        'MeanSlopeVectorAP','SDSlopeVectorAP','NSlopeAP',...
        'ParticleCountAP','APbinArea','OnRatioAP','NEllipsesAP',...
        'ParticleCountProbAP',...
        'EllipsesOnAP','TotalEllipsesAP',...
        'EllipsePos','EllipsesFilteredPos','FilteredParticlesPos',...
        'MeanVectorAllAP','SEVectorAllAP', 'Prefix');
elseif HistoneChannel&strcmp(ExperimentAxis,'DV')
    
    %If we have only one channel get rid of all the cells
    if NChannels==1
        CompiledParticles=CompiledParticles{1};
        MeanVectorAll=MeanVectorAll{1};
        SDVectorAll=SDVectorAll{1};
        NParticlesAll=NParticlesAll{1};
        MaxFrame=MaxFrame{1};
        AllTracesVector=AllTracesVector{1};
    end
    
    
    save([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],...
        'CompiledParticles','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
        'nc12','nc13','nc14','StemLoopEnd','ncFilterID','ncFilter',...
        'MeanVectorAll',...
        'SDVectorAll','NParticlesAll','MaxFrame',...
        'AllTracesVector','MeanCyto','SDCyto','MedianCyto','MaxCyto',...
        'MeanOffsetVector','SDOffsetVector','NOffsetParticles', 'Prefix')
elseif strcmp(ExperimentAxis,'NoAP')
    
    %If we have only one channel get rid of all the cells
    if NChannels==1
        CompiledParticles=CompiledParticles{1};
        MeanVectorAll=MeanVectorAll{1};
        SDVectorAll=SDVectorAll{1};
        NParticlesAll=NParticlesAll{1};
        MaxFrame=MaxFrame{1};
        AllTracesVector=AllTracesVector{1};
    end
    try
        save([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],...
            'CompiledParticles','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
            'nc12','nc13','nc14','StemLoopEnd','ncFilterID','ncFilter',...
            'MeanVectorAll',...
            'SDVectorAll','NParticlesAll','MaxFrame',...
            'AllTracesVector','MeanCyto','SDCyto','MedianCyto','MaxCyto',...
            'MeanOffsetVector','SDOffsetVector','NOffsetParticles')
    catch
        save([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],...
            'CompiledParticles','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
            'nc12','nc13','nc14','StemLoopEnd',...
            'MeanVectorAll',...
            'SDVectorAll','NParticlesAll','MaxFrame',...
            'AllTracesVector','MeanCyto','SDCyto','MedianCyto','MaxCyto',...
            'MeanOffsetVector','SDOffsetVector','NOffsetParticles', 'Prefix')
    end
else
    
    %If we have only one channel get rid of all the cells
    if NChannels==1
        CompiledParticles=CompiledParticles{1};
        APFilter=APFilter{1};
        MeanVectorAP=MeanVectorAP{1};
        SDVectorAP=SDVectorAP{1};
        NParticlesAP=NParticlesAP{1};
        MeanVectorAll=MeanVectorAll{1};
        SDVectorAll=SDVectorAll{1};
        NParticlesAll=NParticlesAll{1};
        MaxFrame=MaxFrame{1};
        AllTracesVector=AllTracesVector{1};
        AllTracesAP=AllTracesAP{1};
        MeanSlopeVectorAP=MeanSlopeVectorAP{1};
        SDSlopeVectorAP=SDSlopeVectorAP{1};
        NSlopeAP=NSlopeAP{1};
        ParticleCountAP={1};
        OnRatioAP=OnRatioAP{1};
        ParticleCountProbAP=ParticleCountProbAP{1};
        EllipsesOnAP=EllipsesOnAP{1};
        MeanVectorAllAP=MeanVectorAllAP{1};
        SEVectorAllAP=SEVectorAllAP{1};    
    end
    
    save([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],...
        'CompiledParticles','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
        'nc12','nc13','nc14','StemLoopEnd','ncFilterID','ncFilter','APbinID','APFilter',...
        'MeanVectorAP','SDVectorAP','NParticlesAP','MeanVectorAll',...
        'SDVectorAll','NParticlesAll','MaxFrame','MinAPIndex','MaxAPIndex',...
        'AllTracesVector','AllTracesAP','MeanCyto','SDCyto','MedianCyto','MaxCyto',...
        'MeanOffsetVector','SDOffsetVector','NOffsetParticles',...
        'MeanSlopeVectorAP','SDSlopeVectorAP','NSlopeAP', 'Prefix')
end