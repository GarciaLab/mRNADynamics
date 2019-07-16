function CompileNuclearProtein(Prefix, varargin)
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

% Default values for Options
ROI=false; % No ROI, as a default
nonROI=false;
xROI = false;
yROI = false;
NameString_ROI='';

%This function will add fluorescence information to each schnitz.

close all

%gives the location of these 5 quantities.
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;


for args=1:length(varargin)
    if strcmp(varargin{args},'ROI')|| strcmp(varargin{args},'nonROI')
        ROI = true;
        % Find the direction of ROI
        if varargin{args+1} == 'x'
            xROI = true;
        elseif varargin{args+1} == 'y'
            yROI = true;
        end
        
        % Find the region of ROI
        if ~isnumeric(varargin{args+2})||~isnumeric(varargin{args+2})
            error('Wrong input parameters. After ''ROI'' you should input the threshold of ROI, [Coord1_ROI Coor2_ROI] ')
        else
            ROI1=varargin{args+2}(1);
            ROI2=varargin{args+2}(2);
        end
    end
end

% Saving Nomenclature
if ROI
    NameString_ROI = 'ROI';
elseif nonROI
    NameString_ROI = 'nonROI';
else
    NameString_ROI = '';
end


[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);




%Load all the information
load([DropboxFolder,filesep,Prefix,'\CompiledParticles.mat'], 'CompiledParticles')
load([DropboxFolder,filesep,Prefix,'\Ellipses.mat'], 'Ellipses')
load([DropboxFolder,filesep,Prefix,'\FrameInfo.mat'], 'FrameInfo')
load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'schnitzcells')

numFrames = length(FrameInfo);

%See if FrameInfo has information about the number of input channels. This
%is not fully implemented yet. If no information is found, then assume we
%have only one input channel.
if isfield(FrameInfo,'NChInput')
    NChannels=FrameInfo(1).NChInput;
else
    NChannels=1;
end

%What type of experiment are we dealing with? Get this out of
%MovieDatabase.xlsx
% [SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
%     Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
%     ] = readMovieDatabase(Prefix);

% refactor in progress, we should replace readMovieDatabase with getExperimentDataFromMovieDatabase
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

%Pre-calculating ExperimentAxis boolean for faster use in later if statements
ExperimentAxisIsNoAP = strcmpi(ExperimentAxis, 'NoAP');

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

NewCyclePos=[nc9,nc10,nc11,nc12,nc13,nc14];
NewCyclePos=NewCyclePos(~(NewCyclePos==0));
NewCyclePos=NewCyclePos(~isnan(NewCyclePos));




%Add the APPosition to Particles if they don't exist yet
if (~isfield(schnitzcells,'APpos'))&&(strcmpi(ExperimentAxis,'AP')||strcmpi(ExperimentAxis,'DV'))
    %First, run this to get the alignment between the zoom-in and zoom-out
    %images:
%     AddParticlePosition(Prefix)
    %Now, add the nuclear position
    AddNuclearPosition(Prefix)
    load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
end



% Put together CompiledNuclei

%Get the actual time corresponding to each frame
if isfield(FrameInfo,'FileMode')
    if strcmp(FrameInfo(1).FileMode,'TIF')%Is this a TIF file if not is it a LSM or LIFE Export
        for j=1:numFrames
            ElapsedTime(j)=etime(datevec(FrameInfo(j).TimeString),datevec(FrameInfo(1).TimeString));
        end
    elseif strcmp(FrameInfo(1).FileMode,'LSM')||strcmp(FrameInfo(1).FileMode,'LIFExport')||strcmp(FrameInfo(1).FileMode, 'LSMExport')% SEANCHANGE If it is a LIFEexport Sum over Fram e info
        for j=1:numFrames
            ElapsedTime(j)=FrameInfo(j).Time-FrameInfo(1).Time;%Finds the elapsed time by subtracting each time point by the initial time point
        end
    else
        error('File mode not supported. Cannot extract time information. Include format in ExportDataForLivemRNA.m')
    end
else
    warning('No FileMode information found. Assuming that this is TIF from the 2-photon.')
    for j=1:numFrames
        ElapsedTime(j)=etime(datevec(FrameInfo(j).TimeString),datevec(FrameInfo(1).TimeString));
    end
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




%Now get the nuclear information for those that were approved
NZSclices=size(schnitzcells(1).Fluo,2);

%CompiledNuclei(length(schnitzcells))=struct;

h=waitbar(0,'Compiling nuclear traces');
k=1;
for i=1:length(schnitzcells)
    

    waitbar(i/length(schnitzcells),h)
    
    if schnitzcells(i).Approved
        %Which frames were approved manually?
        FrameFilter=schnitzcells(i).FrameApproved;
        
        %Check that for the remaining frames we got a good z-profile
        for j=1:length(schnitzcells(i).frames)
            [MaxValue,MaxPos]=max(schnitzcells(i).Fluo(j,:));
            if NZSclices<3
                if (MaxPos==2)||(MaxPos==NZSclices-1)
                    FrameFilter(j)=false;
                end
            end
        end
        
        if sum(FrameFilter)
        
            %Copy the filtered information
            CompiledNuclei(k).P=schnitzcells(i).P;
            CompiledNuclei(k).E=schnitzcells(i).E;
            CompiledNuclei(k).D=schnitzcells(i).D;
            CompiledNuclei(k).Frames=uint16(schnitzcells(i).frames(FrameFilter));
            CompiledNuclei(k).xPos=single(schnitzcells(i).cenx(FrameFilter));
            CompiledNuclei(k).yPos=single(schnitzcells(i).ceny(FrameFilter));
            CompiledNuclei(k).Radius=single(schnitzcells(i).len(FrameFilter));
            CompiledNuclei(k).cellno=uint16(schnitzcells(i).cellno(FrameFilter));
            CompiledNuclei(k).nc=[];
            
            %Save the information about the original schnitz
            CompiledNuclei(k).schnitz=uint16(i);
            
            if ~ExperimentAxisIsNoAP
                CompiledNuclei(k).MeanDV=single(mean(schnitzcells(i).DVpos(FrameFilter)));
                CompiledNuclei(k).MedianDV=single(median(schnitzcells(i).DVpos(FrameFilter)));
                CompiledNuclei(k).MeanAP=single(mean(schnitzcells(i).APpos(FrameFilter)));
                CompiledNuclei(k).MedianAP=single(median(schnitzcells(i).APpos(FrameFilter)));
            end
            
           FluoTimeTrace = ExtractDlFluo(schnitzcells(i).Fluo, .5);
           CompiledNuclei(k).FluoTimeTrace=single(FluoTimeTrace);

            %Copy and extract the fluorescence information
            CompiledNuclei(k).FluoMax=single(squeeze(max(schnitzcells(i).Fluo(FrameFilter,:,:),[],2)));
            %For DV case, need to calculate more
            if strcmpi(ExperimentAxis,'DV')
                DV_Fluo = schnitzcells(i).Fluo(FrameFilter,:,:);
                CompiledNuclei(k).DVFluo = single(DV_Fluo);
                CompiledNuclei(k).FluoMin=single(squeeze(min(DV_Fluo(DV_Fluo>0),[],2)));
                CompiledNuclei(k).FluoMean=single(squeeze(mean(DV_Fluo(DV_Fluo>0),2)));
                %the 'p' field is a parabola fit of DV fluorescence over
                %time
%                 for frame = 1:size(DV_Fluo, 1)
%                     CompiledNuclei(k).parabolaFit{frame} = polyfit(1:size(DV_Fluo,2),DV_Fluo(frame,:),2);
%                 end
                    
            end

            %If there was only one time point and multiple channels,
            %squeeze can lead to a weird shape of the matrix
            if (NChannels>1)&&(size(CompiledNuclei(k).FluoMax,2)==1)
                CompiledNuclei(k).FluoMax=CompiledNuclei(k).FluoMax';                
            end
            
            k=k+1;
        end
    end
end
close(h)     

%% ROI option
% This option is sorting the CompiledNuclei that are within the ROI region
% that is defined by direction, [ROI1 ROI2]
% written by YJK on 6/6/2019

if ROI
    % Sort the CompiledNuclei using the threshold
    t=1;
    
    % Use the ROI1 and ROI2 to split the Nuclei
    for NucleiIndex=1:length(CompiledNuclei)
        if xROI
            if nanmean(CompiledNuclei(NucleiIndex).xPos) > ROI1 &&...
                    nanmean(CompiledNuclei(NucleiIndex).xPos) < ROI2
                CompiledNuclei_ROI(t)=CompiledNuclei(NucleiIndex);
                t=t+1;
            end
        elseif yROI
            if nanmean(CompiledNuclei(NucleiIndex).yPos) > ROI1 &&...
                    nanmean(CompiledNuclei(NucleiIndex).yPos) < ROI2
                CompiledNuclei_ROI(t)=CompiledNuclei(NucleiIndex);
                t=t+1;
            end
        end
    end
    % Redefine the CompiledNuclei sorted out.
    CompiledNuclei = CompileNuclei_ROI;
end


%% Create AP and nc filters

%nc filters:

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
ncFilter=false(length(CompiledNuclei),length(ncFilterID));
for i=1:length(CompiledNuclei)
    if ~isempty(CompiledNuclei(i).Frames)
    
        if ~isempty(CompiledNuclei(i).nc)
            ncFilter(i,find(CompiledNuclei(i).nc==ncFilterID))=true;
        else
            ncsFound=find(CompiledNuclei(i).Frames(1)>=[nc9,nc10,nc11,nc12,nc13,nc14]);
            if ncsFound(end)==1
                CompiledNuclei(i).nc=9;
                ncFilter(i,ncFilterID==9)=true;
            elseif ncsFound(end)==2
                CompiledNuclei(i).nc=10;
                ncFilter(i,ncFilterID==10)=true;
            elseif ncsFound(end)==3
                CompiledNuclei(i).nc=11;
                ncFilter(i,ncFilterID==11)=true;
            elseif ncsFound(end)==4
                CompiledNuclei(i).nc=12;
                ncFilter(i,ncFilterID==12)=true;
            elseif ncsFound(end)==5
                CompiledNuclei(i).nc=13;
                ncFilter(i,ncFilterID==13)=true;
            elseif ncsFound(end)==6
                CompiledNuclei(i).nc=14;
                ncFilter(i,ncFilterID==14)=true;
            end

        end
    end
end



if strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')
    %AP filters:

    %Divide the AP axis into boxes of a certain AP size. We'll see which
    %particle falls where.

    APbinID=0:APResolution:1;

    APFilter=false(length(CompiledNuclei),length(APbinID));
    
    for i=1:length(CompiledNuclei)
        APFilter(i,max(find(APbinID<=CompiledNuclei(i).MeanAP)))=true;
    end

end


%% Information about the cytoplasmic fluroescence
%If the nuclear masks are present then use them. Otherwise just calculate
%the median of the images as a function of time
MeanCyto=[];
SDCyto=[];
MaxCyto=[];
MedianCyto = [];
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

if strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')
    
%Get the data for the individual particles in a matrix that has the frame
%number and the particle number as dimensions. Also, get a vector that
%reports the mean position.
[AllTracesVector,AllTracesAP, AllTracesDV]=AllTracesNuclei(FrameInfo,CompiledNuclei);



    %Mean plot for different AP positions

    %Figure out the AP range to use
    MinAPIndex=1;%min(find(sum(APFilter)));
    MaxAPIndex=size(APFilter,2);%max(find(sum(APFilter)));
    
    %Get the corresponding mean information 
    k=1;
    for ap=MinAPIndex:MaxAPIndex
        [MeanVectorAPTemp,SDVectorAPTemp,NParticlesAPTemp]=AverageTracesNuclei(FrameInfo,...
            CompiledNuclei(APFilter(:,ap)),NChannels);
        MeanVectorAPCell{k}=MeanVectorAPTemp';
        SDVectorAPCell{k}=SDVectorAPTemp';
        NParticlesAPCell{k}=NParticlesAPTemp';
        k=k+1;
    end

    %Turn the information into useful structures
    if NChannels>1
        for ch=1:NChannels
            for ap=MinAPIndex:MaxAPIndex
                MeanVectorAPCell2{ch,ap}=MeanVectorAPCell{ap}{ch};
                SDVectorAPCell2{ch,ap}=SDVectorAPCell{ap}{ch};
                NParticlesAPCell2{ch,ap}=NParticlesAPCell{ap}{ch};
            end
        end

        for ch=1:NChannels
            MeanVectorAP{ch}=cell2mat({MeanVectorAPCell2{ch,:}}')';
            SDVectorAP{ch}=cell2mat({SDVectorAPCell2{ch,:}}')';
            NParticlesAP{ch}=cell2mat({NParticlesAPCell2{ch,:}}')';
        end
    else
        for ap=MinAPIndex:MaxAPIndex
            MeanVectorAPCell2{ap}=MeanVectorAPCell{ap};
            SDVectorAPCell2{ap}=SDVectorAPCell{ap};
            NParticlesAPCell2{ap}=NParticlesAPCell{ap};
        end

        MeanVectorAP=cell2mat(MeanVectorAPCell2);
        SDVectorAP=cell2mat(SDVectorAPCell2);
        NParticlesAP=cell2mat(NParticlesAPCell2);
    end


elseif strcmpi(ExperimentAxis,'NoAP')
    %Get the data for the individual particles in a matrix that has the frame
    %number and the particle number as dimensions. Also, get a vector that
    %reports the mean position.
    [AllTracesVector,AllTracesAP, AllTracesDV]=AllTracesNuclei(FrameInfo,CompiledNuclei,'NoAP');
end
    



%Calculate the mean for all of them
[MeanVectorAll,SDVectorAll,NParticlesAll]=AverageTracesNuclei(FrameInfo,CompiledNuclei);


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

%% Save everything

savedVariables = [savedVariables,...
            'CompiledNuclei','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
            'nc12','nc13','nc14','ncFilterID','ncFilter','APbinID','APFilter',...
            'MeanVectorAP','SDVectorAP','NParticlesAP',...
            'MeanVectorAll','SDVectorAll','NParticlesAll',...
            'MaxFrame',...
            'AllTracesVector','AllTracesAP',...
            'MeanCyto','SDCyto','MedianCyto','MaxCyto',...
            'MeanCytoAPProfile','SDCytoAPProfile','SECytoAPProfile',...
            'IntegrationArea'...
            'DVbinID','DVFilter','MeanVectorDV','SDVectorDV','NParticlesDV',...
                    'MinDVIndex','MaxDVIndex', 'AllTracesDV'];


save([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat',NameString_ROI],...
        savedVariables{:},'-v7.3');

save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells', '-v7.3')
