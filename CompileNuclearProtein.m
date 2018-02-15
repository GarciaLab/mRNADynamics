function CompileNuclearProtein(varargin)
%This code gives me the input
%varargin Variable length input argument list.
%allows any number of arguments to a function.  The variable
%varargin is a cell array containing the optional arguments to the
%function.  varargin must be declared as the last input argument
%and collects all the inputs from that point onwards.

% ROI (Region of Interest, illuminatio) option is added on 10/27/2017 by YJK
% ex ) CompileNuclearProtein(Prefix,'ROI',ROI1,ROI2)
% Assume that the ROI is the top-half of the imaging window,
% ROI1 is the lower boundary of the ROI, ROI2 is the upper boundary of the
% non-ROI 
% Note. The main change is in the 'APbin filter' and 'binning and
% averaging' part of the code.

% Default values for Options
ROI=0; % No ROI, as a default

%This function will add fluorescence information to each schnitz.

close all

%gives the location of these 5 quantities.
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

%Look at the input parameter and use defaults if missing
if isempty(varargin)% returns 1 if it is an empty array and 0 otherwise.an empty array has prod(size(X))==0
    FolderTemp=uigetdir(DropboxFolder,'Select folder with data to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
else
    Prefix=varargin{1};
    for i=2:length(varargin)
        % Adding the option of ROI, and two thresholds like in
        % ComileParticls (YJK : on 10/27/2017)
        if strcmp(varargin{i},'ROI')
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

[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);




%Load all the information
load([DropboxFolder,filesep,Prefix,'\Ellipses.mat'])
load([DropboxFolder,filesep,Prefix,'\FrameInfo.mat'])
load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])

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
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
    Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
    ] = readMovieDatabase(Prefix);

% refactor in progress, we should replace readMovieDatabase with getExperimentDataFromMovieDatabase
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

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
if (~isfield(schnitzcells,'APpos'))&&(strcmpi(ExperimentAxis,'ap'))
    %First, run this to get the alignment between the zoom-in and zoom-out
    %images:
    AddParticlePosition(Prefix)
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
    elseif strcmp(FrameInfo(1).FileMode,'LSM')||strcmp(FrameInfo(1).FileMode,'LIFExport')%If it is a LIFEexport Sum over Fram e info
        for j=1:numFrames
            ElapsedTime(j)=FrameInfo(j).Time-FrameInfo(1).Time;%Finds the elapsed time by subtracting each time point by the initial time point
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
    
    if (schnitzcells(i).Approved==1)
        %Which frames were approved manually?
        FrameFilter=schnitzcells(i).FrameApproved;
        
        %Check that for the remaining frames we got a good z-profile
        for j=1:length(schnitzcells(i).frames)
            [MaxValue,MaxPos]=max(schnitzcells(i).Fluo(j,:));
            if (MaxPos==2)||(MaxPos==NZSclices-1)
                FrameFilter(j)=0;
            end
        end
        
        if sum(FrameFilter)
        
            %Copy the filtered information
            CompiledNuclei(k).P=schnitzcells(i).P;
            CompiledNuclei(k).E=schnitzcells(i).E;
            CompiledNuclei(k).D=schnitzcells(i).D;
            CompiledNuclei(k).Frames=schnitzcells(i).frames(FrameFilter);
            CompiledNuclei(k).xPos=schnitzcells(i).cenx(FrameFilter);
            CompiledNuclei(k).yPos=schnitzcells(i).ceny(FrameFilter);
            CompiledNuclei(k).Radius=schnitzcells(i).len(FrameFilter);
            CompiledNuclei(k).cellno=schnitzcells(i).cellno(FrameFilter);
            CompiledNuclei(k).nc=[];
            
            %Save the information about the original schnitz
            CompiledNuclei(k).schnitz=i;

            
            if strcmpi(ExperimentAxis,'ap')
                %Determine the particles average and median AP position
                CompiledNuclei(k).MeanAP=mean(schnitzcells(i).APpos(FrameFilter));
                CompiledNuclei(k).MedianAP=median(schnitzcells(i).APpos(FrameFilter));
            end

            %Copy and extract the fluorescence information
            CompiledNuclei(k).FluoMax=squeeze(max(schnitzcells(i).Fluo(FrameFilter,:,:),[],2));
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
% This option is separating the CompiledNuclei defined above into
% CompiledNuclei_ROI and COmpiledNuclei_nonROI
% written by YJK on 10/27/2017

if ROI
    % separate the CompiledNuclei into CompiledNuclei_ROI and
    % Particles_nonROI using Threshold
    t=1;
    s=1;

    % Use the ROI1 and ROI2 to split the Particles
    for NucleiIndex=1:length(CompiledNuclei)
        if nanmean(CompiledNuclei(NucleiIndex).yPos) < ROI1
            CompiledNuclei_ROI(t)=CompiledNuclei(NucleiIndex);
            t=t+1;
        elseif nanmean(CompiledNuclei(NucleiIndex).yPos) > ROI2
            CompiledNuclei_nonROI(s)=CompiledNuclei(NucleiIndex);
            s=s+1;
        end
    end
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



if strcmpi(ExperimentAxis,'ap')
    %AP filters:

    %Divide the AP axis into boxes of a certain AP size. We'll see which
    %particle falls where.

    APResolution=APResolution;
    APbinID=0:APResolution:1;

    APFilter=false(length(CompiledNuclei),length(APbinID));
    
    if ROI==1
        %Define two APFilters for ROI and non-ROI respectively
        APFilter_ROI=false(length(CompiledNuclei_ROI),length(APbinID));
        APFilter_nonROI=false(length(CompiledNuclei_nonROI),length(APbinID));
        APFilter=false(length(CompiledNuclei),length(APbinID));

        for i=1:length(CompiledNuclei)
            APFilter(i,max(find(APbinID<=CompiledNuclei(i).MeanAP)))=1;
        end

        for i=1:length(CompiledNuclei_ROI)
            APFilter_ROI(i,max(find(APbinID<=CompiledNuclei_ROI(i).MeanAP)))=1;
        end

        for i=1:length(CompiledNuclei_nonROI)
            APFilter_nonROI(i,max(find(APbinID<=CompiledNuclei_nonROI(i).MeanAP)))=1;
        end
    
    else
        for i=1:length(CompiledNuclei)
            APFilter(i,max(find(APbinID<=CompiledNuclei(i).MeanAP)))=1;
        end
    end
end


%% Information about the cytoplasmic fluroescence
%If the nuclear masks are present then use them. Otherwise just calculate
%the median of the images as a function of time

if strcmp(ExperimentAxis,'AP')
    [MeanCyto,SDCyto,MedianCyto,MaxCyto,...
        MeanCytoAPProfile,SDCytoAPProfile,SECytoAPProfile]=GetCytoMCP(Prefix);
else
    MeanCyto=[];
    SDCyto=[];
    MaxCyto=[];
    if (~isempty(strfind(Channel1{1}, 'Bcd')))
        nameSuffix=['_ch',iIndex(1,2)];
    else
        nameSuffix=['_ch',iIndex(2,2)]; %assuming input channel is channel two. obviously this is going to be wrong in general.
    end
    
    h=waitbar(0,'Calculating the median cyto intentisy');
    for i=1:numFrames
        waitbar(i/numFrames,h)
        for j=1:FrameInfo(1).NumberSlices
            Image(:,:,j)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(i,3),'_z',iIndex(j,2),nameSuffix,'.tif']);
        end
        ImageMax=max(Image,[],3);
        MedianCyto(i)=median(double(ImageMax(:)));
    end
    close(h)
end  


%% Binning and averaging data


if strcmpi(ExperimentAxis,'ap')

   
    %Get the data for the individual particles in a matrix that has the frame
    %number and the particle number as dimensions. Also, get a vector that
    %reports the mean AP position.
    [AllTracesVector,AllTracesAP]=AllTracesNuclei(FrameInfo,CompiledNuclei);

    %Mean plot for different AP positions

    %Figure out the AP range to use
    MinAPIndex=1;%min(find(sum(APFilter)));
    MaxAPIndex=size(APFilter,2);%max(find(sum(APFilter)));
    
    if ROI
        % Mean values for ROI region
        k=1;
        for i=MinAPIndex:MaxAPIndex
            [MeanVectorAPTemp_ROI,SDVectorAPTemp_ROI,NParticlesAPTemp_ROI]=AverageTracesNuclei(FrameInfo,...
                CompiledNuclei_ROI(APFilter_ROI(:,i)),NChannels);
            MeanVectorAPCell_ROI{k}=MeanVectorAPTemp_ROI';
            SDVectorAPCell_ROI{k}=SDVectorAPTemp_ROI';
            NParticlesAPCell_ROI{k}=NParticlesAPTemp_ROI';
            k=k+1;
        end

        %Turn the information into useful structures
        if NChannels>1
            for j=1:NChannels
                for i=MinAPIndex:MaxAPIndex
                    MeanVectorAPCell2_ROI{j,i}=MeanVectorAPCell_ROI{i}{j};
                    SDVectorAPCell2_ROI{j,i}=SDVectorAPCell_ROI{i}{j};
                    NParticlesAPCell2_ROI{j,i}=NParticlesAPCell_ROI{i}{j};
                end
            end

            for j=1:NChannels
                MeanVectorAP_ROI{j}=cell2mat({MeanVectorAPCell2_ROI{j,:}}')';
                SDVectorAP_ROI{j}=cell2mat({SDVectorAPCell2_ROI{j,:}}')';;
                NParticlesAP_ROI{j}=cell2mat({NParticlesAPCell2_ROI{j,:}}')';;
            end
        else
            for i=MinAPIndex:MaxAPIndex
                MeanVectorAPCell2_ROI{j,i}=MeanVectorAPCell_ROI{i};
                SDVectorAPCell2_ROI{j,i}=SDVectorAPCell_ROI{i};
                NParticlesAPCell2_ROI{j,i}=NParticlesAPCell_ROI{i};
            end

                MeanVectorAP_ROI=cell2mat(MeanVectorAPCell2_ROI);
                SDVectorAP_ROI=cell2mat(SDVectorAPCell2_ROI);
                NParticlesAP_ROI=cell2mat(NParticlesAPCell2_ROI);
        end

%       Get the corresponding mean information 
%      (nonROI, CompiledParticles_nonROI)
        k=1;
        for i=MinAPIndex:MaxAPIndex
            [MeanVectorAPTemp_nonROI,SDVectorAPTemp_nonROI,NParticlesAPTemp_nonROI]=AverageTracesNuclei(FrameInfo,...
                CompiledNuclei_nonROI(APFilter_nonROI(:,i)),NChannels);
            MeanVectorAPCell_nonROI{k}=MeanVectorAPTemp_nonROI';
            SDVectorAPCell_nonROI{k}=SDVectorAPTemp_nonROI';
            NParticlesAPCell_nonROI{k}=NParticlesAPTemp_nonROI';
            k=k+1;
        end

        %Turn the information into useful structures
        if NChannels>1
            for j=1:NChannels
                for i=MinAPIndex:MaxAPIndex
                    MeanVectorAPCell2_nonROI{j,i}=MeanVectorAPCell_nonROI{i}{j};
                    SDVectorAPCell2_nonROI{j,i}=SDVectorAPCell_nonROI{i}{j};
                    NParticlesAPCell2_nonROI{j,i}=NParticlesAPCell_nonROI{i}{j};
                end
            end

            for j=1:NChannels
                MeanVectorAP_nonROI{j}=cell2mat({MeanVectorAPCell2_nonROI{j,:}}')';
                SDVectorAP_nonROI{j}=cell2mat({SDVectorAPCell2_nonROI{j,:}}')';;
                NParticlesAP_nonROI{j}=cell2mat({NParticlesAPCell2_nonROI{j,:}}')';;
            end
        else
            for i=MinAPIndex:MaxAPIndex
                MeanVectorAPCell2_nonROI{j,i}=MeanVectorAPCell_nonROI{i};
                SDVectorAPCell2_nonROI{j,i}=SDVectorAPCell_nonROI{i};
                NParticlesAPCell2_nonROI{j,i}=NParticlesAPCell_nonROI{i};
            end

                MeanVectorAP_nonROI=cell2mat(MeanVectorAPCell2_nonROI);
                SDVectorAP_nonROI=cell2mat(SDVectorAPCell2_nonROI);
                NParticlesAP_nonROI=cell2mat(NParticlesAPCell2_nonROI);
        end
        
        %Get the corresponding mean information (for all nuclei, both ROI and non-ROI)
        k=1;
        for i=MinAPIndex:MaxAPIndex
            [MeanVectorAPTemp,SDVectorAPTemp,NParticlesAPTemp]=AverageTracesNuclei(FrameInfo,...
                CompiledNuclei(APFilter(:,i)),NChannels);
            MeanVectorAPCell{k}=MeanVectorAPTemp';
            SDVectorAPCell{k}=SDVectorAPTemp';
            NParticlesAPCell{k}=NParticlesAPTemp';
            k=k+1;
        end

        %Turn the information into useful structures
        if NChannels>1
            for j=1:NChannels
                for i=MinAPIndex:MaxAPIndex
                    MeanVectorAPCell2{j,i}=MeanVectorAPCell{i}{j};
                    SDVectorAPCell2{j,i}=SDVectorAPCell{i}{j};
                    NParticlesAPCell2{j,i}=NParticlesAPCell{i}{j};
                end
            end

            for j=1:NChannels
                MeanVectorAP{j}=cell2mat({MeanVectorAPCell2{j,:}}')';
                SDVectorAP{j}=cell2mat({SDVectorAPCell2{j,:}}')';;
                NParticlesAP{j}=cell2mat({NParticlesAPCell2{j,:}}')';;
            end
        else
            for i=MinAPIndex:MaxAPIndex
                MeanVectorAPCell2{j,i}=MeanVectorAPCell{i};
                SDVectorAPCell2{j,i}=SDVectorAPCell{i};
                NParticlesAPCell2{j,i}=NParticlesAPCell{i};
            end

            MeanVectorAP=cell2mat(MeanVectorAPCell2);
            SDVectorAP=cell2mat(SDVectorAPCell2);
            NParticlesAP=cell2mat(NParticlesAPCell2);
        end

    else % This is the case which we don't use ROI option

        %Get the corresponding mean information
        k=1;
        for i=MinAPIndex:MaxAPIndex
            [MeanVectorAPTemp,SDVectorAPTemp,NParticlesAPTemp]=AverageTracesNuclei(FrameInfo,...
                CompiledNuclei(APFilter(:,i)),NChannels);
            MeanVectorAPCell{k}=MeanVectorAPTemp';
            SDVectorAPCell{k}=SDVectorAPTemp';
            NParticlesAPCell{k}=NParticlesAPTemp';
            k=k+1;
        end

        %Turn the information into useful structures
        if NChannels>1
            for j=1:NChannels
                for i=MinAPIndex:MaxAPIndex
                    MeanVectorAPCell2{j,i}=MeanVectorAPCell{i}{j};
                    SDVectorAPCell2_nonROI{j,i}=SDVectorAPCell{i}{j};
                    NParticlesAPCell2_nonROI{j,i}=NParticlesAPCell{i}{j};
                end
            end

            for j=1:NChannels
                MeanVectorAP{j}=cell2mat({MeanVectorAPCell2{j,:}}')';
                SDVectorAP{j}=cell2mat({SDVectorAPCell2{j,:}}')';;
                NParticlesAP{j}=cell2mat({NParticlesAPCell2{j,:}}')';;
            end
        else
            for i=MinAPIndex:MaxAPIndex
                MeanVectorAPCell2{j,i}=MeanVectorAPCell{i};
                SDVectorAPCell2{j,i}=SDVectorAPCell{i};
                NParticlesAPCell2{j,i}=NParticlesAPCell{i};
            end

            MeanVectorAP=cell2mat(MeanVectorAPCell2);
            SDVectorAP=cell2mat(SDVectorAPCell2);
            NParticlesAP=cell2mat(NParticlesAPCell2);
        end
    end
end
    

%Calculate the mean for all of them
[MeanVectorAll,SDVectorAll,NParticlesAll]=AverageTracesNuclei(FrameInfo,CompiledNuclei);
%Now find the different maxima in each nc

MaxFrame=[];
for i=1:length(NewCyclePos)
    if i==1
        [~,MaxIndex]=max(MeanVectorAll(1:NewCyclePos(1)));
        MaxFrame=[MaxFrame,MaxIndex];
    elseif i<=length(NewCyclePos)
        [~,MaxIndex]=max(MeanVectorAll(NewCyclePos(i-1):NewCyclePos(i)));
        MaxFrame=[MaxFrame,NewCyclePos(i-1)+MaxIndex-1];
    end
end
[~,MaxIndex]=max(MeanVectorAll(NewCyclePos(i):end));
MaxFrame=[MaxFrame,NewCyclePos(i)+MaxIndex-1];





%% Save everything


if strcmpi(ExperimentAxis,'ap')
    
    % consider ROI option
    if ROI
        save([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'],...
            'CompiledNuclei','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
            'nc12','nc13','nc14','ncFilterID','ncFilter','APbinID','APFilter',...
            'MeanVectorAP','SDVectorAP','NParticlesAP','MeanVectorAll',...
            'SDVectorAll','NParticlesAll','MaxFrame','MinAPIndex','MaxAPIndex',...
            'AllTracesVector','AllTracesAP',...
            'MeanCyto','SDCyto','MedianCyto','MaxCyto',...
            'MeanCytoAPProfile','SDCytoAPProfile','SECytoAPProfile',...
            'MeanVectorAP_ROI','SDVectorAP_ROI','NParticlesAP_ROI',...
            'MeanVectorAP_nonROI','SDVectorAP_nonROI','NParticlesAP_nonROI',...
            'CompiledNuclei_ROI','CompiledNuclei_nonROI', '-v7.3')%,...
            %'IntegrationArea')
    else 
        save([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'],...
            'CompiledNuclei','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
            'nc12','nc13','nc14','ncFilterID','ncFilter','APbinID','APFilter',...
            'MeanVectorAP','SDVectorAP','NParticlesAP','MeanVectorAll',...
            'SDVectorAll','NParticlesAll','MaxFrame','MinAPIndex','MaxAPIndex',...
            'AllTracesVector','AllTracesAP',...
            'MeanCyto','SDCyto','MedianCyto','MaxCyto',...
            'MeanCytoAPProfile','SDCytoAPProfile','SECytoAPProfile', '-v7.3')%,...
            %'IntegrationArea')
    end
else
    save([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'],...
        'CompiledNuclei','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
        'nc12','nc13','nc14','ncFilterID','ncFilter','MeanVectorAll',...
        'SDVectorAll','NParticlesAll','MaxFrame',...
        'MeanCyto','SDCyto','MedianCyto','MaxCyto','IntegrationArea', '-v7.3')
end


save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells', '-v7.3')


