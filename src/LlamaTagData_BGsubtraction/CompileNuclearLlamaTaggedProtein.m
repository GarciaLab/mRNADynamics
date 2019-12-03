function CompileNuclearLlamaTaggedProtein(varargin)
%% DESCRIPTION
% Inputs : schnitzcells (Background subtraction done by
% processLlamaTagData_subtractBGFluo.m)
% CompiledNuclei (so that we can update this).

% Outputs : (1) Mean and STD of nuclear fluo (BG subtracted), saved in
% CompiledNuclei.mat

%% Define the folder directory

[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

Prefix=varargin{1};

FilePrefix=[Prefix,'_'];

%% Now get the actual Dropbox folder
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath,...
    configValues, movieDatabasePath, movieDatabaseFolder, movieDatabase]=...
    DetermineLocalFolders(Prefix);

% This is to get the name of the channels from MovieDatabase
[XLSNum,XLSTxt,XLSRaw]=xlsread([DropboxFolder,filesep,'MovieDatabase.csv']);


%Check that FrameInfo exists
if exist([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
end

numFrames = length(FrameInfo);

%See if FrameInfo has information about the number of input channels. This
%is not fully implemented yet. If no information is found, then assume we
%have only one input channel.
if isfield(FrameInfo,'NChInput')
    NChannels=FrameInfo(1).NChInput;
else
    NChannels=1;
end

[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, movieDatabase)

% Load files that needs to be further processed.
% CompiledNuclei
load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'])
% schnitzcells
load([DropboxFolder,filesep,Prefix,filesep,[Prefix '_lin.mat']])
%load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'])

%Check that FrameInfo exists
if exist([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
end

%% Initialize some variables
% NameString_ROI='';
%% Convert the addtional info in schnitzcells to CompiledNuclei
h=waitbar(0,'Compiling nuclear traces');
k=1;
for i=1:length(schnitzcells)
    waitbar(i/length(schnitzcells),h)
    
    if schnitzcells(i).Approved
        %Which frames were approved manually?
        FrameFilter=schnitzcells(i).FrameApproved;
        
        if sum(FrameFilter)
            
            %Copy the filtered information
            CompiledNuclei(k).FluoMax_BGsubtracted=schnitzcells(i).Fluo_BGsubtracted;
            %Save the information about the original schnitz
            % CompiledNuclei(k).schnitz=uint16(i);
            % Indexing
            k=k+1;
        end
    end
end
close(h)     

if strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')
    %AP filters:

    %Divide the AP axis into boxes of a certain AP size. We'll see which
    %particle falls where.

    APbinID=0:APResolution:1;

    APFilter=false(length(CompiledNuclei),length(APbinID));
    
    for i=1:length(CompiledNuclei)
        APFilter(i,max(find(APbinID<=CompiledNuclei(i).MeanAP)))=true;
    end
else
    APbinID = [];
    APFilter = [];
end
%% Binning and averaging data
% Main change here from our current pipeline is the AverageTracesNuclei.m
% Now this AverageTracesNuclei has additional input, "BGsubtracted", in
% which case if it's 1, then it will select the "FluoMax_BGsubtracted" in
% favor of "FluoMax"

% If background is subtracted?
BGsubtracted = true;

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
    k= MinAPIndex;
    for ap=MinAPIndex:MaxAPIndex
        [MeanVectorAPTemp,SDVectorAPTemp,NParticlesAPTemp]=AverageTracesNuclei(FrameInfo,...
            CompiledNuclei(APFilter(:,ap)),NChannels,BGsubtracted);
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
            MeanVectorAPCell2{ap}=double(MeanVectorAPCell{ap});
            SDVectorAPCell2{ap}=double(SDVectorAPCell{ap});
            NParticlesAPCell2{ap}=double(NParticlesAPCell{ap});
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
    
% Redefine the Mean, STD, and N for BG_subtracted nuclei
MeanVectorAP_BGsubtracted = MeanVectorAP;
SDVectorAP_BGsubtracted = SDVectorAP;
NParticlesAP_BGsubtracted = NParticlesAP;

% %Calculate the mean for all of them
% [MeanVectorAll,SDVectorAll,NParticlesAll]=AverageTracesNuclei(FrameInfo,CompiledNuclei);

%% Save additional fields to CompiledNuclei.mat

save([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'],...
       'MeanVectorAP_BGsubtracted', 'SDVectorAP_BGsubtracted', 'NParticlesAP_BGsubtracted',...
       '-append');

% save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells', '-v7.3')
end