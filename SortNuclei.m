function SortNuclei (varargin)
% This function takes the fluorescence of each nucleus at a defined time
% after mitosis. It uses this information to assign a 'Fiduciary
% Fluorescence' value to each nucleus

%% Get Folder info, set up stuff, etc

[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

Prefix=varargin{1};
 
FilePrefix=[Prefix,'_'];

%Now get the actual Dropbox folder
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);

%Load all the information

%load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'])
%load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'])
load([DropboxFolder,filesep,Prefix,filesep,[Prefix '_lin.mat']])

%Check that FrameInfo exists
if exist([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
end

%get frames of each mitosis
[XLSNum,XLSTxt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
nc10Column=find(strcmp(XLSRaw(1,:),'nc10'));
nc11Column=find(strcmp(XLSRaw(1,:),'nc11'));
nc12Column= strcmp(XLSRaw(1,:),'nc12');
nc13Column=find(strcmp(XLSRaw(1,:),'nc13'));
nc14Column=find(strcmp(XLSRaw(1,:),'nc14'));

DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));

% get the mitotic frames from the MovieDatabase excel
nc10=XLSRaw{PrefixRow,nc10Column};
nc11=XLSRaw{PrefixRow,nc11Column};
nc12=XLSRaw{PrefixRow,nc12Column};
nc13=XLSRaw{PrefixRow,nc13Column};
nc14=XLSRaw{PrefixRow,nc14Column};

%%
if exist([DropboxFolder,filesep,Prefix,filesep,'SortedNuclei.mat'])||0
    display('Re-running SortNuclei, former results will be deleted')
    cd ([DropboxFolder,filesep,Prefix])
    delete 'SortedNuclei.mat'
end


%% Get Frame info about mitosis

NcFrames = [nc10,nc11,nc12,nc13,nc14]; %Store the frame where each mitosis occurs in a single array

%Get the absolute time of each mitosis frame
for NCs = 1:length(NcFrames)
    NcFrame = NcFrames(NCs); %get the frame where this nc starts
    if NcFrame == 0 %if the frame was 0 then that mitosis wasn't in the movie
        NcTimes(NCs) = 0;
    else
    AbsTime = FrameInfo(NcFrame).Time; %FrameInfo contains the absolute time for each frame
    NcTimes(NCs) = AbsTime;
    end
end

Ncsi = [10,11,12,13,14]; 
MitoticFrames =[]; %to store the frames where mitosis occur
MitoticTimes=[]; %to store the absolute time where mitosis occur

%MitoticFrames and MitoticTimes matrices:
%first colum contains the frame (or time) where the nc starts, second column is the nc number,

for nc = 1:length(NcFrames)
    ncName = Ncsi(nc);
    MitoticFrames(nc,1) = NcFrames(nc); %start frame
    MitoticFrames(nc,2) = ncName; %identity of nuclear cycle
    MitoticTimes(nc,1) = NcTimes(nc);
    MitoticTimes(nc,2) =  ncName;
end

%% Use Absolute time to define fiduciary fluorescence
% We assign each nucleus a fiduciary fluorescence value based on the value
% at a given constant time after mitosis. We take 'IntWindow' frames before and the same number of
% frames after. Then we take the mean.

schnitzIndex = 1;
IntWindow = 6; %number of frames before and after the 'fiduciary' frame
FiduciaryTime = [160,160,340,650,500]; %Seconds after mitosis, the moment when we meassure the gradient
EffectiveNCs = [max(find(MitoticFrames(:,1)==0)):5]; %find out which NCs we have in the movie

for s = 1:length(schnitzcells) %go over all schnitz
    for nc = EffectiveNCs   %all sorting is done nc-wise. This can handle from nc10-nc14  
        %info about the sorting frame
        Mitosis = MitoticTimes(nc,1); %get the absolute time at which mitosis occurs
        SortTime = Mitosis + FiduciaryTime(nc); %Absolute fiduciary time in seconds
        %find the frame which is the closest to the sorting frame
        ClosestFrame = abs([FrameInfo.Time]-SortTime);
        [dummy SortFrame] = min(ClosestFrame);
        
        %info about nucleus frames
        schnitzFrames = schnitzcells(s).frames;
        StartFrame = FrameInfo(schnitzFrames(1)).Time;
        
        if StartFrame > Mitosis && any(schnitzFrames == SortFrame)  %if nucleus belong to this nc AND is present during arbitrary window
 
            SortFramePos = find(schnitzFrames==SortFrame); %find position of the sorting frame in the frames vector
            schnitzFluo = NucleusConcentration2(schnitzcells(s).Fluo);
            FiduciaryFluo = schnitzFluo(SortFramePos); %get fluorescence value at the sorting frame
            FiduciaryFrames = SortFrame;
           
%           Now take fluorescence around the sorting frame
            for i = 1:IntWindow %number of frames used to integrate fluorescence around the sorting frame                    
                if SortFramePos-i > 0 && SortFramePos+i <= length(schnitzFrames) %are there enough frames around?
                    FiduciaryFluo = [schnitzFluo(SortFramePos-i) FiduciaryFluo schnitzFluo(SortFramePos+i)]; 
                    FiduciaryFrames = [schnitzFrames(SortFramePos-i) FiduciaryFrames schnitzFrames(SortFramePos+i)];
                end
            end
            
            %if the schnitz belonged to the nc and was present at the
            %sorting frame, save it in the 'SortedNuclei' struct
            if ~isnan(nanmean(FiduciaryFluo))
                SortedNuclei(schnitzIndex).SortFrameTime = FrameInfo(SortFrame).Time;
                SortedNuclei(schnitzIndex).nc = nc+9;
                SortedNuclei(schnitzIndex).OriginalSchnitz = s;
                if s==5
                    '5 found'
                    
                end
                SortedNuclei(schnitzIndex).Fluo = schnitzFluo;
                SortedNuclei(schnitzIndex).Frames = schnitzcells(s).frames;
                SortedNuclei(schnitzIndex).FiduciaryFluoT = FiduciaryFluo;
                SortedNuclei(schnitzIndex).FiduciaryFramesT = FiduciaryFrames;
                SortedNuclei(schnitzIndex).MeanFidFluoT = nanmean(FiduciaryFluo);
                SortedNuclei(schnitzIndex).Prefix = Prefix;
                
                SortedNuclei(schnitzIndex).NucleusTime = [];
                MovieFrames = [schnitzcells(s).frames];
                for frm = 1:length(MovieFrames)
                    MovieFrame = MovieFrames(frm);
                    Time = FrameInfo(MovieFrame).Time;
                    SortedNuclei(schnitzIndex).NucleusTime = [SortedNuclei(schnitzIndex).NucleusTime Time];
                end
                schnitzIndex = schnitzIndex+1;
            end
        end
        clear PeakSchnitzFluo
    end
end

%%

%save([DropboxFolder,filesep,Prefix,filesep,'SortedNuclei.mat'],'SortedNuclei')


% /////////////     Add Particle info to SortedNuclei struct     \\\\\\\\\\\\\\\\\

%first, figure out which nuclei are on and off.

%a priori all are OFF
for i = 1:length(SortedNuclei) 
    SortedNuclei(i).On = 0;
end

OnNuclei = [CompiledParticles.Nucleus]; %get the ID of nuclei that contain a particle

%Now mark ON nuclei
for s = 1:length(SortedNuclei)
    OrS = SortedNuclei(s).OriginalSchnitz;
    if any(OnNuclei==OrS)
        SortedNuclei(s).On = 1;
    end
end

%Now let's add particle info to ON nuclei
for s = 1:length(SortedNuclei)
    if SortedNuclei(s).On
        OrS = SortedNuclei(s).OriginalSchnitz;
        Part = find(OnNuclei==OrS);
        SortedNuclei(s).CompiledParticle = Part;
        SortedNuclei(s).mRNA = CompiledParticles(Part).Fluo;
        SortedNuclei(s).mRNAError = CompiledParticles(Part).FluoError;
        SortedNuclei(s).mRNAFrames = CompiledParticles(Part).Frame;
    end
end


%% Add time after mitosis to nuclei and particles
for s = 1:length(SortedNuclei)   
    NC = SortedNuclei(s).nc;
    MitoFrame = MitoticFrames(NC-9,1);
    MitoTime = MitoticTimes(NC-9,1);
    
    NucleusFrames = SortedNuclei(s).Frames;
    for nf = 1:length(NucleusFrames)
        frame = NucleusFrames(nf);
        NucleusTimes(nf) = FrameInfo(frame).Time;
    end
    NucleusT = NucleusTimes - MitoTime;
    SortedNuclei(s).NucleusTimes = NucleusT;
end
%%
for s=1:length(SortedNuclei)
    NC = SortedNuclei(s).nc;
    MitoTime = MitoticTimes(NC-9);
    ParticleFrames = SortedNuclei(s).mRNAFrames;
    
    if ~isempty(ParticleFrames)
        for pf = 1:length(ParticleFrames)
            frame=ParticleFrames(pf);
            ParticleTimes(pf) = FrameInfo(frame).Time;
        end
    else
        ParticleTimes=[];
    end    
    try
    length(ParticleTimes)==length(ParticleFrames)
    catch
       %ERR
    end
    ParticleT = ParticleTimes-MitoTime;
    SortedNuclei(s).ParticleTimes = ParticleT;  
    clear ParticleTimes

end


%% Save the data

save([DropboxFolder,filesep,Prefix,filesep,'SortedNuclei.mat'],'SortedNuclei')
%save([DropboxFolder,filesep,Prefix,filesep,'PeakFrames2'],'PeakFrames')

%% Show output

h(1) = figure
hold on
palette = ['k','k','k','k','k','k','k','k','k','y','m','r','g','b'];
for n = 1:length(SortedNuclei)
    nc = SortedNuclei(n).nc;
    color = palette(nc);
    oriS = SortedNuclei(n).OriginalSchnitz;
    XPos = mean(schnitzcells(oriS).cenx);
    Fluo = SortedNuclei(n).MeanFidFluoT;
    if SortedNuclei(n).On
        plot(XPos,Fluo,'o','Color',color,'MarkerFaceColor',color)
    else
        plot(XPos,Fluo,'o','Color',color)
    end
end
hold off
ylabel('Input fluorescence at sorting frame')
xlabel('Nucleus X position at sorting frame')
%title(strcat('Nuclei lost:',num2str(LostNuclei)))
savefig(h,'SortedNuclei','compact')


