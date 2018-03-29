function AccumulationMovie(Prefix)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initial Setup to get folder names etc %%%%%%%%%%%%%%%%%%%%%

[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, SchnitzcellsFolder,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);

[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)


ncs = [nc9,nc10,nc11,nc12,nc13,nc14];

D = dir([FISHPath,filesep,'Data',filesep,Prefix,filesep,'*His*.tif']);


%%%%%%%%%%% Part 1: Take max nuclei and save in single structure in Data folder %%%%%%%%%%%%%%

MaxNuclei=struct;

TotalTime=length(D);

for i=1:TotalTime

    if i<10
        strii=['00', num2str(i)];
    elseif i<100
        strii=['0', num2str(i)];
    else
        strii=[num2str(i)];
    end

MaxNuclei.(['Time', num2str(i)])= imread([FISHPath,filesep,'Data',filesep,Prefix,filesep,D(i).name]);

end

save([FISHPath,filesep,'Data',filesep,Prefix,filesep,'MaxNuclei.mat'],'MaxNuclei');

%%%%%%%%% Part 2: Perform diff gaussian filtering to segment nuclei 

% Figure out pixel size  

    HyphenPositionR = find(Prefix == '-');
    DateFolderS = Prefix(1 : HyphenPositionR(3)-1);
    LineFolderS = Prefix(HyphenPositionR(3)+1 : end);
    Folder = [SourcePath, filesep, DateFolderS, filesep, LineFolderS];


DTIF=dir([Folder,filesep,'*.tif']);
DLSM=dir([Folder,filesep,'*.lsm']);

if (length(DTIF)>0)&(length(DLSM)==0)
    
        display('2-photon @ Princeton data mode')
        DRaw=DTIF;
        FileMode='TIF';
        
        info = imfinfo([Folder,filesep,DRaw(1).name]);
        
%PixelWidth = info(1).XResolution*2.54e-2;

PixelWidth=200e-9 % Hack !
%%%%%% not reading in directly from header, Xresolution field seems to give same value for all magnifications?
        
        
elseif (length(DTIF)==0)&(length(DLSM)>0)
    
        display('LSM mode')
        DRaw=DLSM;
        FileMode='LSM';
        
load([Folder,filesep,DRaw(1).name(1:end-4)]);

PixelWidth=Datas.LSM_info.VoxelSizeX

    
else
    PixelWidth=200e-9 % Hack !
    %error('File type not recognized')
end
% 
% %%%%%%%%%%%%%%%%%% Use nuclei size and cc times to perform first pass
% %%%%%%%%%%%%%%%%%% filtering segmentation and rough tracking
% 
% RadiidiffCell=round(PixelWidth/2e-7*[25,21,18,15,11,7]);
% 
% OptimalRadius=ones(1,TotalTime);
% 
% FirstDivis = find(ncs,1);
% 
% ncsMOD = [ncs, TotalTime];
% 
% OptimalRadius(1:ncs(FirstDivis)) = RadiidiffCell(FirstDivis-1);
% 
% for i=FirstDivis:6
%     
%     OptimalRadius(ncsMOD(i)+1:ncsMOD(i+1)) = RadiidiffCell(i);
% end    
% 
% LabelNucsCore = SegmentNucleiLive([FISHPath,filesep,'Data',filesep,Prefix,filesep],[],OptimalRadius,round(6*OptimalRadius));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Dilate Nuclei
% 
% DilateNucleiLive([FISHPath,filesep,'Data',filesep,Prefix,filesep])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Make structure with intensity values
 
    DataEve=load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);
    
    Schnitzcells=load([DropboxFolder,filesep,Prefix(1:end),filesep,Prefix(1:end),'_lin.mat']);
    
    %Load Ellipses
    Ellipses=load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat']);

%%%%%%%%%%%%%%%%%%%  Now add the extra information

    DataEve.schnitzcells=Schnitzcells.schnitzcells;
    DataEve.Ellipses=Ellipses.Ellipses;


NN=length(DataEve.ElapsedTime);

Blk = zeros(length(DataEve.CompiledParticles),NN);

PosiX=zeros(length(DataEve.CompiledParticles),1);
PosiY=zeros(length(DataEve.CompiledParticles),1);

jj=1;

    for j=1:length(DataEve.CompiledParticles)
            
            Blk(j,DataEve.CompiledParticles(j).Frame)=DataEve.CompiledParticles(j).Fluo;
       
            PosiX(j) = mean(DataEve.CompiledParticles(j).xPos);
            PosiY(j) = mean(DataEve.CompiledParticles(j).yPos);

    end
%%%%%%%%%%%%%

load([FISHPath,filesep,'Data',filesep,Prefix,filesep,'LabNucDilate.mat']);
load([FISHPath,filesep,'Data',filesep,Prefix,filesep,'LabelNucsCore.mat']);

%%%%%%% First connect the HG Dots to JB nuclei %%%%%%

NucSOFDots = [DataEve(1).CompiledParticles.Nucleus]'; % The nuclei of the different dots

for i=1:length(NucSOFDots)

    frames = DataEve(1).schnitzcells(NucSOFDots(i)).frames;
    cenx = DataEve(1).schnitzcells(NucSOFDots(i)).cenx;
    ceny = DataEve(1).schnitzcells(NucSOFDots(i)).ceny;

    HGDotstoJBNucs(i).JBNucs=[];
    HGDotstoJBNucs(i).HGNuc=NucSOFDots(i);

    for jj=1:length(frames)
        
        JacquesCentroids = cat(1,LabelNucsCore.(['Time',num2str(frames(jj))]).RegionProps.Centroid);
        
        [neighborIds neighborDistances] = kNearestNeighbors(JacquesCentroids, [cenx(jj), ceny(jj)], 1);
        
        JacquesLinearIndex1=LabelNucsCore.(['Time',num2str(frames(jj))]).RegionProps(neighborIds).PixelIdxList;
        
        JacquesIndex=unique(LabelNucsCore.(['Time',num2str(frames(jj))]).ImageZ(JacquesLinearIndex1));
        
        HGDotstoJBNucs(i).JBNucs=[HGDotstoJBNucs(i).JBNucs,JacquesIndex];
        
    end
    
end


%%%%%%% Estimate the amount of cytoplasmic mRNA present ( need to think
%%%%%%% more about this)


 halflife=10;
 ElongationTime =4;
 
 Time=DataEve.ElapsedTime;
 
 Blkk=[];
%  options = odeset('NonNegative',1);
%  for NBlk=1:length(NucSOFDots);
%       [TOUT,YOUT] = ode45(@(t,y)NumericIntegrationmRNA(t,y,halflife,ElongationTime,Blk(NBlk,:),Time),Time,0,options);% Still need to apply shift
%       Blkk=[Blkk;YOUT'];
%  end
% 
% 
% halflife=10;

Blkk=cumtrapz(Time,Blk,2);

%%%%%%%%%%%%%%%%%%%%%%% Make movie %%%%%%%%%%%%%%%%%

maxx=max(Blkk(:));

N=10000; % Number of discrete steps

  writerObj = VideoWriter(['AccumulationMovieInfHalflife.avi']);%,'Uncompressed AVI');
        writerObj.FrameRate = 7;
        open(writerObj);
        
DisplayBAll=struct;
DisplayBOn=struct;
        
for TT=1:TotalTime;
     
     DisplayB=zeros(size(LabNucDilate.Time1.Image));
     DisplayBOnT=zeros(size(LabNucDilate.Time1.Image));
 
 for i=1:length(NucSOFDots)
     
     Ind=find(ismember(LabNucDilate.(['Time', num2str(TT)]).Image,HGDotstoJBNucs(i).JBNucs));
     
     try
        
          DisplayB(Ind)=Blkk(i,TT);
          DisplayBOnT(Ind)=Blk(i,TT);
          
        
     catch
     end
 
 end
 
DisplayBZ=round(DisplayB/maxx*N);

DisplayBZ(1,1)=N;

DisplayBAll.(['Time', num2str(TT)]) = DisplayBZ+1;
DisplayBOn.(['Time', num2str(TT)]) = DisplayBOnT;

imshow(label2rgbBackdropLive(DisplayBZ+1,'jet',[0,0,0],imadjust(im2double(MaxNuclei.(['Time', num2str(TT)])))));

%pause(0.1)
           frame = getframe;
        
        writeVideo(writerObj,frame);
        
        clf
    
end

    close(writerObj);

save([FISHPath,filesep,'Data',filesep,Prefix,filesep,'DisplayBAll.mat'],'DisplayBAll','DisplayBOn');    
   % 
save([DropboxFolder,filesep,Prefix,filesep,'AccumulationData.mat'],'Blk','Blkk','HGDotstoJBNucs');

%    save_to_base(1)

