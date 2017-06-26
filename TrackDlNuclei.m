function TrackDlNuclei (varargin)
%This function creates a new His channel using the Venus channel instead of
%the mCherry channel.
%Pending: Make the function automatically save the former His channel in a
%folder called 'OriginalHis' within the PreProcessedData/Prefix folder

%Information about about folders
[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

%Look at the input parameters and use defaults if missing
%Prefix=[];
ForceAP=0;      %Force AP detection even if it's already there
SkipTraces=0;   %Do not output the individual traces.
SkipFluctuations=0;  %Do not generate the plots of correlations of fluctuations and offset
SkipFits=0;         %Do not generate the fit output (but still does the fit)
SkipMovie=0;        %Do not generate the movie
ApproveAll=0;       %Only use manually approved particles
MinParticles=4;     %Require 4 particles per AP bin or else discard

if isempty(varargin)
    FolderTemp=uigetdir(DefaultDropboxFolder,'Select folder with data to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
else
    Prefix=varargin{1};

    for i=2:length(varargin)
        if strcmp(varargin{i},'InvertChannel')
            display('The inverted dl channel will be used')
            Invert=1;
%         elseif strcmp(varargin{i},'SkipTraces')
%             SkipTraces=1;
%         elseif strcmp(varargin{i},'SkipFluctuations')
%             SkipFluctuations=1;
%         elseif strcmp(varargin{i},'SkipFits')    
%             SkipFits=1;
%         elseif strcmp(varargin{i},'SkipMovie')    
%             SkipMovie=1;
%         elseif strcmp(varargin{i},'SkipAll')        
%             SkipTraces=1;
%             SkipFluctuations=1;
%             SkipFits=1;
%             SkipMovie=1;
%         elseif strcmp(varargin{i},'ApproveAll')    
%             ApproveAll=1;
%         elseif strcmp(varargin{i},'SetMinParticles')
%             MinParticles = input('Set minimum particle threshold:');
        end
    end

end
FilePrefix=[Prefix,'_'];

%Now get the actual Dropbox folder
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);

%What type of experiment are we dealing with? Get this out of
%MovieDatabase.xlsx

[XLSNum,XLSTxt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
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

%//////////////     Load all the information      ///////////////////
%load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
%load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'])
%load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'])
%load([DropboxFolder,filesep,Prefix,filesep,[Prefix '_lin.mat']])

%Check that FrameInfo exists
if exist([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
else
    warning('No FrameInfo.mat found. Trying to continue')
    %Adding frame information
    DHis=dir([PreProcPath,filesep,FilePrefix(1:end-1),filesep,'*His*.tif']);
    FrameInfo(length(DHis)).nc=[];
    %Adding information

    Dz=dir([PreProcPath,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'*001*.tif']);
    NumberSlices=length(Dz)-1;
    
    for i=1:length(FrameInfo)
        FrameInfo(i).NumberSlices=NumberSlices;
    end
    
end

%////////////// END OF TrackmRNA  SCRIPT LINES    \\\\\\\\\\\\\\\\\

%Go to the fake His folder
%open folder. Channel 01 is -usually- the Venus channel. How to automate
%this
D=dir([PreProcPath,filesep,Prefix,filesep,'*ch01.tif']);
H=dir([PreProcPath,filesep,Prefix,filesep,'*His*.tif']);



%we need to know to what frame the histone channel corresponds to
for i=1:length(H)
    Name = H(i).name; %get image name
    SplitName = strsplit(Name,'_');
    FrameName = SplitName{end};
    Frame = str2num(FrameName(1:end-3));
    H(i).Frame = Frame;
end

%we need to know to what frame does it image corresponds to.
%Create field showing to what frame and z position the image belongs to
for i = 1:length(D)
    Name = D(i).name; %get image name
    SplitName = strsplit(Name,'_'); %split name
    FrameName = SplitName{end-2}; %get the frame value in the name..PROBLEM HERE!! how is '001' different to '01'?
    ZposName = SplitName{end-1};
    ChannelName = SplitName{end};
    %Frame = FrameName(2:end); %get rid of the t in the frame value
    Zpos = ZposName(2:end);
    Channel = ChannelName(4);
    D(i).Frame = str2num(FrameName); %assign frame to the current image
    D(i).Zpos = str2num(Zpos);
    D(i).Channel = Channel;
end

%////////////   Z project the Dl-Channel    /////////////
%Create a 4D array where the first two dimensions are the image
%size, the third one is the frame and the fourth one is Z
%TotalFrames = max([D.Frame])
h = waitbar(0,'Please wait...stacking Dl-Venus channel');
%D is the Dl-Venus channel
for frm=1:max([D.Frame])
   for img = 1:length(D)
       if D(img).Frame == frm & D(img).Channel == '1'
           Image = imread(D(img).name);
           if Invert
               Image = imcomplement(Image);
           end
           FrameZ(:,:,frm,D(img).Zpos)=Image;
       end
   end
   waitbar(frm / max([D.Frame]));
end
close(h)

% Re-write histone channel
DlFrames = [1:length(D)];
Zstacks = 23;
MidZ = ceil(Zstacks/2);
MidZ = [MidZ-2,MidZ+2];
for fr = 1:(length(D)/Zstacks);
    frame = num2str(fr);
    Venus = mean(FrameZ(:,:,fr,(MidZ(1):MidZ(2))),4);
    NormV = Venus - min(Venus(:));
    NormV = NormV ./ max(NormV(:));
    Name = H(fr).name;
    Name = ([PreProcPath,filesep,Prefix,filesep,Name]);
    imwrite(NormV,Name);
end








% %%
% %Or we can try to have an automated way to do it
% %Get the relevant folders now:
% [SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
%     DetermineLocalFolders(Prefix);
% %Figure out what type of experiment we have
% [XLSNum,XLSTxt]=xlsread([DropboxFolder,filesep,'MovieDatabase.xlsx']);
% 
% nc10Column=find(strcmp(XLSTxt(1,:),'nc10'));
% nc11Column=find(strcmp(XLSTxt(1,:),'nc11'));
% nc12Column=find(strcmp(XLSTxt(1,:),'nc12'));
% nc13Column=find(strcmp(XLSTxt(1,:),'nc13'));
% nc14Column=find(strcmp(XLSTxt(1,:),'nc14'));
% DataFolderColumn=find(strcmp(XLSTxt(1,:),'DataFolder'));
% 
% % Convert the prefix into the string used in the XLS file
% Dashes = strfind(Prefix, '-');
% PrefixRow = find(strcmp(XLSTxt(:, DataFolderColumn),...
%     [Prefix(1:Dashes(3)-1), '\', Prefix(Dashes(3)+1:end)]));
% if isempty(PrefixRow)
%     PrefixRow = find(strcmp(XLSTxt(:, DataFolderColumn),...
%         [Prefix(1:Dashes(3)-1), '/', Prefix(Dashes(3)+1:end)]));
%     if isempty(PrefixRow)
%         error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
%     end
% end
% 
% nc10 = XLSNum(PrefixRow,10);
% nc11 = XLSNum(PrefixRow,11);
% nc12 = XLSNum(PrefixRow,12);
% nc13 = XLSNum(PrefixRow,13);
% nc14 = XLSNum(PrefixRow,14);
% 
% 
% NCs = [nc10 nc11 nc12 nc13 nc14]
% NCs = NCs(NCs>0);
% MitosisFrames=[];
% for i=1:length(NCs)
%     MitosisFrames=[MitosisFrames NCs(i)-13:NCs(i)+16]
% end
% MitosisFrames(MitosisFrames<=0)=1;
% %%
%Show results
% for fr = ManualMitFrames2
%     frame = num2str(fr)%for labels
%     
%     subplot(3,1,1)
%     Venus = mean(FrameZ(:,:,fr,(8:12)),4);
%     imshow(Venus,[]);
%     title(['Dorsal-Venus channel , Frame' frame])
% 
%     
%     subplot(3,1,2)
%     Name = H(fr).name;
%     Hist = imread(Name);
%     imshow(Hist,[])
%     title(['Fake Histone (mCherry) channel , Frame' frame])
%     
%     subplot(3,1,3)
%     Venus2 = (uint8(Venus))*30; %multiply to match His intensity range
%     VenusHist(:,:,2)=Hist;    
%     VenusHist(:,:,1)=Venus2;
%     MaxVenusHist = max(VenusHist,[],3); 
%     imshow(MaxVenusHist,[])
%     title(['Combined channel, Frame:' frame])
%     
%     imwrite(MaxVenusHist,Name)
%     pause(0.1)
% end
% 
% 
% % 
% % %%
% % %Visualize our new His channel
% % for fr = 1:length(H)
% %     frame = num2str(fr)%for labels
% %     Name = H(fr).name;
% %     Hist = imread(Name);
% %     subplot(1,1,1)
% %     imshow(Hist,[]);
% %     title(['Combined channel, Frame:' frame])
% %     pause(0.1)
% %     
% % end

