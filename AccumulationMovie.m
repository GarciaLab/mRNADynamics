function AccumulationMovie(Prefix)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initial Setup to get folder names etc %%%%%%%%%%%%%%%%%%%%%

%Get the folders, including the default Dropbox one
[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders;
%Now get the actual DropboxFolder
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(Prefix);



%Determine division times
%Load the information about the nc from the XLS file
[Num,Txt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
XLSHeaders=Txt(1,:);
Txt=Txt(2:end,:);

ExperimentTypeColumn=find(strcmp(XLSRaw(1,:),'ExperimentType'));
ExperimentAxisColumn=find(strcmp(XLSRaw(1,:),'ExperimentAxis'));

DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));

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
Channel2Column=find(strcmp(XLSRaw(1,:),'Channel2'));


%Find the corresponding entry in the XLS file
if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))
    warning('This step in CheckParticleTracking will most likely have to be modified to work')
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [Date,'\BcdGFP-HisRFP']));
else
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
    
    if isempty(XLSEntry)
    disp('%%%%%%%%%%%%%%%%%%%%%')
    disp('Dateset could not be found. Check MovieDatabase.xlsx')
    disp('%%%%%%%%%%%%%%%%%%%%%')
    end
end


if strcmp(XLSRaw(XLSEntry,Channel2Column),'His-RFP')
    nc9=XLSRaw{XLSEntry,nc9Column};
    nc10=XLSRaw{XLSEntry,nc10Column};
    nc11=XLSRaw{XLSEntry,nc11Column};
    nc12=XLSRaw{XLSEntry,nc12Column};
    nc13=XLSRaw{XLSEntry,nc13Column};
    nc14=XLSRaw{XLSEntry,nc14Column};
    CF=XLSRaw{XLSEntry,CFColumn};
end
    

%Convert the prefix into the string used in the XLS file
Dashes=findstr(Prefix,'-');

%Find the corresponding entry in the XLS file
if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))
    XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
        [Date,'\BcdGFP-HisRFP']));
else
    XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
        [Prefix(1:Dashes(3)-1),filesep,Prefix(Dashes(3)+1:end)]));
end
ncs=[nc9,nc10,nc11,nc12,nc13,nc14];

D=dir([FISHPath,filesep,'Data',filesep,Prefix,filesep,'*His*.tif']);


%%%%%%%%%%% Part 1: Take max nuclei and save in single structure in Data
%%%%%%%%%%% folder

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

PixelWidth=220e-9 % Hack !
%%%%%% not reading in directly from header, Xresolution field seems to give samevalue for all magnifications?
        
        
elseif (length(DTIF)==0)&(length(DLSM)>0)
    
        display('LSM mode')
        DRaw=DLSM;
        FileMode='LSM';
        
load([Folder,filesep,DRaw(1).name(1:end-4)]);

PixelWidth=Datas.LSM_info.VoxelSizeX

    
else
    error('File type not recognized')
end

%%%%%%%%%%%%%%%%%% Use nuclei size and cc times to perform first pass filtering

save_to_base(1)

% OptimalRadius=[];
% 
% for i=1:TotalTime
% OptimalRadius(i) = DetermineFilterRadius(MaxNuclei.(['Time', num2str(i)]),3,15,50,0.5,20,2);
% end
% 
% OptimalRadius(isnan(OptimalRadius))=10;



