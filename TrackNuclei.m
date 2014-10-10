function TrackNuclei(Prefix)

%This function is just a script that call Laurent's tracking code

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
Channel1Column=find(strcmp(XLSRaw(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSRaw(1,:),'Channel2'));


%Find the corresponding entry in the XLS file
if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))&(isempty(findstr(Prefix,'Bcd1x')))
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


if (strcmp(XLSRaw(XLSEntry,Channel2Column),'His-RFP'))|...
        (strcmp(XLSRaw(XLSEntry,Channel1Column),'His-RFP'))|...
        (strcmp(XLSRaw(XLSEntry,Channel2Column),'MCP-TagRFP(1)'))
    nc9=XLSRaw{XLSEntry,nc9Column};
    nc10=XLSRaw{XLSEntry,nc10Column};
    nc11=XLSRaw{XLSEntry,nc11Column};
    nc12=XLSRaw{XLSEntry,nc12Column};
    nc13=XLSRaw{XLSEntry,nc13Column};
    nc14=XLSRaw{XLSEntry,nc14Column};
    CF=XLSRaw{XLSEntry,CFColumn};
end

%Do we need to convert any NaN chars into doubles?
if strcmp(lower(nc14),'nan')
    nc14=nan;
end
if strcmp(lower(nc13),'nan')
    nc13=nan;
end
if strcmp(lower(nc12),'nan')
    nc12=nan;
end
if strcmp(lower(nc11),'nan')
    nc11=nan;
end
if strcmp(lower(nc10),'nan')
    nc10=nan;
end
if strcmp(lower(nc9),'nan')
    nc9=nan;
end


%Convert the prefix into the string used in the XLS file
Dashes=findstr(Prefix,'-');

%Find the corresponding entry in the XLS file
if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))&(isempty(findstr(Prefix,'Bcd1')))
    XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
        [Date,'\BcdGFP-HisRFP']));
else
    XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
        [Prefix(1:Dashes(3)-1),filesep,Prefix(Dashes(3)+1:end)]));
end
ncs=[nc9,nc10,nc11,nc12,nc13,nc14];

if (length(find(isnan(ncs)))==length(ncs))|(length(ncs)<6)
    error('Have the ncs been defined in MovieDatabase.XLSX?')
end

%Now do the nuclear segmentation and lineage tracking. This should be put
%into an independent function.

%Create the cell array with the names.
D=dir([FISHPath,filesep,'Data',filesep,Prefix,filesep,'*His*.tif']);
for i=1:length(D)
    names{i}=[FISHPath,filesep,'Data',filesep,Prefix,filesep,D(i).name];
end


%Pull the mitosis information from ncs.
ncs=ncs(ncs~=0);

%Note that I'm adding a range of two frames frames before and after the
%determines mitosis
indMit=[ncs'-2,ncs'+2];

%Make sure no indices are negative. This could happen is the nuclear cycle
%started at frame 1, for example.
indMit(indMit<1)=1;

%If we don't have nc14 we'll fool the code into thinking that the last
%frame of the movie was nc14
if isnan(indMit(end,1))
    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
    indMit(end,1)=length(FrameInfo)-6;
    indMit(end,2)=length(FrameInfo)-5;
end    

%If indMit is empty this means that we didn't see any mitosis. Deal with
%this by assigning the first frames to it
if isempty(indMit)
    indMit=[1,2];
end


%Make sure to edit getdefaultParameters.m to change the pixel size
%parameters!!

%Embryo mask
ImageTemp=imread(names{1});
embryo_mask=true(size(ImageTemp));
clear ImageTemp

%Do the tracking for the first time
if ~exist([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
    [nuclei, centers, Dummy, dataStructure] = ...
        mainTracking(names,'indMitosis',indMit,'embryoMask', embryo_mask);
    % names is a cell array containing the names of all frames in the movie in order.
    % indMitosis is an nx2 array containing the first and last frame of mitosis in every row.
    % embryoMask is the possible mask of the embryo. If no embryo edge is visible,
    % true(size(an_image_from_the_movie)) can be given as input.

    % Convert the results to compatible structures and save them
    %Put circles on the nuclei
    [Ellipses] = putCirclesOnNuclei(centers,names,indMit);
    %Convert nuclei structure into schnitzcell structure
    [schnitzcells] = convertNucleiToSchnitzcells(nuclei); 
else
    %Do re-tracking
    display('Re-tracking...')
    warning('Doing re-tracking. This code needs to be able to handle approved schnitz still')
  
    %Load the Ellipses and re-generate the centers
    load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses')
    %centers = updateCentersFromEllipses(Ellipses, centers);
    centers = updateCentersFromEllipses(Ellipses);

    %Load the dataStructure to seed up retracking if it exists
    if exist([FISHPath,filesep,'Analysis',filesep,Prefix,'_',filesep,'dataStructure.mat'])
        load([FISHPath,filesep,'Analysis',filesep,Prefix,'_',filesep,'dataStructure.mat'])
    elseif exist([FISHPath,filesep,'Analysis',filesep,Prefix,'_',filesep,'TrackingDataStructure.mat'])
        load([FISHPath,filesep,'Analysis',filesep,Prefix,'_',filesep,'TrackingDataStructure.mat'])
    end

    % look for a settings file in the Raw Data folder.
    if exist([FISHPath,filesep,'Data',filesep,Prefix,filesep,Prefix,'-AcquisitionSettings.mat'],'file')
        load([FISHPath,filesep,'Data',filesep,Prefix,filesep,Prefix,'-AcquisitionSettings.mat']);
        fields = fieldnames(AcquisitionSettings);
        settingArguments = cell(1,2*length(fields));
        for i=1:length(fields)
            settingArguments{2*i-1} = fields{i}; % field names must match the parameter names in mainTracking
            settingArguments{2*i} = getfield(AcquisitionSettings, fields{i});
            % alter argument names if necessary
            if strcmp(settingArguments{2*i-1}, 'time_resolution')
                settingArguments{2*i-1} = 'time resolution';
            elseif strcmp(settingArguments{2*i-1}, 'space_resolution')
                settingArguments{2*i-1} = 'space resolution';
            end
        end
    else
        settingArguments = {};
    end
    
    
    %Re-run the tracking
    if exist('dataStructure')
        %Edit the names in dataStructure to match the current folder setup
        dataStructure.names=names;
        
        [nuclei, centers, Dummy, dataStructure] = mainTracking(...
            names,'indMitosis',indMit,'embryoMask', embryo_mask,...
            'centers',centers,'dataStructure',dataStructure, settingArguments{:});
    else
        [nuclei, centers, Dummy, dataStructure] = mainTracking(...
            names,'indMitosis',indMit,'embryoMask', embryo_mask,...
            'centers',centers, settingArguments{:});
    end

    %Put circles on the nuclei
    [Ellipses] = putCirclesOnNuclei(centers,names,indMit);
    %Convert nuclei structure into schnitzcell structure
    [schnitzcells] = convertNucleiToSchnitzcells(nuclei); 
end

%Save the information
%Now save
mkdir([DropboxFolder,filesep,Prefix])
save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses')
save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells')
save([FISHPath,filesep,'Analysis',filesep,Prefix,'_',filesep,'dataStructure.mat'],'dataStructure')
