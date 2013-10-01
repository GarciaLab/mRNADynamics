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
[Num,Txt]=xlsread([DefaultDropboxFolder,'\HGMovieDatabaseV2.xlsx']);
XLSHeaders=Txt(1,:);
Txt=Txt(2:end,:);

%Find the different columns.
DataFolderColumn=find(strcmp(XLSHeaders,'DataFolder'));
nc9Column=find(strcmp(XLSHeaders,'nc9'));
nc10Column=find(strcmp(XLSHeaders,'nc10'));
nc11Column=find(strcmp(XLSHeaders,'nc11'));
nc12Column=find(strcmp(XLSHeaders,'nc12'));
nc13Column=find(strcmp(XLSHeaders,'nc13'));
nc14Column=find(strcmp(XLSHeaders,'nc14'));
CFColumn=find(strcmp(XLSHeaders,'CF'));
Channel2Column=find(strcmp(XLSHeaders,'Channel2'));

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


if strcmp(Txt(XLSEntry,Channel2Column),'His-RFP')
    nc9=Num(XLSEntry,nc9Column-6);
    nc10=Num(XLSEntry,nc10Column-6);
    nc11=Num(XLSEntry,nc11Column-6);
    nc12=Num(XLSEntry,nc12Column-6);
    nc13=Num(XLSEntry,nc13Column-6);
    nc14=Num(XLSEntry,nc14Column-6);
    %This is in case the last column for CF is all nan and is not part of
    %the Num matrix
    if size(Num,2)==CFColumn-6    
        CF=Num(XLSEntry,CFColumn-6);
    else
        CF=nan;
    end
else
    error('nc information not define in HGMovieDatabase.xlsx')
end
ncs=[nc9,nc10,nc11,nc12,nc13,nc14];


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
    load([DropboxFolder,filesep,Prefix,'\Ellipses.mat'],'Ellipses')
    %centers = updateCentersFromEllipses(Ellipses, centers);
    centers = updateCentersFromEllipses(Ellipses);

    %Load the dataStructure to seed up retracking if it exists
    if exist([FISHPath,filesep,'Analysis',filesep,Prefix,'_',filesep,'dataStructure.mat'])
        load([FISHPath,filesep,'Analysis',filesep,Prefix,'_',filesep,'dataStructure.mat'])
    elseif exist([FISHPath,filesep,'Analysis',filesep,Prefix,'_',filesep,'TrackingDataStructure.mat'])
        load([FISHPath,filesep,'Analysis',filesep,Prefix,'_',filesep,'TrackingDataStructure.mat'])
    end

    %Re-run the tracking
    if exist('dataStructure')
        %Edit the names in dataStructure to match the current folder setup
        dataStructure.names=names;
        
        [nuclei, centers, Dummy, dataStructure] = mainTracking(...
            names,'indMitosis',indMit,'embryoMask', embryo_mask,...
            'centers',centers,'dataStructure',dataStructure);
    else
        [nuclei, centers, Dummy, dataStructure] = mainTracking(...
            names,'indMitosis',indMit,'embryoMask', embryo_mask,...
            'centers',centers);
    end

    %Put circles on the nuclei
    [Ellipses] = putCirclesOnNuclei(centers,names,indMit);
    %Convert nuclei structure into schnitzcell structure
    [schnitzcells] = convertNucleiToSchnitzcells(nuclei); 
end

%Save the information
%Now save
mkdir([DropboxFolder,filesep,Prefix])
save([DropboxFolder,filesep,Prefix,'\Ellipses.mat'],'Ellipses')
save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells')
save([FISHPath,filesep,'Analysis',filesep,Prefix,'_',filesep,'dataStructure.mat'],'dataStructure')
