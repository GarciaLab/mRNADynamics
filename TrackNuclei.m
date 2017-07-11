function TrackNuclei(Prefix,varargin)

%This function is just a script that call Laurent's tracking code

%Options:
%StitchSchnitz: Run the schnitzcells fixing code by Simon

SkipStitchSchnitz=1;
if ~isempty(varargin)
    if strcmpi(varargin{1},'stitchschnitz')
        SkipStitchSchnitz=0;
    else
        error('Input parameter not recognized')
    end
end



%Get the folders, including the default Dropbox one
[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
%Now get the actual DropboxFolder
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);



%Determine division times
%Load the information about the nc from the XLS file
[Num,Txt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
XLSHeaders=Txt(1,:);
Txt=Txt(2:end,:);

ExperimentTypeColumn=find(strcmp(XLSRaw(1,:),'ExperimentType'));
ExperimentAxisColumn=find(strcmp(XLSRaw(1,:),'ExperimentAxis'));
Channel1Column=find(strcmp(XLSRaw(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSRaw(1,:),'Channel2'));

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
Channel1=XLSRaw(PrefixRow,Channel1Column);
Channel2=XLSRaw(PrefixRow,Channel2Column);
%If Channel2 was left empty, it would contain a NaN, which will cause
%problems below. In that case, replace it by an empty string.
if isnan(Channel2{1})
    Channel2{1}='';
end

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
XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
    [Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));

if isempty(XLSEntry)
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
    if isempty(XLSEntry)
        disp('%%%%%%%%%%%%%%%%%%%%%')
        error('Dateset could not be found. Check MovieDatabase.xlsx')
        disp('%%%%%%%%%%%%%%%%%%%%%')
    end
end


nc9=XLSRaw{XLSEntry,nc9Column};
nc10=XLSRaw{XLSEntry,nc10Column};
nc11=XLSRaw{XLSEntry,nc11Column};
nc12=XLSRaw{XLSEntry,nc12Column};
nc13=XLSRaw{XLSEntry,nc13Column};
nc14=XLSRaw{XLSEntry,nc14Column};
CF=XLSRaw{XLSEntry,CFColumn};


if ~exist('nc9','var')
    error('Cannot find nuclear cycle values. Were they defined in MovieDatabase.XLSX?')
end

%This checks whether all ncs have been defined
ncCheck=[nc9,nc10,nc11,nc12,nc13,nc14];
if length(ncCheck)~=6
    error('Check the nc frames in the MovieDatabase entry. Some might be missing')
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
XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
        [Prefix(1:Dashes(3)-1),filesep,Prefix(Dashes(3)+1:end)]));

%HG: Did we actually need these cases down here?
% if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
%         (isempty(findstr(Prefix,'NoBcd')))&(isempty(findstr(Prefix,'Bcd1')))&(isempty(findstr(Prefix,'Bcd4x')))
%     XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
%         [Date,'\BcdGFP-HisRFP']));
% else
%     XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
%         [Prefix(1:Dashes(3)-1),filesep,Prefix(Dashes(3)+1:end)]));
% end


ncs=[nc9,nc10,nc11,nc12,nc13,nc14];

if (length(find(isnan(ncs)))==length(ncs))|(length(ncs)<6)
    error('Have the ncs been defined in MovieDatabase.XLSX?')
end

%Now do the nuclear segmentation and lineage tracking. This should be put
%into an independent function.

%Create the cell array with the names.
D=dir([PreProcPath,filesep,Prefix,filesep,Prefix,'-His_*.tif']);
names = cell(length(D), 1);
for i=1:length(D)
    names{i}=[PreProcPath,filesep,Prefix,filesep,D(i).name];
end


%Pull the mitosis information from ncs.
ncs=ncs(ncs~=0);

%Note that I'm adding a range of two frames frames before and after the
%determines mitosis
indMit=[ncs'-2,ncs'+2];

%Make sure no indices are negative. This could happen is the nuclear cycle
%started at frame 1, for example.
indMit(indMit<1)=1;

%Check whether nc14 occurred very close to the end of the movie. For those
%frames we'll move the boundary for tracking purposes
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
indMit(indMit>=length(FrameInfo))=indMit(indMit>=length(FrameInfo))-2;


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


  
%Get information about the spatial and temporal resolution
settingArguments{1}='time resolution';
settingArguments{2}=median(diff([FrameInfo.Time]));     %Median separation between frames (in seconds)
settingArguments{3}='space resolution';
settingArguments{4}=FrameInfo(1).PixelSize;




%Do the tracking for the first time
if ~exist([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
    
    
    [nuclei, centers, Dummy, dataStructure] = ...
        mainTracking(names,'indMitosis',indMit,'embryoMask', embryo_mask,...
        settingArguments{:});
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
    if exist([FISHPath,filesep,Prefix,'_',filesep,'dataStructure.mat'])
        load([FISHPath,filesep,Prefix,'_',filesep,'dataStructure.mat'])
    elseif exist([FISHPath,filesep,Prefix,'_',filesep,'TrackingDataStructure.mat'])
        load([FISHPath,filesep,Prefix,'_',filesep,'TrackingDataStructure.mat'])
    end

    % look for a settings file in the Raw Data folder.
    if exist([PreProcPath,filesep,Prefix,filesep,Prefix,'-AcquisitionSettings.mat'],'file')
        load([PreProcPath,filesep,Prefix,filesep,Prefix,'-AcquisitionSettings.mat']);
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

%Add the radius information to the schnitz
for i=1:length(schnitzcells)
    for j=1:length(schnitzcells(i).frames)
        schnitzcells(i).len(:)=...
            mean(Ellipses{schnitzcells(i).frames(j)}(schnitzcells(i).cellno(j),3:4));
    end
end


%Save everything at this point. It will be overwritten later, but it's
%useful for debugging purposes if there's a bug in the code below.
mkdir([DropboxFolder,filesep,Prefix])
save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses')
save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells')
if ~exist([FISHPath,filesep,Prefix,'_'])
    mkdir([FISHPath,filesep,Prefix,'_'])
end
save([FISHPath,filesep,Prefix,'_',filesep,'dataStructure.mat'],'dataStructure')



%Extract the nuclear fluorescence values if we're in the right experiment
%type
if strcmp(lower(ExperimentType),'inputoutput')|strcmp(lower(ExperimentType),'input')
    
    
    %Parse the channel information for the different experiment types
    if strcmp(lower(ExperimentType),'inputoutput')
        InputChannelTemp=strfind({lower(Channel1{1}),lower(Channel2{1})},'mcp');
        if isempty(InputChannelTemp)
            InputChannelTemp=strfind({lower(Channel1{1}),lower(Channel2{1})},'pp7');
                if isempty(InputChannelTemp)
                    InputChannelTemp=strfind({lower(Channel1{1}),lower(Channel2{1})},'lambdan');
                end    
        end
        InputChannelTemp=cellfun(@isempty,InputChannelTemp);
        InputChannel = find(InputChannelTemp);
        
    elseif strcmp(lower(ExperimentType),'input')
        %Parse the information from the different channels
        Channels={Channel1{1},Channel2{1}};

        %Histone channel.
        histoneChannel=find(~cellfun(@isempty,strfind(lower(Channels),'his')));
        if isempty(histoneChannel)
            histoneChannel=0;
        end

        %Input channels
        InputChannel=~cellfun(@isempty,Channels);
        if histoneChannel
            InputChannel(histoneChannel)=0;
        end
        InputChannel=find(InputChannel);
    end
    
    
    %Create the circle that we'll use as the mask
    IntegrationRadius=2;       %Radius of the integration region in um
    IntegrationRadius=floor(IntegrationRadius/FrameInfo(1).PixelSize); %Radius of the integration in pixels
    if ~mod(IntegrationRadius,2)
        IntegrationRadius=IntegrationRadius+1;
    end
    Circle=logical(zeros(3*IntegrationRadius,3*IntegrationRadius));
    Circle=MidpointCircle(Circle,IntegrationRadius,1.5*IntegrationRadius+0.5,...
        1.5*IntegrationRadius+0.5,1);
    
    %Add the Fluo and Mask field to schnitzcells. This is done now so that
    %the parfor doesn't freak out
    schnitzcells(1).Fluo=[];
    schnitzcells(1).Mask=[];
    
    
    if sum(InputChannel)
                
        %Extract the fluroescence of each schnitz, for each channel,
        %at each time point

        %Get the image dimensions
        PixelsPerLine=FrameInfo(1).PixelsPerLine;
        LinesPerFrame=FrameInfo(1).LinesPerFrame;
        %Number of z-slices
        NumberSlices=FrameInfo(1).NumberSlices;
        
        for ChN=1:length(InputChannel)
        
            h=waitbar(0,['Extracting nuclear fluorescence for channel ',num2str(ChN)]);
            for CurrentFrame=1:length(FrameInfo)

                waitbar(CurrentFrame/length(FrameInfo),h);

                %Initialize the image
                Image=zeros(LinesPerFrame,PixelsPerLine,NumberSlices+2);

                %Load the z-stack for this frame
                for CurrentZ=1:(NumberSlices+2)   %Note that I need to add the two extra slices manually
                    if strcmp(lower(ExperimentType),'inputoutput')
                        Image(:,:,CurrentZ)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'_ch',iIndex(InputChannel,2),'.tif']);
                    elseif strcmp(lower(ExperimentType),'input')&(length(InputChannel))==1
                        Image(:,:,CurrentZ)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'.tif']);
                    elseif strcmp(lower(ExperimentType),'input')&(length(InputChannel))>1
                        Image(:,:,CurrentZ)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'_ch',iIndex(InputChannel(ChN),2),'.tif']);
                    end

                end

                TotalSchnitz=length(schnitzcells);

                parfor j=1:TotalSchnitz

                    %HG: Note that I'm calling a function here so that I can
                    %debug the parfor loop above. Ideally, I would have 
                    %parfor loop over images, not schnitzes within an image.
                    %However, I couldn't quite figure out how to do that.
                    schnitzcells(j)=ExtractNuclearFluorescence(schnitzcells(j),...
                        CurrentFrame,...
                        Image,LinesPerFrame,PixelsPerLine,NumberSlices,Circle,IntegrationRadius,ChN);
                end
            end
        close(h)    
        end
    else
        error('Input channel not recognized. Check correct definition in MovieDatabase.XLSX')
    end
end
    

%Save the information
%Now save
mkdir([DropboxFolder,filesep,Prefix])
save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses')
%Change the name of the Circle variable to make it more understandble when
%loaded independently
IntegrationArea=Circle;
save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells','IntegrationArea')
if ~exist([FISHPath,filesep,Prefix,'_'])
    mkdir([FISHPath,filesep,Prefix,'_'])
end
save([FISHPath,filesep,Prefix,'_',filesep,'dataStructure.mat'],'dataStructure')


%Stitch the schnitzcells using Simon's code
if ~SkipStitchSchnitz
    display('Skipping StitchSchnitz')
    StitchSchnitz(Prefix)
end

