function ImportTr2d(Prefix)

%This function grabs the data exported by tr2d and saves it into the
%corresponding Prefix in the format of Ellipses and lineages.

%Get folder and movie length information for Prefix
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);

%JAKE: we need to add this for getting the Dropbox folder
%Get the folders, including the default Dropbox one
[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);

%JAKE: Import the information about experiment
%Load the information about the nc from moviedatabase file
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
%Check whether nc14 occurred very close to the end of the movie. For those
%frames we'll move the boundary for tracking purposes

%JAKE: Added num frames
numFrames = length(FrameInfo);

%File names with exported information from tr2d
ObjectFile='tr2d_objects.csv';
TrackFile='tr2d_tracks.csv';
tr2dExportFolder=[PreProcPath,filesep,Prefix,filesep,'tr2dProject',filesep,...
    'mRNADynamicsExport',filesep];

%Load the tr2d files
Objects=csvread([tr2dExportFolder,ObjectFile],2,0);    %Read starting at the second row
Tracks=csvread([tr2dExportFolder,TrackFile],2,0);

%Create the Ellipses structure. The format is:
%(x, y, a, b, theta, maxcontourvalue, time, particle_id)

%The tr2d object format is (note that time and IDs start from 0):
%# t	 id	 area	 com_x	 com_y	 angle	 r1	 r2

%Note that we're mapping Area in tr2d to maxcontourvalue in Ellipses. We
%were not using the latter at all in the code.

Ellipses=cell(length(FrameInfo),1);
for i=1:length(FrameInfo)
    %Find all objects in this frame
    FrameFilter=(Objects(:,1)==i-1);
    %Copy the information
    Ellipses{i}=Objects(FrameFilter,[4,5,7,8,6,3]);
    Ellipses{i}(:,7)=i; %Frame
    Ellipses{i}(:,8)=Objects(FrameFilter,2)+1;
end

%Create the schnitzcells structure. The format is:
% P
% E
% D
% frames
% cenx
% ceny
% len
% cellno

%The tr2d format is
%# tracklet_id	 parent_tracklet_id	 child_tracklat_id1	 child_tracklat_id2	 (time	 object_id)...
%Note that time and ID are given as a pair of numbers. Also, the csv is
%padded with zeroes at the end of the list of frames and IDs.

for i=1:size(Tracks,1)
    %Find where the list of frames and IDs ends. I need to be careful with
    %to potentially problematic cases
    %1) There is no padding because this is the longest track
    %2) The track exists for only the first frame (denoted by 0 in tr2d)
    
    %Check that this track is not one of the longest ones. To do this,
    %we'll look for the first pair of [0,0] in the frame/index part of
    %Tracks
    
    if ~isempty(findstr(Tracks(i,7:end),[0,0]))
        MaxColumn=min(findstr(Tracks(i,7:end),[0,0]))+5;
        %If MaxColumn is an odd number, it means there was a 0 for the id
        %before the last pair of zeros that mark the end of the track
        if mod(MaxColumn,2)
            MaxColumn=MaxColumn+1;
        end
    else
        MaxColumn=size(Tracks,2);
    end
    
    %Actual track information:
    schnitzcells(i).frames=Tracks(i,5:2:MaxColumn)+1;
    schnitzcells(i).cellno=Tracks(i,6:2:MaxColumn)+1;
    schnitzcells(i).P=Tracks(i,2)+1;
    schnitzcells(i).D=Tracks(i,3)+1;
    schnitzcells(i).E=Tracks(i,4)+1;
    
    %Information coming from Ellipses. Our code should be modified to not
    %need to use duplicated information
    for j=1:length(schnitzcells(i).frames)
        EllipseToPull=...
            find(Ellipses{schnitzcells(i).frames(j)}(:,8)==...
            schnitzcells(i).cellno(j));
        
        schnitzcells(i).cenx(j)=...
            Ellipses{schnitzcells(i).frames(j)}(EllipseToPull,1);
        schnitzcells(i).ceny(j)=...
            Ellipses{schnitzcells(i).frames(j)}(EllipseToPull,2);
        %Calculate the length as the major axis length
        schnitzcells(i).len(j)=...
            max([Ellipses{schnitzcells(i).frames(j)}(EllipseToPull,3),...
                Ellipses{schnitzcells(i).frames(j)}(EllipseToPull,4)]);
    end
end

%JAKE: Added YJK's code for extracting nuclear fluorescence
%Extract the nuclear fluorescence values if we're in the right experiment
%type
if strcmpi(ExperimentType,'inputoutput')||strcmpi(ExperimentType,'input')
    
    
    %Parse the channel information for the different experiment types
    if strcmpi(ExperimentType,'inputoutput')
        InputChannelTemp=strfind({lower(Channel1{1}),lower(Channel2{1})},'mcp');
        if isempty(InputChannelTemp)
            InputChannelTemp=strfind({lower(Channel1{1}),lower(Channel2{1})},'pp7');
                if isempty(InputChannelTemp)
                    InputChannelTemp=strfind({lower(Channel1{1}),lower(Channel2{1})},'lambdan');
                end    
        end
        InputChannelTemp=cellfun(@isempty,InputChannelTemp);
        InputChannel = find(InputChannelTemp);
        
    elseif strcmpi(ExperimentType,'input')
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
    Circle=false(3*IntegrationRadius,3*IntegrationRadius);
    Circle=MidpointCircle(Circle,IntegrationRadius,1.5*IntegrationRadius+0.5,...
        1.5*IntegrationRadius+0.5,1);
    
    %Add the Fluo and Mask field to schnitzcells. This is done now so that
    %the parfor doesn't freak out
    schnitzcells(1).Fluo=[];
    schnitzcells(1).Mask=[];
    
    
    if sum(InputChannel)
                
        %Extract the fluroescence of each schnitz, for each channel,
        %at each time point

        TotalSchnitz=length(schnitzcells);
         
        %Get the image dimensions
        PixelsPerLine=FrameInfo(1).PixelsPerLine;
        LinesPerFrame=FrameInfo(1).LinesPerFrame;
        %Number of z-slices
        NumberSlices=FrameInfo(1).NumberSlices;
        NumberSlices2 = NumberSlices + 2;
        
        %Initialize blank images
        Mask=false(LinesPerFrame,PixelsPerLine);
        Image=zeros(LinesPerFrame,PixelsPerLine,NumberSlices2); 

        
        for ChN=1:length(InputChannel)
        
            h=waitbar(0,['Extracting nuclear fluorescence for channel ',num2str(ChN)]);
            
                if strcmpi(ExperimentType,'inputoutput') || (strcmpi(ExperimentType,'input')&&(length(InputChannel))==1)
                    nameSuffix=['_ch',iIndex(InputChannel,2)];
                elseif strcmpi(ExperimentType,'input')&&(length(InputChannel))>1
                    nameSuffix=['_ch',iIndex(InputChannel(ChN),2)];
                else
                    error('Not sure what happened here. Talk to YJK or SA.')
                end
                
            for CurrentFrame=1:numFrames

                waitbar(CurrentFrame/numFrames,h);
% 
%                 %Initialize the image
%                 Image=zeros(LinesPerFrame,PixelsPerLine,NumberSlices2)
%                 %AR 1/12/18 moved this outside the for loops
                
                %Load the z-stack for this frame

                for CurrentZ=1:NumberSlices2   %Note that I need to add the two extra slices manually
                         Image(:,:,CurrentZ)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),nameSuffix,'.tif']);
                end

                parfor j=1:TotalSchnitz

                    %HG: Note that I'm calling a function here so that I can
                    %debug the parfor loop above. Ideally, I would have 
                    %parfor loop over images, not schnitzes within an image.
                    %However, I couldn't quite figure out how to do that.
                    schnitzcells(j)=ExtractNuclearFluorescence(schnitzcells(j),...
                        CurrentFrame,...
                        Image,LinesPerFrame,PixelsPerLine,NumberSlices2,Circle,IntegrationRadius,InputChannel(ChN));
                end
            end
        close(h)
        end
    else
        error('Input channel not recognized. Check correct definition in MovieDatabase')
    end
end
    

%Save the information
%Now save
mkdir([DropboxFolder,filesep,Prefix])
save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses')

if strcmpi(ExperimentType,'inputoutput')||strcmpi(ExperimentType,'input')
    %Change the name of the Circle variable to make it more understandable when
    %loaded independently
    IntegrationArea=Circle;
    save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells','IntegrationArea')
else
    save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells')
end




