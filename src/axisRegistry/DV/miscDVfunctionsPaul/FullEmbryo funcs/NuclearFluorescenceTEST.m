[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);



%Determine division times
%Load the information about the nc from moviedatabase file
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

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

                for j=1:TotalSchnitz

                    %HG: Note that I'm calling a function here so that I can
                    %debug the parfor loop above. Ideally, I would have 
                    %parfor loop over images, not schnitzes within an image.
                    %However, I couldn't quite figure out how to do that.
                    %schnitzcells(j)=ExtractNuclearFluorescence(schnitzcells(j),...
                    %    CurrentFrame,...
                    %    Image,LinesPerFrame,PixelsPerLine,NumberSlices2,Circle,IntegrationRadius,InputChannel(ChN));
                    schnitzcells(j)=ExtractNuclearFluorescence(schnitzcells(j),...
                        CurrentFrame,...
                        Image,LinesPerFrame,PixelsPerLine,NumberSlices2,Circle,IntegrationRadius,ChN);
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

if ~exist([FISHPath,filesep,Prefix,'_'], 'dir')
    mkdir([FISHPath,filesep,Prefix,'_'])
end
save([FISHPath,filesep,Prefix,'_',filesep,'dataStructure.mat'],'dataStructure')


%Stitch the schnitzcells using Simon's code
if ~SkipStitchSchnitz
    disp 'Skipping StitchSchnitz'
    StitchSchnitz(Prefix)
end