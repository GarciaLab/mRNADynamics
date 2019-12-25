function schnitzcells =...
    integrateSchnitzFluo(Prefix, schnitzcells, FrameInfo,...
    ExperimentType, Channels, PreProcPath)

saveFlag = false;
if nargin == 1
    
    [rawDataPath,ProcPath,DropboxFolder,MS2CodePath, PreProcPath,...
        rawDataFolder, Prefix, ExperimentType,Channel1,Channel2,OutputFolder,...
        Channel3, spotChannels, MovieDataBaseFolder, movieDatabase]...
        = readMovieDatabase(Prefix);
    
    Channels = {Channel1{1},Channel2{1}, Channel3{1}};
    
    [Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution, ...
        Channel1, Channel2, Objective, Power, DataFolderColumnValue, ~, Comments, ...
        nc9, nc10, nc11, nc12, nc13, nc14, CF, Channel3, prophase, metaphase] = getExperimentDataFromMovieDatabase(Prefix, movieDatabase);
  
    DataFolder = [DropboxFolder, filesep, Prefix];

     load([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo');
    
    FilePrefix = [Prefix, '_'];
    schnitzPath = [DropboxFolder, filesep, FilePrefix(1:end - 1), filesep, FilePrefix(1:end - 1), '_lin.mat'];
    
    disp('Loading schnitzcells...')
    load(schnitzPath, 'schnitzcells');
    disp('schnitzcells loaded.')
    
    saveFlag = true;
    
end

numFrames = length(FrameInfo);

%Parse the channel information for the different experiment types
if strcmpi(ExperimentType,'inputoutput')
    
    InputChannelTemp1 = strfind({lower(Channels{1}),lower(Channels{2}), lower(Channels{3})},'input');
    InputChannelTemp2=~cellfun(@isempty,InputChannelTemp1);
    InputChannel = find(InputChannelTemp2);
    
    
elseif strcmpi(ExperimentType,'input')
    %Parse the information from the different channels
    
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

%Initialize fields
schnitzcells(1).Fluo = [];

if sum(InputChannel)
    
    %Extract the fluroescence of each schnitz, for each channel,
    %at each time point
    
    %Get the image dimensions
    PixelsPerLine=FrameInfo(1).PixelsPerLine;
    LinesPerFrame=FrameInfo(1).LinesPerFrame;
    %Number of z-slices
    NumberSlices=FrameInfo(1).NumberSlices;
    NumberSlices2 = NumberSlices + 2;
    
    %Generate reference frame for edge detection
    refFrame = ones(LinesPerFrame,PixelsPerLine,NumberSlices2);
    convRef = convn(refFrame, Circle, 'same');
    edgeMask = convRef~=sum(Circle(:));
    
    for ChN=1:length(InputChannel)
        
        h=waitbar(0,['Extracting nuclear fluorescence for channel ',num2str(ChN)]);
        
        if strcmpi(ExperimentType,'inputoutput') || (strcmpi(ExperimentType,'input')&&(length(InputChannel))==1)
            nameSuffix=['_ch',iIndex(InputChannel,2)];
        elseif strcmpi(ExperimentType,'input')&&(length(InputChannel))>1
            nameSuffix=['_ch',iIndex(InputChannel(ChN),2)];
        else
            error('Not sure what happened here. Talk to YJK or SA.');
        end
        % NL: should parallelize this
        
        tempSchnitz = schnitzcells;
        for CurrentFrame=1:numFrames
            
            waitbar(CurrentFrame/numFrames,h);
            %
            %                 %Initialize the image
            Image=zeros(LinesPerFrame,PixelsPerLine,NumberSlices2);
            
            %Load the z-stack for this frame
            for CurrentZ=1:NumberSlices2   %Note that I need to add the two extra slices manually
                Image(:,:,CurrentZ)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),nameSuffix,'.tif']);
            end
            
            convImage = imfilter(Image, double(Circle), 'same');
            convImage(edgeMask) = NaN;
            for j=1:length(tempSchnitz)
                CurrentIndex=find(tempSchnitz(j).frames==CurrentFrame);
                cenx=min(max(1,round(tempSchnitz(j).cenx(CurrentIndex))),PixelsPerLine);
                ceny=min(max(1,round(tempSchnitz(j).ceny(CurrentIndex))),LinesPerFrame);
                tempSchnitz(j).Fluo(CurrentIndex,1:NumberSlices2,ChN) = single(convImage(ceny,cenx,:));
            end
            
        end
        schnitzcells = tempSchnitz;
        
        
        close(h);
    end
else
    error('Input channel not recognized. Check correct definition in MovieDatabase. Input channels should use the :input notation.');
end

if saveFlag
    save(schnitzPath);
end

end