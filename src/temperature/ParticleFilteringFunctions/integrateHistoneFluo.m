function schnitzcells =...
    integrateHistoneFluo(Prefix, schnitzcells, FrameInfo)

saveFlag = true;

liveExperiment = LiveExperiment(Prefix);

if nargin == 1
    
    DropboxFolder = liveExperiment.userResultsFolder;
    Channels = liveExperiment.Channels;
    
    FrameInfo = getFrameInfo(liveExperiment);
    
    schnitzcells = getSchnitzcells(liveExperiment);
    
    saveFlag = true;
    
end
Channels = liveExperiment.Channels;
Channel1=Channels{1}; Channel2 = Channels{2}; Channel3 = Channels{3};

% Check how many channels have ":Nuclear" in the MovieDatabase.csv
NuclearChannels = [contains(Channel1, 'His', 'IgnoreCase', true),...
    contains(Channel2, 'His', 'IgnoreCase', true),...
    contains(Channel3, 'His', 'IgnoreCase', true)];

HistoneChannel = find(NuclearChannels == 1, 1);

schnitzPath = [liveExperiment.resultsFolder, filesep, Prefix, '_lin.mat'];



hisMat = getHisMat(liveExperiment);

numFrames = length(FrameInfo);


%Create the circle that we'll use as the mask
IntegrationRadius=2;       %Radius of the integration region in um
IntegrationRadius=floor(IntegrationRadius/FrameInfo(1).PixelSize); %Radius of the integration in pixels
if ~mod(IntegrationRadius,2)
    IntegrationRadius=IntegrationRadius+1;
end
Circle=false(3*IntegrationRadius,3*IntegrationRadius);
Circle=double(MidpointCircle(Circle,IntegrationRadius,1.5*IntegrationRadius+0.5,...
    1.5*IntegrationRadius+0.5,1));

%Initialize fields
schnitzcells(1).HistoneFluo = [];

if sum(NuclearChannels)
    
    %Extract the fluroescence of each schnitz, for each channel,
    %at each time point
    
    %Get the image dimensions
    PixelsPerLine=FrameInfo(1).PixelsPerLine;
    LinesPerFrame=FrameInfo(1).LinesPerFrame;
    %Number of z-slices
    nPadding = 2;
    
    nSlices=FrameInfo(1).NumberSlices + nPadding;
    
    
    %Generate reference frame for edge detection
    refFrame = ones(LinesPerFrame,PixelsPerLine, nSlices);
    convRef = convn(refFrame, Circle, 'same');
    edgeMask = convRef~=sum(Circle(:));

    
    tempSchnitz = schnitzcells;
    for CurrentFrame=1:numFrames
        
        try waitbar(CurrentFrame/numFrames,h); catch; end
        % GM: added the uint16 conversion
        
        imStack = double(getMovieFrame(liveExperiment, CurrentFrame, HistoneChannel));

        %
        convImage = imfilter(imStack, Circle, 'same');
        convImage(edgeMask) = NaN;
        
        for j=1:length(tempSchnitz)
            
            cenx=min(max(1,round(tempSchnitz(j).cenx(tempSchnitz(j).frames==CurrentFrame))),PixelsPerLine);
            ceny=min(max(1,round(tempSchnitz(j).ceny(tempSchnitz(j).frames==CurrentFrame))),LinesPerFrame);
            tempSchnitz(j).HistoneFluo(tempSchnitz(j).frames==CurrentFrame,:) = single(convImage(ceny,cenx,:));
            
            
        end %loop of nuclei in a frame
        
    end %loop over frames
    
    schnitzcells = tempSchnitz;
    
    
    try close(h); end
    
end %loop over channels


if saveFlag
    if whos(var2str(schnitzcells)).bytes < 2E9
        save(schnitzPath, 'schnitzcells', '-v6');
    else
        save(schnitzPath, 'schnitzcells', '-v7.3', '-nocompression');
    end
    
end

end