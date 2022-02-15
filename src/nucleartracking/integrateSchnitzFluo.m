function schnitzcells =...
    integrateSchnitzFluo(Prefix, schnitzcells, FrameInfo)

saveFlag = false;

liveExperiment = LiveExperiment(Prefix);

if nargin == 1
        
    DropboxFolder = liveExperiment.userResultsFolder;
    Channels = liveExperiment.Channels;
    
    FrameInfo = getFrameInfo(liveExperiment); 
    
    schnitzcells = getSchnitzcells(liveExperiment);

    saveFlag = true;
    
end

schnitzPath = [liveExperiment.resultsFolder, filesep, Prefix, '_lin.mat'];

Channels = liveExperiment.Channels;

InputChannelIndexes = find(contains(Channels, 'input', 'IgnoreCase', true));

if isempty(InputChannelIndexes)
    warning(['No input channel found. Check correct definition in MovieDatabase.',...
        ' Input channels should use the :input notation.'])
    return;
end

movie = getMovieMat(liveExperiment);
%assume there's just one input channel
movie = double(movie(:, :, :, :, InputChannelIndexes));
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
schnitzcells(1).Fluo = [];

if sum(InputChannelIndexes)
    
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
    chIndex = 0;
    for ChN=InputChannelIndexes
        chIndex = chIndex+1;
        h=waitbar(0,['Extracting nuclear fluorescence for channel ',num2str(ChN)]);
        
%         nameSuffix=['_ch',iIndex(ChN,2)];
        
        % NL: should parallelize this
        
        tempSchnitz = schnitzcells;
        for CurrentFrame=1:numFrames
            
            try waitbar(CurrentFrame/numFrames,h); catch; end
            %
            %                 %Initialize the image
            %             Image=zeros(LinesPerFrame,PixelsPerLine,nSlices);
            
            %             %Load the z-stack for this frame
            %             for CurrentZ=1:nSlices   %Note that I need to add the two extra slices manually
            %                 Image(:,:,CurrentZ)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),nameSuffix,'.tif']);
            %             end
            
            
            convImage = imfilter(movie(:,:,:, CurrentFrame), Circle, 'same');
            convImage(edgeMask) = NaN;
            
            for j=1:length(tempSchnitz)
                
                cenx=min(max(1,round(tempSchnitz(j).cenx(tempSchnitz(j).frames==CurrentFrame))),PixelsPerLine);
                ceny=min(max(1,round(tempSchnitz(j).ceny(tempSchnitz(j).frames==CurrentFrame))),LinesPerFrame);
                tempSchnitz(j).Fluo(tempSchnitz(j).frames==CurrentFrame,:,chIndex) = single(convImage(ceny,cenx,:));
                    
                    
            end %loop of nuclei in a frame
            
        end %loop over frames
        
        schnitzcells = tempSchnitz;
        
        
        try close(h); end
    
    end %loop over channels
    
else 
    
    error('Input channel not recognized. Check correct definition in MovieDatabase.Input channels should use the :input notation.');

end 

if saveFlag
    if whos(var2str(schnitzcells)).bytes < 2E9
        save(schnitzPath, 'schnitzcells', '-v6');
    else
        save(schnitzPath, 'schnitzcells', '-v7.3', '-nocompression');
    end
    
end

end