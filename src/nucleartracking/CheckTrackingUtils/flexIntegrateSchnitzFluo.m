function schnitzcells =...
    flexIntegrateSchnitzFluo(Prefix, radii, schnitzcells, FrameInfo)

saveFlag = true;

liveExperiment = LiveExperiment(Prefix);

if nargin == 1
        
    DropboxFolder = liveExperiment.userResultsFolder;
    Channels = liveExperiment.Channels;
    
    FrameInfo = getFrameInfo(liveExperiment); 
    
    schnitzcells = getSchnitzcells(liveExperiment);

    saveFlag = true;
    
    radii = [2];
    
    
elseif nargin == 2  
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

movieMat = getMovieMat(liveExperiment);

numFrames = length(FrameInfo);

default_radius = radii(1);


%Create the circle that we'll use as the mask
for i=1:length(radii)
    IntegrationRadius = radii(i);
    IntegrationRadius=floor(IntegrationRadius/FrameInfo(1).PixelSize); %Radius of the integration in pixels
    if ~mod(IntegrationRadius,2)
        IntegrationRadius=IntegrationRadius+1;
    end
    Circle=false(3*IntegrationRadius,3*IntegrationRadius);
    Circle=double(MidpointCircle(Circle,IntegrationRadius,1.5*IntegrationRadius+0.5,...
        1.5*IntegrationRadius+0.5,1));
    Fluolabel = ['Fluo' , strrep(num2str(radii(i)), '.', '_'), 'um'];
    %Initialize fields
    if i == 1
        schnitzcells(1).Fluo = [];
    end
    schnitzcells(1).(Fluolabel) = [];
    

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
                % GM: added the uint16 conversion 
                if ~isempty(movieMat)
                    imStack = double(movieMat(:,:,:, CurrentFrame, ChN));
                else
                    imStack = double(getMovieFrame(liveExperiment, CurrentFrame, ChN));
                end
                % 
                convImage = imfilter(imStack, Circle, 'same');
                convImage(edgeMask) = NaN;

                for j=1:length(tempSchnitz)

                    cenx=min(max(1,round(tempSchnitz(j).cenx(tempSchnitz(j).frames==CurrentFrame))),PixelsPerLine);
                    ceny=min(max(1,round(tempSchnitz(j).ceny(tempSchnitz(j).frames==CurrentFrame))),LinesPerFrame);
                    
                    tempSchnitz(j).(Fluolabel)(tempSchnitz(j).frames==CurrentFrame,:,chIndex) = single(convImage(ceny,cenx,:));
                    if i == 1
                        tempSchnitz(j).Fluo(tempSchnitz(j).frames==CurrentFrame,:,chIndex) = ...
                            tempSchnitz(j).(Fluolabel)(tempSchnitz(j).frames==CurrentFrame,:,chIndex);
                    end

                end %loop of nuclei in a frame

            end %loop over frames

            schnitzcells = tempSchnitz;


            try close(h); end

        end %loop over channels

    else 

        error('Input channel not recognized. Check correct definition in MovieDatabase.Input channels should use the :input notation.');

    end 
end
if saveFlag
    if whos(var2str(schnitzcells)).bytes < 2E9
        save(schnitzPath, 'schnitzcells', '-v6');
    else
        save(schnitzPath, 'schnitzcells', '-v7.3', '-nocompression');
    end

end

end