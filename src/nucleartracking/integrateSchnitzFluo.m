function schnitzcells =...
    integrateSchnitzFluo(Prefix, schnitzcells, FrameInfo,...
   PreProcPath)

saveFlag = false;

thisExperiment = liveExperiment(Prefix);

if nargin == 1
    
    [rawDataPath,ProcPath,DropboxFolder,MS2CodePath, PreProcPath,...
        rawDataFolder, Prefix, ~,Channel1,Channel2,OutputFolder,...
        Channel3, spotChannels, MovieDataBaseFolder, movieDatabase]...
        = readMovieDatabase(Prefix);
    
    Channels = {Channel1{1},Channel2{1}, Channel3{1}};
    
    DataFolder = [DropboxFolder, filesep, Prefix];
    
    load([DataFolder, filesep, 'FrameInfo.mat'], 'FrameInfo');
    
    FilePrefix = [Prefix, '_'];
    schnitzPath = [DropboxFolder, filesep, FilePrefix(1:end - 1), filesep, FilePrefix(1:end - 1), '_lin.mat'];
    
    disp('Loading schnitzcells...')
    load(schnitzPath, 'schnitzcells');
    disp('schnitzcells loaded.')
    
    saveFlag = true;
    
end

Channels = thisExperiment.Channels;

InputChannel = find(contains(Channels, 'input', 'IgnoreCase', true));

% load([PreProcPath, filesep, Prefix, filesep, Prefix, '_movieMat.mat'], 'movieMat');

movieMat = double(loadMovieMat(Prefix, 'chRange', InputChannel));
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

if sum(InputChannel)
    
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
    
    for ChN=1:length(InputChannel)
        
        h=waitbar(0,['Extracting nuclear fluorescence for channel ',num2str(ChN)]);
        
        nameSuffix=['_ch',iIndex(InputChannel,2)];
        
        % NL: should parallelize this
        
        tempSchnitz = schnitzcells;
        for CurrentFrame=1:numFrames
            
            waitbar(CurrentFrame/numFrames,h);
            %
            %                 %Initialize the image
            %             Image=zeros(LinesPerFrame,PixelsPerLine,nSlices);
            
            %             %Load the z-stack for this frame
            %             for CurrentZ=1:nSlices   %Note that I need to add the two extra slices manually
            %                 Image(:,:,CurrentZ)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),nameSuffix,'.tif']);
            %             end
            
            
            convImage = imfilter(movieMat(:,:,:, CurrentFrame), Circle, 'same');
            convImage(edgeMask) = NaN;
            
            for j=1:length(tempSchnitz)
%                 CurrentIndex=find(tempSchnitz(j).frames==CurrentFrame);
%                 cenx=min(max(1,round(tempSchnitz(j).cenx(CurrentIndex))),PixelsPerLine);
%                 ceny=min(max(1,round(tempSchnitz(j).ceny(CurrentIndex))),LinesPerFrame);
%                 tempSchnitz(j).Fluo(CurrentIndex,1:nSlices,ChN) = single(convImage(ceny,cenx,:));


                    %trying faster method
                    cenx=min(max(1,round(tempSchnitz(j).cenx(tempSchnitz(j).frames==CurrentFrame))),PixelsPerLine);
                    ceny=min(max(1,round(tempSchnitz(j).ceny(tempSchnitz(j).frames==CurrentFrame))),LinesPerFrame);
                    tempSchnitz(j).Fluo(tempSchnitz(j).frames==CurrentFrame,:,ChN) = single(convImage(ceny,cenx,:));
                    
                    
            end
            
        end
        schnitzcells = tempSchnitz;
        
        
        close(h);
    end
else
    error('Input channel not recognized. Check correct definition in MovieDatabase. Input channels should use the :input notation.');
end

clear movieMat;

if saveFlag
    if whos(var2str(schnitzcells)).bytes < 2E9
        save(schnitzPath, 'schnitzcells', '-v6');
    else
        save(schnitzPath, 'schnitzcells', '-v7.3', '-nocompression');
    end
    
end

end