function maskCytoplasmForWeka(Prefix, radiusScale)

% Prefix = '2019-11-26-2xDl_Venus_snaBAC_MCPmCherry_Leica_Zoom45_21uW14uW_01';
% radiusScale = 1.1;

liveExperiment = LiveExperiment(Prefix);

nCh = numel(liveExperiment.spotChannels);
channels = liveExperiment.Channels;
channelNames = cell(1,nCh);

xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
zDim = liveExperiment.zDim;
nFrames = liveExperiment.nFrames;

dropboxFolder = liveExperiment.userResultsFolder;
resultsFolder = liveExperiment.resultsFolder;
preProcFolder = liveExperiment.preFolder;
maskPreProcFolder = [preProcFolder filesep 'maskedImages\'];
mkdir(maskPreProcFolder);

Ellipses = getEllipses(liveExperiment);

nCh = 1;
% nFrames = 1;
rawImDir = cell(1,nCh);
for channel = 1:nCh
    rawImDir{channel} = dir([preProcFolder filesep '*ch01.tif']);
    
    h = waitbar(1/nFrames, ['Masking raw images for Channel 0' num2str(channel) '...']);
    for currFrame = 1:nFrames
        waitbar(currFrame/nFrames,h);
        %Get the raw image
        currImPath = [rawImDir{1,channel}(currFrame).folder, filesep, rawImDir{1,channel}(currFrame).name];
        %using bfopen is slow, but fits into 2 lines  - figure out a faster
        %way to do this with the Tiff class
        evalc('currIm = bfOpen3DVolume(currImPath);');   %using evalc to suppress fprint statement inside bfopen
        imStack = currIm{1,1}{1,1};     %this is the xDim x yDim x zDim image matrix

        %Make the cytoplasmic mask
        ellipseFrame = Ellipses{currFrame};
        nuclearMask = makeNuclearMask(ellipseFrame, [yDim xDim], radiusScale);
        nuclearMask = nuclearMask >= 1;  %overlapping nuclei are annotated with a value of 2, but I want a binarized mask   
%         figure(1)
%         imshow(nuclearMask,[])
        
        %Create the new file inside the maskedImages folder
        nameSuffix = ['_ch',iIndex(channel,2)];
        imageName = [Prefix, '_', iIndex(currFrame,3), '_mask', nameSuffix, '.tif'];
        paddedZDim = zDim + 2;
        
        for currZSlice = 1:paddedZDim
            imSlice = imStack(:,:,currZSlice);
            maskedSlice = double(imSlice) .* nuclearMask;
            
%             %show the masked z slice
%             figure(2)
%             imshow(maskedSlice,[])
%             pause(2)
            
            %save this masked z slice
            if currZSlice == 1
                imwrite(uint16(maskedSlice), [maskPreProcFolder, filesep, imageName]);
            else
                imwrite(uint16(maskedSlice),[maskPreProcFolder, filesep, imageName], 'WriteMode', 'append');
            end
        end
    end
    close(h)
    
end