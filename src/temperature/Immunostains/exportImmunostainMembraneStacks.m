function exportImmunostainMembraneStacks(imageMat,Prefix)
%%

NEmbryos = size(imageMat, 4);

liveExperiment = LiveExperiment(Prefix);
PreProcFolder = liveExperiment.preFolder;

%Copy the data
waitbarFigure = waitbar(0, ['Extracting ' imagingModality ' images']);


%Counter for number of frames
numberOfFrames = 1;
totalFrames = sum(NReplicates);
NRepsPerEmbryo = max(NReplicates);

ySize = size(AllImages{1}{1,1}, 1);
xSize = size(AllImages{1}{1,1}, 2);
BlankImage = zeros(ySize, xSize, moviePrecision);

hisMat = zeros(ySize, xSize, NRepsPerEmbryo,NEmbryos, hisPrecision);

topZSlice = min(NSlices);

% loop through each series
for embryoIndex = 1:NEmbryos
    for repIndex = 1:NReplicates(embryoIndex)
        waitbar(numberOfFrames/totalFrames, waitbarFigure)
        for channelIndex = 1:NChannels
            
            NameSuffix = ['_ch',iIndex(channelIndex,2)];
            
            NewName = [Prefix, '_Position',iIndex(embryoIndex,3),'_', iIndex(repIndex,3),...
                NameSuffix, '.tif'];
            
            % write bottom slice (always a blank padding image)
            imwrite(BlankImage, [PreProcFolder, filesep, NewName]);
            
            %Copy the rest of the images
            slicesCounter = 1;
            firstImageIndex = (repIndex-1)*NSlices(embryoIndex)*NChannels+channelIndex;
            lastImageIndex = (repIndex-1)*NSlices(embryoIndex)*NChannels+(NSlices(embryoIndex)-1)*NChannels+channelIndex;
            for imageIndex = firstImageIndex:NChannels:lastImageIndex
                
                if imageIndex == firstImageIndex
                    imwrite(AllImages{embryoIndex,1}{imageIndex,1}, [PreProcFolder, filesep, NewName]);
                elseif slicesCounter <= topZSlice
                    % if zPadding, it will process all images (because topZSlice would be max(NSlices)
                    % if no zPadding, it will process images rounding down to the series with least
                    % zSlices, because topZSlice would be min(NSlices)
                    
                    imwrite(AllImages{embryoIndex,1}{imageIndex,1},...
                        [PreProcFolder, filesep, NewName], 'WriteMode', 'append');
                    slicesCounter = slicesCounter + 1;
                    
                end
                
            end
            
        end
        
        
    end
end



