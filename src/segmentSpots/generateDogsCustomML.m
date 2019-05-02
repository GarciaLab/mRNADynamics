% Generates differences of gaussians files that will be later used to run segment spots ML, requires a classifier.
function generateProbMapsCustomML(Prefix, ProcPath, MS2CodePath, PreProcPath, ExperimentType, coatChannel, zSize, numFrames, nCh, initialFrame, ignoreMemoryCheck, classifierPathCh1, classifierFolder)

cd(MS2CodePath);
OutputFolder1 = [ProcPath, filesep, Prefix, '_', filesep, 'dogs'];
mkdir(OutputFolder1)

stacksPath = [PreProcPath, filesep, Prefix, filesep, 'stacks'];

if isempty(classifierPathCh1)
    [classifierPathCh1, classifierFolder] = uigetfile([MS2CodePath, '/src/classifiers/*.model']);
end

if nCh == 2
    [classifierPathCh2, ~] = uigetfile([MS2CodePath, '/src/classifiers/*.model']);
end


filterMovieWaitbarFigure = waitbar(0, 'Running CustomML Classifier');

for channelIndex = 1:nCh
    
    if strcmpi(ExperimentType, 'inputoutput')
        nameSuffix = ['_ch', iIndex(coatChannel, 2)];
    else
        nameSuffix = ['_ch', iIndex(channelIndex, 2)];
    end
    
    currentFrameWaitbar = waitbar(0, ['Making ch0', num2str(channelIndex), ' .tif stacks for Weka classification']);
    
    for currentFrame = initialFrame:numFrames
        currentFrameWaitbar = waitbar(currentFrame / numFrames, currentFrameWaitbar);
        set(currentFrameWaitbar, 'units', 'normalized', 'position', [0.4, .15, .25, .1]);
        
        %generate instead of load
        rawStackName = [stacksPath, filesep, iIndex(currentFrame, 3), nameSuffix, '.tif'];
        
    end
    
    %
    pMap = [];
    
    for m = 1:2:zSize*2
        pMap(:, :, ceil(m / 2)) = pMapTemp(:, :, m); %the even images in the original array are negatives of the odds
    end
    
    clear pMapTemp;
    
    for i = 1:size(pMap, 3)
        p_name = ['prob', Prefix, '_', iIndex(currentFrame, 3), '_z', iIndex(i, 2), nameSuffix, '.tif'];
        imwrite(uint16(pMap(:, :, i)), [OutputFolder1, filesep, p_name])
    end
    
    clear pMap;
    
end

close(currentFrameWaitbar);
waitbar(channelIndex / nCh, filterMovieWaitbarFigure);
end

close(filterMovieWaitbarFigure);

end
