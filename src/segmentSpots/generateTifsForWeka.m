% Generates the TIF stacks necessary for doing Weka classification.
% Recommended to run this before making a new classifier.
function generateTifsForWeka(Prefix, PreProcPath, numFrames,...
    nCh,coatChannel, zSize, initialFrame)

mm = false;


movieChFile = [PreProcPath, filesep, Prefix, filesep, Prefix,'_movieMatCh', num2str(coatChannel), '.mat'];
movieFile = [PreProcPath, filesep, Prefix, filesep, Prefix,'_movieMat.mat'];
if exist(movieFile, 'file') || exist(movieChFile, 'file')
    movieMat = loadMovieMat(Prefix, 'chRange', coatChannel);
    mm=true;
end

% Create stacks subfolder
stacksPath = [PreProcPath, filesep, Prefix, filesep, 'stacks'];
mkdir(stacksPath);

% Sets waitbar for the whole tif generation process
tifStacksWaitbar = waitbar(0, 'Making .tif stacks for Weka classification');

% Make requisite TIF stacks for classification
for channelIndex = 1:nCh
    
    nameSuffix = ['_ch', iIndex(coatChannel, 2)];
    
    currentFrameWaitbar = waitbar(0, ['Making ch0', num2str(channelIndex), ' .tif stacks for Weka classification']);
    
    for currentFrame = initialFrame:numFrames
        currentFrameWaitbar = waitbar(currentFrame / numFrames, currentFrameWaitbar);
        set(currentFrameWaitbar, 'units', 'normalized', 'position', [0.4, .15, .25, .1]);
        rawStackName = [stacksPath, filesep, iIndex(currentFrame, 3), nameSuffix, '.tif'];
        
        %Don't write new stacks if they're already made.
        %2018-08-22 MT: Now takes into account 1 vs 2 spot channels
        %when determining if you've already made the stacks
        if length(dir([stacksPath, filesep, '*.tif'])) ~= numFrames * nCh
            rawStackArray = [];
            
            if ~mm
                
                for i = 1:zSize
                    fileName = [PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(currentFrame, 3), '_z', iIndex(i, 2), ...
                        nameSuffix, '.tif'];
                    rawStackArray(:, :, i) = imread(fileName);
                end
                
            else
                
                if numel(size(movieMat))==5
                    rawStackArray = squeeze(movieMat(:, :, :, currentFrame, coatChannel));
                elseif numel(size(movieMat))==4
                    rawStackArray = movieMat(:, :, :, currentFrame);
                end
                
            end
            
            if max(rawStackArray(:)) < 255 %max uint8 value
                rawStackArray = uint8(rawStackArray);
            else
                rawStackArray = uint16(rawStackArray);
            end
            
            imwrite(rawStackArray(:, :, 1), rawStackName);
            
            for z = 2:size(rawStackArray, 3)
                
                imwrite(rawStackArray(:, :, z), rawStackName, 'WriteMode', 'append');
                
            end
            
            clear rawStackArray;

        end
        
    end
    
    close(currentFrameWaitbar);
    waitbar(channelIndex / nCh, tifStacksWaitbar);
    
end

close(tifStacksWaitbar);

end
