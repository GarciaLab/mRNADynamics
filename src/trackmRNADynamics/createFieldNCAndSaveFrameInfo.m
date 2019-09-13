
function createFieldNCAndSaveFrameInfo(FrameInfo, OutputFolder, nc9, nc10, nc11, nc12, nc13, nc14)
% creating the field nc for FrameInfo
if exist([OutputFolder, filesep, 'FrameInfo.mat'], 'file')
    numberOfFrames = length(FrameInfo);
    
    for currentFrame = 1:numberOfFrames
        
        if currentFrame < nc9
            FrameInfo(currentFrame).nc = 8;
        elseif (currentFrame >= nc9) & (currentFrame < nc10)
            FrameInfo(currentFrame).nc = 9;
        elseif (currentFrame >= nc10) & (currentFrame < nc11)
            FrameInfo(currentFrame).nc = 10;
        elseif (currentFrame >= nc11) & (currentFrame <= nc12)
            FrameInfo(currentFrame).nc = 11;
        elseif (currentFrame >= nc12) & ( currentFrame <= nc13 | isnan(nc13) )
            FrameInfo(currentFrame).nc = 12;
        elseif (currentFrame >= nc13) &  ( currentFrame <= nc14 | isnan(nc14) )
            FrameInfo(currentFrame).nc = 13;
        elseif currentFrame >= nc14
            FrameInfo(currentFrame).nc = 14;
        end
        
    end
    
    save([OutputFolder, filesep, 'FrameInfo.mat'], 'FrameInfo')
else
    warning('Tried to save nc frame information, but could not since there is no FrameInfo.mat')
end

end