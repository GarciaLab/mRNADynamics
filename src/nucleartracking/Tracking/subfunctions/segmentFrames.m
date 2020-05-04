function xy = segmentFrames(FrameInfo,hisMat,firstFrame,lastFrame,nucleusDiameter, embryoMask, varargin)


    nFrames = lastFrame-firstFrame+1;
    frameNum = firstFrame:lastFrame;
    xy = cell(lastFrame-firstFrame+1,1);
    
    
    parfor j = 1:nFrames %for j = 1:nFrames
        [xy{j}, ~] = findNuclei(FrameInfo,hisMat, frameNum(j), nucleusDiameter, embryoMask);

    end
    
end