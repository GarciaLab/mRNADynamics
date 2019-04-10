function xy = segmentFrames(FrameInfo,names,firstFrame,lastFrame,nucleusDiameter, embryoMask, varargin)

    update_waitbar = false;
    if nargin > 6
        try
            h_waitbar_segmentation = varargin{1};
            update_waitbar = true;
        catch
            %do nothing
        end
    end
    nFrames = lastFrame-firstFrame+1;
    frameNum = firstFrame:lastFrame;
    xy = cell(lastFrame-firstFrame+1,1);
    for j = 1:nFrames

        [xy{j}, ~] = findNuclei(FrameInfo,names, frameNum(j), nucleusDiameter, embryoMask);

        if update_waitbar
            waitbar(j/nFrames,h_waitbar_segmentation, ['Segmentation progress : ', num2str(j), ' frames processed out of ', num2str(nFrames)]);
        end
    end
    
end