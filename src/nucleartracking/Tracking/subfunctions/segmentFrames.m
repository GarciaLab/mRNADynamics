function xy = segmentFrames(FrameInfo,hisMat,firstFrame,lastFrame,nucleusDiameter,...
    embryoMask, varargin)

    nFrames = lastFrame-firstFrame+1;
    frameNum = firstFrame:lastFrame;
    xy = cell(lastFrame-firstFrame+1,1);

    if nargin > 6
        % Added/adapted by GM 1/9/19 to include mutlithresh option 
        useMultithresh=false;
        k=1;

        while k<=length(varargin)
            if strcmpi(varargin{k},'useMultithresh')
                useMultithresh = true;

                warning('Using multithresh function to threshold nuclei during segmentation (segment frames).')
            end
            k=k+1;
        end
         % Added/adapted by GM 1/9/19 to include mutlithresh option 
        
    end 
    
    for j = 1:nFrames 
         if ~useMultithresh
            [xy{j}, ~] = findNuclei(FrameInfo, hisMat(:,:, frameNum(j)), nucleusDiameter, embryoMask);
         else
            [xy{j}, ~] = findNuclei(FrameInfo, hisMat(:,:, frameNum(j)), nucleusDiameter, embryoMask, 'useMultithresh');
        end
    end
    
end