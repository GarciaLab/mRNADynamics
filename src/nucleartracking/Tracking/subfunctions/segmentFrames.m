function xy = segmentFrames(FrameInfo,names,firstFrame,lastFrame,nucleusDiameter, embryoMask, varargin)

    update_waitbar = false;
    if nargin > 6
        % Added/adapted by GM 1/9/19 to include mutlithresh option 
        useMultithresh=false;
        k=1;

        while k<=length(varargin)
            if strcmpi(varargin{k},'useMultithresh')
                useMultithresh = true;

                warning('Using multithresh function to threshold nuclei during segmentation (segment frames).')
            elseif k == 1 % This just leaves what was done before edits to include Multithresh option 
                try
                    h_waitbar_segmentation = varargin{1};
                    update_waitbar = true;
                catch
                    %do nothing
                end
            end
            k=k+1;
        end
         % Added/adapted by GM 1/9/19 to include mutlithresh option 
        
    end
    nFrames = lastFrame-firstFrame+1;
    frameNum = firstFrame:lastFrame;
    xy = cell(lastFrame-firstFrame+1,1);
    
    update_waitbar = false;
    
    if useMultithresh
        
    end
    %parfor j = 1:nFrames
    for j = 1:nFrames   
        if ~useMultithresh
            [xy{j}, ~] = findNuclei(FrameInfo,names, frameNum(j), nucleusDiameter, embryoMask);
        else
            [xy{j}, ~] = findNuclei(FrameInfo,names, frameNum(j), nucleusDiameter, embryoMask, 'useMultithresh');
        end

        if update_waitbar
            waitbar(j/nFrames,h_waitbar_segmentation, ['Segmentation progress : ', num2str(j), ' frames processed out of ', num2str(nFrames)]);
        end
    end
    
end