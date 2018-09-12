function xy = segmentFrames(Prefix,names,firstFrame,lastFrame,nucleusDiameter, embryoMask, varargin)

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

        [xy{j}, ~] = findNuclei(Prefix,names, frameNum(j), nucleusDiameter, embryoMask);

        if update_waitbar
            progress = findall(h_waitbar_segmentation,'type','patch');
            progress = get(progress,'XData');
            try
                progress = progress(2)/100;
                waitbar((progress*numel(names)+1)/numel(names),h_waitbar_segmentation,['Segmentation progress : ' num2str((progress*numel(names)+1)) ' processed out of ' num2str(numel(names))])
            catch
%                 warning('There is a problem with calling waitbar. Is this something with the Matlab version?')
            end
        end
    end
    
end