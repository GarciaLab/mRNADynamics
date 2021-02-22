function [MeanCyto, SDCyto, MaxCyto, MedianCyto] =...
    ...
    getCytoplasmStatistics(...
    ...
    APExperiment, HistoneChannel, Prefix, numFrames, PreProcPath,...
    FrameInfo, NChannels)

 %Information about the cytoplasm
    %If the nuclear masks are present then use them. Otherwise just calculate
    %the median of the images as a function of time
    
    MeanCyto=[];
    SDCyto=[];
    MaxCyto=[];
    MedianCyto=[];
    
    liveExperiment = LiveExperiment(Prefix);

    if NChannels==1
    
        if HistoneChannel && APExperiment
            [MeanCyto,SDCyto,MedianCyto,MaxCyto]=GetCytoMCP(Prefix);
        else
            h=waitbar(0,'Calculating the median cyto intentisy');
            for frame=1:numFrames
                waitbar(frame/numFrames,h)
                Image = getMovieFrame(liveExperiment, frame, spotChannels);
                ImageMax=max(Image,[],3);
                MedianCyto(frame)=median(double(ImageMax(:)));
            end
            close(h)
        end
    else
   %?
    end


end