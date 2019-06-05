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

    %HG on 8/6/16: Why was this commented out? Did I do this?


    if NChannels==1
    
        if HistoneChannel && APExperiment
            [MeanCyto,SDCyto,MedianCyto,MaxCyto]=GetCytoMCP(Prefix);
        else
            h=waitbar(0,'Calculating the median cyto intentisy');
            for i=1:numFrames
                waitbar(i/numFrames,h)
                for j=1:FrameInfo(1).NumberSlices
                    % JP on 6/3/2019: I'm adding a hardcoded _ch01 to
                    % comply with naming convention, otherwise it fails. I
                    % hardcoded 01 since we are inside the if NChannels ==
                    % 1 block.
                    Image(:,:,j)=imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(i,3),'_z',iIndex(j,2),'_ch01','.tif']);
                end
                ImageMax=max(Image,[],3);
                MedianCyto(i)=median(double(ImageMax(:)));
            end
            close(h)
        end
    else
   %?
    end


end