function doFrameSkipping(SkipFrames, FrameInfo, OutputFolder)
  
%Skipping frames?
  if ~isempty(SkipFrames)
    %Filter FrameInfo
    FrameFilter=ones(size(FrameInfo));
    FrameFilter(SkipFrames)=0;
    FrameInfo=FrameInfo(logical(FrameFilter));
    
    %Rename all Histone channel images
    %Find all the Histone files
    D=dir([OutputFolder,filesep,'*-His*.tif']);
    %Delete the skipped files
    for i=SkipFrames
        delete([OutputFolder,filesep,D(i).name]);
    end
    %Rename all remaining files
    D=dir([OutputFolder,filesep,'*-His*.tif']);
    for i=1:length(D)
        if ~strcmp([OutputFolder,filesep,D(i).name],...
                [OutputFolder,filesep,D(i).name(1:end-7),iIndex(i,3),'.tif'])
            movefile([OutputFolder,filesep,D(i).name],...
                [OutputFolder,filesep,D(i).name(1:end-7),iIndex(i,3),'.tif'])
        end
    end
    
    %Rename all coat protein channel images
    %Find all the coat protein files
    D=dir([OutputFolder,filesep,'*_z01.tif']);
    %Delete the skipped files
    for i=SkipFrames
        D2=dir([OutputFolder,filesep,D(i).name(1:end-6),'*.tif']);
        for j=1:length(D2)
            delete([OutputFolder,filesep,D2(j).name])
        end
    end
    %Rename all remaining files
    D=dir([OutputFolder,filesep,'*_z01.tif']);
    for i=1:length(D)
        if ~strcmp([OutputFolder,filesep,D(i).name],...
                    [OutputFolder,filesep,D(i).name(1:end-11),iIndex(i,3),'_z01.tif'])
            D2=dir([OutputFolder,filesep,D(i).name(1:end-6),'*.tif']);
            for j=1:length(D2)
                movefile([OutputFolder,filesep,D2(j).name],...
                    [OutputFolder,filesep,D2(j).name(1:end-11),iIndex(i,3),D2(j).name(end-7:end)])
            end
        end
    end
  end
end