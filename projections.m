function [maxProj, medianProj] = projections(Prefix, CurrentFrame, ZSlices, NDigits)
% This function will return the max and median projection of the frame
% (along the z axis). It must be given Prefic, CurrentFrame, TotalFrames,
% ZSlixes, NDigits (which are all variables from CheckParticleTracking). 

Images = []; % This is to store all the z stacks into one 3D matrix. 

[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);
DataFolder=[DropboxFolder,filesep,Prefix];

FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];

%Only use this if you aren't using maxZProjectionTest as a function...
% CurrentFrame = 22;
% NDigits = 3;
%Now get the actual folders
% [SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
%     DetermineLocalFolders(FilePrefix(1:end-1));

% if exist([DataFolder,filesep,'FrameInfo.mat'])
%     load([DataFolder,filesep,'FrameInfo.mat'])
% else
%     disp('Error')
% end
%
% ZSlices=FrameInfo(1).NumberSlices+2; %Note that the blank slices are included

for CurrentZ = 1:ZSlices
    Images(:,:,CurrentZ) = imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
        FilePrefix,iIndex(CurrentFrame,NDigits),'_z',iIndex(CurrentZ,2),'.tif']);
end
maxProj = max(Images,[],3);
medianProj = median(Images,3);
end
            
            
            