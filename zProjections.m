function [maxProj, medianProj] = zProjections(Prefix, currentFrame, zSlices, NDigits)
% This function will return the max and median projection of the frame
% (along the z axis). It must be given Prefix, CurrentFrame,
% ZSlixes, NDigits (which are all variables from CheckParticleTracking). 

Images = []; % This is to store all the z stacks into one 3D matrix. 

[~,~,DropboxFolder,~,PreProcPath]=...
    DetermineLocalFolders(Prefix);
DataFolder=[DropboxFolder,filesep,Prefix];
f = load([DataFolder, filsep, FrameInfo]);

FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];

[PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
        FilePrefix,iIndex(currentFrame,NDigits),'_z',iIndex(currentZ,2),'.tif']);
    
Images = zeros(f.FrameInfo.LinesPerFrame, f.FrameInfo.PixelsPerLine, zSlices);
for currentZ = 1:zSlices
    Images(:,:,currentZ) = imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
        FilePrefix,iIndex(currentFrame,NDigits),'_z',iIndex(currentZ,2),'.tif']);
end
maxProj = max(Images,[],3);
medianProj = median(Images,3);
end            