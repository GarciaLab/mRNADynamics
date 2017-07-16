function [maxProj, medianProj] = zProjections(Prefix, currentFrame, zSlices, nDigits,DropboxFolder,PreProcPath)
% zProjections(Prefix, currentFrame, zSlices, nDigits,DropboxFolder,PreProcPath)
%
% DESCRIPTION
% This function will calculate the max time projection of the movie.
% By default it will make the max time projection using the max z
% projections.
%
% ARGUEMENTS
% Prefix: Prefix of the data set to analyze
% currentFrame: The index of the z stack for the projection
% zSlices: The number of zSlices in the stack (currentZ goes from
%          1:zSlices)
% nDigits: This is 3 or 4 depending on whether totalFrames < 1E3 or 1E4
%          respectively.
% DropboxFolder: Dropbox folder path
% PreProcPath: Preprocessed data path
%
% OPTIONS
% none
%
% OUTPUT
% This will return the max and median z projections of the current z stack.
%
% Author (contact): Emma Luu (emma_luu@berkeley.edu)
% Created: 06/06/2017
% Last Updated: 06/22/2017
%
% Documented by: Emma Luu (emma_luu@berkeley.edu)

%% Information about the Folders and Frames 
% [~,~,DropboxFolder,~,PreProcPath]=...
%     DetermineLocalFolders(Prefix);
DataFolder=[DropboxFolder,filesep,Prefix];
load([DataFolder, filesep, 'FrameInfo.mat']);
FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];

%% Making Projection
% This is to store all the z stacks into a 3D matrix. 
Images =zeros(FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine, zSlices);

for currentZ = 1:zSlices
    Images(:,:,currentZ) = imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
        FilePrefix,iIndex(currentFrame,nDigits),'_z',iIndex(currentZ,2),'.tif']);
end
maxProj = max(Images,[],3);
medianProj = median(Images,3);
end            