function proj = zProjections(Prefix, currentChannel, currentFrame, zSlices, nDigits,DropboxFolder,PreProcPath, FrameInfo, projType)
% zProjections(Prefix, currentFrame, zSlices, nDigits,DropboxFolder,PreProcPath)
%
% DESCRIPTION
% This function will calculate the max time projection of the movie.
% By default it will make the max time projection using the max z
% projections.
%
% ARGUEMENTS
% Prefix: Prefix of the data set to analyze
% currentChannel: The current channel of the spots. It is assumed that this
%                 is a number value.
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
% Last Updated: 11/12/2017
%
% Documented by: Emma Luu (emma_luu@berkeley.edu)

%% Making Projection
% This is to store all the z stacks into a 3D matrix. 

if strcmpi(projType, 'median')
    Images =zeros(FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine, zSlices);

    for currentZ = 2:zSlices-1
        Images(:,:,currentZ) = imread([PreProcPath,filesep,Prefix,filesep,...
            Prefix,'_',iIndex(currentFrame,nDigits),'_z',iIndex(currentZ,2),'_ch0', num2str(currentChannel) ,'.tif']);
    end
    proj = median(Images,3);
elseif strcmpi(projType, 'max')
    im = zeros(FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine, 'uint16');
    proj = im;
    for currentZ = 2:zSlices-1
        im = imread([PreProcPath,filesep,Prefix,filesep,...
                Prefix,'_',iIndex(currentFrame,nDigits),'_z',iIndex(currentZ,2),'_ch0', num2str(currentChannel) ,'.tif']);
        proj = max(proj, im);
    end
end

end            