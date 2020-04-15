function exportTifStackToAvi(fileIn, fileOut, varargin)

% AR 4/2020

%fileIn can be either a path to a file or a .mat  imagestack

load('ReferenceHist.mat', 'ReferenceHist');

if ischar(fileIn)
    stack = imreadStack(fileIn);
else
    stack = fileIn;
end

g = @(y) arrayfun(@(x) histeq(mat2gray(y), ReferenceHist), 3, 'UniformOutput',false);
tempStack = g(stack); 
stack = uint8(round(256*tempStack{1}));
nFrames = size(stack, 3);

outputVideo = VideoWriter(fileOut, 'Grayscale AVI');
outputVideo.FrameRate = 30;
open(outputVideo)
%Loop through the image sequence, load each image, and then write it to the video.

newStack = zeros(size(stack, 1), size(stack, 2), 1, size(stack, 3), 'uint8');


for f= 1:nFrames
    newStack(:, :, 1, f) = stack(:, :, f);
end



newStack = histeq(newStack, ReferenceHist);
writeVideo(outputVideo, newStack);
% Finalize the video file.
close(outputVideo)


end