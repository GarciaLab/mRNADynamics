function exportTifStackToAvi(fileIn, fileOut)


fileIn = 'X:\Armando\tr2dProject\RAW.tif';
fileOut = ['X:\Armando', filesep, 'hisVideo.avi'];

load('ReferenceHist.mat', 'ReferenceHist');
stack = imreadStack(fileIn);
stack = uint8(stack);
nFrames = size(stack, 3); 

file
outputVideo = VideoWriter(fileOut, 'Grayscale AVI');
outputVideo.FrameRate = 30;
open(outputVideo)
%Loop through the image sequence, load each image, and then write it to the video.
newStack = zeros(size(stack, 1), size(stack, 2), 1, size(stack, 3), 'uint8');
for ii = 1:nFrames
   newStack(:, :, 1, ii) = stack(:, :, ii);
end

newStack = histeq(newStack, ReferenceHist);
writeVideo(outputVideo, newStack);
% Finalize the video file.
close(outputVideo)


end