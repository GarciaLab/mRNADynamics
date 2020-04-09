function obj = setupSystemObjects()

stack = imreadStack('X:\Armando\tr2dProject\RAW.tif');
stack = uint8(stack);
nFrames = size(stack, 3); 

outputVideo = VideoWriter(['X:\Armando', filesep, 'hisVideo.avi']);
outputVideo.FrameRate = 15;
open(outputVideo)
%Loop through the image sequence, load each image, and then write it to the video.
for ii = 1:nFrames
   writeVideo(outputVideo, stack(:, :, ii));
end
% Finalize the video file.
close(outputVideo)

obj.reader = VideoReader(['X:\Armando', filesep, 'hisVideo.avi']);
obj.maskPlayer = vision.VideoPlayer('Position', [740, 400, 700, 400]);
obj.videoPlayer = vision.VideoPlayer('Position', [20, 400, 700, 400]);

obj.detector = vision.ForegroundDetector('NumGaussians', 8, ...
            'NumTrainingFrames', 40, 'MinimumBackgroundRatio', 0.7);

obj.blobAnalyser = vision.BlobAnalysis('BoundingBoxOutputPort', true, ...
            'AreaOutputPort', true, 'CentroidOutputPort', true, ...
            'MinimumBlobArea', 400); 
        
        


