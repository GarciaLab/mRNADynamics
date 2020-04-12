function obj = setupSystemObjects(hisVideoFile)
%
% This is a subfunction for tracknuclei_computervision that
%creates the computer vision object

%AR 4/2020

obj.reader = VideoReader(hisVideoFile);
obj.maskPlayer = vision.VideoPlayer('Position', [255 205 415 300]);
obj.videoPlayer = vision.VideoPlayer('Position', [726 205 415 300]);

obj.detector = vision.ForegroundDetector('NumGaussians', 8, ...
    'NumTrainingFrames', 40, 'MinimumBackgroundRatio', 0.4);

obj.blobAnalyser = vision.BlobAnalysis('BoundingBoxOutputPort', true, ...
    'AreaOutputPort', true, 'CentroidOutputPort', true, ...
    'MinimumBlobArea', 400);

end