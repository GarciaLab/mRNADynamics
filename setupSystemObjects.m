function obj = setupSystemObjects()
% 
obj.reader = VideoReader(['X:\Armando', filesep, 'hisVideo.avi']);
obj.maskPlayer = vision.VideoPlayer('Position', [255 205 415 300]);
obj.videoPlayer = vision.VideoPlayer('Position', [726 205 415 300]);

obj.detector = vision.ForegroundDetector('NumGaussians', 8, ...
            'NumTrainingFrames', 40, 'MinimumBackgroundRatio', 0.4);

obj.blobAnalyser = vision.BlobAnalysis('BoundingBoxOutputPort', true, ...
            'AreaOutputPort', true, 'CentroidOutputPort', true, ...
            'MinimumBlobArea', 400); 
        
end