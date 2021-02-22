function obj = setupSystemObjects(hisVideoFile)
%
% This is a subfunction for tracknuclei_computervision that
%creates the computer vision object

%AR 4/2020

obj.reader = VideoReader(hisVideoFile);
obj.maskPlayer = vision.VideoPlayer('Position', [255 205 415 300]);
obj.videoPlayer = vision.VideoPlayer('Position', [726 205 415 300]);


end