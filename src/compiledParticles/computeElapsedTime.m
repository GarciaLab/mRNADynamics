function elapsedTime_min = computeElapsedTime(FrameInfo)

%Get the actual time corresponding to each frame

elapsedTime_min = [];
numFrames = length(FrameInfo);

for f=1:numFrames
     if ~isPrinceton2Photon
         elapsedTime_min(f)=FrameInfo(f).Time-FrameInfo(1).Time;
     else
        elapsedTime_min(f)=etime(datevec(FrameInfo(f).TimeString),datevec(FrameInfo(1).TimeString));
     end
end

 %time was stored in seconds in FrameInfo
elapsedTime_min=elapsedTime_min/60; 



end