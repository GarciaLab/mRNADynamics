function [CurrentFrame, ManualZFlag] = changeFrame(NewFrame, MaxFrame)
%CHANGEFRAME 
%   Change the frame of a particle in CheckParticleTracking

    MinFrame = 1;
    CurrentFrame = max(NewFrame, MinFrame);
    CurrentFrame = min(CurrentFrame, MaxFrame);
    ManualZFlag = 0;
    
end

