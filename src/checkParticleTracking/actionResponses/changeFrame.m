function [CurrentFrame, ManualZFlag] = changeFrame(NewFrame, MaxFrame)
%CHANGEFRAME Summary of this function goes here
%   Detailed explanation goes here
    MinFrame = 1;
    CurrentFrame = max(NewFrame, MinFrame);
    CurrentFrame = min(CurrentFrame, MaxFrame);
    ManualZFlag = 0;
end

