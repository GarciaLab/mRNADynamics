function [CurrentZ,ManualZFlag] = changeZSlice(NewZ, ZSlices)
%CHANGESLICE Summary of this function goes here
%   Detailed explanation goes here
    CurrentZ = min(max(1, NewZ), ZSlices);
    ManualZFlag = 1;
end

