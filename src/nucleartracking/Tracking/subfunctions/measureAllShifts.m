function [ shifts ] = measureAllShifts( hisMat, varargin )
%MEASUREALLSHIFTS Summary of this function goes here
%   Detailed explanation goes here

if nargin > 1
    h_waitbar = varargin{1};
end
numberOfFrames = size(hisMat, 3);
shifts = zeros(numberOfFrames-1,2);

im2 = hisMat(:,:, 1);

for i = 2:numberOfFrames
    
    im1 = im2;
    im2 = hisMat(:, :, i);
    
    [shifts(i-1,1), shifts(i-1,2)] = measure_shift(im1,im2);
    waitbar(0.35+(i-1)/numberOfFrames*0.45,h_waitbar)

end

