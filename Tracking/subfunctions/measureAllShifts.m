function [ shifts ] = measureAllShifts( names, varargin )
%MEASUREALLSHIFTS Summary of this function goes here
%   Detailed explanation goes here

if nargin > 1
    h_waitbar = varargin{1};
end

shifts = zeros(numel(names)-1,2);

im2 = imread(names{1});

for i = 2:numel(names)
    
    im1 = im2;
    im2 = imread(names{i});
    
    [shifts(i-1,1), shifts(i-1,2)] = measure_shift(im1,im2);
    waitbar(0.35+(i-1)/numel(names)*0.45,h_waitbar)

end

