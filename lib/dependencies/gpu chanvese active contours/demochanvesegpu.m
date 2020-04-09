% demo for chenvese function
% Copyright (c) 2009, 
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved  
% http://sites.google.com/site/rexstribeofimageprocessing/
%%
%-- Chan & Vese method on gray and color image
%   Find contours of objects
close all
clear all

tic
% I = imread('saturngray.png');
I = imread('Sample9SGS.gif');
size(I)
% Customized Mask
I = gsingle(I);
m = zeros(size(I,1),size(I,2),gsingle);
m(20:120,20:120) = 1;

% accelereyes:
% Change nIter to cover more of image (will take more time)
nIter = 500;

% call compute function (same for CPU and GPU)
seg = chenvese(I,m,nIter,0.1,'chan'); % ability on gray image
geval(seg);
toc
