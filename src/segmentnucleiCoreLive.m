
%%                              segmentnucleiCoreLive.m
% Overview:
% This code segments the nuclei in a confocal projection of an early
% drosophila embryo to find the cores of the nuclei. It finds the cores of the nuclei using a laplacian of
% gaussian filter and then segments these using a watershed algorithm to
% split joined cores. Output is a binary image of the cores of the nuclei.
%
% Input:
%
% I -  fluorescent image of nuclei
% sizefilt -  size of the initial filter to roughly find nuclei (10- 30)
% radiusfilt -  radius of the initial filter to roughly find nuclei  (5-10)
% coursefilterthreshold - threshold for segmenting course filter (if choose
% 0 threshold will be chosen automatically using Otsu's method) 
%
% Output:
%
% NucCourseFiltT - bw image of core nuclei
%
% Comments:
% The output is used as input for other functions that do a more complete
% segmentation of the nuclei


function NucCourseFiltT=segmentnucleiCoreLive(I,sizefilt,radiusfilt,coursefilterthreshold)

troubleshooting=0;

%%%%%% Finding core nuclear pixels %%%%%%%%%%%

H = -fspecial('log',sizefilt,radiusfilt);                      % Defining the first filter kernel
II=single(I);                                                  % Using single format speeds things up
NucCourseFilt = imfilter(II,H,'replicate');                    % Apply first coarse filter

NucCourseFilt = (NucCourseFilt./max(NucCourseFilt(:)));        % Scale filtered result
NucCourseFilt=NucCourseFilt.*(NucCourseFilt>0);                % Remove negative pixles.

if troubleshooting
imshowbig(NucCourseFilt)
end

if coursefilterthreshold==0
NucCourseFiltT=im2bw(NucCourseFilt,graythresh(NucCourseFilt)); % Threshold the filtered nuclear image with automatic thresholding
else
NucCourseFiltT=im2bw(NucCourseFilt,coursefilterthreshold);     % Threshold the filtered nuclear image user chosen threshold
end

NucSplit = watershed(-(NucCourseFilt));                        % Watershed the filtered image
NucCourseFiltT=logical(NucCourseFiltT).*logical(NucSplit);     % Segmented core pixels



%%%%%%%%