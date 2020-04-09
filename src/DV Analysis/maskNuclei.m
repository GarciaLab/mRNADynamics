function nuclearMask = maskNuclei(im, sigma)

sigma = round(sigma); 

% 1) smooth and amplify edges
h = fspecial('log', 3*sigma, sigma);
imlog = (imfilter(im, h, 'symmetric'));
% 2) otsu threshold the filtered image to binarize
imBin = ~imbinarize(imlog);
figure(1); imagesc(imlog);title('original'); colormap gray;

[darkCenters, darkRadii, metric] = imfindcircles(imlog, [20, 40], 'ObjectPolarity', 'dark',  'Method', 'TwoStage', 'Sensitivity', .85);
[brightCenters, brightRadii, metric] = imfindcircles(imlog, [20, 40], 'ObjectPolarity', 'bright',  'Method', 'TwoStage', 'Sensitivity', .85);
centers = vertcat(darkCenters, brightCenters); radii = vertcat(darkRadii, brightRadii);
viscircles(centers, radii)

figure(2); imagesc(~imBin); title('inverted');colormap gray;

%2.5) decide polarization and complement if dark
%logic- if the polarization is dark, then there will be 2 regions. if
%right, the number of regions will be roughly the number of nuclei. also, erosion will remove nuclei with the dark polarity.  
% 
ses = strel('disk', round(sigma/2));
imEr = imerode(imBin, ses); %imEr = bwmorph(imEr, 'thin');
imErComp = imerode(~imBin, ses);% imErComp = bwmorph(imErComp, 'thin');
% imEr = bwmorph(imBin, 'thin');
% imErComp = bwmorph(~imBin, 'thick');

figure(3); imagesc(imEr);title('original eroded'); colormap gray;

figure(4); imagesc(~imErComp); title('inverted eroded');colormap gray;
% figure
% imshow(brightEr, [])
% figure
% imshow(brightErComp, [])
% stats = regionprops(brightEr);
% statsComp = regionprops(brightErComp);
% 
% area = sum([stats.Area]);
% areaComp = sum([statsComp.Area]);
% nRegions = length(regionprops(brightEr));
% nRegionsComp = length(regionprops(brightErComp));
euler = regionprops(imEr, 'EulerNumber');
euler = [euler.EulerNumber];
eulerComp = regionprops(~imEr, 'EulerNumber');
eulerComp = [eulerComp.EulerNumber];
% imshow(imerode(brightBin, ses), [])

% if area > areaComp & nRegions < nRegionsComp
if max(abs(euler)) > max(abs(eulerComp)) %lots of holes if the image has dark polarization
    imBin = ~imBin;
end
% 
% figure
% imshow(brightBin, [])

% 3) clean up the binary image 
se = strel('disk',sigma);
brightopen = imopen(imBin,se);
% 4) compute the distance xform in preparation for watershed
D = -bwdist(~brightopen);
% 5) clean up the distance xform to prevent overwatershedding
mask = imextendedmin(D,3); %3 pixel gaps
D2 = imimposemin(D,mask);
% 6) last step- watershed to get nice, disconnected, convex nuclei
L = watershed(D2);
L(~brightopen) = 0;

nuclearMask = ~~L;

figure()
imshow(nuclearMask, [])

end