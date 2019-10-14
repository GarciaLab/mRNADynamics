function  [nuclearFluo,numNuclei] = integrateNuclearFluo(image)

% next steps is to add options for rmin rmax used here or have that be the
% options used in the extractNuclearFluoDosage orrrrr have this code run
% through different ranges of radii and pick the best range base on the
% integration....?
warning('off')
imageSize = size(image);
threshold1 = 50;
rmin = 4;
rmax = 20;


% figure('Position',[19 56 1217 562])
% subplot(1,2,1)
% imshow(image,[]);
% %imageTemp = imsharpen(image);
% [centersDark,radiiDark] = imfindcircles(image,[rmin rmax],'ObjectPolarity','dark');
% hold on 
% viscircles(centersDark,radiiDark,'Color',[0.9290 0.6940 0.1250])
% title(num2str(length(radiiDark)))
% 
% subplot(1,2,2)
tempImage = threshold1.*uint16(image<=threshold1) + image;
%tempImage2 = threshold2.*uint16(image>=threshold2) + tempImage;
%tempImageTemp = imsharpen(tempImage);
% imshow(tempImage,[]);

[centersDark,radiiDark] = imfindcircles(tempImage,[rmin rmax],'ObjectPolarity','dark');
% hold on 
% viscircles(centersDark,radiiDark,'Color',[0.9290 0.6940 0.1250])
% title(num2str(length(radiiDark)))

%histogram(image);
%xlabel('Intensity')
%ylabel('Frequency')
%imshow(image>=50,[])
%finalImage = image.*nuclearMask;
%imshow(finalImage,[]);

nuclearMask = makeMask(centersDark,radiiDark,imageSize);

nuclearFluo = sum(sum(image.*nuclearMask));
numNuclei = length(radiiDark);
end