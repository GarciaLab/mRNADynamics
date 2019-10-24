function  [nuclearFluo,numNuclei] = integrateNuclearFluo(image,displayFigures)

% next steps is to add options for rmin rmax used here or have that be the
% options used in the extractNuclearFluoDosage orrrrr have this code run
% through different ranges of radii and pick the best range base on the
% integration....?
warning('off')
imageSize = size(image);
threshold1 = 50;
rmin = 4;
rmax = 30;
colorOfCirlces = [0.9290 0.6940 0.1250];


tempImage = threshold1.*uint16(image<=threshold1) + image;
%tempImage = imgaussfilt(image,'FilterSize',round(mean([rmin,rmax])));
[centersDark,radiiDark] = imfindcircles(tempImage,[rmin rmax],'ObjectPolarity','dark');

if displayFigures == 1
    figure('Position',[19 56 1217 562])
    subplot(1,2,1)
    imshow(image,[]);
%     imageTemp = imsharpen(image);
    [centersDark,radiiDark] = imfindcircles(image,[rmin rmax],'ObjectPolarity','dark');
    hold on
    viscircles(centersDark,radiiDark,'Color',colorOfCirlces);
    title(['Original : ',num2str(length(radiiDark))])
    
    subplot(1,2,2)
%     tempImageTemp = imsharpen(tempImage);
    imshow(tempImage,[]);
    hold on
    viscircles(centersDark,radiiDark,'Color',colorOfCirlces);
    title(['Threshold : ',num2str(length(radiiDark))])
    
    % code for histogram of intensities
elseif isequal(displayFigures,'histogram')
    figure('Position',[19 56 1217 562])
    subplot(1,2,1)
    imshow(image,[]);
%     imageTemp = imsharpen(image);
    [centersDark,radiiDark] = imfindcircles(image,[rmin rmax],'ObjectPolarity','dark');
    hold on
    viscircles(centersDark,radiiDark,'Color',colorOfCirlces);
    title(['Original : ',num2str(length(radiiDark))])
    
    subplot(1,2,2)
%     h = histogram(image);
%     bins = h.Values;
%     binCenters = h.BinEdges + h.BinWidth/2;
%     if length(binCenters) > 2
%         [peaks,locations]= findpeaks(bins,binCenters(2:length(binCenters)));
%         hold on
%         for i = 1:length(locations)
%             xValue = locations(i);
%             xValues = [xValue xValue];
%             plot(xValues,[0 10^5])
%         end
%     end
%     xlabel('Intensity')
%     ylabel('Frequency')
    [B,L] = bwboundaries(tempImage);%,'noholes');
    imshow(label2rgb(L, @jet, [.5 .5 .5]))
    hold on
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
    end
    
    %sigh what did past Emma do? :'( 
    
end

nuclearMask = makeMask(centersDark,radiiDark,imageSize);
% imshow(nuclearMask,[]);
nuclearFluo = sum(sum(image.*nuclearMask));
numNuclei = length(radiiDark);
end