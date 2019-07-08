function makeSurfImForDVCorrection(Prefix)

%give this function the path to the surface image max projection as input.
%the output is a classified image 


[SourcePath, ~, DefaultDropboxFolder, DropboxFolder, ~, ~,...
~, ~] = DetermineAllLocalFolders(Prefix);

%this is generated around line 400 of addparticleposition
surf=imread([DropboxFolder,filesep,Prefix,filesep,'DV',filesep,'surf_max.tif']);

%%
%let hough transform do its best
if size(surf, 1) == 1024
    surfbig = imresize(surf, 2);
end
surfbigfilt = imgaussfilt(surfbig, 2);
rows = size(surfbigfilt, 1);
cols = size(surfbigfilt, 2);
close all
imlog = surfbigfilt > 6;
Overlay = figure;
overlayAxes = axes(Overlay);
imshow(surfbigfilt, []);
% [centers,radii] = imfindcircles(surfbigfilt,[6 11], 'Sensitivity', .96);%
[centers,radii] = imfindcircles(imlog,[6 11], 'Sensitivity', .96);
imin = surfbigfilt < 3;
[centers2,radii2] = imfindcircles(imin,[6 11], 'Sensitivity', .94);
centers = vertcat(centers, centers2);
radii = vertcat(radii, radii2);
viscircles(centers, radii,'Color','b', 'LineWidth', .4);

rad = 9;

%%make corrections
cc=1;
% while (cc~='x')
%     
%     ct=waitforbuttonpress;
%     cc=get(Overlay,'currentcharacter');
%     cm=get(overlayAxes,'CurrentPoint');
%     
%     if (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
%         cc=1;
%         if cm(1,2)>0 & cm(1,1)>0 & cm(1,2) <= size(surfbigfilt, 1) & cm(1,1) <= size(surfbigfilt, 2)
%             centers = vertcat(centers, [cm(1, 1), cm(1, 2)]);
%             radii = vertcat(radii, 9);
%         end
%         imshow(surfbigfilt, [])
%         viscircles(centers, radii,'Color','b', 'LineWidth', .4);
%         
%     elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'alt'))
%         cc=1;
%         point = [cm(1,1), cm(1,2)];
%         if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=size(surfbigfilt, 1))& (cm(1,1)<=size(surfbigfilt, 2))
%             dif = centers - point;
%             [~, nearestNuc] = min(abs(dif));
%             centers(nearestNuc(2), :) = [];
%             radii(nearestNuc(2), :) = [];
%             imshow(surfbigfilt, [])
%             viscircles(centers, radii,'Color','b', 'LineWidth', .4);
%         end
%     elseif (ct~=0)&(cc=='9')    %Debug mode
%         keyboard
%     end
% end

%%
%make a binary image

classSurf = false(rows, cols);
for i = 1:size(centers, 1)
    rad = floor(radii(i));
    ctr =flip(round(centers(i, :)));
    for y = ctr(1) - rad:ctr(1) + rad
        for x = ctr(2) - rad:ctr(2) + rad
            if norm([y, x] - ctr) < rad
                classSurf(y, x) = true;
            end
        end
    end
end

classSurf = imresize(classSurf, .5);
classFig = figure();
imshow(classSurf);

imwrite(classSurf, [DropboxFolder,filesep,Prefix,filesep,'DV',filesep,'classSurf.tif']);
